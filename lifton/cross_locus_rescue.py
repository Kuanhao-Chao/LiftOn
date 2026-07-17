"""Phase D — duplicate-safe, gene-level, CROSS-LOCUS miniprot rescue (opt-in).

Figure-4 deep-dive (benchmarks/compare/figure4_outliers.py + GFF3 forensics)
found the very-distant per-transcript deficit vs miniprot is TWO mechanisms:

  * SAME-LOCUS (chicken ~83% of regressed): Liftoff lifts a truncated/broken CDS
    but miniprot's clean model is at the SAME target locus -> captured by the
    3-way miniprot-only merge candidate (run_liftoff candidate-3).
  * CROSS-LOCUS paralog disagreement (zebrafish ~86%, arab->rice ~88%): Liftoff
    lifts to the SYNTENIC copy (a near-garbage stub for some variants) while
    miniprot aligns the protein cleanly to a DIFFERENT chromosome (the paralog /
    duplicate). ``has_valid_miniprot`` is False there (miniprot does not overlap
    the Liftoff locus), so candidate-3 never fires -- this residual is unaddressed.

Iteration 15 (Gap #2) tried to recover these by ADDING the miniprot model, which
produced DUPLICATE gene models (the weak Liftoff gene stayed). This pass instead
REPLACES: it drops every block of a weak reference gene from the (already
written) output GFF3 and appends a single strictly-better cross-locus miniprot
model, so the gene id appears exactly once -> duplicate-safe by construction.

Mechanics (runs AFTER the Step-7/8 writer + score.txt have closed, so it never
perturbs a default decision; the whole pass is gated behind
``args.cross_locus_rescue`` -> the default path never imports this module):

  1. Index score.txt -> per reference gene: best emitted protein identity + the
     set of target seqids its emitted model(s) sit on.
  2. Walk miniprot mRNAs; for a gene whose best emitted identity is WEAK
     (< MAX_LIFTOFF) and whose miniprot hit is on a DIFFERENT seqid (cross-locus),
     build + ORF-rescue + score the miniprot model. Keep it iff its identity beats
     the emitted best by >= MIN_GAIN. Dedup per reference gene (keep the best).
  3. Rewrite the output GFF3: stream it, dropping every block whose top-level
     feature resolves (via get_ID_base) to a replaced reference gene, then append
     the replacement models (tagged ``lifton_rescue=cross_locus``).

SYNTENY CAVEAT: a replacement moves the gene to miniprot's (possibly paralogous)
locus -- it trades synteny for protein identity. This is why the pass is opt-in
(``--miniprot-cross-locus-rescue`` / env ``LIFTON_CROSS_LOCUS_RESCUE``) and the
gain must clear a strict, regression-free A/B before any promotion.
"""
import copy
import io
import os
import re
import sys
from collections import ChainMap

from lifton import run_miniprot, lifton_utils, logger

_COPY_SUFFIX = re.compile(r"_\d+$")


def _isolated_ref_features(ref_features_dict, ref_gene_id):
    """Give an experimental replacement private copy-counter state."""
    if ref_gene_id not in ref_features_dict:
        return ref_features_dict
    return ChainMap(
        {ref_gene_id: copy.copy(ref_features_dict[ref_gene_id])},
        ref_features_dict,
    )


def _enabled(args):
    env = os.environ.get("LIFTON_CROSS_LOCUS_RESCUE")
    if env is not None:
        return env.strip().lower() not in ("", "0", "false", "no", "off")
    return bool(getattr(args, "cross_locus_rescue", False))


def _max_liftoff(args):
    """A gene is 'weak' (eligible for replacement) if its best emitted protein
    identity is below this. Env LIFTON_CROSS_LOCUS_MAX_LIFTOFF overrides."""
    try:
        return float(os.environ.get("LIFTON_CROSS_LOCUS_MAX_LIFTOFF",
                                    getattr(args, "cross_locus_max_liftoff", 0.5)))
    except (TypeError, ValueError):
        return 0.5


def _min_gain(args):
    """Replace only if the cross-locus miniprot identity beats the emitted best by
    at least this margin. Env LIFTON_CROSS_LOCUS_MIN_GAIN overrides."""
    try:
        return float(os.environ.get("LIFTON_CROSS_LOCUS_MIN_GAIN",
                                    getattr(args, "cross_locus_min_gain", 0.3)))
    except (TypeError, ValueError):
        return 0.3


def _parse_loc(loc):
    """'seqid:start-end' -> (seqid, start, end) ints, or None on parse failure."""
    try:
        seqid, span = loc.rsplit(":", 1)
        start, end = span.split("-", 1)
        return seqid, int(start), int(end)
    except (ValueError, IndexError):
        return None


def _index_emitted(score_path, ref_features_dict, ref_features_reverse_dict):
    """score.txt -> {ref_gene_id: {"best_pi": float, "intervals": [(seqid,start,end), ...]}}.

    score.txt cols: 0=transcript_id 1=liftoff 2=miniprot 3=lifton_dna
                    4=lifton_aa 5=annotation 6=status 7=seqid:start-end
    The lifted transcript id resolves to its reference gene via get_ID_base
    (strip a LiftOn copy suffix) -> ref_features_reverse_dict (ref_trans->ref_gene),
    the SAME namespace get_ref_ids_miniprot returns, so the join is exact.

    ``intervals`` (not a bare seqid set) is what lets the caller distinguish a
    genuinely SAME-locus miniprot hit (overlaps an already-emitted interval —
    candidate-3's job) from a cross-locus one that merely happens to sit on the
    SAME chromosome as a weak lift (two different genes on one big chromosome
    are NOT "the same locus" — a same-seqid-only check would wrongly treat any
    same-chromosome paralog as already-covered and never rescue it).
    """
    index = {}
    if not os.path.exists(score_path):
        return index
    with open(score_path) as fh:
        for line in fh:
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 8:
                continue
            tid, loc = cols[0], cols[7]
            try:
                pi = float(cols[4])
            except (ValueError, IndexError):
                continue
            # Resolve the lifted transcript id to its reference gene. Try the id
            # as-is, then strip a trailing LiftOn copy suffix (_N) and retry --
            # ref_features_reverse_dict is keyed by transcript, so get_ID_base
            # (which validates against the gene-keyed ref_features_dict) cannot
            # strip a transcript-level copy suffix here.
            ref_gene = ref_features_reverse_dict.get(tid)
            if ref_gene is None:
                ref_gene = ref_features_reverse_dict.get(_COPY_SUFFIX.sub("", tid))
            if ref_gene is None:
                continue
            parsed = _parse_loc(loc)
            slot = index.setdefault(ref_gene, {"best_pi": pi, "intervals": []})
            slot["best_pi"] = max(slot["best_pi"], pi)
            if parsed is not None:
                slot["intervals"].append(parsed)
    return index


def _overlaps_any(seqid, start, end, intervals):
    for iv_seqid, iv_start, iv_end in intervals:
        if iv_seqid == seqid and start <= iv_end and iv_start <= end:
            return True
    return False


def _build_replacement_text(mtrans, m_feature_db, ref_db, ref_gene_id, ref_trans_id,
                            tgt_fai, ref_proteins, ref_trans, tree_dict,
                            ref_features_dict, args):
    """Build + ORF-rescue + score a miniprot-only model; return
    (protein_identity, gff3_text) or (None, None) on failure / no protein.

    Mirrors the Iter-23 rescue build (run_miniprot.lifton_miniprot_with_ref_protein
    + orf_search_protein) but captures the GFF3 lines to a string (the writer is
    already closed) instead of streaming to fw.
    """
    if ref_trans_id not in ref_proteins or ref_trans_id not in ref_trans:
        return None, None
    candidate_ref_features = _isolated_ref_features(
        ref_features_dict, ref_gene_id,
    )
    lifton_gene, lifton_trans, transcript_id, lifton_status = \
        run_miniprot.lifton_miniprot_with_ref_protein(
            mtrans, m_feature_db, ref_db.db_connection, ref_gene_id, ref_trans_id,
            tgt_fai, ref_proteins, ref_trans, {}, candidate_ref_features, args)
    # Lifton_GENE.__init__ auto-assigns a copy-numbered id (e.g. "gene1_1") when
    # ref_features_dict[ref_gene_id].copy_num was already bumped by the WEAK
    # gene this replaces (Step 7 constructed one Lifton_GENE for it already).
    # This is a REPLACE, not an additional copy -- force the canonical
    # ref_gene_id back onto the gene AND its transcript's Parent linkage (set
    # during add_miniprot_transcript from the pre-override id), so the
    # appended block truly replaces the dropped one (same gene id) instead of
    # silently renaming it (which would make the weak gene id vanish and a
    # new "_1" id appear -- still duplicate-safe, but not a replace).
    lifton_gene.entry.attributes["ID"] = [ref_gene_id]
    lifton_gene.entry.id = ref_gene_id
    lifton_gene.entry.attributes.pop("extra_copy_number", None)
    for trans in lifton_gene.transcripts.values():
        if "Parent" in trans.entry.attributes:
            trans.entry.attributes["Parent"] = [ref_gene_id]
    lifton_gene.transcripts[transcript_id].entry.attributes[
        "lifton_rescue"] = ["cross_locus"]
    lifton_gene.orf_search_protein(lifton_trans.entry.id, ref_trans_id, tgt_fai,
                                   ref_proteins, ref_trans, lifton_status)
    pi = lifton_status.lifton_aa
    # Mirror the Iter-23 rescue's write-side tail (miniprot_rescue.py) so the
    # emitted mRNA carries protein_identity/dna_identity/status -- without
    # this, write_entry emits a bare mRNA row with no identity attributes.
    lifton_gene.add_lifton_gene_status_attrs("miniprot")
    lifton_gene.add_lifton_trans_status_attrs(transcript_id, lifton_status)
    buf = io.StringIO()
    # write_entry's coding/non-coding stats-counter branch indexes this dict by
    # a fixed key set ("coding"/"non-coding"/"other", lifton/io/feature_serializer
    # .py:22) -- a bare {} raises KeyError('coding'). Mirror lifton.py:757's shape
    # (a throwaway copy since we only need the write side effect, not the counts).
    _throwaway_stats = {"coding": {}, "non-coding": {}, "other": {}}
    written = lifton_gene.write_entry(buf, _throwaway_stats)
    if written is False or getattr(
            lifton_gene, "_serialization_failures", []):
        return None, None
    return pi, buf.getvalue()


def _rewrite_output(out_path, replace_set, replacement_texts, ref_features_dict):
    """Stream out_path -> tmp, dropping every block whose top-level feature
    resolves to a ref gene in replace_set, append the replacement texts, atomic
    rename. A 'block' = a top-level line (no Parent= attribute) + the following
    child lines until the next top-level line (LiftOn writes contiguous blocks).
    Returns (n_blocks_dropped, n_lines_dropped).
    """
    tmp = out_path + ".tmp_xlocus"
    dropping = False
    n_blocks = 0
    n_lines = 0
    with open(out_path) as src, open(tmp, "w") as dst:
        for line in src:
            if line.startswith("#") or not line.strip():
                dst.write(line)
                continue
            cols = line.split("\t")
            if len(cols) < 9:
                dst.write(line)
                continue
            attrs = cols[8]
            is_top_level = "Parent=" not in attrs
            if is_top_level:
                # decide this block's fate from its own (top-level) id
                fid = ""
                for kv in attrs.split(";"):
                    if kv.startswith("ID="):
                        fid = kv[3:].strip()
                        break
                ref_gene = lifton_utils.get_ID_base(fid, ref_features_dict)
                dropping = ref_gene in replace_set
                if dropping:
                    n_blocks += 1
            if dropping:
                n_lines += 1
                continue
            dst.write(line)
        # append the replacement models at the end
        for ref_gene in sorted(replacement_texts):
            dst.write(replacement_texts[ref_gene])
    os.replace(tmp, out_path)
    return n_blocks, n_lines


def cross_locus_rescue_pass(out_path, score_path, m_feature_db, ref_db, tree_dict,
                            tgt_fai, ref_proteins, ref_trans, ref_features_dict,
                            m_id_2_ref_id_trans_dict, ref_features_len_dict,
                            ref_trans_exon_num_dict, ref_features_reverse_dict,
                            args):
    """Replace weak DNA-lifted genes with strictly-better cross-locus miniprot
    models. Returns the number of genes replaced; inert (returns 0) when the
    flag is off or there is no miniprot evidence.
    """
    if not _enabled(args) or m_feature_db is None:
        return 0

    max_liftoff = _max_liftoff(args)
    min_gain = _min_gain(args)
    index = _index_emitted(score_path, ref_features_dict, ref_features_reverse_dict)
    if not index:
        return 0

    # ref_gene -> (best miniprot pi seen, gff3 text). Dedup keeps the best hit.
    chosen = {}
    chosen_pi = {}

    try:
        mtranscripts = sorted(
            m_feature_db.features_of_type('mRNA'),
            key=lambda m: (m.seqid, m.start, m.end, m.attributes["ID"][0]))
    except Exception as e:
        logger.log_error(f"cross-locus rescue: failed to enumerate mRNAs: {e}")
        return 0

    for mtrans in mtranscripts:
        try:
            mtrans_id = mtrans.attributes["ID"][0]
            ref_gene_id, ref_trans_id = lifton_utils.get_ref_ids_miniprot(
                ref_features_reverse_dict, mtrans_id, m_id_2_ref_id_trans_dict)
            if ref_gene_id is None or ref_trans_id is None:
                continue

            slot = index.get(ref_gene_id)
            if slot is None:                       # gene was never emitted -> Iter-23's job
                continue
            if slot["best_pi"] >= max_liftoff:     # emitted model is good enough
                continue
            if _overlaps_any(mtrans.seqid, mtrans.start, mtrans.end, slot["intervals"]):
                continue                            # genuinely overlaps -> candidate-3's job
            # processed-pseudogene filter (mirror process_miniprot / Iter-23)
            if (len(list(m_feature_db.children(mtrans, featuretype='CDS'))) == 1
                    and ref_trans_exon_num_dict.get(ref_trans_id, 0) > 1):
                continue

            pi, text = _build_replacement_text(
                mtrans, m_feature_db, ref_db, ref_gene_id, ref_trans_id, tgt_fai,
                ref_proteins, ref_trans, tree_dict, ref_features_dict, args)
            if pi is None or text is None:
                continue
            if pi <= slot["best_pi"] + min_gain:   # not strictly-enough better
                continue
            if pi > chosen_pi.get(ref_gene_id, -1.0):
                chosen_pi[ref_gene_id] = pi
                chosen[ref_gene_id] = text
        except Exception as e:
            logger.log_error(f"cross-locus rescue error ({mtrans.id}): {e}")

    if not chosen:
        return 0

    replace_set = set(chosen)
    n_blocks, n_lines = _rewrite_output(out_path, replace_set, chosen,
                                        ref_features_dict)
    sys.stderr.write(
        f"\n[LiftOn] cross-locus rescue: {len(chosen)} gene(s) replaced with a "
        f"better miniprot model (dropped {n_blocks} weak block(s); "
        f"weak<{max_liftoff:.2f}, gain>={min_gain:.2f}).\n")
    sys.stderr.flush()
    return len(chosen)
