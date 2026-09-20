"""Iteration 23 — clean separate-pass miniprot-only rescue.

The default LiftOn pipeline never falls back to a miniprot-only gene model for a
coding gene the DNA lift (Liftoff / Step 7) missed ENTIRELY and that Step 8 then
dropped for falling outside the tight ``-min_miniprot``/``-max_miniprot`` length
band. The Iteration-22 in-loop rescue (an ``elif`` inside
``run_miniprot.process_miniprot``) recovered those genes but was NOT
default-ready: building a ``Lifton_GENE`` adds it to the SHARED Step-8
suppression ``tree_dict`` (``lifton_class.Lifton_GENE.__init__``), so a rescue
emitted mid-loop suppressed a LATER overlapping candidate the default would have
emitted (the ``off ⊄ on`` swap), and a multi-hit ref transcript got a redundant
second model.

This module runs the rescue as a SEPARATE POST-LIFT PASS, AFTER the Step-7 and
Step-8 loops have fully closed. Because no rescue mutates ``tree_dict`` during
any default decision, the default Step-7+8 output is byte-identical between
flag-OFF and flag-ON, so ``on = default ∪ rescues ⊇ default = off`` (off ⊆ on by
construction). A ref-gene-id dedup set (the genes Step 7 + Step 8 already
emitted) plus the now-final suppression tree give 0-redundant emission. Mirrors
the shipped lifton2 design (``lifton2/lifton2/miniprot_rescue.py``).

The whole pass is gated behind ``args.miniprot_rescue``; the caller never imports
this module when the flag is OFF, so the default path is provably inert.
"""
import io
import os
import sys
import time
from collections import defaultdict

from intervaltree import Interval

from lifton import (align, coreutils, lifton_class, orf_completion,
                    run_miniprot, lifton_utils, logger)
from lifton.intervals import _make_interval
from lifton.locus_pipeline import DeferredStateJournal, commit_locus_delta


# Divergence-adaptive protein-identity floor (PROMOTED to default ON). The
# miniprot-only rescue (Iteration 23) used a FIXED floor (default 0.5). When the
# DNA lift placed only a small fraction of the genes miniprot ALSO found, the pair
# is divergent and the rescue is the dominant recall lever, so lower the floor
# toward ADAPT_FLOOR_MIN to admit more genuinely-missing genes; above ADAPT_R_HIGH
# recall (same/close-species) the floor stays at the user's base (rescue is
# marginal there -> keep precision). Mirrors lifton2's shipped design
# (lifton2/lifton2/miniprot_rescue.py). Thresholds env-overridable for the A/B.
ADAPT_FLOOR_MIN = float(os.environ.get("LIFTON_RESCUE_FLOOR_MIN") or 0.30)
ADAPT_R_LOW = float(os.environ.get("LIFTON_RESCUE_R_LOW") or 0.10)
ADAPT_R_HIGH = float(os.environ.get("LIFTON_RESCUE_R_HIGH") or 0.50)

# Protein-coverage gate (v1.0.12, sub-pass B). The length band above compares a
# miniprot model's genomic span with the reference gene's CDS span, and both
# spans include introns. Intron length scales with genome size, so across
# species the ratio says little about whether the hit is complete: on human to
# zebrafish, chicken and xenopus, 75-83 % of the missed coding genes miniprot
# finds at protein identity >= 0.5 fail that band, although their alignments
# cover the whole reference protein (median coverage 1.00). Sub-pass B
# reconsiders exactly the candidates the band rejected, gated instead on the
# fraction of the reference protein the hit aligns (miniprot's
# ``Target=<id> <start> <end>``), which does not depend on intron length.
# On by default since v1.0.12: 12/12 A/B cells (8-cell ladder, whole genomes)
# added genes with 0 lost, 0 duplicates and 0 regressions
# (benchmarks/compare/rescue_extension_ab.coverage_gate.md).
COVERAGE_GATE_DEFAULT = True
COVERAGE_MIN_DEFAULT = 0.8
# Coverage bounds the hit from below; this bounds it from above. A model whose
# CDS is more than 1.5 times the reference coding length is mostly inserted
# sequence, which the adaptive PI floor (as low as 0.30) would not reject alone.
CDS_LENGTH_MAX_RATIO = 1.5


# Isoform-aware rescue (v1.0.12). The rescue deduplicates by reference gene, so
# it used to emit one transcript per rescued gene; a human gene carries about six.
# With this on, once every rescue gene is placed, each one also receives the
# other reference transcripts of the same gene whose miniprot hits sit at its
# locus (same sequence and strand, overlapping the placed hit), each held to the
# same protein-identity floor. Placement decisions are made before any isoform
# is attached, so the set of rescued genes and their order do not change.
# On by default since v1.0.12 (benchmarks/compare/rescue_extension_ab.isoforms.md).
ISOFORMS_DEFAULT = True


def _rescue_isoforms_on(args):
    """Whether the isoform pass runs. Env ``LIFTON_RESCUE_ISOFORMS`` wins (1/0),
    then ``args.rescue_isoforms`` when set, else the module default."""
    env = os.environ.get("LIFTON_RESCUE_ISOFORMS")
    if env is not None:
        return env.strip().lower() not in ("0", "", "false", "no")
    value = getattr(args, "rescue_isoforms", None)
    return ISOFORMS_DEFAULT if value is None else bool(value)


class _Accepted:
    """A rescue gene that passed every gate, staged but not yet written."""

    __slots__ = ("lifton_gene", "mtrans", "ref_gene_id", "ref_trans_id",
                 "block", "staged_stats", "score_text")

    def __init__(self, lifton_gene, mtrans, ref_gene_id, ref_trans_id, block,
                 staged_stats, score_text):
        self.lifton_gene = lifton_gene
        self.mtrans = mtrans
        self.ref_gene_id = ref_gene_id
        self.ref_trans_id = ref_trans_id
        self.block = block
        self.staged_stats = staged_stats
        self.score_text = score_text


class _RescuePublisher:
    """Writes accepted genes at once, or holds them until :meth:`flush`.

    Holding changes only when the text is written, not what is decided: every
    placement effect (copy number, suppression interval, dedup) is committed at
    acceptance either way, and held genes are written in acceptance order.
    """

    def __init__(self, fw, fw_score, transcripts_stats_dict, defer):
        self.fw = fw
        self.fw_score = fw_score
        self.transcripts_stats_dict = transcripts_stats_dict
        self.defer = defer
        self.accepted = []

    def publish(self, accepted):
        if self.defer:
            self.accepted.append(accepted)
        else:
            self._write(accepted)

    def _write(self, accepted):
        self.fw.write(accepted.block)
        _merge_stats(self.transcripts_stats_dict, accepted.staged_stats)
        if accepted.score_text:
            self.fw_score.write(accepted.score_text)

    def flush(self):
        for accepted in self.accepted:
            self._write(accepted)
        self.accepted = []


def _coverage_gate_on(args):
    """Whether sub-pass B runs. Env ``LIFTON_RESCUE_COVERAGE_GATE`` wins (1/0),
    then ``args.coverage_rescue_gate`` when set, else the module default."""
    env = os.environ.get("LIFTON_RESCUE_COVERAGE_GATE")
    if env is not None:
        return env.strip().lower() not in ("0", "", "false", "no")
    value = getattr(args, "coverage_rescue_gate", None)
    return COVERAGE_GATE_DEFAULT if value is None else bool(value)


def _coverage_min(args):
    value = getattr(args, "miniprot_rescue_coverage_min", None)
    if value is not None:
        return float(value)
    try:
        return float(os.environ.get("LIFTON_RESCUE_COVERAGE_MIN")
                     or COVERAGE_MIN_DEFAULT)
    except ValueError:
        return COVERAGE_MIN_DEFAULT


def _first_attribute(feature, key):
    try:
        return feature.attributes[key][0]
    except (KeyError, IndexError, TypeError):
        return None


def _protein_length(ref_proteins, ref_trans_id):
    """Reference protein length without its terminal stop, or None."""
    try:
        protein = ref_proteins[ref_trans_id]
        length = len(protein)
        if length and str(protein[length - 1:length]) == "*":
            length -= 1
    except (KeyError, TypeError):
        return None
    return length if length > 0 else None


def miniprot_protein_coverage(mtrans, ref_proteins, ref_trans_id):
    """Fraction of the reference protein a miniprot hit aligns, in [0, 1].

    Reads ``Target=<id> <start> <end>`` (1-based protein coordinates). The
    reference protein's terminal stop, when present, is not counted. Returns
    None when the attribute or the protein is unusable.
    """
    target = _first_attribute(mtrans, "Target")
    fields = target.split() if target else []
    if len(fields) < 3:
        return None
    try:
        qstart, qend = int(fields[1]), int(fields[2])
    except ValueError:
        return None
    length = _protein_length(ref_proteins, ref_trans_id)
    if length is None or qend < qstart:
        return None
    return min(1.0, (qend - qstart + 1) / length)


def _candidate_quality_key(mtrans, coverage):
    """Order candidates best first: miniprot's per-protein ``Rank``, then its
    alignment ``Identity``, then coverage, then position and ID so the order is
    total and run-to-run deterministic.

    ``coverage`` is None for a sub-pass that has no coverage signal, which
    orders on the rest. It has to be tolerated rather than assumed: negating
    None raises, and the caller's ``except Exception`` turned that into a
    silently dropped candidate -- a whole sub-pass doing nothing, reported as
    per-candidate errors nobody reads.
    """
    try:
        rank = int(_first_attribute(mtrans, "Rank"))
    except (TypeError, ValueError):
        rank = sys.maxsize
    try:
        identity = float(_first_attribute(mtrans, "Identity"))
    except (TypeError, ValueError):
        identity = 0.0
    return (rank, -identity, -(coverage or 0.0), mtrans.seqid,
            int(mtrans.start), int(mtrans.end), mtrans.attributes["ID"][0])


def _adaptive_floor_on(args):
    """Whether the divergence-adaptive floor is active. Env
    ``LIFTON_RESCUE_ADAPTIVE_FLOOR`` wins (1/0), else ``args.adaptive_rescue_floor``
    (default True -- PROMOTED to default ON). Opt out with
    ``--no-adaptive-rescue-floor`` / ``LIFTON_RESCUE_ADAPTIVE_FLOOR=0`` (restores
    the fixed ``-miniprot_rescue_min_id`` floor)."""
    env = os.environ.get("LIFTON_RESCUE_ADAPTIVE_FLOOR")
    if env is not None:
        return env not in ("0", "", "false", "False")
    return bool(getattr(args, "adaptive_rescue_floor", True))


def _adaptive_floor(recall, base):
    """Map DNA-lift gene recall -> PI floor in ``[min(ADAPT_FLOOR_MIN, base), base]``:
    ``recall<=ADAPT_R_LOW`` -> the low floor (very-distant, recall is the priority);
    ``recall>=ADAPT_R_HIGH`` -> ``base`` (well-synteny, keep precision); linear in
    between. Never rises above ``base``."""
    lo = min(ADAPT_FLOOR_MIN, base)
    if recall <= ADAPT_R_LOW:
        return lo
    if recall >= ADAPT_R_HIGH:
        return base
    frac = (recall - ADAPT_R_LOW) / (ADAPT_R_HIGH - ADAPT_R_LOW)
    return lo + frac * (base - lo)


def _dna_lift_recall(universe, emitted):
    """Fraction of the miniprot-found ref genes (``universe``) that the DNA lift +
    default Step 8 already placed (``emitted``) -- a cheap, divergence-sensitive
    proxy: ~1 on same-species, low on very-distant. Returns 1.0 (=> base floor) if
    the universe is empty (no miniprot-found ref genes -> the rescue is marginal)."""
    universe = set(universe)
    if not universe:
        return 1.0
    return len(universe & set(emitted)) / len(universe)


def _stage_gene(lifton_gene, cds_id_allocator=None):
    """Serialize to a private buffer and collect stats before committing state."""
    buffer = io.StringIO()
    staged_stats = {'coding': {}, 'non-coding': {}, 'other': {}}
    if cds_id_allocator is None:
        written = lifton_gene.write_entry(buffer, staged_stats)
    else:
        # Claim CDS IDs provisionally: a candidate rejected below must not burn the
        # identifiers, or the transcript that legitimately owns a stable cds- ID gets
        # pushed to a "-1" variant purely because some earlier candidate was discarded.
        with cds_id_allocator.tentative() as scope:
            written = lifton_gene.write_entry(
                buffer,
                staged_stats,
                cds_id_allocator=cds_id_allocator,
            )
            if written is not False:
                scope.commit()
    failures = getattr(lifton_gene, "_serialization_failures", [])
    if written is False:
        return None, None, failures
    return buffer.getvalue(), staged_stats, failures


def _merge_stats(target, staged):
    for feature_type, entries in staged.items():
        bucket = target[feature_type]
        for feature_id, count in entries.items():
            bucket[feature_id] = bucket.get(feature_id, 0) + count


def _record_failure(args, mtrans, error, *, feature_id=None):
    records = getattr(args, "_rescue_failure_records", None)
    if records is None:
        records = []
        args._rescue_failure_records = records
    record = {
        "mRNA": getattr(mtrans, "id", "<unknown>"),
        "message": str(error),
    }
    if feature_id is not None:
        record["feature_id"] = feature_id
    records.append(record)


def rescue_miniprot_only_pass(m_feature_db, ref_db, tree_dict, tgt_fai,
                              ref_proteins, ref_trans, ref_features_dict,
                              m_id_2_ref_id_trans_dict, ref_features_len_dict,
                              ref_trans_exon_num_dict, ref_features_reverse_dict,
                              emitted_ref_gene_ids, fw, fw_score,
                              transcripts_stats_dict, args):
    """Emit a miniprot-only gene model for each reference coding gene the DNA
    lift + default Step 8 missed but miniprot found at a non-overlapping locus,
    protein-identity >= floor, non-redundant.

    Runs strictly after the Step-8 loop closes, so it never perturbs a default
    decision. Returns the number of genes added; inert (returns 0) when
    ``m_feature_db is None``.
    """
    if m_feature_db is None:
        return 0

    # Sub-phase wall clocks, published on args and recorded in the run manifest.
    # The rescue is the largest phase of a distant-species run, and which of its
    # passes the time goes to is not otherwise visible.
    timings = {}
    args._rescue_timings = timings
    started = time.perf_counter()

    def _mark(name, since):
        timings[name] = round(time.perf_counter() - since, 3)
        return time.perf_counter()

    floor = float(getattr(args, "miniprot_rescue_min_id", 0.5))
    added = 0
    isoforms = _rescue_isoforms_on(args)
    publisher = _RescuePublisher(fw, fw_score, transcripts_stats_dict,
                                 defer=isoforms)

    # Deterministic iteration: sort the miniprot mRNAs by (seqid, start, end,
    # ID) so -t8 == -t1 and the cascading suppression (each emit mutates
    # tree_dict, so a later overlapping candidate is suppressed) is fully
    # reproducible run-to-run.
    try:
        mtranscripts = sorted(
            m_feature_db.features_of_type('mRNA'),
            key=lambda m: (m.seqid, m.start, m.end, m.attributes["ID"][0]),
        )
    except Exception as e:
        logger.log_error(f"miniprot-only rescue: failed to enumerate mRNAs: {e}")
        return 0

    # Divergence-adaptive PI floor (default ON): the DNA-lift gene recall is the
    # fraction of the miniprot-found ref genes already emitted by Step 7 + Step 8
    # (~1 same-species, low very-distant); lower the floor toward ADAPT_FLOOR_MIN
    # as recall drops. off subset of on is preserved: this only LOWERS the floor
    # within the already-additive separate pass, so it can only ADD rescues, never
    # remove or change a default-emitted feature.
    if _adaptive_floor_on(args):
        universe = set()
        for _m in mtranscripts:
            try:
                _gid, _tid = lifton_utils.get_ref_ids_miniprot(
                    ref_features_reverse_dict, _m.attributes["ID"][0],
                    m_id_2_ref_id_trans_dict)
            except Exception:
                _gid = None
            if _gid is not None:
                universe.add(_gid)
        recall = _dna_lift_recall(universe, emitted_ref_gene_ids)
        adapted = _adaptive_floor(recall, floor)
        if abs(adapted - floor) > 1e-9:
            sys.stderr.write(
                f"[LiftOn] miniprot-only rescue: adaptive PI floor {adapted:.2f} "
                f"(DNA-lift gene recall {recall:.3f}, base {floor:.2f}).\n")
            sys.stderr.flush()
        floor = adapted
    started = _mark("enumerate_and_floor", started)

    for mtrans in mtranscripts:
        try:
            mtrans_id = mtrans.attributes["ID"][0]

            # (1) resolve canonical ref ids (same namespace as Step 7's
            #     get_ref_ids_liftoff harvest -> dedup is comparable).
            ref_gene_id, ref_trans_id = lifton_utils.get_ref_ids_miniprot(
                ref_features_reverse_dict, mtrans_id, m_id_2_ref_id_trans_dict)
            if ref_gene_id is None or ref_trans_id is None:
                continue

            # (2) DEDUP: skip a ref gene already emitted by Step 7 (DNA lift),
            #     Step 8 (default miniprot), or an earlier rescue in this pass.
            #     This is what makes the pass 0-redundant + off ⊆ on.
            #
            #     A reference gene that belongs at a SECOND target locus is
            #     handled by its own sub-pass after this one, not here: letting
            #     it through mid-walk occupies a locus a later default
            #     candidate would have taken. See _second_locus_subpass.
            if ref_gene_id in emitted_ref_gene_ids:
                continue

            # (3) protein availability (mirror process_miniprot:381)
            if ref_trans_id not in ref_proteins or ref_trans_id not in ref_trans:
                continue

            # (4) overlap suppression vs the FINAL tree_dict: a gene already
            #     sits at this target locus (DNA lift, default Step 8, or a
            #     prior rescue). check_ovps_ratio reads tree_dict per-seqid.
            mtrans_interval = Interval(mtrans.start, mtrans.end, mtrans_id)
            if lifton_utils.check_ovps_ratio(mtrans, mtrans_interval,
                                             args.overlap, tree_dict):
                continue

            # (5) processed-pseudogene filter (mirror process_miniprot:386):
            #     1 CDS in miniprot but >1 exon in the reference.
            if (len(list(m_feature_db.children(mtrans, featuretype='CDS'))) == 1
                    and ref_trans_exon_num_dict.get(ref_trans_id, 0) > 1):
                continue

            # (6) length-ratio sanity band (the WIDE rescue band; the PI floor
            #     at step 8 is the real quality gate). Guards against a
            #     catastrophically mis-scaled miniprot hit.
            ref_len = ref_features_len_dict.get(ref_gene_id)
            if not ref_len:
                continue
            ratio = (mtrans.end - mtrans.start + 1) / ref_len
            if not run_miniprot._miniprot_rescue_band_ok(ratio, args):
                continue

            # (7)-(11) Build, score against the PI floor, ORF-search, stage,
            #          commit, and publish (see _build_and_accept).
            if _build_and_accept(
                    mtrans, ref_gene_id, ref_trans_id, ratio, floor,
                    m_feature_db, ref_db, tree_dict, tgt_fai, ref_proteins,
                    ref_trans, ref_features_dict, emitted_ref_gene_ids,
                    publisher, args):
                added += 1
        except Exception as e:
            logger.log_error(f"miniprot-only rescue error ({mtrans.id}): {e}")
            _record_failure(args, mtrans, e)

    if added:
        sys.stderr.write(
            f"\n[LiftOn] miniprot-only rescue: {added} gene(s) added "
            f"(protein-identity floor {floor:.2f}).\n")
        sys.stderr.flush()

    # Sub-pass B runs only after sub-pass A has finished, and only adds genes
    # at loci still free in the final tree, so everything above is unchanged
    # whether or not it runs.
    started = _mark("subpass_a", started)
    coverage_added = 0
    if _coverage_gate_on(args):
        coverage_added = _coverage_gate_subpass(
            mtranscripts, floor, m_feature_db, ref_db, tree_dict, tgt_fai,
            ref_proteins, ref_trans, ref_features_dict,
            m_id_2_ref_id_trans_dict, ref_features_len_dict,
            ref_trans_exon_num_dict, ref_features_reverse_dict,
            emitted_ref_gene_ids, publisher, args)
    args._rescue_coverage_gate_added = coverage_added
    started = _mark("subpass_b_coverage", started)

    # Isoforms attach only after every gene is placed, so no placement decision
    # can depend on them.
    isoforms_added = 0
    if isoforms:
        isoforms_added = _isoform_pass(
            publisher.accepted, mtranscripts, floor, m_feature_db, ref_db,
            tree_dict, tgt_fai, ref_proteins, ref_trans, ref_features_dict,
            m_id_2_ref_id_trans_dict, ref_features_len_dict,
            ref_features_reverse_dict, args)
        started = _mark("isoform_pass", started)

    # Sub-pass C runs LAST, after the isoform pass, and that ordering is the
    # whole point. Placed before it, a second-locus gene occupies ground a
    # default gene needed to widen into, and `_extension_collides` then refuses
    # the isoform: on human to zebrafish that silently dropped 21 transcripts
    # of mean identity 0.662 -- `gene-CCND3` lost four isoforms scoring
    # 0.584-0.628 to make room for a second-locus model scoring 0.356. The
    # gene-level gate did not see it, because the gene survived and only its
    # transcripts vanished.
    #
    # Running last costs some second-locus genes, since a widened default gene
    # now holds ground they wanted. That is the right trade: it makes the
    # default output -- every gene, every isoform, every span -- provably
    # untouched, which is what "additive" has to mean.
    second_locus_added = 0
    if _second_locus_on(args):
        second_locus_added = _second_locus_subpass(
            mtranscripts, floor, m_feature_db, ref_db, tree_dict, tgt_fai,
            ref_proteins, ref_trans, ref_features_dict,
            m_id_2_ref_id_trans_dict, ref_features_len_dict,
            ref_trans_exon_num_dict, ref_features_reverse_dict,
            emitted_ref_gene_ids, publisher, args)
        started = _mark("subpass_c_second_locus", started)
    args._rescue_second_locus_added = second_locus_added
    publisher.flush()
    _mark("publish", started)
    args._rescue_isoforms_added = isoforms_added
    return added + coverage_added


#: PROMOTED to default-on. The gate was re-scored against the target's own
#: annotation, because source-annotation recall structurally cannot see what
#: this recovers: on human to zebrafish it covers 690 more real GRCz11
#: protein-coding genes, with 0 genes and 0 transcripts lost, 0 default genes
#: moved and 0 regressions, and the off arm byte-identical to the default.
#: The divergence ladder passes its safety gate on 8 of 8 cells, the rate
#: tracks known duplication history (rice to sorghum highest, at 9.06 genes per
#: 1,000 reference proteins), and the same-species control places exactly
#: nothing. See notes/second_locus_rescue_gate.md.
SECOND_LOCUS_DEFAULT = True

#: How many EXTRA loci one reference gene may be given. A duplicated genome
#: wants one; a repeat family must not be allowed to spray copies.
SECOND_LOCUS_MAX_DEFAULT = 1


def _second_locus_on(args):
    """Is the second-locus rescue on? ``LIFTON_RESCUE_SECOND_LOCUS`` wins over
    the resolved flag, as every other rescue switch does."""
    env = os.environ.get("LIFTON_RESCUE_SECOND_LOCUS")
    if env is not None:
        return env.strip().lower() not in ("", "0", "false", "no", "off")
    resolved = getattr(args, "rescue_second_locus", None)
    return SECOND_LOCUS_DEFAULT if resolved is None else bool(resolved)


def _second_locus_max(args):
    env = os.environ.get("LIFTON_RESCUE_SECOND_LOCUS_MAX")
    if env is not None:
        try:
            return max(0, int(env))
        except ValueError:
            pass
    value = getattr(args, "rescue_second_locus_max", None)
    return SECOND_LOCUS_MAX_DEFAULT if value is None else max(0, int(value))


def _second_locus_allowed(ref_gene_id, counts, args):
    """May this already-emitted reference gene be placed at one more locus?"""
    if not _second_locus_on(args):
        return False
    return counts.get(ref_gene_id, 0) < _second_locus_max(args)


def _build_and_accept(mtrans, ref_gene_id, ref_trans_id, ratio, floor,
                      m_feature_db, ref_db, tree_dict, tgt_fai, ref_proteins,
                      ref_trans, ref_features_dict, emitted_ref_gene_ids,
                      publisher, args, extra_attrs=()):
    """Build one miniprot-only gene; if it passes the PI floor and serializes,
    commit its placement and hand it to ``publisher``.

    Returns True when the gene was accepted. Shared by both sub-passes.
    """
    # (7) Build + score against isolated mutable state. Construction must not
    #     consume a copy number or add a suppression interval until every
    #     quality/serialization gate has passed.
    state_journal = DeferredStateJournal(
        ref_features_dict, buffer_score=True,
    )
    lifton_gene, lifton_trans, transcript_id, lifton_status = \
        run_miniprot.lifton_miniprot_with_ref_protein(
            mtrans, m_feature_db, ref_db.db_connection, ref_gene_id,
            ref_trans_id, tgt_fai, ref_proteins, ref_trans, tree_dict,
            ref_features_dict, args, state_journal=state_journal)

    # (8) PI floor -- the quality gate (mirror the in-loop pre-ORF check at
    #     process_miniprot:401, on the miniprot alignment identity set by
    #     lifton_miniprot_with_ref_protein).
    if lifton_status.lifton_aa < floor:
        return False

    # (9) tag attrs BEFORE the ORF/status tail (parity with the in-loop rescue
    #     ordering at process_miniprot:403-404).
    attributes = lifton_gene.transcripts[transcript_id].entry.attributes
    attributes["miniprot_annotation_ratio"] = [f"{ratio:.3f}"]
    attributes["lifton_rescue"] = ["miniprot_only"]
    for key, value in extra_attrs:
        attributes[key] = [value]

    # (10) ORF search + status tail (mirror process_miniprot:409-413).
    lifton_gene.orf_search_protein(lifton_trans.entry.id, ref_trans_id,
                                   tgt_fai, ref_proteins, ref_trans,
                                   lifton_status)
    # miniprot's CDS ends at the last aligned codon, so the model lacks the stop
    # the reference protein has. Complete it AFTER the ORF search, which
    # therefore sees exactly the sequence it would have seen.
    orf_completion.complete_and_rescore(
        lifton_trans, mtrans, tgt_fai, ref_proteins, ref_trans_id,
        lifton_status, args)
    lifton_utils.print_lifton_status(transcript_id, mtrans,
                                     lifton_status, DEBUG=args.debug)
    lifton_gene.add_lifton_gene_status_attrs("miniprot")
    lifton_gene.add_lifton_trans_status_attrs(transcript_id, lifton_status)
    # (11) Stage the complete hierarchy, commit the copy/interval journal, then
    #      publish. A rejected or unserializable model is fully inert, making
    #      fixed-floor output a strict subset of adaptive-floor output.
    block, staged_stats, serialization_failures = _stage_gene(
        lifton_gene,
        cds_id_allocator=getattr(args, "_cds_id_allocator", None),
    )
    for feature_id, detail in serialization_failures:
        _record_failure(args, mtrans, detail, feature_id=feature_id)
    if block is None:
        return False
    lifton_utils.write_lifton_status(
        state_journal.score_handle, transcript_id, mtrans, lifton_status,
    )
    delta = state_journal.finish()
    commit_locus_delta(delta, ref_features_dict, tree_dict)
    emitted_ref_gene_ids.add(ref_gene_id)
    publisher.publish(_Accepted(lifton_gene, mtrans, ref_gene_id, ref_trans_id,
                                block, staged_stats, delta.score_text))
    return True


def _coverage_gate_subpass(mtranscripts, floor, m_feature_db, ref_db,
                           tree_dict, tgt_fai, ref_proteins, ref_trans,
                           ref_features_dict, m_id_2_ref_id_trans_dict,
                           ref_features_len_dict, ref_trans_exon_num_dict,
                           ref_features_reverse_dict, emitted_ref_gene_ids,
                           publisher, args):
    """Sub-pass B: rescue candidates the span band rejected, by coverage.

    Every other gate is the same as sub-pass A (dedup, protein availability,
    overlap with the now-final tree, the processed-pseudogene filter, the PI
    floor). Candidates are tried best first rather than in coordinate order,
    so when two hits compete for one gene or one locus the better one wins.
    """
    coverage_min = _coverage_min(args)
    candidates = []
    for mtrans in mtranscripts:
        try:
            mtrans_id = mtrans.attributes["ID"][0]
            ref_gene_id, ref_trans_id = lifton_utils.get_ref_ids_miniprot(
                ref_features_reverse_dict, mtrans_id, m_id_2_ref_id_trans_dict)
            if ref_gene_id is None or ref_trans_id is None:
                continue
            if ref_gene_id in emitted_ref_gene_ids:
                continue
            if ref_trans_id not in ref_proteins or ref_trans_id not in ref_trans:
                continue
            ref_len = ref_features_len_dict.get(ref_gene_id)
            if not ref_len:
                continue
            ratio = (mtrans.end - mtrans.start + 1) / ref_len
            # Inside the band, sub-pass A already made the decision.
            if run_miniprot._miniprot_rescue_band_ok(ratio, args):
                continue
            coverage = miniprot_protein_coverage(mtrans, ref_proteins,
                                                 ref_trans_id)
            if coverage is None or coverage < coverage_min:
                continue
            cds_children = list(m_feature_db.children(mtrans, featuretype='CDS'))
            cds_length = sum(int(c.end) - int(c.start) + 1 for c in cds_children)
            protein_length = _protein_length(ref_proteins, ref_trans_id)
            if cds_length > CDS_LENGTH_MAX_RATIO * 3 * protein_length:
                continue
            candidates.append((_candidate_quality_key(mtrans, coverage), mtrans,
                               ref_gene_id, ref_trans_id, ratio, coverage,
                               len(cds_children)))
        except Exception as e:
            logger.log_error(
                f"miniprot-only rescue (coverage gate) error ({mtrans.id}): {e}")
            _record_failure(args, mtrans, e)
    candidates.sort(key=lambda candidate: candidate[0])

    added = 0
    for (_, mtrans, ref_gene_id, ref_trans_id, ratio, coverage,
         n_cds) in candidates:
        try:
            if ref_gene_id in emitted_ref_gene_ids:
                continue
            mtrans_interval = Interval(mtrans.start, mtrans.end,
                                       mtrans.attributes["ID"][0])
            if lifton_utils.check_ovps_ratio(mtrans, mtrans_interval,
                                             args.overlap, tree_dict):
                continue
            # Processed-pseudogene filter, as in sub-pass A.
            if n_cds == 1 and ref_trans_exon_num_dict.get(ref_trans_id, 0) > 1:
                continue
            if _build_and_accept(
                    mtrans, ref_gene_id, ref_trans_id, ratio, floor,
                    m_feature_db, ref_db, tree_dict, tgt_fai, ref_proteins,
                    ref_trans, ref_features_dict, emitted_ref_gene_ids,
                    publisher, args,
                    extra_attrs=(("rescue_gate", "protein_coverage"),
                                 ("miniprot_protein_coverage",
                                  f"{coverage:.3f}"))):
                added += 1
        except Exception as e:
            logger.log_error(
                f"miniprot-only rescue (coverage gate) error ({mtrans.id}): {e}")
            _record_failure(args, mtrans, e)

    if added:
        sys.stderr.write(
            f"[LiftOn] miniprot-only rescue, protein-coverage gate: {added} "
            f"gene(s) added (coverage >= {coverage_min:.2f}, protein-identity "
            f"floor {floor:.2f}).\n")
        sys.stderr.flush()
    return added


def _second_locus_subpass(mtranscripts, floor, m_feature_db, ref_db,
                          tree_dict, tgt_fai, ref_proteins, ref_trans,
                          ref_features_dict, m_id_2_ref_id_trans_dict,
                          ref_features_len_dict, ref_trans_exon_num_dict,
                          ref_features_reverse_dict, emitted_ref_gene_ids,
                          publisher, args):
    """Sub-pass C: place a reference gene at a SECOND target locus.

    A whole-genome duplication gives the target two genes where the reference
    has one, and the rescue's reference-gene dedup refuses the second. Measured
    against zebrafish's own GRCz11 annotation, that refusal hides 2,051 real
    target genes on human to zebrafish
    (notes/coortholog_recall_measurement_2026-09.md).

    This runs AFTER sub-passes A and B, for the reason sub-pass B runs after A.
    A first cut relaxed the dedup inside sub-pass A instead, reasoning that the
    overlap gate would keep it additive. It does not: an acceptance commits its
    interval into the shared tree, so an extra gene placed early in the walk
    suppresses a later default candidate that wanted the same free locus. The
    A/B measured the cost -- 119 genes lost, 165 models overlapping one already
    emitted, 33 transcripts regressed against 7 improved -- which is the
    Iteration-22 swap in a new place. Deferring the whole sub-pass restores the
    property the design claimed: every default decision is already committed,
    so nothing here can change one.

    ``emitted_ref_gene_ids`` is therefore inverted relative to its use above:
    here a gene must ALREADY be in it. The per-gene cap bounds how many extra
    loci one reference gene may take, so a repeat family cannot spray copies.
    """
    maximum = _second_locus_max(args)
    if maximum <= 0:
        return 0
    candidates = []
    for mtrans in mtranscripts:
        try:
            mtrans_id = mtrans.attributes["ID"][0]
            ref_gene_id, ref_trans_id = lifton_utils.get_ref_ids_miniprot(
                ref_features_reverse_dict, mtrans_id, m_id_2_ref_id_trans_dict)
            if ref_gene_id is None or ref_trans_id is None:
                continue
            # The inversion: only a gene already placed is a candidate here.
            if ref_gene_id not in emitted_ref_gene_ids:
                continue
            if ref_trans_id not in ref_proteins or ref_trans_id not in ref_trans:
                continue
            ref_len = ref_features_len_dict.get(ref_gene_id)
            if not ref_len:
                continue
            ratio = (mtrans.end - mtrans.start + 1) / ref_len
            if not run_miniprot._miniprot_rescue_band_ok(ratio, args):
                continue
            cds_children = list(m_feature_db.children(mtrans, featuretype='CDS'))
            candidates.append((_candidate_quality_key(mtrans, None), mtrans,
                               ref_gene_id, ref_trans_id, ratio,
                               len(cds_children)))
        except Exception as e:
            logger.log_error(
                f"miniprot-only rescue (second locus) error ({mtrans.id}): {e}")
            _record_failure(args, mtrans, e)
    # Best first, so when two hits compete for one free locus the better wins.
    candidates.sort(key=lambda candidate: candidate[0])

    added = 0
    counts = {}
    for _, mtrans, ref_gene_id, ref_trans_id, ratio, n_cds in candidates:
        try:
            if counts.get(ref_gene_id, 0) >= maximum:
                continue
            mtrans_interval = Interval(mtrans.start, mtrans.end,
                                       mtrans.attributes["ID"][0])
            if lifton_utils.check_ovps_ratio(mtrans, mtrans_interval,
                                             args.overlap, tree_dict):
                continue
            if n_cds == 1 and ref_trans_exon_num_dict.get(ref_trans_id, 0) > 1:
                continue
            if _build_and_accept(
                    mtrans, ref_gene_id, ref_trans_id, ratio, floor,
                    m_feature_db, ref_db, tree_dict, tgt_fai, ref_proteins,
                    ref_trans, ref_features_dict, emitted_ref_gene_ids,
                    publisher, args,
                    extra_attrs=(("lifton_rescue_second_locus", "true"),)):
                added += 1
                counts[ref_gene_id] = counts.get(ref_gene_id, 0) + 1
        except Exception as e:
            logger.log_error(
                f"miniprot-only rescue (second locus) error ({mtrans.id}): {e}")
            _record_failure(args, mtrans, e)

    if added:
        sys.stderr.write(
            f"[LiftOn] miniprot-only rescue, second locus: {added} gene(s) "
            f"placed at a locus no emitted model reaches (at most {maximum} "
            f"per reference gene).\n")
        sys.stderr.flush()
    return added


def _extension_collides(tree_dict, seqid, gene_id, old_start, old_end,
                        new_start, new_end):
    """True if widening [old_start, old_end] to [new_start, new_end] would
    reach into any other gene's interval in the (final) suppression tree."""
    tree = tree_dict.get(seqid)
    if tree is None:
        return False
    for low, high in ((new_start, old_start - 1), (old_end + 1, new_end)):
        if low > high:
            continue
        for interval in tree.overlap(low, high + 1):
            if interval.data != gene_id:
                return True
    return False


class _GeneView:
    """The immutable fields of a placed gene that an isoform build needs, so
    worker threads never touch the live gene object."""

    __slots__ = ("ref_gene_id", "gene_id", "copy_num", "is_non_coding",
                 "stop_completion")

    def __init__(self, lifton_gene, stop_completion=True):
        self.ref_gene_id = lifton_gene.ref_gene_id
        self.gene_id = lifton_gene.entry.id
        self.copy_num = lifton_gene.copy_num
        self.is_non_coding = lifton_gene.is_non_coding
        # Resolved on the parent: a forked worker has no ``args``.
        self.stop_completion = stop_completion


def _isoform_candidates(accepted, hits, ref_proteins, ref_trans):
    """Best hit per other transcript of the gene, co-located with the placed
    hit, in the order the serial pass tries them."""
    primary = accepted.mtrans
    best = {}
    for mtrans, ref_trans_id in hits:
        if ref_trans_id == accepted.ref_trans_id:
            continue
        if mtrans.seqid != primary.seqid or mtrans.strand != primary.strand:
            continue
        if mtrans.end < primary.start or mtrans.start > primary.end:
            continue
        if ref_trans_id not in ref_proteins or ref_trans_id not in ref_trans:
            continue
        key = _candidate_quality_key(mtrans, 0.0)
        if ref_trans_id not in best or key < best[ref_trans_id][0]:
            best[ref_trans_id] = (key, mtrans)
    return [(ref_trans_id, mtrans) for ref_trans_id, (_, mtrans)
            in sorted(best.items(), key=lambda item: item[1][0])]


def _score_isoform(view, mtrans, m_entry, cds_children, ref_trans_attrs,
                   ref_trans_id, floor, ref_len, tgt_fai, ref_proteins, ref_trans):
    """Build and score one detached isoform transcript.

    Mirrors ``run_miniprot.lifton_miniprot_with_ref_protein`` followed by the
    rescue's floor, tags and ORF search. Every database row arrives
    prefetched, and nothing shared is written, so it is safe on a worker
    thread. Returns ``(transcript, status, passed_floor)``.
    """
    transcript = lifton_class.Lifton_TRANS(
        ref_trans_id, view.ref_gene_id, view.gene_id, view.copy_num,
        coreutils.clone_feature(mtrans), ref_trans_attrs)
    transcript.entry.seqid = mtrans.seqid
    transcript.entry.start = mtrans.start
    transcript.entry.end = mtrans.end
    for cds in cds_children:
        transcript.add_exon(cds)
        transcript.add_cds(coreutils.clone_feature(cds))
    status = lifton_class.Lifton_Status()
    alignment = align.lifton_parasail_align(transcript, m_entry, tgt_fai,
                                            ref_proteins, ref_trans_id)
    status.annotation = "miniprot"
    status.lifton_aa = alignment.identity
    if status.lifton_aa < floor:
        return transcript, status, False
    attributes = transcript.entry.attributes
    if ref_len:
        attributes["miniprot_annotation_ratio"] = [
            f"{(mtrans.end - mtrans.start + 1) / ref_len:.3f}"]
    attributes["lifton_rescue"] = ["miniprot_only"]
    attributes["rescue_isoform"] = ["true"]
    ref_protein_seq = (str(ref_proteins[ref_trans_id])
                       if ref_trans_id in ref_proteins.keys() else None)
    ref_trans_seq = (str(ref_trans[ref_trans_id])
                     if ref_trans_id in ref_trans.keys() else None)
    transcript.orf_search_protein(tgt_fai, ref_protein_seq, ref_trans_seq,
                                  status, is_non_coding=view.is_non_coding,
                                  eval_only=False)
    orf_completion.complete_and_rescore(
        transcript, m_entry, tgt_fai, ref_proteins, ref_trans_id, status,
        enabled_override=view.stop_completion)
    transcript.add_lifton_trans_status_attrs(status)
    return transcript, status, True


def _attach_scored_isoforms(accepted, scored, tree_dict, args):
    """Attach one gene's scored isoforms in candidate order; return how many.

    The decision per candidate is the serial pass's: skip it if widening the
    gene to cover it would reach another gene's interval in the (final)
    suppression tree, otherwise keep it if it cleared the identity floor. A
    widened gene's interval joins the tree so two genes cannot widen into
    each other. The gene is re-staged once; if that fails, the placement-time
    model is kept unchanged.
    """
    lifton_gene = accepted.lifton_gene
    gene_id = lifton_gene.entry.id
    seqid = lifton_gene.entry.seqid
    original_span = (lifton_gene.entry.start, lifton_gene.entry.end)
    attached, widened = [], []
    for mtrans, outcome in scored:
        start, end = int(lifton_gene.entry.start), int(lifton_gene.entry.end)
        new_start = min(start, int(mtrans.start))
        new_end = max(end, int(mtrans.end))
        if _extension_collides(tree_dict, seqid, gene_id, start, end,
                               new_start, new_end):
            continue
        if isinstance(outcome, BaseException):
            # The serial pass built this candidate here and would have
            # stopped the gene with this error.
            raise outcome
        transcript, status, passed = outcome
        if not passed:
            continue
        lifton_utils.print_lifton_status(transcript.entry.id, mtrans, status,
                                         DEBUG=args.debug)
        lifton_gene.transcripts[transcript.entry.id] = transcript
        attached.append((transcript.entry.id, mtrans, status))
        if (new_start, new_end) != (start, end):
            lifton_gene.entry.start, lifton_gene.entry.end = new_start, new_end
            interval = _make_interval(new_start, new_end, gene_id)
            tree_dict.setdefault(seqid, _new_tree()).add(interval)
            widened.append(interval)
    if not attached:
        lifton_gene.entry.start, lifton_gene.entry.end = original_span
        return 0

    block, staged_stats, serialization_failures = _stage_gene(
        lifton_gene, cds_id_allocator=getattr(args, "_cds_id_allocator", None))
    if block is None:
        for feature_id, detail in serialization_failures:
            _record_failure(args, accepted.mtrans, detail, feature_id=feature_id)
        for transcript_id, _, _ in attached:
            del lifton_gene.transcripts[transcript_id]
        for interval in widened:
            tree_dict[seqid].discard(interval)
        lifton_gene.entry.start, lifton_gene.entry.end = original_span
        return 0
    score = io.StringIO()
    for transcript_id, mtrans, status in attached:
        lifton_utils.write_lifton_status(score, transcript_id, mtrans, status)
    accepted.block = block
    accepted.staged_stats = staged_stats
    accepted.score_text = (accepted.score_text or "") + score.getvalue()
    return len(attached)


def _new_tree():
    from intervaltree import IntervalTree
    return IntervalTree()


def _isoform_pass(accepted, mtranscripts, floor, m_feature_db, ref_db,
                  tree_dict, tgt_fai, ref_proteins, ref_trans,
                  ref_features_dict, m_id_2_ref_id_trans_dict,
                  ref_features_len_dict, ref_features_reverse_dict, args):
    """Attach co-located isoforms to every accepted gene, in acceptance order.

    Three phases. The main thread gathers each gene's candidates and prefetches
    their database rows. The detached transcripts are then built and scored,
    by forked worker processes when there are threads and work to share (see
    ``_score_isoform_jobs``). The main thread finally attaches them gene by
    gene with the serial decision rule, so the output does not depend on
    ``--threads``.
    """
    if not accepted:
        return 0
    hits_by_gene = defaultdict(list)
    for mtrans in mtranscripts:
        try:
            ref_gene_id, ref_trans_id = lifton_utils.get_ref_ids_miniprot(
                ref_features_reverse_dict, mtrans.attributes["ID"][0],
                m_id_2_ref_id_trans_dict)
        except Exception:
            continue
        if ref_gene_id is not None and ref_trans_id is not None:
            hits_by_gene[ref_gene_id].append((mtrans, ref_trans_id))

    max_inflight = _rescue_max_inflight(args)
    timings = getattr(args, "_rescue_timings", None)
    elapsed = {"prefetch": 0.0, "score": 0.0, "attach": 0.0}
    counters = {"jobs": 0, "high_water": 0}
    added = 0

    def _drain(plans, jobs):
        """Score and attach one bounded batch, then let it go.

        Scoring reads nothing that attachment writes -- ``tree_dict`` is
        touched only by ``_attach_scored_isoforms`` -- and records are attached
        in acceptance order either way, so draining in batches gives the same
        output as draining once at the end. That is what makes the bound free.
        """
        nonlocal added
        if not plans:
            return
        counters["jobs"] += len(jobs)
        counters["high_water"] = max(counters["high_water"], len(jobs))
        started = time.perf_counter()
        results = _score_isoform_jobs(jobs, tgt_fai, ref_proteins, ref_trans,
                                      args)
        elapsed["score"] += time.perf_counter() - started
        started = time.perf_counter()
        for record, slots in plans:
            if not slots:
                continue
            scored = [(mtrans,
                       value if isinstance(value, BaseException)
                       else results[value])
                      for mtrans, value in slots]
            try:
                added += _attach_scored_isoforms(record, scored, tree_dict, args)
            except Exception as e:
                logger.log_error(
                    f"miniprot-only rescue (isoforms) error "
                    f"({record.mtrans.id}): {e}")
                _record_failure(args, record.mtrans, e)
        elapsed["attach"] += time.perf_counter() - started

    plans, jobs = [], []
    phase_started = time.perf_counter()
    stop_completion = orf_completion.enabled(args)
    for record in accepted:
        candidates = _isoform_candidates(
            record, hits_by_gene.get(record.ref_gene_id, ()), ref_proteins,
            ref_trans)
        if not candidates:
            plans.append((record, []))
            continue
        view = _GeneView(record.lifton_gene, stop_completion)
        ref_len = ref_features_len_dict.get(record.ref_gene_id) or 0
        slots = []
        for ref_trans_id, mtrans in candidates:
            try:
                # ``mtrans`` came from the same features_of_type('mRNA') sweep,
                # so it IS the row a lookup by its own ID would return; the
                # scorer only ever reads the entry's strand and stores it.
                prefetched = (
                    view, mtrans, mtrans,
                    list(m_feature_db.children(mtrans, featuretype='CDS')),
                    ref_db.db_connection[ref_trans_id].attributes,
                    ref_trans_id, floor, ref_len)
            except Exception as error:   # surfaces where the serial pass would
                slots.append((mtrans, error))
                continue
            slots.append((mtrans, len(jobs)))
            jobs.append(prefetched)
        plans.append((record, slots))
        # Hold at most `max_inflight` prefetched jobs at once. Unbounded, this
        # list held every isoform job in the genome: 55,852 of them on human to
        # zebrafish, each keeping an mRNA row, its CDS rows and an attribute
        # map alive, and the pass accounted for 3.72 GiB of a 6.26 GiB peak.
        if max_inflight and len(jobs) >= max_inflight:
            elapsed["prefetch"] += time.perf_counter() - phase_started
            _drain(plans, jobs)
            plans, jobs = [], []
            phase_started = time.perf_counter()

    elapsed["prefetch"] += time.perf_counter() - phase_started
    _drain(plans, jobs)
    plans, jobs = [], []
    if timings is not None:
        timings["isoform_prefetch"] = round(elapsed["prefetch"], 3)
        timings["isoform_score"] = round(elapsed["score"], 3)
        timings["isoform_attach"] = round(elapsed["attach"], 3)
        timings["isoform_jobs"] = counters["jobs"]
        timings["isoform_jobs_high_water"] = counters["high_water"]
    if added:
        sys.stderr.write(
            f"[LiftOn] miniprot-only rescue: {added} additional isoform(s) "
            f"attached to rescued genes.\n")
        sys.stderr.flush()
    return added


# Scoring an isoform is mostly Python (ORF search, variant calls, the windowed
# aligner's anchoring), so threads do not overlap it; forked processes do. The
# jobs are shared copy-on-write through this module global, set just before
# the pool forks and cleared after; only index ranges go out and scored
# transcripts come back.
_ISOFORM_SHARED = {}
ISOFORM_PROCESS_MIN_JOBS = 32


def _reset_inherited_signal_handlers():
    """A forked worker must not run the parent's SIGTERM/SIGHUP handlers
    (``Pool.terminate`` signals its workers); restore the defaults."""
    import signal
    for signum in (signal.SIGTERM, getattr(signal, "SIGHUP", None)):
        if signum is not None:
            signal.signal(signum, signal.SIG_DFL)


def _isoform_worker_init():
    from pyfaidx import Fasta
    _reset_inherited_signal_handlers()
    # A forked child shares the parent's file offsets, so reopen every FASTA
    # (LiftOn opens them all as plain ``Fasta(path)``).
    _ISOFORM_SHARED["fastas"] = tuple(
        Fasta(path) for path in _ISOFORM_SHARED["paths"])


def _score_isoform_range(bounds):
    start, stop = bounds
    tgt_fai, ref_proteins, ref_trans = _ISOFORM_SHARED["fastas"]
    out = []
    for job in _ISOFORM_SHARED["jobs"][start:stop]:
        try:
            out.append(_score_isoform(*job, tgt_fai, ref_proteins, ref_trans))
        except Exception as error:
            # Any exception must survive pickling back to the parent.
            out.append(RuntimeError(f"{type(error).__name__}: {error}"))
    return out


#: How many prefetched isoform jobs may be held at once. Output is identical at
#: any value -- scoring reads nothing attachment writes -- so this trades peak
#: memory against pool restarts. 0 restores the unbounded single batch.
RESCUE_MAX_INFLIGHT_DEFAULT = 8192


def _rescue_max_inflight(args):
    env = os.environ.get("LIFTON_RESCUE_MAX_INFLIGHT")
    if env is not None:
        try:
            return max(0, int(env))
        except ValueError:
            pass
    value = getattr(args, "rescue_max_inflight", None)
    return RESCUE_MAX_INFLIGHT_DEFAULT if value is None else max(0, int(value))


def _isoform_workers(args, n_jobs):
    env = os.environ.get("LIFTON_RESCUE_ISOFORM_WORKERS")
    if env is not None:
        return max(0, int(env))
    if n_jobs < ISOFORM_PROCESS_MIN_JOBS:
        return 0
    return max(1, int(getattr(args, "threads", 1) or 1))


def _score_isoform_jobs_serial(jobs, tgt_fai, ref_proteins, ref_trans):
    """The in-process scoring loop the forked pool parallelises."""
    results = []
    for job in jobs:
        try:
            results.append(_score_isoform(*job, tgt_fai, ref_proteins,
                                          ref_trans))
        except Exception as error:
            results.append(error)
    return results


def _score_isoform_jobs(jobs, tgt_fai, ref_proteins, ref_trans, args):
    """Score every prefetched isoform job, in order; an exception is returned
    in place of its result. Output is identical whichever way it runs."""
    workers = _isoform_workers(args, len(jobs))
    paths = tuple(getattr(fasta, "filename", None)
                  for fasta in (tgt_fai, ref_proteins, ref_trans))
    if workers <= 1 or None in paths:
        return _score_isoform_jobs_serial(jobs, tgt_fai, ref_proteins,
                                          ref_trans)
    import multiprocessing
    chunks = max(1, min(len(jobs), workers * 4))
    step = -(-len(jobs) // chunks)
    ranges = [(start, min(start + step, len(jobs)))
              for start in range(0, len(jobs), step)]
    _ISOFORM_SHARED.update(jobs=jobs, paths=paths)
    try:
        try:
            with multiprocessing.get_context("fork").Pool(
                    workers, initializer=_isoform_worker_init) as pool:
                parts = pool.map(_score_isoform_range, ranges, chunksize=1)
        except OSError as error:
            # Same hazard as the parallel lift: under strict overcommit the
            # kernel charges each fork the parent's whole address space, so a
            # large parent can fail to start workers with hundreds of GB free.
            # Scoring is identical in-process, so degrade instead of aborting.
            logger.log_warning(
                f"Isoform rescue could not start {workers} worker(s) ({error}); "
                f"scoring {len(jobs)} isoform job(s) in-process. Set "
                f"LIFTON_RESCUE_ISOFORM_WORKERS to a smaller value to keep the "
                f"parallel path."
            )
            return _score_isoform_jobs_serial(jobs, tgt_fai, ref_proteins,
                                              ref_trans)
    finally:
        _ISOFORM_SHARED.clear()
    return [result for part in parts for result in part]
