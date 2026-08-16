"""Custom evaluator (PRIMARY metric) reusing LiftOn's own parasail kernel.

For each tool output GFF3 + target genome, resolve every lifted transcript
container (``mRNA`` or ``transcript``) to its reference transcript, extract the
lifted transcript DNA + translated protein from the
TARGET genome with LiftOn's own ``extract_sequence`` helpers, and align them
against the reference transcript/protein with LiftOn's ``align`` +
``get_id_fraction`` kernel. Emits a per-transcript TSV and a summary JSON with
completeness, protein-identity, and DNA-identity statistics.

DNA basis: miniprot has no exon/UTR (CDS only), so for miniprot both reference
and lifted DNA are CDS-spliced (``dna_basis="cds"``); for Liftoff/LiftOn both are
exon-spliced full transcripts (``dna_basis="transcript"``). Protein identity is
always CDS-translated and directly comparable across all three tools.
"""
from __future__ import annotations

import json
import os
import re
import statistics
import subprocess
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import gffutils
import pyfaidx

from lifton import align, annotation as _lifton_annotation
from lifton.exceptions import LiftOnAlignmentError
from lifton.extract_sequence import get_dna_sequence, get_protein_sequence

from . import id_mapping

_ERR_SUMMARY_RE = re.compile(r"Errors\s*:\s*(\d+)")


def count_gff3_validator_errors(gff_path, python_exe, env) -> int:
    """Structured error count from ``python -m lifton.gff3_validator``.

    Parses the single authoritative ``Errors   : N`` summary line the validator
    prints (gff3_validator.py:336). This replaces the crude ``grep -c "error"``
    line count the A/B harnesses used, which conflated the summary line, the
    ``[ERROR] Check: …`` detail lines, and any incidental "error" substrings —
    and consequently *undercounted* (candidate-3 drosophila read as 1->2 when the
    true structured delta was 0->4 strand_consistency errors). Falls back to
    counting ``[ERROR]`` detail lines if the summary line is ever absent.
    """
    r = subprocess.run(
        [python_exe, "-m", "lifton.gff3_validator", str(gff_path)],
        env=env, capture_output=True, text=True)
    out = (r.stdout or "") + (r.stderr or "")
    m = _ERR_SUMMARY_RE.findall(out)
    if m:
        return int(m[-1])
    return sum(1 for ln in out.splitlines() if "[ERROR]" in ln)


def _build_db(gff_path: str):
    """Build a gffutils FeatureDB over a GFF, reusing LiftOn's Annotation loader
    which has a 3-strategy fallback for messy real-world GFFs (e.g. Liftoff
    output with thousands of duplicate ``ID=1.0`` features and multi-segment CDS
    sharing IDs — which plain ``create_unique`` chokes on)."""
    # Always rebuild from the CURRENT GFF: unlink any stale gffutils sidecar
    # first. Without this, re-scoring an OVERWRITTEN output GFF silently reused
    # the previous run's `<gff>_db` (the LiftOn Annotation loader's force= did
    # not rebuild the on-disk file), so the eval scored the OLD output — e.g. a
    # corrected full devel.gff3 (7856 mRNA) was scored against a leftover
    # 2730-mRNA DB, fabricating a spurious under-recovery. Mirrors
    # tool_runners._clean_input_dbs for input GFFs.
    for sidecar in (str(gff_path) + "_db", str(gff_path) + ".eval_db"):
        try:
            os.unlink(sidecar)
        except FileNotFoundError:
            pass
    try:
        return _lifton_annotation.Annotation(
            str(gff_path), infer_genes=False, infer_transcripts=False,
            merge_strategy="create_unique", id_spec=None, force=True,
            verbose=False, auto_convert_gtf=False,
        ).db_connection
    except Exception:
        # last-ditch: gffutils with the most tolerant merge strategy
        return gffutils.create_db(
            str(gff_path), dbfn=str(gff_path) + ".eval_db", force=True,
            keep_order=True, merge_strategy="create_unique",
            id_spec=None, sort_attribute_values=False,
            disable_infer_genes=True, disable_infer_transcripts=True,
        )


def _children(db, mrna, ftype):
    try:
        return list(db.children(mrna, featuretype=ftype, order_by="start"))
    except Exception:
        return []


def _transcript_ftype(db):
    """Choose the preferred transcript-container traversal order.

    RefSeq/NCBI annotate coding transcripts as ``mRNA``; Ensembl/gffread GTF use
    ``transcript``. Both types are scored by :func:`_transcript_features`; putting
    the dominant type first preserves the historical RefSeq traversal order and
    gives deterministic precedence when malformed inputs reuse one logical ID
    across both types."""
    try:
        n_m = db.count_features_of_type("mRNA")
        n_t = db.count_features_of_type("transcript")
    except Exception:
        n_m = sum(1 for _ in db.features_of_type("mRNA"))
        n_t = sum(1 for _ in db.features_of_type("transcript"))
    return "transcript" if n_t > n_m else "mRNA"


def _transcript_ftypes(db):
    """Return present coding-transcript types with the dominant type first."""
    preferred = _transcript_ftype(db)
    other = "transcript" if preferred == "mRNA" else "mRNA"
    counts = {}
    for feature_type in (preferred, other):
        try:
            counts[feature_type] = db.count_features_of_type(feature_type)
        except Exception:
            counts[feature_type] = sum(
                1 for _ in db.features_of_type(feature_type)
            )
    return tuple(
        feature_type for feature_type in (preferred, other)
        if counts[feature_type] > 0
    )


def _transcript_features(db):
    """Yield mRNA/transcript rows once, treating their types as equivalent."""
    seen_ids = set()
    for feature_type in _transcript_ftypes(db):
        for feature in db.features_of_type(feature_type):
            if feature.id in seen_ids:
                continue
            seen_ids.add(feature.id)
            yield feature


# sub-feature types excluded from the "overall feature completeness" aggregate
# (we count genes / transcripts / ncRNAs / pseudogenes as features, not their
# exon/CDS/UTR parts).
_SUBFEATURE_TYPES = frozenset({
    "exon", "CDS", "start_codon", "stop_codon",
    "five_prime_UTR", "three_prime_UTR", "intron",
})


def feature_index(db):
    """One pass over a FeatureDB -> (ids_by_type, all_ids, census).

    ids_by_type: {featuretype: set(feature_id)}; all_ids: every feature id (used
    as the copy-suffix base check); census: {featuretype: count}.
    """
    ids_by_type = defaultdict(set)
    all_ids = set()
    census = Counter()
    for f in db.all_features():
        ids_by_type[f.featuretype].add(f.id)
        all_ids.add(f.id)
        census[f.featuretype] += 1
    return ids_by_type, all_ids, census


def completeness_by_type(tool_db, ref_ids_by_type, ref_all_ids):
    """Per reference feature type, how many reference features survived into the
    tool output (matched by id, copy-suffix-aware via id_mapping.strip_copy_suffix).

    Also records the tool's own per-type feature census and an ``_overall_``
    aggregate over non-sub-feature types. NOTE: miniprot emits its own MP* ids
    (resolved to the reference only via the Target attribute, captured separately
    by completeness_coding), so its id-match recovery here is meaningful only for
    the mRNA/CDS it produces — the per-type recovery view is primarily a
    full-annotation (Liftoff/LiftOn) metric.
    """
    tool_base_counts = Counter()
    tool_census = Counter()
    for f in tool_db.all_features():
        tool_census[f.featuretype] += 1
        base, _ = id_mapping.strip_copy_suffix(str(f.id), ref_all_ids)
        tool_base_counts[base] += 1

    out = {}
    agg_ref = agg_rec = 0
    for ftype, rids in ref_ids_by_type.items():
        n_ref = len(rids)
        if not n_ref:
            continue
        n_rec = sum(1 for rid in rids if tool_base_counts.get(rid, 0) >= 1)
        n_extra = sum(max(0, tool_base_counts.get(rid, 0) - 1) for rid in rids)
        tool_feature_count = tool_census.get(ftype, 0)
        if ftype in {"mRNA", "transcript"}:
            tool_feature_count = sum(
                tool_census.get(transcript_type, 0)
                for transcript_type in ("mRNA", "transcript")
            )
        out[ftype] = {
            "n_reference": n_ref,
            "n_recovered": n_rec,
            "pct_recovered": round(n_rec / n_ref, 5),
            "n_extra_copies": n_extra,
            "n_tool_features": tool_feature_count,
        }
        if ftype not in _SUBFEATURE_TYPES:
            agg_ref += n_ref
            agg_rec += n_rec
    out["_overall_"] = {
        "n_reference": agg_ref,
        "n_recovered": agg_rec,
        "pct_recovered": round(agg_rec / agg_ref, 5) if agg_ref else None,
        "n_extra_copies": 0,
        "n_tool_features": sum(v for k, v in tool_census.items()
                               if k not in _SUBFEATURE_TYPES),
    }
    return out


def _safe_prot_identity(lifted_prot, ref_prot):
    if not lifted_prot or not ref_prot:
        return None
    try:
        return align.protein_align(lifted_prot, ref_prot).identity
    except (LiftOnAlignmentError, Exception):
        return None


def _safe_dna_identity(lifted_dna, ref_dna):
    if not lifted_dna or not ref_dna:
        return None
    try:
        return align.trans_align(lifted_dna, ref_dna).identity
    except (LiftOnAlignmentError, Exception):
        return None


# ---------------------------------------------------------------------------
# reference materialization
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# Structural metrics (Phase 3). Coordinate-INDEPENDENT: the reference and the
# lifted transcript live on different genomes, so every comparison uses
# spliced-transcript boundary POSITIONS (cumulative segment lengths), never
# absolute genomic coordinates. These catch structural regressions that
# per-residue protein identity hides -- a rearranged intron chain or an ORF that
# lost its start/stop but still scores high identity on the aligned core.
# ---------------------------------------------------------------------------

def _internal_boundaries(features, strand):
    """Internal splice-junction positions in spliced 5'->3' coordinates: the
    cumulative segment length at each internal boundary. Coordinate-independent,
    so a reference transcript (one genome) and its lift (another) are directly
    comparable. <2 segments -> [] (no internal junctions)."""
    segs = sorted((min(f.start, f.end), max(f.start, f.end)) for f in features)
    if len(segs) < 2:
        return []
    lengths = [e - s + 1 for (s, e) in segs]
    if strand == "-":
        lengths = lengths[::-1]
    bounds, acc = [], 0
    for L in lengths[:-1]:
        acc += L
        bounds.append(acc)
    return bounds


def _boundary_snsp(ref_bounds, lifted_bounds):
    """(sn, sp, exact) for two internal-boundary lists. sn = fraction of the
    reference junctions recovered; sp = fraction of the lifted junctions that are
    correct; exact = the full intron chain matches (order + positions). Empty
    reference boundary set -> sn = 1.0 (vacuously recovered)."""
    rb, lb = set(ref_bounds), set(lifted_bounds)
    inter = len(rb & lb)
    sn = inter / len(rb) if rb else 1.0
    sp = inter / len(lb) if lb else 1.0
    return sn, sp, int(list(ref_bounds) == list(lifted_bounds))


def _orf_validity(prot):
    """(start_ok, stop_ok, no_internal_stop, valid) for a translated protein.
    start_ok = begins with M; stop_ok = ends with a stop codon (trailing '*');
    no_internal_stop = no '*' before the terminal; valid = all three. Empty
    protein -> all 0."""
    if not prot:
        return (0, 0, 0, 0)
    start_ok = int(prot.startswith("M"))
    stop_ok = int(prot.endswith("*"))
    body = prot[:-1] if prot.endswith("*") else prot
    no_internal = int("*" not in body)
    return (start_ok, stop_ok, no_internal,
            int(bool(start_ok and stop_ok and no_internal)))


def build_reference(ref_gff: str, ref_fa: str, log=print) -> tuple:
    """Return (ref, ref_index).

    ref: {transcript_id: {dna_exon, dna_cds, prot, is_coding}} keyed by the
    reference transcript-container ID.
    ref_index: {"ids_by_type", "all_ids", "census"} for all-feature-type
    completeness (see feature_index)."""
    db = _build_db(ref_gff)
    fa = pyfaidx.Fasta(ref_fa)
    ref: dict[str, dict] = {}
    for mrna in _transcript_features(db):
        exons = _children(db, mrna, "exon")
        cds = _children(db, mrna, ("CDS", "stop_codon"))
        cds_only = [c for c in cds]
        dna_exon = get_dna_sequence(mrna, fa, exons) if exons else ""
        dna_cds = get_dna_sequence(mrna, fa, cds_only) if cds_only else ""
        prot = get_protein_sequence(mrna, fa, cds_only) if cds_only else ""
        ref[mrna.id] = {
            "dna_exon": dna_exon or "",
            "dna_cds": dna_cds or "",
            "prot": (prot or "").rstrip("*") + ("*" if prot and prot.endswith("*") else ""),
            "is_coding": bool(cds_only) and bool(prot),
            # Phase 3 structural refs (coordinate-independent boundary positions).
            "exon_bounds": _internal_boundaries(exons, mrna.strand) if exons else [],
            "cds_bounds": _internal_boundaries(cds_only, mrna.strand) if cds_only else [],
            "n_exons": len(exons), "n_cds": len(cds_only),
        }
    n_coding = sum(1 for v in ref.values() if v["is_coding"])
    ref_ids_by_type, ref_all_ids, ref_census = feature_index(db)
    log(f"  [eval] reference: {len(ref)} transcript containers ({n_coding} coding); "
        f"{len(ref_census)} feature types, {sum(ref_census.values())} features total")
    ref_index = {
        "ids_by_type": ref_ids_by_type,
        "all_ids": ref_all_ids,
        "census": dict(ref_census),
    }
    return ref, ref_index


# ---------------------------------------------------------------------------
# per-tool evaluation
# ---------------------------------------------------------------------------

def _eval_key(d):
    """Ranking key for picking the best copy per reference id (max protein
    identity, then DNA identity)."""
    return (d.get("protein_identity") or -1, d.get("dna_identity") or -1)


def _eval_one_mrna(mrna, exons, cds, ref_id, ref, fa, is_miniprot):
    """Worker: extract lifted protein/DNA (pyfaidx + pre-fetched CDS/exon lists)
    and align them against the reference. No gffutils DB access here, so it is
    safe to run concurrently (parasail releases the GIL). Returns the per-mRNA
    record WITHOUT copy_count (assigned later, in submission order)."""
    r = ref[ref_id]
    lifted_prot = get_protein_sequence(mrna, fa, cds) if cds else ""
    if is_miniprot:
        lifted_dna = get_dna_sequence(mrna, fa, cds) if cds else ""
        ref_dna = r["dna_cds"]
        dna_basis = "cds"
    else:
        lifted_dna = get_dna_sequence(mrna, fa, exons) if exons else (
            get_dna_sequence(mrna, fa, cds) if cds else "")
        ref_dna = r["dna_exon"] or r["dna_cds"]
        dna_basis = "transcript" if exons else "cds"
    prot_id = _safe_prot_identity(lifted_prot, r["prot"]) if r["is_coding"] else None
    dna_id = _safe_dna_identity(lifted_dna, ref_dna)
    # Phase 3 structural metrics (coordinate-independent). For coding transcripts
    # the intron chain is compared on CDS boundaries (the protein-relevant
    # structure, present for all three tools incl. CDS-only miniprot); non-coding
    # falls back to exon boundaries. ORF-validity is coding-only.
    use_cds = bool(r["is_coding"])
    lifted_bounds = (_internal_boundaries(cds, mrna.strand) if (use_cds and cds)
                     else (_internal_boundaries(exons, mrna.strand) if exons else []))
    ref_bounds = r["cds_bounds"] if use_cds else r["exon_bounds"]
    sn, sp, exact = _boundary_snsp(ref_bounds, lifted_bounds)
    start_ok, stop_ok, no_internal, orf_valid = (
        _orf_validity(lifted_prot) if use_cds else (0, 0, 0, 0))
    return {
        "ref_mrna_id": ref_id,
        "tool_feature_id": mrna.id,
        "is_coding": int(r["is_coding"]),
        "protein_identity": prot_id,
        "dna_identity": dna_id,
        "ref_prot_len": len((r["prot"] or "").rstrip("*")),
        "lifted_prot_len": len((lifted_prot or "").rstrip("*")),
        "ref_dna_len": len(ref_dna or ""),
        "lifted_dna_len": len(lifted_dna or ""),
        "dna_basis": dna_basis,
        "seqid": mrna.seqid,
        # structural (Phase 3)
        "n_cds_lifted": len(cds), "n_cds_ref": r["n_cds"],
        "intron_chain_exact": exact,
        "exon_sn": round(sn, 5), "exon_sp": round(sp, 5),
        "orf_start_ok": start_ok, "orf_stop_ok": stop_ok, "orf_valid": orf_valid,
    }


def evaluate_tool(tool: str, tool_gff: str, tgt_fa: str, ref: dict,
                  manifest: dict, out_dir: Path, profile: dict | None,
                  log=print, ref_index: dict | None = None, threads: int = 1) -> dict:
    ref_ids = set(ref.keys())
    space = manifest.get("miniprot_target_space", "protein")
    acc_to_mrna = manifest.get("protein_acc_to_mrna", {})
    is_miniprot = tool == "miniprot"

    db = _build_db(tool_gff)
    fa = pyfaidx.Fasta(tgt_fa)

    target_map = id_mapping.build_miniprot_target_map(db) if is_miniprot else {}

    # all-feature-type completeness (independent of transcript-container scoring)
    cbt = {}
    if ref_index is not None:
        cbt = completeness_by_type(db, ref_index["ids_by_type"], ref_index["all_ids"])

    # --- phase 1: materialize per-mRNA gffutils reads on the parent thread ---
    # gffutils SQLite is not concurrent-safe; the alignment workers below only
    # touch the thread-safe pyfaidx Fasta + the pre-fetched exon/CDS lists.
    payloads = []   # (idx, mrna, exons, cds, ref_id) for valid ref_id only
    n_features = 0
    n_unmapped_id = 0
    for mrna in _transcript_features(db):
        n_features += 1
        if is_miniprot:
            ref_id = id_mapping.resolve_miniprot(mrna, target_map, space, acc_to_mrna, ref_ids)
        else:
            ref_id, _copy = id_mapping.resolve_liftoff_lifton(mrna, ref_ids)
        if ref_id is None or ref_id not in ref:
            n_unmapped_id += 1
            continue
        exons = _children(db, mrna, "exon")
        cds = _children(db, mrna, ("CDS", "stop_codon"))
        payloads.append((len(payloads), mrna, exons, cds, ref_id))

    # --- phase 2: extract + align (parallel; parasail releases the GIL) ---
    recs = [None] * len(payloads)
    if threads and threads > 1 and len(payloads) > 1:
        with ThreadPoolExecutor(max_workers=int(threads)) as ex:
            futs = {ex.submit(_eval_one_mrna, mrna, exons, cds, ref_id, ref, fa,
                              is_miniprot): idx
                    for (idx, mrna, exons, cds, ref_id) in payloads}
            for fut in as_completed(futs):
                recs[futs[fut]] = fut.result()
    else:
        for (idx, mrna, exons, cds, ref_id) in payloads:
            recs[idx] = _eval_one_mrna(mrna, exons, cds, ref_id, ref, fa, is_miniprot)

    # --- phase 3: assemble `best` per reference id in submission order ---
    # (deterministic dedup: identical to the serial tie-break, independent of
    #  worker completion order, so threads=1 and threads=N agree byte-for-byte.)
    best: dict[str, dict] = {}
    for rec in recs:
        ref_id = rec["ref_mrna_id"]
        prev = best.get(ref_id)
        if prev is None:
            rec["copy_count"] = 1
            best[ref_id] = rec
        else:
            rec["copy_count"] = prev["copy_count"] + 1
            if _eval_key(rec) > _eval_key(prev):
                best[ref_id] = rec
            else:
                prev["copy_count"] += 1

    # assemble per-reference rows (left join on the reference set)
    rows = []
    for ref_id, r in ref.items():
        rec = best.get(ref_id)
        if rec is None:
            rows.append({
                "ref_mrna_id": ref_id, "tool_feature_id": "", "recovered": 0,
                "is_coding": int(r["is_coding"]), "copy_index": 0,
                "protein_identity": "", "dna_identity": "",
                "ref_prot_len": len((r["prot"] or "").rstrip("*")), "lifted_prot_len": 0,
                "ref_dna_len": len(r["dna_exon"] or r["dna_cds"] or ""), "lifted_dna_len": 0,
                "dna_basis": "", "status": "lost", "seqid": "",
                "n_cds_lifted": 0, "n_cds_ref": r["n_cds"],
                "intron_chain_exact": "", "exon_sn": "", "exon_sp": "",
                "orf_start_ok": "", "orf_stop_ok": "", "orf_valid": "",
            })
        else:
            status = "recovered"
            if rec["is_coding"] and rec["protein_identity"] is None:
                status = "map_failed"
            rows.append({
                "ref_mrna_id": ref_id, "tool_feature_id": rec["tool_feature_id"],
                "recovered": 1, "is_coding": rec["is_coding"],
                "copy_index": rec["copy_count"] - 1,
                "protein_identity": rec["protein_identity"] if rec["protein_identity"] is not None else "",
                "dna_identity": rec["dna_identity"] if rec["dna_identity"] is not None else "",
                "ref_prot_len": rec["ref_prot_len"], "lifted_prot_len": rec["lifted_prot_len"],
                "ref_dna_len": rec["ref_dna_len"], "lifted_dna_len": rec["lifted_dna_len"],
                "dna_basis": rec["dna_basis"], "status": status, "seqid": rec["seqid"],
                "n_cds_lifted": rec["n_cds_lifted"], "n_cds_ref": rec["n_cds_ref"],
                "intron_chain_exact": rec["intron_chain_exact"],
                "exon_sn": rec["exon_sn"], "exon_sp": rec["exon_sp"],
                "orf_start_ok": rec["orf_start_ok"], "orf_stop_ok": rec["orf_stop_ok"],
                "orf_valid": rec["orf_valid"],
            })

    # write TSV
    out_dir.mkdir(parents=True, exist_ok=True)
    cols = ["ref_mrna_id", "tool_feature_id", "recovered", "is_coding", "copy_index",
            "protein_identity", "dna_identity", "ref_prot_len", "lifted_prot_len",
            "ref_dna_len", "lifted_dna_len", "dna_basis", "status", "seqid",
            "n_cds_lifted", "n_cds_ref", "intron_chain_exact", "exon_sn", "exon_sp",
            "orf_start_ok", "orf_stop_ok", "orf_valid"]
    tsv = out_dir / f"{tool}.transcripts.tsv"
    with open(tsv, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for row in rows:
            fh.write("\t".join(str(row[c]) for c in cols) + "\n")

    # per-feature-type completeness side TSV (does not disturb transcripts.tsv,
    # which reporter.py indexes positionally)
    if cbt:
        ft_tsv = out_dir / f"{tool}.feature_types.tsv"
        with open(ft_tsv, "w") as fh:
            fh.write("feature_type\tn_reference\tn_recovered\tpct_recovered"
                     "\tn_extra_copies\tn_tool_features\n")
            for ftype in sorted(cbt):
                d = cbt[ftype]
                fh.write(f"{ftype}\t{d['n_reference']}\t{d['n_recovered']}\t"
                         f"{d['pct_recovered']}\t{d['n_extra_copies']}\t"
                         f"{d['n_tool_features']}\n")

    # aggregate
    n_ref_total = len(ref)
    n_ref_coding = sum(1 for v in ref.values() if v["is_coding"])
    recovered = [r for r in rows if r["recovered"]]
    recovered_coding = [r for r in recovered if r["is_coding"]]
    extra_copies = sum(max(0, r["copy_index"]) for r in rows)
    prot_ids = [float(r["protein_identity"]) for r in recovered_coding
                if r["protein_identity"] != ""]
    dna_ids = [float(r["dna_identity"]) for r in recovered if r["dna_identity"] != ""]

    # Phase 3 structural aggregates over recovered coding transcripts. These are
    # eval-only additive summary keys (protein/dna identity + completeness are
    # unchanged), so downstream JSON consumers that read by key are unaffected.
    def _frac(key):
        vals = [r[key] for r in recovered_coding if r.get(key) not in ("", None)]
        return round(sum(int(v) for v in vals) / len(vals), 5) if vals else None

    def _mean_key(key):
        vals = [float(r[key]) for r in recovered_coding if r.get(key) not in ("", None)]
        return round(sum(vals) / len(vals), 5) if vals else None

    structural = {
        "orf_valid_frac": _frac("orf_valid"),
        "orf_start_frac": _frac("orf_start_ok"),
        "orf_stop_frac": _frac("orf_stop_ok"),
        "intron_chain_exact_frac": _frac("intron_chain_exact"),
        "exon_sn_mean": _mean_key("exon_sn"),
        "exon_sp_mean": _mean_key("exon_sp"),
        "n_scored": len(recovered_coding),
    }

    def _stats(vals):
        if not vals:
            return {"mean": None, "median": None, "pct_identical": None, "n": 0}
        return {
            "mean": round(statistics.fmean(vals), 5),
            "median": round(statistics.median(vals), 5),
            "pct_identical": round(sum(1 for v in vals if v >= 0.999999) / len(vals), 5),
            "n": len(vals),
        }

    def _hist(vals, bins=10):
        h = [0] * bins
        for v in vals:
            idx = min(bins - 1, int(v * bins))
            h[idx] += 1
        return h

    summary = {
        "tool": tool, "benchmark": manifest["id"], "species": manifest["species"],
        "cross_species": manifest["cross_species"],
        "n_reference_total": n_ref_total, "n_reference_coding": n_ref_coding,
        "n_tool_features": n_features, "n_unmapped_id": n_unmapped_id,
        "n_recovered_any": len(recovered), "n_recovered_coding": len(recovered_coding),
        "completeness_all": round(len(recovered) / n_ref_total, 5) if n_ref_total else None,
        "completeness_coding": round(len(recovered_coding) / n_ref_coding, 5) if n_ref_coding else None,
        "completeness_by_type": cbt,
        # headline all-feature completeness is id-recovery based, so it is n/a for
        # miniprot (its MP* ids never match reference ids — see completeness_coding)
        "completeness_feature_total": (
            None if is_miniprot else cbt.get("_overall_", {}).get("pct_recovered")),
        "n_extra_copies": extra_copies,
        "dna_basis": "cds" if is_miniprot else "transcript",
        "protein_identity": _stats(prot_ids),
        "dna_identity": _stats(dna_ids),
        "protein_identity_hist": _hist(prot_ids),
        "structural": structural,
        "profile": {
            "wall_clock_seconds": (profile or {}).get("wall_clock_seconds"),
            "peak_rss_mb": (profile or {}).get("peak_rss_mb"),
        },
        "transcripts_tsv": str(tsv),
    }
    (out_dir / f"{tool}.summary.json").write_text(json.dumps(summary, indent=2))
    log(f"  [eval:{tool}] completeness_coding="
        f"{summary['completeness_coding']} mean_prot_id={summary['protein_identity']['mean']} "
        f"mean_dna_id={summary['dna_identity']['mean']} (mapped {len(recovered_coding)}/{n_ref_coding})")
    return summary


def evaluate_all(manifest: dict, work_dir: Path, profiles: dict, force: bool,
                 log=print, mode="genes", threads: int = 1) -> dict:
    suf = "" if mode == "genes" else f"__{mode}"
    tools_sub = "tools" if mode == "genes" else f"tools_{mode}"
    eval_sub = "eval" if mode == "genes" else f"eval_{mode}"
    done = work_dir / ".done" / f"eval{suf}.done"
    eval_dir = work_dir / eval_sub
    out_path = eval_dir / "summaries.json"
    if done.exists() and not force and out_path.exists():
        log(f"  [eval:{mode}] cached")
        return json.loads(out_path.read_text())

    p = manifest["paths"]
    ref, ref_index = build_reference(p["ref_gff"], p["ref_fa"], log)
    tool_gffs = {
        "liftoff": work_dir / tools_sub / "liftoff" / "liftoff.gff3",
        # miniprot is mode-independent: always the shared genes-mode output
        "miniprot": work_dir / "tools" / "miniprot" / "miniprot.gff3",
        "lifton": work_dir / tools_sub / "lifton" / "lifton.gff3",
    }
    # 4th variant (genes-mode only): LiftOn with --optimize. evaluate_tool keys
    # purely on the tool-name string, so it scores identically to "lifton"
    # (resolve_liftoff_lifton, transcript DNA basis).
    if mode == "genes":
        tool_gffs["lifton_optimize"] = work_dir / "tools" / "lifton_optimize" / "lifton.gff3"
    # present tools, preserving the liftoff/miniprot/lifton order
    present = []
    for tool, gff in tool_gffs.items():
        if gff.exists():
            present.append((tool, gff))
        else:
            log(f"  [eval:{mode}] SKIP {tool}: {gff} missing")

    def _run_tool(tool, gff):
        return evaluate_tool(
            tool, str(gff), p["tgt_fa"], ref, manifest, eval_dir,
            (profiles or {}).get(tool), log, ref_index=ref_index, threads=threads)

    # The 3 tools are independent — separate gffutils DB files (``<gff>_db``),
    # separate output files, and shared READ-ONLY ref/ref_index — so run them
    # concurrently to overlap their serial DB builds + alignment. Each still uses
    # ``threads`` align workers (3*threads stays within the core budget).
    # Determinism: assemble ``summaries`` in fixed tool order, not completion order.
    results = {}
    _concurrent_tools = (threads and threads > 1 and len(present) > 1
                         and not os.environ.get("LIFTON_EVAL_SERIAL_TOOLS"))
    if _concurrent_tools:
        with ThreadPoolExecutor(max_workers=len(present)) as ex:
            futs = {ex.submit(_run_tool, tool, gff): tool for tool, gff in present}
            for fut in as_completed(futs):
                results[futs[fut]] = fut.result()
    else:
        for tool, gff in present:
            results[tool] = _run_tool(tool, gff)

    summaries = {tool: results[tool] for tool, _ in present}
    out_path.write_text(json.dumps(summaries, indent=2))
    done.parent.mkdir(parents=True, exist_ok=True)
    done.write_text("ok\n")
    return summaries
