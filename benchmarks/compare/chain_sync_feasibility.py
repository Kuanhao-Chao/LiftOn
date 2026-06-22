#!/usr/bin/env python
"""Iteration-4 #1 feasibility — measure recoverable identity from relaxing the
exact-genomic-end sync constraint (protein_maximization.py:259), READ-ONLY.

The headroom diagnostic showed ~33% of drosophila merge-fired transcripts are
"single-chunk": no sync point ever fires, so the chaining awards the WHOLE
protein all-or-nothing to one source. This prototype answers — BEFORE we build
the risky tiling reconciliation — *how much protein identity finer mixing could
actually recover*.

It runs LiftOn IN-PROCESS (so it can monkeypatch `chaining_algorithm`) on the
cached `-L`/`-M`, serial (`-t 1`). For every merge-fired transcript it records
the two input alignments and computes, all in REFERENCE-AA coordinates:
  - whole_l / whole_m : fraction of ref AAs each source matches alone.
  - current_est       : max(whole_l, whole_m) — what a single-chunk merge gets
                        (exact for single-chunk; a lower bound for multi-chunk).
  - oracle_either     : fraction of ref AAs matched by Liftoff OR miniprot — the
                        per-position UPPER BOUND of any mixing strategy.
  - oracle_region     : take the better source per Liftoff-CDS ref-AA region —
                        the realistic ceiling for CDS-granular sync-relaxation.
  - gap_either / gap_region vs current_est = recoverable headroom.

Go/no-go: if the single-chunk gap is ~0, even perfect mixing can't help (both
sources match the same positions) → STOP, #1 isn't worth the tiling risk. If it
is material, proceed to design the tiling-safe relaxation.

Usage (repo root, lifton_devel env):  one bid per process (fresh monkeypatch):
    python -m benchmarks.compare.chain_sync_feasibility <bid>
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path

from lifton import protein_maximization as pm
from lifton import get_id_fraction
from lifton import lifton as lifton_main

from .tool_runners import _clean_input_dbs, _compose_env  # noqa: F401 (env parity)

HERE = Path(__file__).resolve().parent
WORK = HERE / "work"
REG = json.loads((HERE / "benchmarks.json").read_text())
TOOLS = REG["tools"]

RECORDS = []
_ORIG = pm.chaining_algorithm


def _matched_ref_positions(ref_aln, query_aln):
    """Set of reference-AA indices (over non-gap ref columns) the query matches,
    plus the reference protein length (non-gap ref count)."""
    matched = set()
    ref_i = -1
    for r, q in zip(ref_aln, query_aln):
        if r != "-":
            ref_i += 1
            if q == r:
                matched.add(ref_i)
    return matched, ref_i + 1


def _aln_col_to_refaa(ref_aln, col):
    """Convert an alignment-column index to a reference-AA index (non-gap count)."""
    c = int(math.ceil(col))
    return sum(1 for ch in ref_aln[:c] if ch != "-")


def _boundary_items(bounds):
    """cdss_protein_aln_boundaries is a dict {idx: (start, end)} in the real
    pipeline (align.get_cdss_protein_boundary) but a list of (start, end) in the
    unit-test fixtures. Normalise to a list of (start, end) tuples."""
    if isinstance(bounds, dict):
        return [bounds[k] for k in sorted(bounds)]
    return list(bounds)


def _region_oracle(l_aln, l_matched, m_matched, ref_len):
    """Take the better source per Liftoff-CDS reference-AA region (the realistic
    ceiling for CDS-granular sync-relaxation)."""
    bounds = _boundary_items(getattr(l_aln, "cdss_protein_aln_boundaries", None) or [])
    if not bounds or ref_len <= 0:
        return (max(len(l_matched), len(m_matched)) / ref_len) if ref_len else 0.0
    edges = sorted({0} | {_aln_col_to_refaa(l_aln.ref_aln, b) for _, b in bounds} | {ref_len})
    total = 0
    for lo, hi in zip(edges, edges[1:]):
        lc = sum(1 for i in l_matched if lo <= i < hi)
        mc = sum(1 for i in m_matched if lo <= i < hi)
        total += max(lc, mc)
    return total / ref_len if ref_len else 0.0


def _analyze(l_aln, m_aln, chains):
    if not l_aln.ref_aln or not m_aln.ref_aln:
        return None
    l_matched, l_len = _matched_ref_positions(l_aln.ref_aln, l_aln.query_aln)
    m_matched, m_len = _matched_ref_positions(m_aln.ref_aln, m_aln.query_aln)
    ref_len = max(l_len, m_len)
    if ref_len <= 0:
        return None
    whole_l = len(l_matched) / ref_len
    whole_m = len(m_matched) / ref_len
    current = max(whole_l, whole_m)
    oracle_either = len(l_matched | m_matched) / ref_len
    oracle_region = _region_oracle(l_aln, l_matched, m_matched, ref_len)
    return {
        "n_chunks": len(chains),
        "ref_len": ref_len,
        "whole_l": round(whole_l, 4),
        "whole_m": round(whole_m, 4),
        "current_est": round(current, 4),
        "oracle_either": round(oracle_either, 4),
        "oracle_region": round(oracle_region, 4),
        "gap_either": round(oracle_either - current, 4),
        "gap_region": round(oracle_region - current, 4),
    }


def _instrumented(l_aln, m_aln, fai, DEBUG):
    cds_list, chains = _ORIG(l_aln, m_aln, fai, DEBUG)
    try:
        rec = _analyze(l_aln, m_aln, chains)
        if rec is not None:
            RECORDS.append(rec)
    except Exception as exc:  # never break the pipeline for instrumentation
        RECORDS.append({"error": repr(exc)})
    return cds_list, chains


def _ann_db(bid):
    for b in REG["benchmarks"]:
        if b["id"] == bid:
            return b.get("annotation_database", "RefSeq")
    return "RefSeq"


def _summarize(records, label):
    errors = [r for r in records if "error" in r]
    recs = [r for r in records if "error" not in r]
    single = [r for r in recs if r["n_chunks"] <= 1]
    multi = [r for r in recs if r["n_chunks"] > 1]

    def agg(rs, key):
        vals = [r[key] for r in rs]
        return round(sum(vals) / len(vals), 5) if vals else 0.0

    def n_gt(rs, key, thr):
        return sum(1 for r in rs if r[key] > thr)

    out = {
        "benchmark": label,
        "n_merge_fired": len(recs),
        "n_errors": len(errors),
        "sample_error": errors[0]["error"] if errors else None,
        "n_single_chunk": len(single),
        "n_multi_chunk": len(multi),
        "single_chunk": {
            "mean_current_est": agg(single, "current_est"),
            "mean_oracle_either": agg(single, "oracle_either"),
            "mean_oracle_region": agg(single, "oracle_region"),
            "mean_gap_either": agg(single, "gap_either"),
            "mean_gap_region": agg(single, "gap_region"),
            "n_gap_region_gt_1e-3": n_gt(single, "gap_region", 1e-3),
            "n_gap_region_gt_1e-2": n_gt(single, "gap_region", 1e-2),
            "n_gap_either_gt_1e-2": n_gt(single, "gap_either", 1e-2),
        },
        "all_merge_fired": {
            "mean_gap_region": agg(recs, "gap_region"),
            "n_gap_region_gt_1e-2": n_gt(recs, "gap_region", 1e-2),
        },
    }
    return out, single


def main(argv=None):
    args_in = (argv if argv is not None else sys.argv[1:])
    if not args_in:
        print("usage: chain_sync_feasibility <bid>", file=sys.stderr)
        return 2
    bid = args_in[0]
    pm.chaining_algorithm = _instrumented   # monkeypatch BEFORE the run

    W = WORK / bid
    man = json.loads((W / "subset" / "subset.manifest.json").read_text())
    p = man["paths"]
    root = W / "_chain_sync_feasibility"
    root.mkdir(parents=True, exist_ok=True)
    liftoff = W / "tools" / "liftoff" / "liftoff.gff3"
    miniprot = W / "tools" / "miniprot" / "miniprot.gff3"
    _clean_input_dbs(p["ref_gff"], liftoff, miniprot)

    argv_lifton = ["-t", "1", "-copies", "-ad", _ann_db(bid),
                   "-g", p["ref_gff"], "-L", str(liftoff), "-M", str(miniprot),
                   "-o", str(root / "out.gff3"), p["tgt_fa"], p["ref_fa"]]
    print(f"=== {bid}: in-process LiftOn (instrumented chaining) ===", flush=True)
    parsed = lifton_main.parse_args(argv_lifton)
    lifton_main.run_all_lifton_steps(parsed)

    summary, single = _summarize(RECORDS, bid)
    # top recoverable single-chunk transcripts (by region gap)
    summary["top_recoverable_single_chunk"] = sorted(
        ({"ref_len": r["ref_len"], "current_est": r["current_est"],
          "oracle_region": r["oracle_region"], "gap_region": r["gap_region"]}
         for r in single), key=lambda r: -r["gap_region"])[:10]

    (HERE / f"chain_sync_feasibility.{bid}.json").write_text(json.dumps(summary, indent=2))
    s = summary["single_chunk"]
    print(f"\n  [{bid}] merge-fired={summary['n_merge_fired']} "
          f"single-chunk={summary['n_single_chunk']}", flush=True)
    print(f"    single-chunk current={s['mean_current_est']} "
          f"oracle_either={s['mean_oracle_either']} oracle_region={s['mean_oracle_region']}",
          flush=True)
    print(f"    single-chunk mean gap_region={s['mean_gap_region']} "
          f"gap_either={s['mean_gap_either']}; "
          f"#gap_region>1e-2={s['n_gap_region_gt_1e-2']} "
          f"#gap_region>1e-3={s['n_gap_region_gt_1e-3']}", flush=True)
    print(f"\nDONE: wrote chain_sync_feasibility.{bid}.json", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
