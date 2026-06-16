#!/usr/bin/env python
"""Step-3 gffutils query-collapse micro-bench A/B (Iteration 18).

Step 3 (`lifton/extract_sequence.py`) extracted transcript DNA + protein
FASTAs by issuing up to **three** `gffutils.children()` queries per feature:
(a) `featuretype='exon'`, (b) `featuretype=('start_codon','CDS','stop_codon')`,
and (c) — only when a+b are both empty — the recursive-descent
`level=1, order_by='start'`. Iteration 18 collapses these to **one** ordered
query + an in-Python featuretype partition (`_stream_inner` /
`__inner_extract_feature`). The change is byte-neutral (the exon/CDS lists are
re-sorted by start in `merge_children_intervals` before concatenation, so SQL
row order is discarded for those branches; the recursive-descent branch keeps
`order_by='start'`).

This micro-bench times **only Step 3** (the rest of the pipeline would dilute
the small signal). It compares, on a real benchmark reference DB built ONCE and
reused across all repeats:

  * baseline  = the pre-Iteration-18 3-query extractor, loaded verbatim from
                `git show HEAD:lifton/extract_sequence.py` (HEAD = the Iter-17
                commit, which still carries the 3-query code) — no
                re-implementation, so the baseline is byte-exact.
  * branch    = the working-tree 1-query extractor (the current import).

Two metrics:
  1. children() call count — DETERMINISTIC (noise-free) proof of the
     round-trip reduction (the optimization's mechanism).
  2. wall-clock median ± IQR over N repeats (alternated to cancel drift) — the
     GO/NO-GO measurement. The payoff is small, so the gate is: byte-identical
     FASTAs AND a median wall reduction clearly above the run-to-run noise band.

Decision rule (per the Iteration-18 plan):
  GO   iff byte-identical AND median(baseline) - median(branch) > 2 * noise_band
  NO-GO otherwise  -> revert lifton/extract_sequence.py, keep this script +
                      results as the audit trail (Iter-9/11/13/15 precedent).

Usage (repo root, lifton_devel env):
    python -m benchmarks.compare.extract_query_collapse_ab [IDS...]
    EXTRACT_AB_REPEATS=25 python -m benchmarks.compare.extract_query_collapse_ab drosophila
"""
from __future__ import annotations

import importlib.util
import json
import os
import statistics
import subprocess
import sys
import tempfile
import time
from pathlib import Path

from pyfaidx import Fasta

from lifton import annotation, extract_sequence
from .tool_runners import _clean_input_dbs

HERE = Path(__file__).resolve().parent
WORK = HERE / "work"
REPO = HERE.parent.parent  # repo root (…/LiftOn)

DEFAULT_IDS = ["drosophila", "mouse_to_rat", "rice"]
REPEATS = int(os.environ.get("EXTRACT_AB_REPEATS", "15"))
WARMUP = 2
FEATURES = ["gene"]  # the -f-less default parent type (get_parent_features_to_lift(None))


def _load_baseline():
    """Import the pre-Iteration-18 3-query extractor verbatim from git HEAD.

    HEAD is the Iteration-17 commit (the Step-3 collapse is uncommitted in the
    working tree), so HEAD:lifton/extract_sequence.py is the exact 3-query
    baseline — no hand re-implementation, no mismatch risk.
    """
    src = subprocess.check_output(
        ["git", "show", "HEAD:lifton/extract_sequence.py"], cwd=str(REPO)
    ).decode("utf-8")
    if "all_children = list(ref_db.db_connection.children(" in src:
        raise RuntimeError(
            "HEAD already contains the collapsed (1-query) extractor; this "
            "A/B needs HEAD to be the pre-collapse baseline. Run before "
            "committing Iteration 18, or point at the right baseline ref."
        )
    tmp = tempfile.NamedTemporaryFile(
        suffix="_extract_baseline.py", delete=False, mode="w")
    tmp.write(src)
    tmp.close()
    spec = importlib.util.spec_from_file_location(
        "extract_sequence_baseline", tmp.name)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod, tmp.name


def _count_children_calls(extract_fn, ref_db, ref_fai, out_dir):
    """Run `extract_fn` once, counting how many times db_connection.children
    is invoked (= SQLite round-trips on the per-feature path)."""
    conn = ref_db.db_connection
    orig = conn.children
    n = [0]

    def wrapped(*a, **k):
        n[0] += 1
        return orig(*a, **k)

    conn.children = wrapped
    try:
        extract_fn(ref_db, FEATURES, ref_fai, out_dir)
    finally:
        conn.children = orig
    return n[0]


def _time_once(extract_fn, ref_db, ref_fai, out_dir):
    t0 = time.perf_counter()
    extract_fn(ref_db, FEATURES, ref_fai, out_dir)
    return time.perf_counter() - t0


def _iqr(xs):
    if len(xs) < 4:
        return max(xs) - min(xs) if xs else 0.0
    q = statistics.quantiles(xs, n=4)  # [p25, p50, p75]
    return q[2] - q[0]


def _fasta_bytes(out_dir):
    tp = Path(out_dir) / "transcripts.fa"
    pp = Path(out_dir) / "proteins.fa"
    return tp.read_bytes(), pp.read_bytes()


def _bench_one(bid, baseline_extract, branch_extract):
    man = json.loads(
        (WORK / bid / "subset" / "subset.manifest.json").read_text())
    p = man["paths"]
    ref_gff, ref_fa = p["ref_gff"], p["ref_fa"]

    # Build the reference DB ONCE (excluded from timing) and reuse for every
    # repeat of BOTH arms, so we time only the extraction queries, not the
    # gffutils create_db. Clean any stale sibling _db first for a clean build.
    _clean_input_dbs(ref_gff)
    print(f"  [{bid}] building reference Annotation (once)…", flush=True)
    ref_db = annotation.Annotation(ref_gff, False, False)
    ref_fai = Fasta(ref_fa)

    with tempfile.TemporaryDirectory(prefix=f"extract_ab_{bid}_") as td:
        td = Path(td)
        dir_base = td / "baseline"
        dir_branch = td / "branch"
        scratch = td / "scratch"

        # 1) Byte-identity (also primes the page/SQLite cache).
        baseline_extract(ref_db, FEATURES, ref_fai, str(dir_base))
        branch_extract(ref_db, FEATURES, ref_fai, str(dir_branch))
        b_tx, b_pr = _fasta_bytes(dir_base)
        n_tx, n_pr = _fasta_bytes(dir_branch)
        byte_identical = (b_tx == n_tx) and (b_pr == n_pr)

        # 2) Deterministic children() call-count (the mechanism proof).
        q_base = _count_children_calls(baseline_extract, ref_db, ref_fai, str(scratch))
        q_branch = _count_children_calls(branch_extract, ref_db, ref_fai, str(scratch))

        # 3) Wall-clock, N repeats, arms ALTERNATED to cancel any drift.
        for _ in range(WARMUP):
            _time_once(baseline_extract, ref_db, ref_fai, str(scratch))
            _time_once(branch_extract, ref_db, ref_fai, str(scratch))
        t_base, t_branch = [], []
        for i in range(REPEATS):
            if i % 2 == 0:
                t_base.append(_time_once(baseline_extract, ref_db, ref_fai, str(scratch)))
                t_branch.append(_time_once(branch_extract, ref_db, ref_fai, str(scratch)))
            else:
                t_branch.append(_time_once(branch_extract, ref_db, ref_fai, str(scratch)))
                t_base.append(_time_once(baseline_extract, ref_db, ref_fai, str(scratch)))

    med_base = statistics.median(t_base)
    med_branch = statistics.median(t_branch)
    iqr_base, iqr_branch = _iqr(t_base), _iqr(t_branch)
    noise = max(iqr_base, iqr_branch)
    delta = med_base - med_branch
    speedup = med_base / med_branch if med_branch > 0 else None
    saved_pct = 100.0 * delta / med_base if med_base > 0 else None
    wall_go = byte_identical and delta > 2 * noise

    rec = {
        "benchmark": bid,
        "repeats": REPEATS,
        "byte_identical": byte_identical,
        "children_calls": {"baseline": q_base, "branch": q_branch,
                           "reduction_x": round(q_base / q_branch, 2) if q_branch else None},
        "wall_ms": {
            "baseline_median": round(med_base * 1e3, 3),
            "branch_median": round(med_branch * 1e3, 3),
            "baseline_iqr": round(iqr_base * 1e3, 3),
            "branch_iqr": round(iqr_branch * 1e3, 3),
            "noise_band": round(noise * 1e3, 3),
            "delta_median": round(delta * 1e3, 3),
        },
        "speedup": round(speedup, 3) if speedup else None,
        "saved_pct": round(saved_pct, 2) if saved_pct is not None else None,
        "wall_gate": {
            "rule": "byte_identical AND delta_median > 2 * noise_band",
            "delta_gt_2x_noise": bool(delta > 2 * noise),
            "GO": bool(wall_go),
        },
    }
    return rec


def main(argv=None):
    ids = (argv if argv is not None else sys.argv[1:]) or DEFAULT_IDS
    baseline_mod, tmp_path = _load_baseline()
    baseline_extract = baseline_mod.extract_features_to_fasta
    branch_extract = extract_sequence.extract_features_to_fasta
    try:
        results = []
        for bid in ids:
            print(f"=== {bid}: Step-3 query-collapse micro-bench "
                  f"(N={REPEATS} repeats) ===", flush=True)
            rec = _bench_one(bid, baseline_extract, branch_extract)
            results.append(rec)
            cc = rec["children_calls"]
            w = rec["wall_ms"]
            print(f"  [{bid}] byte-identical={rec['byte_identical']} | "
                  f"children() calls {cc['baseline']}->{cc['branch']} "
                  f"({cc['reduction_x']}x fewer)", flush=True)
            print(f"  [{bid}] Step-3 wall median {w['baseline_median']}ms -> "
                  f"{w['branch_median']}ms (Δ{w['delta_median']}ms, "
                  f"{rec['speedup']}x / saved {rec['saved_pct']}%) | "
                  f"noise band ±{w['noise_band']}ms | "
                  f"wall-GO={rec['wall_gate']['GO']}", flush=True)
            if not rec["byte_identical"]:
                print(f"  !! WARNING {bid}: FASTAs differ — byte-neutrality "
                      f"FAILED; the collapse is NOT safe, investigate.",
                      flush=True)
    finally:
        try:
            os.unlink(tmp_path)
        except OSError:
            pass

    (HERE / "extract_query_collapse_ab.json").write_text(
        json.dumps(results, indent=2))
    _write_md(results, HERE / "extract_query_collapse_ab.md")
    print("\nDONE: wrote extract_query_collapse_ab.json + "
          "extract_query_collapse_ab.md", flush=True)
    all_byte = all(r["byte_identical"] for r in results)
    any_wall_go = any(r["wall_gate"]["GO"] for r in results)
    print(f"SUMMARY: byte-identical all={all_byte}; any wall-GO={any_wall_go}",
          flush=True)
    return 0


def _write_md(results, path):
    lines = [
        "## Step-3 gffutils query-collapse micro-bench — Iteration 18\n",
        "Step 3 issued up to **3** `gffutils.children()` queries per feature "
        "(exon; CDS-tuple; recursive-descent fallback). Iteration 18 collapses "
        "them to **one** `order_by='start'` query + an in-Python featuretype "
        "partition (byte-neutral: `merge_children_intervals` re-sorts the "
        "exon/CDS lists before concatenation, so SQL order is discarded for "
        "those branches; the recursive-descent branch keeps `order_by='start'`).\n",
        "This times **only Step 3** on a real reference DB built ONCE and "
        "reused. `baseline` = the pre-Iter-18 3-query extractor loaded verbatim "
        "from `git show HEAD:lifton/extract_sequence.py`; `branch` = the "
        "working-tree 1-query extractor. The `children()` call count is the "
        "deterministic (noise-free) mechanism proof; the wall-clock median±IQR "
        "is the GO/NO-GO measurement.\n",
        "**Gate:** byte-identical FASTAs AND median wall reduction > 2× the "
        "run-to-run noise band.\n",
        "| Dataset | byte-ident | children() calls | reduction | Step-3 wall median | speedup | wall-GO |",
        "|---|---|---|---|---|---|---|",
    ]
    for r in results:
        cc = r["children_calls"]
        w = r["wall_ms"]
        lines.append(
            "| {} | {} | {}→{} | {}× | {}ms→{}ms (Δ{}ms, ±{} noise) | {}× ({}%) | {} |"
            .format(
                r["benchmark"],
                "yes" if r["byte_identical"] else "**NO**",
                cc["baseline"], cc["branch"], cc["reduction_x"],
                w["baseline_median"], w["branch_median"],
                w["delta_median"], w["noise_band"],
                r["speedup"], r["saved_pct"],
                "**yes**" if r["wall_gate"]["GO"] else "no",
            ))
    lines.append(
        "\n**Interpretation:** the `children()` call-count reduction is exact "
        "and deterministic (the optimization's mechanism). The wall-clock gate "
        "decides GO/NO-GO: if the median Step-3 reduction does not clearly "
        "exceed the noise band on any dataset, the change is byte-correct but "
        "unmeasurable → NO-GO + revert, keeping this script + results as the "
        "audit trail (the Iter-9/11/13/15 precedent).")
    path.write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    sys.exit(main())
