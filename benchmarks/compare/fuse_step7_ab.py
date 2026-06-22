#!/usr/bin/env python
"""A/B for the FUSED Step-7 pool (Iteration 10).

Pre-Iteration-10, the parallel Step-7 path ran two sequential phases with a
hard barrier: a 4-thread prefetcher pool built ALL `MaterialisedLocus`
payloads (SQLite-bound), then a separate N-thread worker pool ran ALL the
parasail processing — wall = `T_materialise(4) + T_process(N)`, never
overlapping. Iteration 10 FUSES them into one pool: each worker materialises
its own locus (thread-local DB) then immediately processes it, so SQLite I/O
overlaps GIL-released parasail → wall → ~`max(T_mat, T_proc)`.

The change is BYTE-NEUTRAL (same payloads → same `process_locus_native` →
same `LocusResult`, ordered writer unchanged). This A/B proves both halves:
  - **byte-identity** (the HARD gate): fused vs two-phase are identical;
  - **speed**: fused is materially faster at Step 7 than two-phase.

Both states run on the same cached `-L`/`-M` (Step 4 skipped → isolates the
Step-7 dispatch + dodges the Liftoff `-copies` fresh-run non-determinism) with
the SAME argv (`-t N --locus-pipeline`); only the env gate differs:

  state "two_phase"  LIFTON_FUSE_STEP7=0   prefetcher pool + barrier (baseline)
  state "fused"      (default)             fused single pool        ← the iteration

`LIFTON_FUSE_STEP7=0` restores the exact pre-fusion path IN THE SAME BUILD, so
no second checkout is needed. Step-7 wall is probed via `LIFTON_PERF_STEP7`.
If `fused` ever regresses on a materialise-heavy dataset (rice), rerun that id
with `LIFTON_FUSE_MAT_CONCURRENCY=4` to test capping the materialise half.

Usage (repo root, lifton_devel env):
    python -m benchmarks.compare.fuse_step7_ab [IDS...]
    LIFTON_AB_THREADS=4 python -m benchmarks.compare.fuse_step7_ab drosophila
"""
from __future__ import annotations

import json
import os
import sys
from pathlib import Path

from .profiling import run_profiled
from .tool_runners import _clean_input_dbs, _compose_env
from .align_window_ab import _ann_db, _TRANSCRIPT_TYPES

HERE = Path(__file__).resolve().parent
WORK = HERE / "work"
REG = json.loads((HERE / "benchmarks.json").read_text())
TOOLS = REG["tools"]

# Step-7-dominated subsets (enough loci for the fan-out + materialise overlap
# to matter). rice is the most materialise-heavy → the contention canary.
DEFAULT_IDS = ["drosophila", "mouse_to_rat", "rice"]
THREADS = os.environ.get("LIFTON_AB_THREADS", "8")

STATES = {
    "two_phase": {"LIFTON_FUSE_STEP7": "0"},   # pre-Iteration-10 barrier baseline
    "fused":     {},                            # default = fused single pool
}


def _run(bid, state, env_extra, p, root):
    liftoff_gff = WORK / bid / "tools" / "liftoff" / "liftoff.gff3"
    miniprot_gff = WORK / bid / "tools" / "miniprot" / "miniprot.gff3"
    out_dir = root / state
    out_dir.mkdir(parents=True, exist_ok=True)
    out = out_dir / f"{state}.gff3"
    # Rebuild the input gffutils sibling DBs fresh per state (inputs unchanged
    # → identical rebuild) so the comparison is deterministic.
    _clean_input_dbs(p["ref_gff"], liftoff_gff, miniprot_gff)
    # SAME argv for every state: only the LIFTON_FUSE_STEP7 env gate differs.
    argv = [TOOLS["lifton_bin"], "-t", THREADS, "-copies", "-ad", _ann_db(bid),
            "-g", p["ref_gff"], "-L", str(liftoff_gff), "-M", str(miniprot_gff),
            "-o", str(out), "--locus-pipeline", p["tgt_fa"], p["ref_fa"]]
    env = _compose_env(TOOLS)
    env["LIFTON_PERF_STEP7"] = "1"   # emit the Step-7 wall probe to stderr
    env.update(env_extra)
    pr = run_profiled(argv, label=f"fuse_{state}",
                      log_dir=root / "logs", env=env, log=print)
    if pr.exit_code != 0 or not out.exists() or out.stat().st_size == 0:
        raise RuntimeError(f"{bid}: state {state} failed (exit {pr.exit_code}); "
                           f"see {pr.stderr_path}")
    return out, pr


def _parse_step7_wall(stderr_path):
    """Pull the LIFTON_PERF_STEP7 probe line ('[LiftOn][perf] Step7 wall
    (mode): N.NNs') — the Step-4-free measure of the fusion speedup."""
    try:
        for ln in Path(stderr_path).read_text().splitlines():
            if "[LiftOn][perf] Step7 wall" in ln:
                return float(ln.rsplit(":", 1)[1].strip().rstrip("s"))
    except (OSError, ValueError):
        pass
    return None


def _feature_counts(gff_path):
    n_gene = n_tx = n_cds = n_total = 0
    for ln in Path(gff_path).read_text().splitlines():
        if not ln or ln.startswith("#"):
            continue
        c = ln.split("\t")
        if len(c) < 9:
            continue
        n_total += 1
        t = c[2]
        if t == "gene" or t.endswith("_gene") or t == "pseudogene":
            n_gene += 1
        elif t in _TRANSCRIPT_TYPES:
            n_tx += 1
        elif t == "CDS":
            n_cds += 1
    return {"gene": n_gene, "transcript": n_tx, "CDS": n_cds, "total": n_total}


def main(argv=None):
    ids = (argv if argv is not None else sys.argv[1:]) or DEFAULT_IDS
    results = []
    for bid in ids:
        print(f"=== {bid}: fused Step 7 vs two-phase A/B "
              f"(-t {THREADS}, cached -L/-M) ===", flush=True)
        man = json.loads((WORK / bid / "subset" / "subset.manifest.json").read_text())
        p = man["paths"]
        root = WORK / bid / "_fuse_step7_ab"
        root.mkdir(parents=True, exist_ok=True)

        outs, prs = {}, {}
        for state, env_extra in STATES.items():
            print(f"--- {bid} state {state} ---", flush=True)
            outs[state], prs[state] = _run(bid, state, env_extra, p, root)

        identical = outs["two_phase"].read_bytes() == outs["fused"].read_bytes()
        counts = {s: _feature_counts(outs[s]) for s in STATES}
        completeness_ok = counts["two_phase"] == counts["fused"]

        wall = {s: prs[s].wall_clock_seconds for s in STATES}
        s7 = {s: _parse_step7_wall(prs[s].stderr_path) for s in STATES}
        rss = {s: getattr(prs[s], "peak_rss_mb", 0) or 0 for s in STATES}

        s7_speedup = (round(s7["two_phase"] / max(s7["fused"], 1e-9), 2)
                      if (s7["two_phase"] and s7["fused"]) else None)
        s7_saved_pct = (round(100.0 * (s7["two_phase"] - s7["fused"]) /
                              max(s7["two_phase"], 1e-9), 1)
                        if (s7["two_phase"] and s7["fused"]) else None)
        total_speedup = round(wall["two_phase"] / max(wall["fused"], 1e-9), 2)
        total_saved_pct = round(100.0 * (wall["two_phase"] - wall["fused"]) /
                                max(wall["two_phase"], 1e-9), 1)
        rss_overhead_pct = round(100.0 * (rss["fused"] - rss["two_phase"]) /
                                 max(rss["two_phase"], 1e-9), 1)

        rec = {
            "benchmark": bid,
            "threads": THREADS,
            "byte_identical": identical,
            "step7_wall_s": s7,
            "step7_speedup": s7_speedup,
            "step7_saved_pct": s7_saved_pct,
            "total_wall_s": {s: round(wall[s], 1) for s in STATES},
            "total_speedup": total_speedup,
            "total_saved_pct": total_saved_pct,
            "peak_rss_mb": {s: round(rss[s], 0) for s in STATES},
            "rss_overhead_pct": rss_overhead_pct,
            "completeness": {"counts": counts, "unchanged": completeness_ok},
            # Gate: byte-identity (hard) AND a Step-7 wall win.
            "gate": {"byte_identical": identical,
                     "step7_win": bool(s7_saved_pct and s7_saved_pct > 0),
                     "completeness_unchanged": completeness_ok},
        }
        results.append(rec)
        print(f"\n  [{bid}] byte-identical two_phase==fused = {identical}", flush=True)
        print(f"  STEP 7 wall: two_phase={s7['two_phase']}s fused={s7['fused']}s "
              f"(speedup {s7_speedup}× / saved {s7_saved_pct}%)", flush=True)
        print(f"  total wall: two_phase={rec['total_wall_s']['two_phase']}s "
              f"fused={rec['total_wall_s']['fused']}s "
              f"({total_speedup}× / {total_saved_pct}%) | RSS "
              f"{rec['peak_rss_mb']['two_phase']}→{rec['peak_rss_mb']['fused']}MB "
              f"({rss_overhead_pct:+.1f}%)", flush=True)
        print(f"  completeness unchanged = {completeness_ok}", flush=True)
        if not identical:
            print(f"  !! WARNING: outputs differ — byte-identity gate FAILED; "
                  f"the fusion is NOT byte-neutral, do not ship", flush=True)
        if s7_saved_pct is not None and s7_saved_pct <= 0:
            print(f"  !! NOTE: fused Step-7 wall did not beat two-phase on {bid} "
                  f"— rerun with LIFTON_FUSE_MAT_CONCURRENCY=4 to test capping "
                  f"the materialise half (contention canary).", flush=True)

    (HERE / "fuse_step7_ab.json").write_text(json.dumps(results, indent=2))
    _write_md(results, HERE / "fuse_step7_ab.md")
    print("\nDONE: wrote fuse_step7_ab.json + fuse_step7_ab.md", flush=True)
    return 0


def _write_md(results, path):
    lines = [
        "## Fused Step-7 pool vs two-phase A/B — Iteration 10\n",
        "Both states run on cached `-L`/`-M` (Step 4 skipped → Step-7 dispatch "
        "isolated, byte-identity clean) with the SAME argv (`-t N "
        "--locus-pipeline`), differing only by the `LIFTON_FUSE_STEP7` env gate. "
        "**two_phase** = `LIFTON_FUSE_STEP7=0` (pre-Iteration-10 prefetcher pool + "
        "barrier); **fused** = default (materialise + process in one pool). The "
        "change is pure scheduling. Gate: **byte-identical (hard) AND a Step-7 "
        "wall win**.\n",
        f"Threads: `-t {THREADS}`. **Step-7 wall** (probe `LIFTON_PERF_STEP7`) is "
        "the direct measure; **total wall** is end-to-end.\n",
        "| Dataset | byte-ident | Step-7 wall two_phase→fused | Step-7 speedup | total wall two_phase→fused | total speedup | peak RSS Δ | complete |",
        "|---|---|---|---|---|---|---|---|",
    ]
    for r in results:
        s7, w = r["step7_wall_s"], r["total_wall_s"]
        lines.append(
            "| {} | {} | {}s→{}s | {}× ({}%) | {}s→{}s | {}× ({}%) | {:+.1f}% | {} |".format(
                r["benchmark"], "yes" if r["byte_identical"] else "**NO**",
                s7["two_phase"], s7["fused"], r["step7_speedup"], r["step7_saved_pct"],
                w["two_phase"], w["fused"], r["total_speedup"], r["total_saved_pct"],
                r["rss_overhead_pct"],
                "yes" if r["completeness"]["unchanged"] else "**NO**"))
    lines.append(
        "\n**Interpretation:** byte-identical is the hard gate (the fusion is "
        "pure scheduling). A positive Step-7 speedup confirms the SQLite-bound "
        "materialise now overlaps the GIL-released parasail instead of running "
        "in a separate barrier'd phase. The rigorous byte-identity proof is the "
        "unit suite (`tests/test_fuse_step7.py` + the 24-cell matrix "
        "`tests/test_native_matrix.py`, whose on-disk cells take the fused path).")
    path.write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    sys.exit(main())
