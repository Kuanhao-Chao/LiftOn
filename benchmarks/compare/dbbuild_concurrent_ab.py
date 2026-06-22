#!/usr/bin/env python
"""A/B for the CONCURRENT Step-5 gffutils DB build (Iteration 11).

AUDIT RECORD (Iteration-11 NO-GO). The `lifton/` change this A/B drove
(threaded concurrent Step-5 DB build, gated by `LIFTON_CONCURRENT_DBBUILD` +
`LIFTON_PERF_STEP5`) was implemented, measured, and **REVERTED**. It was
byte-identical on all three datasets BUT consistently SLOWER at the Step-5
build itself (mouse 28.4->30.8s -8.4%, drosophila 11.9->13.2s -11.6%, rice
14.0->15.1s -8.3%) because `gffutils.create_db` is GIL-bound Python (line-by-
line parsing + Feature construction), so two builds on two THREADS time-slice
on one core with GIL contention instead of running in parallel. Threading
helps only GIL-RELEASING work (the miniprot SUBPROCESS in concurrent Step 4,
parasail in Step-7 fusion) — not a pure-Python gffutils build. True
parallelism would need a separate PROCESS (own GIL), but the ceiling is only
~2.5% total wall (mammalian-only) and it forks mid-pipeline. This script +
its committed `dbbuild_concurrent_ab.{json,md}` are the audit trail; re-running
requires restoring the instrumentation. See the CLAUDE.md Iteration-11 note.

Step 5 builds two INDEPENDENT gffutils DBs serially: `liftoff.gff3` →
`liftoff.gff3_db`, then `miniprot.gff3` → `miniprot.gff3_db`. Iteration 11
builds them concurrently — the miniprot DB on a worker thread, liftoff on the
main thread — then reopens the miniprot DB on the main thread by its on-disk
`.dbfn` (sidestepping SQLite thread-affinity). Wall collapses from
`t_liftoff_db + t_miniprot_db` toward `max(...)`, saving ≈ the smaller
(miniprot) build. The change is pure scheduling of *where* the deterministic
`.db` file is built → BYTE-NEUTRAL.

This A/B proves both halves:
  - **byte-identity** (the HARD gate): serial vs concurrent are identical;
  - **speed**: the concurrent Step-5 build is faster than serial.

Both states reuse cached `-L`/`-M` so Step 4 is a fast file-load and the
Step-5 DB build is a large fraction of total wall (measurable). `_clean_input_dbs`
removes the stale `_db` siblings so each state rebuilds the DBs from scratch.
Same argv for every state; only the `LIFTON_CONCURRENT_DBBUILD` env gate differs:

  state "serial"      LIFTON_CONCURRENT_DBBUILD=0   serial build (baseline)
  state "concurrent"  (default)                     concurrent build  ← the iteration

Step-5 wall is probed via `LIFTON_PERF_STEP5`. The win scales with GFF size, so
mouse (liftoff 255K / miniprot 120K lines) is the canary.

Usage (repo root, lifton_devel env):
    python -m benchmarks.compare.dbbuild_concurrent_ab [IDS...]
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

# Win scales with DB-build cost (∝ GFF size). mouse is the mammalian canary
# (liftoff 255K / miniprot 120K lines); drosophila/rice are ~92K/46K.
DEFAULT_IDS = ["mouse", "drosophila", "rice"]

STATES = {
    "serial":     {"LIFTON_CONCURRENT_DBBUILD": "0"},  # serial build baseline
    "concurrent": {},                                   # default = concurrent build
}


def _run(bid, state, env_extra, p, root):
    liftoff_gff = WORK / bid / "tools" / "liftoff" / "liftoff.gff3"
    miniprot_gff = WORK / bid / "tools" / "miniprot" / "miniprot.gff3"
    out_dir = root / state
    out_dir.mkdir(parents=True, exist_ok=True)
    out = out_dir / f"{state}.gff3"
    # Force a fresh _db rebuild per state (inputs unchanged → identical rebuild).
    _clean_input_dbs(p["ref_gff"], liftoff_gff, miniprot_gff)
    # SAME argv for every state; only the LIFTON_CONCURRENT_DBBUILD env differs.
    # No -t / --locus-pipeline: Step 7 is irrelevant here, so keep it serial to
    # isolate the Step-5 DB build and reduce noise.
    argv = [TOOLS["lifton_bin"], "-copies", "-ad", _ann_db(bid),
            "-g", p["ref_gff"], "-L", str(liftoff_gff), "-M", str(miniprot_gff),
            "-o", str(out), p["tgt_fa"], p["ref_fa"]]
    env = _compose_env(TOOLS)
    env["LIFTON_PERF_STEP5"] = "1"   # emit the Step-5 wall probe to stderr
    env.update(env_extra)
    pr = run_profiled(argv, label=f"dbbuild_{state}",
                      log_dir=root / "logs", env=env, log=print)
    if pr.exit_code != 0 or not out.exists() or out.stat().st_size == 0:
        raise RuntimeError(f"{bid}: state {state} failed (exit {pr.exit_code}); "
                           f"see {pr.stderr_path}")
    return out, pr


def _parse_step5_wall(stderr_path):
    """Pull the LIFTON_PERF_STEP5 probe line ('[LiftOn][perf] Step5 wall
    (mode): N.NNs') — the direct measure of the concurrent-build saving."""
    try:
        for ln in Path(stderr_path).read_text().splitlines():
            if "[LiftOn][perf] Step5 wall" in ln:
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
        print(f"=== {bid}: concurrent Step-5 DB build vs serial A/B "
              f"(cached -L/-M) ===", flush=True)
        man = json.loads((WORK / bid / "subset" / "subset.manifest.json").read_text())
        p = man["paths"]
        root = WORK / bid / "_dbbuild_concurrent_ab"
        root.mkdir(parents=True, exist_ok=True)

        outs, prs = {}, {}
        for state, env_extra in STATES.items():
            print(f"--- {bid} state {state} ---", flush=True)
            outs[state], prs[state] = _run(bid, state, env_extra, p, root)

        identical = outs["serial"].read_bytes() == outs["concurrent"].read_bytes()
        counts = {s: _feature_counts(outs[s]) for s in STATES}
        completeness_ok = counts["serial"] == counts["concurrent"]

        wall = {s: prs[s].wall_clock_seconds for s in STATES}
        s5 = {s: _parse_step5_wall(prs[s].stderr_path) for s in STATES}
        rss = {s: getattr(prs[s], "peak_rss_mb", 0) or 0 for s in STATES}

        s5_speedup = (round(s5["serial"] / max(s5["concurrent"], 1e-9), 2)
                      if (s5["serial"] and s5["concurrent"]) else None)
        s5_saved_pct = (round(100.0 * (s5["serial"] - s5["concurrent"]) /
                              max(s5["serial"], 1e-9), 1)
                        if (s5["serial"] and s5["concurrent"]) else None)
        total_speedup = round(wall["serial"] / max(wall["concurrent"], 1e-9), 2)
        total_saved_pct = round(100.0 * (wall["serial"] - wall["concurrent"]) /
                                max(wall["serial"], 1e-9), 1)
        rss_overhead_pct = round(100.0 * (rss["concurrent"] - rss["serial"]) /
                                 max(rss["serial"], 1e-9), 1)

        rec = {
            "benchmark": bid,
            "byte_identical": identical,
            "step5_wall_s": s5,
            "step5_speedup": s5_speedup,
            "step5_saved_pct": s5_saved_pct,
            "total_wall_s": {s: round(wall[s], 1) for s in STATES},
            "total_speedup": total_speedup,
            "total_saved_pct": total_saved_pct,
            "peak_rss_mb": {s: round(rss[s], 0) for s in STATES},
            "rss_overhead_pct": rss_overhead_pct,
            "completeness": {"counts": counts, "unchanged": completeness_ok},
            "gate": {"byte_identical": identical,
                     "step5_win": bool(s5_saved_pct and s5_saved_pct > 0),
                     "completeness_unchanged": completeness_ok},
        }
        results.append(rec)
        print(f"\n  [{bid}] byte-identical serial==concurrent = {identical}", flush=True)
        print(f"  STEP 5 wall: serial={s5['serial']}s concurrent={s5['concurrent']}s "
              f"(speedup {s5_speedup}× / saved {s5_saved_pct}%)", flush=True)
        print(f"  total wall: serial={rec['total_wall_s']['serial']}s "
              f"concurrent={rec['total_wall_s']['concurrent']}s "
              f"({total_speedup}× / {total_saved_pct}%) | RSS "
              f"{rec['peak_rss_mb']['serial']}→{rec['peak_rss_mb']['concurrent']}MB "
              f"({rss_overhead_pct:+.1f}%)", flush=True)
        print(f"  completeness unchanged = {completeness_ok}", flush=True)
        if not identical:
            print(f"  !! WARNING: outputs differ — byte-identity gate FAILED; "
                  f"the concurrent build is NOT byte-neutral, do not ship", flush=True)

    (HERE / "dbbuild_concurrent_ab.json").write_text(json.dumps(results, indent=2))
    _write_md(results, HERE / "dbbuild_concurrent_ab.md")
    print("\nDONE: wrote dbbuild_concurrent_ab.json + dbbuild_concurrent_ab.md",
          flush=True)
    return 0


def _write_md(results, path):
    lines = [
        "## Concurrent Step-5 DB build vs serial A/B — Iteration 11\n",
        "Both states reuse cached `-L`/`-M` (Step 4 a fast load → Step-5 DB "
        "build a large fraction of wall) with the SAME argv, differing only by "
        "the `LIFTON_CONCURRENT_DBBUILD` env gate. **serial** = "
        "`LIFTON_CONCURRENT_DBBUILD=0` (build liftoff then miniprot DB serially); "
        "**concurrent** = default (miniprot DB on a worker, liftoff on main, "
        "reopen on main). Pure scheduling of where the `.db` file is built. "
        "Gate: **byte-identical (hard) AND a Step-5 wall win**.\n",
        "**Step-5 wall** (probe `LIFTON_PERF_STEP5`) is the direct measure; "
        "**total wall** is end-to-end.\n",
        "| Dataset | byte-ident | Step-5 wall serial→concurrent | Step-5 speedup | total wall serial→concurrent | total speedup | peak RSS Δ | complete |",
        "|---|---|---|---|---|---|---|---|",
    ]
    for r in results:
        s5, w = r["step5_wall_s"], r["total_wall_s"]
        lines.append(
            "| {} | {} | {}s→{}s | {}× ({}%) | {}s→{}s | {}× ({}%) | {:+.1f}% | {} |".format(
                r["benchmark"], "yes" if r["byte_identical"] else "**NO**",
                s5["serial"], s5["concurrent"], r["step5_speedup"], r["step5_saved_pct"],
                w["serial"], w["concurrent"], r["total_speedup"], r["total_saved_pct"],
                r["rss_overhead_pct"],
                "yes" if r["completeness"]["unchanged"] else "**NO**"))
    lines.append(
        "\n**Interpretation:** byte-identical is the hard gate (the change is "
        "pure scheduling — the `.db` file contents are identical, only *where* "
        "it is built moves). A positive Step-5 speedup confirms the miniprot DB "
        "build now overlaps the liftoff DB build. The win scales with DB-build "
        "cost (∝ GFF size), so mouse (255K/120K lines) shows the largest saving. "
        "The rigorous byte-identity proof is the unit suite "
        "(`tests/test_dbbuild_concurrent.py` + the 24-cell matrix "
        "`tests/test_native_matrix.py`, whose cells all take the concurrent path).")
    path.write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    sys.exit(main())
