#!/usr/bin/env python
"""A/B for fresh parallel Step 7 WITHOUT --native (Iteration 8).

Before Iteration 8, ``--threads N --locus-pipeline`` only fanned Step 7 out
across worker threads when ``--native`` (or ``LIFTON_PARALLEL_FORCE``) was
set; on a plain default (gffutils SQLite) run it silently downgraded to
serial. Iteration 8 routes the parallel path through the materialise +
proxy-DB machinery unconditionally, so it parallelises on any backend
without ``--native`` — and stays byte-identical.

This A/B proves both halves of the claim on real data:
  - **byte-identity** (the HARD gate): serial vs parallel-no-native vs
    parallel-native are all identical;
  - **speed**: parallel-no-native is materially faster at Step 7 than
    serial, and ≈ parallel-native (the win is now available WITHOUT
    --native).

Unlike `parallel_aligners_ab.py` (which must run Step 4 fresh), this A/B
**reuses cached `-L`/`-M`** so Step 4 is skipped. That (a) isolates the
Step-7 dispatch — the only thing this iteration changes — and (b) makes
the byte-identity comparison clean (Liftoff's `-copies` is non-
deterministic across FRESH runs; reusing fixed `-L`/`-M` removes that
noise). All three states run at the same ``-t`` so only the Step-7 fan-out
differs.

  state "serial"           -t N                          Step 7 serial
  state "parallel"         -t N --locus-pipeline         Step 7 parallel, NO --native  ← the iteration
  state "parallel_native"  -t N --locus-pipeline --native Step 7 parallel + --native (parity check)

Step-7 wall is probed via ``LIFTON_PERF_STEP7``; total wall is end-to-end
(dominated by Step 7 here since Step 4 is a cached load).

Usage (repo root, lifton_devel env):
    python -m benchmarks.compare.fresh_parallel_step7_ab [IDS...]
    LIFTON_AB_THREADS=4 python -m benchmarks.compare.fresh_parallel_step7_ab drosophila
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

# Step-7-dominated subsets (enough loci for fan-out to pay off). bee is
# omitted by default — at bee scale Step 7 is a small fraction.
DEFAULT_IDS = ["drosophila", "mouse_to_rat", "rice"]
THREADS = os.environ.get("LIFTON_AB_THREADS", "8")

STATES = {
    "serial":          [],                              # -t N, no fan-out
    "parallel":        ["--locus-pipeline"],            # the Iteration-8 path (NO --native)
    "parallel_native": ["--locus-pipeline", "--native"],  # parity check
}


def _run(bid, state, flag, p, root):
    liftoff_gff = WORK / bid / "tools" / "liftoff" / "liftoff.gff3"
    miniprot_gff = WORK / bid / "tools" / "miniprot" / "miniprot.gff3"
    out_dir = root / state
    out_dir.mkdir(parents=True, exist_ok=True)
    out = out_dir / f"{state}.gff3"
    # Rebuild the input gffutils sibling DBs fresh for each state so the
    # run is deterministic (inputs unchanged → identical rebuild).
    _clean_input_dbs(p["ref_gff"], liftoff_gff, miniprot_gff)
    # Same -t for every state: only the Step-7 fan-out (--locus-pipeline /
    # --native) differs, so the cached Step-4 load is a constant.
    argv = [TOOLS["lifton_bin"], "-t", THREADS, "-copies", "-ad", _ann_db(bid),
            "-g", p["ref_gff"], "-L", str(liftoff_gff), "-M", str(miniprot_gff),
            "-o", str(out), *flag, p["tgt_fa"], p["ref_fa"]]
    env = _compose_env(TOOLS)
    env["LIFTON_PERF_STEP7"] = "1"   # emit the Step-7 wall probe line to stderr
    pr = run_profiled(argv, label=f"fresh_p7_{state}",
                      log_dir=root / "logs", env=env, log=print)
    if pr.exit_code != 0 or not out.exists() or out.stat().st_size == 0:
        raise RuntimeError(f"{bid}: state {state} failed (exit {pr.exit_code}); "
                           f"see {pr.stderr_path}")
    return out, pr


def _parse_step7_wall(stderr_path):
    """Pull the LIFTON_PERF_STEP7 probe line ('[LiftOn][perf] Step7 wall
    (mode): N.NNs') — the Step-4-free measure of the fan-out speedup."""
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
        print(f"=== {bid}: fresh parallel Step 7 (no --native) A/B "
              f"(-t {THREADS}, cached -L/-M) ===", flush=True)
        man = json.loads((WORK / bid / "subset" / "subset.manifest.json").read_text())
        p = man["paths"]
        root = WORK / bid / "_fresh_parallel_step7_ab"
        root.mkdir(parents=True, exist_ok=True)

        outs, prs = {}, {}
        for state, flag in STATES.items():
            print(f"--- {bid} state {state} ---", flush=True)
            outs[state], prs[state] = _run(bid, state, flag, p, root)

        b_serial = outs["serial"].read_bytes()
        b_par = outs["parallel"].read_bytes()
        b_par_n = outs["parallel_native"].read_bytes()
        ident_par = (b_serial == b_par)               # the iteration's gate
        ident_native_parity = (b_par == b_par_n)      # no-native == native
        all_identical = ident_par and ident_native_parity

        counts = {s: _feature_counts(outs[s]) for s in STATES}
        completeness_ok = (counts["serial"] == counts["parallel"] ==
                           counts["parallel_native"])

        wall = {s: prs[s].wall_clock_seconds for s in STATES}
        s7 = {s: _parse_step7_wall(prs[s].stderr_path) for s in STATES}

        s7_speedup = (round(s7["serial"] / max(s7["parallel"], 1e-9), 2)
                      if (s7["serial"] and s7["parallel"]) else None)
        s7_saved_pct = (round(100.0 * (s7["serial"] - s7["parallel"]) /
                              max(s7["serial"], 1e-9), 1)
                        if (s7["serial"] and s7["parallel"]) else None)
        total_speedup = round(wall["serial"] / max(wall["parallel"], 1e-9), 2)
        total_saved_pct = round(100.0 * (wall["serial"] - wall["parallel"]) /
                                max(wall["serial"], 1e-9), 1)

        rec = {
            "benchmark": bid,
            "threads": THREADS,
            "byte_identical_serial_vs_parallel": ident_par,
            "byte_identical_native_parity": ident_native_parity,
            "step7_wall_s": s7,
            "step7_speedup": s7_speedup,
            "step7_saved_pct": s7_saved_pct,
            "total_wall_s": {s: round(wall[s], 1) for s in STATES},
            "total_speedup": total_speedup,
            "total_saved_pct": total_saved_pct,
            "completeness": {"counts": counts, "unchanged": completeness_ok},
            # Gate: byte-identity (hard) AND a Step-7 wall win.
            "gate": {"byte_identical": all_identical,
                     "step7_win": bool(s7_saved_pct and s7_saved_pct > 0),
                     "completeness_unchanged": completeness_ok},
        }
        results.append(rec)
        print(f"\n  [{bid}] byte-identical serial==parallel(no-native)={ident_par}  "
              f"parallel==native={ident_native_parity}", flush=True)
        print(f"  STEP 7 wall: serial={s7['serial']}s "
              f"parallel(no-native)={s7['parallel']}s "
              f"native={s7['parallel_native']}s "
              f"(speedup {s7_speedup}× / saved {s7_saved_pct}%)", flush=True)
        print(f"  total wall: serial={rec['total_wall_s']['serial']}s "
              f"parallel={rec['total_wall_s']['parallel']}s "
              f"({total_speedup}× / {total_saved_pct}%)", flush=True)
        print(f"  completeness unchanged={completeness_ok}", flush=True)
        if not all_identical:
            print(f"  !! WARNING: outputs differ — byte-identity gate FAILED; "
                  f"investigate before trusting the speed numbers", flush=True)

    (HERE / "fresh_parallel_step7_ab.json").write_text(json.dumps(results, indent=2))
    _write_md(results, HERE / "fresh_parallel_step7_ab.md")
    print("\nDONE: wrote fresh_parallel_step7_ab.json + fresh_parallel_step7_ab.md",
          flush=True)
    return 0


def _write_md(results, path):
    lines = [
        "## Fresh parallel Step 7 without `--native` — Iteration 8\n",
        "Each state reuses cached `-L`/`-M` (Step 4 skipped → Step-7 dispatch "
        "isolated, byte-identity clean) and runs at the same `-t`, differing "
        "only by `--locus-pipeline` / `--native`. The change routes parallel "
        "Step 7 through the materialise + proxy-DB path on **any** backend, so "
        "`--locus-pipeline` now fans out WITHOUT `--native`. Gate: **byte-"
        "identical (hard) AND a Step-7 wall win**.\n",
        f"Threads: `-t {THREADS}`. **Step-7 wall** (probe `LIFTON_PERF_STEP7`) "
        "is the direct measure; **total wall** is end-to-end.\n",
        "| Dataset | byte-ident serial=parallel | native parity | Step-7 wall serial→parallel | Step-7 speedup | total wall serial→parallel | total speedup | complete |",
        "|---|---|---|---|---|---|---|---|",
    ]
    for r in results:
        s7, w = r["step7_wall_s"], r["total_wall_s"]
        lines.append(
            "| {} | {} | {} | {}s→{}s | {}× ({}%) | {}s→{}s | {}× ({}%) | {} |".format(
                r["benchmark"],
                "yes" if r["byte_identical_serial_vs_parallel"] else "**NO**",
                "yes" if r["byte_identical_native_parity"] else "**NO**",
                s7["serial"], s7["parallel"], r["step7_speedup"], r["step7_saved_pct"],
                w["serial"], w["parallel"], r["total_speedup"], r["total_saved_pct"],
                "yes" if r["completeness"]["unchanged"] else "**NO**"))
    lines.append(
        "\n**Interpretation:** byte-identical across all three states is the "
        "hard gate (the change is pure scheduling). `parallel ≈ "
        "parallel_native` Step-7 wall confirms the speedup no longer requires "
        "`--native`. The rigorous byte-identity proof is the unit suite "
        "(`tests/test_fresh_parallel_step7.py` + the now-non-native-exercising "
        "24-cell matrix `tests/test_native_matrix.py`).")
    path.write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    sys.exit(main())
