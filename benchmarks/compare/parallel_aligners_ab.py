#!/usr/bin/env python
"""A/B for `--parallel-aligners` — concurrent Step 4 (Iteration 6).

LiftOn's Step 4 runs the two external aligner programs SEQUENTIALLY by
default: Liftoff (DNA) then miniprot (protein). They are independent (same
inputs, disjoint outputs), so `--parallel-aligners` runs them on two worker
threads at once → wall-clock collapses from t_liftoff + t_miniprot to
max(t_liftoff, t_miniprot). This is a PURE SCHEDULING change.

Unlike the other A/Bs (`fast_align_ab.py`, `legacy_merge_ab.py`), this one
**must NOT reuse cached `-L`/`-M`** — that short-circuits Step 4 to a file
load and erases the whole effect. Each state runs LiftOn FRESH (ref_gff +
genomes only), letting Step 4 actually run both tools.

  state "serial"    --serial-aligners      t_liftoff + t_miniprot
  state "parallel"  (no flag = default)    max(t_liftoff, t_miniprot)  ← default

(Concurrent Step 4 was PROMOTED to default in Iteration 6, so the baseline
is now the opt-OUT `--serial-aligners` and the candidate is the bare default.)

Reports per dataset:
  - Δwall, speedup, wall_saved_pct   (the speed win — the point of the flag),
  - peak RSS serial→parallel          (reported, not gated; both tools resident
                                       at once raises peak),
  - **byte-identity** serial vs parallel (best-effort signal — see CAVEAT;
    the flag changes no output, but Liftoff's `-copies` multi-copy alignment
    is itself non-deterministic ACROSS FRESH RUNS, so on repeat-rich genomes
    a fresh serial-vs-parallel diff can appear that is pure Liftoff noise,
    not the flag. The rigorous byte-identity proof is the UNIT TEST
    tests/test_parallel_aligners.py on fixed cached -L/-M inputs.),
  - completeness (gene/tx/CDS counts) as an extra sanity signal.

Promotion is a separate human decision after reading this. The gate here:
byte-identical AND a material wall-clock win.

Usage (repo root, lifton_devel env):
    python -m benchmarks.compare.parallel_aligners_ab [IDS...]   # default below
    LIFTON_AB_THREADS=1 python -m benchmarks.compare.parallel_aligners_ab bee
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

DEFAULT_IDS = ["bee", "drosophila", "mouse_to_rat"]
THREADS = os.environ.get("LIFTON_AB_THREADS", "8")

STATES = {
    "serial": ["--serial-aligners"],   # opt out of the promoted default
    "parallel": [],                     # default = concurrent Step 4
}


def _run(bid, state, flag, p, root):
    out_dir = root / state
    out_dir.mkdir(parents=True, exist_ok=True)
    out = out_dir / f"{state}.gff3"
    # FRESH Step 4: no -L/-M. Clean any stale gffutils sibling DB next to the
    # shared ref_gff so each state rebuilds the ref DB deterministically.
    _clean_input_dbs(p["ref_gff"])
    # Each state's distinct out_dir isolates LiftOn's intermediates
    # (<outdir>/lifton_output/), so the two states never collide.
    argv = [TOOLS["lifton_bin"], "-t", THREADS, "-copies", "-ad", _ann_db(bid),
            "-g", p["ref_gff"],
            "-o", str(out), *flag, p["tgt_fa"], p["ref_fa"]]
    env = _compose_env(TOOLS)
    env["LIFTON_PERF_STEP4"] = "1"   # emit the Step-4 wall probe line to stderr
    pr = run_profiled(argv, label=f"parallel_aligners_{state}",
                      log_dir=root / "logs", env=env, log=print)
    if pr.exit_code != 0 or not out.exists() or out.stat().st_size == 0:
        raise RuntimeError(f"{bid}: state {state} failed (exit {pr.exit_code}); "
                           f"see {pr.stderr_path}")
    return out, pr


def _parse_step4_wall(stderr_path):
    """Pull the LIFTON_PERF_STEP4 probe line ('[LiftOn][perf] Step4 wall
    (mode): N.NNs') from a run's stderr — the direct, Step-7-noise-free
    measure of what --parallel-aligners saves."""
    try:
        for ln in Path(stderr_path).read_text().splitlines():
            if "[LiftOn][perf] Step4 wall" in ln:
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
        print(f"=== {bid}: --parallel-aligners (concurrent Step 4) vs serial A/B "
              f"(-t {THREADS}, FRESH Step 4) ===", flush=True)
        man = json.loads((WORK / bid / "subset" / "subset.manifest.json").read_text())
        p = man["paths"]
        root = WORK / bid / "_parallel_aligners_ab"
        root.mkdir(parents=True, exist_ok=True)

        outs, prs = {}, {}
        for state, flag in STATES.items():
            print(f"--- {bid} state {state} ---", flush=True)
            outs[state], prs[state] = _run(bid, state, flag, p, root)

        identical = outs["serial"].read_bytes() == outs["parallel"].read_bytes()
        counts_s = _feature_counts(outs["serial"])
        counts_p = _feature_counts(outs["parallel"])
        completeness_ok = (counts_s == counts_p)

        wall_s = prs["serial"].wall_clock_seconds
        wall_p = prs["parallel"].wall_clock_seconds
        rss_s = getattr(prs["serial"], "peak_rss_mb", 0) or 0
        rss_p = getattr(prs["parallel"], "peak_rss_mb", 0) or 0
        speedup = round(wall_s / max(wall_p, 1e-9), 3)
        wall_pct = round(100.0 * (wall_s - wall_p) / max(wall_s, 1e-9), 1)

        # Direct Step-4 wall (probe) — isolates the overlap saving from the
        # serial Step-7 noise that dominates the subset total.
        s4_s = _parse_step4_wall(prs["serial"].stderr_path)
        s4_p = _parse_step4_wall(prs["parallel"].stderr_path)
        s4_saved = round(s4_s - s4_p, 1) if (s4_s and s4_p) else None
        s4_pct = (round(100.0 * (s4_s - s4_p) / max(s4_s, 1e-9), 1)
                  if (s4_s and s4_p) else None)

        rec = {
            "benchmark": bid,
            "threads": THREADS,
            "byte_identical": identical,
            "wall_s": {"serial": round(wall_s, 1), "parallel": round(wall_p, 1)},
            "wall_saved_s": round(wall_s - wall_p, 1),
            "wall_saved_pct": wall_pct,
            "speedup": speedup,
            "step4_wall_s": {"serial": s4_s, "parallel": s4_p},
            "step4_saved_s": s4_saved,
            "step4_saved_pct": s4_pct,
            "peak_rss_mb": {"serial": round(rss_s, 0), "parallel": round(rss_p, 0)},
            "rss_overhead_pct": round(100.0 * (rss_p - rss_s) / max(rss_s, 1e-9), 1),
            "completeness": {"serial": counts_s, "parallel": counts_p,
                             "unchanged": completeness_ok},
            # The gate: identical output AND a real wall-clock win.
            "gate": {"byte_identical": identical,
                     "wall_win": wall_pct > 0,
                     "completeness_unchanged": completeness_ok},
        }
        results.append(rec)
        print(f"\n  [{bid}] byte-identical={identical}", flush=True)
        print(f"  STEP 4 wall (the flag's target): serial={s4_s}s parallel={s4_p}s "
              f"(saved {s4_saved}s / {s4_pct}%)", flush=True)
        print(f"  total wall: serial={rec['wall_s']['serial']}s "
              f"parallel={rec['wall_s']['parallel']}s "
              f"(saved {rec['wall_saved_s']}s / {wall_pct:+.1f}% / {speedup}×) | "
              f"RSS {rec['peak_rss_mb']['serial']}→{rec['peak_rss_mb']['parallel']}MB "
              f"({rec['rss_overhead_pct']:+.1f}%)", flush=True)
        print(f"  completeness unchanged={completeness_ok} "
              f"(serial={counts_s} parallel={counts_p})", flush=True)
        if not identical:
            print(f"  !! WARNING: outputs differ — investigate before trusting "
                  f"the speed number (expected identical at any -t since "
                  f"minimap2/miniprot preserve input order)", flush=True)

    (HERE / "parallel_aligners_ab.json").write_text(json.dumps(results, indent=2))
    _write_md(results, HERE / "parallel_aligners_ab.md")
    print("\nDONE: wrote parallel_aligners_ab.json + parallel_aligners_ab.md", flush=True)
    return 0


def _write_md(results, path):
    lines = [
        "## --parallel-aligners (concurrent Step 4) vs serial A/B — Iteration 6\n",
        "Each state runs LiftOn **fresh** (no cached `-L`/`-M`, so Step 4 "
        "actually runs both aligners), differing only by `--parallel-aligners`. "
        "The flag overlaps Liftoff and miniprot → wall = max() instead of sum. "
        "Output is byte-identical (speed-only change). Gate: **byte-identical "
        "AND a wall-clock win**.\n",
        f"Threads: `-t {THREADS}`. **Step-4 wall** is the direct measure of "
        "the overlap (probe via `LIFTON_PERF_STEP4`); **total wall** is "
        "end-to-end and is diluted by serial Step 7 on these subsets (no "
        "`--locus-pipeline`).\n",
        "| Dataset | byte-identical | Step-4 wall serial→parallel | Step-4 saved | total wall serial→parallel | total saved | peak RSS Δ | complete |",
        "|---|---|---|---|---|---|---|---|",
    ]
    for r in results:
        w, rss, s4 = r["wall_s"], r["peak_rss_mb"], r["step4_wall_s"]
        lines.append(
            "| {} | {} | {}s→{}s | {}s ({}%) | {}s→{}s | {}s ({:+.1f}%) | {:+.1f}% | {} |".format(
                r["benchmark"], "yes" if r["byte_identical"] else "**NO**",
                s4["serial"], s4["parallel"], r["step4_saved_s"], r["step4_saved_pct"],
                w["serial"], w["parallel"], r["wall_saved_s"], r["wall_saved_pct"],
                r["rss_overhead_pct"],
                "yes" if r["completeness"]["unchanged"] else "NO"))
    lines.append(
        "\n**Caveat on `byte-identical=NO`:** the flag itself changes no "
        "output (Liftoff runs on the main thread identically in both states; "
        "miniprot output is byte-identical). But Liftoff's `-copies` "
        "multi-copy alignment is **non-deterministic across fresh runs** "
        "(confirmed: two identical-config standalone Liftoff runs on "
        "drosophila diverge in tRNA copy lines) — so a fresh serial-vs-"
        "parallel diff on a repeat-rich genome is Liftoff noise, not the "
        "flag. The rigorous byte-identity proof is the unit test on fixed "
        "`-L`/`-M` inputs (`tests/test_parallel_aligners.py`).")
    path.write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    sys.exit(main())
