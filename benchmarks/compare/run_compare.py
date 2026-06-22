"""Driver for the three-way Liftoff/miniprot/LiftOn benchmark comparison.

Usage (from the repo root, under the lifton_devel env):
    python -m benchmarks.compare.run_compare --benchmarks human_mane
    python -m benchmarks.compare.run_compare --all -t 8
    python -m benchmarks.compare.run_compare --all --only-eval     # re-run eval+report only
    python -m benchmarks.compare.run_compare --benchmarks bee --force

Stages per benchmark (each resumable via work/<id>/.done/<stage>.done):
    subset → tools (liftoff, miniprot, lifton) → eval → liftofftools → report
"""
from __future__ import annotations

import argparse
import json
import sys
import traceback
from pathlib import Path

from . import subset_builder, tool_runners, evaluator, liftofftools_wrap, reporter

HERE = Path(__file__).resolve().parent
REGISTRY = HERE / "benchmarks.json"
WORK_ROOT = HERE / "work"
PLANS_DIR = HERE.parents[1] / "plans"


def _log(msg=""):
    print(msg, flush=True)


def load_registry() -> dict:
    return json.loads(REGISTRY.read_text())


def _load_cached_profiles(wd: Path, mode: str) -> dict:
    """Load tool profiles for a feature mode without re-running. liftoff/lifton
    come from the mode-specific tools dir; miniprot is always the shared genes
    profile (it is mode-independent)."""
    tools_sub = "tools" if mode == "genes" else f"tools_{mode}"
    profiles = {}
    for t in ("liftoff", "lifton"):
        pp = wd / tools_sub / "logs" / f"{t}.profile.json"
        if pp.exists():
            profiles[t] = json.loads(pp.read_text())
    mp = wd / "tools" / "logs" / "miniprot.profile.json"
    if mp.exists():
        profiles["miniprot"] = json.loads(mp.read_text())
    # lifton_optimize is genes-mode only
    if mode == "genes":
        op = wd / "tools" / "logs" / "lifton_optimize.profile.json"
        if op.exists():
            profiles["lifton_optimize"] = json.loads(op.read_text())
    return profiles


def run_benchmark(bench: dict, tools: dict, threads: int, stages: set, force: bool,
                  log=_log, modes=("genes",)) -> dict:
    bid = bench["id"]
    wd = WORK_ROOT / bid
    log(f"\n=== Benchmark: {bid} ({bench['species']}) ===")
    result = {"id": bid, "ok": False}
    manifest = None
    try:
        # subset is mode-independent: build/load it once and share across modes
        if "subset" in stages:
            manifest = subset_builder.build_subset(bench, wd, tools, threads, force, log)
        else:
            manifest = json.loads((wd / "subset" / "subset.manifest.json").read_text())
        for mode in modes:   # "genes" first so allfeat can reuse the genes miniprot
            log(f"  --- feature mode: {mode} ---")
            if "tools" in stages:
                profiles = tool_runners.run_all_tools(bench, manifest, wd, tools,
                                                      threads, force, log, mode=mode)
            else:
                profiles = _load_cached_profiles(wd, mode)
            if "eval" in stages:
                evaluator.evaluate_all(manifest, wd, profiles, force, log, mode=mode,
                                       threads=threads)
        if "liftofftools" in stages:   # genes-only cross-check (variants scores coding transcripts)
            liftofftools_wrap.run_crosscheck(manifest, wd, tools, force, log)
        result["ok"] = True
    except Exception as exc:
        log(f"  !! {bid} FAILED: {exc}")
        log(traceback.format_exc())
        result["error"] = str(exc)
    return result


def _run_selected(selected, tools, threads, stages, force, jobs, modes=("genes",)):
    """Run benchmarks sequentially (jobs==1) or concurrently (jobs>1).

    Benchmarks are fully independent (separate work/<id>/ dirs + .done
    sentinels), so concurrency is safe. In parallel mode each benchmark logs to
    its own work/<id>/run.log to keep output readable; a one-line start/finish
    is printed to the console.
    """
    if jobs <= 1 or len(selected) <= 1:
        return [run_benchmark(b, tools, threads, stages, force, modes=modes) for b in selected]

    from concurrent.futures import ThreadPoolExecutor, as_completed

    def _task(b):
        bid = b["id"]
        logf = WORK_ROOT / bid / "run.log"
        logf.parent.mkdir(parents=True, exist_ok=True)
        fh = open(logf, "w")

        def _blog(msg=""):
            fh.write(str(msg) + "\n")
            fh.flush()
        _log(f"[start] {bid}")
        try:
            r = run_benchmark(b, tools, threads, stages, force, log=_blog, modes=modes)
        finally:
            fh.close()
        _log(f"[{'done' if r['ok'] else 'FAIL'}] {bid}"
             + ("" if r["ok"] else f" — {r.get('error','')}") + f"  (log: {logf})")
        return r

    results = []
    with ThreadPoolExecutor(max_workers=jobs) as ex:
        futs = {ex.submit(_task, b): b["id"] for b in selected}
        for fut in as_completed(futs):
            results.append(fut.result())
    return results


def main(argv=None):
    reg = load_registry()
    ids = [b["id"] for b in reg["benchmarks"]]
    ap = argparse.ArgumentParser(description="Liftoff/miniprot/LiftOn benchmark comparison")
    ap.add_argument("--benchmarks", nargs="+", choices=ids, help="benchmark id(s) to run")
    ap.add_argument("--all", action="store_true", help="run all benchmarks")
    ap.add_argument("-t", "--threads", type=int, default=8)
    ap.add_argument("-j", "--jobs", type=int, default=1,
                    help="number of benchmarks to run concurrently (each uses --threads "
                         "threads; benchmarks are independent so this is safe). Default 1.")
    ap.add_argument("--force", action="store_true", help="ignore .done sentinels")
    ap.add_argument("--only-subset", action="store_true")
    ap.add_argument("--only-tools", action="store_true")
    ap.add_argument("--only-eval", action="store_true", help="eval + liftofftools + report only")
    ap.add_argument("--no-report", action="store_true")
    ap.add_argument("--no-liftofftools", action="store_true")
    ap.add_argument("--feature-modes", nargs="+", choices=["genes", "allfeat"],
                    default=["genes", "allfeat"],
                    help="feature mode set(s) to build/maintain: genes (Liftoff lifts only the "
                         "gene hierarchy) and/or allfeat (Liftoff lifts all top-level annotation "
                         "types). Default: both.")
    args = ap.parse_args(argv)
    # always run genes before allfeat (allfeat reuses the genes-mode miniprot)
    modes = tuple(m for m in ("genes", "allfeat") if m in args.feature_modes)

    if args.all:
        selected = reg["benchmarks"]
    elif args.benchmarks:
        selected = [b for b in reg["benchmarks"] if b["id"] in args.benchmarks]
    else:
        ap.error("specify --all or --benchmarks ID [ID ...]")

    all_stages = {"subset", "tools", "eval", "liftofftools"}
    if args.only_subset:
        stages = {"subset"}
    elif args.only_tools:
        stages = {"subset", "tools"}
    elif args.only_eval:
        stages = {"eval", "liftofftools"}
    else:
        stages = set(all_stages)
    if args.no_liftofftools:
        stages.discard("liftofftools")

    _log(f"Running {len(selected)} benchmark(s) with jobs={args.jobs}, threads={args.threads} "
         f"each, feature modes={list(modes)}")
    results = _run_selected(selected, reg["tools"], args.threads, stages, args.force, args.jobs,
                            modes=modes)

    # report aggregates ALL benchmarks that have eval output (not just this
    # run's selection), so re-running a single benchmark still yields a complete
    # cross-benchmark report.
    done_ids = [b["id"] for b in reg["benchmarks"]
                if (WORK_ROOT / b["id"] / "eval" / "summaries.json").exists()]
    if done_ids and not args.no_report and not args.only_subset and not args.only_tools:
        _log("\n=== Building comparison report ===")
        data = reporter.load(done_ids, WORK_ROOT, modes=modes)
        PLANS_DIR.mkdir(parents=True, exist_ok=True)
        flat = reporter.build_master_json_tsv(data, PLANS_DIR / "benchmark_plots")
        # also drop a flat tsv at the plans root
        (PLANS_DIR / "benchmark_comparison.tsv").write_text(
            (PLANS_DIR / "benchmark_plots" / "comparison.tsv").read_text()
            if (PLANS_DIR / "benchmark_plots" / "comparison.tsv").exists() else "")
        main_id = "human_mane" if "human_mane" in data else (done_ids[0] if done_ids else None)
        md = reporter.write_report(data, PLANS_DIR / "benchmark_comparison_report.md",
                                   main_id=main_id, log=_log)
        reporter.render_html(md, PLANS_DIR / "benchmark_comparison_report.html", _log)
        _log(f"\nReport: {PLANS_DIR / 'benchmark_comparison_report.html'}")

    n_ok = sum(1 for r in results if r["ok"])
    _log(f"\n=== Done: {n_ok}/{len(results)} benchmarks ok ===")
    for r in results:
        _log(f"  {r['id']}: {'OK' if r['ok'] else 'FAILED — ' + r.get('error','')}")
    return 0 if n_ok == len(results) else 1


if __name__ == "__main__":
    sys.exit(main())
