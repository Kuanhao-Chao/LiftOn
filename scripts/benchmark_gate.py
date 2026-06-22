#!/usr/bin/env python
"""Regression + benchmark gate for LiftOn algorithmic changes.

Runs two checks and exits non-zero if either fails:

  1. pytest, including the 24-cell byte-identity matrix
     (tests/test_native_matrix.py) and the hermetic integration pipeline.
     This guards that the default output path stays byte-frozen and that all
     fast-path configs remain mutually byte-identical.
  2. A FAST single-chromosome benchmark: re-runs ONLY the LiftOn step (reusing
     the cached Liftoff + miniprot outputs via -L/-M) plus the custom evaluator,
     then compares LiftOn's protein identity / completeness / wall-clock against
     a committed baseline. Any regression beyond the thresholds below fails.

This is the gate a `.git/hooks/pre-push` (see scripts/hooks/pre-push) and the
`make benchmark-gate` target invoke after every algorithmic change. The full
8-benchmark suite (`make benchmark-suite`) runs on demand / nightly.

Usage (under the lifton_devel env):
    python scripts/benchmark_gate.py                  # run gate, fail on regression
    python scripts/benchmark_gate.py --update-baseline # reseed baseline, exit 0
    python scripts/benchmark_gate.py --skip-pytest     # benchmark check only
    python scripts/benchmark_gate.py --benchmark rice  # use a different fast dataset
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
COMPARE = REPO / "benchmarks" / "compare"
WORK = COMPARE / "work"
BASELINE_DIR = COMPARE / "baseline"
PY = sys.executable

# Fastest cached benchmark and the MAIN one — regressions here matter most. The
# full all-feature-type coverage is validated by the on-demand 8-benchmark suite.
DEFAULT_FAST_ID = "human_mane"

# Accuracy metrics: max allowed ABSOLUTE drop vs baseline. Wall clock: max
# allowed RELATIVE increase. Tune here if the fast benchmark proves noisy.
THRESHOLDS = {
    "protein_identity_mean": 0.005,
    "completeness_coding": 0.01,
    "completeness_feature_total": 0.01,
    "wall_clock_regress_frac": 0.25,
}

PYTEST_TARGETS = ["tests/test_native_matrix.py", "tests/test_integration_pipeline.py"]

# repo root must be importable for `benchmarks.compare.*`
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))


def _run(cmd):
    print("+ " + " ".join(str(c) for c in cmd), flush=True)
    return subprocess.run(cmd, cwd=str(REPO))


def run_pytest_gate() -> bool:
    r = _run([PY, "-m", "pytest", *PYTEST_TARGETS, "-q", "-p", "no:cacheprovider"])
    return r.returncode == 0


def run_fast_benchmark(fast_id: str) -> bool:
    """Refresh ONLY the LiftOn step + eval (cached subset/Liftoff/miniprot are
    reused) so the gate measures the current algorithm without paying for the
    external aligners."""
    for sentinel in ("lifton.done", "eval.done"):
        p = WORK / fast_id / ".done" / sentinel
        if p.exists():
            p.unlink()
    from benchmarks.compare import run_compare
    rc = run_compare.main(["--benchmarks", fast_id, "--feature-modes", "genes",
                           "--no-liftofftools", "--no-report"])
    return rc == 0


def load_metrics(fast_id: str) -> dict:
    sm = json.loads((WORK / fast_id / "eval" / "summaries.json").read_text())
    s = sm["lifton"]
    return {
        "benchmark": fast_id,
        "protein_identity_mean": (s.get("protein_identity") or {}).get("mean"),
        "completeness_coding": s.get("completeness_coding"),
        "completeness_feature_total": s.get("completeness_feature_total"),
        "wall_clock_seconds": (s.get("profile") or {}).get("wall_clock_seconds"),
    }


def baseline_path(fast_id: str) -> Path:
    return BASELINE_DIR / f"{fast_id}.baseline.json"


def check(cur: dict, base: dict):
    rows = []
    ok = True
    for key in ("protein_identity_mean", "completeness_coding",
                "completeness_feature_total"):
        b, c, thr = base.get(key), cur.get(key), THRESHOLDS[key]
        if b is None or c is None:
            rows.append((key, b, c, "—", "skip (n/a)"))
            continue
        drop = b - c
        passed = drop <= thr
        ok = ok and passed
        rows.append((key, f"{b:.5f}", f"{c:.5f}", f"{c - b:+.5f}",
                     "PASS" if passed else f"FAIL (drop {drop:.5f} > {thr})"))
    b, c = base.get("wall_clock_seconds"), cur.get("wall_clock_seconds")
    thr = THRESHOLDS["wall_clock_regress_frac"]
    if b and c:
        frac = (c - b) / b
        passed = frac <= thr
        ok = ok and passed
        rows.append(("wall_clock_seconds", f"{b:.1f}", f"{c:.1f}", f"{frac * 100:+.0f}%",
                     "PASS" if passed else f"FAIL (+{frac * 100:.0f}% > {thr * 100:.0f}%)"))
    else:
        rows.append(("wall_clock_seconds", b, c, "—", "skip (n/a)"))
    return ok, rows


def _print_table(rows):
    w = max(len(r[0]) for r in rows)
    print(f"\n  {'metric'.ljust(w)}  {'baseline':>10}  {'current':>10}  {'delta':>9}  verdict")
    print("  " + "-" * (w + 46))
    for name, b, c, d, v in rows:
        print(f"  {str(name).ljust(w)}  {str(b):>10}  {str(c):>10}  {str(d):>9}  {v}")


def main(argv=None):
    ap = argparse.ArgumentParser(description="LiftOn regression + benchmark gate")
    ap.add_argument("--update-baseline", action="store_true",
                    help="reseed the baseline from a fresh run and exit 0")
    ap.add_argument("--skip-pytest", action="store_true",
                    help="skip the pytest/byte-identity gate (benchmark check only)")
    ap.add_argument("--benchmark", default=DEFAULT_FAST_ID,
                    help=f"fast benchmark id (default {DEFAULT_FAST_ID})")
    a = ap.parse_args(argv)
    fast_id = a.benchmark

    if not a.update_baseline and not a.skip_pytest:
        print("=== [1/2] pytest + 24-cell byte-identity gate ===")
        if not run_pytest_gate():
            print("\nGATE FAIL: pytest (incl. 24-cell byte-identity) failed.")
            return 1

    print(f"\n=== {'reseed' if a.update_baseline else '[2/2]'} fast benchmark: {fast_id} ===")
    if not run_fast_benchmark(fast_id):
        print("\nGATE FAIL: fast benchmark run failed.")
        return 1
    cur = load_metrics(fast_id)

    if a.update_baseline:
        BASELINE_DIR.mkdir(parents=True, exist_ok=True)
        baseline_path(fast_id).write_text(json.dumps(cur, indent=2) + "\n")
        print(f"\nbaseline updated: {baseline_path(fast_id)}")
        print(json.dumps(cur, indent=2))
        return 0

    bp = baseline_path(fast_id)
    if not bp.exists():
        print(f"\nGATE: no baseline at {bp}.\n"
              f"Seed it with: python scripts/benchmark_gate.py --update-baseline"
              f"{'' if fast_id == DEFAULT_FAST_ID else f' --benchmark {fast_id}'}")
        return 1
    base = json.loads(bp.read_text())
    ok, rows = check(cur, base)
    _print_table(rows)
    print(f"\nGATE {'PASS' if ok else 'FAIL'} ({fast_id})")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
