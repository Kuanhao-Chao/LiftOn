#!/usr/bin/env bash
# Phase 16 follow-up — re-run the bee benchmark after Tier 1-5 fixes.
#
# Designed to be launched in a *detached* tmux session, e.g.:
#   tmux new-session -d -s phase16-bee \
#       'bash /ccb/salz3/kh.chao/LiftOn/benchmarks/phase16_rerun_bee.sh'
#
# Then monitor with:
#   tmux attach -t phase16-bee     # interactive view
#   tail -f /ccb/salz3/kh.chao/LiftOn/benchmarks/results/phase16_rerun_*.log
#
# Why bee only: the diagnostic ran 5/5 datasets to ERROR with the same
# root cause, so a single dataset (bee, ~2.5 min wall-clock) is enough
# to validate the fixes. Once green, the harness can re-run all 5.
#
# This script is idempotent: re-running it is safe; the bee inputs
# already on disk will not be re-downloaded, and a successful prior
# run is gated by .lifton.done (use `--force` via PHASE16_FORCE=1 to
# override).

set -euo pipefail

ENV_BIN=/home/kh.chao/miniconda3/envs/lifton_devel/bin
REPO=/ccb/salz3/kh.chao/LiftOn

if [[ ! -x "$ENV_BIN/python" ]]; then
    echo "[phase16] ERROR: env binary not found: $ENV_BIN/python" >&2
    echo "[phase16]        edit ENV_BIN at the top of this script." >&2
    exit 2
fi

# Critical: prepend the lifton_devel bin to PATH so the harness's
# subprocess.run(["lifton", ...]) resolves to THIS env's lifton (which
# has the Phase 16 fixes), not the base conda env's stale shim
# (which dies at import time on a numpy/parasail incompatibility).
export PATH="$ENV_BIN:$PATH"
export PYTHONNOUSERSITE=1

# Phase 17a-1 attempted to opt-into gffbase here so that
# `lifton/parallel.py:_backend_supports_threads` would return True and
# the parallel Step 7 would engage on bee. That experiment failed in a
# discoverable way: gffbase's per-feature children() path is far slower
# than gffutils' SQLite B-tree on real RefSeq input — Step 3
# (extract_features_to_fasta) progressed by 0 features in 13 minutes
# during the verification run, vs. seconds under gffutils. The perf
# regression dwarfs the parallelism gain on bee. A future phase
# (Phase 12 follow-up) needs to either (a) fix the gffbase children()
# regression so it's a viable default, or (b) complete the
# materialisation refactor in lifton/locus_pipeline.py so workers no
# longer touch ctx.ref_db at all (then the guard can be loosened).
# Until one of those lands, the benchmark uses the gffutils default and
# accepts the serial Step 7 fallback.
# (No LIFTON_USE_GFFBASE export here.)

TS=$(date -u +%Y%m%dT%H%M%SZ)
LOG_DIR="$REPO/benchmarks/results"
LOG="$LOG_DIR/phase16_rerun_${TS}.log"
mkdir -p "$LOG_DIR"

# All output (stdout AND stderr) goes to both console and log file.
exec > >(tee -a "$LOG") 2>&1

cd "$REPO"

echo "============================================================"
echo " Phase 16 re-run — bee dataset"
echo " UTC start:    $TS"
echo " Repo:         $REPO"
echo " Env:          $ENV_BIN"
echo " Combined log: $LOG"
echo "============================================================"
echo

echo "[1/4] Installing mappy (Phase 16 Tier 5 dependency)"
echo "      Required by --native; without it Liftoff falls back to"
echo "      the legacy subprocess path with a stderr warning per call."
"$ENV_BIN/pip" install --no-input --quiet mappy
"$ENV_BIN/python" -c "import mappy; print(f'      mappy {mappy.__version__} OK ({mappy.__file__})')"
echo

echo "[2/4] Sanity check: harness + run_liftoff regression tests"
"$ENV_BIN/python" -m pytest \
    tests/test_benchmark_harness.py \
    tests/test_run_liftoff_traceback.py \
    -q
echo

echo "[3/4] Re-running bee benchmark"
echo "      Expected wall-clock: ~2.5 min on a typical Slurm node."
echo "      The Tier 1 parser fix means the summary table will show"
echo "      real numeric values (mapped/lost/identity) instead of"
echo "      'ERROR: could not convert...'"
echo
FORCE_ARG=""
if [[ "${PHASE16_FORCE:-0}" == "1" ]]; then
    FORCE_ARG="--force"
fi
"$ENV_BIN/python" benchmarks/run_benchmarks.py --datasets bee $FORCE_ARG

echo
echo "[4/4] Post-run inventory"
echo "------------------------------------------------------------"
ls -la "$REPO/benchmarks/results/bee/" 2>/dev/null || true
echo
echo "Output GFF3 (head -5):"
head -n 5 "$REPO/benchmarks/results/bee/lifton.gff3" 2>/dev/null \
    || echo "  (lifton.gff3 not produced — check stderr.log tail above)"
echo
echo "stderr.log line count (should be << 100K with Tier 4 fix):"
wc -l "$REPO/benchmarks/results/bee/logs/lift.stderr.log" 2>/dev/null \
    || echo "  (no stderr.log)"
echo
echo "Last 60 lines of stderr.log:"
tail -n 60 "$REPO/benchmarks/results/bee/logs/lift.stderr.log" 2>/dev/null \
    || echo "  (no stderr.log)"
echo
echo "============================================================"
echo " Phase 16 re-run finished at $(date -u +%Y%m%dT%H%M%SZ)"
echo " Combined log: $LOG"
echo "============================================================"
