#!/usr/bin/env bash
# Phase 17c follow-on — generic re-run wrapper, parameterised by $DATASET.
#
# Designed to be launched in a *detached* tmux session, e.g.:
#   tmux new-session -d -s phase17c-rice \
#       'DATASET=rice PHASE16_FORCE=1 \
#        bash /ccb/salz3/kh.chao/LiftOn/benchmarks/phase17_rerun.sh'
#
#   tmux new-session -d -s phase17c-bee \
#       'DATASET=bee PHASE16_FORCE=1 \
#        bash /ccb/salz3/kh.chao/LiftOn/benchmarks/phase17_rerun.sh'
#
# Then monitor with:
#   tmux attach -t phase17c-${DATASET}     # interactive view; detach Ctrl-b d
#   tail -f /ccb/salz3/kh.chao/LiftOn/benchmarks/results/phase17_rerun_*.log
#
# Environment variables:
#   DATASET         dataset id from benchmarks/datasets.json (default: bee).
#                   Recognised values: bee, rice, mouse, human, arabidopsis.
#   PHASE16_FORCE   set to 1 to bypass the .lifton.done short-circuit
#                   (forces a re-run even if a prior run completed).
#
# Inherits the Phase 16 wrapper's PATH + PYTHONNOUSERSITE invariants so
# subprocess.run(["lifton", ...]) resolves to the lifton_devel env's
# lifton (which carries the Phase 16/17 fixes), not the base conda
# env's stale shim.

set -euo pipefail

ENV_BIN=/home/kh.chao/miniconda3/envs/lifton_devel/bin
REPO=/ccb/salz3/kh.chao/LiftOn
DATASET="${DATASET:-bee}"

if [[ ! -x "$ENV_BIN/python" ]]; then
    echo "[phase17] ERROR: env binary not found: $ENV_BIN/python" >&2
    echo "[phase17]        edit ENV_BIN at the top of this script." >&2
    exit 2
fi

# Critical (Phase 16): prepend the lifton_devel bin to PATH so the
# harness's subprocess.run(["lifton", ...]) resolves to THIS env's
# lifton, not the base conda env's stale shim (which dies at import
# time on a numpy/parasail incompatibility).
export PATH="$ENV_BIN:$PATH"
export PYTHONNOUSERSITE=1

# Phase 17a-1 attempted to opt-into gffbase here so that
# `lifton/parallel.py:_backend_supports_threads` would return True.
# That experiment failed in a discoverable way: gffbase's per-feature
# children() path is far slower than gffutils' SQLite B-tree on real
# RefSeq input. Phase 17b solved the problem differently — workers go
# through the materialised-payload + proxy-DB path, so the guard
# accepts gffutils + native without forcing gffbase. Therefore: NO
# LIFTON_USE_GFFBASE export here. Default backend remains gffutils.

TS=$(date -u +%Y%m%dT%H%M%SZ)
LOG_DIR="$REPO/benchmarks/results"
LOG="$LOG_DIR/phase17_rerun_${DATASET}_${TS}.log"
mkdir -p "$LOG_DIR"

# All output (stdout AND stderr) goes to both console and log file.
exec > >(tee -a "$LOG") 2>&1

cd "$REPO"

echo "============================================================"
echo " Phase 17 re-run — ${DATASET} dataset"
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

echo "[2/4] Sanity check: harness + run_liftoff + locus-pipeline regression tests"
"$ENV_BIN/python" -m pytest \
    tests/test_benchmark_harness.py \
    tests/test_run_liftoff_traceback.py \
    tests/test_locus_materialise.py \
    -q
echo

echo "[3/4] Re-running ${DATASET} benchmark"
echo "      Phase 17b parallel path is now active by default under"
echo "      --locus-pipeline -t 8 --native (the Phase 9 guard accepts"
echo "      gffutils + native via the materialised-payload + proxy-DB"
echo "      route; no LIFTON_USE_GFFBASE opt-in required)."
echo
FORCE_ARG=""
if [[ "${PHASE16_FORCE:-0}" == "1" ]]; then
    FORCE_ARG="--force"
fi
"$ENV_BIN/python" benchmarks/run_benchmarks.py --datasets "$DATASET" $FORCE_ARG

DS_DIR="$REPO/benchmarks/results/$DATASET"
echo
echo "[4/4] Post-run inventory"
echo "------------------------------------------------------------"
ls -la "$DS_DIR/" 2>/dev/null || true
echo
echo "Output GFF3 (head -5):"
head -n 5 "$DS_DIR/lifton.gff3" 2>/dev/null \
    || echo "  (lifton.gff3 not produced — check stderr.log tail above)"
echo
echo "Output GFF3 (size + line count):"
ls -la "$DS_DIR/lifton.gff3" 2>/dev/null \
    && wc -l "$DS_DIR/lifton.gff3" 2>/dev/null \
    || echo "  (no lifton.gff3)"
echo
echo "stderr.log line count (should be << 100K with Phase 16 Tier 4 fix):"
wc -l "$DS_DIR/logs/lift.stderr.log" 2>/dev/null \
    || echo "  (no stderr.log)"
echo
echo "Serial-fallback warning count (Phase 17b expects 0):"
grep -c "Falling back to serial" "$DS_DIR/logs/lift.stderr.log" 2>/dev/null \
    || echo "0 (file absent or no match)"
echo
echo "Wall-clock from /usr/bin/time -v:"
grep -E "Elapsed|Maximum resident set size|User time|System time|Percent of CPU" \
    "$DS_DIR/logs/lift.time.log" 2>/dev/null \
    || echo "  (no time.log)"
echo
echo "Last 60 lines of stderr.log:"
tail -n 60 "$DS_DIR/logs/lift.stderr.log" 2>/dev/null \
    || echo "  (no stderr.log)"
echo
echo "============================================================"
echo " Phase 17 re-run (${DATASET}) finished at $(date -u +%Y%m%dT%H%M%SZ)"
echo " Combined log: $LOG"
echo "============================================================"
