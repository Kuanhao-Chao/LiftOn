#!/usr/bin/env bash
# Phase 16 — convenience wrapper around run_benchmarks.py.
#
# Activates the lifton-test conda env (best effort), runs the
# benchmark harness, and tees stdout/stderr to a timestamped log.
#
# Usage examples:
#   ./benchmarks/run_benchmarks.sh                       # all 5 datasets
#   ./benchmarks/run_benchmarks.sh --datasets bee rice   # subset
#   ./benchmarks/run_benchmarks.sh --download-only       # pre-fetch only
#
# Any extra args are forwarded verbatim to run_benchmarks.py.

set -euo pipefail

HERE="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
REPO_ROOT="$( cd "$HERE/.." >/dev/null 2>&1 && pwd )"

# Best-effort: activate the lifton-test conda env if it exists.
CONDA_SH="/opt/anaconda3/etc/profile.d/conda.sh"
if [[ -f "$CONDA_SH" ]]; then
    # shellcheck disable=SC1090
    source "$CONDA_SH"
    if conda env list | grep -qE "^lifton-test\b"; then
        conda activate lifton-test
    fi
fi

mkdir -p "$HERE/results"
LOG="$HERE/results/run_$(date -u +%Y%m%dT%H%M%SZ).log"

echo "[bench] writing combined log to $LOG"
cd "$REPO_ROOT"
python "$HERE/run_benchmarks.py" "$@" 2>&1 | tee "$LOG"
