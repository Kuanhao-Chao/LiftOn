#!/usr/bin/env bash
# Per-genome devel re-lift + revalidate (one tmux session per genome).
# Usage: run_validity_single.sh <benchmark_id> [threads]
set -u
bid="$1"; threads="${2:-4}"
cd /ccb/salz3/kh.chao/LiftOn
export PYTHONNOUSERSITE=1 PYTHONHASHSEED=0 PATH="/ccb/sw/bin:/home/kh.chao/bin:$PATH"
OUT=benchmarks/compare/_validity_relift
echo "[$bid] $(date) starting single re-lift (threads=$threads)"
/home/kh.chao/miniconda3/envs/lifton_devel/bin/python -m benchmarks.compare.validity_relift \
    --single "$bid" --threads "$threads" --out-dir "$OUT"
echo "[$bid] $(date) finished single re-lift (rc=$?)"
