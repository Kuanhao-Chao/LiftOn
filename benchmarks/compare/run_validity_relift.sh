#!/usr/bin/env bash
# Track-7 validity re-lift launcher — runs the 17-genome devel re-lift + revalidate
# inside tmux so it survives session disconnects. The driver writes the JSON
# incrementally (atomic per genome), so a restart only re-does in-flight genomes.
set -uo pipefail
cd /ccb/salz3/kh.chao/LiftOn
export PYTHONNOUSERSITE=1 PYTHONHASHSEED=0 PATH="/ccb/sw/bin:/home/kh.chao/bin:$PATH"
LOG="benchmarks/logs/validity_relift_tmux_$(date +%Y%m%d_%H%M%S).log"
echo "$LOG" > /tmp/lifton_relift_logpath.txt
echo "[run_validity_relift] $(date) logging to $LOG"
/home/kh.chao/miniconda3/envs/lifton_devel/bin/python -m benchmarks.compare.validity_relift \
    --workers 4 --threads 8 2>&1 | tee "$LOG"
rc=${PIPESTATUS[0]}
echo "[run_validity_relift] $(date) EXIT $rc" | tee -a "$LOG"
echo "$rc" > /tmp/lifton_relift_exit.txt
