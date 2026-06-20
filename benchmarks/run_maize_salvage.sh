#!/usr/bin/env bash
set -u; cd /ccb/salz3/kh.chao/LiftOn
export PYTHONNOUSERSITE=1 PYTHONHASHSEED=0
export PATH=/home/kh.chao/miniconda3/envs/lifton_devel/bin:/ccb/sw/bin:/home/kh.chao/bin:$PATH
PY=/home/kh.chao/miniconda3/envs/lifton_devel/bin/python
TS=$(date +%Y%m%d_%H%M%S); LOG=benchmarks/logs/maize_salvage_${TS}.log
exec > >(tee -a "$LOG") 2>&1
echo "=== MAIZE SALVAGE (pinned tgt CM039150.1) start=$(date) ==="
$PY -m benchmarks.compare.fourway_compare t1_maize_b73_to_mo17
echo "=== maize scored rc=$? $(date) ==="
$PY -m benchmarks.compare.master_report; echo "master_report rc=$?"
$PY -m benchmarks.compare.blog_figures;  echo "blog_figures rc=$?"
touch benchmarks/logs/MAIZE_SALVAGE_DONE_${TS}
echo "=== MAIZE SALVAGE DONE $(date) marker=MAIZE_SALVAGE_DONE_${TS} ==="
