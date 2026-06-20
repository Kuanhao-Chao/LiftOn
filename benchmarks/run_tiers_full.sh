#!/usr/bin/env bash
# benchmarks/run_tiers_full.sh
# ---------------------------------------------------------------------------
# FULL-GENOME (end-to-end, whole-genome) 4-way comparison for all 12 tier pairs.
# Each cell runs `fourway_compare --full <id>` (Liftoff / miniprot / LiftOn v1.0.8 /
# LiftOn devel=v2.0.0) on the WHOLE ref+tgt genomes from benchmarks.json.
#
# Parallel-safe: each cell writes its OWN results file via FOURWAY_RESULTS_JSON
# (benchmarks/compare/_full_runs/<id>.json) -- NO race on the shared
# fourway_results.json. Merge happens in a separate step after completion.
#
# Concurrency cap = 4 (mix of genome sizes; RAM up to ~30 GB/cell, ~600 GB free).
# Crash-tolerant: a cell that fails (e.g. soybean's standalone-Liftoff KeyError,
# or a v1.0.8 crash) is logged and skipped; siblings continue.
#
# Disconnect-safe: run inside a detached tmux session:
#   tmux new-session -d -s lifton_full "bash /ccb/salz3/kh.chao/LiftOn/benchmarks/run_tiers_full.sh"
#   tmux attach -t lifton_full
# ---------------------------------------------------------------------------
set -u
cd /ccb/salz3/kh.chao/LiftOn || exit 2
export PYTHONNOUSERSITE=1 PYTHONHASHSEED=0
export PATH=/home/kh.chao/miniconda3/envs/lifton_devel/bin:/ccb/sw/bin:/home/kh.chao/bin:$PATH
PY=/home/kh.chao/miniconda3/envs/lifton_devel/bin/python
mkdir -p benchmarks/logs benchmarks/compare/_full_runs
TS=$(date +%Y%m%d_%H%M%S)
MAIN=benchmarks/logs/tier_full_pipeline_${TS}.log
exec > >(tee -a "$MAIN") 2>&1
MAXJOBS=4

echo "############################################################"
echo "### TIER FULL-GENOME 4-WAY   start=$(date)   ts=$TS   maxjobs=$MAXJOBS"
echo "### main log: $MAIN"
echo "############################################################"

# cheap -> expensive (1 Gb first, 3 Gb human pairs last)
IDS="t4_drosophila_to_bee t2_tomato_to_potato t1_tomato_microtom_to_heinz \
t1_soybean_w82_to_lee t1_maize_b73_to_mo17 t2_mouse_to_caroli t3_dog_to_cat \
t2_human_to_gorilla t3_human_to_macaque t3_human_to_marmoset \
t4_human_to_chicken t4_human_to_xenopus"

run_cell(){
  local id=$1
  local rj="benchmarks/compare/_full_runs/${id}.json"
  local lg="benchmarks/logs/full_${id}_${TS}.log"
  echo "[launch] $id  $(date)  -> $lg"
  FOURWAY_RESULTS_JSON="$rj" $PY -m benchmarks.compare.fourway_compare --full "$id" > "$lg" 2>&1
  echo "[done]   $id  rc=$?  $(date)"
}

for id in $IDS; do
  while [ "$(jobs -rp | wc -l)" -ge "$MAXJOBS" ]; do sleep 30; done
  run_cell "$id" &
  sleep 5   # stagger launches so provenance/index setup don't thunder
done
wait
echo "=== all full-genome cells finished $(date) ==="

# tally what produced a result
echo "--- per-cell result files ---"
for id in $IDS; do
  rj="benchmarks/compare/_full_runs/${id}.json"
  if [ -s "$rj" ]; then
    $PY -c "import json;d=json.load(open('$rj'));k=[x for x in d if x.startswith('full:')];print('  OK  $id ->',k, {t:(round(d[k[0]]['mean_pi'].get(t),4) if d[k[0]]['mean_pi'].get(t) else None) for t in d[k[0]]['tools']} if k else {})" 2>/dev/null || echo "  ??  $id (unreadable json)"
  else
    echo "  FAIL $id (no result file -- check benchmarks/logs/full_${id}_${TS}.log)"
  fi
done

touch benchmarks/logs/TIER_FULL_DONE_${TS}
echo "############################################################"
echo "### TIER FULL-GENOME DONE   end=$(date)   marker=TIER_FULL_DONE_${TS}"
echo "### NEXT: merge benchmarks/compare/_full_runs/*.json into fourway_results.json"
echo "############################################################"
