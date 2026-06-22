#!/usr/bin/env bash
# Track-7: launch each of the 17 full-genome devel re-lifts in its OWN tmux
# session (relift_<bid>), all in parallel. Each writes a per-genome result file
# to _validity_relift/<bid>.result.json + a <bid>.done sentinel (no shared-JSON
# race). Once all 17 .done land, run the merge (validity_relift --merge) to fold
# them into fourway_results.json.
set -u
cd /ccb/salz3/kh.chao/LiftOn
THREADS=${1:-4}
OUT=benchmarks/compare/_validity_relift
SINGLE=/ccb/salz3/kh.chao/LiftOn/benchmarks/compare/run_validity_single.sh
mkdir -p "$OUT"
rm -f "$OUT"/*.result.json "$OUT"/*.done "$OUT"/*.log 2>/dev/null   # fresh start

IDS="t4_drosophila_to_bee drosophila arabidopsis_to_rice human_to_zebrafish \
t4_human_to_xenopus t4_human_to_chicken t2_tomato_to_potato \
t1_tomato_microtom_to_heinz t1_maize_b73_to_mo17 t3_dog_to_cat \
t2_mouse_to_caroli t3_human_to_macaque t3_human_to_marmoset t2_human_to_gorilla \
bee arabidopsis rice"

n=0
for bid in $IDS; do
  sess="relift_${bid}"
  tmux kill-session -t "$sess" 2>/dev/null
  tmux new-session -d -s "$sess" "bash $SINGLE $bid $THREADS > $OUT/${bid}.log 2>&1"
  n=$((n + 1))
  echo "launched tmux $sess (threads=$THREADS)"
done
echo "LAUNCHED $n per-genome tmux sessions (threads=$THREADS each). Merge with:"
echo "  python -m benchmarks.compare.validity_relift --merge --out-dir $OUT"
