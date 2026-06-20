#!/usr/bin/env bash
# benchmarks/run_tiers_4way.sh  (v2 — maize decoupled)
# ---------------------------------------------------------------------------
# Resumable, disconnect-safe driver for the 12-cell tier 4-way comparison
# (Liftoff / miniprot / LiftOn v1.0.8 stable / LiftOn devel = v2.0.0).
#
# t1_maize_b73_to_mo17 is the long pole: its largest chromosome is ~313 Mb and
# its target (Mo17 T2T) is ~2.2 Gb, so its subset minimap2 + liftoff dominate.
# We therefore DECOUPLE maize so it never blocks the other 11:
#
#   bg       build maize inputs in parallel (writes only work/maize/, no race)
#   Phase 1  fourway_compare <11 light cells>   -> report + figures  (~1 h)
#   Phase 2  (wait maize build) fourway_compare <maize> -> full report
#
# Exactly ONE fourway_compare runs at a time (Phase 1 fully finishes before
# Phase 2 starts), so the unlocked fourway_results.json write is race-free.
# Idempotent: re-running reuses .done build stages and re-scores cells.
#
# Run detached so it survives disconnect:
#   tmux new-session -d -s lifton_tiers "bash /ccb/salz3/kh.chao/LiftOn/benchmarks/run_tiers_4way.sh"
#   tmux attach -t lifton_tiers
# ---------------------------------------------------------------------------
set -u
cd /ccb/salz3/kh.chao/LiftOn || exit 2
export PYTHONNOUSERSITE=1 PYTHONHASHSEED=0
export PATH=/home/kh.chao/miniconda3/envs/lifton_devel/bin:/ccb/sw/bin:/home/kh.chao/bin:$PATH
PY=/home/kh.chao/miniconda3/envs/lifton_devel/bin/python
mkdir -p benchmarks/logs
TS=$(date +%Y%m%d_%H%M%S)
MAIN=benchmarks/logs/tier_pipeline_${TS}.log
exec > >(tee -a "$MAIN") 2>&1

echo "############################################################"
echo "### TIER 4-WAY PIPELINE v2 (maize decoupled)  start=$(date)  ts=$TS"
echo "### main log : $MAIN"
echo "############################################################"

# 11 light cells; soybean (only un-warmed light cell) ordered LAST so it builds
# inline via fourway's _ensure_inputs after the 10 cached cells have scored.
LIGHT11="t1_tomato_microtom_to_heinz t2_tomato_to_potato \
t2_human_to_gorilla t2_mouse_to_caroli t3_human_to_macaque t3_dog_to_cat \
t3_human_to_marmoset t4_human_to_chicken t4_human_to_xenopus t4_drosophila_to_bee \
t1_soybean_w82_to_lee"
MAIZE="t1_maize_b73_to_mo17"

# --- background: build maize inputs (long pole; overlaps Phase 1; no shared file) ---
echo "=== bg: build maize inputs (parallel, long pole) $(date) ==="
$PY -m benchmarks.compare.build_inputs $MAIZE -t 8 --no-force \
    > benchmarks/logs/prewarm_maize_${TS}.log 2>&1 &
MZ=$!

# --- Phase 1: score the 11 light cells (single writer) ---
echo; echo "=== PHASE 1: fourway_compare over 11 light cells $(date) ==="
$PY -m benchmarks.compare.fourway_compare $LIGHT11
echo "=== PHASE 1 scoring done rc=$? $(date) ==="
$PY -m benchmarks.compare.master_report; echo "  master_report(P1) rc=$?"
$PY -m benchmarks.compare.blog_figures;  echo "  blog_figures(P1)  rc=$?"
touch benchmarks/logs/TIER_PHASE1_DONE_${TS}
echo "=== PHASE 1 COMPLETE (11 light cells scored + report) $(date) ==="
$PY - <<'PYEOF'
import json
d = json.load(open('benchmarks/compare/fourway_results.json'))
t = sorted(k for k in d if k.split(':',1)[1][:3] in ('t1_','t2_','t3_','t4_'))
print(f"  [P1] fourway_results tier_keys={len(t)}")
for k in t:
    r=d[k]; mp=r.get('mean_pi',{}); cc=r.get('completeness_coding',{})
    print(f"    {k}: devel_meanpi={mp.get('lifton_devel')} devel_compl={cc.get('lifton_devel')} tools={r.get('tools')}")
PYEOF

# --- Phase 2: maize (after its parallel build; still single writer) ---
echo; echo "=== PHASE 2: wait maize build then score $(date) ==="
wait $MZ; echo "  maize build_inputs rc=$? $(date)"
$PY -m benchmarks.compare.fourway_compare $MAIZE
echo "=== PHASE 2 scoring done rc=$? $(date) ==="
$PY -m benchmarks.compare.master_report; echo "  master_report(P2) rc=$?"
$PY -m benchmarks.compare.blog_figures;  echo "  blog_figures(P2)  rc=$?"

$PY - <<'PYEOF'
import json
d = json.load(open('benchmarks/compare/fourway_results.json'))
t = sorted(k for k in d if k.split(':',1)[1][:3] in ('t1_','t2_','t3_','t4_'))
print(f"  [FINAL] fourway_results total={len(d)} tier_keys={len(t)}")
PYEOF

touch benchmarks/logs/TIER_PIPELINE_DONE_${TS}
echo "############################################################"
echo "### PIPELINE COMPLETE (all 12)  end=$(date)  marker=TIER_PIPELINE_DONE_${TS}"
echo "############################################################"
