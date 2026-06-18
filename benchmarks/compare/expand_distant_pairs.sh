#!/usr/bin/env bash
# =============================================================================
# Expand the 4-way benchmark with 7 distantly-related cross-species pairs.
#
#   plant : arabidopsis_to_lyrata (close) | rice_to_sorghum (distant) |
#           arabidopsis_to_rice  (very-distant, eudicot->monocot)   [+ --full]
#   fish  : zebrafish_to_medaka  (distant) | human_to_zebrafish     [+ --full]
#           (very-distant, mammal->fish; human chr20, controlled vs human_to_mouse)
#   insect: drosophila_to_anopheles (very-distant, fly->mosquito)
#   fungi : cerevisiae_to_pombe      (very-distant, ~500 My)
#
# Refs reuse on-disk genomes; only TARGETS are downloaded. Very-distant pairs use
# tgt_chrom=WHOLE (subset_builder WHOLE mode) because genome-wide DNA synteny is
# gone — the reference is still subset to one chromosome, the whole target genome
# goes to Liftoff + miniprot.
#
# Resilient: per-pair failures are logged, not fatal (fourway_compare also wraps
# each id in try/except). A subset-only report is written before the heavy --full
# flagships so partial results are always captured.
#
# Run detached so it survives disconnect:
#   tmux new -d -s lifton_bench \
#     'bash benchmarks/compare/expand_distant_pairs.sh 2>&1 | tee benchmarks/compare/expand_distant_pairs.log'
# =============================================================================
set -uo pipefail

REPO=/ccb/salz3/kh.chao/LiftOn
PY=/home/kh.chao/miniconda3/envs/lifton_devel/bin/python
export PYTHONNOUSERSITE=1 PYTHONHASHSEED=0
cd "$REPO"

# pairs whose target genome must be downloaded (REUSE-only pairs are skipped)
FETCH_IDS="arabidopsis_to_lyrata rice_to_sorghum zebrafish_to_medaka drosophila_to_anopheles cerevisiae_to_pombe"
# all 7 new pairs, subset 4-way
ALL_IDS="arabidopsis_to_lyrata rice_to_sorghum arabidopsis_to_rice zebrafish_to_medaka human_to_zebrafish drosophila_to_anopheles cerevisiae_to_pombe"
# full-genome flagships (lighter plant first, heavy human/zebrafish last)
FULL_IDS="arabidopsis_to_rice human_to_zebrafish"

ts(){ date "+%Y-%m-%d %H:%M:%S"; }
banner(){ echo; echo "######## $(ts)  $* ########"; echo; }

banner "START expand_distant_pairs  (host=$(hostname), pid=$$)"

# --- Phase 1: fetch target genomes -----------------------------------------
banner "PHASE 1/6: fetch target genomes ($FETCH_IDS)"
"$PY" -m benchmarks.compare.fetch_new_pairs $FETCH_IDS \
  || echo "!! fetch_new_pairs returned nonzero (continuing; missing inputs surface per-pair below)"

# --- Phase 2: build inputs (subset + liftoff + miniprot) per pair -----------
banner "PHASE 2/6: build inputs per pair"
for id in $ALL_IDS; do
  banner "build_inputs $id"
  if "$PY" -m benchmarks.compare.build_inputs "$id" -t 8 --no-force; then
    echo ">> build_inputs OK: $id"
  else
    echo "!! build_inputs FAILED: $id (continuing with remaining pairs)"
  fi
done

# --- Phase 3: 4-way subset comparison (all 7) ------------------------------
banner "PHASE 3/6: 4-way subset comparison"
"$PY" -m benchmarks.compare.fourway_compare $ALL_IDS \
  || echo "!! fourway_compare (subset) returned nonzero"

# --- Phase 4: intermediate report (subset results captured) ----------------
banner "PHASE 4/6: regenerate report + figures (subset checkpoint)"
"$PY" -m benchmarks.compare.master_report  || echo "!! master_report returned nonzero"
"$PY" -m benchmarks.compare.report_figures || echo "!! report_figures returned nonzero"

# --- Phase 5: full-genome flagships ----------------------------------------
banner "PHASE 5/6: 4-way FULL flagships ($FULL_IDS)"
"$PY" -m benchmarks.compare.fourway_compare --full $FULL_IDS \
  || echo "!! fourway_compare (--full) returned nonzero"

# --- Phase 6: final report -------------------------------------------------
banner "PHASE 6/6: regenerate report + figures (final, incl. fulls)"
"$PY" -m benchmarks.compare.master_report  || echo "!! master_report returned nonzero"
"$PY" -m benchmarks.compare.report_figures || echo "!! report_figures returned nonzero"

banner "ALL DONE"
echo "results : $REPO/benchmarks/compare/fourway_results.json"
echo "report  : $REPO/benchmarks/compare/MASTER_COMPARISON_REPORT.md"
echo "figures : $REPO/benchmarks/compare/figures/  (+ report rfig_*.png)"
