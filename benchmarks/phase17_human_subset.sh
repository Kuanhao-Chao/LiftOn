#!/usr/bin/env bash
# Phase 17d Step 2 — human chromosome subset benchmark.
#
# Direct lifton invocation with `-chroms <file>` so the run scopes to a
# tractable chromosome subset of the full GRCh38 → CHM13 lift-over.
# Bypasses run_benchmarks.py (which has no per-dataset flag injection)
# and writes its own results dir + log under benchmarks/results/.
#
# Designed to be launched in a *detached* tmux session, e.g.:
#   tmux new-session -d -s phase17d-human-subset \
#       'bash /ccb/salz3/kh.chao/LiftOn/benchmarks/phase17_human_subset.sh'
#
# The chromosome subset (chr19, chr22, chrY, chrM — see
# `benchmarks/human_subset.chroms`) is ~180 MB / ~5 % of the GRCh38
# genome but exercises gene-dense regions + sex chromosome +
# mitochondrial logic + long human protein isoforms (which test the
# per-locus parasail O(P²) hot-spot more heavily than bee/rice).

set -euo pipefail

ENV_BIN=/home/kh.chao/miniconda3/envs/lifton_devel/bin
REPO=/ccb/salz3/kh.chao/LiftOn
DATA="$REPO/benchmarks/data/human"
OUTDIR="$REPO/benchmarks/results/human_subset"
CHROMS="$REPO/benchmarks/human_subset.chroms"

if [[ ! -x "$ENV_BIN/python" ]]; then
    echo "[phase17d] ERROR: env binary not found: $ENV_BIN/python" >&2
    exit 2
fi
if [[ ! -f "$CHROMS" ]]; then
    echo "[phase17d] ERROR: chroms file not found: $CHROMS" >&2
    exit 2
fi
for f in "$DATA/NCBI_RefSeq_no_rRNA.gff" \
         "$DATA/chm13v2.0.fa" \
         "$DATA/GCF_000001405.40_GRCh38.p14_genomic.fna" ; do
    if [[ ! -f "$f" ]]; then
        echo "[phase17d] ERROR: input not found: $f" >&2
        exit 2
    fi
done

# Phase 16 PATH discipline — resolve `lifton` to lifton_devel's binary,
# not the base conda env's stale shim.
export PATH="$ENV_BIN:$PATH"
export PYTHONNOUSERSITE=1

TS=$(date -u +%Y%m%dT%H%M%SZ)
LOG_DIR="$REPO/benchmarks/results"
LOG="$LOG_DIR/phase17d_human_subset_${TS}.log"
mkdir -p "$LOG_DIR" "$OUTDIR" "$OUTDIR/logs"

# All output goes to both console AND the combined log.
exec > >(tee -a "$LOG") 2>&1

cd "$REPO"

echo "============================================================"
echo " Phase 17d Step 2 — human chr19+22+Y+M subset"
echo " UTC start:    $TS"
echo " Repo:         $REPO"
echo " Env:          $ENV_BIN"
echo " Output dir:   $OUTDIR"
echo " Combined log: $LOG"
echo " Chroms file:"
sed 's/^/   /' "$CHROMS"
echo "============================================================"
echo

echo "[1/4] Sanity check: harness + run_liftoff + locus_pipeline regression tests"
"$ENV_BIN/python" -m pytest \
    tests/test_benchmark_harness.py \
    tests/test_run_liftoff_traceback.py \
    tests/test_locus_materialise.py \
    -q
echo

echo "[2/4] Confirming mappy availability (Phase 16 Tier 5 dep)"
"$ENV_BIN/python" -c "import mappy; print(f'      mappy {mappy.__version__} OK')"
echo

echo "[3/4] Running lifton with -chroms (subset)"
echo "      Phase 17b parallel path is active under --locus-pipeline"
echo "      -t 8 --native (the materialised-payload + proxy-DB route"
echo "      means workers are DB-free; no LIFTON_USE_GFFBASE needed)."
echo
echo "      Command: lifton --stream --inmemory-liftoff --locus-pipeline"
echo "               -t 8 --native -copies -chroms <subset>"
echo "               -g <ref.gff> -o <out.gff3> <target.fa> <ref.fa>"
echo

# /usr/bin/time -v writes its summary to stderr; the -o option directs
# it to a file instead, leaving the lift's own stderr clean.
/usr/bin/time -v -o "$OUTDIR/logs/lift.time.log" \
  "$ENV_BIN/lifton" \
    --stream --inmemory-liftoff --locus-pipeline -t 8 --native -copies \
    -chroms "$CHROMS" \
    -g "$DATA/NCBI_RefSeq_no_rRNA.gff" \
    -o "$OUTDIR/lifton.gff3" \
    "$DATA/chm13v2.0.fa" \
    "$DATA/GCF_000001405.40_GRCh38.p14_genomic.fna" \
    > "$OUTDIR/logs/lift.stdout.log" \
    2> "$OUTDIR/logs/lift.stderr.log"

echo
echo "[4/4] Post-run inventory"
echo "------------------------------------------------------------"
ls -la "$OUTDIR/" 2>/dev/null || true
echo
echo "Output GFF3 size + line count:"
ls -la "$OUTDIR/lifton.gff3" 2>/dev/null \
    && wc -l "$OUTDIR/lifton.gff3" 2>/dev/null \
    || echo "  (no lifton.gff3 — check stderr.log tail)"
echo
echo "Output GFF3 (head -5):"
head -n 5 "$OUTDIR/lifton.gff3" 2>/dev/null \
    || echo "  (no lifton.gff3)"
echo
echo "Wall-clock from /usr/bin/time -v:"
grep -E "Elapsed|Maximum resident set size|User time|System time|Percent of CPU|Exit status" \
    "$OUTDIR/logs/lift.time.log" 2>/dev/null \
    || echo "  (no time.log)"
echo
echo "stderr.log line count (Phase 16 Tier 4 expects << 100K):"
wc -l "$OUTDIR/logs/lift.stderr.log" 2>/dev/null \
    || echo "  (no stderr.log)"
echo
echo "Serial-fallback warning count (Phase 17b expects 0):"
grep -c "Falling back to serial" "$OUTDIR/logs/lift.stderr.log" 2>/dev/null \
    || echo "0 (file absent or no match)"
echo
echo "Last 60 lines of stderr.log:"
tail -n 60 "$OUTDIR/logs/lift.stderr.log" 2>/dev/null \
    || echo "  (no stderr.log)"
echo
echo "============================================================"
echo " Phase 17d Step 2 finished at $(date -u +%Y%m%dT%H%M%SZ)"
echo " Combined log: $LOG"
echo "============================================================"
