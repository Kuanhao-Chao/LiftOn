#!/bin/bash
# The pre-release gates that need real data and a real install:
#   * the chr22 example at -t 1 and -t 8, byte-identical to each other
#   * gff3-validate on that output
#   * wheel + sdist build
#   * a --no-cache-dir install into a fresh venv (the empty-cache path that
#     caught the broken `cigar` dependency before v1.0.10 reached PyPI)
#   * the installed wheel lifting chr22 to the SAME bytes as the tree -- the
#     smoke test alone only proves the wheel imports and answers -h
#
# Usage: bash benchmarks/release_gates.sh [output_dir]
set -u
ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
PY=${LIFTON_PY:-$(command -v python)}
BIN=$(dirname "$PY")
OUT=${1:-$ROOT/benchmarks/compare/_v1012/release_gates}
rm -rf "$OUT"; mkdir -p "$OUT"

echo "== chr22 example, -t 1 and -t 8 =="
cd "$ROOT/test"
for t in 1 8; do
  rm -rf "$OUT/chr22_t$t"; mkdir -p "$OUT/chr22_t$t"
  /usr/bin/time -v -o "$OUT/chr22_t$t.time" \
    "$BIN/lifton" -t $t -g GRCh38_chr22.gff3 -o "$OUT/chr22_t$t/lifton.gff3" \
      -copies -sc 0.95 chm13_chr22.fa GRCh38_chr22.fa > "$OUT/chr22_t$t.log" 2>&1
  echo "[chr22 t=$t] exit=$?  $(md5sum "$OUT/chr22_t$t/lifton.gff3" 2>/dev/null | cut -d' ' -f1)"
done
"$BIN/gff3-validate" "$OUT/chr22_t1/lifton.gff3" > "$OUT/chr22_validate.txt" 2>&1
echo "[gff3-validate] exit=$?"
grep -iE "^(errors|warnings)|is_valid|VALID" "$OUT/chr22_validate.txt" | head -5

echo "== wheel + sdist build =="
cd "$ROOT"
rm -rf "$OUT/dist"
"$PY" -m build --outdir "$OUT/dist" > "$OUT/build.log" 2>&1
echo "[build] exit=$?"
ls "$OUT/dist"

echo "== clean-venv smoke =="
"$PY" -m venv "$OUT/venv" > "$OUT/venv.log" 2>&1
WHL=$(ls "$OUT"/dist/*.whl 2>/dev/null | head -1)
"$OUT/venv/bin/pip" install --no-cache-dir "$WHL" >> "$OUT/venv.log" 2>&1
echo "[pip install] exit=$?"
"$OUT/venv/bin/lifton" --version;      echo "[--version] exit=$?"
"$OUT/venv/bin/lifton" -h > /dev/null; echo "[-h] exit=$?"
"$OUT/venv/bin/gff3-validate" -h > /dev/null; echo "[gff3-validate -h] exit=$?"
"$OUT/venv/bin/lifton" -h 2>&1 | grep -cE "no-orf-stop-completion|intermediate-dir" | sed 's/^/[new flags in -h] /'

# Importing and answering -h does not prove the packaged artifact LIFTS the same.
# Run the chr22 example from the venv and compare with the tree's own output.
echo "== the wheel lifts identically =="
cd "$ROOT/test"
rm -rf "$OUT/wheel_chr22"; mkdir -p "$OUT/wheel_chr22"
"$OUT/venv/bin/lifton" -t 8 -g GRCh38_chr22.gff3 -o "$OUT/wheel_chr22/lifton.gff3" \
  -copies -sc 0.95 chm13_chr22.fa GRCh38_chr22.fa > "$OUT/wheel_chr22.log" 2>&1
echo "[wheel chr22] exit=$?"
WHEEL_MD5=$(md5sum "$OUT/wheel_chr22/lifton.gff3" 2>/dev/null | cut -d" " -f1)
TREE_MD5=$(md5sum "$OUT/chr22_t8/lifton.gff3" 2>/dev/null | cut -d" " -f1)
echo "[wheel md5] $WHEEL_MD5"
echo "[tree  md5] $TREE_MD5"
[ -n "$WHEEL_MD5" ] && [ "$WHEEL_MD5" = "$TREE_MD5" ] \
  && echo "[wheel == tree] yes" || echo "[wheel == tree] NO -- investigate"
echo "[DONE] exit=$?"
