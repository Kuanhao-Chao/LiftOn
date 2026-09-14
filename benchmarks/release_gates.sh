#!/bin/bash
# Run real chr22 lifts, compare serial/threaded output, build distributions,
# install the wheel without a pip cache, and prove the installed wheel lifts
# identically. Every failed stage is a failed gate; prior evidence is retained.
# Usage: LIFTON_PY=/path/to/python bash benchmarks/release_gates.sh [new_output_dir]
set -euo pipefail
ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
PY=${LIFTON_PY:-$(command -v python)}
PY=$(command -v "$PY")
BIN=$(dirname "$PY")
if [[ $# -gt 0 ]]; then
  # mkdir must fail on an existing directory, including an empty one. Never
  # erase earlier evidence or accidentally reuse a stale successful output.
  mkdir -p "$(dirname "$1")"
  mkdir "$1"
  OUT=$(cd "$1" && pwd)
else
  mkdir -p "$ROOT/benchmarks/compare/_runs"
  OUT=$(mktemp -d "$ROOT/benchmarks/compare/_runs/release-gates.XXXXXXXX")
fi
trap 'status=$?; if (( status != 0 )); then echo "[FAILED] exit=$status; evidence: $OUT" >&2; fi' EXIT
printf 'Evidence: %s\n' "$OUT"
export PYTHONNOUSERSITE=1

run_chr22() {
  local executable=$1
  local threads=$2
  local destination=$3
  mkdir "$destination"
  "$executable" -t "$threads" -g "$ROOT/test/GRCh38_chr22.gff3" \
    -o "$destination/lifton.gff3" -copies -sc 0.95 \
    "$ROOT/test/chm13_chr22.fa" "$ROOT/test/GRCh38_chr22.fa"
  test -s "$destination/lifton.gff3"
}

echo '== chr22 example, -t 1 and -t 8 =='
cd "$OUT"
for t in 1 8; do
  run_chr22 "$BIN/lifton" "$t" "$OUT/chr22_t$t" > "$OUT/chr22_t$t.log" 2>&1
  sha256sum "$OUT/chr22_t$t/lifton.gff3"
done
cmp "$OUT/chr22_t1/lifton.gff3" "$OUT/chr22_t8/lifton.gff3"
"$BIN/gff3-validate" "$OUT/chr22_t1/lifton.gff3" > "$OUT/chr22_validate.txt" 2>&1

echo '== wheel + sdist build =='
cd "$ROOT"
"$PY" -m build --outdir "$OUT/dist" > "$OUT/build.log" 2>&1
shopt -s nullglob
wheels=("$OUT"/dist/*.whl)
sdists=("$OUT"/dist/*.tar.gz)
[[ ${#wheels[@]} -eq 1 && ${#sdists[@]} -eq 1 ]]
sha256sum "${wheels[0]}" "${sdists[0]}" > "$OUT/distributions.sha256"

echo '== clean-venv smoke =='
# Run outside the source tree with no PYTHONPATH/user-site leakage.
cd "$OUT"
unset PYTHONPATH
"$PY" -m venv "$OUT/venv" > "$OUT/venv.log" 2>&1
"$OUT/venv/bin/pip" install --no-cache-dir "${wheels[0]}" >> "$OUT/venv.log" 2>&1
"$OUT/venv/bin/lifton" --version
"$OUT/venv/bin/lifton" -h > "$OUT/wheel_help.txt"
"$OUT/venv/bin/gff3-validate" -h > "$OUT/validator_help.txt"
grep -q -- '--no-orf-stop-completion' "$OUT/wheel_help.txt"
grep -q -- '--intermediate-dir' "$OUT/wheel_help.txt"

echo '== installed wheel lifts identically =='
run_chr22 "$OUT/venv/bin/lifton" 8 "$OUT/wheel_chr22" > "$OUT/wheel_chr22.log" 2>&1
"$OUT/venv/bin/gff3-validate" "$OUT/wheel_chr22/lifton.gff3" > "$OUT/wheel_validate.txt" 2>&1
cmp "$OUT/chr22_t8/lifton.gff3" "$OUT/wheel_chr22/lifton.gff3"
sha256sum "$OUT/wheel_chr22/lifton.gff3"
echo '[DONE] exit=0'
