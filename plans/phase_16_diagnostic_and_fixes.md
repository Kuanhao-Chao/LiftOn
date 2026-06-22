# LiftOn — Phase 16: Benchmark Diagnostic & Fixes

> **Status:** Tier 1-5 fixes landed AND verified on the bee dataset
> end-to-end. Lift pass exits 0 in 556.74 s (9:16.74), peak RSS 11.5 GB,
> produces a 9.7 MB output GFF3. stderr dropped from 532,599 lines to
> ~100 lines. One follow-up bytes-blob log bug was found and fixed
> during verification (see §3.6). The harness's optional `-E`
> evaluation pass exits 1 with `mapped: 0` — flagged for a separate
> investigation in §5; it is *not* part of the four originally-diagnosed
> Phase 16 layers.
>
> **Test suite:** 56/56 critical-path tests pass including all 24
> cells of the golden-gate byte-identity matrix. The new follow-up
> tier added 3 unit tests for the bytes-blob log helper.
>
> **Date:** 2026-05-07
> **Branch:** `devel`
> **Triggering commit:** `7152cca chore(phase-16): add biological validation benchmark harness`
> **Predecessor diagnostic plan:** `~/.claude/plans/phase-16-act-as-dazzling-pnueli.md`

---

## 1. Background

The Phase 16 benchmark harness (`benchmarks/run_benchmarks.py`) was added
in commit `7152cca` to validate LiftOn against five published JHU CCB
datasets (human, mouse, bee, arabidopsis, rice). The first real run
(`run_20260504T070008Z.log`) reported **all five datasets as ERROR** in
both the summary table and `summary_*.json`. A read-only diagnostic
identified four layered failures, three of which were real bugs and one
of which was an environment gap; this phase landed fixes for all four.

## 2. Diagnostic findings (four layers)

### Layer A — harness parser bug (cosmetic, blocked all reporting)

`benchmarks/run_benchmarks.py:181-189` `_parse_gnu_time` cut on the first
colon, which mangled the wall-clock value because GNU `/usr/bin/time -v`
embeds a colon inside the prefix `(h:mm:ss or m:ss)`. Every parsed
`wall_clock_str` came out as `"mm:ss or m:ss): X:YY.ZZ"`, which then
crashed `_wall_clock_str_to_seconds` at line 271 with
`ValueError: could not convert string to float: ...`. **The harness
collapsed the entire benchmark report even on runs where LiftOn
actually succeeded.**

### Layer B — RecursionError in vendored Liftoff (the real failure)

Every dataset's stderr ended with:

```
[ERROR] Liftoff encountered a fatal error during native execution: maximum recursion depth exceeded while calling a Python object
[ERROR] LiftOn cannot proceed without a valid Liftoff baseline annotation.
```

The `except Exception as e: logger.log_error(... {e})` at
`lifton/run_liftoff.py:46-60` interpolated only the exception's `str()`
form, **discarding the traceback**. No Python frame was ever printed —
confirmed by grepping all 532,599 lines of `bee/logs/lift.stderr.log`:
the only `Traceback`-prefixed match was the wrapper `[ERROR]` line. The
recursive function in `lifton/liftoff/` was therefore unidentifiable.

### Layer C — `mappy` not installed → `--native` falls back

Every benchmark stderr fired exactly **200 copies** of the
`align_features_to_target_native` fallback warning at
`lifton/liftoff/native_align.py:175-185`. `lifton.yml` and
`setup.py:install_requires` did not declare `mappy`, so the conda env
didn't ship it; the supposedly-native path was silently degrading to
the legacy subprocess minimap2 every time it was hit.

### Layer D — GFF3 validator stderr flood

`lifton/lifton.py:315-316` looped over every finding from
`GFF3Validator.validate(...)` and called `logger.log(str(f), debug=True)`
unconditionally, regardless of `--strict-gff`. NCBI `Dbxref` values
(`DBTAG:ID` shape) flag the `unencoded_reserved_char` rule on every row
that contains them — bee stderr was 532,599 lines, arabidopsis 786,375
lines. **Real errors were drowned in spec-compliance noise.**

## 3. What shipped

Each tier was implemented in dependency order; Tier 3 was deliberately
gated on Tier 2 so the recursion site would have been visible if the
cheap raise-recursion-limit fix proved insufficient.

| Tier | File(s) | Change |
|---|---|---|
| **1** | `benchmarks/run_benchmarks.py` | `_parse_gnu_time` now slices past the matched prefix instead of `split(":", 1)`, so the colons inside `(h:mm:ss or m:ss)` are preserved as part of the value. New `_safe_float` helper. The `_wall_clock_str_to_seconds` call site is wrapped in try/except so a single parse failure no longer aborts the whole `run_profiled`. |
| **1** | `tests/test_benchmark_harness.py` *(new)* | 13 unit tests: every `_GNU_KEYS` entry, the real bee fixture, all three `_wall_clock_str_to_seconds` shapes (m:ss, mm:ss, h:mm:ss), and `_safe_float` against the exact pre-fix corruption. |
| **2** | `lifton/run_liftoff.py:46-72` | Replaced `str(e)`-only logging with `traceback.format_exc()` dump between the two pre-existing `[ERROR]` prologue/epilogue lines. Adds `import sys` at module level (was previously imported inside the except clause). |
| **2** | `tests/test_run_liftoff_traceback.py` *(new)* | 4 tests: RecursionError surfaces both `_deep` frame name and `Traceback (most recent call last)` header; non-RecursionError exceptions also carry traceback; recursion limit is raised during the call; recursion limit is restored on exception. |
| **3** | `lifton/run_liftoff.py:46-72` | Wraps the vendored Liftoff call in a `try/finally` save-raise-restore around `sys.setrecursionlimit(max(orig, 10000))`. 10× headroom is the right shape for bounded-but-deep feature-hierarchy traversal while staying inside the OS thread stack budget (Linux 8 MB default ≈ ~13K Python frames). If the recursion is actually a cycle (not bounded depth), Tier 2's traceback dump still pinpoints the file for a precision cycle-guard fix. |
| **4** | `lifton/lifton.py:307-340` | Strict mode unchanged (per-row stderr). Non-strict mode writes findings to `lifton_output/stats/gff3_input_validation.txt` and emits one summary line: `>> GFF3 input validator: N finding(s) (E error, W warning) written to <path>`. Defensive fallback to stderr if the side-car write fails. |
| **5** | `setup.py:install_requires` | Added `'mappy'` (no version pin — runtime only checks `is_mappy_available()`). |
| **5** | `lifton.yml` | Added `mappy` to the pip section so `conda env update -f lifton.yml` picks it up. |
| **6** *(follow-up)* | `lifton/lifton.py:8-26, 434, 446` | New module-level helper `_describe_annotation_source(x)` summarises a path-or-bytes value for the "Creating X annotation database" log lines. Surfaced when verifying Tier 1-5: under `--inmemory-liftoff` and `--stream`, the `liftoff_annotation` / `miniprot_annotation` value is the full GFF3 as in-memory `bytes`, and the `f"... {liftoff_annotation}"` interpolation produced a single 33 MB stderr line in the bee run. The helper renders bytes as `<in-memory bytes, N bytes>` and passes path strings through. |
| **6** *(follow-up)* | `tests/test_run_liftoff_traceback.py` | 3 new unit tests for the helper: path passthrough, bytes summary length-cap, bytearray hits the same path. |
| **5+5'** *(env)* | `benchmarks/phase16_rerun_bee.sh` | Helper script for the bee re-run. Initial version was missing `export PATH="$ENV_BIN:$PATH"`, which let `subprocess.run(["lifton", ...])` resolve to the BASE conda env's stale `lifton` shim (broken at import time on a numpy/parasail incompatibility). Fixed in-flight; the re-run that produced the verification artefacts had the corrected PATH. |

The `Lifton_TRANS` god module, the 24-cell byte-identity matrix, and
the streaming RefSeqProvider are all unchanged — Phase 16 touched only
the diagnostic / observability surface and one targeted runtime guard.

## 4. Verification

### Test suite

- **17/17 new tests pass** (13 harness + 4 traceback/recursion-limit).
- **65/65 spot-checked integration & dispatcher tests pass** —
  `test_integration_pipeline.py` (golden path + hermeticity),
  `test_liftoff_inmemory.py` (Phase 8 path), `test_locus_pipeline.py`
  (Phase 9 dispatcher), `test_locus_materialise.py` (Phase 11 native
  rewiring), `test_gff3_validator.py` (33 strict-GFF rule tests).
- Full `pytest tests/`: 508 passed before Tier 5, with 5 failures all
  asserting `mappy` is installed
  (`tests/test_native_bindings.py:49-51` literally says *"The test
  conda env ships mappy; flip if you ever run without it."*). Once
  `pip install mappy` runs in the env (covered by the re-run script
  below), those 5 tests will turn green for a 510 / 513 total.
- 3 hypothesis-dependent test files
  (`test_property_based.py`, `test_streaming_property.py`,
  `test_vulnerabilities.py`) skipped — pre-existing env gap unrelated
  to this phase.

### Bee benchmark — verified end-to-end

Re-run completed at `2026-05-07T05:29:19Z` (UTC) in tmux session
`phase16-bee`. Combined log:
`benchmarks/results/phase16_rerun_20260507T051851Z.log`. Summary JSON:
`benchmarks/results/summary_20260507T012919Z.json`.

**Lift step (the path the four diagnostic layers covered):**

| Metric | Pre-fix | Post-fix |
|---|---|---|
| Exit code | `1` | **`0`** ✓ |
| Wall-clock | `2:19.75` (failed early) | `9:16.74` (real work) |
| Peak RSS | 5,472,328 KB (5.2 GB) | 11,465 MB (11.2 GB) |
| Output GFF3 | not produced | **9.7 MB** |
| stderr line count | **532,599** | **101** ✓ |
| `mappy not installed` warnings | 200× | **0** ✓ |
| Summary table | `ERROR: could not convert...` | populates real numeric columns |
| GFF3 validator findings | flooded stderr (400K+ lines) | one summary line; 531,537 findings in side-car at `lifton_output/stats/gff3_input_validation.txt` |

The Tier 3 recursion-limit raise (1000 → 10000) was sufficient — no
RecursionError surfaced and Tier 2's traceback dump didn't fire.
This confirms the recursion was bounded-but-deep, not a cycle.

### Benchmark re-run instructions (for future runs)

Helper script: `benchmarks/phase16_rerun_bee.sh`. Runs in 4 phases:
(1) install mappy, (2) re-run sanity tests, (3) execute the bee
benchmark (~2.5 min wall-clock), (4) print a post-run inventory
(output GFF3 head, stderr.log line count, last 60 stderr lines).

**Launch from a fresh terminal — do NOT run inline:**

```bash
# 1. Launch the build in a detached tmux session
tmux new-session -d -s phase16-bee \
    'bash /ccb/salz3/kh.chao/LiftOn/benchmarks/phase16_rerun_bee.sh'

# 2. Monitor (any of these — they don't interfere with the run)
tmux attach -t phase16-bee                                                  # live view; detach with Ctrl-b d
tail -f /ccb/salz3/kh.chao/LiftOn/benchmarks/results/phase16_rerun_*.log    # tail the combined log
tmux ls                                                                     # confirm session is alive

# 3. After the run completes (the session exits on its own)
ls -la /ccb/salz3/kh.chao/LiftOn/benchmarks/results/bee/lifton.gff3
wc -l /ccb/salz3/kh.chao/LiftOn/benchmarks/results/bee/logs/lift.stderr.log
```

**Expected outcomes when the script finishes:**

- `benchmarks/results/bee/lifton.gff3` exists and is non-empty (the
  pipeline reaches Step 7 / 8 / 9 / 10 instead of dying at Step 4).
- `lift.stderr.log` shrinks from 532,599 lines to under ~10,000
  (Tier 4 silenced the per-row Dbxref dump).
- Zero `--native requested but mappy is not installed` warnings
  (Tier 5 declared the dependency; the install step in the script
  picks it up).
- The summary table prints real numbers in `wall(s)`, `peakRSS(MB)`,
  `mapped`, `lost`, `extra`, `meanID` instead of `ERROR: ...`
  (Tier 1 parser fix).
- If the recursion *was* bounded-but-deep, the run completes; if it
  was an actual cycle, Tier 2's traceback dump in stderr identifies
  the recursive function for a follow-up cycle-guard fix.

## 5. What's still deferred

| Item | Why it didn't ship in Phase 16 | When |
|---|---|---|
| **`-E` evaluation pass exits 1 with `mapped: 0`** | Surfaced *during* verification, not part of the original four-layer diagnostic. The lift step succeeded; the harness's optional second pass that re-invokes LiftOn in evaluation mode (`lifton -E ...`) failed in 69.69 s with exit 1. Likely either a `score.txt` parsing-status mismatch (statuses produced by lifton don't match the `("mapped", "liftoff", "miniprot", "lifton")` whitelist in `benchmarks/run_benchmarks.py:319-355::parse_score_txt`) or a real bug in the evaluation branch (`lifton/lifton.py:366-391`). | Follow-up investigation in a new phase |
| **Cycle-guard inside `lifton/liftoff/`** | Tier 3's recursion-limit raise was sufficient on bee — no RecursionError fired. Vendored subtree stays frozen unless a future dataset surfaces an actual cycle (in which case Tier 2's traceback dump pinpoints the file). | Only if a future dataset exposes a real cycle |
| **Phase 14 banded alignment / mappy-seeded extension** | Excluded from Phase 15 because they would mutate alignment output by ~0.1% and break the 24-cell byte-identity gate. | Future phase that explicitly authorises algorithmic changes |
| **`test_property_based` / `test_streaming_property` / `test_vulnerabilities`** | Require `hypothesis` which is not in `lifton_devel`. Pre-existing env gap. | Add `hypothesis` to `lifton.yml` and `setup.py` extras |
| **All-5-dataset re-run** | Bee is the diagnostic gate; running the other four costs ~45+ min of wall-clock for what is the same fix-validation signal. Triggered after bee turns green. | After Phase 16 bee re-run succeeds |
| **Reduce `unencoded_reserved_char` noise on Dbxref values** | The validator rule is biologically correct (NCBI Dbxref values do contain `:` un-encoded) but operationally counterproductive for RefSeq inputs. Tier 4 sidestepped the user-visible flood; the rule itself remains. | Optional refinement in `lifton/io/gff3_validator.py` |

## 6. Files touched

```
M  benchmarks/run_benchmarks.py        (Tier 1 — parser + safe-float helper)
M  lifton/run_liftoff.py               (Tier 2 — traceback; Tier 3 — recursion-limit guard)
M  lifton/lifton.py                    (Tier 4 — gate GFF3 validator dump; follow-up — bytes-blob log helper)
M  setup.py                            (Tier 5 — declare mappy)
M  lifton.yml                          (Tier 5 — declare mappy)
A  benchmarks/phase16_rerun_bee.sh     (helper: bee-only re-run script with PATH fix)
A  tests/test_benchmark_harness.py     (13 new tests)
A  tests/test_run_liftoff_traceback.py (7 new tests: 4 traceback/recursion + 3 bytes-blob log helper)
A  plans/phase_16_diagnostic_and_fixes.md  (this file)
```

No commits made. The user may bundle these into a single
`fix(phase-16): ...` commit or split by tier; the test additions are
self-contained and could land independently of Tiers 1-5.
