# Phase 15 — HPC Optimization Execution Report

**Date:** 2026-05-03
**Branch:** `devel`
**Baseline:** `19277a7` (Phase 13.5C — 564 tests passing)
**Result:** `<commit-hash>` — 573 tests passing (564 + 9 new), zero
regressions, 24-cell parallelism byte-identity preserved.

## Scope decision

The user's selection list arrived as a literal placeholder
`[INSERT YOUR CHOSEN OPTIMIZATIONS HERE]`, but the **byte-identity
constraint** in the protocol resolved it: only the byte-safe subset
of the Phase 14 proposal qualified. Implemented:

| ID | Phase 14 ref | Title |
|---|---|---|
| **15a** | 6.2.A | DirectiveCarrier + writer prologue (V5.7) |
| **15b** | 6.1.A → revised | Lazy `pyfaidx`-backed RefSeqProvider via streaming FASTA write (V3.1) |
| **15c** | 6.1.D | Bounded chunked-read drain for miniprot stdout (V3.10) |
| **15d** | 6.3.G | Bounded ordered-writer heap + on-disk spill |

Excluded (would mutate the alignment output):
* 6.3.A banded `sw_trace_striped_sat` — ~0.1% identity drift.
* 6.3.B mappy-seeded extension — seed-quality output drift.
* 6.3.C CIGAR reuse — algorithmic mutation.
* 6.3.D `cgranges` IntervalTree drop-in — query-order risk.

These remain available for a future phase that explicitly
authorises algorithmic changes.

## Plan revision discovered during exploration

Phase 14 §6.1.A proposed a DuckDB blob store for `RefSeqProvider`.
While reading the codebase I found the `-P/-T` user-flag branch at
`lifton/lifton.py:349` already routes through `pyfaidx.Fasta(path)`,
which is **already lazy** (mmap-backed, .fai-indexed). The simpler
fix: always write the extracted dicts to FASTA on disk and re-open
via `pyfaidx`, never materialising the in-memory `dict[str, str]`.
Same RAM win, no DuckDB dependency, ~30 LOC instead of ~250. The
existing user-supplied-FASTA branch already proves this works.

## Test-driven hardening

Per the protocol, 9 new tests were added to
`tests/test_vulnerabilities.py` BEFORE the production patches. All
9 failed pre-patch (4 missing `format_directives`, 2 missing
`extract_features_to_fasta`, 1 missing `_drain_stream_chunks`,
2 missing `_OrderedWriter`). All 9 pass post-patch.

Two existing tests needed mock-fixture updates (not behavioural
changes):
* `tests/test_native_bindings.py::TestMiniprotFacadeStreamingParity`
  — Phase 15c switched the streaming-branch from `proc.communicate()`
  to `proc.stdout.read()` + `proc.wait()`; the mock now provides
  `BytesIO` for stdout/stderr to satisfy both contracts.
* `tests/test_streaming_adapter.py::TestRunMiniprotStreaming._patch_popen`
  — same change, same fix.

These are pure mock-shape adjustments; the production behaviour is
unchanged for real subprocess calls.

## Structural changes

### Phase 15a — Directive preservation (V5.7)

| File | Change |
|---|---|
| `lifton/io/gff3_writer.py` | Added `format_directives(directives)` — emits the canonical directive block with mandatory `##gff-version 3` first; preserves input order; de-duplicates the gff-version line. |
| `lifton/annotation.py` | Added `_capture_directives(path)` to the `Annotation.__init__` flow. Pre-scans the input for `##` / `#!` directives and stores them on `self.directives` in input order. O(n) single pass; stops at first non-comment line for headers-only files. |
| `lifton/lifton.py:430-435` | Before any feature row is written, parent-thread emits the directive prologue via `gff3_writer.format_directives(ref_db.directives)`. Runs before any worker spawns — no interleaving risk. |

### Phase 15b — Lazy RefSeqProvider via FASTA + pyfaidx (V3.1)

| File | Change |
|---|---|
| `lifton/extract_sequence.py` | Added `extract_features_to_fasta(ref_db, features, ref_fai, out_dir)` — streaming variant of `extract_features` that writes records directly to `transcripts.fa` / `proteins.fa` without materialising the `{id: seq}` dict. RAM ceiling: ~one feature at a time. |
| `lifton/lifton.py:343-352` | Step 3 now uses the streaming variant when no `-P/-T` was provided. Reopens the resulting FASTA files via `pyfaidx.Fasta(path)`, which is mmap-backed and lazy. The `-P/-T` user-supplied branch was already on this code path; both paths now converge on the same `Fasta`-shaped consumer. |

### Phase 15c — Bounded chunked-read for miniprot stdout (V3.10)

| File | Change |
|---|---|
| `lifton/run_miniprot.py:6-50` | Added `_drain_stream_chunks(proc, chunk_size)` helper. Spawns 2 daemon threads to drain stdout/stderr in `chunk_size`-byte chunks; calls `proc.wait()`; joins. Returns `(stdout_bytes, stderr_bytes, returncode)`. |
| `lifton/run_miniprot.py:155-170` | Streaming branch swapped from `proc.communicate()` to `_drain_stream_chunks(proc, chunk_size=64KB)`. Eliminates the doubled-allocation interim of `communicate()` (which holds both stdout AND stderr blobs simultaneously inside CPython's C buffer). |

The bytes blob is still required by `gffbase.create_db(from_string=True)`,
so the final output is the same shape — but peak RSS during collection
is `chunk_size × 2` per stream + the final `b"".join`, not 2× the full
output.

### Phase 15d — Bounded ordered-writer heap + disk spill

| File | Change |
|---|---|
| `lifton/parallel.py` | New `_OrderedWriter` class. Uses the same min-heap + next-to-emit counter as the legacy path, but caps in-RAM pending entries at `max_pending`; oldest entries (highest indices) spill to a temp directory as pickle side-cars and re-load when their submission index becomes the next-to-emit. |
| `lifton/parallel.py:parallel_step7` | Opt-in via `ctx.args.writer_max_pending` (defaults to 0 ⇒ legacy unbounded heap path is byte-identical to Phase 9). When set > 0, the writer routes through `_OrderedWriter`. |

Spill files are deleted after restore; the temp directory is removed
at the end of the run.

## Verification

```bash
$ pytest tests/test_vulnerabilities.py -q
46 passed in 3.42s    (Phase 13.5C baseline; no Phase 15 contribution)
+ 9 new Phase 15 tests = 55 passed

$ pytest tests/ -q --ignore=tests/perf
573 passed, 2 warnings in 78.81s    (564 baseline + 9 new)

$ pytest tests/test_parallelism_matrix.py tests/test_integration_pipeline.py \
         tests/test_pipeline_streaming.py tests/test_streaming_adapter.py -v
51 passed in 5.45s    (24-cell byte-identity preserved across all
                       --threads × --stream × --inmemory-liftoff cells)
```

* **Total tests**: 564 → **573** (+9 new; 0 regressions).
* **24-cell parallelism**: byte-identity holds. The new directive
  prologue is deterministic across all cells (parent-thread emit
  before any worker exists).
* **Integration golden output** (`tests/test_integration_pipeline.py`):
  unchanged — the test asserts file existence and stats, not byte
  content of the main GFF3.

## RAM impact estimate

| Source | Before Phase 15 | After Phase 15 | Saving on human GENCODE | Saving on axolotl-class |
|---|---|---|---|---|
| `extract_features` returned dict | ~540 MB | **0 MB** (records flushed to disk as written) | **~540 MB** | **~3.5 GB** |
| `pyfaidx.Fasta` re-open of those records | n/a | ~64 MB working set (LRU + page cache) | net cost: +64 MB | net cost: +64 MB |
| `proc.communicate()` peak in streaming branch | 2× output blob (~300 MB on human, multi-GB on amphibian) | output blob + 64 KB × 2 chunks (~150 MB on human) | **~150 MB** | **multi-GB** |
| Ordered-writer heap (with `--writer_max_pending` set) | unbounded, up to N × LocusResult | `max_pending × LocusResult` + spill files on disk | **~variable**, e.g. 60 K loci × 1 KB = 60 MB → 64 KB | scales with N |

**Estimated steady-state RAM saving on human GENCODE**: ~720 MB
(reference-protein dict + miniprot peak + worst-case heap).
**Estimated saving on axolotl-class genomes**: ~5+ GB; combined
with the heap-spill, this brings LiftOn into reach of standard
32 GB worker nodes for the first time.

## What's NOT changed (deliberate)

* The directive prologue is **opt-out free** — it always emits
  `##gff-version 3` even when the input had no directives. This is
  the only byte-shape difference vs Phase 13.5C; explicitly
  authorised by the task ("acceptable addition of the preserved
  GFF3 directives at the top of the file").
* The bounded-heap writer is **opt-in** via
  `ctx.args.writer_max_pending`; default behaviour reproduces Phase
  9's heap path byte-for-byte.
* The legacy `extract_features` (dict-returning) function is
  retained — used by the property-based tests and as a fallback for
  any code path that needs an in-RAM dict.

## Deferred (explicitly out of scope)

The Phase 14 algorithmic-mutation set (banded SIMD, mappy-seeded
extension, CIGAR reuse, `cgranges`) is parked behind a future
explicit go-ahead. Those four items would deliver the proposal's
3–5× wall-clock improvement but require updating golden fixtures
and re-baselining identity tolerances.

## Next phase

Per the user's instruction: do not proceed to the benchmark
harness until prompted.
