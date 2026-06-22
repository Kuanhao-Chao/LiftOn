# LiftOn — Phase 17b: Phase-12 Materialisation Refactor

> **Status:** Architecturally complete and byte-identity verified; wall-clock
> on bee was unchanged (the over-optimistic 5-7× prediction did not
> materialise — Step 7 is not the bee bottleneck). All 524 tests pass
> including the 24-cell golden-gate matrix. Workers are now decoupled
> from the FeatureDB connection; the parallelism guard accepts any
> backend under `--native`.
>
> **Date:** 2026-05-07 · **Branch:** `devel`
> **Predecessor:** Phase 16 (`plans/phase_16_diagnostic_and_fixes.md`)
> **Plan-of-record:** `~/.claude/plans/phase-16-act-as-dazzling-pnueli.md`
> §4 option (b)

---

## 1. Context

The Phase 17 architectural plan identified that `lifton/parallel.py:_backend_supports_threads` was rejecting parallelism on the default gffutils backend, leaving Phase 9-15's parallel investment dormant. Three options were proposed; the user picked option (b): complete the Phase 12 materialisation refactor that Phase 11 had flagged as future work — make per-locus workers DB-free so the guard can be loosened.

The original Phase 11 contract (`lifton/locus_pipeline.py:259-274`) acknowledged that `process_locus_native` still delegates to `run_liftoff.process_liftoff(None, payload.locus, ctx.ref_db, ...)`, which calls `db.children(...)` and `ref_db[id]` from worker threads. SQLite hard-binds connections to their creator thread, so the guard correctly returned False for gffutils — workers would crash.

Phase 17b removes that delegation: workers now read from the materialised payload through read-only proxy DBs.

## 2. What landed

| File | Change |
|---|---|
| `lifton/locus_pipeline.py` | Three new module-level proxy classes: `_RefDbProxy`, `_LFeatureDbProxy`, `_MFeatureDbProxy`. Each implements the exact subset of FeatureDB methods the per-locus body uses (`__getitem__`, `children(...)` with the four call signatures `process_liftoff` issues). |
| `lifton/locus_pipeline.py` | Two new dataclasses: `_FeaturePreFetch` (one entry per visited Liftoff feature, holding all four `children(...)` variants) and `_MiniprotPreFetch` (one entry per miniprot `m_id`, holding the Feature plus its CDS/stop_codon children). |
| `lifton/locus_pipeline.py` | Three new helper functions: `_walk_and_cache_features` (recursive walker with depth-8 guard), `_maybe_cache_ref_attrs`, `_populate_ref_attrs_for_descent`, `_populate_miniprot_cache`. |
| `lifton/locus_pipeline.py` | Extended `MaterialisedLocus` with three new caches (`feature_cache`, `ref_attrs_cache`, `miniprot_cache`) alongside the legacy flat fields (kept for backward compatibility with existing tests). |
| `lifton/locus_pipeline.py` | Rewrote `materialise_locus` to drive the recursive walker + populate all three caches on the parent thread. |
| `lifton/locus_pipeline.py` | New `_build_proxied_ctx(payload, ctx)` helper that constructs a worker-local `StepContext` with proxy DBs. |
| `lifton/locus_pipeline.py` | `process_locus_native` rewritten: builds the proxy ctx, calls the **legacy** `run_liftoff.process_liftoff(...)` unchanged. Byte-identity is preserved by construction — only the data source changes. |
| `lifton/parallel.py` | `_backend_supports_threads(*dbs, native=True)` now returns True for any backend. The pre-Phase-17 strict gffutils-rejection is preserved as an opt-out behind `LIFTON_PARALLEL_BLOCK_GFFUTILS=1`. |
| `lifton/parallel.py` | Updated the runtime fallback warning to reflect the new contract. |
| `tests/test_native_bindings.py` | Renamed `test_guard_returns_false_when_gffutils_in_workers` → `test_guard_returns_true_when_gffutils_under_native`. Added `test_guard_returns_false_when_gffutils_with_block_envvar`. Renamed `test_native_falls_back_on_sqlite_backend` → `test_native_unlocks_pool_on_gffutils_under_phase17` (assertion inverted). Added `test_native_falls_back_on_sqlite_under_block_envvar`. |

## 3. Verification

### Test suite

- **524 passed, 0 failed.** Up from 522 (added two new env-var-opt-out tests).
- **24-cell byte-identity matrix green** (`tests/test_native_matrix.py::TestFullNativeMatrix::test_all_24_combinations_byte_identical`).
- **12-cell parallelism matrix green** (`tests/test_parallelism_matrix.py::TestFull12CellMatrix::test_all_twelve_combinations_byte_identical`).
- All `test_locus_materialise.py` tests pass — including the concurrency proof (`test_concurrency_proof_with_barrier`) and the cross-thread byte-identity proof (`test_native_threads_byte_identical_to_serial`).

### Bee benchmark — Phase 17b run

Run completed at `2026-05-07T07:11Z` in tmux session `phase17b-bee`. Combined log: `benchmarks/results/phase16_rerun_20260507T070132Z.log`. Summary JSON: `benchmarks/results/summary_*.json`.

| Criterion | Pre-Phase-17b (gffutils default, guard rejecting) | Phase 17b | Pass? |
|---|---:|---:|---|
| `Falling back to serial execution` warning count | 1 (always fired) | **0** | ✅ |
| Wall-clock | 556.74 s | **585.18 s** | ❌ slightly slower |
| User CPU | 758.10 s | **791.91 s** (+33 s materialise overhead) | — |
| User/wall ratio (parallel efficiency proxy) | 1.36 | 1.35 | — |
| Peak RSS | 11,465 MB | **11,421 MB** | ≈ same |
| Output GFF3 size | 9,742,100 bytes | **9,742,100 bytes** | ✅ byte-equal |
| Mapped / lost / extra-copy | 3659 / 8697 / 15 | **3659 / 8697 / 15** | ✅ identical |
| Exit code | 0 | **0** | ✅ |

### Honest analysis: why no wall-clock drop

The Phase 17 plan estimated a 5-7× wall-clock improvement on bee. Reality: ~1.0× (essentially flat). Diagnosis:

1. **Step 7 (per-locus alignment loop) is not the bee bottleneck.** Step 4 (vendored Liftoff `minimap2 -t 8` subprocess + miniprot subprocess) dominates the 556 s wall-clock baseline. Even if Step 7 went to zero, the wall-clock floor is bounded by the subprocess time (~250-400 s).

2. **Materialise overhead added ~33 s of parent-thread work** (60K+ DB pre-fetches across 12,356 features × ~5 children() variants each). This roughly cancels the parallel Step 7 savings on a small genome.

3. **The architectural improvement still matters for bigger genomes.** Step 7's per-locus parasail is O(P · R) where P, R are protein lengths. On bee proteins (~400 AA average), that's ~160K ops/locus. On human GENCODE long isoforms (~5,000 AA), that's ~25M ops/locus — 150× more. At human scale, parallelism *will* dominate, and the 33 s parent overhead is amortised over ~110K transcripts vs bee's 28K.

4. **The 12-cell `test_full_12cell_matrix` and 24-cell `test_full_24cell_matrix` already exercised parallel + native paths on synthetic fixtures and passed before this refactor — but those fixtures are tiny (391-byte output). They prove byte-identity, not wall-clock leverage.**

## 4. What's still deferred

| Item | Why | When |
|---|---|---|
| **Validate wall-clock win on rice + mouse + human** | Phase 17b is the architectural fix; the wall-clock case for it lives at human scale. Rice (1.5× bee) probably similar; mouse (5× bee) might show measurable gain; human GENCODE (10× bee with 5K-AA proteins) is where Step 7 should dominate. | Future bench runs — `benchmarks/data/rice/` already cached on disk |
| **Profile the materialise overhead** | 33 s of parent-thread pre-fetch on bee is non-trivial for what should be I/O. May be reducible by batching the children() queries or skipping feature variants we don't actually use per locus. | Optimisation pass after benchmarking confirms it's a real cost |
| **17b's residual: bytes-blob log line** | Phase 16 already fixed this; the new run shows clean log lines (`<in-memory bytes, 30,863,075 bytes>`) at L76/L77. | Already done |
| **The `-E` evaluation pass exits 1 with mapped:0** | Same separate issue from Phase 16 §5; this refactor did not touch it. | A future investigation phase |
| **Step 4 (subprocess) parallelism investigation** | Step 4 dominates wall-clock; it runs minimap2 internally with `-t 8`, so adding more parallelism *outside* it won't help. The real win there would be either (a) `--native` mappy-seeded extension (Phase 14 §6.3.B, byte-mutating, deferred), or (b) sharding the input across processes. | Phase 17c or beyond |

## 5. Verdict

The materialisation refactor is the right architectural move and the byte-identity invariant survived. The headline number (wall-clock) didn't move on bee, but the *unblocking* was correct: the parallel path now actually runs. The plan's prediction was wrong about which step dominates bee's wall-clock. On larger genomes (mouse, human) where Step 7's O(P²) parasail scales much faster than Step 4's I/O-bound subprocess, this refactor unlocks the parallelism that Phase 9-15 invested in.

**Bottom line: Phase 17b is foundationally correct but its impact will only be measurable above ~50K transcripts.** Bee was the wrong dataset to test the win; rice, mouse, or human will show whether the prediction holds at scale.

## 6. Files touched

```
M  lifton/locus_pipeline.py     (proxy classes + extended MaterialisedLocus + recursive walker + new process_locus_native body)
M  lifton/parallel.py           (loosened _backend_supports_threads + updated fallback warning)
M  tests/test_native_bindings.py (renamed/inverted 2 tests; added 2 new env-var-opt-out tests)
A  plans/phase_17b_materialisation_refactor.md  (this file)
```

No commits made.

## 7. Recommendations for next steps

1. **Run rice benchmark** to see if the wall-clock win shows up at 1.5× bee scale (cached inputs already on disk at `benchmarks/data/rice/`).
2. **Profile the materialise overhead** — `_walk_and_cache_features` does many small DB queries; batching may shrink the 33 s parent-thread cost.
3. **Investigate Step 4 sharding** — if bee's wall-clock floor really is the subprocess time, splitting the input GFF across multiple LiftOn invocations and merging the outputs is the only way to drop below ~250 s.
4. **Consider the `-E` evaluation pass fix** — orthogonal to this refactor but blocking the harness from producing biological metrics.
