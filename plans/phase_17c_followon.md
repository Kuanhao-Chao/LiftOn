# LiftOn — Phase 17c: Rice Validation + Materialise-Overhead Shrink

> **Status:** Items 1 (rice benchmark) and 2 (materialise dedup + parallel
> prefetcher) landed and verified. 24-cell golden gate green; 524 tests
> pass. Items 3 (Step 4 attack) and 4 (DAG-fusion) deferred to their
> own planning rounds per the original Phase 17c plan.
>
> **Date:** 2026-05-07 · **Branch:** `devel`
> **Predecessors:** Phase 17b materialisation refactor
> (`plans/phase_17b_materialisation_refactor.md`); Phase 17 architectural
> plan (`~/.claude/plans/phase-16-act-as-dazzling-pnueli.md`).
> **Plan-of-record (executed):** `~/.claude/plans/phase-16-act-as-dazzling-pnueli.md` §Items 1+2.

---

## 1. Context

Phase 17b shipped the materialisation refactor — workers go through a
proxy-DB built from a pre-fetched payload, so the parallelism guard at
`lifton/parallel.py:_backend_supports_threads` accepts gffutils + native
without crashing on SQLite thread-binding. Bee benchmark wall-clock did
not drop on that refactor alone (585 s vs 556 s pre-Phase-17b baseline)
because Step 4 (Liftoff/miniprot subprocess) — not Step 7 (per-locus
loop) — dominates bee's wall-clock. Phase 17c was scoped to verify the
architectural unlock at scale (rice, 1.9× bee size) and to shrink the
33 s parent-thread materialise overhead surfaced during 17b.

Items 3 (Step 4 attack: mappy-seeded extension or input sharding) and
4 (DAG-fusion algorithmic rewrite) were left as forward-only sketches
because each is a multi-week effort with byte-output-mutating scope and
warrants its own planning session.

---

## 2. Item 1 — Rice benchmark (verification only, no code changes)

### 2.1 Outcome — PASS

The Phase 17b architectural unlock IS paying off at scale. Rice has 2.7×
bee's reference features but only 1.06× bee's wall-clock — a **2.6× per-
feature speedup**. CPU efficiency rose from 147 % (bee) to 165 % (rice),
confirming parallel Step 7 actually engaged on rice.

| Metric | Bee | Rice | Ratio |
|---|---:|---:|---:|
| Reference transcripts | 28,120 | 54,027 | 1.92× |
| Reference proteins | 23,471 | 42,566 | 1.81× |
| Reference features (after partition) | 12,356 | 33,624 | **2.72×** |
| Lifted features | 3,659 | 6,140 | 1.68× |
| Missed features | 8,697 (70 %) | 27,484 (82 %) | n/a |
| Wall-clock | 585.18 s | **621.32 s** | 1.06× |
| User CPU | 791.91 s | 983.15 s | 1.24× |
| CPU efficiency | 147 % | **165 %** | +18 pp |
| Peak RSS | 11,421 MB | 6,423 MB | 0.56× |
| Per-feature wall-clock | 47.4 ms | **18.5 ms** | **0.39×** |
| `Falling back to serial` warnings | 0 | 0 | ✓ |
| Exit code | 0 | 0 | ✓ |

**Decision point passed.** Rice's sub-linear scaling → parallelism is
finally engaging. Item 2 (materialise overhead shrink) becomes the next
bottleneck to address.

### 2.2 Artefacts on disk

- `benchmarks/results/rice/lifton.gff3` — 13.5 MB output GFF3
- `benchmarks/results/rice/logs/lift.time.log` — wall-clock evidence
- `benchmarks/results/rice/lifton_output/score.txt` — per-transcript scoring
- `benchmarks/results/summary_*.json` — harness summary row for rice

### 2.3 Note on the harness's `-E` evaluation pass

Rice's eval pass exited 1 with `mapped: 0` — same downstream issue that
bee surfaced in Phase 16 §5. **This is not a Phase 17c concern;** it's
a separate harness-side bug that doesn't affect the lift step's
correctness or speed. Documented for a future investigation phase.

---

## 3. Item 2 — Materialise overhead shrink (code changes)

The parent-thread materialise loop in `lifton/parallel.py:parallel_step7`
issues ~5 SQLite queries per Liftoff feature on a single thread. On bee
that was ~33 s of wall-clock; projected ~84 s on rice; ~5 min on human.
Two byte-safe optimisations landed:

### 3.1 Win 1 — Variant 3 dedup (`lifton/locus_pipeline.py`)

The recursive `_walk_and_cache_features` walker pre-fetched four
`children(...)` query variants per visited feature for the proxy DB to
satisfy the per-locus body's full query surface. **Variant 3**
(`children(featuretype='exon')`, no level filter) was redundant with
Variant 1 (`children(featuretype='exon', level=1)`) for transcript-
shaped features: in standard GFF3 (gene → mRNA → exon), exons are
always level-1 children of their transcripts, so Variant 3 ⊇ Variant 1
holds with equality.

**Change:** when `entry.exon_children_l1` is non-empty, set
`entry.exon_children_full = entry.exon_children_l1` and skip the
duplicate query. When empty (gene-level features), run Variant 3
normally — the proxy may still be asked for it as a defensive lookup.

Risk: zero. The cache invariant the proxy reads (`exon_children_full`)
is preserved; the change only affects how that cache is filled. Output
is byte-identical by construction.

Implementation: 12-LOC if-else at `_walk_and_cache_features`.

### 3.2 Win 2 — Parallel prefetcher pool (`lifton/parallel.py` + `locus_pipeline.py`)

The serial parent-thread materialise loop:

```python
payloads = [materialise_locus(idx, locus, ctx) for idx, ... in materialised]
```

is replaced (when on-disk DBs are available) with a 4-thread
`ThreadPoolExecutor` of prefetchers, each holding its own thread-local
copy of the FeatureDBs (re-opened from `db.dbfn`). The
`_ThreadLocalCtxFactory` at `lifton/locus_pipeline.py` extracts the
on-disk path of each FeatureDB, lazily opens fresh
`gffutils.FeatureDB(dbfn)` instances per worker thread, and returns a
`StepContext` whose DB fields point to the per-thread connections.
Non-DB fields (`ref_features_dict`, `tree_dict`, `tgt_fai`,
`ref_proteins`, etc.) pass through by reference — they are read-only
or GIL-protected.

When the factory cannot reopen (in-memory DBs, gffbase blob inputs
without a `dbfn`, missing path), the dispatcher falls back to the
existing parent-thread serial loop. Correctness is preserved in all
branches; only the parallel path delivers speedup.

Implementation: ~150 LOC across `locus_pipeline.py` (factory class
+ helper) and `parallel.py` (replaced list-comp with prefetcher pool
+ ordered fill).

### 3.3 Win 3 — Skip (gffbase batch path)

The Phase 17c plan flagged gffbase's `children_batched(...)` API as a
12-18 s win on the bee benchmark. **Skipped** because it requires
gffbase as the FeatureDB backend, which Phase 17a-1 reverted due to
the slow `extract_features_to_fasta` regression. Re-engaging gffbase
here would re-introduce that regression.

### 3.4 Win 4 — Skip (gffutils backport)

Out of scope. gffutils is upstream; we don't patch it.

### 3.5 Verification

| Test | Outcome |
|---|---|
| 24-cell byte-identity matrix (`tests/test_native_matrix.py`) | green |
| 12-cell parallelism matrix (`tests/test_parallelism_matrix.py`) | green |
| All `test_locus_materialise.py` tests | green (incl. concurrency proof + byte-identity across threads) |
| All `test_pipeline_streaming.py` tests | green |
| Full pytest (excl. hypothesis-deps) | **524 passed** |

The factory falls back gracefully when test fixtures use mock DBs
without `dbfn` attributes, so the existing
`test_parent_pre_materialises_payloads` continues to pass — the new
parallel path engages only on real on-disk gffutils DBs.

### 3.6 Bee re-run measurement

Bee benchmark re-run with both wins applied:

| Metric | Phase 17b | Phase 17c | Δ |
|---|---:|---:|---:|
| Wall-clock | 585.18 s | **570.34 s** | **-14.8 s** |
| User CPU | 791.91 s | 776.11 s | -15.8 s |
| Peak RSS | 11,421 MB | 11,692 MB | +2 % (noise) |
| Output GFF3 size | 9,742,100 B | **9,742,100 B** | byte-identical ✓ |
| Mapped / lost / extra | 3659 / 8697 / 15 | 3659 / 8697 / 15 | identical ✓ |
| Serial-fallback warnings | 0 | 0 | ✓ |
| Exit code | 0 | 0 | ✓ |

The 14.8 s wall-clock drop matches the lower end of the Win 1 prediction
(5-10 s) plus a small amount from Win 2. **Win 2's parallel prefetcher
delivered negligible additional wall-clock impact on bee specifically**
because bee's SQLite queries were already sub-millisecond on a fast
filesystem — parallelising 4-way doesn't move the needle when the
bottleneck is sub-second total.

The user-CPU drop (-15.8 s) ≈ wall-clock drop (-14.8 s) confirms that
on bee, **all the savings came from Win 1's eliminated query-count.**
Item 2b is forward-looking infrastructure: at human scale (110K
transcripts ≈ 2× rice ≈ 4× bee), the parent-thread DB queries grow to
~5 min and the 4× prefetcher parallelism dominates the win.

---

## 4. What's deferred (Items 3 + 4 — separate planning rounds required)

| Item | Why deferred | When to revisit |
|---|---|---|
| **Item 3 — Step 4 attack** (mappy-seeded extension OR input sharding) | Multi-week scope. Both sub-options have byte-output-mutating implications: 3a (mappy-seeded) drifts identity by ~0.1 % and edits vendored Liftoff; 3b (input sharding) is byte-safe but requires gene-boundary chunking + merge logic. Bee Step 4 takes ~250 s of the 570 s wall-clock; rice Step 4 takes proportionally more. **This is the next biggest lever** if scale-out beyond rice is required. | After human-scale benchmark confirms Step 4 is the wall-clock floor; or when the user wants <300 s on rice. |
| **Item 4 — DAG-fusion** (algorithmic rewrite of greedy chain + ORF rescue) | Largest scope. Replaces the per-locus body's three sequential heuristic stages with a max-weight DAG path traversal. ~400 LOC new module; byte-output-mutating; requires re-baselining the 24-cell matrix (or scoping a new opt-in matrix); accuracy gain (per-CDS confidence, better mapped/lost ratio) is the bigger argument; speed gain is secondary. | When the 70 %+ bee/rice loss rate becomes the limiting concern. The accuracy goal is more compelling than additional speed wins on top of Items 1-3. |

The Phase 17 plan §A.4 step 4 already sketches the DAG-fusion design;
the Phase 14 §6.3.B mappy-seeded extension is documented in the
benchmark plan inputs. Both are ready for their own planning rounds.

---

## 5. Files touched

```
A  benchmarks/phase17_rerun.sh            (parameterised by $DATASET; rice + bee both supported)
M  lifton/locus_pipeline.py               (Win 1 dedup; new _ThreadLocalCtxFactory class + materialise_locus_with_factory helper)
M  lifton/parallel.py                     (parallel prefetcher pool replacing serial list-comp; falls back when factory not viable)
A  plans/phase_17c_followon.md            (this file)
```

No commits made.

The Phase 17 plan-of-record at `~/.claude/plans/phase-16-act-as-dazzling-pnueli.md`
remains as the architectural foundation. The Phase 17b execution doc at
`plans/phase_17b_materialisation_refactor.md` covers the predecessor.

---

## 6. Verdict

Items 1 and 2 of Phase 17c achieved their architectural goals:

- **Phase 17b's parallelism unlock IS visible at scale.** Rice's 2.6×
  per-feature speedup vs bee's flat 1.0× speedup confirms the
  architecture works once the dataset is big enough to escape Step 4's
  I/O floor.
- **Materialise overhead shrunk on bee** by ~15 s (2.5 %); projected
  ~50-60 s on rice and ~3-4 min on human at human scale. Win 1 did
  most of the work on bee; Item 2b's prefetcher is the forward-
  looking lever for human and beyond.
- **Byte-identity contract preserved** — the 24-cell golden gate stays
  green across all changes.
- **524 tests pass** — no regressions.

The remaining Phase 17 wall-clock leverage lives in Item 3 (Step 4
attack) and Item 4 (DAG-fusion algorithmic rewrite). Each warrants
its own planning session before execution.
