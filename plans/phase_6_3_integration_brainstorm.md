# LiftOn — Phase 6.3: Deep Integration Brainstorm

> **Mandate:** brainstorm 3-4 architectures that bypass the
> sequential disk + DB round-trip in the current pipeline. No code is
> written this phase; the user picks the winning architecture from
> the ranked list below.

---

## 1. Today's pipeline — where the time actually goes

Tracing `run_all_lifton_steps` (`lifton/lifton.py:208-525`) for a
typical chr22 run, the I/O cliff is concentrated in five places:

| Stage | Lines | I/O cost (qualitative) | Disk artifacts produced |
|---|---|---|---|
| **S1** Ref FASTA + ref GFF DB build | `lifton.py:233-238` | medium (gffutils SQLite build) | `<ref.gff3>_db` SQLite |
| **S2** Eager `extract_features` of all transcripts/proteins | `lifton.py:255` → `extract_sequence.py:71-94` | medium (full ref-protein dict in RAM) | `intermediate_files/{proteins,transcripts}.fa` |
| **S3a** Vendored Liftoff → minimap2 SAM → re-parse → write GFF3 | `run_liftoff.py:33-58` → `liftoff/align_features.py:53-66` → `liftoff/write_new_gff.py:13-28` | **high** (minimap2 writes per-chr SAMs to `args.directory`; pysam re-parses; write_new_gff dumps a fresh GFF3) | per-chr SAM files; `liftoff_output/liftoff/liftoff.gff3` |
| **S3b** miniprot subprocess → GFF3 file | `run_miniprot.py:62-67` | **high** (full miniprot stdout captured to disk) | `liftoff_output/miniprot/miniprot.gff3` |
| **S4** Re-ingest both GFFs as gffutils DBs | `lifton.py:362, 375` | **high** (two more gffutils SQLite builds) | `liftoff.gff3_db`, `miniprot.gff3_db` SQLite caches |
| **S5** Sequential walk over Liftoff genes + miniprot mRNAs | `lifton.py:418-444` | medium (per-feature SQLite seeks) | output `lifton.gff3` |

The **double materialisation** is the architectural smell:
- Liftoff hands LiftOn a finished `lifton_output/liftoff/liftoff.gff3`
  on disk. LiftOn then **re-parses that file** into a SECOND gffutils
  DB it didn't need to write in the first place.
- Same pattern with miniprot.
- Even within Liftoff, each `align_single_chroms` writes per-chr SAMs
  to `args.directory` and the parent re-opens them with pysam.

Two-pass disk round-trips dominate wall time on real annotations. The
Phase 6.1 audit (and gffbase's own MIGRATION.md) measured the gffutils
SQLite ingest at O(seconds) on a chr22 GFF, and **two of these run per
LiftOn invocation** in S4 alone.

---

## 2. The four architectural strategies, ranked by `(performance × confidence) / risk`

| Rank | Strategy | Expected speedup vs today | Implementation effort | Schedule risk | Net score |
|---|---|---|---|---|---|
| **#1** | **Streaming Adapter Layer** — keep the external `minimap2` / `miniprot` binaries; capture their stdout, push records straight into in-memory gffbase / Python objects without writing intermediate GFF3 files | **3-5×** end-to-end on chr22; ~10-15× on the S4 re-ingest step | Low-medium | **Low** | **Highest** |
| **#2** | **Vendored-Liftoff in-memory refactor** — change Liftoff's public entrypoint to return generators of `Lifton_GENE`-shaped objects directly, eliminating S3a's `write_new_gff` + S4 re-ingest | 2-3× on top of #1 | Medium | Medium | High |
| **#3** | **Algorithmic Fusion (per-locus pipeline)** — interleave Liftoff alignment, miniprot alignment, parasail, protein-maximisation chaining, and ORF rescue PER LOCUS as a single coroutine; eliminate the global GFF3 materialisation entirely; output emitted incrementally as loci complete | 5-8× wall-clock on a 64-core box; major memory ceiling drop | High | High | Medium |
| **#4** | **Native-binding alignment (mappy + pyminiprot)** — drop the `subprocess` calls altogether, use `mappy` (the official Python bindings to minimap2) and a `pyminiprot` equivalent in-process. Combines naturally with #3 | 1.5-2× on top of #3; eliminates fork+exec overhead | Very high (pyminiprot bindings don't exist yet — would need to build them) | Very high (depends on a binding we don't yet have) | Lowest |

The four strategies stack: #1 is independently shippable, #2 builds on
#1, #3 builds on #2, #4 is the long-term endgame. Each phase below can
land independently while keeping the test suite green.

---

## 3. Strategy #1 — Streaming Adapter Layer (recommended first)

### What it changes
- `lifton/run_miniprot.py:62-67` — replace `subprocess.run(...,
  stdout=fw)` (writes `miniprot.gff3`) with `subprocess.Popen(...,
  stdout=PIPE)` and stream the stdout pipe **directly** into
  `gffbase.parse_bytes(...)` (already exists at
  `lifton/gffbase/parser.py:127`) or into an in-memory
  `lifton.gffbase.create_db(..., from_string=True)` (line 39 of
  `lifton/gffbase/create_db.py` already supports `from_string`).
- `lifton/run_liftoff.py:33-58` — keep the vendored Liftoff body as-is
  but tee its `lifted_feature_list` (already an in-memory dict in
  Liftoff!) into an in-memory gffbase `FeatureDB(":memory:")` instead
  of writing `liftoff.gff3` then re-parsing it in S4.
- `lifton/lifton.py:362, 375` — instead of `Annotation(<gff3 path>)`,
  accept a live `gffbase.FeatureDB` constructed in-memory.

### Why it's #1
- **Zero algorithmic risk.** The external aligners produce identical
  bytes; we just read them from a pipe instead of a file.
- **No vendored-Liftoff surgery.** Strategy #1 *consumes* Liftoff's
  existing in-memory `lifted_feature_list` instead of refactoring it.
- **Maximum return on minimum code surface.** Eliminates two of the
  three SQLite materialisations (S3 output + S4 re-ingest), which is
  60-70 % of today's I/O cliff per the Phase 6.1 numbers.
- **Phase 6.2 adapter pattern reusable.** The `gffbase_adapter`
  module designed in `phase_6_2_migration_plan.md` is the exact
  surface where the streaming variant lands.

### Concrete steps (no code in this phase — these are the implementation hooks)
1. **gffbase**: confirm `create_db(data, dbfn=':memory:',
   from_string=True)` works on a multi-megabyte string. Verified
   against `python/gffbase/create_db.py:53-62`: yes, materialises to
   a tempfile, ingests, deletes.
2. **gffbase improvement (out of scope for this phase)**: file a
   feature request to make `from_string` accept an iterator-of-bytes
   so we can avoid even the temp file. Without it, #1 still wins by
   removing S4's *read*, just not the *write*.
3. **LiftOn `run_miniprot`**: switch to `Popen` with
   `stdout=PIPE`; pass the bytes object straight to
   `gffbase.create_db(..., from_string=True)`. Total LOC change ≈ 30.
4. **LiftOn `run_liftoff`**: tap into Liftoff's
   `lifted_feature_list` (returned by `liftoff_main.run_all_liftoff_steps`
   as the in-memory dict; written to GFF3 by `write_new_gff` only as
   a side effect at lines 36, 42 of `liftoff_main.py`). Skip the
   `write_new_gff` call when LiftOn is the consumer; instead
   serialise it into an in-memory `FeatureDB` directly.
5. **LiftOn `Annotation`**: extend the Phase 6.2 backend resolver to
   accept a pre-built `FeatureDB` object as a constructor argument
   so `lifton.py:362, 375` can pass it through without
   re-instantiating from a file path.

### Estimated speedup
- S3b (miniprot write): saved entirely → **~5-15s on chr22**.
- S4 (re-ingest both DBs): saved → **~10-30s on chr22**.
- S3a (Liftoff write): saved when Strategy #2 lands; until then we
  still pay the GFF3 write but skip the re-ingest.
- Net: **3-5× wall-clock reduction** end-to-end on chr22.

### Risks
- `subprocess.Popen` deadlock if we don't drain stderr in parallel.
  Use `communicate()` or threaded stderr drainer.
- Memory pressure: full miniprot output may be 100-500 MB on a
  mammalian genome. Acceptable; gffbase ingest is streaming
  internally so peak RSS scales O(records) not O(file size).

### Test gate (Phase 5 baseline preserved)
- All 297 tests green with `--backend=gffbase --stream` flag off
  (default).
- New `tests/test_streaming_pipeline.py` parametrises over
  `stream={off,on}` and asserts byte-identical output GFF3.

---

## 4. Strategy #2 — Vendored-Liftoff in-memory refactor

### What it changes
- Add a new public entrypoint in `lifton/liftoff/liftoff_main.py`:
  `run_all_liftoff_steps_inmemory(args, ref_db) -> Iterator[LiftedFeature]`
  that yields the items in `lifted_feature_list` instead of routing
  them through `write_new_gff.write_new_gff`.
- Replace `lifton/run_liftoff.py:33-58` to consume the iterator,
  pushing each lifted feature directly into LiftOn's in-memory
  `gffbase.FeatureDB(':memory:')`.
- The vendored Liftoff continues to exist on disk and gets the same
  three retry strategies as before.

### Why it's #2
- Builds on Strategy #1's adapter layer; can land as a follow-up PR
  without re-doing #1 work.
- Removes the **last** GFF3 materialisation (Liftoff's own output).
- Makes the vendored Liftoff an in-process library, finally
  collapsing the "fork → child writes file → parent reads file" anti-
  pattern that Phase 1 audit flagged.

### Concrete hooks
1. `liftoff_main.py:42` — `write_new_gff.write_new_gff(...)` is the
   sole serialisation point. Already a no-op for the in-memory caller
   if we expose `lifted_feature_list` directly to `run_liftoff.py`.
2. `liftoff/lift_features.py` and `liftoff/find_best_mapping.py` —
   already work with in-memory `aligned_seg`/`new_feature` objects;
   no change needed.
3. The retry-on-failure logic (Liftoff has none today) is inherited
   from Strategy #1's adapter.

### Risks
- **Vendored fork drift.** Adding a new entrypoint means the next
  upstream Liftoff merge will need conflict resolution. Acceptable
  cost for the integration win; document the entrypoint in the
  vendored README.
- **Memory.** `lifted_feature_list` is a dict keyed by gene id with
  the entire chr22 worth of lifted features in RAM. Already fits
  today (Liftoff itself does this); we're just plumbing it forward.

### Test gate
- 297 + new tests still green.
- Wall-clock benchmark: Strategy #2 should remove an additional
  **5-10s** on top of Strategy #1's chr22 time.

---

## 5. Strategy #3 — Algorithmic Fusion (per-locus pipeline)

### What it changes
The current pipeline is GFF-major: parse the entire reference GFF →
align everything → re-parse → walk genes. Strategy #3 inverts to
**locus-major**: for each gene, run the small slice of work that gene
needs, then emit and discard.

```
for gene in stream_reference_genes(ref_gff):           # generator
    minimap2_hits   = align_dna_chunk(gene)            # mappy or popen-of-1
    miniprot_hits   = align_protein(gene.protein)      # ditto
    lifton_gene     = lift(gene, minimap2_hits)
    miniprot_trans  = lift_protein(gene, miniprot_hits)
    if both:
        chained = protein_maximization.chaining_algorithm(...)
    if mutations:
        rescue = orf.rescue.find_best_orf(...)
    write_entry(lifton_gene, fw)                       # emit immediately
    del lifton_gene                                    # free RAM
```

### Why it's #3 (not #1)
- **Highest theoretical speedup.** No global materialisation, ever.
  Memory ceiling drops to O(single-locus). Pipeline becomes trivially
  parallel at the locus level — combines naturally with the Phase 4
  Step 9 `ProcessPoolExecutor` plan (`phase_4_roadmap.md` §9).
- **But it requires a deep refactor.** Rewires
  `run_all_lifton_steps` end-to-end. Touches every module that
  currently assumes "the reference DB is already built."
- **Algorithmic correctness risk.** The current global-then-iterate
  model has subtle dependencies — e.g.,
  `lifton_utils.miniprot_id_mapping` (`lifton_utils.py:431-460`)
  builds a corpus-wide ref↔miniprot map before any per-gene work
  starts. Going locus-major requires re-deriving those globals
  per-locus or building them lazily.

### Concrete hooks
1. New module `lifton/locus_pipeline.py` that owns the per-gene
   coroutine.
2. Move `protein_maximization.chaining_algorithm` and
   `Lifton_TRANS.orf_search_protein` so they can be invoked without a
   pre-built `ref_proteins` dict — pass the single ref protein the
   locus needs.
3. Use a small in-memory IntervalTree of "alignments seen so far" to
   detect overlapping miniprot mRNAs from neighbouring genes (the
   only cross-gene state that survives in the locus-major model).
4. Existing `tree_dict` (`intervals.py:6-13`) becomes a streaming
   per-chromosome tree built incrementally as Liftoff genes are
   lifted.

### Risks
- **Determinism.** Today's output gene ordering is deterministic via
  the Liftoff feature order. With locus-major + parallel workers, we
  need an ordered-writer buffer (the same one Phase 4 §9 designed for
  parallelism).
- **Cross-gene merging.** Step 8 of the current pipeline
  (`lifton.py:435-441`) handles miniprot mRNAs that don't overlap any
  Liftoff gene. In locus-major mode this becomes a "deferred" pass at
  the end; or we treat miniprot-only mRNAs as their own pseudo-loci.
- **Vendored-Liftoff incompatibility.** Liftoff today is GFF-major:
  it builds the global `lifted_feature_list` before returning. Going
  locus-major **for free** requires Liftoff to expose a per-locus
  yield, which it doesn't. Workarounds: (a) keep Liftoff GFF-major
  and only fuse Steps 7+8; (b) refactor Liftoff's
  `lift_original_annotation` to yield. (a) is much safer for Phase 6.

### Test gate
- All 297 tests still green; integration test asserts byte-identical
  output GFF3.
- New benchmark: peak RSS drops by ≥ 50 % on chr22 vs Strategy #2.

---

## 6. Strategy #4 — Native-binding alignment (mappy + pyminiprot)

### What it changes
- Replace `subprocess` calls to `minimap2` with the `mappy` Python
  module (official; ships with minimap2 source). `mappy.Aligner(idx)`
  + `aligner.map(seq)` returns hits as Python objects with no
  fork/exec/SAM-write overhead.
- Replace `subprocess` calls to `miniprot` with… **whatever does not
  exist today.** miniprot has no `pymiriprot` binding in pip; the C
  source under `lh3/miniprot` exposes a `mp_align` API but no
  PyO3/pybind11 wrapper.

### Why it's #4 (lowest practical priority)
- **`mappy` ready today.** Drop-in for `align_features.py` once
  Strategy #2 has eliminated the SAM-file round-trip.
- **`pyminiprot` does NOT exist.** We would need to build it
  ourselves — about 1-2 weeks of PyO3 / pybind11 work to wrap
  miniprot's `mp_align` and yield Python tuples. This is real
  engineering, not glue.
- **Marginal additional speedup.** Subprocess fork overhead is
  ~30-50 ms per invocation; minimap2's actual alignment is the
  dominant cost. Switching to `mappy` saves the fork tax once per
  chromosome batch, which is ~seconds, not minutes.
- **Useful only AFTER Strategies #1-#3 land.** Without locus-major
  fusion, the per-call alignment volume isn't high enough to make
  the binding overhead worth it.

### Concrete hooks (deferred)
1. Phase 6.4 (or later): write a `pyminiprot` PyO3 binding upstream
   in the miniprot repo (or a sidecar repo), publish to PyPI.
2. Phase 6.5 (or later): swap `subprocess` for `mappy` /
   `pyminiprot` in `align_features.py` and `run_miniprot.py`.

### Risks
- New binding maintenance burden — adds a Rust/PyO3 dependency
  separate from the existing gffbase Rust extension.
- Version skew: every new minimap2 / miniprot release would need a
  binding rebuild.

### Test gate
- Same as Strategy #1-#3; byte-identical output.

---

## 7. Recommended schedule

| Phase | Strategy | Wall-clock target | Memory target |
|---|---|---|---|
| **6.3.1** | Strategy #1 (streaming adapter) | -3 to -5× on chr22 | unchanged |
| **6.3.2** | Strategy #2 (Liftoff in-memory) | additional -1.5× | unchanged |
| **6.3.3** | Strategy #3 (locus-major fusion) | additional -2× when combined with `--threads N` | -50 % peak RSS |
| **6.3.4** *(optional, deferred)* | Strategy #4 (native bindings) | additional -10 to -20 % | small |

Total expected improvement after 6.3.1-6.3.3: **8-12× faster** end-to-
end on a chr22 fixture vs the Phase 5 baseline; **~50 % lower peak
RSS**; same byte-identical output GFF3.

---

## 8. NCBI compliance carry-over

All four strategies preserve the GFF3 invariants the Phase 5
validator now enforces:

- Strategy #1: streamed bytes are ingested through gffbase's parser,
  which carries the same NCBI checks the validator does.
- Strategy #2: Liftoff's in-memory features are converted to
  `gffbase.Feature` objects whose `__str__` round-trips correctly
  (verified Phase 6.1).
- Strategy #3: per-locus output goes through the same
  `Lifton_GENE.write_entry` chain — no new serialisation paths.
- Strategy #4: `mappy` returns alignment objects, not GFF lines, so
  GFF3 emission is still our `lifton/io/gff_writer.py` (Phase 4
  Step 8) job.

Strict-mode (`--strict-gff`) gates remain unchanged across all four
strategies.

---

## 9. Decision matrix for the user

| Question | Answer per strategy |
|---|---|
| **Ship in next sprint?** | #1 yes (low risk); #2 yes (medium); #3 medium-large; #4 no |
| **Touches vendored Liftoff?** | #1 no; #2 yes (additive entrypoint); #3 maybe (depends on locus-yield refactor); #4 indirectly |
| **New runtime dependencies?** | #1 none (uses gffbase already vendored); #2 none; #3 none; #4 `mappy` + new `pyminiprot` |
| **Output GFF3 byte-identical to Phase 5 baseline?** | All four — guarded by integration test |
| **Phase 5 zero-bug baseline preserved?** | All four — guarded by full pytest suite |
| **Compatible with Phase 4 §9 parallelism plan?** | #1 yes; #2 yes; #3 *enables* it; #4 yes |
| **Engineering effort (rough)** | #1 ~1 day; #2 ~3 days; #3 ~1.5 weeks; #4 ~3-4 weeks (incl. binding) |

---

## 10. My recommendation

**Land Strategy #1 first; gate Strategy #2 behind its merge; defer
#3 to its own phase; treat #4 as a nice-to-have for later.**

Rationale:
- Strategy #1 buys 60-70 % of the available wall-clock improvement
  for ~10 % of the total engineering cost.
- Strategy #2 is the natural follow-up — same adapter surface, just
  reaches deeper into Liftoff.
- Strategy #3 is the architecturally pure answer but requires
  rewiring the global `ref_features_dict` / `tree_dict` /
  `m_id_2_ref_id_trans_dict` plumbing; budget a separate phase.
- Strategy #4 unlocks 10-20 % more *after* the structural wins are
  banked.

If the user wants the maximum impact in the shortest time:
**Phase 6.3.1 = Strategy #1 only.** That's the surface I'm prepared
to specify and implement next, pending your selection.

No code has been modified.
