# LiftOn — Phase 6.4: Step-by-Step Roadmap to Full Native Integration

> **Goal:** reach Strategy #4 (in-process `mappy` + `pyminiprot` with
> locus-major fusion) without ever breaking the Phase 5 zero-bug
> baseline. Each step is independently shippable, each step has its
> own coverage / regression / golden-output gate, and each step is
> reversible by a single `git revert`.
>
> **Reference standards:**
> - NCBI GFF3 — every step's output is gated by `--strict-gff`.
> - Phase 5 baseline — **297 tests, 0 failures, 0 xfails.** This
>   number only ever goes up; commits that drop it get reverted.
>
> **Total schedule estimate:** ~6 weeks of focused work for the four
> sub-phases (6.4.A through 6.4.D), excluding any unforeseen `pyminiprot`
> binding work.

---

## North-star architecture (Phase 6.4.D end-state)

```
                ┌──────────────────────────────────────────────────────────┐
                │   for each gene g in stream_reference_genes(ref_gff):    │
                │       hits_dna  ← mappy.Aligner.map(g.dna)               │  in-process
                │       hits_prot ← pyminiprot.align(g.protein)            │  in-process
                │       lifted    ← lift_in_memory(g, hits_dna)            │
                │       chained   ← chain(lifted, hits_prot)               │
                │       rescued   ← orf_rescue_if_needed(chained)          │
                │       gff_writer.write(rescued)                          │  ordered emit
                │       del lifted, chained, rescued                       │  RAM bound
                └──────────────────────────────────────────────────────────┘
                                          │
                  ProcessPoolExecutor distributes loci across workers,
                  ordered-writer buffer preserves deterministic output.
```

No SAM files. No intermediate GFFs. No SQLite. No fork-per-aligner.
Output identical to the Phase 5 baseline.

---

## Sub-phase map

| Sub-phase | Strategy this lands | Net wall-clock gain (cumulative) | New deps | Risk | Effort |
|---|---|---|---|---|---|
| **6.4.A** | Strategy #1 — streaming adapter | 3-5× | none | Low | ~1 week |
| **6.4.B** | Strategy #2 — vendored-Liftoff in-memory | 4-7× | none | Medium | ~1 week |
| **6.4.C** | Strategy #3 — locus-major fusion + parallelism | 8-12× + 50 % RSS drop | none | High | ~2 weeks |
| **6.4.D** | Strategy #4 — `mappy` + `pyminiprot` native bindings | 9-14× + lower latency | `mappy`, build new `pyminiprot` | Very high | ~2 weeks |

Each sub-phase ends only when **all** of the following hold:
- Full pytest suite green (count never decreases).
- Coverage gate (`lifton_class.py` ≥ 90 %, `lifton_utils.py` ≥ 90 %)
  preserved.
- Byte-identical golden GFF3 vs Phase 5 baseline (or a documented,
  reviewed normalisation).
- Wall-clock benchmark on the chr22 fixture meets the sub-phase target.
- Phase report written to `plans/phase_6_4_<letter>_<topic>.md`.

---

# Sub-phase 6.4.A — Streaming Adapter Layer

> **Strategy #1 of the brainstorm.** Lowest-risk, biggest-bang-for-buck
> first move. Eliminates two of the three SQLite materialisations.

## Why this comes first
- No vendored-Liftoff surgery; no new dependencies.
- Reuses the Phase 6.2 adapter surface that's already approved.
- Buys ~60-70 % of the total available speedup at ~10 % of the
  total engineering cost.
- If anything else slips, this single landing already pays for itself.

## Design

### A.1 — Add `from_string=` plumbing in `lifton/gffbase_adapter.py`
Phase 6.2's adapter already wraps `gffbase.create_db`. Extend it:

```python
def build_database_from_string(
    *, gff_text: str | bytes, dbfn: str = ":memory:",
    infer_genes: bool, infer_transcripts: bool, force: bool,
) -> _gffbase.FeatureDB:
    """Ingest a GFF3 blob held in RAM, no intermediate file."""
    return _gffbase.create_db(
        gff_text, dbfn=dbfn,
        from_string=True, force=force,
        disable_infer_genes=not infer_genes,
        disable_infer_transcripts=not infer_transcripts,
    )
```

`gffbase.create_db` already supports `from_string=True`
(`lifton/gffbase/create_db.py:53-62`). Phase 6.4.A consumers pass the
captured stdout of `minimap2`/`miniprot` straight in.

### A.2 — Stream miniprot into RAM
**File:** `lifton/run_miniprot.py:62-67` (today writes a file).

Replace `subprocess.run(..., stdout=fw)` with:
```python
proc = subprocess.Popen(
    [miniprot_path, "--gff-only", tgt_genome, ref_proteins_file]
    + args.mp_options.split(),
    stdout=subprocess.PIPE, stderr=subprocess.PIPE,
    bufsize=1 << 20,
)
gff_bytes, err_bytes = proc.communicate()       # drains both pipes
if proc.returncode != 0:
    return None
return gff_bytes                                 # caller passes to gffbase
```

The caller in `lifton/lifton.py:375` switches from
`Annotation(<file_path>)` to `Annotation(gff_bytes_or_path,
backend="gffbase")`. The Phase 6.2 `Annotation` class learns to accept
either a path string or a `bytes` blob.

### A.3 — Stream Liftoff's already-in-memory output
The vendored Liftoff already builds `lifted_feature_list` in RAM
(`liftoff/liftoff_main.py:36, 42`); the only reason it touches disk is
that `write_new_gff.write_new_gff(...)` serialises to
`liftoff_outdir/liftoff.gff3` for the parent process to re-parse.

**File:** `lifton/run_liftoff.py:33-58`. Two-step refactor inside this
sub-phase:
1. **A.3.a:** keep the GFF3 write but skip the SQLite re-ingest by
   reading the GFF3 file's bytes back into `gffbase.create_db(...,
   from_string=True)`. Saves S4's ingest cost without touching
   Liftoff's API yet.
2. **A.3.b:** *(deferred to 6.4.B)* eliminate the GFF3 write entirely.

### A.4 — Make `Annotation` polymorphic on input type
**File:** `lifton/annotation.py:Annotation.__init__`.

```python
def __init__(self, source, ..., *, backend=None):
    if isinstance(source, (bytes, bytearray)):
        # In-memory GFF3 blob
        self._mode = "blob"
        ...
    elif os.path.exists(source):
        self._mode = "path"
        ...
    else:
        raise FileNotFoundError(source)
```

The gffutils branch refuses `bytes` input (raises a clean error); the
gffbase branch routes through `build_database_from_string`. Default
backend stays `gffutils` so the Phase 5 baseline is unchanged.

## Tests

### Unit tests (~15 new)
`tests/test_streaming_adapter.py`:
- `build_database_from_string` round-trips a synthetic GFF3 blob
  and the result of `db.features_of_type` matches the file-based
  ingest byte-for-byte.
- `Annotation(<bytes>, backend="gffbase")` works.
- `Annotation(<bytes>, backend="gffutils")` raises a clear error.
- `Popen + communicate` deadlock guard: oversized fake-miniprot
  stdout (~10 MB) ingests cleanly without hanging.
- Failure modes:
  - miniprot exits non-zero → return `None`, pipeline continues with
    Liftoff-only.
  - miniprot stdout is empty → return `None`.
  - miniprot stderr contains `"ERROR"` even with exit 0 → return
    `None` (mirrors today's logic).

### Integration tests (~5 new)
`tests/test_pipeline_streaming.py` parametrises
`tests/test_integration_pipeline.py` over `--stream={off,on}`:
- Asserts byte-identical output GFF3.
- Asserts `score.txt`, `unmapped_features.txt`, `extra_copy_features.txt`,
  `mapped_feature.txt` are byte-identical.
- Asserts no intermediate `.gff3_db`, `.duckdb`, or `miniprot.gff3`
  files are written when `--stream` is on (verified by recursive
  directory diff on `lifton_outdir/`).

### Property tests (~3 new)
`tests/test_streaming_property.py` (Hypothesis):
- For any random GFF3 blob the validator accepts, `from_string`
  ingest produces the same `db.features_of_type` count as
  file-based ingest.
- For any random `Popen` stdout containing valid GFF3, no record is
  truncated by `communicate()`.

### Performance gate (informational)
`tests/perf/test_streaming_speedup.py` (skipped by default,
`pytest -m perf`):
- chr22 reference + chr22 target.
- Asserts `streaming_total_walltime <= 0.5 * baseline_total_walltime`.
- Asserts no `.gff3_db` SQLite file appears in `lifton_outdir/`.

## Acceptance gates for 6.4.A
- [x] All 297 + ~23 new tests green.
- [x] Coverage on `lifton_class.py` and `lifton_utils.py` still ≥ 90 %.
- [x] Byte-identical golden GFF3 across `--stream={off,on}`.
- [x] chr22 wall-clock ≥ 2× faster with `--stream` on (target: 3-5×).
- [x] Phase report written to `plans/phase_6_4_A_streaming.md`.

---

# Sub-phase 6.4.B — Vendored-Liftoff In-Memory Refactor

> **Strategy #2.** Builds on the streaming adapter to remove the
> last GFF3 materialisation (Liftoff's own output).

## Why this comes after 6.4.A
- The streaming-from-bytes plumbing already exists in
  `Annotation`, so this sub-phase only needs to teach Liftoff to
  *yield* its in-memory output instead of writing it.
- Touching the vendored Liftoff is risky; doing it after 6.4.A means
  we have a regression net specifically for Liftoff output.

## Design

### B.1 — New entrypoint in vendored Liftoff
**File:** `lifton/liftoff/liftoff_main.py`.

```python
def run_all_liftoff_steps_inmemory(args, ref_db):
    """Same as run_all_liftoff_steps, but instead of calling
    write_new_gff.write_new_gff(...), returns the
    (lifted_feature_list, feature_db, ref_parent_order) tuple
    so the parent process can serialise / consume in RAM."""
    # ... same body as run_all_liftoff_steps minus the two
    # write_new_gff calls at lines 36, 42 ...
    return lifted_feature_list, feature_db, ref_parent_order
```

The original `run_all_liftoff_steps` becomes a thin wrapper that
calls `run_all_liftoff_steps_inmemory` and then writes the GFF3 — so
existing callers (CI smoke tests, manual users) keep working.

### B.2 — In-memory GFF3 emitter
**New file:** `lifton/liftoff/inmemory_emitter.py` (or extend
`write_new_gff.py`).

```python
def lifted_features_to_gff3_bytes(lifted_features, feature_db, args
                                  ) -> bytes:
    """Same logic as write_new_gff.write_new_gff but emits to a
    BytesIO buffer instead of a file."""
```

This guarantees the byte-identical output gate against the file path.

### B.3 — Wire the new path into `lifton/run_liftoff.py`
**File:** `lifton/run_liftoff.py`. New branch when `args.stream` and
`args.backend == "gffbase"`:

```python
lifted, db, order = liftoff_main.run_all_liftoff_steps_inmemory(
    liftoff_args, ref_db,
)
gff_bytes = inmemory_emitter.lifted_features_to_gff3_bytes(
    lifted, db, liftoff_args,
)
return gff_bytes              # pass to Annotation(<bytes>)
```

The fallback path (write GFF3 + re-parse) remains for the gffutils
backend or when `--stream` is off.

## Tests

### Unit tests (~10 new)
`tests/test_liftoff_inmemory.py`:
- `run_all_liftoff_steps_inmemory(args, ref_db)` returns a 3-tuple.
- The tuple's `lifted_feature_list` matches the dict that
  `run_all_liftoff_steps` populates (same keys, same coords).
- `lifted_features_to_gff3_bytes(...)` is byte-identical to
  `write_new_gff` writing to a tempfile and reading it back.

### Integration tests (~3 new)
`tests/test_pipeline_streaming.py` extends with a third parameter
value `--inmemory-liftoff={off,on}`:
- All 4 combinations (`stream`×`inmemory-liftoff` ∈ 2×2) produce
  byte-identical output GFF3.

### Differential tests (~5 new)
`tests/test_liftoff_differential.py`:
- Run vendored Liftoff against the same input via both code paths
  (file write + read vs in-memory yield).
- Assert: identical lifted-feature counts, identical coordinate
  ranges, identical `extra_copy_number` values, identical attribute
  ordering on every feature.

### Performance gate
- `tests/perf/test_liftoff_inmemory_speedup.py`: chr22 wall-clock
  shaved by ≥ 5 s vs sub-phase 6.4.A alone.
- Disk-write count (via `strace -c -e write` or the equivalent on
  macOS via `dtrace`) drops to zero for the Liftoff stage.

## Acceptance gates for 6.4.B
- [x] 297 + 6.4.A + 6.4.B tests green.
- [x] Coverage ≥ 90 % preserved.
- [x] Byte-identical output across all four 2×2 combinations.
- [x] chr22 wall-clock ≥ 4× faster than Phase 5 baseline (cumulative
      with 6.4.A).
- [x] Phase report at `plans/phase_6_4_B_liftoff_inmemory.md`.

---

# Sub-phase 6.4.C — Locus-major Fusion + Parallelism

> **Strategy #3.** The structural pivot. Per-locus pipeline; the
> reference DB stops being a global oracle and becomes a streaming
> source. Combines naturally with `--threads N` from the Phase 4
> Step 9 plan.

## Why this comes after 6.4.B
- 6.4.A and 6.4.B remove disk I/O. 6.4.C removes **global
  materialisation**, which is the next bottleneck once disk is
  cheap.
- Locus-major fusion is impossible while the pipeline still depends
  on the post-Liftoff write-then-re-ingest pattern that 6.4.B
  eliminates.
- Parallelism becomes "free" in the locus-major model — every
  locus is independent except for the IntervalTree cross-gene
  overlap check, which we already know how to share-as-snapshot.

## Design

### C.1 — Stream reference features
**New file:** `lifton/locus_pipeline.py`.

```python
def stream_reference_genes(ref_db) -> Iterator[Lifton_GeneSpec]:
    """Yield one Lifton_GeneSpec per gene in chromosome order.
    Each spec carries: ref_gene_id, attrs, coordinates, ordered
    transcript ids, ordered exon coords, the reference protein
    string (loaded lazily via RefSeqProvider), the parent's
    miniprot Target id (looked up in m_id_2_ref_id_trans_dict)."""
```

Uses the lazy `RefSeqProvider` from Phase 4 Step 4 (already planned)
so we never realise the full `ref_proteins` dict.

### C.2 — Per-locus pipeline coroutine
**File:** `lifton/locus_pipeline.py`.

```python
def process_locus(spec, *, target_fai, mp_index, args):
    """Run alignment + chaining + ORF rescue for ONE gene. Returns
    a fully-built Lifton_GENE ready to write."""
    # Today these calls are spread across run_liftoff.process_liftoff,
    # run_miniprot.process_miniprot, lifton_utils.LiftOn_miniprot_alignment,
    # protein_maximization.chaining_algorithm, Lifton_TRANS.orf_search_protein.
    # In 6.4.C they all happen in this one function on this one gene.
    ...
    return lifton_gene
```

### C.3 — Ordered, parallel writer
**File:** `lifton/parallel/locus_pool.py` (new).

Mirrors the design from `phase_4_roadmap.md` §9: a
`ProcessPoolExecutor` whose workers run `process_locus`, and a
deterministic ordered-writer buffer in the parent that drains
completed loci in submission order to keep output deterministic.

Default `--threads 1` keeps the serial path so 6.4.B → 6.4.C is a
zero-diff transition until the user opts in.

### C.4 — Cross-locus state
The only legitimately-shared state today is:
- `tree_dict` (per-chromosome IntervalTree of accepted Liftoff
  genes). Built incrementally; in the parallel path each worker
  gets a snapshot, and the parent re-merges.
- `m_id_2_ref_id_trans_dict` and `ref_id_2_m_id_trans_dict`. Built
  once before the locus loop starts; immutable from there. Workers
  receive a frozen copy via `initializer=`.
- `ref_features_dict`. Same — frozen once.

## Tests

### Unit tests (~25 new)
`tests/test_locus_pipeline.py`:
- `stream_reference_genes` yields genes in deterministic order.
- `process_locus` on a synthetic single-gene fixture produces an
  output equal to the legacy pipeline's output for that gene.
- Snapshot semantics: a worker that reads `tree_dict` after the
  parent has updated it gets the snapshot, not the live tree
  (tested via a fake-IPC harness).
- All five mutation rescue branches in `find_variants` reachable
  per-locus.

### Differential tests (~10 new)
`tests/test_locus_pipeline_differential.py`:
- For every gene in the chr22 fixture, assert the locus-major
  output equals the legacy whole-pipeline output for that gene
  (byte-identical lines after sorting both by `(seqid, start, type)`).

### Parallelism tests (~8 new)
`tests/test_locus_pipeline_parallel.py`:
- `--threads 1` produces output byte-identical to the serial path.
- `--threads 4` produces output byte-identical to `--threads 1`
  (determinism).
- `--threads N` with synthetic CPU-bound load shows ≥ 2× speedup
  on a 4-core box.
- No deadlock when a worker raises an exception (test via deliberate
  RuntimeError in a poison-pill locus).

### Memory gate
`tests/perf/test_locus_memory.py`:
- `mprof` peak RSS on chr22 with `--threads 1` is ≤ 50 % of the
  Phase 5 baseline.

## Acceptance gates for 6.4.C
- [x] 297 + 6.4.A + 6.4.B + 6.4.C tests green.
- [x] Coverage ≥ 90 % preserved on `lifton_class.py`, `lifton_utils.py`;
      new `locus_pipeline.py` lands at ≥ 95 % coverage.
- [x] Byte-identical output across `--threads ∈ {1, 2, 4, 8}`.
- [x] chr22 wall-clock ≥ 8× faster than Phase 5 baseline at
      `--threads 4`.
- [x] Peak RSS ≤ 50 % of Phase 5 baseline.
- [x] Phase report at `plans/phase_6_4_C_locus_fusion.md`.

---

# Sub-phase 6.4.D — Native bindings (`mappy` + `pyminiprot`)

> **Strategy #4.** The endgame. Drop the `subprocess` calls entirely;
> alignment becomes a function call inside the worker.

## Why this comes last
- Without 6.4.A/B/C the per-call alignment volume isn't high enough
  to make the binding tax pay off.
- `mappy` is on PyPI today; `pyminiprot` is **not**. Phase 6.4.D
  has a pre-step (D.0 below) that builds the binding.
- Switching aligners last means everything else is already proven
  on the slower-but-correct subprocess path. We only have to prove
  the binding produces the same hits as the subprocess; the rest of
  the pipeline is unchanged.

## Sub-phase 6.4.D structure

### D.0 — Build `pyminiprot` (separate repo / sidecar package)
- Fork the upstream `lh3/miniprot` C source.
- Wrap `mp_align` and the alignment iterator with PyO3 (matches the
  toolchain we already have for `gffbase._native`) or pybind11.
- Public API target:
  ```python
  import pyminiprot
  idx = pyminiprot.Index(target_fa, threads=4)
  for hit in idx.align(protein_seq, opts="-S 0.5 ..."):
      yield hit.to_gff3_record()        # same shape as the CLI's GFF3 output
  ```
- Ship a wheel for arm64-mac, x86_64-linux, x86_64-mac.
- **This is a parallel project**, not in-tree to LiftOn. Phase 6.4.D
  cannot start in LiftOn until this lands.

### D.1 — Replace `minimap2` subprocess in vendored Liftoff
**Files:** `lifton/liftoff/align_features.py:53-66, 109`.

```python
# before
subprocess.run([minimap2_path, '-o', output_file, ...])

# after
import mappy
aligner = mappy.Aligner(target_file, preset="splice", n_threads=...)
for query_name, query_seq in iter_features_fasta(features_file):
    for hit in aligner.map(query_seq):
        yield hit                         # consumed by parse_alignment in-process
```

Eliminates the per-chr SAM file write **and** the pysam re-parse at
`align_features.py:120-130`.

### D.2 — Replace `miniprot` subprocess
**File:** `lifton/run_miniprot.py`. Once `pyminiprot` exists:

```python
import pyminiprot
mp = pyminiprot.Index(target_fa, threads=args.threads)
hits = list(mp.align(reference_protein_seq, opts=args.mp_options))
# hits → in-memory feature objects, no GFF3 round-trip
```

This dovetails with sub-phase 6.4.C's per-locus model: the index is
built once at `process_pool_initializer` time and queried per-locus.

### D.3 — Locus-pipeline integration
**File:** `lifton/locus_pipeline.py`. The `process_locus` coroutine
gains:
```python
def process_locus(spec, ..., minimap2_aligner, miniprot_index):
    dna_hits  = list(minimap2_aligner.map(spec.transcript_seq))
    prot_hits = list(miniprot_index.align(spec.protein_seq))
    ...
```

The `ProcessPoolExecutor` initializer hands every worker its own
`mappy.Aligner` and `pyminiprot.Index` instance (both are
thread/process-safe in their respective bindings).

## Tests

### Differential correctness tests (~30 new — biggest test budget)
This is the highest-risk sub-phase; we need exhaustive parity proof.

`tests/test_mappy_subprocess_parity.py`:
- For every chromosome in the chr22 fixture, run alignment via
  both subprocess `minimap2` and `mappy` Python.
- Assert: identical hit counts, identical CIGAR strings, identical
  query/target coordinates, identical mapping quality.

`tests/test_pyminiprot_subprocess_parity.py`:
- Same exhaustive parity for `miniprot`.
- For 1000 random reference proteins: same hit set per protein.
- For 100 deliberately-tricky proteins (frameshift, premature stop,
  selenocysteine, ribosomal slippage) the hit set is identical.

### Integration tests (~10 new)
`tests/test_pipeline_native.py`:
- chr22 end-to-end with `--native` flag on (toggles binding usage).
- Byte-identical output GFF3 vs the Phase 6.4.C subprocess path.
- Falls back gracefully to subprocess when `mappy` or `pyminiprot`
  is not installed.

### Robustness tests (~5 new)
- `mappy.Aligner` build on a malformed FASTA → user-friendly error.
- `pyminiprot.Index` constructor on a missing target → clean fail.
- Worker-process crash isolation: a binding fault in one worker
  does not bring down sibling workers (test via segfault injection
  on a fake binding).

### Performance gate
`tests/perf/test_native_speedup.py`:
- chr22 wall-clock with `--native` ≤ 80 % of `--threads 4` from
  6.4.C (i.e., 20 % savings on top of locus fusion).
- Per-call latency: `aligner.map(seq)` ≤ 1 ms median for a 5 kb
  transcript on the chr22 reference; subprocess minimap2 takes
  ≥ 30 ms per fork+exec on the same input.

## Acceptance gates for 6.4.D
- [x] 297 + 6.4.A + 6.4.B + 6.4.C + 6.4.D tests green (~430+ total).
- [x] Coverage ≥ 90 % on core modules.
- [x] Byte-identical output GFF3 across all four sub-phases.
- [x] chr22 wall-clock ≥ 9× faster than Phase 5 baseline at
      `--threads 4 --native`.
- [x] `pyminiprot` published to PyPI as a sibling package.
- [x] Phase report at `plans/phase_6_4_D_native_bindings.md`.

---

## Cross-cutting test hierarchy

The whole roadmap leans on a layered test pyramid:

```
   ┌─────────────────────────────────────────────────┐
   │                Performance gates                │   informational
   │   (perf marker; only enforced on milestone PRs) │
   ├─────────────────────────────────────────────────┤
   │              Integration / golden                │   each sub-phase
   │   (byte-identical chr22 output GFF3 across all   │   adds variants
   │     flag combinations)                           │
   ├─────────────────────────────────────────────────┤
   │            Differential correctness              │   each new path
   │   (new path produces same answer as old path     │   gets a parity
   │     for every input in the corpus)               │   test suite
   ├─────────────────────────────────────────────────┤
   │               Property-based                     │   Hypothesis;
   │   (algebraic invariants on overlap math, ID      │   200 examples
   │     stripping, frame computation)                │   per strategy
   ├─────────────────────────────────────────────────┤
   │                 Unit tests                       │   per-function;
   │   (covers every branch reachable by the suite)   │   ≥ 90 % gate
   └─────────────────────────────────────────────────┘
```

The Phase 5 hermetic-pipeline contract carries forward: **no tests
spawn `minimap2` / `miniprot` / pre-built native binaries.** Every
new test layer monkey-patches the binary at the boundary where the
Phase 5 suite already does.

---

## Cross-cutting CI matrix

| Sub-phase | Default flag values | New CI dimensions |
|---|---|---|
| 6.4.A | `--stream=off` | `pytest -m "not perf"` only — perf gate runs nightly |
| 6.4.B | `--stream=off, --inmemory-liftoff=off` | + `--inmemory-liftoff=on` matrix |
| 6.4.C | `--threads=1, --locus-pipeline=off` | + `--threads={1,4}` × `--locus-pipeline={off,on}` |
| 6.4.D | `--native=off` | + `--native={off,on}` (only if `mappy`/`pyminiprot` installed; skipped otherwise) |

**Default flags stay off across all sub-phases** until a final cleanup
phase (Phase 6.5 — separate approval) flips them on and removes the
legacy code paths. This guarantees zero regression for users who
upgrade through 6.4.A → 6.4.D without changing their command line.

---

## Risk register and mitigations

| Risk | Likelihood | Impact | Mitigation |
|---|---|---|---|
| Streaming Popen deadlocks on full pipe buffer | Medium | High | Always use `communicate()` or thread-drained `stderr`; add a regression test with a 100 MB synthetic miniprot stdout |
| Vendored-Liftoff `lifted_feature_list` non-determinism | Low | High | Pin gene-order via the existing `ref_parent_order` list; differential test against the file-write path |
| Locus-major changes the cross-gene IntervalTree semantics | Medium | High | The IntervalTree is built incrementally either way; freeze a snapshot per worker, re-merge in the parent, regression-test against single-thread |
| `mappy` produces SAM with subtly different optional tags | Medium | Medium | Differential test asserts only the tags LiftOn actually reads; document the diff in `plans/phase_6_4_D_aligner_diff.md` |
| `pyminiprot` binding hits a memory leak on long runs | High | High | 24-hour soak test on the GENCODE v49 corpus; worker-recycle policy in the pool to bound RSS |
| User upgrades and sees a different output GFF3 | Low (we gate on byte-identity) | Critical | Each sub-phase commits a golden chr22 output; CI diffs against it; force-changes require a reviewed `plans/phase_6_4_*_diff_log.md` |

---

## Detailed unit-test budget

Cumulative test counts at each gate, rolled up from the per-sub-phase
tables above:

| Sub-phase | New unit | New integration | New differential | New property | New perf | Cumulative total |
|---|---:|---:|---:|---:|---:|---:|
| 6.4.A | 15 | 5 | — | 3 | 2 | **~322** |
| 6.4.B | 10 | 3 | 5 | — | 2 | **~342** |
| 6.4.C | 25 | — | 10 | — | 5 | **~382** |
| 6.4.D | 5 | 10 | 30 | — | 3 | **~430** |

Property tests are the smallest layer because gffbase already brings
its own Hypothesis suite; we only re-run our existing strategies
under the new code paths.

---

## Hand-off check before each sub-phase begins

Before starting any sub-phase X, confirm:

1. The previous sub-phase's `plans/phase_6_4_*.md` report is merged.
2. The Phase 5 baseline test count + this-roadmap's cumulative
   adds are all green on `devel`.
3. The byte-identical chr22 golden GFF3 is committed at
   `tests/golden/chr22_phase_<previous>.gff3` and CI diffs against
   it.
4. Coverage gates for `lifton_class.py` and `lifton_utils.py` are
   each ≥ 90 %.
5. The user has explicitly approved the sub-phase plan in chat.

If any of those is missing, the next sub-phase doesn't start.

---

## What this phase produced (this document only)

- The four-stage roadmap (6.4.A → 6.4.D) above.
- No code changes.
- No test changes.
- No source-tree changes.

**Awaiting your approval to begin sub-phase 6.4.A.**
