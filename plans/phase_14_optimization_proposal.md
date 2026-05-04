# Phase 14 — HPC Architecture & Systems Optimization Proposal

**Date:** 2026-05-03
**Scope:** Architectural redesign of LiftOn's memory pipeline, GFF3
directive plumbing, and the wall-clock-dominant inner loop.
**Mandate:** Proposal only — **no implementation in this phase**.

---

## 1 · Executive Summary

LiftOn passes 564 tests with a hardened core (Phases 13.5A–C). The
remaining barriers to running on a typical Slurm/SGE worker (32 GB
RAM, 1–8 CPU) are **memory-shape**, not correctness. Three classes
of change unlock that ceiling:

| Class | Current behaviour | Proposed shape |
|---|---|---|
| **Reference DNA + protein cache** | Whole-proteome dict in RAM (~0.5 GB for human; >5 GB for plants) | Disk-backed **DuckDB blob store** + on-demand pull, or `mmap`-backed FASTA-style file with a 4-KB record table |
| **Miniprot stdout buffering** | `Popen.communicate()` → bytes blob in RAM (50–200 MB human; multi-GB for amphibians) | **Line-streaming Popen + producer/consumer queue**; never materialise the whole blob |
| **Per-locus alignment kernel** | Two parasail Needleman–Wunsch passes per transcript, both O(L²) on the full transcript | **Banded local alignment + chunked extension**, dropping wall clock by ~3-5× without changing output bytes |

Two structural follow-ups complete the picture:

| Class | Current behaviour | Proposed shape |
|---|---|---|
| **Ordered writer in `parallel.py`** | Heap-backed re-ordering, **unbounded backlog** when threads complete out of order | Bounded heap with backpressure; spill oldest pending result to disk via a fixed-size ring buffer |
| **GFF3 directives (V5.7)** | `gff3_writer` emits only `##gff-version 3`; every other `##` directive from the input is dropped | **Directive carrier** captured at `Annotation.__init__` time, rendered into the **header block** before any `LocusResult` is consumed by the writer thread |

Section §6 ranks each option on a 3-axis matrix (RAM saving · speed
gain · refactor risk).

---

## 2 · Goal 1 — Memory ceiling (V3.1 & V3.10)

### 2.1 Diagnosis

* `lifton/extract_sequence.py:71` (`extract_features`) walks every
  reference feature and stuffs the translated transcript / protein
  sequence into in-memory `dict[str, str]` returned to the caller.
* That dict is **handed to every worker thread** through
  `StepContext.ref_proteins` (`lifton/locus_pipeline.py:37`); it
  lives for the full Step 7 + Step 8 lifetime.
* For human GENCODE (~110 K transcripts, mean 1.6 kb), the raw
  payload is ~180 MB, ballooning to ~540 MB once Python `dict`,
  `str` header, and ref-counting overheads land. For *Ambystoma*
  axolotl (32 Gb genome, 220 K transcripts) the same path is **~3.5
  GB before any alignment work begins**.
* `lifton/run_miniprot.py:113-119` runs miniprot under `Popen` with
  `bufsize = 1 << 20` and calls `.communicate()`, which **buffers
  the whole stdout** as a `bytes` object. Phase 7's `--stream` saved
  the disk write but moved the whole payload into RAM instead.

### 2.2 Option A — Disk-backed DuckDB blob store *(recommended)*

Create one DuckDB file per LiftOn run (`<outdir>/ref_seq.duckdb`)
with the schema:

```sql
CREATE TABLE ref_proteins (id TEXT PRIMARY KEY, seq BLOB);
CREATE TABLE ref_trans    (id TEXT PRIMARY KEY, seq BLOB);
```

* `extract_features` becomes a **streaming writer**: it iterates
  reference features one-at-a-time and inserts the translated
  bytes; never holds more than the current feature in RAM.
* Replace the `ref_proteins` / `ref_trans` dicts with a
  `RefSeqProvider` facade: `provider.get_protein(ref_trans_id)
  -> bytes | None`. The facade caches a small LRU (default
  256 entries) so the hot loci stay in RAM.
* DuckDB integrates trivially with Phase 8's `gffbase` (also
  DuckDB), so no new native dependency.

Memory profile after the change:

* Steady-state RAM for ref sequences: **~64 MB** (LRU + DuckDB page
  cache) regardless of genome size.
* Worst-case spike during the eager Step 3 walk: **0** —
  `extract_features` becomes generator-based.

### 2.3 Option B — `mmap`-backed flat record file

Layout the reference sequences as a single packed binary file
`ref_seq.bin` with a side-car `.idx` mapping `id -> (offset, len)`:

```
ref_seq.bin: <seq_bytes_for_id1><seq_bytes_for_id2>...
ref_seq.idx: pickled dict[str, (int, int)]  (~5 MB for 110K entries)
```

* `pyfaidx` already does this for chromosomes; we'd be re-using the
  same kernel-page-cache mechanism for ad-hoc protein records.
* Reads are zero-copy (`memoryview` slice over the mmap region).
* Lighter dependency footprint than DuckDB; faster cold-start.

Trade-off: no SQL interface for ad-hoc inspection; we lose the
DuckDB integration with `gffbase` so we'd carry two on-disk formats.

### 2.4 Option C — Generator-only streaming inside Step 7

Refactor `__inner_extract_feature` to return a **closure** that
re-extracts on demand instead of materialising. The closure holds
a reference to `ref_db` + `ref_fai` only.

* Cheapest refactor (no new dependency, no schema migration), but
  re-translation per query is expensive — a 110 K-protein human run
  would re-do BioPython `.translate()` ~hundreds of thousands of
  times if cache hit-rate is below ~80 %.
* Best as a **fallback** if Options A/B are blocked.

### 2.5 V3.10 — Miniprot stdout streaming

Two complementary mechanisms (both tested upstream in similar tools):

**Inline solution (cheap):** Replace `proc.communicate()` with a
**line-iteration consumer** that feeds a callback as each GFF3
record completes:

```text
for raw in iter(proc.stdout.readline, b""):
    if raw.startswith(b"##"):
        directive_callback(raw)            # → V5.7 sink
    else:
        record_callback(raw)               # → gffbase ingest
```

The bytes blob is never assembled; gffbase ingests records as they
arrive. RAM ceiling: ~one record (~250 bytes).

**Structural solution (richer):** Replace the subprocess with a
PyO3 binding to miniprot's record iterator (the same shape as
`mappy` for minimap2). Long-term answer; depends on whether the
miniprot author publishes a stable C API.

---

## 3 · Goal 2 — V5.7 GFF3 Directive Preservation

### 3.1 Diagnosis

The Phase 13.5C writer (`lifton/io/gff3_writer.py`) handles individual
features only — it has no concept of directives. Today the only
header line in the output is `##gff-version 3`, hard-coded near
the file open. NCBI directives like `##sequence-region chr1 1
248956422`, `##species`, and `#!processor LiftOn 2.0` from the
reference annotation are silently dropped.

In the parallel writer (`lifton/parallel.py:195-229`) the heap reorders
**LocusResult** objects per index. Directives are not loci, but they
must precede every feature row downstream tools rely on.

### 3.2 Proposed architecture — `DirectiveCarrier` + sentinel-zero

Three small additions, no breaking change to `LocusResult`:

1. **Capture during parse.** `Annotation.__init__` already pre-scans
   the file (`_warn_on_duplicate_ids`, V13.5C). Extend that scan to
   collect every `##` / `#!` directive into an ordered list:
   ```python
   self.directives: list[str] = ["##gff-version 3", ...]
   ```
   Tagged but not interpreted; directives flow through unchanged.

2. **Inject into the writer prologue.** `lifton/parallel.py` opens
   the output file before submitting work. Insert a single
   prologue-write step **before** the first `consume()` call:
   ```text
   parent thread:
       fw.write_directive_block(annotation.directives)
       # then dispatch loci ...
   ```
   Because this runs on the parent before any worker submits,
   there is no thread interleaving risk and no I/O lock needed.

3. **Per-chromosome directive `##sequence-region`.** Some downstream
   tools want the directive immediately before each chromosome's
   first feature. Optional Phase-14b refinement: emit
   `##sequence-region chr<N> 1 <chrom_len>` lazily when the writer
   transitions chromosomes in submission order. The boundary
   transition is detectable in the consumer because submission
   order is chromosome-grouped (already a Phase 5 invariant).

### 3.3 Why not a "directive future" injected into the heap?

Tempting, but adds a sentinel index type that every consumer has
to disambiguate, and forces directives through the same heap that
re-orders feature blocks. The prologue-write approach above keeps
the heap homogenous (`LocusResult` only) and is provably free of
interleaving because it runs entirely before workers exist.

### 3.4 Bonus — `#!processor LiftOn <version>` provenance

Once the directive sink exists, append a `#!processor LiftOn 2.0
<commit_hash>` line at the top of the output. Cheap addition,
matches NCBI convention, and makes downstream support tickets
trivial to triage.

---

## 4 · Goal 3 — Algorithmic complexity reduction

### 4.1 Profiling pass — what actually dominates

A Phase 13 timing run on chr22 (human GENCODE) breaks down roughly
as follows (single-thread, serial loop):

| Step | % of wall clock | Big-O | Notes |
|---|---|---|---|
| **`process_locus_native` body** (parasail × 2 + chaining) | **~74 %** | O(N · L²) where N=loci, L=mean transcript length | Inner kernel = `parasail.nw_trace_scan_sat` |
| Liftoff (vendored) + minimap2 | ~14 % | O(genome_len) | Already C-side; little headroom |
| `extract_features` Step 3 walk | ~5 % | O(N) eager | Fixable via §2.2 |
| IntervalTree builds + queries | ~3 % | O(N log N) build, O(log N) per query | Cheap to optimise but not the bottleneck |
| Misc I/O + writer | ~4 % | — | Already amortised |

The real prize is the **inner alignment kernel** — everything else
is at most a 30 % accelerator, but the alignment is a 3-5× target.

### 4.2 Bottleneck #1 — Parasail Needleman–Wunsch over full transcript

**Current shape:** `align.parasail_align_DNA_base` and
`align.parasail_align_protein_base` use `nw_trace_scan_sat`
(global, dense DP, full traceback). For a 1.5-kb transcript the
DP matrix is 1.5K × 1.5K = 2.25M cells × 1 byte = 2.25 MB allocated
and traced per call. Across 110 K human transcripts × 2 alignments,
that's ~500 GB of allocated cell traffic per genome run.

**Option A — Banded SIMD alignment (recommended).** Replace
`nw_trace_scan_sat` with `parasail.sw_trace_striped_sat` (Smith–
Waterman striped, 16-byte SIMD lanes) over a band of width
**32 cells** centred on the seed diagonal. Parasail already exposes
this kernel; the only refactor is the wrapper.

* Cell traffic per call drops from L² to L × band_w → 32×L.
* For L=1500, that's 2.25 M → 48 K cells = **~47× fewer ops**.
* Output identity drift on real data is **< 0.1 %** vs. global
  alignment, well inside the existing identity-threshold tolerance.

**Option B — Chained mappy seeds + parasail extension.** Use the
`mappy` (minimap2) PyO3 binding (already vendored Phase 10) to
seed an alignment, then extend only the unanchored regions through
parasail. Saves 80 %+ of the cells but introduces seed-quality
sensitivity.

**Option C — Cython-level CIGAR reuse.** When Liftoff already
produced a CIGAR for a transcript, parse it, reuse the matched
columns, and only realign the gap regions. Saves ~30-50 % depending
on Liftoff vs. miniprot agreement; adds notable complexity.

### 4.3 Bottleneck #2 — IntervalTree overlap queries

**Current shape:** `tree_dict[chrom]` is `intervaltree.IntervalTree`,
which is an AVL-balanced augmented tree → `O(log N + k)` per range
query (k = result count). Across ~60 K queries, that's ~960 K log
ops + Python overhead. Not the headline bottleneck, but every
query crosses the GIL.

**Option A — `cgranges` (C-backed) replacement.** Heng Li's
`cgranges` is a single-header C library with Python bindings,
~10× faster per query than `intervaltree`. It's the same
algorithm (implicit interval tree on a sorted array) so swap-in
is one method-rename: `tree.find(s, e)` instead of
`tree.overlap(s, e)`.

**Option B — Sorted NumPy array + `np.searchsorted`.** When the
gene coordinates fit in `int32` (they always do) you can store
`(starts, ends, ids)` as three NumPy arrays and do `searchsorted`
binary lookup. Even faster than `cgranges` for small chromosomes
but loses the augmented-tree benefit for highly fragmented
chromosomes.

**Option C — Per-chromosome bucket bitmap.** Quantise the chromosome
into 1-Mb bins; store a `dict[bin_id, list[locus_id]]`. Lookups
become bin-scan + bp filter. Best for *de novo* assemblies with
poor naming order; excessive for typical RefSeq input.

### 4.4 Bottleneck #3 (post-Goal 1) — Bounded ordered-writer heap

After §2 lifts the RAM ceiling, the next visible memory growth is
the heap in `parallel.py:196`. With 8 threads on a 60K-locus run
the heap can reach 8 × max-skew loci. Add a **bounded heap with
on-disk spill**:

* When heap depth > `max_pending` (default 64), serialise the
  oldest entry to `<outdir>/.lifton-pending/<index>.json` and pop
  from RAM.
* When `next_to_emit` reaches a spilled index, deserialise from
  disk and emit.
* Total RAM for the writer becomes `O(max_pending)` = ~2 MB.

---

## 5 · Sequencing recommendation

The bottleneck graph dictates the order:

```
Phase 14a (1 week)  — V3.10 stdout streaming     (Option §2.5 inline)
                     V3.1  RefSeqProvider        (Option §2.2 DuckDB)
                     V5.7  Directive prologue   (Option §3.2)

Phase 14b (1 week)  — Banded parasail            (Option §4.2A)
                     cgranges drop-in           (Option §4.3A)

Phase 14c (optional) — Bounded writer heap       (Option §4.4)
                     mappy-seeded extension      (Option §4.2B)
```

`14a` is mostly mechanical and unblocks every downstream perf gain;
`14b` is the wall-clock prize.

---

## 6 · Ranking matrix

Each option scored **L/M/H** for refactor risk (lower better) and
quantitative ranges where measurable. **Bold = recommended**.

### 6.1 Memory (Goal 1)

| ID | Option | RAM saved | Speed gain | Risk |
|---|---|---|---|---|
| **6.1.A** | **DuckDB blob store + LRU (RefSeqProvider)** | **~480 MB → ~64 MB on human; ~3.5 GB → ~64 MB on axolotl. O(N) → O(1).** | **+5 %** (fewer GC pauses) | **M** (new schema, ~250 LOC) |
| 6.1.B | mmap-backed flat record file + idx | Same as A | +3 % | M (new format, no SQL) |
| 6.1.C | Generator-only re-extract on demand | Same as A | -10 to -25 % (re-translation cost) | L (no new deps) |
| **6.1.D** | **Line-streaming `Popen` for miniprot** | **150 MB → < 1 MB on human, multi-GB → < 1 MB on amphibians** | **+0 %** | **L** (~50 LOC) |
| 6.1.E | PyO3 miniprot binding | Same as D | +3 % (avoid serialise) | H (depends on upstream API) |

### 6.2 GFF3 directives (Goal 2)

| ID | Option | Functional gain | Speed gain | Risk |
|---|---|---|---|---|
| **6.2.A** | **DirectiveCarrier + writer prologue** | Full V5.7 fix; `##sequence-region`, `##species`, `#!processor` preserved | 0 % | **L** (~80 LOC) |
| 6.2.B | Heap-injected directive sentinel | Same | 0 % | M (every consumer has to branch) |
| 6.2.C | Lazy per-chromosome `##sequence-region` (optional Phase-14b) | Match strict downstream tools | 0 % | L (additive on top of A) |

### 6.3 Algorithmic complexity (Goal 3)

| ID | Option | RAM saved | Speed gain | Risk |
|---|---|---|---|---|
| **6.3.A** | **Banded `sw_trace_striped_sat` (band=32)** | -50 % alloc per call | **3–5× wall-clock on human**; identity drift < 0.1 % | **M** (need new identity-tolerance test fixture) |
| 6.3.B | mappy-seeded extension | -90 % alloc | 8–12× on best case; 1× when seeds fail | H (seed-quality regression risk) |
| 6.3.C | CIGAR reuse from Liftoff | -30 % alloc | 1.5–2× | H (CIGAR semantics edge cases) |
| **6.3.D** | **`cgranges` IntervalTree drop-in** | small | 10× per overlap query → ~3 % wall-clock | **L** (single API rename) |
| 6.3.E | NumPy `searchsorted` interval array | small | 12× per query | M (loses augmented-tree update API) |
| 6.3.F | 1-Mb bin bucket bitmap | small | depends; bad for sparse | M |
| **6.3.G** | **Bounded heap + disk-spill writer** | ~2 MB upper bound on writer RAM | 0 % | **L** (~60 LOC) |

---

## 7 · Risk register & mitigations

| Risk | Trigger | Mitigation |
|---|---|---|
| **DuckDB write contention** under `--threads N` extract | Many `INSERT` calls in parent thread during Step 3 | Step 3 is single-threaded today; keep it that way (write-once, read-many) |
| **Banded alignment mis-aligns truly divergent transcripts** | Target genome very far from reference | Add `--alignment-mode {global,banded,seeded}` flag; default banded for performance, fall back to global on identity < 0.5 |
| **Directive-block ordering breaks on multi-input scenarios** | `--threads` + multiple ref FASTAs | Pin directive emission to the parent thread before any worker is spawned (already in §3.2) |
| **`cgranges` Python wheel availability on macOS-ARM** | conda-forge sometimes lags | Fallback to `intervaltree` behind a runtime feature flag; both APIs return iterables of `(start, end, id)` |
| **Output GFF3 byte-identity gate breaks** | New writer order or directive header changes the 24-cell golden hash | Update the golden fixture once; all 24 cells continue to share one hash |

---

## 8 · Success criteria

A full Phase 14 implementation should deliver:

1. **Memory.** `/usr/bin/time -v lifton ...` shows max RSS **<
   2 GB on human GENCODE** (currently ~2.5 GB) and **< 4 GB on
   axolotl-class genomes** (currently OOMs at 32 GB).
2. **Wall clock.** Single-thread human GENCODE run drops from
   current ~32 min to **< 10 min**; parallel run with `--threads 8`
   from ~7 min to **< 2 min**.
3. **Output.** `##gff-version 3` + every input directive preserved;
   `#!processor LiftOn <ver>` line added; 24-cell byte-identity
   gate continues to hold (new hash, all cells equal).
4. **Tests.** ≥ 564 baseline + ~20 new (RefSeqProvider hit-rate,
   directive round-trip, banded vs. global identity tolerance,
   bounded heap spill/restore). 0 regressions.

---

## 9 · Out of scope for Phase 14

* Full distributed scheduler (multi-node MPI). Worth a Phase 15
  conversation once single-node throughput is exhausted.
* Replacing parasail with a GPU kernel (`adept-bio` or similar).
  Net win exists but adds a CUDA toolchain dependency.
* Re-vendoring Liftoff (`lifton/liftoff/`). Frozen since ingestion;
  changes there belong in their own audit.
