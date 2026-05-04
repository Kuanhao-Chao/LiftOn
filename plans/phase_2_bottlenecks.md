# LiftOn — Phase 2: Deep Bottleneck & Concurrency Analysis

> Read-only analysis. No code rewritten yet. All file:line references
> are anchors for Phase 3 implementation.

---

## 1. Memory & I/O Overhead

### 1a. `extract_sequence.extract_features` — eager dict materialisation
**Location:** `lifton/extract_sequence.py:25-35`

```
def extract_features(ref_db, features, ref_fai):
    ref_trans = {}
    ref_proteins = {}
    for feature in features:
        for locus in ref_db.db_connection.features_of_type(feature):
            __inner_extract_feature(...)   # mutates the two dicts
    return ref_trans, ref_proteins
```

**Why this hurts**

- Two dicts (`ref_trans`, `ref_proteins`) hold one entry per top-level
  feature × all transcripts. On a mammalian GFF (~60k genes,
  ~200k transcripts) this is several hundred MB of strings — and
  they live for the entire run because they are returned to
  `lifton.py:256` and then passed to `process_liftoff` /
  `process_miniprot` for every gene.
- The recursion in `__inner_extract_feature`
  (`extract_sequence.py:38-55`) walks the gffutils SQLite DB but
  re-reads every child via `ref_db.db_connection.children()` — each
  call is a SQL query.
- `get_dna_sequence` (line 73) reads exon slices via
  `pyfaidx`, but builds the per-transcript sequence as a Python
  `str` via repeated `+=` (O(n²) per transcript).
- The downstream code path that *does* use `ref_trans[id]` and
  `ref_proteins[id]` only reads each entry **once per matched
  transcript** (e.g. `lifton_utils.LiftOn_miniprot_alignment` and
  `align.protein_align` look up `ref_proteins[ref_id]` then never
  again). So the dicts are pure resident state for one-shot
  consumers.

**Refactor target — switch to a lazy, callable lookup**

Replace the two dicts with a single class that resolves on demand
and (optionally) caches with an LRU bound:

```
class RefSeqProvider:
    def __init__(self, ref_db, ref_fai, *, cache_size=4096): ...
    @lru_cache(maxsize=cache_size)
    def trans(self, ref_id: str) -> str | None: ...
    @lru_cache(maxsize=cache_size)
    def protein(self, ref_id: str) -> str | None: ...
```

Implementation notes:
- Backed by gffutils: a single `db_connection[ref_id]` lookup +
  one `children(..., featuretype=("exon"|"CDS"))` query per call.
- Drop the recursive double-read in `__inner_extract_feature`; call
  sites already know whether they want trans or protein.
- Use `io.StringIO` (or `''.join(parts)`) instead of `+=` for the
  exon concatenation — kills the O(n²) growth in
  `get_dna_sequence:73-85`.
- Replace `pyfaidx.Fasta(...).keys()` membership check
  (line 78) with a precomputed `frozenset(fasta.keys())` once at
  init — `keys()` on `pyfaidx` is cheap but called per transcript.
- Preserve the existing on-disk fallback at `lifton.py:259-263`
  (when the user supplies pre-extracted FASTAs): the provider
  should accept either a gffutils source or a `pyfaidx.Fasta` and
  expose the same `.trans(id)` / `.protein(id)` API. This lets the
  rest of the pipeline stay identical while we swap the producer.

**Order-of-magnitude impact**

For a human-scale annotation: peak resident size of
`ref_trans + ref_proteins` drops from ~400 MB (full materialisation)
to a bounded LRU at ~30 MB (cache_size=4096 × ~7 KB avg per pair),
plus the gffutils SQLite page cache (already present, not double-counted).

### 1b. `Annotation.get_feature_dict` — full re-materialisation
**Location:** `lifton/annotation.py:150-155`

```
def get_feature_dict(self, feature_types):
    id_to_feature = {}
    features = self.get_features_of_type(feature_types)
    for feature in features:
        id_to_feature[feature.id] = feature
    return id_to_feature
```

**Why this hurts**

- `get_features_of_type` (line 143) materialises the full list
  before the dict is built; both live in memory simultaneously.
- The returned dict is essentially a second copy of the SQLite
  index — gffutils already supports O(1) lookup via
  `self.db_connection[feature_id]`.
- Call sites in `lifton_utils` (search for `get_feature_dict`)
  use `dict[id]` lookups that are 1:1 replaceable by
  `db_connection[id]` (which hits the SQLite index directly,
  prepared-statement-cached).

**Refactor target**

- Replace `get_feature_dict` callers with a thin
  `FeatureLookup` wrapper that delegates to `db_connection[id]`
  on miss and (optionally) memoises hot keys.
- For places that *need* an ID set (e.g. for "is this in the
  reference?"), expose a separate `feature_id_set(types) -> set[str]`
  helper that builds only a `set` of IDs (10× lighter than a dict
  of full features).
- Same idea for `get_features_of_type` (line 143): change the
  return type to a generator / iterator and let the caller decide
  whether to listify.

### 1c. Other I/O hotspots worth fixing in the same pass

| File | Lines | Issue | Fix |
|---|---|---|---|
| `extract_sequence.py:73-85` | sequence concatenation `+=` | O(n²) | use list-and-join |
| `extract_sequence.py:78` | `if chrom not in fasta.keys()` per transcript | repeated dict materialisation by `pyfaidx` | precompute `chrom_set = set(fasta.keys())` once |
| `annotation.py:80-88` | `get_protein_coding_features` materialises a `[child for child in db.children(feature) if ...]` per gene | full list build just to test `len > 0` | use `any(c.featuretype=='CDS' for c in db.children(feature))` |
| `annotation.py:182-188` | `get_num_levels` walks every level via list comprehension | recursive list build | inline counter; break early |
| `lifton.py:262-263` | `Fasta(ref_trans_file)` / `Fasta(ref_proteins_file)` reload from disk every run | rebuild `.fai` index every invocation | add `lifton_utils.cached_fasta()` helper that respects mtime |

---

## 2. The "God Module" — `Lifton_TRANS` (`lifton_class.py:218-680`)

### Mixed responsibilities (current state)

`Lifton_TRANS` currently owns four orthogonal concerns:

| Concern | Lines | Methods |
|---|---|---|
| **Construction & exon/CDS bookkeeping** | 218-264 | `__init__`, `add_exon`, `add_cds`, `update_gffutil_entry_trans`, `mv_exon_idx` |
| **Coordinate reconciliation** (5-case CDS↔exon merge) | 266-451 | `update_cds_list` (the 185-line god method), `update_boundaries` |
| **Sequence assembly + alignment dispatch** | 453-524 | `get_coding_seq`, `get_coding_trans_seq`, `translate_coding_seq`, `align_coding_seq`, `align_trans_seq` |
| **ORF rescue + CDS boundary patching** | 526-660 | `orf_search_protein`, `__find_orfs`, `__update_cds_boundary`, `__iterate_exons_update_cds`, `__get_cds_frame` |
| **GFF serialisation** | 662-679 | `write_entry`, `print_transcript` |

### Decoupling target — three sibling modules, one slim model

Proposed split (no behaviour change in Phase 3):

```
lifton/
├── model/
│   └── transcript.py        # data only: Lifton_TRANS = dataclass-ish
├── reconcile/
│   └── cds_exon.py          # pure functions on (exons, cds_list)
├── orf/
│   └── rescue.py            # ORF finding + CDS boundary patching
├── seq/
│   └── assembly.py          # get_coding_seq, get_coding_trans_seq
├── align/
│   └── transcript_align.py  # align_coding_seq, align_trans_seq
└── io/
    └── gff_writer.py        # write_entry for GENE/TRANS/EXON/CDS
```

### Concrete extraction map

#### a) `reconcile/cds_exon.py`
Move `update_cds_list` (266-451) and the helper `update_boundaries`
(672-674) into pure functions:

```
def reconcile_cds_with_exons(
    exons: list[Lifton_EXON],
    cds_list: list[Lifton_CDS],
    strand: str,
) -> list[Lifton_EXON]: ...
```

Why pure: the method only reads `self.entry.strand` and mutates
`self.exons`. Make strand an argument, return the new list, let the
caller assign. This kills the implicit `self` coupling and makes the
five branches independently unit-testable. The five cases become five
named functions (`_reconcile_single_cds`, `_reconcile_single_exon`,
`_reconcile_init_head`, `_reconcile_inner`, `_reconcile_tail`).

Side effect: `lifton_utils.segments_overlap_length` (used 4× inside
`update_cds_list`) becomes the only cross-module dependency, and it
is already pure — so this module imports nothing from
`lifton_class`.

#### b) `orf/rescue.py`
Move `__find_orfs` (553-599), `__update_cds_boundary` (600-604),
`__iterate_exons_update_cds` (606-657), `__get_cds_frame` (659-660).
Recast as:

```
@dataclass
class ORFCandidate:
    start: int
    end: int
    identity: float

def find_best_orf(
    trans_seq: str,
    ref_protein_seq: str,
) -> ORFCandidate | None: ...

def patch_cds_boundary(
    exons: list[Lifton_EXON],
    orf: ORFCandidate,
    strand: str,
) -> None: ...
```

The `orf_search_protein` orchestrator (526-551) stays as the
top-level method on `Lifton_TRANS` but becomes a thin dispatcher
that calls into `seq.assembly`, `align.transcript_align`,
`variants.find_variants`, then conditionally `orf.rescue.find_best_orf`.

This also breaks the circular-feeling import:
`lifton_class.__find_orfs` currently does
`lifton_class.Lifton_ORF(...)` (line 578) — once `Lifton_ORF` moves
into `model/transcript.py`, the rescue module imports it cleanly.

#### c) `seq/assembly.py`
Move `get_coding_seq` (453-468) and `get_coding_trans_seq`
(470-494). These already take `(self, fai)` and only touch
`self.exons` — convert to:

```
def coding_seq(exons, fai) -> tuple[str, list, list[int]]: ...
def coding_and_trans_seq(exons, fai, strand: str) -> tuple[str, str]: ...
```

Bonus fix: `get_coding_trans_seq:480-482` builds `trans_seq` via
`Seq(...).upper() + Seq(...).upper()` then re-converts to `str` —
two unnecessary BioPython object constructions per exon. Pull the
upper-casing out of the loop.

#### d) `align/transcript_align.py`
Move `translate_coding_seq` (496-500), `align_coding_seq` (502-514),
`align_trans_seq` (516-524). All three are stateless wrappers around
`align.protein_align` / `align.trans_align` plus
`Lifton_Status` mutation — convert to free functions that
**return** the mutated status delta instead of mutating in place
(makes them safe to call from worker processes — see §3).

#### e) `io/gff_writer.py`
Move `write_entry` (`Lifton_GENE`, `Lifton_TRANS`, `Lifton_EXON`,
`Lifton_CDS` — `lifton_class.py:158, 662, 718, 739`). They are
currently tightly entangled with the four classes via `self.entry`
attribute access; lift them into a single visitor function so the
output format lives in one file and can grow alternative writers
(GTF, BED12) later.

### Net effect

`lifton_class.py` shrinks from 742 LOC to ≈ 230 LOC of pure data
classes; the five new modules average ~120 LOC each, all
unit-testable in isolation. The package-level cycle `lifton_utils ↔
lifton_class` (Phase 1, §3) collapses because the new modules
import only from `model/` and `lifton_utils`, never the reverse.

---

## 3. Concurrency Potential

### Current bottleneck profile

The Step 7 loop (`lifton.py:338-346`) is the hot path:

```
for feature in features:
    for locus in l_feature_db.features_of_type(feature):
        lifton_gene = run_liftoff.process_liftoff(...)  # parasail+ORF
        ...
        lifton_gene.write_entry(fw, transcripts_stats_dict)
```

`process_liftoff` per gene runs:
1. `Lifton_GENE` build (cheap)
2. **parasail protein alignment** (`align.py:50-67`) — CPU-bound,
   no shared state once `ref_proteins[ref_id]` is fetched
3. **parasail DNA alignment** (`align.py:91-105`) — same
4. Optional **protein-maximization chaining**
   (`protein_maximization.py:68-128`) — pure CPU
5. **ORF rescue** if mutations found (`__find_orfs:553-599`) — pure
   CPU, currently O(L) per frame × 3 frames per transcript

Steps 2-5 are embarrassingly parallel **per gene**. The only
serial dependency is the *output* — GFF entries should appear in a
deterministic order to keep diffs sane.

### Proposed parallelisation

**Use `concurrent.futures.ProcessPoolExecutor` with a producer/
ordered-consumer pattern**, not `multiprocessing.Pool.imap` — the
ordered-consumer pattern preserves output ordering with a small
reordering buffer instead of forcing chunk-level synchronisation.

```
from concurrent.futures import ProcessPoolExecutor, as_completed

def parallel_step7(features, l_feature_db, ..., n_workers):
    with ProcessPoolExecutor(max_workers=n_workers,
                             initializer=_init_worker,
                             initargs=(args, ref_genome_path,
                                       ref_db_path, ref_proteins_path,
                                       ref_trans_path)) as ex:
        futures = {}
        for feature in features:
            for locus in l_feature_db.features_of_type(feature):
                fut = ex.submit(_process_liftoff_remote,
                                locus.id, feature)
                futures[fut] = (locus.seqid, locus.start, locus.id)

        # Buffered, deterministic-order writer
        pending = {}
        next_idx = 0
        for fut in as_completed(futures):
            entry = fut.result()
            ...
```

### Why this works

- **Inputs are picklable bytes**: pass `(locus_id, feature_type)` to
  the worker, not the gffutils Feature object. The worker rehydrates
  via `db_connection[locus_id]` (cheap with the SQLite cache).
- **Reference data is read-only**: the worker's `_init_worker`
  opens its own `pyfaidx.Fasta(ref_genome)`, its own
  `pyfaidx.Fasta(target_genome)`, its own `gffutils.FeatureDB`, and
  (after §1 refactor) its own `RefSeqProvider`. All four are
  read-only file handles — perfect for forked-or-spawned workers
  on macOS (use `mp.set_start_method('spawn')` for safety).
- **No GIL contention**: `parasail` releases the GIL during the
  alignment kernel — but worker isolation buys us linear scaling
  on modern multicore boxes regardless. We pick processes not
  threads to also parallelise the pure-Python ORF search.
- **Determinism**: futures are kept in submission order; results
  are written in that order via a small heap-buffered consumer
  (`heapq.heappush(buf, (idx, entry))`, then drain prefix). Worst-
  case buffer depth ≈ `n_workers × 2`.

### What CANNOT be parallelised

- **Step 4** (`exec_liftoff` + `exec_miniprot`, `lifton.py:302-304`).
  Both already use their own thread/process pools internally
  (`-t` flag forwarded). Don't double-parallelise.
- **Step 5** (database creation, `lifton.py:310-314`). gffutils SQLite
  build is single-writer; can't be sped up here.
- **Step 8** (`lifton.py:352-359`, miniprot-only loop).
  Conceptually parallelisable identically to Step 7, but it
  reads from the *output* of Step 7 (the IntervalTree of accepted
  Liftoff genes via `tree_dict`). Either run after Step 7's writer
  drains, or share the IntervalTree as a `mp.Manager` object
  (slow). Recommendation: keep Step 8 sequential in v1; revisit
  once Step 7 is verified correct.
- **Final writer** (`fw.write` in `write_entry`). Single file
  handle, single thread.

### CLI integration

Add `--threads N` (alias `-t N`) to `lifton.parse_args` (currently
no such flag — verify in `lifton.py:136-205`). Default to
`min(8, os.cpu_count())`. `-t 1` must take the existing serial path
to keep regression diffs trivial during the rollout.

### Expected speedup

Step 7 dominates wall time on real GFFs (Phase 1 timing harness
exists at `lifton.py:375-414`). On a 10-core box we should see
~6-8× on Step 7, translating to ~3-5× end-to-end depending on the
liftoff/miniprot share of the run.

---

## 4. Environment & Dependency Debt

### 4a. Prune dead deps
- **`interlap`** — declared in `setup.py:13`, never imported. Drop.
- **`pytest`** — declared as runtime; move to
  `[project.optional-dependencies] test = [...]`.
- **`ujson`** — verify with `grep -rn "import ujson\|from ujson" lifton/`
  before removal; if unused, drop.
- **`networkx`** — used inside the vendored `liftoff/`
  (`find_best_mapping.py`). Keep.

### 4b. Unify Python to 3.10+
Single source of truth is `pyproject.toml` after migration. The
four current sources must converge:

| File | Current | After |
|---|---|---|
| `setup.py:14` | `>=3.6` | (file removed) |
| `lifton.yml:18` | `3.11.9` | `3.11` |
| `Dockerfile:2` | `python:3.8-slim` | `python:3.11-slim` |
| `.github/workflows/tests.yml` | `3.12` | matrix `[3.10, 3.11, 3.12]` |
| `pyproject.toml` (new) | — | `requires-python = ">=3.10"` |

Why 3.10+ (not 3.11+): `match`/structural pattern matching is
attractive for the 5-case `update_cds_list` rewrite; `int | None`
PEP-604 unions land cleanly; large user base still on 3.10.

### 4c. Migrate to `pyproject.toml` (PEP 621 + setuptools backend)
Minimum viable structure:

```
[build-system]
requires = ["setuptools>=68", "wheel"]
build-backend = "setuptools.build_meta"

[project]
name = "lifton"
version = "1.0.5"
requires-python = ">=3.10"
dependencies = [
    "biopython>=1.81,<2",
    "cigar>=0.1.3",
    "gffutils>=0.13,<0.14",
    "intervaltree>=3.1.0,<4",
    "networkx>=3.3,<4",
    "numpy>=1.22,<3",
    "parasail>=1.3.4,<2",
    "pyfaidx>=0.8,<1",
    "pysam>=0.22,<1",
]

[project.optional-dependencies]
test = ["pytest>=7", "pytest-xdist"]

[project.scripts]
lifton = "lifton.lifton:main"

[tool.setuptools.packages.find]
include = ["lifton*"]
```

Add upper bounds (every dep) — Phase 1 §3 flagged the lower-bound-
only policy as High risk. Generate a `requirements-lock.txt` via
`uv pip compile` for reproducible builds; the conda `lifton.yml`
becomes a thin wrapper that `pip install -e .` after creating the
3.11 env.

### 4d. Add the missing `minimap2` startup check
**Where:** `lifton/lifton.py:433` — currently only calls
`run_miniprot.check_miniprot_installation()` (defined at
`run_miniprot.py:6-24`).

`minimap2` is invoked transitively inside the vendored liftoff
(`lifton/liftoff/align_features.py:59,64,109`). It must be checked
**before** any heavy work begins.

```
# new helper, mirroring check_miniprot_installation
def check_minimap2_installation():
    if shutil.which("minimap2") is None:
        sys.exit("ERROR: minimap2 not found on PATH. "
                 "Install via `conda install -c bioconda minimap2`.")
```

Wire into `main()` immediately after `check_miniprot_installation()`.
While in there, also check `samtools` if any vendored liftoff path
uses it (grep `lifton/liftoff/` for `samtools`).

---

## Prioritised Refactor Targets

Ranked by `(impact × confidence) ÷ effort`. Top items deliver
biggest user-visible win for least implementation risk.

| # | Target | Why now | Files to touch | Effort | Risk |
|---|---|---|---|---|---|
| **1** | **Add `--threads` + parallel Step 7 via `ProcessPoolExecutor`** | Biggest wall-time win; isolated change; ordered-writer keeps output identical | `lifton/lifton.py:338-346`, new `lifton/parallel.py` | M | M (must validate determinism on a chr22 golden run) |
| **2** | **Replace `extract_sequence.extract_features` dicts with `RefSeqProvider` (LRU)** | Cuts peak RSS by ~hundreds of MB; required precondition for §1 worker isolation | `lifton/extract_sequence.py:25-94`, `lifton/lifton.py:256-263`, callers in `lifton_utils.py` | M | L |
| **3** | **Migrate to `pyproject.toml`, unify Python to 3.10+, prune `interlap`/move `pytest`, add upper bounds** | Unlocks modern toolchain; fixes four-way version drift | `pyproject.toml` (new), delete `setup.py`, `Dockerfile:2`, `lifton.yml`, `.github/workflows/tests.yml` | S | L |
| **4** | **Add `check_minimap2_installation()` startup guard** | One-line UX fix; prevents wasted compute on misconfigured systems | `lifton/lifton.py:417-435`, `lifton/run_liftoff.py` | XS | XS |
| **5** | **Extract `reconcile/cds_exon.py` from `Lifton_TRANS.update_cds_list`** | Largest single hot method (185 LOC, 5 cases); pure function rewrite is straight lift | `lifton/lifton_class.py:266-451` → new `lifton/reconcile/cds_exon.py` | M | M (high test debt — write golden tests first) |
| **6** | **Extract `orf/rescue.py`** | Removes 130 LOC of CPU-bound logic from the god module; cleanly callable by parallel workers | `lifton/lifton_class.py:553-660` → new `lifton/orf/rescue.py` | S | L |
| **7** | **Replace `Annotation.get_feature_dict` callers with direct `db_connection[id]` lookups** | Removes a redundant in-memory copy of the GFF index | `lifton/annotation.py:150-155` + grep callers | S | L |
| **8** | **Replace per-transcript `+=` string growth with list-and-join in `get_dna_sequence`** | O(n²) → O(n); cheap | `lifton/extract_sequence.py:73-85` | XS | XS |
| **9** | **Extract `seq/assembly.py` and `align/transcript_align.py`** | Makes alignment helpers picklable for workers in §1 | `lifton/lifton_class.py:453-524` | S | L |
| **10** | **Extract `io/gff_writer.py` (visitor for GENE/TRANS/EXON/CDS write_entry)** | Centralises output format; opens door for GTF/BED12 emitters later | `lifton/lifton_class.py:158, 662, 718, 739` | S | L |
| **11** | **Add a chr22 pytest smoke test (golden GFF comparison)** | Mandatory safety net before any of #1, #5, #6 land | `test/` (new) | S | XS |

**Recommended implementation order:**
`#11 → #4 → #3 → #8 → #2 → #7 → #6 → #5 → #9 → #10 → #1`

i.e. land tests + cheap wins + memory fix first; do the big
god-module split while the suite is green; finish with
parallelism on top of the now-decoupled, picklable units.

---

## Verification strategy (for Phase 3)

1. **Golden output**: capture today's `lifton_output/lifton.gff3`,
   `score.txt`, `unmapped_features.txt` for a chr22 reference vs
   chr22 target run. Each refactor PR must produce byte-identical
   outputs (or a documented, reviewed diff).
2. **Memory benchmark**: `mprof run lifton ...` before/after #2 on
   a full mammalian GFF.
3. **Wall-time benchmark**: existing `--measure_time` harness
   (`lifton.py:375-414`) — record baseline now; report Step 7 delta
   per refactor.
4. **Test suite**: extend `test/` with unit tests against the new
   `reconcile/`, `orf/`, `seq/`, `align/transcript_align.py`
   modules using small synthetic exon/CDS fixtures.
