# LiftOn — Phase 1: Codebase Audit & Reconnaissance

## Context

This document captures the results of a read-only architectural audit of the
`LiftOn` repository (`/Users/chaokuan-hao/Documents/Projects/LiftOn`,
upstream `https://github.com/Kuanhao-Chao/LiftOn`) commissioned as the
foundation for a forthcoming refactoring effort. No functional code has been
modified. The audit answers four questions: (1) the end-to-end execution
flow, (2) the dominant in-memory data structures, (3) the dependency and
build surface, and (4) the structural red flags that should drive Phase 2
design choices.

LiftOn is a genome-annotation lift-over tool that combines DNA-level
alignment (a vendored fork of `liftoff` → `minimap2`) with protein-level
alignment (`miniprot`) and an in-house "protein-maximization" chaining
algorithm to produce a merged GFF3 annotation on a target genome.

---

## 1. Execution Flow (CLI → Output)

### Entry point
- `setup.py:16` registers `lifton = lifton.lifton:main`.
- `lifton/lifton.py:417-435` — `main()` prints a banner, calls
  `parse_args()` (lines 136-205), validates that `miniprot` is on `PATH`
  (`run_miniprot.py:6-24`), then dispatches to
  `run_all_lifton_steps()`.

### `run_all_lifton_steps()` — `lifton/lifton.py:208-415`
Nine strictly sequential phases:

| Step | Lines | Action |
|---|---|---|
| 0–1 | 213-238 | Load target & reference FASTA via `pyfaidx`; build reference `Annotation` (gffutils SQLite). |
| 2 | 241-245 | `lifton_utils.get_parent_features_to_lift()` + `get_ref_liffover_features()` partition features into protein-coding vs non-coding. |
| 3 | 248-267 | `extract_sequence` materialises reference transcript DNA & protein dicts; flags truncated proteins. |
| 4 | 300-304 | Run liftoff (in-process, vendored) and miniprot (subprocess). |
| 5 | 310-322 | Re-load liftoff & miniprot GFFs as gffutils DBs; open output writers (GFF, score.txt, unmapped, extra-copy, optional chain). |
| 6 | 327-328 | Build ref-ID ↔ miniprot-ID map (`lifton_utils.miniprot_id_mapping`, lines 384-405); seed per-chromosome `IntervalTree`s. |
| 7 | 338-346 | **Liftoff branch**: per gene → `Lifton_GENE` → align → protein-maximization chaining → ORF rescue → `write_entry()`. |
| 8 | 352-359 | **Miniprot-only branch**: for miniprot mRNAs not covered by liftoff, build a `Lifton_GENE` from miniprot coords, align, write. |
| 9 | 365-372 | `stats.print_report()`; close handles. |

### Vendored-liftoff invocation
`lifton/run_liftoff.py:8-30` does **not** shell out — it
`copy.deepcopy(args)`, redirects output to
`{outdir}/liftoff/liftoff.gff3`, then calls
`liftoff_main.run_all_liftoff_steps(args)` directly
(`lifton/liftoff/liftoff_main.py:10-42`). Liftoff itself shells out to
`minimap2` (`lifton/liftoff/align_features.py:59,64,109`).

### Miniprot invocation
`lifton/run_miniprot.py:27-54` — true subprocess (`subprocess.run` on
line 49) writing to a temp GFF.

### Output writing
- Primary GFF3: handle opened at `lifton.py:315`, written by
  `Lifton_GENE.write_entry` → `Lifton_TRANS.write_entry` →
  `Lifton_EXON.write_entry` → `Lifton_CDS.write_entry`
  (`lifton_class.py:158, 662, 718, 739`). Each leaf delegates to
  gffutils' `__str__()`.
- Side-channel files: `score.txt`, `unmapped_features.txt`,
  `extra_copy_features.txt`, optional `chain.txt`
  (`lifton.py:316-321`).

### Concurrency
Fully sequential; the only parallelism is the thread count forwarded to
`minimap2` and `miniprot`. No multiprocessing of genes.

---

## 2. Primary Data Structures

### Annotation object model (`lifton/lifton_class.py`, 742 LOC)
- `Lifton_GENE` (lines 57-186) — root; holds the gffutils gene entry,
  `is_protein_coding`, an ID→`Lifton_TRANS` dict, and copy-number.
- `Lifton_TRANS` (lines 218-680) — by far the largest class; owns an
  ordered exon list + dict, the 5-case CDS/exon overlap reconciler
  (`update_cds_list`, lines 266-451), and the in-house ORF finder
  (`__find_orfs`, 553-599).
- `Lifton_EXON` (682-725) and `Lifton_CDS` (727-743) — thin gffutils
  wrappers; an exon points to ≤ 1 CDS.
- `Lifton_Status` (12-22) and `Lifton_Alignment` (24-46) — alignment
  metadata; the latter stores aligned strings, CIGAR, identity, and the
  per-CDS protein-coordinate boundary dict.

### GFF storage — `gffutils` SQLite
- `lifton/annotation.py:20` tries to attach an existing `<gff>_db`
  cache; line 37 falls back to `gffutils.create_db(...)`.
- `get_protein_coding_features()` (83-88) and `get_feature_dict()`
  (150-155) iterate features and materialise an in-memory ID→record
  dict — i.e. the SQLite is also realised as a Python dict during a
  run.

### Reference sequence dicts — full upfront load
`lifton/extract_sequence.py:25-55` walks every gene in the reference
DB, extracts merged CDS DNA (`get_dna_sequence`, 73-85), translates
(`get_protein_sequence`, 92-94), and stores **all** transcripts in
`ref_trans` and **all** proteins in `ref_proteins` Python dicts. These
remain resident for the whole run.

### Coordinate / interval queries
- Per-chromosome `IntervalTree` built in `lifton/intervals.py:1-13`
  and seeded at `lifton_class.py:76-79` for O(log n) overlap queries
  during the miniprot-merge phase.
- `interlap` is declared in `setup.py` but **not actually imported**
  anywhere in `/lifton/` (dead dependency).

### Sequence access
- FASTA: `pyfaidx.Fasta` — lazy, index-backed; e.g.
  `lifton_class.py:461` (`exon.cds.entry.sequence(fai)`).
- Protein alignment: `parasail.nw_trace_scan_sat` with BLOSUM62
  (`align.py:50-67`).
- DNA alignment: `parasail` with a custom `ACGT*` 1/-3 matrix
  (`align.py:91-105`).
- Identity: gap-compressed BLAST identity in
  `get_id_fraction.get_AA_id_fraction` (23-40); region-restricted in
  `get_partial_id_fraction` (1-19), used by the chaining algorithm
  (`protein_maximization.py:37-50`).

### Memory hot-spots (rank-ordered)
1. **gffutils SQLite + materialised feature dicts** — full GFF tree.
2. **`ref_trans` / `ref_proteins`** — every reference transcript &
   protein resident for the whole run (~hundreds of MB on a mammalian
   genome).
3. **`Lifton_GENE/TRANS/EXON/CDS` tree** — mirror of the GFF
   hierarchy in Python objects.
4. **Per-chromosome `IntervalTree`s** — modest.
5. **Transient `Lifton_Alignment` matrices** — per-transcript-pair,
   discarded.

The "load everything into RAM" pattern in `extract_sequence.py:25-55`
combined with the gffutils dict materialisation is the single biggest
scalability ceiling.

---

## 3. Dependencies & Build System

### Declared Python dependencies (`setup.py:13`)
`numpy>=1.22.0`, `biopython>=1.76`, `cigar>=0.1.3`, `parasail>=1.2.4`,
`intervaltree>=3.1.0`, `interlap>=0.2.6` *(unused)*, `networkx>=3.3`,
`pyfaidx>=0.5.8`, `pysam>=0.19.1`, `gffutils>=0.10.1`, `ujson>=3.2.0`,
`pytest>=7.0.0`. **All lower-bound only — no upper caps.**

### Three-way Python version drift
| Surface | Python |
|---|---|
| `setup.py:14` | `>=3.6` (EOL) |
| `lifton.yml:18` | `3.11.9` |
| `Dockerfile:2` | `3.8-slim` (EOL) |
| `.github/workflows/tests.yml` | `3.12` |

Four sources of truth, all different.

### `setup.py` ↔ `lifton.yml` divergence
`lifton.yml` pins concrete versions (e.g. `pysam==0.22.1`,
`numpy==2.1.1`, `gffutils==0.13`) and adds `argcomplete`, `argh`,
`importlib-metadata`, `simplejson`, `sortedcontainers`, `zipp` that
are not declared in `setup.py`. Conversely `setup.py` declares
`interlap` and `pytest` runtime-deps; `pytest` belongs in extras, and
`interlap` is dead.

### External binaries (not vendored, not version-checked beyond presence)
- `minimap2` — invoked via `liftoff.align_features` at lines 59, 64,
  109.
- `miniprot` — `run_miniprot.py:20, 49`.
- The vendored `lifton/liftoff/` (~3 KLOC, 18 files) is itself a
  fork of upstream Liftoff with no recorded provenance commit.

### Build / packaging
- No `pyproject.toml`, no `setup.cfg` — legacy `setup.py` only;
  PEP 517/518 non-compliant, blocks modern resolvers (`uv`, `pdm`,
  `pip-compile`).
- No C/C++ extensions — pure Python (heavy lifting deferred to
  `parasail`, `pysam`, `pyfaidx`, external binaries).
- `Dockerfile` does `pip install --no-cache-dir .` against the
  outdated 3.8 base.
- CI (`.github/workflows`) runs a single shell example, not pytest.

### Tests
98 LOC total, both inside `lifton/liftoff/tests/`
(`test_basic.py` 48 LOC, `test_advanced.py` 49 LOC). **Zero tests** for
`lifton.py`, `lifton_class.py`, `lifton_utils.py`, `align.py`,
`protein_maximization.py`, `annotation.py`, `extract_sequence.py`,
`stats.py`, `variants.py`, or the `run_*.py` wrappers.

---

## 4. Structural Red Flags

### God modules
| File | Bytes | LOC | Notes |
|---|---|---|---|
| `lifton_class.py` | 34 KB | 742 | 5 classes, 57 methods; `Lifton_TRANS` alone owns ~460 lines including the 5-case CDS reconciler and a hand-rolled ORF finder. |
| `lifton.py` | 25 KB | 435 | Argparse + 9-step orchestrator interleaved. |
| `lifton_utils.py` | 20 KB | 529 | Kitchen-sink: arg massaging, ID mapping, alignment dispatch, overlap math. |
| `liftoff/find_best_mapping.py` | — | 392 | Vendored. |
| `liftoff/polish.py` | — | 331 | Vendored. |

### Coupling
Top-of-file imports show a star-shaped dependency on
`lifton_utils` and `lifton_class`:
- `lifton.py` imports `intervals, lifton_utils, annotation,
  extract_sequence, stats, logger, run_liftoff, run_miniprot,
  run_evaluation`.
- `lifton_utils.py` imports `align, lifton_class, run_liftoff,
  run_miniprot, logger`.
- `lifton_class.py` imports `align, lifton_utils, get_id_fraction,
  variants, logger`.

This forms a cycle at the package level (`lifton_utils ↔
lifton_class`), currently kept benign only by deferred attribute
access. Any import-order change risks breaking the package.

### Vendored fork debt
`lifton/liftoff/` is a full copy of upstream Liftoff (18 files,
~3 KLOC). No `UPSTREAM` marker, no diff log, no version pin. Future
`liftoff` bugfixes will not propagate.

### Memory & scaling
- `extract_sequence.py:25-55` — eager full-reference load.
- `annotation.get_feature_dict` — full materialisation of GFF.
- No streaming, no per-chromosome chunking, no parallelism over genes.

### Engineering hygiene
- `grep "-> "` across the package: **0 type hints**.
- Docstrings present but inconsistent; no module-level docstrings in
  the core files; class docstrings missing on the largest classes.
- `grep -E "TODO|FIXME|XXX|HACK"`: **0 markers** — the absence here
  is informational, not necessarily positive (it may indicate latent
  debt was never tagged).
- No hardcoded user paths in the production tree (positive).
- No `global` declarations / mutable module state (positive).

---

## Detailed Phase 1 Summary Report

### 1. Execution Flow
A single `main()` in `lifton/lifton.py` parses arguments, then runs a
linear nine-step pipeline (`run_all_lifton_steps`, lines 208-415):
load FASTAs → build reference gffutils DB → partition coding vs
non-coding → extract reference transcripts/proteins into RAM → run
**vendored** Liftoff in-process and **external** Miniprot via
`subprocess.run` → re-ingest both GFFs → walk every Liftoff gene
(align with parasail, run protein-maximization chaining against
the matching Miniprot transcript when present, ORF-rescue, write
entry) → walk Miniprot-only mRNAs and emit any that pass the
overlap and length-ratio filters → emit `score.txt`,
`unmapped_features.txt`, `extra_copy_features.txt`, optional
`chain.txt`, and the final GFF3. Output is serialised by a four-level
chain of `write_entry` calls in `lifton_class.py`. The whole
pipeline is sequential; the only parallelism is the thread count
forwarded to the external aligners.

### 2. Data Structures
The dominant Python object model is a four-level
`Lifton_GENE → Lifton_TRANS → Lifton_EXON → Lifton_CDS` tree defined
in `lifton/lifton_class.py` (742 LOC, 57 methods). Reference
annotation is held twice: as a `gffutils` SQLite database
(`annotation.py`) and as Python dicts produced by
`get_feature_dict`. Reference DNA and protein sequences are loaded
**eagerly** for every transcript into `ref_trans` and `ref_proteins`
dicts (`extract_sequence.py:25-55`) and stay resident for the whole
run. Coordinate math uses one `intervaltree.IntervalTree` per
chromosome (`intervals.py`); `interlap` is declared but unused.
Sequence I/O uses `pyfaidx` (lazy, indexed); pairwise alignment uses
`parasail` (BLOSUM62 for protein, custom 1/-3 ACGT matrix for DNA);
protein-maximization chains alignments via gap-compressed BLAST
identity in `get_id_fraction.py`. Memory pressure is dominated by the
gffutils DB + the eager reference-sequence dicts.

### 3. Dependencies
Pure-Python project, no compiled extensions. Runtime stack:
`numpy, biopython, cigar, parasail, intervaltree, networkx,
pyfaidx, pysam, gffutils, ujson` — all declared in `setup.py:13`
with **lower bounds only**. `interlap` is declared but unused;
`pytest` is misclassified as a runtime dep. The conda env
(`lifton.yml`) pins concrete versions that diverge sharply from
`setup.py`. Python version is specified four times with four
different answers: `setup.py>=3.6`, `lifton.yml=3.11.9`,
`Dockerfile=3.8-slim`, CI workflow=3.12. The build system is a
single legacy `setup.py` (no `pyproject.toml`, no `setup.cfg`).
External binaries — `minimap2` (used inside the vendored Liftoff)
and `miniprot` (subprocess) — must be on `PATH`; only `miniprot`'s
presence is checked. The `lifton/liftoff/` directory is a complete
in-tree fork of upstream Liftoff (~3 KLOC, 18 files) with no
provenance marker.

### 4. Structural Red Flags
**Critical:** four-way Python-version drift across `setup.py`,
`lifton.yml`, `Dockerfile`, and CI; `lifton_class.py` is a 742-LOC,
57-method god module whose `Lifton_TRANS` mixes coordinate
reconciliation, ORF discovery, and serialisation. **High:** all
dependencies are lower-bound-only (reproducibility / supply-chain
risk); a star-shaped coupling around `lifton_utils` and
`lifton_class` produces a benign-but-fragile package-level cycle;
the vendored Liftoff fork is unmaintained drift. **Medium:** the
`extract_sequence.py:25-55` upfront load of every reference
transcript/protein into Python dicts is the main memory ceiling;
`gffutils` features are also fully materialised; there is no
parallelism over genes. **Engineering hygiene:** zero type hints
across the package, inconsistent docstrings, ~98 LOC of tests
covering only the vendored Liftoff (the LiftOn core is untested),
no `pyproject.toml`, CI runs a shell example rather than a test
suite. Positives: no hardcoded user paths, no mutable module-level
globals, no `subprocess shell=True` injection vectors.

---

## Recommended Phase 2 Focus Areas (preview, non-prescriptive)
1. **Packaging & version unification** — single source of truth via
   `pyproject.toml`; align Python on 3.11+ across Dockerfile, conda,
   and CI; pin upper bounds.
2. **Decompose `lifton_class.py`** — split `Lifton_TRANS` into
   coordinate-reconciler, ORF-finder, and GFF-serialiser modules;
   introduce dataclasses + type hints.
3. **Decouple `lifton_utils`** — break the `lifton_utils ↔
   lifton_class` cycle by moving alignment-dispatch helpers behind a
   small interface module.
4. **Bound memory** — replace eager `ref_trans`/`ref_proteins` dicts
   with a per-chromosome streaming pass keyed off the gffutils DB;
   evaluate `pysam.FastaFile` vs `pyfaidx`.
5. **Vendored-liftoff strategy** — either re-base on upstream and
   record the commit, or carve the fork out into a clearly-owned
   sibling package.
6. **Test scaffolding** — add a chr22 smoke test runnable under
   pytest, with golden-output comparison; backfill unit tests for
   `lifton_class`, `align`, `protein_maximization`,
   `get_id_fraction`.

---

## Verification
This is a read-only audit; no code paths were exercised. Findings
were derived from direct reads of the listed files and `grep`/`wc`
across the package. The execution-flow trace can be cross-checked by
running:

```
lifton -g <ref.gff3> <ref.fa> <target.fa> -o /tmp/lifton_out.gff3
```

and comparing the order of writes to `lifton_output/score.txt`,
`stats/unmapped_features.txt`, and the final GFF3 against the step
table above.
