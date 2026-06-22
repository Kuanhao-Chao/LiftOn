# LiftOn — Reimplementation-Grade Technical Specification

> **Subject:** LiftOn v1.0.8 (homology-based genome-annotation lift-over).
> **Source of truth:** the working tree at branch `devel`, commit `7152cca` (+ uncommitted
> Phase-16/17 changes), package root `lifton/`.
> **Purpose:** a complete, granular specification from which a developer can rebuild LiftOn
> — including its two vendored subtrees — *from scratch, without reading the original source.*

## How to read this document

Every claim is grounded in the real code and cited as `path:line` (relative to the repo root).
Algorithms are given as numbered steps with exact arithmetic, constants, thresholds, and the
verbatim boolean conditions of each branch. Data structures are enumerated field-by-field. The
spec is **language-agnostic in intent** but uses the actual Python signatures and names as the
authoritative reference; a faithful re-implementation in any language must reproduce the same
observable behaviour and, critically, the same output bytes (see *Byte-identity contract* below).

**Scope.** Full depth across all three bodies of code:

| Body | Location | LOC | Treatment |
|---|---|---|---|
| First-party LiftOn | `lifton/*.py`, `lifton/io/`, `lifton/native_bindings/` | ~8.6 K | Full |
| Vendored Liftoff fork | `lifton/liftoff/` (drives `minimap2`/`mappy`) | ~2.9 K | Full (§6.5) |
| Vendored gffbase | `lifton/gffbase/` (DuckDB FeatureDB + Rust `_native` + pure-Python fallback) | ~4.1 K | Full (§5.6) |

External binaries (`minimap2`, `miniprot`) and Python wheels (`parasail`, `pysam`, `pyfaidx`,
`gffutils`, `duckdb`, `pyarrow`, `mappy`, `biopython`) are specified at their **invocation
contract** — exact argv, flags, parameters, and the shape of data exchanged — since a
re-implementation reuses these real tools.

## Document map

- **§1 Architecture** — what LiftOn is, the module/dependency graph, the import cycle, fast-path flags.
- **§2 Execution lifecycle** — `main` → 11 pipeline steps → the per-locus engine → termination.
- **§3 Core algorithms** — sequence assembly, parasail alignment, identity math, variant
  classification, the protein-maximization chaining algorithm, ORF rescue, CDS↔exon reconciliation.
- **§4 Data structures** — the `Lifton_*` object model and all pipeline state objects.
- **§5 Annotation backend** — the `Annotation` class, feature partitioning, and the full gffbase internals.
- **§6 External tools** — Liftoff/miniprot/mappy integration and the full vendored Liftoff internals.
- **§7 I/O formats** — every input, intermediate, and output file format.
- **§8 Validation** — the input-side NCBI validator and the output-side GFF3 validator.
- **§9 Parallelization & performance** — the deterministic threading model and Phase-17 materialisation.
- **§10 Reimplementation guidance** — build order, determinism invariants, test strategy, appendices.

## Conventions

- **Coordinates** are **1-based, fully closed** (GFF3 convention). A feature spanning bases
  `start..end` has length `end - start + 1`. Strand is `+` / `-` (occasionally `.`).
- **Identity** values are fractions in `[0, 1]`, gap-collapsed BLAST-style (see §3.3).
- A `Gotcha:` callout marks behaviour that is subtle, order-dependent, or load-bearing for the
  byte-identity contract — a naive re-implementation will get these wrong.

## Byte-identity contract (the governing invariant)

LiftOn ships five orthogonal "fast-path" flags (`--stream`, `--inmemory-liftoff`,
`--threads N`/`--locus-pipeline`, `--native`) that change *how* work is scheduled and where data
lives, **never what is computed**. The defining regression gate is the **24-cell matrix**: every
combination of `--stream × --inmemory-liftoff × --threads∈{1,2,4} × --native` (= 24 cells) must
emit a **byte-identical** output GFF3. A re-implementation that reorders loci, changes attribute
serialization, or perturbs alignment by even ~0.1 % breaks this contract. The invariants that
preserve it are catalogued in §10.2.

## Quick-reference vocabulary

**Per-transcript `status` (annotation provenance)** — set on `lifton_status.annotation`:
`Liftoff` (DNA-lift kept as-is), `LiftOn_chaining_algorithm` (protein-maximized merge of Liftoff
+ miniprot), `no_ref_protein` (no reference protein to align against). Miniprot-only transcripts
(§2.4) carry their own provenance.

**Mutation vocabulary (9 types)** — recorded on the `mutation` attribute and driving ORF rescue:
`identical`, `synonymous`, `nonsynonymous`, `inframe_insertion`, `inframe_deletion`,
`frameshift`, `start_lost`, `stop_missing`, `stop_codon_gain` (plus the non-coding /
loss sentinels `non_coding`, `full_transcript_loss`, `no_protein`). Exact triggering conditions
are in §3.4; the consolidated table is Appendix C.

---
