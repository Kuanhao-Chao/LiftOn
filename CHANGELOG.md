# Changelog

All notable changes to **LiftOn** are documented here. This project follows
[Keep a Changelog](https://keepachangelog.com/en/1.0.0/) conventions and
[Semantic Versioning](https://semver.org/).

## [1.0.9] - 2026-06-21

This is an incremental release that turns on several accuracy- and
completeness-improving defaults, hardens LiftOn against whole-genome-abort
crashes seen on full RefSeq/cross-species genomes, and adds a family of
byte-identical performance fast-paths and validation tools. **Some new defaults
change the output annotation relative to v1.0.8** — every such change ships with
an opt-out flag that restores the previous behaviour (see *Changed* below).

### Changed (output-affecting defaults — results differ vs v1.0.8)

- **Gene-like lift is now the default.** LiftOn now auto-detects *every*
  reference top-level parent type that has a transcript/exon hierarchy and lifts
  them all — pseudogenes, `ncRNA_gene`s, structured mobile elements, etc. — not
  just `gene`. This *adds* features to the output. Pass `--gene-only` to restore
  the old `gene`-only lift. (`--lift-gene-like` is a kept no-op alias, since this
  is now the default.) An explicit `-f/--features` always overrides the
  auto-detection.
- **Best-of-outcome Liftoff↔miniprot merge is now the default.** Per transcript,
  LiftOn keeps whichever of {chained merge + ORF-rescue, Liftoff + ORF-rescue}
  yields the higher emitted protein identity, instead of applying the chained
  CDS unconditionally. This avoids merges that could silently frameshift
  downstream CDS on divergent inputs. Pass `--legacy-merge` to restore the
  pre-promotion unconditional merge (the published-manuscript behaviour).
  (`--optimize` is a kept no-op alias.)
- **Banded / windowed alignment is now the default for all gene sizes.** The
  protein/DNA aligner uses anchor-windowed alignment above ~2500 aa / 8000 nt
  (giant genes are always memory-bounded, so titin-scale transcripts no longer
  OOM). This is identity-exact on same-species lifts and mean-neutral
  cross-species, while being substantially faster and far lighter on memory.
  Pass `--full-dp-align` to restore the exact giant-only full-DP path.
  (`--fast-align` is a kept no-op alias.)
- **Miniprot-only rescue is now default-ON.** When the DNA lift misses a
  reference coding gene *entirely* (its miniprot mRNA overlaps no lifted gene
  locus), LiftOn now emits the miniprot-only model, tagged
  `lifton_rescue=miniprot_only`. The rescue runs as a separate pass after the
  main lift closes, gated by a protein-identity floor rather than the tight
  miniprot length-ratio band, with a dedup guard so it never produces redundant
  or lost models. It recovers genuinely-missing genes at large evolutionary
  distance. Pass `--no-miniprot-rescue` (or set `LIFTON_MINIPROT_RESCUE=0`) to
  restore the pre-1.0.9 lift. (`--miniprot-rescue` is a kept no-op alias.)

### Fixed (robustness / crashes)

- **Gene-like child double-lift crash on full RefSeq genomes.** A gene-like
  feature (e.g. an `ncRNA`/`pseudogene`) that is a *child* of a gene was being
  enumerated a second time as a top-level locus, producing a duplicate FASTA key
  and crashing the run. Full **Arabidopsis** now completes (~99.9% of coding
  transcripts recovered, vs ~28% before the crash) and full **rice** likewise
  (~77% → ~99.9%). Top-level-only annotations are byte-identical.
- **Inverted-coordinate write crash.** A single malformed transcript with
  `start > end` (seen e.g. on the dog→cat lift) used to abort the entire
  ~60k-transcript genome during the write phase. Such a feature is now skipped
  and logged, and the rest of the genome completes.
- **A malformed feature no longer aborts a whole genome.** The transcript writer
  now catches the project's validation exception so one bad feature is dropped
  and logged instead of propagating out of the parent write phase.
- **Whole-genome-abort hardening.** Robustness fixes around `consume()` /
  `__str__` and a recursion-limit guard (`sys.setrecursionlimit`) plus a full
  traceback dump in the vendored-Liftoff call path, so deep/odd inputs fail
  loudly per-feature instead of silently killing the run.
- **GFF3 parent-child containment + unique-ID guarantees (output now validates
  clean).** On frameshift-corrected RefSeq models the chaining / ORF-boundary
  patching could leave a CDS extending a few bp past its exon (so the
  CDS/exon fell outside the mRNA span) and could reuse stale reference
  `exon-<acc>-N` IDs (duplicate exon IDs across the rebuilt structure) — invalid
  GFF3 the reference itself does not have. LiftOn now normalises containment at
  write time: each exon is widened to cover its CDS, exons are sorted and (only
  on a real collision) re-numbered 5′→3′, and the transcript/gene span is set to
  the child envelope. This is **boundary/ID-only — coding sequence and protein
  identity are unchanged** — and a no-op on already-valid transcripts. The
  bundled `gff3-validate` / `--validate-output` now reports **zero errors** on
  the benchmark genomes (e.g. full dog→cat went from ~190 containment + ~28k
  duplicate-ID errors to a clean *VALID* verdict). Set
  `LIFTON_NO_CONTAINMENT_NORMALIZE=1` to reproduce the pre-1.0.9 bytes. The
  validator was also corrected to accept a spec-valid discontinuous CDS (the
  multiple segments of a multi-exon CDS legitimately share one ID).

### Added (performance — byte-identical fast-paths, same output)

These flags optimise wall-clock or memory and are **byte-identical to the
default output** (pinned by the 24-cell byte-identity matrix test):

- `--stream` — pipe miniprot output straight into an in-memory database,
  skipping the `miniprot.gff3` disk round-trip and SQLite re-ingest.
- `--inmemory-liftoff` — feed Liftoff's lifted features to the database
  in-process, skipping the `liftoff.gff3` disk write and re-ingest.
- `--threads N --locus-pipeline` — fan out the per-locus alignment work across a
  thread pool; output is emitted in submission order so `--threads N` is
  byte-identical to `--threads 1`. Now works on the default backend without
  `--native`.
- `--native` — route miniprot through the in-process native facade and unlock
  in-process threading; falls back gracefully if `mappy` is not installed.
- **Concurrent aligner step is now the default** — miniprot (subprocess) and
  Liftoff (DNA) now overlap, collapsing that step's wall-clock to the larger of
  the two. Pass `--serial-aligners` to opt out. (`--parallel-aligners` is a kept
  no-op alias.)
- **Fused parallel Step 7** — the per-locus materialise and process phases are
  fused into one pool, lowering both wall-clock and peak memory.
- **Sequence-extraction query collapse** — collapses the per-feature database
  queries during sequence extraction, reducing round-trips ~2.2–2.4× and
  speeding up that step ~25–34%.
- **miniprot `-t` now scales with `--threads`** — miniprot was previously pinned
  to its built-in default of 4 threads regardless of `-t/--threads`; it now
  scales with LiftOn's `--threads` (the default `-t 1` is byte-identical, since
  no `-t` is emitted there).

### Added (validation)

- `--strict-gff` — run the NCBI GFF3 input-side validator on the reference
  annotation and exit non-zero on any spec violation (missing
  `##gff-version 3`, `start>end`, negative coordinates, unencoded reserved
  characters, dangling `Parent`, etc.).
- `--validate-output` (and `--validate-verbose`) — re-validate the just-written
  output GFF3 (hierarchy / CDS phase / containment / LiftOn-attribute checks)
  and print a structured report.
- **`gff3-validate` console script** — a standalone GFF3 validator installed
  alongside `lifton` (`gff3-validate path/to/out.gff3`).

### Packaging

- **Python floor raised to `>=3.9`** (drops EOL Pythons 3.6–3.8; the code uses
  PEP-585 builtin generics, native in 3.9).
- **`mappy` is an optional dependency** that enables the in-process `--native`
  path; the runtime falls back gracefully when it is absent.
- Added a `MANIFEST.in`, a `pyproject.toml` (PEP 517/518 build configuration),
  and PyPI trove classifiers / project URLs.
- The vendored `gffbase` ships its pure-Python fallback parser (no pre-built
  `.so`), so installs work without a Rust toolchain.

## [1.0.8]

Prior release. See the project documentation and the git history for details.

[1.0.9]: https://github.com/Kuanhao-Chao/LiftOn/releases/tag/v1.0.9
