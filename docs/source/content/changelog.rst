
|

Changelog
===========

v1.0.9
-------

Incremental release (2026-06-21). Turns on several accuracy- and
completeness-improving defaults, hardens LiftOn against whole-genome-abort
crashes on full RefSeq / cross-species genomes, and adds byte-identical
performance fast-paths plus validation tools. **Some new defaults change the
output annotation relative to v1.0.8** — each ships with an opt-out flag that
restores the previous behaviour. See the project ``CHANGELOG.md`` for the full
list.

**Changed (output-affecting defaults — results differ vs v1.0.8):**

- **Gene-like lift is now the default.** LiftOn auto-detects every reference
  top-level parent type with a transcript/exon hierarchy and lifts them all
  (pseudogenes, ``ncRNA_gene``, structured mobile elements, ...), not just
  ``gene``. This *adds* features. Pass ``--gene-only`` to restore the old
  ``gene``-only lift (``--lift-gene-like`` is a kept no-op alias).

- **Best-of-outcome Liftoff/miniprot merge is now the default.** Per transcript,
  LiftOn keeps whichever of {chained merge + ORF-rescue, Liftoff + ORF-rescue}
  yields the higher emitted protein identity, avoiding merges that could
  silently frameshift downstream CDS. Pass ``--legacy-merge`` to restore the
  pre-promotion unconditional merge (``--optimize`` is a kept no-op alias).

- **Banded / windowed alignment is now the default for all gene sizes.** The
  aligner uses anchor-windowed alignment above ~2500 aa / 8000 nt (giant genes
  are always memory-bounded, so titin-scale transcripts no longer OOM):
  identity-exact on same-species lifts, mean-neutral cross-species, much faster
  and lighter on memory. Pass ``--full-dp-align`` to restore the exact
  giant-only full-DP path (``--fast-align`` is a kept no-op alias).

- **Miniprot-only rescue is now default-ON.** When the DNA lift misses a
  reference coding gene entirely (its miniprot mRNA overlaps no lifted gene
  locus), LiftOn emits the miniprot-only model, tagged
  ``lifton_rescue=miniprot_only``. It runs as a separate pass after the main
  lift closes, gated by a protein-identity floor with a dedup guard, recovering
  genuinely-missing genes at large evolutionary distance. Pass
  ``--no-miniprot-rescue`` (or ``LIFTON_MINIPROT_RESCUE=0``) to opt out
  (``--miniprot-rescue`` is a kept no-op alias).

**Fixed (robustness / crashes):**

- **Gene-like child double-lift crash on full RefSeq genomes.** A gene-like
  feature that is a *child* of a gene was being enumerated again as a top-level
  locus, producing a duplicate FASTA key and crashing the run. Full Arabidopsis
  now completes (~99.9% of coding transcripts recovered, vs ~28% before the
  crash) and full rice likewise (~77% → ~99.9%). Top-level-only annotations are
  byte-identical.

- **Inverted-coordinate write crash.** A single malformed transcript with
  ``start > end`` (e.g. on the dog→cat lift) used to abort an entire
  ~60k-transcript genome during the write phase; such a feature is now skipped
  and logged, and the rest of the genome completes.

- **A malformed feature no longer aborts a whole genome.** The transcript writer
  now catches the project's validation exception so one bad feature is dropped
  and logged instead of propagating out of the parent write phase.

- **Whole-genome-abort hardening** around ``consume()`` / ``__str__`` plus a
  recursion-limit guard and a full traceback dump in the vendored-Liftoff call
  path, so deep/odd inputs fail loudly per-feature instead of silently killing
  the run.

**Added (performance — byte-identical fast-paths, same output, faster/lighter):**

- ``--stream`` — pipe miniprot output straight into an in-memory database,
  skipping the ``miniprot.gff3`` disk round-trip and SQLite re-ingest.
- ``--inmemory-liftoff`` — feed Liftoff's lifted features to the database
  in-process, skipping the ``liftoff.gff3`` disk write and re-ingest.
- ``--threads N --locus-pipeline`` — fan out per-locus work across a thread pool;
  output is emitted in submission order so ``--threads N`` is byte-identical to
  ``--threads 1`` (now works on the default backend without ``--native``).
- ``--native`` — route miniprot through the in-process native facade and unlock
  in-process threading; falls back gracefully if ``mappy`` is not installed.
- **Concurrent aligner step is now the default** (miniprot and Liftoff overlap);
  pass ``--serial-aligners`` to opt out (``--parallel-aligners`` is a kept no-op
  alias).
- **Fused parallel Step 7** — per-locus materialise and process phases are fused
  into one pool, lowering both wall-clock and peak memory.
- **Sequence-extraction query collapse** — collapses per-feature database
  queries (~2.2–2.4× fewer round-trips, ~25–34% faster extraction).
- **miniprot ``-t`` now scales with ``--threads``** (was pinned to its built-in
  default of 4; the default ``-t 1`` is byte-identical).

**Added (validation):**

- ``--strict-gff`` — run the NCBI GFF3 input-side validator on the reference
  annotation and exit non-zero on any spec violation.
- ``--validate-output`` / ``--validate-verbose`` — re-validate the just-written
  output GFF3 and print a structured report.
- **``gff3-validate`` console script** — a standalone GFF3 validator installed
  alongside ``lifton``.

**Packaging:**

- Python floor raised to ``>=3.10`` (the ``networkx>=3.3`` dependency requires Python ≥3.10; 3.9 is EOL).
- ``mappy`` is an optional dependency enabling the in-process ``--native`` path;
  the runtime falls back gracefully when it is absent.
- Added ``MANIFEST.in``, ``pyproject.toml`` (PEP 517/518 build config), and PyPI
  trove classifiers / project URLs.
- The vendored ``gffbase`` ships its pure-Python fallback parser (no pre-built
  ``.so``), so installs work without a Rust toolchain.

v1.0.7
-------

**Bug Fixes:**

- **Fixed gffutils UNIQUE constraint errors**: Enhanced duplicate feature ID handling with automatic recovery strategies. The system now automatically handles duplicate IDs in polished liftoff output and RefSeq annotations with non-overlapping CDS features, preventing silent failures during database creation.

- **Fixed GTF file format processing**: Added automatic GTF format detection and proper handling. GTF files are now correctly processed with automatic gene/transcript inference, and optional automatic conversion to GFF3 format using gffread or agat tools.

- **Fixed ID parsing for IDs ending with numbers**: Improved `get_ID_base()` function to safely handle feature IDs that naturally end with underscore and number (e.g., `FMUND_1`). The function now only removes suffixes when confirmed to be copy numbers, preventing silent failures.

- **Fixed CDS ID preservation**: CDS features now preserve their IDs in output files, complying with GFF3 specification. CDS features from the same mRNA can now share the same ID as required by the standard.

**Improvements:**

- **Enhanced biotype attribute support**: Added support for generic `biotype` attribute as fallback when `gene_biotype` (RefSeq) or `gene_type` (GENCODE/ENSEMBL/CHESS) are not present. This ensures protein-coding features are correctly identified regardless of annotation source.

- **Automatic GTF to GFF3 conversion**: Added optional automatic conversion of GTF files to GFF3 format for better compatibility. Conversion uses gffread (preferred) or agat tools if available, with graceful fallback to direct GTF processing.

- **Improved error handling**: Enhanced error messages and recovery strategies for database creation failures, providing better user guidance and automatic problem resolution.

- **Better format detection**: Improved file format detection logic that checks multiple lines and patterns to reliably distinguish between GTF and GFF3 formats, with GFF3 as the safe default.

v1.0.0
-------

- Initial release of LiftOn
- Release via the documentation (http://ccb.jhu.edu/lifton)
- Released via the paper (bioRxiv coming soon!)


|
|
|
|
|



.. image:: ../_images/jhu-logo-dark.png
   :alt: My Logo
   :class: logo, header-image only-light
   :align: center

.. image:: ../_images/jhu-logo-white.png
   :alt: My Logo
   :class: logo, header-image only-dark
   :align: center

