
|

Changelog
===========

v1.0.13
--------

Packaging release candidate (2026-09-16), addressing GH #78.

- Standard pip installs exclude mappy; ``lifton[native]`` or prebuilt Conda mappy
  supplies the optional, explicitly enabled experimental binding.
- Missing-aligner messages give an installation command. Pip installs Python
  packages, while fresh standard lifts require minimap2 and miniprot on PATH.
- Correct pip/macOS/source instructions and document complete Seqera environments.
- Qualify built wheel/sdist in compiler-free containers and execute fresh native
  lifts before publishing. Genuine optional binding tests remain separately tested.
- Updated Bioconda recipe preparation removes cigar and the obsolete setuptools
  cap and retains the DuckDB exclusions, prebuilt mappy and both aligners.

v1.0.12
--------

Resource-safety and diagnostics release (2026-09-14), motivated by GH #71's two
native ``SIGSEGV`` failures on an approximately 20-Gb target. It also recovers
far more genes between distantly related species and runs multi-threaded by
default.

**Accuracy between distant species:**

- A second miniprot-only rescue sub-pass reconsiders hits the rescue's length
  band rejected. That band compares genomic spans, which include introns and
  therefore follow genome size: lifting from human into fish, bird, and frog
  genomes it discarded 75–83 % of the missed genes miniprot found, although
  their alignments covered the whole reference protein. The sub-pass gates on
  protein coverage (≥ 0.8) with the same identity floor. Whole-genome
  primary-assembly gene recall rises from 0.315 to 0.593 (human → zebrafish),
  0.364 to 0.615 (human → chicken) and 0.361 to 0.638 (human → xenopus), with
  no gene lost, no duplicate model, and 99.7 % of added models overlapping an
  annotated CDS of the released target annotation. ``--no-coverage-rescue-gate``
  opts out.
- Rescued genes now carry the other reference transcripts whose miniprot hits
  lie at their locus, each held to the same identity floor, without changing
  gene placement. ``--no-rescue-isoforms`` opts out.

**Speed:**

- ``--threads N`` with N > 1 now fans out Steps 7 and 8 without
  ``--locus-pipeline``; output is byte-identical to the serial path.
  ``--no-locus-pipeline`` restores serial processing.
- Liftoff's alignment parsing and GFF3 writing are faster, and with
  ``--threads N`` its lift loop runs one forked worker per reference
  chromosome (``--no-parallel-lift`` opts out). Output is identical; the
  fresh-Liftoff aligner phase at ``-t 8`` halves on drosophila
  (261 s to 134 s).

- Miniprot-derived models now carry their stop codon. miniprot's CDS ends at
  the last aligned codon, while the reference convention -- and every other
  model LiftOn emits -- includes the stop, and the ORF search could not add it
  because such a model has no UTR to search. Only 39-59 % of rescued models on
  the distant whole genomes ended in a stop. The terminal CDS and its exon now
  grow by those three bases when the genome has them and the reference protein
  ends in a stop, after the ORF search and only when the model does not get
  worse. ``--no-orf-stop-completion`` opts out.
- ``-dir/--intermediate-dir`` puts a run's intermediate files, statistics,
  score table and manifest where you choose, so concurrent runs sharing an
  output directory no longer share one ``lifton_output/`` (GH #14).

**Changed and fixed:**

- ``-copies`` extra gene copies keep their transcripts, exons and CDS. Liftoff
  suffixes every feature of an extra copy with ``_<extra_copy_number>``, so a
  copy arrives as ``gene-X_1`` / ``rna-X_1``. LiftOn resolved the gene id back
  to the reference correctly but looked the transcript id up verbatim, so the
  lookup failed and the transcript was skipped -- taking every exon and CDS
  with it and leaving a gene line with no children. Across the 17-genome
  benchmark set this affected about 4,400 genes (rice 539 of 815 copy genes,
  human to zebrafish 1,178 of 1,881), identically in v1.0.11, so it was present
  in every release. The reference transcript is now resolved the way the gene
  already was: the exact id first, and only on a miss the ``_<N>`` copy base,
  accepted only when that base really is a transcript of the reference gene the
  copy belongs to -- so reference ids that genuinely end in ``_<int>`` are left
  alone.
- A gene emitted with no child features is counted and reported at the end of
  the run and in ``run_manifest.json`` as ``genes_emitted_without_children``.
- A flat annotation -- a prokaryotic GFF3 with top-level ``CDS`` rows and no
  ``gene`` -- now lifts those rows instead of selecting nothing and failing
  several steps later inside Liftoff, and an empty selection stops the run at
  once, naming the feature types the annotation contains (GH #37).
- ``run_manifest.json`` records where the miniprot-only rescue spends its time,
  split across candidate placement, the coverage sub-pass, and the isoform
  pass's prefetch, scoring and attachment.
- Targets above 4,000,000,000 bases run Liftoff/minimap2 before miniprot by
  default, preventing the two index-memory peaks from overlapping. Smaller
  targets retain concurrent execution. ``--serial-aligners`` and
  ``--parallel-aligners`` are mutually exclusive force overrides.
- Liftoff creates no more Python workers than alignment tasks and divides the
  configured thread budget across those tasks. One target with ``--threads 40``
  therefore creates one worker whose minimap2 command receives 40 threads.
- ``run_manifest.json`` records target statistics, schedule/reason, exact
  aligner commands, last completed stages, bounded stderr tails, return codes,
  and POSIX signal names. ``-11`` is reported as ``SIGSEGV``; OOM is not
  asserted without OS or scheduler evidence.
- miniprot older than 0.14 is rejected for a target sequence at least 2^31
  bases, covering its known large-sequence limitation.
- Trans-spliced copies with repeated logical parent IDs now bind children to
  the matching parent on the same sequence.
- Intermediate Liftoff and miniprot files are written to
  ``lifton_output/liftoff/`` and ``lifton_output/miniprot/`` again; v1.0.10
  and v1.0.11 wrote them to ``lifton_outputliftoff/`` and
  ``lifton_outputminiprot/``.

See :doc:`large_genome_resource_failures` for diagnosis and recovery guidance.

v1.0.11
--------

Single-fix release (2026-08-01). v1.0.10 could not build a database from some
annotations that earlier versions handled, aborting the run outright. If you lift
plant genomes, or any annotation produced with ``-copies``, upgrade.

**Fixed:**

- **``UNIQUE constraint failed: features.id`` on annotations containing copy
  features.** ``gffutils`` disambiguates a repeated ``cds-X`` by renaming it to
  ``cds-X_1`` — exactly the suffix Liftoff's ``-copies`` mode gives extra gene
  copies. Where that generated name already belonged to a real feature the insert
  failed and LiftOn exited. LiftOn now renames only ids that are *both* repeated
  and whose ``<id>_<n>`` already exists, leaving legitimate discontinuous CDS —
  which share an id by design — untouched. Measured across the benchmark corpus:
  v1.0.8 ingests 34 of 34 inputs, v1.0.10 ingested 31, v1.0.11 ingests 34. Same
  error class as GH #47, #12 and #7.

- **A defect in the database-build fallback is no longer reported as a bad input
  file.** ``NameError``, ``AttributeError``, ``TypeError`` and ``ImportError``
  raised inside the fallback now propagate instead of being reported as a failed
  build.

v1.0.10
--------

Incremental release (2026-07-30). Consolidates everything that accumulated
behind the v1.0.10 tag: the originally staged v1.0.10 work, a batch of
GitHub-issue fixes, a full-repository correctness audit, and a CDS-attribute
parity fix. Neither v1.0.10 nor the intermediate v1.0.10.1 was ever published,
so **v1.0.9 is the baseline for every comparison below**. As in v1.0.9, every
output-affecting default ships with an opt-out flag. See the project
``CHANGELOG.md`` for the full list.

**Changed (output-affecting defaults — results differ vs v1.0.9):**

- **A third best-of-outcome merge candidate.** Alongside {chained merge +
  ORF-rescue, Liftoff + ORF-rescue}, LiftOn now scores miniprot's *native*
  CDS-only model and adopts it only when its ORF-rescued protein identity is
  *strictly* higher than the two-way winner — so per-transcript identity never
  decreases — and never for an antisense hit. Pass ``--no-miniprot-candidate``
  to restore the two-way merge.

- **A divergence-adaptive miniprot-only rescue floor.** The rescue pass used a
  fixed protein-identity floor of 0.50; that floor is now lowered toward 0.30 as
  the DNA lift's gene recall drops, recovering more genuinely missing genes at
  large evolutionary distance while staying inert on same- and close-species
  lifts. Pass ``--no-adaptive-rescue-floor`` to restore the fixed floor.

- **Richer, spec-valid CDS rows.** Rebuilt CDS lines lost the reference's
  descriptive attributes (``Dbxref``, ``product``, ``protein_id``, ``gene``,
  ``locus_tag``, ...) and carried no ``ID`` at all. They now inherit those
  attributes and share one ``ID=cds-<transcript>`` — the spec's discontinuous-CDS
  form. Coordinates and the encoded protein are untouched; output grows 12–43%.
  ``LIFTON_NO_CDS_ATTR_CARRY=1`` and ``LIFTON_NO_CONTAINMENT_NORMALIZE=1``
  reproduce the previous bytes.

- **Coding transcripts are harmonized to ``mRNA``.** The DNA-lift path preserved
  the reference featuretype — often the generic ``transcript`` — while the
  miniprot path emitted ``mRNA``, so one output labelled the same kind of feature
  two ways. ``LIFTON_NO_MRNA_HARMONIZE=1`` preserves the reference type.

**Fixed (reported issues):**

- ``no such column: <id>`` from the DNA lift on strict SQLite builds (GH #35),
  which is why the same input worked on one machine and failed on another.
- A single unresolvable alignment killing an entire chromosome (GH #39).
- Mixed-strand output from an antisense merge (GH #33).
- UTR indels reported as frameshifts on close pairs (GH #46).
- CDS emitted with no ``ID`` (GH #32, GH #8) and coding transcripts typed
  inconsistently (GH #28).
- Header-only Liftoff output from an invalid or LFS-pointer minimap2 index
  (GH #57), and a ``--stream`` ingest crash on DuckDB 1.5.3/1.5.4 (GH #56).
- minimap2 is now preflight-checked at startup and documented as a requirement
  (GH #43, GH #11); the run summary is also written to ``stats/summary.txt``
  (GH #50).

**Fixed (robustness):**

- **A skipped locus no longer withholds the whole annotation.** One bad gene in
  a 60,000-gene genome used to produce exit 2 and only a ``*.partial.gff3``.
  Per-locus failures now publish, report ``partial_success``, and are recorded
  in ``run_manifest.json``; ``--strict-completeness`` restores the old behaviour.
- Write-funnel crashes on gene-like/organellar children and on inverted
  (``start > end``) coordinates, which aborted a whole-genome write.
- Three-level hierarchies (``gene → primary_transcript → miRNA → exon``) are
  lifted instead of silently dropped, and LiftOn no longer rejects its own valid
  Ensembl/GENCODE-shaped output.

**Added:**

- **Auditable run manifests and transactional output** —
  ``lifton_output/run_manifest.json`` records sanitized arguments, SHA-256 input
  fingerprints, tool versions, phase timings, counts, validation and failures;
  the GFF3 is staged and published atomically only on success.
- **Always-on structural output validation** before publication, recorded in the
  run manifest and not bypassable by ``--allow-partial-output``.
- **Content-addressed annotation caches** that rebuild when their manifest no
  longer matches the source, parser settings, backend, or LiftOn version.

**Performance (byte-neutral):**

- **GFF3 validation is bounded** — one contiguous top-level block at a time
  instead of the whole file: 5.64 GB → 267 MB on a full dog→cat lift (21× less
  memory, 1.9× faster), with identical reports.
- **Step 7 is ~21% faster** on *Drosophila* at ``-t 8`` (redundant leaf queries
  removed, vectorized identity counters, a pruning ordered gate), with
  byte-identical output and flat peak memory.
- **Features are cloned rather than deep-copied** (33.06 µs → 5.00 µs per rich
  CDS), and per-stage in-flight bounds (``--step7-max-inflight`` and siblings)
  cap the memory a parallel run can hold.

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
- ``--native`` — enable experimental native compatibility hooks. The mappy
  Liftoff route also requires ``LIFTON_NATIVE_LIFTOFF_ALIGN=1``; miniprot and
  bounded locus workers keep their proven paths.
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
- ``mappy`` is an optional dependency for the explicitly enabled native
  Liftoff route; the runtime falls back gracefully when it is absent.
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
- Release via the documentation (https://khchao.com/LiftOn)
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
