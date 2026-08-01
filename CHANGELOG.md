# Changelog

All notable changes to **LiftOn** are documented here. This project follows
[Keep a Changelog](https://keepachangelog.com/en/1.0.0/) conventions and
[Semantic Versioning](https://semver.org/).

## [Unreleased]

Nothing yet.

## [1.0.11] - 2026-08-01

A single-fix release. v1.0.10 could not ingest some annotations that v1.0.8
handled, aborting the run outright; if you lift plant genomes, or any annotation
produced with `-copies`, upgrade.

### Fixed

- **`UNIQUE constraint failed: features.id` on annotations containing copy
  features.** `gffutils` disambiguates a repeated `cds-X` by renaming it to
  `cds-X_1` — which is exactly the suffix Liftoff's `-copies` mode gives extra
  gene copies. When that generated name already belongs to a real feature the
  insert fails and LiftOn exits. Two ids in the rice benchmark trigger it, and
  two is enough to lose a whole genome.

  v1.0.8 recovered because its fallback renamed repeats into a `_dup<n>`
  namespace that cannot collide. v1.0.10 replaced that with a no-op to stop
  rewriting logical `ID`s — a sound goal, since renaming the segments of one
  discontinuous CDS splits a single feature in two — but it left the fallback
  identical to the strategy that had just failed.

  LiftOn now renames only ids that are *both* repeated and whose `<id>_<n>`
  already exists in the file, using facts the preflight scan already collects.
  On the rice input that is 2 ids; the other 4,677 shared CDS ids, which are
  legitimate discontinuous CDS, are left untouched. Measured across the
  benchmark corpus: v1.0.8 ingests 34 of 34 inputs, v1.0.10 ingested 31, and
  v1.0.11 is back to 34. Same error class as GH #47, #12 and #7.

- **A defect in the database-build fallback is no longer reported as a bad
  input file.** The strategy handlers catch `Exception` so a malformed
  annotation falls through to the next strategy; that breadth also swallowed
  programming errors raised inside the fallback itself, which then surfaced as
  "DB build failed" with no way to tell the two apart. `NameError`,
  `AttributeError`, `TypeError` and `ImportError` now propagate.

## [1.0.10] - 2026-07-30

This release consolidates everything that accumulated behind the v1.0.10 tag: the
originally staged v1.0.10 work, the GitHub-issue batch that had been prepared as
v1.0.10.1, a full-repository correctness audit, and a CDS-attribute parity fix.
Neither v1.0.10 nor v1.0.10.1 was ever published, so **v1.0.9 is the baseline for
every comparison below**.

The themes are dual evidence used more aggressively (a third merge candidate, a
divergence-adaptive rescue floor), output that is complete and spec-valid enough to
hand to downstream tools unedited, whole-genome runs that finish instead of aborting
on one bad locus, and a large byte-neutral performance and memory reduction. As in
v1.0.9, every output-affecting default ships with an opt-out flag, and the 24-cell
byte-identity matrix stays green.

### Changed (output-affecting defaults — results differ vs v1.0.9)

- **Third best-of-outcome merge candidate ("candidate-3"), default-ON.** In
  addition to {chained merge + ORF-rescue, Liftoff + ORF-rescue}, LiftOn now
  evaluates a third candidate that emits miniprot's *native* CDS-only model, and
  adopts it only when its ORF-rescued protein identity is *strictly* higher than
  the two-way winner. The candidate is built cleanly (fresh exon/CDS scaffold,
  not routed through `update_cds_list`) and is gated by a strand-consistency
  guard so an antisense miniprot hit can never replace a lifted transcript. This
  is output-additive on the transcripts it improves and is a strict
  per-transcript non-regression. Pass `--no-miniprot-candidate` (or set
  `LIFTON_MINIPROT_CANDIDATE=0`) to restore the two-way merge.
  (`--miniprot-candidate` is a kept no-op alias.)
- **Divergence-adaptive miniprot-only rescue floor, default-ON.** The
  miniprot-only rescue pass (v1.0.9) used a fixed protein-identity floor of 0.50.
  v1.0.10 lowers that floor toward 0.30 as the DNA lift's gene recall drops
  (inert on high-recall same/close-species lifts), recovering more genuinely
  missing genes at large evolutionary distance while remaining redundancy-free
  (`off ⊆ on`, 0 lost / 0 redundant by construction). Pass
  `--no-adaptive-rescue-floor` (or `LIFTON_RESCUE_ADAPTIVE_FLOOR=0`) to restore
  the fixed 0.50 floor. (`--adaptive-rescue-floor` is a kept no-op alias.)
- **Every CDS row now carries the reference's descriptive attributes.** The
  chaining merge and the ORF-rescue boundary patch rebuild a transcript's CDS
  list and used to reset the emitted attributes to `{Parent}`, so an output
  described some CDS richly and most not at all: on the full drosophila lift,
  7,127 of 7,142 pure-Liftoff CDS rows carried
  `Dbxref`/`Name`/`gbkey`/`gene`/`locus_tag`/`product`/`protein_id`/
  `orig_transcript_id` while only 21 of 151,277 merged rows did — even though
  all 28,051 merged mRNA rows kept those same attributes. A rebuilt CDS now
  inherits its transcript's reference CDS attributes. `ID` remains owned by the
  write-funnel allocator and `Parent` is still rewritten, so the change only
  adds attributes. Across the eight-dataset divergence ladder, CDS coverage
  rises from as low as 0.06 to 1.0 with columns 1–8 byte-identical on every
  row, no attribute lost or rewritten, `protein_identity` unchanged on every
  transcript, and validity flat; output grows 12–43%.
  `LIFTON_NO_CDS_ATTR_CARRY=1` reproduces the previous bytes.

  Confirmed on a full dog→cat genome (1,815,199 rows): **692,520 CDS rows gain
  their descriptive attributes**, CDS coverage reaches 739,862/739,862, and
  columns 1–8 are byte-identical on every row — so `protein_identity` is
  unchanged across all 66,276 transcripts and `gff3-validate` reports the same
  0 errors / 176 warnings in both arms. Output grows 398 MB → 524 MB (+31.5 %)
  and peak RSS by 2.5 % (2.90 → 2.98 GB).
- **Every CDS now carries an `ID` (GH #32, GH #8).** Rebuilt CDS lines were
  emitted with no `ID=` at all. A transcript's CDS segments now share one
  `ID=cds-<transcript>` — the spec's discontinuous-CDS form (same ID, same
  Parent, same type). Only the attribute is added; coordinates and the encoded
  protein are untouched. `LIFTON_NO_CONTAINMENT_NORMALIZE=1` reproduces the
  previous bytes.
- **Coding transcripts are harmonized to `mRNA` (GH #28).** The Liftoff and
  chaining paths preserved the reference featuretype — often the generic
  `transcript` — while the miniprot path emitted `mRNA`, so one output labelled
  the same kind of feature two different ways depending on which evidence won.
  A transcript that has a CDS is an mRNA by SO definition, so coding transcripts
  are now canonicalized to `mRNA`; non-coding transcripts keep their type. This
  also corrects the coding/other split in the statistics report, which keys on
  that featuretype. `LIFTON_NO_MRNA_HARMONIZE=1` preserves the reference type.

### Added

- **Auditable run manifests and transactional output.** Every completed CLI
  run writes `lifton_output/run_manifest.json` with sanitized arguments,
  SHA-256 input fingerprints, software/tool versions, backend and cache
  choices, phase timings, resource usage, counts, validation, and failures.
  GFF3 output is staged and atomically published only after success; failed
  data is preserved as `*.partial.gff3`. Partial publication requires the
  explicit `--allow-partial-output` option.
- **Always-on structural output validation.** Before publication, LiftOn now
  streams over every staged GFF3 and rejects malformed columns, coordinates,
  strands, CDS phases, IDs, or parent relationships. This mandatory gate is
  recorded as `gff3_structural` in the run manifest and cannot be bypassed by
  `--allow-partial-output`; `--validate-output` retains the deeper optional
  project-specific checks.
- **Content-addressed annotation caches.** gffutils and gffbase databases now
  require a matching manifest covering source content, parser/inference
  settings, backend/schema, and LiftOn version. Stale or corrupt caches rebuild
  under a temporary name and publish atomically.
- **A minimap2 preflight (GH #43, GH #11).** minimap2 is a hard runtime
  dependency of the default DNA lift but, unlike miniprot, was neither
  documented as one nor checked, so a missing binary surfaced as a deep
  `FileNotFoundError` mid-run. It is now checked at startup with a clear
  message and listed in the installation requirements.
- **The run summary is written to `stats/summary.txt` (GH #50).** The
  per-category side-car files existed but the summary block itself went only to
  stderr, so it was lost from a redirected or scripted run.
- **Coordinate-independent structural evaluation metrics** in the benchmark
  evaluator — per-transcript intron-chain exactness, exon sensitivity/precision
  (Sn·Sp), and ORF validity (start `M` / stop `*` / no internal stop) — scored on
  spliced 5′→3′ boundary positions, catching structural regressions that
  protein-identity alone hides.
- **Guarded benchmark build controller.** Maintainers can run resumable gate,
  subset, full, and end-to-end stages in isolated tmux cells with frozen
  provenance, resource admission checks, artifact validation, and explicit
  retry/reconciliation commands. Controller runs never update the canonical
  benchmark baseline.

### Fixed (reported issues)

- **`pip install lifton` failed in any clean environment.** LiftOn depended on
  `cigar` 0.1.3 — unmaintained since 2015, distributed as an sdist only, whose
  `setup.py` calls `ez_setup.use_setuptools()` to *download* a setuptools egg
  from a URL that no longer serves one. Every fresh install without a cached
  wheel therefore died with
  `tarfile.ReadError: file could not be opened successfully: not a gzip file`.
  LiftOn used exactly one call from that package; it is now parsed in-tree
  (`coreutils.parse_cigar_items`, byte-for-byte equivalent and without the
  upstream generator's PEP 479 `raise StopIteration` bug) and the dependency is
  removed. This also retires the `setuptools<81` build cap, which existed only
  to work around the same package. Affected every prior release, not just this
  one.
- **`no such column: <id>` from the DNA lift on some machines (GH #35).** The
  vendored Liftoff built its `id IN (...)` clauses by double-quoting each value,
  which SQL reads as an *identifier*, not a string. Modern SQLite rejects it;
  older builds silently accepted it through the DQS misfeature — hence the same
  input working on one machine and failing on another. IDs are now bound as
  parameters, which is correct on every SQLite build, injection-safe, and
  handles IDs containing quotes.
- **One bad alignment no longer kills a whole chromosome (GH #39).** On GRCz11
  zebrafish, a gene-family copy whose query name converted to an id absent from
  the parent hierarchy raised `KeyError` from an unguarded lookup, taking the
  entire chromosome with it while every other chromosome lifted normally. That
  alignment is now dropped with a warning naming the feature, and the rest of
  the chromosome lifts.
- **Strand-inconsistent output from an antisense merge (GH #33).** The
  candidate-3 path already refused an opposite-strand miniprot hit, but the
  best-of-outcome chaining merge did not, so miniprot CDS could be grafted onto
  a Liftoff transcript on the other strand — producing exon and CDS rows with
  opposite strands whose better-looking ORF score was a wrong-frame artifact.
  The guard now covers both paths; such transcripts fall back to pure
  Liftoff + ORF-rescue.
- **UTR indels are no longer reported as frameshifts (GH #46).** Frameshift
  detection ran on the full-transcript alignment, so a 1–2 bp indel in a UTR —
  ubiquitous on close pairs — was scored as a frameshift in `score.txt` and
  triggered a needless ORF re-search that could perturb the emitted model.
  Detection is now scoped to the coding sub-alignment.
- **Header-only Liftoff output from an invalid minimap2 index (GH #57).**
  LiftOn now validates cached `.mmi` files, detects unresolved Git LFS
  pointers and corrupt or stale indexes, and rebuilds them in the run's
  `lifton_output/intermediate_files/` directory. Minimap2 failures and
  malformed SAM headers are reported at the alignment stage instead of later
  appearing as an annotation-format error. Generated minimap2 indexes are no
  longer stored in Git LFS.
- **`--stream` miniprot ingest crash (GH #56).** DuckDB 1.5.3 and 1.5.4 can
  fail while appending gffbase's internal spatial `GEOMETRY` column once a
  large GFF3 crosses a row-group boundary. Those releases are excluded, and
  existing environments automatically use the equivalent B-tree query path.
  Streamed blobs also retry once on a fresh connection without geometry if
  another DuckDB build raises the same internal error. Users can force this
  safe path with `LIFTON_DISABLE_RTREE=1`.

### Fixed (robustness / correctness)

- **A skipped locus no longer withholds the whole annotation.** Any per-locus
  processing or serialization failure previously blocked publication, so one
  bad gene in a 60,000-gene genome produced exit 2 and only a
  `*.partial.gff3`. Per-locus failures now publish, report the run as
  `partial_success`, and are recorded in `run_manifest.json`; publication is
  blocked only by failures that mean the output itself cannot be trusted
  (structural validation, report writing, or a non-empty reference that
  produced no features). `--strict-completeness` restores the previous
  behaviour and `--allow-partial-output` still overrides.
- **Gene-like-feature write-funnel crash** on genomes with organellar / gene-like
  children (`LiftOn_FEATURE.normalize_containment` added), and **inverted-coordinate
  write crash** — a malformed (`start > end`) transcript now drops with a logged
  warning instead of aborting the whole-genome write phase.
- **Write-time containment normalization** eliminates the duplicate-exon-ID and
  CDS-past-exon errors LiftOn's chaining could introduce (the reference has none);
  CDS coordinates are never touched, so protein identity is unchanged
  (`LIFTON_NO_CONTAINMENT_NORMALIZE=1` restores the pre-normalization bytes).
- **LiftOn no longer rejects its own valid output.** The mandatory structural
  validator recognised transcript levels through a closed allowlist while the
  lift auto-detects gene-like types and preserves the reference's middle level
  verbatim, so an Ensembl/GENCODE `pseudogene → pseudogenic_transcript → exon`
  (and RefSeq `misc_RNA`, `pre_miRNA`, `pseudogenic_tRNA`) was lifted
  correctly and then failed validation. Transcript levels are now recognised
  by shape; genuinely wrong hierarchies are still rejected.
- **Three-level hierarchies are lifted instead of dropped.** A reference shaped
  `gene → primary_transcript → miRNA → exon` silently lost the nested
  transcript and all of its exons and CDS.
- **An unusable miniprot alignment no longer drops a whole gene.** A miniprot
  mRNA with no CDS children produced a `None` alignment that was dereferenced,
  removing the entire Liftoff locus even though its DNA lift was sound. The
  chaining whole-side fallbacks now also emit CDS in the order
  `update_cds_list` expects.
- **Gene-like child double-lift fix** (a `ncRNA`/`pseudogene` that is also a child
  of a coding gene is no longer enumerated twice), **deterministic `-copies`
  numbering**, **symmetric empty-Liftoff chaining**, **`update_cds_list` root exon
  sort**, and a **loud failure** when Step-3 extraction produces an empty
  `transcripts.fa` (reference-genome seqids not matching the annotation).
- **Copy genes no longer get doubly-suffixed CDS identifiers.** Every chained
  CDS of a transcript shared one attribute dict, and `Lifton_CDS.__init__`
  reads `extra_copy_number` and then pops it — so only the first segment had its
  copy suffix stripped, the rest kept it, and the write funnel then appended the
  suffix again, emitting `cds-NP_9.1_1_1`. Annotations without extra gene copies
  are unaffected.
- **`--legacy-merge` no longer duplicates exons.** A single-CDS transcript
  re-emitted the CDS-bearing exon for every downstream exon (a 7-exon
  transcript emitted 12).
- **Rejected rescue candidates no longer burn stable CDS IDs**, which made
  output depend on how many candidates happened to be rejected.
- **A CDS spanning several exons is reported instead of silently duplicated,**
  and each exon now receives its own copy rather than all of them sharing one
  object that a later edit could rewrite.
- **The mandatory output validator now reports the true error count.**
  `validate_gff3_structure` never populated `issue_totals`/`severity_totals`,
  so the run manifest and the failure message reported at most 50 errors per
  check however broken the file was — and in counts-only mode
  (`max_issues_per_check=0`) `is_valid` returned True for a file with real
  errors. The deep `gff3-validate` CLI keeps its documented per-check cap.
- **Staged output is fsynced before it is published.** The pipeline closed the
  staged GFF3 with a bare `close()`, so `OutputTransaction.close()`
  short-circuited and skipped its flush and `os.fsync`; `commit()` could then
  publish, via `os.replace`, data that had never reached disk.
- **Output validation exit status and directive provenance.** A validation
  failure now exits non-zero and cannot replace an existing output. Reference
  assembly/species/coordinate directives are no longer copied onto target
  coordinates.
- **Deterministic, isolated locus state.** Parallel Step 7 buffers copy/tree,
  score, and chain side effects per locus and commits them in submission order;
  miniprot candidates do not consume copy IDs or suppression intervals before
  acceptance. GFF3 hierarchies are buffered so malformed children cannot leave
  half-written transcript blocks, and statistics count emitted models only.
- **gffbase root scans are totally ordered.** `_order_clause(None)` sorted by
  `file_order` alone, which rows synthesised by parent inference do not carry,
  so the serial and parallel passes could assign different submission indices
  to the same input.
- **The miniprot child proxy honours the requested ordering.** It served the one
  ordering its cache was built with whatever the caller asked for, so a request
  for file order — what the real backends return — would have been answered with
  start-sorted rows. Latent: every current caller asks for `start`.
- **A read-only reference directory no longer aborts the run.**
  `annotation_cache.cache_lock` created `<db>.lock` next to the database and
  raised an uncaught `PermissionError` on a shared read-only mount; it now
  falls back to a temp-directory lock and, failing that, proceeds unlocked with
  a warning.
- **Cached `-L`/`-M` runs work without miniprot installed** — the preflight ran
  before the short-circuit. `check_miniprot_installed` no longer leaks a
  version banner into `-o stdout` output and now checks the exit status.
- **A missing binary or typo'd path reports immediately** instead of appearing
  to hang through a multi-gigabyte input fingerprint.
- **`gff3-validate`** false positives corrected: a spec-valid discontinuous CDS
  (multiple segments sharing one ID) is exempted, and `strand '?'` (stranded but
  unknown) is accepted.

### Performance (byte-neutral)

- **GFF3 validation is bounded instead of whole-file.** `validate_gff3_file`,
  reachable from `--validate-output` and the `gff3-validate` CLI, materialised a
  record per row and walked that list five times: **406 MB** for a 92,675-row
  output and **5.64 GB** for a full dog→cat lift. It now validates one
  contiguous top-level block at a time — **39 MB** and **267 MB** respectively,
  21× less memory at mammalian scale and 1.9× faster, with reports identical
  down to issue order and message text. Files that are not block-ordered, or
  that contain a genuine duplicate ID, fall back to the previous path unchanged.
- **Faster per-locus processing.** The scalar materialiser no longer walks
  exon/CDS leaves (millions of redundant SQLite queries on large genomes), the
  per-alignment identity counters are vectorized (2.4–2.8× protein, 1.6–1.8×
  DNA), and the Step-7 ordered gate prunes committed intervals instead of
  rescanning an ever-growing list under a single lock. Measured on
  *Drosophila* at `-t 8` over two alternated repeats: **Step 7 −21 %**
  (161.6 s → 127.4 s), peak RSS flat, byte-identical output. Verified
  identity-neutral on a mammalian dog→cat lift (mean protein-identity
  delta 0.00000000 over 4,288 transcripts).
- **Features are cloned, not deep-copied.** `copy.deepcopy` ranked second in
  both Step-7 profiles (36.2 M calls on *Drosophila*), and 43 % of a single
  `Feature` copy recursed through attribute strings while 24 % re-copied the
  `dialect` dict — one object every feature of a database already shares and
  nothing in LiftOn writes to. `coreutils.clone_feature` / `clone_attributes`
  copy exactly what is mutable and share the rest: **33.06 µs → 5.00 µs** per
  rich CDS. `Lifton_EXON` and `Lifton_CDS` gain `__deepcopy__` so the hottest
  caller speeds up without changing its call site. Byte-identical output.
- **Less wasted work in the Step-7 hot loop.** A new `LIFTON_PROFILE_STEP7`
  cProfile gate ranked the per-locus phase, and four provably equivalent
  rewrites followed: the parasail alphabet check (71.3 M generator calls for a
  test that almost always answers "clean") now decides by set membership;
  `_coding_subalignment` locates two boundaries and slices instead of building
  two character lists; `__find_orfs` drops a write-only accumulator and a
  quadratic re-scan; and `LiftOn_translate` no longer fetches a coding sequence
  it discards on the next line. `get_truncated_protein` gains a counting
  companion so a dictionary of pyfaidx records is not built for its `len()`.
  Output is **byte-identical** on the drosophila and dog→cat subsets. Step-7
  wall over six alternated repeats against the previous commit: dog→cat
  **198.5 s → 187.7 s (−5.4 %, faster in 6 of 6 repeats)**; on drosophila the
  −1.9 % mean sits inside that benchmark's own 8.9 % run-to-run spread, so no
  claim is made there.
- **One fewer database query per miniprot candidate.** The same
  `children(..., ('CDS','stop_codon'))` statement was issued twice per
  candidate — 9.6 % of every SQL statement Step 7 ran on *Drosophila*
  (72,913 → 65,901). Output is byte-identical.
- **Bounded locus scheduling and direct miniprot ingest.** Parallel Step 7
  now keeps at most `2 * --threads` submitted-but-not-emitted loci by default
  (`--step7-max-inflight` overrides it) and lazily materializes full payloads.
  Parallel Step 8 and evaluation use the same ordered bound, configurable with
  `--step8-max-inflight` and `--evaluation-max-inflight`. `--stream` now parses
  miniprot stdout incrementally into a staged, checkpointed DuckDB database,
  avoiding both `miniprot.gff3` and a second in-memory hit graph.
- **Vendored-Liftoff speed backports.** `seperate_parents_and_children`
  (hierarchy build) now selects only the parent/child/feature id columns rather
  than `SELECT *`, and `convert_all_children_coords` (child-coordinate conversion)
  computes relative coordinates once, hoists invariants, early-breaks, and avoids
  a copy. Both are byte-neutral on a fresh Liftoff lift (`LIFTON_LEGACY_HIERARCHY=1`
  / `LIFTON_LEGACY_CONVERT=1` restore the prior code paths).

### Development

- `make test-fault` runs deterministic injected-failure coverage;
  `make test-fault-stress` opts into the longer repeated stress cases.
- The `cigar` dependency is gone, so the `setuptools<81` build cap it required
  is gone with it, and `tests/test_packaging_metadata.py` now asserts `cigar`
  stays out of the requirements instead.

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
- `--native` — enable experimental native compatibility hooks. The mappy
  Liftoff route additionally requires `LIFTON_NATIVE_LIFTOFF_ALIGN=1`; the
  guarded miniprot and bounded-locus paths do not require this flag.
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

- **Python floor raised to `>=3.10`** (the `networkx>=3.3` dependency requires
  Python ≥3.10, and 3.9 reached EOL in 2025-10).
- **`mappy` is an optional dependency** for the explicitly enabled native
  Liftoff route; the runtime falls back gracefully when it is absent.
- Added a `MANIFEST.in`, a `pyproject.toml` (PEP 517/518 build configuration),
  and PyPI trove classifiers / project URLs.
- The vendored `gffbase` ships its pure-Python fallback parser (no pre-built
  `.so`), so installs work without a Rust toolchain.

## [1.0.8]

Prior release. See the project documentation and the git history for details.

[1.0.10]: https://github.com/Kuanhao-Chao/LiftOn/releases/tag/v1.0.10
[1.0.9]: https://github.com/Kuanhao-Chao/LiftOn/releases/tag/v1.0.9
