## 5. Annotation Backend & Feature Database

This section specifies how LiftOn ingests a reference (or intermediate) GFF3/GTF annotation into a queryable feature database, how it partitions features into protein-coding vs non-coding work lists, and the small but load-bearing library of ID-parsing, ID-mapping, overlap, and per-locus-alignment helpers that the rest of the pipeline relies on. Three source files are covered: `lifton/annotation.py` (the `Annotation` class and its DB-build machinery), `lifton/gffbase_adapter.py` (the DuckDB shim), and `lifton/lifton_utils.py` (the helper library).

The reader must be able to reimplement all of §5.1–§5.5 from this text alone. The internals of the gffbase Rust/DuckDB backend itself are deferred to §5.6.

---

### 5.1 The `Annotation` class

`Annotation` (`lifton/annotation.py:33`) wraps one feature database built from one annotation file (or one in-memory GFF3 blob). It is the only object LiftOn uses to query the reference annotation, the lifted-Liftoff GFF, and the miniprot GFF. The exact same class serves all three, distinguished only by what was passed in.

#### 5.1.1 Constructor signature and the polymorphic `source` argument

```python
def __init__(self, source, infer_genes, infer_transcripts,
             merge_strategy="create_unique", id_spec=None,
             force=False, verbose=False, auto_convert_gtf=True,
             *, backend=None):
```
(`lifton/annotation.py:44-56`)

| Parameter | Type | Default | Meaning |
|---|---|---|---|
| `source` | `str` \| `os.PathLike` \| `bytes` \| `bytearray` | — | Either a filesystem path to a GFF3/GTF/SQLite-DB file, **or** an in-memory GFF3 blob (bytes, or a str that "looks like" GFF3). |
| `infer_genes` | `bool` | — | Passed to the DB builder as `disable_infer_genes = not infer_genes`. GTF input forces this to `True` (see §5.2). |
| `infer_transcripts` | `bool` | — | Passed as `disable_infer_transcripts = not infer_transcripts`. GTF input forces `True`. |
| `merge_strategy` | `str` | `"create_unique"` | gffutils/gffbase merge strategy for duplicate IDs. |
| `id_spec` | any | `None` | gffutils `id_spec` (which attribute to use as the feature ID). |
| `force` | `bool` | `False` | Overwrite an existing DB file when building. Note: the fallback strategies (§5.1.5) always pass `force=True` regardless of this. |
| `verbose` | `bool` | `False` | Emit success/warning log lines. |
| `auto_convert_gtf` | `bool` | `True` | Attempt GTF→GFF3 conversion via gffread/agat (see §5.2). |
| `backend` (kw-only) | `Optional[str]` | `None` | Force `"gffutils"` or `"gffbase"`; `None` ⇒ env/default resolution. |

**Blob detection (`annotation.py:63-68`).** Immediately after entry:
1. Import `gffbase_adapter` as `_adapter`.
2. `self._is_blob = isinstance(source, (bytes, bytearray)) or _adapter.looks_like_gff3_blob(source)`.
   - `looks_like_gff3_blob` (`gffbase_adapter.py:114-122`) returns `True` only for `bytes`/`bytearray` whose first 4096 bytes contain a tab (`b"\t"`) **or** start with `b"#"`. For any non-bytes input (e.g. a `str` path) it returns `False`. **Gotcha:** a `str` GFF3 blob is therefore *not* detected by `looks_like_gff3_blob`; the only way a `str` reaches the blob branch is if `isinstance(source, (bytes, bytearray))` — which it isn't — so a `str` source is always treated as a path. The comment at `annotation.py:66-68` ("str-shaped GFF3 blob — coerce to bytes") is defensive but currently unreachable for plain `str`; in practice the streaming caller passes `bytes`.
3. If `self._is_blob` and `source` is not already bytes, `source = source.encode("utf-8")`.

**Backend resolution (`annotation.py:70-75`).** `self.backend = self._resolve_backend(backend, is_blob=self._is_blob)`. If the input is a blob but the resolved backend is not `"gffbase"`, raise `ValueError("In-memory GFF3 blob input requires backend='gffbase'; the legacy gffutils backend can only accept a file path.")`.

#### 5.1.2 Backend selection — `_resolve_backend`

`@staticmethod _resolve_backend(explicit, *, is_blob)` (`annotation.py:149-181`). Priority order, first match wins:

1. If `explicit in ("gffutils", "gffbase")` → return `explicit`.
2. If `is_blob` → return `"gffbase"` (gffutils has no string-ingest path).
3. If the environment variable `LIFTON_USE_GFFBASE` is set to any truthy value (`os.environ.get("LIFTON_USE_GFFBASE")` is non-empty/non-`None`) → return `"gffbase"`.
4. Otherwise → return `"gffutils"` (the default).

**Gotcha (performance, not byte-output):** The default deliberately stays `gffutils`. A prior attempt (Phase 17a-1) to flip the global default to gffbase caused a ~30–50× slowdown of the Step-3 sequence extractor on real RefSeq-scale inputs because gffbase's per-feature `children(...)` query path is far slower than gffutils' SQLite B-tree for that access pattern. gffbase is opt-in (env var or kwarg) to unlock real threading; the benchmark wrapper `benchmarks/phase16_rerun_bee.sh` sets `LIFTON_USE_GFFBASE=1` so the parallel path is exercised.

#### 5.1.3 Directive capture

`self.directives: list[str] = []` is initialised before any DB work. For the path (non-blob) branch only, `self._capture_directives(self.file_name)` is called (`annotation.py:90-92`).

`_capture_directives(path)` (`annotation.py:716-738`) — Phase 15a, preserves GFF3 header directives in input order so the output writer can re-emit them:
1. Open `path` with `encoding="utf-8", errors="replace"`. On `OSError`, return silently (no directives).
2. Iterate lines; strip trailing `\r\n` to `line`.
3. Skip empty lines (`continue`).
4. If `line.startswith("##")` or `line.startswith("#!")`: append the full line to `self.directives`, continue.
5. If `line.startswith("#")` (a bare comment, single `#`): `continue` (keep scanning, do **not** store).
6. The first non-`#` line is a data line — `return` immediately (header is over).

**Gotcha:** Directives appearing *after* the first feature row are not captured (the scan stops at the first data line). For the blob branch `self.directives` stays empty.

#### 5.1.4 Constructor control flow after setup

Recorded fields (path branch, `annotation.py:79-87`): `self.file_name` (= `source` for paths, or the literal string `"<in-memory blob>"` for blobs), `self.merge_strategy`, `self.id_spec`, `self.force`, `self.verbose`, `self.auto_convert_gtf`.

**Blob branch (`annotation.py:95-107`).** If `self._is_blob`:
1. Set `self.infer_genes`, `self.infer_transcripts`.
2. `self._db_connection = _adapter.build_database_from_string(gff_text=source, infer_genes=..., infer_transcripts=..., merge_strategy=..., id_spec=..., force=True, verbose=...)`.
3. `return` — no preflight validation, no format detection, no directive capture.

**Path branch (`annotation.py:109-142`).** Otherwise, in order:

0. **Existing-DB short-circuit (`annotation.py:110-118`).** If `os.path.exists(self.file_name)` and `os.path.getsize(self.file_name) > 0`: `try: self._db_connection = gffutils.FeatureDB(self.file_name); return`. If `gffutils.FeatureDB(...)` raises (the file isn't a valid SQLite gffutils DB), swallow the exception and fall through. **Gotcha:** this short-circuit always uses `gffutils.FeatureDB` even when `self.backend == "gffbase"`; a path that *is* already a SQLite gffutils DB is opened as gffutils regardless of backend.
1. **Format detection (`annotation.py:121`).** `file_format = self._detect_file_format()` (see §5.2). Branches:
   - `"GTF format"` → `self._handle_gtf_input(...)` (may convert, may raise; see §5.2).
   - `"Unknown format"` → warn if verbose ("Assuming GFF3"), set `infer_genes`/`infer_transcripts` from the constructor args.
   - anything else (GFF3) → set `infer_genes`/`infer_transcripts` from args.
2. **Pre-flight validation (`annotation.py:139`).** `self._run_preflight_validation()` (see §5.1.6).
3. **DB build (`annotation.py:142`).** `self._get_db_connection()` (see §5.1.5).

#### 5.1.5 Database building — `_get_db_connection` and the 3-strategy fallback

`_get_db_connection` (`annotation.py:302-334`) dispatches on backend:

- **gffbase branch (`annotation.py:311-325`):** `db = _adapter.open_existing_db(self.file_name)`; if `None`, `db = _adapter.build_database(file_name=..., infer_genes=..., infer_transcripts=..., merge_strategy=..., id_spec=..., force=True, verbose=...)`. Set `self._db_connection = db`; return. There is **no** 3-strategy fallback for gffbase (its ingest deduplicates internally, per the adapter doc).
- **gffutils branch (`annotation.py:327-334`):** `db_path = self.file_name + "_db"`. `try: self._db_connection = gffutils.FeatureDB(db_path)` (re-use a previously built DB cache); on any exception, `self._db_connection = self._build_database()`.

`_build_database()` (`annotation.py:336-442`) builds a gffutils SQLite DB with up to three progressively more permissive strategies. Precondition: `disable_genes = not self.infer_genes`, `disable_transcripts = not self.infer_transcripts`. Before strategy 1, it calls `self._warn_on_duplicate_ids(self.file_name)` (§5.1.7).

All three strategies call `gffutils.create_db(self.file_name, self.file_name + "_db", ...)` and differ only in `merge_strategy`/`force`/`transform`:

| Strategy | Line | `merge_strategy` | `force` | `transform` | On success |
|---|---|---|---|---|---|
| 1 | 366-378 | `self.merge_strategy` (default `"create_unique"`) | `self.force` | `self._get_transform_func()` | return `db`, print success |
| 2 | 396-413 | `"create_unique"` | `True` | `self._get_unique_id_transform()` | return `db`, warn "renamed _dup1…" |
| 3 | 423-439 | `"merge"` | `True` | `self._get_unique_id_transform()` | return `db`, warn "merged" |

Control flow between strategies:
1. Try strategy 1. On exception `exc1`: call `print_db_build_error(...)`, then check `_is_duplicate_id_error(str(exc1))`. If **not** a duplicate-ID error → `self._fatal_db_error(exc1, "Strategy 1 ...")` (prints a boxed error block and `sys.exit(1)`). Only a recognised duplicate-ID error falls through.
2. Try strategy 2. On exception `exc2`: `print_db_build_error(...)`, fall through (no early exit).
3. Try strategy 3. On exception `exc3`: `print_db_build_error(...)`, then `self._fatal_db_error(exc3, "all 3 strategies")` → `sys.exit(1)`.

`_is_duplicate_id_error(exc_str)` (`annotation.py:780-790`): returns `True` iff the lowercased exception string contains any of `"unique constraint failed"`, `"uniqueconstraintviolation"`, or `"duplicate"`.

`_fatal_db_error(exc, strategy_desc)` (`annotation.py:444-467`): logs a structured error section (file, strategy, error, possible causes, suggested fixes including a `gffread -E` command and a `grep` duplicate-ID hunt), then `sys.exit(1)`.

**Transform functions.**
- `_get_transform_func()` (`annotation.py:474-477`): returns `None` if `not self.infer_genes`; otherwise returns the module-level `_transform_func_gtf`.
- `_transform_func_gtf(x)` (`annotation.py:793-797`): if `"transcript_id" in x.attributes`, append `"_transcript"` to `x.attributes["transcript_id"][0]` (avoids GTF transcript_id/gene_id ID collisions); return `x`. Module alias `transform_func(x)` (`annotation.py:805-807`) just forwards to it.
- `_get_unique_id_transform()` (`annotation.py:479-516`): returns a closure `unique_id_transform(feature)` with a private `seen_ids: dict` counter:
  1. If `self.infer_genes`, first apply `_transform_func_gtf(feature)`.
  2. Determine `original_id`: `feature.id` if truthy, else `feature.attributes["ID"][0]` if `"ID"` present and truthy, else `None`.
  3. If `original_id` is falsy → return `feature` unchanged.
  4. If `original_id` already in `seen_ids`: increment its count, set `new_id = f"{original_id}_dup{seen_ids[original_id]}"`, write it back into `feature.attributes["ID"] = [new_id]` (if `"ID"` present) and `feature.id = new_id` (if `feature` has `.id`).
  5. Else: `seen_ids[original_id] = 0` (first sighting gets count 0; the first *duplicate* becomes `_dup1`).

#### 5.1.6 Pre-flight validation — `_run_preflight_validation`

(`annotation.py:262-294`) Runs `validate_annotation_file(self.file_name, max_duplicate_examples=20, check_orphan_parents=True)` (from `lifton.annotation_validator`, out of scope here). Then:
1. If `result.errors` or `result.warnings` is non-empty → `print_validation_report(result)`.
2. Compute `truly_fatal`: `True` iff any error string (lowercased) contains any of the keywords `"not found"`, `"not readable"`, `"empty"`, `"no valid 9-column"`, `"not a regular file"`.
3. If `truly_fatal` → `logger.log_error(...)` then `sys.exit(1)`.
4. Otherwise (duplicate IDs / orphan parents) → continue; these are handled by the DB-build fallback.

#### 5.1.7 Duplicate-ID pre-scan — `_warn_on_duplicate_ids`

(`annotation.py:740-773`) Called once before strategy 1. Single O(n) pass over the file:
1. Open `path` (`errors="replace"`). On `OSError`, return.
2. For each line: skip if empty or starts with `#`. Split on `\t`; skip if not exactly 9 columns. Take column 9 (`attrs`), split on `;`, find the first pair starting with `"ID="`, extract `ident = pair[3:]`, increment `seen[ident]`, then `break` (one ID per row).
3. `dups = sorted(k for k, v in seen.items() if v > 1)`. If non-empty, log a warning naming the first 5 IDs (with `(+N more)`), stating that gffutils will auto-rename collisions and that downstream lookups against the original IDs may miss the renamed rows.

#### 5.1.8 Public property and query methods

`db_connection` is exposed both as a `@property` returning `self._db_connection` (`annotation.py:187-189`) and, for safety during early `__init__`, via `__getattr__` (`annotation.py:709-714`) which raises `AttributeError(name)` if `_db_connection` is not yet set, and raises a generic `AttributeError` for any other unknown attribute.

The following query methods all delegate to `self._db_connection` (a gffutils `FeatureDB` or a gffbase `FeatureDB`):

| Method | Line | Behaviour |
|---|---|---|
| `get_protein_coding_features(feature_types)` | 599-610 | For each type, for each top-level feature (`is_highest_parent`), collect its `CDS` children; if ≥1 CDS, append `feature.id`. Returns list of IDs. |
| `get_noncoding_features(feature_types)` | 612-623 | Same, but appends `feature.id` only when the CDS child count is `== 0`. |
| `get_novel_protein_coding_features(ref_genes, feature_types)` | 625-631 | `get_protein_coding_features` minus IDs present in `ref_genes`. |
| `get_novel_noncoding_features(ref_genes, feature_types)` | 633-639 | `get_noncoding_features` minus IDs in `ref_genes`. |
| `get_all_parent_feature_ids(feature_types)` | 641-649 | IDs of all top-level features of the given types. |
| `make_parent_to_child_dict(protein_coding, gene)` | 651-664 | Build `defaultdict(list)` mapping each level-1 parent id → its qualifying children. A child qualifies when `protein_coding and child.featuretype == "CDS"`, OR `protein_coding is False and self.is_lowest_child(child)`. If the gene has **no** children at all, set `child_dict[gene] = [self._db_connection[gene]]`. |
| `get_features_of_type(feature_types)` | 666-670 | Flat list of features for all given types. |
| `get_feature_dict(feature_types)` | 672-677 | `{feature.id: feature}` for all features of the given types. |
| `get_source_name(feature_name)` | 679-684 | If the feature has `extra_copy_number` and its value equals the last `_`-split component of `feature.id`, return the id with that last component stripped; else return `feature.id`. |
| `is_lowest_child(feature_name)` | 686-687 | `len(list(children(feature_name))) == 0`. |
| `is_highest_parent(feature_name)` | 689-690 | `len(list(parents(feature_name))) == 0`. |
| `get_paralog_name(feature_name)` | 692-697 | If `extra_copy_number` in attributes, return `"_".join(attributes["ID"][0].split("_")[:-1])`; else `""`. |
| `get_num_levels(feature_name)` | 699-705 | Count hierarchy depth by repeatedly querying `children(..., level=level)` until empty; returns the deepest level reached. |

Module-level standalone helpers in `annotation.py` (kept for back-compat imports):
- `get_perc_id(feature)` (810-813): `float(feature.attributes["sequence_ID"][0])` if present, else `0.0`.
- `merge_children_intervals(children)` (816-828): collect `[start, end]` per child, sort by start, then linearly merge overlapping/abutting intervals (merge when `current[0] <= previous[1]`, updating `previous[1] = max(previous[1], current[1])`). Returns the merged interval list.

---

### 5.2 Format detection and GTF→GFF3 auto-conversion

#### 5.2.1 Format detection

`_detect_file_format()` (`annotation.py:196-203`) calls `extract_sequence.determine_file_format(self.file_name)`. On any exception it logs (if verbose) and returns the safe default `"GFF format"`. Return strings are `"GTF format"`, `"GFF format"`, or `"Unknown format"`.

`_detect_file_format_for(path)` (`annotation.py:251-255`) is the same but for an arbitrary path; on exception returns `"GFF format"`.

#### 5.2.2 GTF handling — `_handle_gtf_input`

(`annotation.py:205-249`) Triggered when `_detect_file_format()` returns `"GTF format"`:
1. Print an informational note to `stderr` ("Detected GTF format … GTF support is experimental …") suggesting manual `gffread -E` or `agat_sp_gtf2gff.pl` conversion.
2. **Force inference:** `self.infer_genes = True` and `self.infer_transcripts = True` (GTF lacks explicit gene/transcript rows). This overrides the constructor args.
3. If `self.auto_convert_gtf` (default `True`):
   - `converted = self._convert_gtf_to_gff3()`.
   - If `converted` is a non-empty existing file **and** `self._detect_file_format_for(converted) == "GFF format"` → log success, set `self.file_name = converted` (subsequent DB build uses the converted file).
   - If converted but not recognised as GFF → warn, keep the original GTF.
   - **If conversion failed entirely** (`converted` is `None`/empty/missing) → raise `LiftOnInputError(...)` with an actionable message (V1.9 fix: fail loudly rather than feeding unparseable GTF into gffutils minutes later).
4. If `self.auto_convert_gtf` is `False` (i.e. `--no-auto-convert-gtf`): restore `self.infer_genes`/`self.infer_transcripts` from the **constructor** arguments (`infer_genes`, `infer_transcripts`) and use the GTF directly.

#### 5.2.3 Conversion — `_convert_gtf_to_gff3`

(`annotation.py:523-579`) Returns the output path or `None`:
1. `base = os.path.splitext(self.file_name)[0]`; `output_file = base + "_converted.gff3"`.
2. `gffread_ok = self._tool_available("gffread")`; `agat_cmd = self._find_agat_command()`.
3. If neither tool available → warn (if verbose) with install hints, return `None`.
4. **gffread first** (preferred): run `["gffread", "-E", self.file_name, "-o", output_file]` with `capture_output=True, text=True, check=True`. If `output_file` exists and is non-empty → return it. On `CalledProcessError`/`FileNotFoundError` → warn and fall through. Empty output → warn.
5. **agat fallback:** build the command depending on the resolved agat command name: if `"gtf2gff" in agat_cmd` use `[agat_cmd, "--gtf", self.file_name, "-o", output_file]`, else `[agat_cmd, "--gff", self.file_name, "-o", output_file]`. Run with `check=True`. Return `output_file` if non-empty; on error warn and continue.
6. If both fail → return `None`.

Helpers: `_tool_available(name)` (`annotation.py:581-586`) runs `["which", name]` with `check=True` and returns whether it succeeds. `_find_agat_command()` (`annotation.py:588-592`) probes, in order, `agat_sp_gtf2gff.pl`, `agat_convert_sp_gff2gtf.pl`, `agat`, returning the first available or `None`.

---

### 5.3 Feature partitioning

Two functions in `lifton_utils.py` decide which reference features become lift-over work and how they are classified coding vs non-coding.

#### 5.3.1 `get_parent_features_to_lift(feature_types_file)`

(`lifton_utils.py:230-246`) Reads the list of top-level feature types LiftOn walks:
1. Default `feature_types = ["gene"]`.
2. If `feature_types_file` (the `-f` / `--features` argument) is not `None`: reset to `[]` and append each `line.rstrip()` read from that file.

Returns `feature_types` (e.g. `["gene"]`, or `["gene", "pseudogene", ...]`).

#### 5.3.2 `get_ref_liffover_features(features, ref_db, intermediate_dir, args)`

(`lifton_utils.py:334-423`) Walks the reference DB and produces the master dictionaries the rest of the pipeline keys against. It also writes two side files: `{intermediate_dir}/ref_feature.txt` (gene-level coding classification) and `{intermediate_dir}/ref_transcript.txt` (transcript-level classification).

**Outputs (returned as a 4-tuple, in order):**

| Return | Type | Keyed by | Value |
|---|---|---|---|
| `ref_features_dict` | `dict` | gene/locus id | a `Lifton_feature` object (see field table below) |
| `ref_features_len_dict` | `dict` | gene/locus id | CDS-span length: `last_CDS.end - first_CDS.start + 1`, or `0` if no CDS |
| `ref_features_reverse_dict` | `dict` | transcript id | parent gene/locus id |
| `ref_trans_exon_num_dict` | `dict` | transcript id | number of CDS children of that transcript (`0` if none) |

A sentinel entry is seeded first (`lifton_utils.py:353-354`): `ref_features_dict["LiftOn-gene"] = Lifton_feature("Lifton-gene")` (a placeholder root for miniprot-only/novel genes). **Gotcha:** the constructor is called with the string `"Lifton-gene"` but the dict key is `"LiftOn-gene"` (capital O) — these are intentionally different spellings; downstream code keys against `"LiftOn-gene"`.

**Algorithm, per top-level feature type `f_itr` in `features`, per `locus` in `ref_db.db_connection.features_of_type(f_itr)`:**
1. `CDS_children = list(children(locus, featuretype='CDS'))` (all descendants, any level).
2. `feature = Lifton_feature(locus.id)`.
3. **Biotype-key precedence** to find `gene_type_key` (`lifton_utils.py:362-382`), driven by `args.annotation_database.upper()`:

   | `annotation_database` (uppercased) | Key precedence (first present wins) |
   |---|---|
   | `"REFSEQ"` | `gene_biotype` → `biotype` |
   | `"GENCODE"`, `"ENSEMBL"`, `"CHESS"` | `gene_type` → `biotype` |
   | anything else | `gene_biotype` → `gene_type` → `biotype` |

   If none of the candidate keys is present in `locus.attributes`, `gene_type_key` stays `None`.
4. **Coding classification (`lifton_utils.py:384-394`):**
   - If `gene_type_key is not None`:
     - If `locus.attributes[gene_type_key][0] == "protein_coding"` **and** `len(CDS_children) > 0` → `feature.is_protein_coding = True`; write `"{locus.id}\tcoding\n"` to `ref_feature.txt`.
     - elif the value is `"lncRNA"` or `"ncRNA"` → `feature.is_non_coding = True`; write `"…\tnon-coding\n"`.
     - else → write `"…\tother\n"`.
   - else (`gene_type_key is None`) → write `"…\tother\n"`.

   **Gotcha:** coding requires *both* the `protein_coding` biotype tag and ≥1 CDS child. A feature tagged `protein_coding` with zero CDS children is classified `other`, not coding.
5. **Children walk (`lifton_utils.py:395-414`):**
   - `exon_children = list(children(locus, featuretype='exon', level=1, order_by='start'))`.
   - If `len(exon_children) > 0` (the locus *directly* owns exons — i.e. a gene row that itself carries exons, no transcript layer): call `__process_ref_liffover_features(locus, ref_db, None)`. With `feature=None`, this helper (`lifton_utils.py:426-428`) is a no-op, so **no** reverse-dict / exon-count entries are written for this branch. This is the single-level (transcript-less) case.
   - Else (the normal gene→transcript→exon hierarchy): for each `transcript` in `children(locus, level=1)`:
     a. `__process_ref_liffover_features(transcript, ref_db, feature)` → adds `transcript.id` into `feature.children` (a set).
     b. Populate `ref_features_reverse_dict[<trans-key>] = locus.id`.
     c. `all_CDS_in_trans = list(children(transcript, featuretype='CDS', order_by='start'))`; set `ref_trans_exon_num_dict[<trans-key>] = len(all_CDS_in_trans)` (≥0).
     d. Transcript-level classification written to `ref_transcript.txt`:
        - if `feature.is_protein_coding and transcript.featuretype == "mRNA"` → `"…\tcoding\n"`;
        - elif `feature.is_non_coding and transcript.featuretype in {"ncRNA","nc_RNA","lncRNA","lnc_RNA"}` → `"…\tnon-coding\n"`;
        - else → `"…\tother\n"`.
6. **Gene-level dict population (`lifton_utils.py:415-420`):**
   - `ref_features_dict[<gene-key>] = feature`.
   - `all_CDS_children = list(children(locus, featuretype='CDS', order_by='start'))`; if non-empty, `ref_features_len_dict[<gene-key>] = all_CDS_children[-1].end - all_CDS_children[0].start + 1`, else `0`.

**The `evaluation_liftoff_chm13` key-rewrite (Gotcha).** Throughout, the dictionary *keys* depend on the boolean `args.evaluation_liftoff_chm13`:
- transcript keys: `transcript.id` normally, else `locus.id[4:]` (strip first 4 chars of the gene id).
- gene keys: `locus.id` normally, else `locus.id[5:]` (strip first 5 chars).

This special evaluation mode rewrites the keys to match a CHM13 naming convention; the default path (`False`) uses the raw IDs. The slice lengths (4 for transcripts, 5 for genes) are asymmetric and load-bearing for that mode.

**`Lifton_feature` fields** (constructed via `lifton_class.Lifton_feature(id)`; relevant fields touched here):

| Field | Type | Set where | Meaning |
|---|---|---|---|
| `id` | `str` | constructor arg | gene/locus id |
| `is_protein_coding` | `bool` | `lifton_utils.py:386` | classified protein-coding |
| `is_non_coding` | `bool` | `lifton_utils.py:389` | classified lncRNA/ncRNA |
| `children` | `set` | `__process_ref_liffover_features` (`:428`) | set of child transcript ids |

---

### 5.4 ID parsing/formation, ID mapping, and overlap/ratio + per-locus alignment helpers

#### 5.4.1 `get_ID_base(id, ref_features_dict=None)` — copy-number suffix stripping

(`lifton_utils.py:163-210`) LiftOn appends a numeric copy-number suffix (`_1`, `_2`, …) to lifted features that have extra copies. This function recovers the original ("base") ID **only when it is safe to do so**.

Algorithm:
1. `splits = id.split("_")`. If `len(splits) < 2` → return `id` unchanged (guards single-component numeric ids like `"0"` from collapsing to `""` — Phase 5 bug #5).
2. `last_part = splits[-1]`. Attempt `copy_num = int(last_part)`. On `ValueError`/`IndexError` (last part not an integer) → return `id` (no copy-number suffix).
3. `id_base = "_".join(splits[:-1])`.
4. If `ref_features_dict is not None`:
   - If `id_base in ref_features_dict.keys()` → return `id_base` (confirmed: the base exists, so the suffix was LiftOn-added).
   - Else → return `id` (the trailing number is part of the original ID, e.g. `FMUND_1`).
5. If `ref_features_dict is None` → return `id` (conservative: cannot verify, so never strip).

**Gotcha:** Without a reference dict you can *never* strip the suffix. The only safe strip is when the stripped base is itself a known reference id. This prevents corrupting natural IDs that legitimately end in `_<number>`.

`get_ID(feature)` (`lifton_utils.py:213-227`) returns `(feature.id, get_ID_base(feature.id))` — note it calls `get_ID_base` with `ref_features_dict=None`, so it returns the id and (effectively) the unmodified id unless the id is single-component.

#### 5.4.2 Liftoff-ID → reference-ID resolution

`__extract_ref_ids(ref_features_dict, liftoff_id)` (`lifton_utils.py:493-502`):
1. If `liftoff_id` is directly a key of `ref_features_dict` → return it.
2. Else `ref_id = get_ID_base(liftoff_id, ref_features_dict)`; if `ref_id in ref_features_dict` → return `ref_id`, else `None`.

`get_ref_ids_liftoff(ref_features_dict, liftoff_gene_id, liftoff_trans_id)` (`lifton_utils.py:458-490`) handles three call shapes:
- `(gene_id, None)` → `(None, __extract_ref_ids(..., gene_id))` returned as `(None, ref_trans_id)` — note the asymmetry: the result of resolving `liftoff_trans_id` is placed in the second slot only when the gene id is None; here gene is provided so this path is `(gene_id, None)` mapped via the gene id. (See exact branches below.)
- Branch 1 (`liftoff_gene_id is None`): `ref_trans_id = __extract_ref_ids(ref_features_dict, liftoff_trans_id)`; return `(None, ref_trans_id)`.
- Branch 2 (`liftoff_trans_id is None`): `ref_gene_id = __extract_ref_ids(ref_features_dict, liftoff_gene_id)`; return `(ref_gene_id, None)`.
- Branch 3 (both present): `ref_gene_id = __extract_ref_ids(ref_features_dict, liftoff_gene_id)`. If `None` → return `(None, None)`. Else `ref_trans_id = get_ID_base(liftoff_trans_id, None)` (note: `None` dict — conservative, only strips single-digit suffix-free single-component change, effectively unchanged) and return `(ref_gene_id, ref_trans_id)`.

**Gotcha:** in branch 3 the transcript id is resolved with `ref_features_dict=None` because `ref_features_dict` contains *gene* ids, not transcript ids — so passing it would never find a match and would never strip. The conservative `None` call is deliberate.

#### 5.4.3 `miniprot_id_mapping(m_feature_db)` — ref_id ↔ miniprot ids

(`lifton_utils.py:431-455`) Builds the bidirectional map between reference transcript ids and miniprot-assigned mRNA ids.
1. If `m_feature_db is None` → return two empty dicts.
2. For each `feature` in `m_feature_db.features_of_type("mRNA")`:
   - `miniprot_id = feature["ID"][0]`.
   - `aa_trans_id = str(feature.attributes["Target"][0]).split(" ")[0]` — the reference protein/transcript id is the first whitespace-delimited token of miniprot's `Target` attribute.
   - Append `miniprot_id` to `ref_id_2_m_id_trans_dict[aa_trans_id]` (a list; created on first sighting).
   - `m_id_2_ref_id_trans_dict[miniprot_id] = aa_trans_id`.
3. Return `(ref_id_2_m_id_trans_dict, m_id_2_ref_id_trans_dict)`.

So one reference id can map to **several** miniprot ids (multiple alignments), but each miniprot id maps to exactly one reference id.

`get_ref_ids_miniprot(ref_features_reverse_dict, miniprot_trans_id, m_id_2_ref_id_trans_dict)` (`lifton_utils.py:505-511`):
1. If `miniprot_trans_id` not in `m_id_2_ref_id_trans_dict` → `(None, None)`.
2. `ref_trans_id = m_id_2_ref_id_trans_dict[miniprot_trans_id]`.
3. If `ref_trans_id` not in `ref_features_reverse_dict` → `(None, ref_trans_id)`.
4. Else → `(ref_features_reverse_dict[ref_trans_id], ref_trans_id)` (gene id, transcript id).

#### 5.4.4 Overlap helpers — `segments_overlap_length` and `check_ovps_ratio`

`segments_overlap_length(segment1, segment2)` (`lifton_utils.py:539-567`):
1. Each segment must be a 2-tuple `(start, end)`; else `raise ValueError("Segments must have exactly 2 endpoints")`.
2. `s1, e1 = segment1`; `s2, e2 = segment2`.
3. `ovp_len = min(e1, e2) - max(s1, s2) + 1` (inclusive-coordinate overlap length).
4. If `ovp_len < 0` → clamp to `0` (V2.9 fix: disjoint segments must not yield negative bp).
5. `ovp = ovp_len > 0`. Return `(ovp_len, ovp)`.

**Gotcha (symmetry):** the overlap length is computed directly from sorted endpoints (`min`/`max`), making the function fully symmetric and order-independent (Phase 5 bug #4). The invariant guaranteed is `ovp_len >= 0` and `ovp_len > 0 ⟺ ovp`.

`check_ovps_ratio(mtrans, mtrans_interval, overlap_ratio, tree_dict)` (`lifton_utils.py:570-602`): tests whether a miniprot transcript overlaps an existing reference feature by more than a fractional threshold.
1. If `mtrans.seqid` not in `tree_dict` → return `False`.
2. Extract `(begin, end)` from `mtrans_interval`: if it has a `.begin` attribute (an `intervaltree.Interval`), use `.begin`/`.end`; else treat as a `(begin, end)` tuple (Phase 5 bug #3 — accept either shape).
3. `ovps = tree_dict[mtrans.seqid].overlap(begin, end)`.
4. For each overlapping interval `ovp`:
   - `ovp_len, _ = segments_overlap_length((mtrans_interval[0], mtrans_interval[1]), (ovp[0], ovp[1]))`.
   - `ref_len = ovp[1] - ovp[0] + 1`; `target_len = mtrans_interval[1] - mtrans_interval[0] + 1`.
   - **Threshold test:** if `(ovp_len / min(ref_len, target_len)) > overlap_ratio` → set `is_overlapped = True` and `break`.
5. Return `is_overlapped`.

The ratio is normalised by the **smaller** of the two interval lengths, so even a small feature fully contained in a large one (or vice versa) counts as a high-ratio overlap. `overlap_ratio` is supplied by the caller (the pipeline default lives in the argument parser, not here).

#### 5.4.5 Per-locus alignment drivers (used by §2.3)

`LiftOn_eval_alignment(eval_trans, locus, tgt_fai, ref_proteins, ref_trans_id, lifton_status)` (`lifton_utils.py:249-253`): runs `align.lifton_parasail_align(...)`; if a non-`None` alignment is returned, set `lifton_status.lifton_aa = eval_aln.identity`; return the alignment.

`LiftOn_liftoff_alignment(lifton_trans, locus, tgt_fai, ref_proteins, ref_trans_id, lifton_status)` (`lifton_utils.py:256-260`): runs `align.lifton_parasail_align(...)`; if non-`None`, set `lifton_status.liftoff = liftoff_aln.identity`; return it. This produces the protein-level identity of the *Liftoff-derived* transcript.

`LiftOn_miniprot_alignment(chromosome, transcript, m_id_dict, m_feature_db, tree_dict, fai, ref_proteins, ref_trans_id, lifton_status)` (`lifton_utils.py:263-331`): finds and aligns the best miniprot transcript for a reference transcript. Returns `(m_lifton_aln, has_valid_miniprot)`.

Initial state: `m_lifton_aln = None`, `has_valid_miniprot = False`.

**Entry guard (`lifton_utils.py:284`):** proceed only if `ref_trans_id in m_id_dict.keys()` **and** `ref_trans_id in ref_proteins.keys()`. Otherwise return `(None, False)`.

For each candidate `m_id` in `m_id_dict[ref_trans_id]` (the list of miniprot ids for this reference id), with `m_entry = m_feature_db[m_id]`:

1. **Check 1 — locus overlap (`:290-294`):** `_, overlap = segments_overlap_length((m_entry.start, m_entry.end), (transcript.start, transcript.end))`. If `not overlap` **or** `m_entry.seqid != transcript.seqid` → `continue` (skip this candidate). The miniprot hit must overlap the current Liftoff transcript's span on the same chromosome.
2. **Check 2 — cross-gene guard (`:303-315`):**
   - `ovps_liftoff = tree_dict[chromosome].overlap(transcript.start, transcript.end)`.
   - `ovps_miniprot = tree_dict[chromosome].overlap(m_entry.start, m_entry.end)`.
   - Build `liftoff_set = {ovp[2] for ovp in ovps_liftoff}` (the gene ids the Liftoff transcript overlaps; `ovp[2]` is the interval's stored data).
   - Set `miniprot_cross_gene_loci = True` and break if any `ovp_miniprot[2]` (gene id) is **not** in `liftoff_set` — i.e. the miniprot hit reaches into a gene the Liftoff transcript does not.
   - If `miniprot_cross_gene_loci` → `continue`. **Gotcha:** this rejects miniprot hits that span more genes than the corresponding Liftoff transcript, preventing miniprot from gluing adjacent genes together.
3. **Valid candidate exists (`:317`):** set `has_valid_miniprot = True`.
4. **Build a `Lifton_TRANS` for the candidate (`:318-326`):**
   - `miniprot_trans = lifton_class.Lifton_TRANS(m_id, "", "", 0, m_entry, {})`.
   - Add exons: iterate `m_feature_db.children(m_entry, featuretype=('CDS','stop_codon'), order_by='start')`, calling `miniprot_trans.add_exon(exon)` for each.
   - Add CDS: iterate the same query again (a fresh generator), calling `miniprot_trans.add_cds(cds)` for each (`cds_num` is counted but unused downstream here). **Gotcha:** both exon and CDS lists are built from the **same** `('CDS','stop_codon')` query — miniprot output carries CDS+stop_codon features, and LiftOn treats them as both the exon structure and the CDS structure.
5. **Align and keep the best (`:327-330`):**
   - `tmp_m_lifton_aln = align.lifton_parasail_align(miniprot_trans, m_entry, fai, ref_proteins, ref_trans_id)`.
   - If `m_lifton_aln is None` **or** `tmp_m_lifton_aln.identity > lifton_status.miniprot`: set `m_lifton_aln = tmp_m_lifton_aln` and `lifton_status.miniprot = m_lifton_aln.identity`. So the best miniprot alignment (highest protein identity) wins, and `lifton_status.miniprot` records that best identity.

Return `(m_lifton_aln, has_valid_miniprot)`. `has_valid_miniprot` is `True` if *any* candidate passed Checks 1 and 2 (even if its alignment was poor); the caller uses it to decide whether miniprot evidence is available for the protein-maximization chaining step.

---

### 5.5 The gffbase adapter shim (`gffbase_adapter.py`)

`gffbase_adapter.py` is a pure translation layer: it maps the handful of gffutils-shaped calls `Annotation` makes onto the gffbase (`lifton.gffbase`, imported as `_gffbase`) API. No LiftOn algorithmic logic lives here. It exposes five functions.

| Function | Line | Purpose |
|---|---|---|
| `db_path_for(file_name)` | 22-27 | Returns the DuckDB cache path: `file_name + ".duckdb"`. **Distinct suffix** from gffutils' `<file>_db` SQLite cache so the two backends never clobber each other on disk. |
| `open_existing_db(file_name)` | 30-43 | If `db_path_for(file_name)` exists, return `_gffbase.FeatureDB(path)`; on any exception or missing file return `None`. Mirrors the legacy `try: gffutils.FeatureDB(<file>_db)` short-circuit. |
| `build_database(*, file_name, infer_genes, infer_transcripts, merge_strategy="create_unique", id_spec=None, force=True, verbose=False, transform=None)` | 46-73 | Build a DuckDB FeatureDB from a file path via `_gffbase.create_db(file_name, dbfn=db_path_for(file_name), force, verbose, merge_strategy, id_spec, transform, disable_infer_genes=not infer_genes, disable_infer_transcripts=not infer_transcripts)`. The 3-strategy retry collapses to **one** call because gffbase deduplicates IDs internally. |
| `build_database_from_string(*, gff_text, dbfn=":memory:", infer_genes=False, infer_transcripts=False, merge_strategy="create_unique", id_spec=None, force=True, verbose=False, transform=None)` | 76-111 | Ingest an in-memory GFF3 blob. If `gff_text` is bytes, `.decode("utf-8")` first. Call `_gffbase.create_db(gff_text, dbfn, from_string=True, force, verbose, merge_strategy, id_spec, transform, disable_infer_genes=not infer_genes, disable_infer_transcripts=not infer_transcripts)`. `from_string=True` makes gffbase materialise a tempfile internally so the caller never touches disk. Default `dbfn=":memory:"` ⇒ purely in-memory DB. This is the `--stream` miniprot fast path. |
| `looks_like_gff3_blob(value)` | 114-122 | Cheap blob/path discriminator (see §5.1.1): `True` only for `bytes`/`bytearray` whose first 4096 bytes contain a tab or start with `b"#"`; `False` for anything else (including `str`). |

**Translation contract.** The shim normalises three gffutils idioms onto gffbase:
1. **Cache discovery:** gffutils opens `<file>_db`; gffbase opens `<file>.duckdb`. `open_existing_db` performs this discovery so `_get_db_connection` (gffbase branch) can short-circuit identically to the gffutils branch.
2. **Inference flags inverted:** LiftOn carries `infer_genes`/`infer_transcripts` as positive booleans; gffbase (like gffutils) takes the negated `disable_infer_*`. The shim does the negation at the boundary (`not infer_genes`, `not infer_transcripts`).
3. **Build collapse:** the gffutils backend's three-strategy duplicate-ID fallback is unnecessary because gffbase's ingest deduplicates internally — hence `build_database` is a single `create_db` call with `force=True`.

**Gotcha (no transform plumbing from `Annotation`):** `Annotation._get_db_connection`'s gffbase branch (`annotation.py:315-323`) does **not** pass a `transform=` argument to `_adapter.build_database`, so the `transform` parameter of `build_database`/`build_database_from_string` defaults to `None` on the LiftOn path. The GTF `_transform_func_gtf`/`_get_unique_id_transform` machinery is therefore gffutils-only; gffbase relies on its internal dedup instead. Any reimplementation must ensure gffbase ingest produces the same effective IDs without those Python-side transforms.
