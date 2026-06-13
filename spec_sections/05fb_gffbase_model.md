### 5.6.B gffbase data model, ingestion & schema

This section specifies the vendored **gffbase** library's data model, on-disk DuckDB schema, the GFF3/GTF parser (native-Rust-or-Python-fallback), the ingestion pipeline that fills the schema, the public `Feature` object, the iterator/writer surface, merge predicates, and the legacy-SQLite export path. gffbase is a first-party DuckDB-backed successor to `gffutils`. A reimplementer should be able to reproduce the byte-level behaviour of every object and SQL pass from this text alone.

All file/line citations are to `lifton/gffbase/`. Coordinates throughout are **1-based, inclusive** (GFF convention).

---

#### 5.6.B.1 The `ParsedFeature` slotted dataclass (`feature.py:22-86`)

`ParsedFeature` is the lightweight record the parser emits, one per data row. It is a `@dataclass(slots=True)`. It carries the raw col-9 bytes for byte-faithful round-trip plus a pre-decoded list of attribute triples.

| Field | Type | Meaning |
|---|---|---|
| `seqid` | `str` | Column 1 (sequence/chromosome name). |
| `source` | `str` | Column 2. |
| `featuretype` | `str` | Column 3. |
| `start` | `Optional[int]` | Column 4, 1-based; `None` if `.`/empty. |
| `end` | `Optional[int]` | Column 5, inclusive; `None` if `.`/empty. |
| `score` | `str` | Column 6, kept as string (`.` or float text). |
| `strand` | `str` | Column 7. |
| `frame` | `str` | Column 8 (phase): one of `.`,`0`,`1`,`2`. |
| `attributes_blob` | `bytes` | Raw column-9 bytes, UTF-8. Preserved verbatim. |
| `attributes_pairs` | `List[Tuple[str,str,int]]` | Long-form `(key, value, multivalue_index)` triples; `idx` preserves multi-value order, e.g. `Parent=a,b,c` → `(Parent,a,0),(Parent,b,1),(Parent,c,2)`. Default empty list. |
| `extra` | `List[str]` | Tab-separated columns past column 9, if any. Default empty list. |

Properties / methods:
- `chrom` (property, `:41`) → alias for `seqid`.
- `stop` (property, `:45`) → alias for `end`.
- `attributes_dict()` (`:49`) → materializes `{key: [values…]}`, preserving first-seen key order and multi-value order. Iterates `attributes_pairs` and `out.setdefault(k, []).append(v)`.
- `from_tuple(tup)` (classmethod, `:58`) → builds from the **11-tuple** the Rust extension yields, in this exact order: `(seqid, source, featuretype, start, end, score, strand, frame, blob, pairs, extra)`. `blob` is coerced to `bytes` if not already; each pair is coerced to `(k, v, int(i))`; `extra` is `list(extra)`.

---

#### 5.6.B.2 The public `Feature` class (`feature.py:218-457`)

`Feature` is the backward-compatible public object mirroring legacy `gffutils.Feature`. It uses `__slots__` (`:226-234`): `seqid, source, featuretype, start, end, score, strand, frame, attributes, extra, bin, id, dialect, file_order, keep_order, sort_attribute_values, _attributes_blob, children`.

##### 5.6.B.2.1 Constructor (`__init__`, `feature.py:236-291`)

Signature defaults: `seqid="."`, `source="."`, `featuretype="."`, `start="."`, `end="."`, `score="."`, `strand="."`, `frame="."`, `attributes=None`, `extra=None`, `bin=None`, `id=None`, `dialect=None`, `file_order=None`, `keep_order=False`, `sort_attribute_values=False`.

Construction algorithm:
1. Assign `seqid, source, featuretype` directly.
2. `start = _coord_to_int(start)`, `end = _coord_to_int(end)`. `_coord_to_int(v)` (`:101`): returns `None` if `v` is `None`/`""`/`"."`; returns `v` if already `int`; else `int(v)`.
3. `score`/`strand`/`frame`: if `None` substitute `"."`, else keep as given.
4. `bin = bin`, `id = id`, `dialect = dialect or {}`, `file_order`, `keep_order`, `sort_attribute_values` assigned. `_attributes_blob = None`, `children = None`.
5. `fmt = dialect.get("fmt","gff3")`.
6. Resolve `attributes` (the field is always a `_LazyAttributes`, see §5.6.B.2.2):
   - if already `_LazyAttributes` → store as-is;
   - elif `bytes`/`bytearray` → set `_attributes_blob = bytes(attributes)` and wrap as `_LazyAttributes(blob=…, dialect_fmt=fmt)` (parsing deferred);
   - elif `None` → `_LazyAttributes(initial={}, dialect_fmt=fmt)`;
   - else → `_LazyAttributes(initial=attributes, dialect_fmt=fmt)`.
7. Resolve `extra`:
   - `None` → `[]`;
   - `bytes`/`bytearray` → decode UTF-8 (errors="replace"), `text.split("\t")` if non-empty else `[]`;
   - `str` → `text.split("\t")` if non-empty else `[]`;
   - else → `list(extra)`.

##### 5.6.B.2.2 `_LazyAttributes` (`feature.py:110-215`)

A `MutableMapping` whose values are **always lists**. `__slots__ = ("_d","_blob","_dialect_fmt","_parsed")`.

- `__init__(initial=None, blob=None, dialect_fmt="gff3")`: sets `_d={}`, `_blob`, `_dialect_fmt`, `_parsed=False`. If `initial` is given, `_ingest(initial)`, then `_parsed=True`, `_blob=None`.
- `_ingest(initial)` (`:137`): another `_LazyAttributes` → materialize it, copy each `k → list(v)`; `Mapping` → store `v if list else [v]`; `list` of triples/pairs → `setdefault(k,[]).append(v)`; otherwise `TypeError`.
- `_materialize()` (`:157`): no-op if `_parsed`. If `_blob` is `None`/`b""` → set `_parsed=True`. Otherwise decode blob UTF-8 (errors="replace"), call `_pyfallback.attributes.parse_attributes(text)`, append each parsed pair to `_d`, set `_parsed=True`, `_blob=None`. **Gotcha:** materialization is the trigger that abandons byte-faithful round-trip (see §5.6.B.2.4).
- All `MutableMapping` dunders (`__getitem__`, `__setitem__`, `__delitem__`, `__iter__`, `__len__`, `__contains__`, `items`, `keys`, `values`) call `_materialize()` first. `__setitem__` wraps non-list values as `[value]`.

##### 5.6.B.2.3 Aliases & dunders

| Member | Line | Behaviour |
|---|---|---|
| `chrom` get/set | `295`/`299` | get/set `seqid`. |
| `stop` get/set | `303`/`307` | get returns `end`; set runs `_coord_to_int`. |
| `__len__` | `313` | `0` if start or end is `None`, else `end - start + 1` (inclusive length). |
| `__repr__` | `318` | `<Feature {type} ({seqid}:{start}-{end}[{strand}]) at {hex(id)}>`. |
| `__str__` / `__unicode__` | `324`/`327` | `_format_line()` (the GFF3 serialization). |
| `__hash__` | `330` | `hash(self._format_line())`. |
| `__eq__` | `333` | `True` iff both are `Feature` and `_format_line()` strings are equal; else `NotImplemented`. |
| `__ne__` | `338` | negation of `__eq__`. |
| `__getitem__` | `344` | `int` key indexes `_FIELD_ORDER`; else indexes `attributes`. |
| `__setitem__` | `349` | `int` key sets field via `_FIELD_ORDER`; else sets attribute. |

`_FIELD_ORDER` (`feature.py:95-98`) = `("seqid","source","featuretype","start","end","score","strand","frame","attributes")`. So `feat[0]` is seqid, `feat[8]` is the attributes mapping.

**Gotcha:** equality and hashing are defined purely on the serialized GFF3 line, **not** on object identity or `id`. Two features with identical column text are equal.

##### 5.6.B.2.4 GFF3/GTF serialization (`_format_attributes` `:357`, `_format_line` `:397`)

`_format_line()`:
1. `start_s = "." if start is None else str(start)`; same for `end_s`.
2. Build columns list: `[seqid, source, featuretype, start_s, end_s, score, strand, frame, _format_attributes()]` then `extend(self.extra)`.
3. Return `"\t".join(cols)`.

`_format_attributes()`:
1. **Byte-faithful fast path:** if `_attributes_blob is not None` AND `attributes` is a `_LazyAttributes` AND `not attributes._parsed` → return `_attributes_blob.decode("utf-8", errors="replace")` verbatim. **Gotcha:** any read/mutation of `attributes` flips `_parsed=True` and abandons this path, so the re-emitted col-9 may differ in key order/encoding from the original.
2. Otherwise reconstruct from the mapping. Resolve from dialect: `fmt = dialect.get("fmt","gff3")`; `sep = "; " if dialect["field separator"]=="; " else ";"`; `kv_sep = dialect.get("keyval separator") or ("=" if fmt=="gff3" else " ")`; `multival = dialect.get("multival separator", ",")`.
3. `items = list(attributes.items())`. If `sort_attribute_values`, replace each value list with `sorted(v)`.
4. For each `(k, vs)`: ensure `vs_list` is a list.
   - GFF3: `joined = multival.join(str(v) for v in vs_list)`; append `f"{k}{kv_sep}{joined}"`.
   - GTF: `quoted = dialect.get("quoted GFF2 values", True)`; for each value append `f'{k}{kv_sep}"{v}"'` if quoted else `f"{k}{kv_sep}{v}"` (one part **per value**, not joined).
5. `s = sep.join(parts)`. If `dialect["trailing semicolon"]` and `s` doesn't end with `;`, append `;`. If `dialect["leading semicolon"]`, prepend `;`.

##### 5.6.B.2.5 Legacy methods

- `astuple(encoding=None)` (`:410`): returns the legacy **12-tuple** for SQLite export: `(id, seqid, source, featuretype, start, end, score, strand, frame, attributes_json, extra_json, bin)`. `attributes_json = json.dumps({k:list(v) …}, separators=(",",":"))`; `extra_json = json.dumps(extra, …)` or `"[]"`. `bin` is `self.bin` if set, else `self.calc_bin()`.
- `calc_bin(_bin=None)` (`:431`): if `_bin` given, set and return it; if start/end `None` return `None`; else `bin = bin_from_coords(start, end)` from `._bins` (UCSC binning, §5.6.B.9), store and return.
- `sequence(fasta, use_strand=True)` (`:441`): if `fasta` is a `str` path, import `pyfaidx` (raise `ImportError` with guidance if missing) and open `pyfaidx.Fasta(fasta)`; else treat `fasta` as a pyfaidx-style mapping. Extract `str(fa[self.seqid][self.start-1 : self.end])` (**0-based slice from 1-based start**). If `use_strand` and `strand=="-"`, apply `_revcomp`. Returns the string.
  - `_revcomp` (`:462`) uses `_COMPLEMENT = str.maketrans("ACGTNacgtn","TGCANtgcan")`: `seq.translate(_COMPLEMENT)[::-1]`.

##### 5.6.B.2.6 `feature_from_row` (`feature.py:477-497`)

Builds a `Feature` from a DuckDB result row. Row column order `_DB_ROW_FIELDS` (`:471`) = `(id, seqid, source, featuretype, start, end, score, strand, frame, attributes_blob, extra_blob, file_order)` — note **12 columns**, this is the shape `FeatureDB._yield_features` returns. Construction: `score`/`strand`/`frame` default to `"."` if `None`; `attributes=blob` (so the lazy blob path is taken, deferring parse); `extra=extra_blob if extra_blob else None`; `id=fid`; `dialect=dialect or {"fmt":"gff3"}`; `file_order` passed through. **Gotcha:** features loaded from the DB but never inspected pay zero attribute-decode cost (lazy blob).

---

#### 5.6.B.3 The DuckDB schema (`schema.py`)

`SCHEMA_VERSION = "1"` (`schema.py:14`). The DDL (`schema.py:16-83`, executed via one `con.execute(DDL)`) creates these tables (all `CREATE TABLE IF NOT EXISTS`):

**`features`** — one row per feature (authored + synthetic):

| Column | Type | Notes |
|---|---|---|
| `id` | `VARCHAR PRIMARY KEY` | derived per §5.6.B.5.2; unique. |
| `seqid` | `VARCHAR NOT NULL` | |
| `source` | `VARCHAR` | |
| `featuretype` | `VARCHAR NOT NULL` | |
| `start` | `BIGINT NOT NULL` | 1-based. `0` substituted when parser start was `None` (ingest `:141`). |
| `"end"` | `BIGINT NOT NULL` | quoted because `end` is reserved. `0` substituted for `None`. |
| `score` | `VARCHAR` | |
| `strand` | `VARCHAR` | |
| `frame` | `VARCHAR` | |
| `attributes_blob` | `BLOB` | raw col-9 bytes. |
| `extra_blob` | `BLOB` | tab-joined extra columns, UTF-8 bytes, or `b""`. |
| `file_order` | `BIGINT` | 1-based encounter order; `NULL` for synthesized rows. |
| `is_synthetic` | `BOOLEAN DEFAULT FALSE` | `TRUE` for GTF-synthesized gene/transcript rows. |
| `seqid_y` | `BIGINT` | R-tree y-band (see §5.6.B.5.4). |

The spatial path additionally adds `bbox GEOMETRY` via `ALTER TABLE features ADD COLUMN IF NOT EXISTS bbox GEOMETRY` (ingest `:287`).

**`seqid_map`** (`:39`): `seqid VARCHAR PRIMARY KEY`, `seqid_y BIGINT NOT NULL`. One row per distinct seqid → its R-tree y-band; consulted at query time to translate region bounds into band space.

**`attributes`** (`:44`): `feature_id VARCHAR NOT NULL`, `key VARCHAR NOT NULL`, `value VARCHAR NOT NULL`, `idx SMALLINT NOT NULL DEFAULT 0`. Long-form attribute storage; `idx` preserves multi-value order.

**`edges`** (`:51`): `parent VARCHAR NOT NULL`, `child VARCHAR NOT NULL`. Direct parent→child relations.

**`closure`** (`:56`): `ancestor VARCHAR NOT NULL`, `descendant VARCHAR NOT NULL`, `depth SMALLINT NOT NULL`. Transitive closure of `edges`; `depth=1` is a direct edge.

**`meta`** (`:62`): `key VARCHAR PRIMARY KEY`, `value VARCHAR`. Key/value config (see §5.6.B.5.6).

**`directives`** (`:69`): `seq BIGINT PRIMARY KEY DEFAULT nextval('directive_seq')`, `directive VARCHAR NOT NULL`. Backed by `CREATE SEQUENCE directive_seq START 1` (`:67`). Preserves `##`-directive lines in encounter order.

**`autoincrements`** (`:74`): `base VARCHAR PRIMARY KEY`, `n BIGINT`. Legacy compat table (typically empty).

**`duplicates`** (`:79`): `original_id VARCHAR NOT NULL`, `new_id VARCHAR PRIMARY KEY`. Records duplicate-id remappings (§5.6.B.5.2).

##### Post-load indexes (`POST_LOAD_INDEXES`, `schema.py:90-99`)

Built **only after** bulk insert + closure materialization. All `CREATE INDEX IF NOT EXISTS`:

| Index | Table(columns) |
|---|---|
| `features_type` | `features(featuretype)` |
| `features_seqstart` | `features(seqid, start, "end")` — the universal region-query B-tree fallback. |
| `attributes_kv` | `attributes(key, value)` |
| `attributes_fid` | `attributes(feature_id)` |
| `edges_parent` | `edges(parent)` |
| `edges_child` | `edges(child)` |
| `closure_ancestor` | `closure(ancestor, depth)` |
| `closure_descend` | `closure(descendant, depth)` |

The R-tree index `features_rtree ON features USING RTREE (bbox)` is created separately by `_finalize_rtree` (§5.6.B.5.4) when the spatial extension loaded. **Gotcha (Phase 19):** the redundant `features_seqid` index was dropped — any `(seqid)`-only predicate is served by the leading prefix of `features_seqstart`.

##### SQLite-compat views (`COMPAT_VIEWS_SQL`, `schema.py:242-254`)

Run after closure is populated. Let `FeatureDB.execute()` accept legacy-schema SQL:
- `features_compat`: projects features with `CAST(attributes_blob AS VARCHAR) AS attributes`, `CAST(extra_blob AS VARCHAR) AS extra`, and a literal `0 AS bin`. **Documented break:** `attributes` here is the raw col-9 string, **not** legacy JSON.
- `relations_compat`: `SELECT ancestor AS parent, descendant AS child, depth AS level FROM closure`.

---

#### 5.6.B.4 The parser (`parser.py` + `_pyfallback/parser.py` + `dialect.py` + `_pyfallback/attributes.py`)

##### 5.6.B.4.1 Engine dispatch (`parser.py`)

At import time (`:16-22`) it tries `from . import _native as _rust`; on `ImportError`, `_rust=None`, `_NATIVE=False`. `native_available()` (`:25`) reports `_NATIVE`.

`_resolve_engine(engine)` (`:74`):
1. `None`/`"auto"` → `"rust"` if `_NATIVE` else `"python"`.
2. `"rust"` but not native → `RuntimeError("Rust extension not built…")`.
3. value not in `{"rust","python"}` → `ValueError`.

Public functions:
- `parse_gff(path, *, checklines=10, force_dialect_check=False, force_gff=False, strict=True, engine="auto")` (`:86`) → resolves engine, calls `_rust.parse_file(...)` or `_pyparser.parse_file(...)`, wraps in `_Iterator(it, native=…)`. Handles `.gz` (pyfallback `_open` uses `gzip` if path ends `.gz`).
- `parse_bytes(data, …)` (`:127`) → same but from a `bytes` buffer; pyfallback decodes UTF-8 (errors="replace") into a `StringIO`.
- `detect_dialect(path, *, checklines=10, engine="auto")` (`:156`) → returns the dialect dict only.

`_Iterator` (`:30`) is the adapter that always yields `ParsedFeature`: `__next__` calls `next(self._inner)`; if native, runs `ParsedFeature.from_tuple(item)`; else passes through. It exposes `dialect()`, `directives()` (as a list), and a `warnings` property (`:62`) returning the inner iterator's structured-error dicts (`line_no`,`kind`,`message`) collected when `strict=False`, else `[]`.

`strict` semantics: `strict=True` (default) raises `GFFFormatError` on the first malformed line; `strict=False` skips malformed lines silently and records them in `.warnings`.

##### 5.6.B.4.2 Line streaming & two-pass dialect peek (`_pyfallback/parser.py`)

The pure-Python fallback is the behavioural oracle; the Rust crate must match it. `_FallbackIterator` (`:307`) wraps the `_stream_features` generator (`:222`).

`_iter_lines(stream)` (`:76`) strips a single trailing `\n` then a single trailing `\r` per line (handles CRLF).

`_stream_features(stream, checklines, force_dialect_check, force_gff, strict, warnings, directives)`:
1. **Peek phase** — iterate lines, `line_no += 1` each:
   - empty line → skip.
   - line starts `##`: if `##FASTA` set `fasta_reached=True` and `break`; else append to `directives`, continue.
   - line starts `#` (single, not `##`) → skip (comment).
   - else parse via `_parse_line_into_feature(line, line_no)`. On `GFFFormatError`: if `_maybe_handle(e)` returns `True` (non-strict) append a warning dict and `continue`; else re-raise.
   - compute the line's dialect observation `obs = parse_attributes(blob)[1]`, append to `samples`; append the feature to `buffered`.
   - if `not force_dialect_check and len(buffered) >= checklines`, `break` (default `checklines=10`).
2. `dialect = merge_dialects(samples)`. If `force_gff`: force `dialect["fmt"]="gff3"`, `dialect["keyval separator"]="="`.
3. Yield each buffered feature as `(feat, directives, dialect)`.
4. If `fasta_reached`, `return`.
5. **Stream phase** — continue iterating remaining lines with the same directive/comment/error handling, yielding `(feat, directives, dialect)`. A `##FASTA` line ends iteration.

`_maybe_handle(err)` (`:240`): if `strict`, returns `False` (caller raises); else appends `{line_no, kind, message}` to `warnings` and returns `True` (caller skips).

`_FallbackIterator.dialect()`/`directives()` lazily `_drain_for_metadata()` (`:335`) — forces the generator to its first yield (or exhaustion) so directives + dialect are populated even on directives-only/empty inputs, stashing the first record back via `_chain` so `__next__` still sees it. `detect_dialect` (`:395`) drives the iterator up to `checklines` records non-strictly, then returns `it.dialect()`.

##### 5.6.B.4.3 Per-line field parsing & validation (`_pyfallback/parser.py:184-219`)

`_parse_line_into_feature(line, line_no)`:
1. `fields = line.split("\t")`. If `len < 9` → `GFFFormatError("expected at least 9 tab-separated fields…", kind="TooFewFields")`.
2. Unpack `seqid, source, featuretype, start_s, end_s, score, strand, frame = fields[:8]`; `blob = fields[8]`; `extra = fields[9:]`.
3. `pairs, _obs = parse_attributes(blob)`.
4. `start = _coord_or_error(start_s, line_no, "start")`, same for `end` (`:171`: `.`/`""`→`None`; non-int → `GFFFormatError kind="InvalidCoordinate"`).
5. `_validate(...)` — return value (an error) raises if non-`None`.
6. Construct `ParsedFeature` with `attributes_blob=blob.encode("utf-8")`.

`_validate` (`:91-168`) returns an error (or `None`) checking, **in order**:

| Check | Condition that errors | `kind` |
|---|---|---|
| empty seqid | `not seqid` | `EmptySeqid` |
| empty featuretype | `not featuretype` | `EmptyFeaturetype` |
| whitespace in type | `any(ch.isspace() for ch in featuretype)` | `InvalidFeaturetype` |
| start ≥ 1 | `start is not None and start < 1` | `InvalidCoordinate` |
| end ≥ start | `start is not None and end is not None and end < start` | `InvalidCoordinate` |
| strand | `strand not in ("+","-","?",".")` | `InvalidStrand` |
| phase domain | `frame not in (".","0","1","2")` | `InvalidPhase` |
| CDS phase required | `featuretype=="CDS" and frame=="."` | `InvalidPhase` |
| score | `score not in ("",".")` and `float(score)` raises | `InvalidScore` |
| attribute structure | for `trimmed = blob.strip()` non-empty and `!= "."`: error if `(not has_eq and not has_quote) or n_pairs==0` where `has_eq = "=" in trimmed`, `has_quote = '"' in trimmed` | `InvalidAttribute` |

**Gotcha:** the CDS-requires-phase rule and the attribute-structure rule are stricter than vanilla GFF3; a CDS row with `.` phase is rejected in strict mode.

##### 5.6.B.4.4 Column-9 attribute parsing (`_pyfallback/attributes.py`)

`parse_attributes(blob)` (`:17`) returns `(pairs, dialect_observation)` where `pairs = [(key, value, multivalue_idx)…]`:
1. Empty blob → `([], default_dialect())`.
2. Observe `leading semicolon` (lstripped startswith `;`), `trailing semicolon` (rstripped endswith `;`), and `field separator` = `"; "` if present, elif `" ; "`, else `";"`.
3. `_split_top_level_semicolons` (`:76`) splits on `;` **outside** double quotes (toggling `in_quotes` on each `"`); a `;` inside quotes sets `obs["semicolon in quotes"]=True`.
4. For each non-empty stripped segment:
   - `_split_keyval` (`:93`): if `"="` present, split on **first** `=` → `(key, val, "=")`; else split on first whitespace → `(key, val, " ")`; else `(seg, "", "=")`.
   - `local_fmt = "gff3" if kv_sep=="=" else "gtf"`. First seen sets `detected_fmt`.
   - `_strip_quotes` (`:104`): strip one pair of surrounding `"`; sets `obs["quoted GFF2 values"]=True` if quoted.
   - GFF3: `multi_values = _split_unquoted_commas(clean_val)` (split on `,` outside quotes); GTF: `[clean_val]` (no comma splitting).
   - track first-appearance `order`; if a key recurs set `obs["repeated keys"]=True`.
   - for each value: GFF3 applies `urllib.parse.unquote` (percent-decode); GTF leaves verbatim. Append `(key, decoded, counter)` with running `counter`.
5. `obs["fmt"] = detected_fmt or "gff3"`; `obs["keyval separator"] = " " if fmt=="gtf" else "="`; `obs["order"] = order`.

**Gotcha:** GFF3 values are percent-decoded on parse, but the byte-faithful round-trip path (§5.6.B.2.4) re-emits the original blob unchanged — decoding only affects the materialized mapping.

##### 5.6.B.4.5 Dialect reconciliation (`dialect.py`)

`default_dialect()` (`:13`):

| Key | Default |
|---|---|
| `fmt` | `"gff3"` |
| `field separator` | `";"` |
| `keyval separator` | `"="` |
| `multival separator` | `","` |
| `leading semicolon` | `False` |
| `trailing semicolon` | `False` |
| `quoted GFF2 values` | `False` |
| `repeated keys` | `False` |
| `semicolon in quotes` | `False` |
| `order` | `[]` |

`merge_dialects(samples)` (`:28`): empty → default. `fmt = "gtf"` iff `n_gtf > (len - n_gtf)` (strict majority of GTF observations; ties → gff3). `keyval = " "` for gtf else `"="`. `field separator` = plurality vote (`max(set, key=count)`). The five boolean keys are OR-ed across all samples. `order` = union of per-sample `order` lists in first-appearance order.

---

#### 5.6.B.5 The ingestion pipeline (`ingest.py`)

`from_file(path, dbfn=":memory:", *, force=False, batch_size=50_000, max_depth=8, disable_infer_genes=False, disable_infer_transcripts=False, gtf_subfeature="exon", engine="auto", build_rtree=True)` (`ingest.py:246`) returns `(con, IngestStats)`. Defaults: `DEFAULT_BATCH_SIZE=50_000` (`:39`), `DEFAULT_MAX_DEPTH=8` (`:40`), `SEQID_Y_BAND=1_000_000` (`:496`).

##### 5.6.B.5.1 Setup
1. If `dbfn != ":memory:"`: if it exists and not `force` → `ValueError`; if exists and `force` → `os.unlink`.
2. `con = duckdb.connect(dbfn)`; `_apply_pragmas(con)`; `con.execute(DDL)`.
3. `_apply_pragmas` (`:441`): if env `GFFUTILS2_THREADS` set → `PRAGMA threads = N`; then `PRAGMA disable_progress_bar` (swallow `duckdb.Error`).
4. **Spatial up-front (Phase 19):** `has_spatial = build_rtree and not _rtree_disabled_by_env() and _try_load_spatial(con)`. `_rtree_disabled_by_env` (`:499`) returns True if env `GFFBASE_TEST_DISABLE_RTREE` ∈ `{"1","true","yes"}`. `_try_load_spatial` (`:508`) runs `INSTALL spatial` + `LOAD spatial`, returns True iff both succeed (else `duckdb.Error` → False). If `has_spatial`, `ALTER TABLE features ADD COLUMN IF NOT EXISTS bbox GEOMETRY`.

##### 5.6.B.5.2 Bulk parse → Arrow → DuckDB
1. `it = _parser.parse_gff(path, engine=engine)`. `seqid_to_y = {}`. `builder = _ArrowBatchBuilder(seqid_to_y, has_spatial)`. `autoinc={}`, `id_counts={}`, `duplicate_pairs=[]`, `file_order=0`, `n_raw=0`, `_fmt_cache=None`.
2. For each parsed `feat`:
   - `file_order += 1`, `n_raw += 1`.
   - On the first feature, `_fmt_cache = _dialect_fmt_safe(it)` (`:428`: `it.dialect().get("fmt","gff3")`, default `"gff3"` on exception). Resolved once to avoid per-record Rust↔Python cost.
   - `fid = _derive_id(feat, _fmt_cache, autoinc)` (see ID inference below).
   - **Duplicate-id handling (mimics gffutils `create_unique`):** `seen = id_counts.get(fid,0)`. If `seen` truthy: `new_fid = f"{fid}__{seen+1}"`, append `(fid,new_fid)` to `duplicate_pairs`, `id_counts[fid]=seen+1`, `fid=new_fid`. Else `id_counts[fid]=1`. **Gotcha:** the **first** occurrence keeps the bare id; duplicates get `__2`, `__3`, … suffixes (double underscore). This is essential for NCBI RefSeq where one `ID=cds-…` spans multiple exon-segment rows.
   - `builder.append(fid, feat, file_order)`; flush when `len(builder) >= batch_size`.
3. Final `builder.flush_into(con)`.
4. If `duplicate_pairs`: register an Arrow table `{original_id, new_id}` and `INSERT INTO duplicates`.

**ID inference — `_derive_id(feat, dialect_fmt, autoincrement)` (`:228`):**
- GFF3 (`dialect_fmt=="gff3"`): scan `attributes_pairs` for `key=="ID"`, return its value.
- Otherwise (GTF, or GFF3 with no `ID`): `n = autoincrement.get(featuretype,0)+1`; `autoincrement[featuretype]=n`; return `f"{featuretype}_{n}"`. **Gotcha:** GTF leaf rows (exon/CDS) get synthetic ids like `exon_1`; gene/transcript ids are filled by the synthesis passes.

##### 5.6.B.5.3 `_ArrowBatchBuilder` (`ingest.py:64-220`)

Column-oriented accumulator. Two explicit Arrow schemas (no inference):
- `FEATURES_SCHEMA` (`:75`): `id`(string), `seqid`(string), `source`(string), `featuretype`(string), `start`(int64), `end`(int64), `score`(string), `strand`(string), `frame`(string), `attributes_blob`(binary), `extra_blob`(binary), `file_order`(int64), `is_synthetic`(bool), `seqid_y`(int64).
- `ATTRIBUTES_SCHEMA` (`:92`): `feature_id`(string), `key`(string), `value`(string), `idx`(int16).

`append(feat_id, feat, file_order)` (`:128`):
1. **Lazy y-band:** `y = seqid_to_y.get(seqid)`; if `None`, `y = len(seqid_to_y) * SEQID_Y_BAND` and store. So the first distinct seqid gets `seqid_y=0`, second `1_000_000`, etc. (encounter order). The dict is shared across batches so bands are stable for the whole file.
2. Append all feature columns. `start`/`end` substitute `0` when the parser value is `None`. `extra_blob = ("\t".join(extra)).encode("utf-8")` if extra else `b""`. `is_synthetic = False`.
3. For each `(k,v,idx)` in `feat.attributes_pairs`, append to the attribute columns.

`flush_into(con)` (`:186`): no-op if empty. Build the two Arrow tables, `con.register("__staging_features", …)` and `__staging_attributes`. Then:
- If `has_spatial`: `INSERT INTO features (… , bbox) SELECT …, ST_MakeEnvelope(start, seqid_y, "end", seqid_y+1) FROM __staging_features` — the bbox envelope is `[start..end] × [seqid_y .. seqid_y+1]` so each seqid occupies a 1-unit-tall band.
- Else: `INSERT INTO features (… columns without bbox …) SELECT * FROM __staging_features`.
- `INSERT INTO attributes SELECT * FROM __staging_attributes`. Unregister both staging tables and `_reset()`.

##### 5.6.B.5.4 Normalization passes (set-based SQL)

After bulk load, `dialect = it.dialect()`, `directives = list(it.directives())`, `fmt = dialect.get("fmt","gff3")`. Directives are inserted in one shot via an Arrow `{directive}` table → `INSERT INTO directives`.

**Branch on `fmt`:**

If `fmt == "gtf"`:
1. If not `disable_infer_transcripts`: `_synthesize_transcripts(con, gtf_subfeature)`.
2. If not `disable_infer_genes`: `_synthesize_genes(con, gtf_subfeature)`.
3. `con.execute(EDGES_FROM_GTF)`.
4. Patch `seqid_y`/`bbox` for synthesized rows (they were inserted with `NULL` seqid_y) via a single targeted `UPDATE … FROM seqid_map m WHERE features.seqid = m.seqid AND features.seqid_y IS NULL` (with or without the `ST_MakeEnvelope` bbox term depending on `has_spatial`).

Else (GFF3): `con.execute(EDGES_FROM_PARENT)`.

**`EDGES_FROM_PARENT`** (`schema.py:109`): `INSERT INTO edges SELECT a.value AS parent, a.feature_id AS child FROM attributes a WHERE a.key='Parent'`. **Gotcha:** a feature with `Parent=p1,p2` produces **two** edge rows (one per multi-value).

**`EDGES_FROM_GTF`** (`schema.py:120`): a `UNION ALL` of two inserts:
- transcript→child: rows where `key='transcript_id'`, `f.featuretype <> 'transcript'`, `a.value <> f.id` (no self-loop), and the parent id exists.
- gene→transcript: rows where `key='gene_id'`, `f.featuretype = 'transcript'`, `a.value <> f.id`, and the parent exists.

**`_synthesize_transcripts(con, subfeature)`** (`:455`): counts existing `transcript` rows, runs `GTF_SYNTHESIZE_TRANSCRIPTS` with param `subfeature` (default `"exon"`), counts after, `n = after-before`, then `GTF_SYNTHESIZE_TRANSCRIPT_ATTRS` and `GTF_PROPAGATE_GENE_ID`. Returns `n`.
- `GTF_SYNTHESIZE_TRANSCRIPTS` (`schema.py:143`): `INSERT INTO features … SELECT a.value AS id, ANY_VALUE(f.seqid), source='gffbase_derived', featuretype='transcript', MIN(f.start), MAX(f."end"), score='.', ANY_VALUE(f.strand), frame='.', NULL attrs/extra/file_order, is_synthetic=TRUE FROM features f JOIN attributes a ON a.feature_id=f.id AND a.key='transcript_id' WHERE f.featuretype = ? AND a.value NOT IN (SELECT id FROM features) GROUP BY a.value`. So one synthetic transcript per distinct `transcript_id` over the subfeature rows, spanning min-start..max-end.
- `GTF_SYNTHESIZE_TRANSCRIPT_ATTRS` (`:170`): inserts a self-referential `transcript_id=<id>` attribute (idx 0) for each synthesized transcript.
- `GTF_PROPAGATE_GENE_ID` (`:183`): for each synthesized transcript, picks the most common `gene_id` among its subfeatures' `(transcript_id, gene_id)` pairs (CTE `pairs` counts, `ranked` uses `ROW_NUMBER() OVER (PARTITION BY transcript_id ORDER BY cnt DESC, gene_id)`, takes `rn=1`) and inserts that `gene_id` attribute (idx 0). **Gotcha:** ties broken by lexicographic `gene_id`. This deliberately avoids the `edges` table so edge population can be a single later pass (no duplicate edges).

**`_synthesize_genes(con, subfeature)`** (`:474`): counts existing `gene` rows, runs `GTF_SYNTHESIZE_GENES` with param `subfeature`, computes `n`, then mirrors `gene_id` into attributes: `INSERT INTO attributes SELECT f.id, 'gene_id', f.id, 0 FROM features f WHERE featuretype='gene' AND is_synthetic=TRUE`. Returns `n`.
- `GTF_SYNTHESIZE_GENES` (`schema.py:211`): same shape as transcript synthesis but `featuretype='gene'`, grouping rows where `f.featuretype IN ('transcript', ?)` (transcript + subfeature) carrying a `gene_id`, `a.value NOT IN (SELECT id FROM features)`, `GROUP BY a.value`.

##### 5.6.B.5.5 Closure, indexes, R-tree, compat views
1. **Closure:** `con.execute(CLOSURE_RECURSIVE_CTE, [max_depth])`. `CLOSURE_RECURSIVE_CTE` (`schema.py:261`): `INSERT INTO closure WITH RECURSIVE walk(ancestor, descendant, depth) AS (SELECT parent, child, 1 FROM edges UNION ALL SELECT w.ancestor, e.child, w.depth+1 FROM walk w JOIN edges e ON e.parent=w.descendant WHERE w.depth < ?) SELECT * FROM walk`. So depth-1 rows are direct edges; recursion bounded by `max_depth` (default 8). **Gotcha:** depth is capped — relations deeper than `max_depth` are not materialized.
2. **Indexes:** `con.execute(POST_LOAD_INDEXES)` (§5.6.B.3).
3. **R-tree finalize:** if `has_spatial`, `rtree_built = _finalize_rtree(con, seqid_to_y)` (`:519`): `DELETE FROM seqid_map`, `executemany INSERT INTO seqid_map` from `list(seqid_to_y.items())` (stable encounter order, first seqid → y=0), then `CREATE INDEX IF NOT EXISTS features_rtree ON features USING RTREE (bbox)`. Returns True, or False on `duckdb.Error`. Since bbox was populated inline during INSERT, this involves **no UPDATE passes**.
4. **Compat views:** `con.execute(COMPAT_VIEWS_SQL)` (must run after closure exists).

##### 5.6.B.5.6 Stats & meta
- `n_attributes`, `n_edges`, `n_closure` from `COUNT(*)` queries.
- `_write_meta(con, dialect, fmt, *, rtree_built, max_depth)` (`:546`): computes `closure_max_depth = MAX(depth) FROM closure` (0 if empty); `INSERT OR REPLACE INTO meta` the rows: `schema_version=SCHEMA_VERSION`, `dialect=json.dumps(dialect or {})`, `fmt`, `rtree_built="true"/"false"`, `max_depth`, `closure_max_depth`. The relational dispatcher reads `closure_max_depth` to choose cached-closure vs. dynamic CTE without a per-call query.
- Returns `IngestStats(...)`.

**`IngestStats`** (`ingest.py:43-55`) fields: `n_features_raw`, `n_features_synthetic_transcripts`, `n_features_synthetic_genes`, `n_attributes`, `n_edges`, `n_closure_rows`, `rtree_built` (bool), `fmt` (str), `dialect` (dict), `directives` (list).

---

#### 5.6.B.6 `create_db` legacy facade (`create_db.py`)

`create_db(data, dbfn, *, …)` (`create_db.py:22`) is the drop-in successor to `gffutils.create_db`. It accepts the full legacy kwarg list; **actively honored** kwargs (Phase 5): `data`, `dbfn`, `force`, `checklines`, `from_string`, `disable_infer_genes`, `disable_infer_transcripts`, `gtf_subfeature`. All other kwargs (`id_spec`, `merge_strategy="error"`, `transform`, `gtf_transcript_key`, `gtf_gene_key`, `infer_gene_extent`, `keep_order`, `text_factory=str`, `pragmas`, `sort_attribute_values`, `dialect`, `force_gff`, `force_dialect_check`, `verbose`, `force_merge_fields`) are accepted-but-no-op for signature compatibility.

Algorithm:
1. If `from_string`: write `data` to a `NamedTemporaryFile(suffix=".gff3", delete=False)` and set `path` to its name, `cleanup_path = path`. Else `path = data`.
2. `con, stats = _ingest.from_file(path, dbfn=dbfn, force=force, disable_infer_genes=…, disable_infer_transcripts=…, gtf_subfeature=…)`.
3. `finally`: if a temp file was created and not `_keep_tempfiles`, `os.unlink` (swallow `OSError`).
4. Return `FeatureDB((con, stats), keep_order=…, sort_attribute_values=…, text_factory=…, pragmas=…)`.

---

#### 5.6.B.7 `DataIterator` (`iterators.py`)

`DataIterator(data, checklines=10, transform=None, force_dialect_check=False, from_string=False, **kwargs)` (`iterators.py:79`) is a factory returning `_DataIterator` (`:18`). Unlike the raw parser, it yields full **`Feature`** objects.

- Construction: if `from_string`, `parser.parse_bytes(data-as-bytes, …)`; else `parser.parse_gff(data, …)`. Stores `_transform`.
- `__next__` (`:47`): `pf = next(self._inner)`; build `Feature(seqid=pf.seqid, …, attributes=pf.attributes_blob, extra="\t".join(pf.extra) or None, dialect=self._inner.dialect() or {"fmt":"gff3"})`. If `_transform` set: `out = transform(feat)`; if `out is False` skip this record (recurse to next); if `out is not None` replace `feat`; return `feat`. **Gotcha:** `attributes` is passed as the raw blob → lazy, byte-faithful.
- Properties `dialect` (`:70`) and `directives` (`:74`) proxy the inner parser.

---

#### 5.6.B.8 `GFFWriter` (`gffwriter.py`)

`GFFWriter(out, with_header=True, in_place=False)` (`gffwriter.py:16`) writes `Feature` records back to GFF/GTF.

Constructor (`:19`):
- If `out` has a `.write` attribute → use it directly as `_fh`.
- Elif `in_place` → atomic write: open a `NamedTemporaryFile` in the **target's directory** with suffix `.gffbase.tmp`, then reopen as a UTF-8 text file; `_target_path` is the final path, `_opened_path` the temp.
- Else → open `out` directly for write.
- If `with_header`, immediately write `##gff-version 3\n`.

Methods:
- `write_rec(rec)` (`:50`): if `rec` is `str`, `rec.rstrip("\n")`; else `str(rec)` (uses `Feature.__str__`). Write line + `\n`.
- `write_recs(recs)` (`:57`): loop `write_rec`.
- `write_gene_recs(db, gene_id)` (`:61`): resolve gene (`db[gene_id]` if str), write it, then write `db.children(gene, level=None, order_by="start")` (all descendants).
- `write_mRNA_children(db, mrna_id)` (`:67`): write mrna then `db.children(mrna, level=1, order_by="start")` (direct children only).
- `write_exon_children(db, exon_id)` (`:73`): write exon then its `level=1` children.
- `close()` (`:79`): flush + close; if `in_place`, `shutil.move(_opened_path, _target_path)` (atomic swap).
- Context-manager protocol (`__enter__`/`__exit__`) closes on exit.

---

#### 5.6.B.9 UCSC binning (`_bins.py`)

`bin_from_coords(start, end)` (`_bins.py:16`) — used only by the SQLite export's `bin` column; the runtime query layer uses the R-tree / seqstart B-tree and never reads `bin`. Constants: `_BINOFFSETS=(512+64+8+1, 64+8+1, 8+1, 1, 0)`, `_BINFIRSTSHIFT=17`, `_BINNEXTSHIFT=3`. Algorithm: convert to 0-based half-open (`start0=start-1`, `end0=end`); `start_bin = start0 >> 17`, `end_bin = (end0-1) >> 17`; for each offset, if `start_bin==end_bin` return `offset+start_bin`, else shift both `>> 3`; fallback return `0`.

---

#### 5.6.B.10 Merge predicates (`merge_criteria.py`) and `FeatureDB.merge`

`merge_criteria.py` is a set of pure predicates with signature `(acc: Feature, cur: Feature, components: list[Feature]) -> bool`. `acc` is the running accumulator; `cur` the candidate; `components` the already-folded features. Each returns True if `cur` should be merged into `acc`.

| Predicate | Line | Returns True when |
|---|---|---|
| `seqid(acc,cur,comp)` | `17` | `acc.seqid == cur.seqid` |
| `strand(acc,cur,comp)` | `21` | `acc.strand == cur.strand` |
| `feature_type(acc,cur,comp)` | `25` | `acc.featuretype == cur.featuretype` |
| `exact_coordinates_only` | `29` | `acc.start==cur.start and acc.end==cur.end` |
| `overlap_end_inclusive` | `33` | `acc.start <= cur.start <= acc.end + 1` (overlap or immediately abutting) |
| `overlap_start_inclusive` | `38` | `acc.start <= cur.end + 1 <= acc.end + 1` |
| `overlap_any_inclusive` | `42` | `overlap_end_inclusive(...) or overlap_start_inclusive(...)` |
| `overlap_end_threshold(t)` | `48` | factory → predicate `abs(acc.end - cur.start) <= t` |
| `overlap_start_threshold(t)` | `54` | factory → predicate `abs(acc.start - cur.end) <= t` |
| `overlap_any_threshold(t)` | `60` | factory → `end_thr OR start_thr` |

**Gotcha:** the `+1` in the inclusive overlap predicates means features that are exactly adjacent (gap of 0 bases between them) merge — this is the "immediately after" allowance.

`FeatureDB.merge(features, merge_criteria=None, multiline=False)` (`interface.py:1216`):
1. Default criteria (`:1219`) = `(seqid, overlap_end_inclusive, strand, feature_type)`.
2. `feats = sorted(features, key=lambda f: (f.seqid, f.start, f.end))`. Empty → return.
3. Iterate; `accum=None`, `components=[]`. For each `f`: if `accum is None`, `accum = _clone_for_merge(f)`, `components=[f]`. Else if **all** predicates pass, `accum.end = max(accum.end, f.end)` and append to `components`. Else: set `accum.children = list(components)`, `yield accum`, then start a new accumulator from `f`.
4. After the loop, yield the final accumulator (with its `children`).
- `_clone_for_merge(f)` (`:1242`) makes a fresh `Feature` copying seqid/source/type/start/end/score/strand/frame, a deep-copied attributes mapping, and the dialect.

**Gotcha:** merging extends only the **end** (`max(accum.end, f.end)`); since input is sorted by `(seqid,start,end)`, `accum.start` is the minimum start of the run and never changes.

---

#### 5.6.B.11 Legacy SQLite export (`sqlite_export.py`)

`export_sqlite(con, path, force=False)` (`sqlite_export.py:72`) writes a legacy-`gffutils`-format SQLite `.db` from a gffbase DuckDB connection. Returns `os.path.abspath(path)`.

The target SQLite schema `_LEGACY_SCHEMA` (`:25-69`) creates tables `features(id,seqid,source,featuretype,start,end,score,strand,frame,attributes,extra,bin, PK id)`, `relations(parent,child,level, PK(parent,child,level))`, `meta(dialect,version)`, `directives(directive)`, `autoincrements(base,n, PK base)`, `duplicates(idspecid,newid, PK newid)`, plus indexes `featuretype`, `seqidstartend` on `(seqid,start,end)`, `relationsparent`, `relationschild`, `binindex` on `(bin)`.

Algorithm:
1. If `path` exists: not `force` → `ValueError`; else `os.unlink`.
2. `sqlite_con = sqlite3.connect(path)`; `executescript(_LEGACY_SCHEMA)`.
3. **Features:** query `SELECT id, seqid, source, featuretype, start, "end", score, strand, frame, CAST(attributes_blob AS VARCHAR) AS attributes, CAST(extra_blob AS VARCHAR) AS extra FROM features ORDER BY file_order NULLS LAST, id`. **Gotcha:** synthesized rows (NULL file_order) sort last, then by id. For each row compute `ucsc_bin = bin_from_coords(start,end) if start and end else None`, substitute `""` for null attributes/extra, and `executemany INSERT INTO features` (12 columns including bin).
4. **Relations:** `SELECT ancestor, descendant, depth FROM closure` → `INSERT INTO relations` (`level=depth`).
5. **Meta:** read all `meta` key/values; insert `(meta.get("dialect", '{"fmt":"gff3"}'), "gffbase-export")`.
6. **Directives:** `SELECT directive FROM directives ORDER BY seq` → insert.
7. **Autoincrements:** best-effort `SELECT base,n FROM autoincrements`; insert if any (swallow `duckdb.Error`).
8. `commit()` in `try`, `sqlite_con.close()` in `finally`.

**Documented break:** the exported `attributes` column is the raw col-9 UTF-8 string, **not** legacy JSON. Code reading attributes via `Feature.attributes` continues to work; raw-SQL consumers expecting JSON-decoded attributes must migrate.
