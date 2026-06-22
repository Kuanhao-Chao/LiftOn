### 5.6 Vendored gffbase internals — Part A: public API & query interface

`gffbase` is LiftOn's first-party, DuckDB-backed successor to `gffutils`. It is imported under the `lifton.gffbase` namespace and is a **drop-in replacement** for the subset of the `gffutils` API that LiftOn consumes (`FeatureDB`, `Feature`, `DataIterator`, `GFFWriter`, `create_db`, plus the exception classes). This Part A documents the **read/query API surface a reimplementer must replicate**: the package's public exports, the `FeatureDB` class and every method LiftOn (and its proxies) call against it, the DuckDB connection-management contract, the exact SQL each query variant emits, the directive/header storage model, the UCSC bin computation, and the exception hierarchy. Ingestion (`create_db`/`ingest`) and the on-disk schema are documented in Part B; this section assumes the schema tables exist and describes only how they are queried.

The single most important reimplementation contract: **the row ordering, return types, and exhaustion semantics of every generator must match exactly**, because the proxy DBs in `lifton/locus_pipeline.py` rebuild this API from in-memory Python dicts and must be byte-for-byte interchangeable with the real `FeatureDB` (the 24-cell byte-identity matrix pins this). Where ordering is load-bearing, it is flagged "Gotcha:".

#### 5.6.1 Package surface (`lifton/gffbase/__init__.py`)

The package re-exports a flat legacy-compatible surface. `__all__` (`__init__.py:38-61`) is:

| Export | Source module | Role |
|---|---|---|
| `create_db` | `.create_db` | Build a `FeatureDB` from a GFF3/GTF file or string (Part B). |
| `FeatureDB` | `.interface` | The query class documented here. |
| `Feature` | `.feature` | Public 1-based-inclusive feature object (§5.6.10). |
| `DataIterator` | `.iterators` | Streaming line→`Feature` iterator (Part B). |
| `GFFWriter` | `.gffwriter` | Serialise `Feature`s back to GFF3 (Part B). |
| `example_filename` | `.helpers` | Locate a packaged test fixture (§5.6.11). |
| `FeatureNotFoundError`, `DuplicateIDError`, `AttributeStringError`, `EmptyInputError`, `GFFFormatError` | `.exceptions` (or `._native`) | Exception hierarchy (§5.6.12). |
| `merge_criteria` | `.merge_criteria` | Predicate functions for `FeatureDB.merge`. |
| `ParsedFeature` | `.feature` | Slotted parser-output dataclass. |
| `parse_gff`, `parse_bytes`, `detect_dialect`, `native_available` | `.parser` | Low-level parse entry points (Part B). |
| `ingest`, `export_sqlite` | `.ingest`, `.sqlite_export` | Ingestion engine + legacy SQLite dump. |
| `__version__` | literal | `"0.1.0"` (`__init__.py:63`). |

**`GFFFormatError` rebinding (`__init__.py:23-26`).** The package tries `from ._native import GFFFormatError` first; if the Rust extension is loaded, that PyO3 `create_exception!` type is the canonical class (it is what Rust actually raises). On `ImportError`, it falls back to the pure-Python `GFFFormatError` from `.exceptions`. Both inherit `ValueError`, so `except ValueError` and `except GFFFormatError` keep working regardless of which path raised. A reimplementer must ensure whichever class is publicly bound is the one raised by the parser, or `isinstance`/`except` checks break.

#### 5.6.2 `FeatureDB` — construction & connection management

```python
def __init__(self, dbfn, default_encoding="utf-8", keep_order=False,
             pragmas=None, sort_attribute_values=False, text_factory=str)
```
(`interface.py:46-138`). The `dbfn` argument is polymorphic; resolution (`interface.py:62-77`):

| `dbfn` type | Branch | Effect |
|---|---|---|
| `duckdb.DuckDBPyConnection` | `interface.py:62-64` | Adopt the connection directly; `self.dbfn = ":existing-connection:"`. |
| `(con, IngestStats)` 2-tuple where `con` is a connection | `interface.py:65-70` | Internal path used by `create_db`; adopt `dbfn[0]`; `self.dbfn = ":existing-connection:"`. |
| `str` (path) | `interface.py:71-73` | `self.conn = duckdb.connect(dbfn, read_only=False)`; `self.dbfn` is the path string. |
| anything else | `interface.py:74-77` | `raise TypeError`. |

**Gotcha:** `self.dbfn` is the load-bearing field that `_ThreadLocalCtxFactory._extract_dbfn` (`locus_pipeline.py:625-639`) reads to decide whether a DB can be reopened per-thread. A path string ⇒ reopen-able; `":existing-connection:"` / `":memory:"` ⇒ not reopen-able (falls back to the parent connection).

After resolving the connection, construction does the following, **in order**:

1. **Apply pragmas** (`interface.py:79-80`): if `pragmas` was passed, call `self.set_pragmas(pragmas)` (§5.6.9).
2. **Read provenance from `meta`** (`interface.py:83-87`): `meta = self._read_meta()` → `dict`. Then:
   - `self.dialect = self._parse_dialect(meta.get("dialect"))` — JSON-decode the stored dialect string; on missing/invalid value return `{"fmt": "gff3"}` (`interface.py:140-147`).
   - `self.fmt = meta.get("fmt", "gff3")`.
   - `self._rtree_built = meta.get("rtree_built", "false").lower() == "true"`.
   - `self._max_depth = int(meta.get("max_depth", "8"))` — the materialised-closure depth bound; default **8**.
3. **Confirm R-tree existence** (`interface.py:88-90`): if `_rtree_built`, re-check via `_has_rtree_index()` (queries `duckdb_indexes()` for an index named `features_rtree` on table `features`; `interface.py:156-164`). A stale meta flag is overridden by the actual index presence.
4. **Load the spatial extension** (`interface.py:96-104`): if `_rtree_built`, attempt `LOAD spatial`. The R-tree's functions (`ST_Intersects`, `ST_MakeEnvelope`) are **not** auto-loaded just because the index exists. On `duckdb.Error`, retry `INSTALL spatial; LOAD spatial`; if that also fails (offline HPC node), set `_rtree_built = False` and fall back to the B-tree path. **Gotcha:** the R-tree-vs-B-tree choice is purely a planner detail — both paths return identical rows (`interface.py:10-13` docstring), so output is unaffected by the fallback.
5. **Read closure max-depth** (`interface.py:109-117`): `self._closure_max_depth = int(meta["closure_max_depth"])` if present, else `SELECT MAX(depth) FROM closure` (0 when the table is empty or absent). This is the corpus's *actual* hierarchy depth, distinct from `_max_depth` (the cache's *capacity*).
6. **Load the seqid→y-band map** (`interface.py:122-131`): `self._seqid_y_map = {}`; if `_rtree_built`, populate from `SELECT seqid, seqid_y FROM seqid_map` as `{seqid: int(seqid_y)}`. On error (older DB), set `_rtree_built = False`.
7. **Load directives** (`interface.py:133-138`): `self.directives = [row[0] for row in conn.execute("SELECT directive FROM directives ORDER BY seq")]` — the ordered list of `##` header lines preserved from the input (see §5.6.8).

Instance fields set by the constructor:

| Field | Type | Meaning |
|---|---|---|
| `conn` | `duckdb.DuckDBPyConnection` | The live connection; all queries route through it. |
| `dbfn` | `str` | Path, or `":existing-connection:"`. |
| `default_encoding` | `str` | `"utf-8"`; carried, not used in query path. |
| `keep_order`, `sort_attribute_values`, `text_factory` | mixed | gffutils-compat passthroughs. |
| `_analyzed_flag` | `bool` | Set by `analyze()` (§5.6.9). |
| `dialect` | `dict` | Used by `feature_from_row` to format col-9 on `Feature.__str__`. |
| `fmt` | `str` | `"gff3"` / `"gtf"`. |
| `_rtree_built` | `bool` | Whether `region()` may use the spatial path. |
| `_max_depth` | `int` | Closure-cache capacity (default 8). |
| `_closure_max_depth` | `int` | Actual corpus depth. |
| `_seqid_y_map` | `dict[str,int]` | seqid → R-tree y-band. |
| `directives` | `list[str]` | Ordered `##` header lines. |

#### 5.6.3 `__getitem__` / `__contains__` — point lookup by ID

```python
def __getitem__(self, key) -> Feature        # interface.py:187-194
def __contains__(self, key) -> bool           # interface.py:196-201
```
- `target_id = key.id if isinstance(key, Feature) else key` — both accept a `Feature` or a raw ID string.
- `__getitem__` runs `SELECT {_SELECT_FEATURE} FROM features WHERE id = ?` with `[target_id]`, `fetchone()`. If `None`, **raise `FeatureNotFoundError(target_id)`** (`interface.py:193`). Otherwise return `feature_from_row(row, dialect=self.dialect)`.
- `__contains__` runs `SELECT 1 FROM features WHERE id = ? LIMIT 1`; returns `row is not None`.

`_SELECT_FEATURE` (`interface.py:33-36`) is the exact 12-column projection, in this order, that `feature_from_row` expects:
```
id, seqid, source, featuretype, start, "end", score, strand, frame,
attributes_blob, extra_blob, file_order
```
**Gotcha:** the `end` column is a DuckDB reserved word and **must** be quoted as `"end"` everywhere. Every SQL string in the module that references it quotes it.

**Gotcha (proxy contract):** `_RefDbProxy.__getitem__` and `_LFeatureDbProxy.__getitem__`/`_MFeatureDbProxy.__getitem__` (`locus_pipeline.py:173-176, 204-207, 250-253`) **raise `KeyError`** on a miss, *not* `FeatureNotFoundError`. The legacy consumer at `run_liftoff.py:252-253` catches `(KeyError, FeatureNotFoundError)`, so both are acceptable; a reimplementer of either side must keep the union catchable.

#### 5.6.4 Scans — `all_features` / `features_of_type`

```python
def all_features(self, limit=None, strand=None, featuretype=None,
                 order_by=None, reverse=False, completely_within=False)
                 -> Iterator[Feature]                    # interface.py:230-245
def features_of_type(self, featuretype, limit=None, strand=None,
                     order_by=None, reverse=False, completely_within=False)
                     -> Iterator[Feature]                # interface.py:247-260
```
`features_of_type` is a thin wrapper that delegates to `all_features` with `featuretype` forwarded (`interface.py:256-259`). Both build SQL via `_build_scan_sql` (`interface.py:262-296`) then stream through `_yield_features`.

**`_build_scan_sql` WHERE-clause construction**, in order:
1. **`limit` as a region triple** (`interface.py:268-279`): if `limit is not None`, parse it through `_normalize_region_args(limit, None, None, None)` → `(rseqid, rstart, rend)`. (Legacy `limit` is a *region restrictor*, not a row cap.) If `rseqid` present: `seqid = ?`. If both `rstart` and `rend` present:
   - `completely_within=True` → `start >= ? AND "end" <= ?` with params `[rstart, rend]`.
   - else (overlap) → `start <= ? AND "end" >= ?` with params `[rend, rstart]`. **Note the param swap** — overlap binds region-end to `start <= ?` and region-start to `"end" >= ?`.
2. **`strand`** (`interface.py:280-282`): `strand = ?`.
3. **`featuretype`** (`interface.py:283-290`): if a `list`/`tuple`/`set` → `featuretype IN (?,?,…)`; else `featuretype = ?`.
4. SQL: `SELECT {_SELECT_FEATURE} FROM features` + `WHERE ` + clauses joined by `" AND "` (only if any) + `ORDER BY ` + `_order_clause(order_by, reverse)`.

**`_order_clause(order_by, reverse)`** (`interface.py:298-313`):

| `order_by` | ORDER column |
|---|---|
| `None` | `file_order` (default — preserves input order) |
| `"length"` | `("end" - start)` |
| one of `seqid, source, featuretype, start, end, score, strand, frame, file_order` | that column (with `end` → `"end"`) |
| anything else | passed through verbatim (power-user literal SQL) |

Direction is `DESC` if `reverse` else `ASC`. **Gotcha:** the *default* ORDER BY is `file_order ASC` — this is what makes scans replay the input file order, which downstream byte-identity depends on.

#### 5.6.5 `region` — spatial overlap query with R-tree/B-tree dispatch

```python
def region(self, region=None, seqid=None, start=None, end=None,
           strand=None, featuretype=None, completely_within=False)
           -> Iterator[Feature]                          # interface.py:319-345
```
Steps:
1. `rseqid, rstart, rend = self._normalize_region_args(region, seqid, start, end)`.
2. `use_rtree = self._rtree_built and rseqid is not None and rstart is not None and rend is not None` (`interface.py:331-336`). The R-tree path is used **only** when all three bounds are known *and* the R-tree is live.
3. Build SQL via `_region_sql_rtree` or `_region_sql_btree`; stream via `_yield_features`.

**`_normalize_region_args(region, seqid, start, end)`** (static, `interface.py:414-443`) — returns `(seqid, start, end)` from four input shapes:

| Input shape | Result |
|---|---|
| `region` and any of `seqid/start/end` both given | `raise ValueError("pass either region or seqid/start/end, not both")` |
| `region is None` | `(seqid, start, end)` straight through |
| `str` `"chr"` | `(chr, None, None)` |
| `str` `"chr:start-end"` or `"chr:start-end:strand"` | split on first `:`, then split rest on first `-`; strip optional `:strand` suffix from end; `(chr, int(start), int(end))` |
| `tuple` length 3 | `(t[0], int(t[1]), int(t[2]))` |
| `tuple` length 2 | `(t[0], None, None)` |
| `tuple` other length | `raise ValueError` |
| `Feature` | `(f.seqid, f.start, f.end)` |
| other | `raise TypeError` |

**`_region_sql_btree`** (`interface.py:377-412`) — the universal path:
1. `seqid = ?` if `seqid` given.
2. If both `start` and `end` given: `completely_within` → `start >= ? AND "end" <= ?` (`[start, end]`); else overlap → `start <= ? AND "end" >= ?` (`[end, start]`).
3. Else if only `start` given → `start >= ?`; else if only `end` → `"end" <= ?`.
4. `strand = ?`, `featuretype = ?`/`IN (…)` as in scans.
5. `SELECT {_SELECT_FEATURE} FROM features` + optional `WHERE` + **`ORDER BY start`**.

**`_region_sql_rtree`** (`interface.py:347-375`) — spatial path:
1. `seqid_y = self._seqid_y_map.get(seqid)`. If `None` (seqid not in this DB), **delegate to `_region_sql_btree`** (returns a valid empty-row query; `interface.py:351-354`).
2. WHERE: `seqid = ?` **and** `ST_Intersects(bbox, ST_MakeEnvelope(?, ?, ?, ?))`, with params `[seqid, start, seqid_y, end, seqid_y + 1]`. **Gotcha:** each seqid occupies its own 1-unit-tall y-band `[seqid_y, seqid_y+1]`, so the spatial envelope is tight on both axes and never returns cross-chromosome candidates.
3. `completely_within` adds `start >= ? AND "end" <= ?` (`[start, end]`).
4. `strand`/`featuretype` filters as above.
5. `ORDER BY start`.

**Gotcha:** both `region()` paths order by `start ASC` (not `file_order`), unlike `all_features`. A reimplementer must not unify the two — `region()`'s ordering is positional.

There is also a vectorised `region_batched(regions, featuretype=None, completely_within=False, format="arrow")` (`interface.py:449-560`) that performs a single spatial JOIN between a registered staging Arrow table `__staging_regions` and `features`, returning a column-oriented PyArrow/pandas/polars result keyed by `query_idx`. It is **not used by the LiftOn pipeline** (LiftOn uses the row-by-row generators), but its result schema is the 14-column shape in `_empty_region_batched` (`interface.py:562-592`). Documented here for completeness; reimplementers targeting only LiftOn may omit it.

#### 5.6.6 `children` / `parents` — closure-cache vs dynamic-CTE traversal

```python
def children(self, id, level=None, featuretype=None, order_by=None,
             reverse=False, limit=None, completely_within=False)
             -> Iterator[Feature]                        # interface.py:598-612
def parents(self, id, level=None, featuretype=None, order_by=None,
            reverse=False, completely_within=False, limit=None)
            -> Iterator[Feature]                         # interface.py:614-628
```
**Gotcha:** the parameter order differs — `children` has `(…, reverse, limit, completely_within)` while `parents` has `(…, reverse, completely_within, limit)`. This matches gffutils. Both accept a `Feature` or raw ID (`target_id = id.id if isinstance(id, Feature) else id`), then call `_relation_query(...)` with `direction="children"` or `"parents"`.

**`_relation_query`** (`interface.py:837-860`):
1. `use_dynamic = self._dispatch_relation(level, target_id, direction)`.
2. Build SQL via `_relation_sql_dynamic` or `_relation_sql_cached`.
3. Stream via `_yield_features`.

**`_dispatch_relation(level, target_id, direction)`** (`interface.py:862-895`) — returns `True` ⇒ use the dynamic recursive CTE; `False` ⇒ use the materialised closure table. Decision matrix (verbatim, `interface.py:867-877`):

| `level` | `closure_max_depth` | `_has_overflow?` | Path |
|---|---|---|---|
| `None` | `== 0` (no edges materialised) | — | dynamic CTE |
| `None` | `>= 1` | no | closure cache |
| `None` | `>= 1` | yes | dynamic CTE |
| `int` | `<= _max_depth` | — | closure cache |
| `int` | `> _max_depth` | — | dynamic CTE |

Concretely (`interface.py:889-895`):
- If `level is not None`: return `level > self._max_depth`.
- Else if `self._closure_max_depth == 0`: return `True`.
- Else: return `self._has_overflow(target_id, direction)`.

**`_has_overflow(target_id, direction)`** (`interface.py:897-916`): checks whether the requested anchor has descendants/ancestors *past* `_max_depth`, i.e. the closure cache is truncated for this anchor. For `children`:
```sql
SELECT EXISTS (SELECT 1 FROM edges e JOIN closure c ON c.descendant = e.parent
               WHERE c.ancestor = ? AND c.depth = ?)
```
For `parents`, the symmetric query joins `c.ancestor = e.child WHERE c.descendant = ?`. Params `[target_id, self._max_depth]`. **Gotcha:** the cache is *preferred whenever it can serve the request* — Phase 7's measurement (GENCODE v45) found the cache 3.85× faster than forced-dynamic; dynamic is only the correctness fallback for traversals past `_max_depth`.

**`_relation_sql_cached`** (`interface.py:918-944`) — set-based JOIN on `closure`:
- For `children`: `join_col = c.descendant`, `anchor_col = c.ancestor`. For `parents`: reversed.
- WHERE: `{anchor_col} = ?` (param `target_id`); if `level is not None`: `c.depth = ?`; then featuretype filter; then limit filter.
- SQL: `SELECT {_select_feature_aliased('f')} FROM closure c JOIN features f ON f.id = {join_col} WHERE … ORDER BY {_order_clause_qualified(order_by, reverse, 'f')}`.

**`_relation_sql_dynamic`** (`interface.py:946-992`) — recursive CTE over `edges`:
- For `children`: seed `SELECT child, 1 FROM edges WHERE parent = ?`, recurse `JOIN edges e ON e.parent = w.id`. For `parents`: seed on `child = ?`, recurse on `e.child = w.id`, select `parent`.
- Walk bound `max_walk = level if level is not None else max(64, self._max_depth * 4)` (`interface.py:961`). The recursion stops at `WHERE w.depth < ?` bound by `max_walk`.
- CTE shape:
```sql
WITH RECURSIVE walk(id, depth) AS (
    SELECT <child|parent>, 1 FROM edges WHERE <parent|child> = ?
    UNION ALL
    SELECT e.<child|parent>, w.depth + 1
    FROM walk w JOIN edges e ON e.<parent|child> = w.id
    WHERE w.depth < ?)
SELECT {f-aliased cols} FROM walk w JOIN features f ON f.id = w.id
WHERE 1=1 [AND w.depth = ?] [featuretype] [limit] ORDER BY {order}
```
Params: `[target_id, max_walk]`, then `level` (if `level is not None`, appended as `w.depth = ?`), then featuretype params, then limit params.

**Shared filter helpers** (used by both cached and dynamic paths and by `region_batched`):
- `_featuretype_filter(featuretype, qualifier="f")` (`interface.py:1002-1009`): `None` → `([],[])`; list/tuple/set → `[f"{q}.featuretype IN (?,?,…)"]` + list of values; scalar → `[f"{q}.featuretype = ?"]` + `[value]`.
- `_limit_filter(limit, completely_within, qualifier="f")` (`interface.py:1011-1027`): parse `limit` through `_normalize_region_args`; emit `{q}.seqid = ?`, and the same overlap/within `start`/`"end"` predicates as the scan path (param swap on overlap).
- `_select_feature_aliased(alias)` (`interface.py:994-1000`): the 12-column projection, each column prefixed with `alias.` (with `alias."end"`).
- `_order_clause_qualified(order_by, reverse, qualifier)` (`interface.py:1029-1043`): identical table to `_order_clause` but every column is prefixed with `{qualifier}.` (default `{qualifier}.file_order`). **Gotcha:** `children`/`parents` default to `f.file_order ASC` — descendant rows come back in the *file order they were ingested*, not positionally. The proxy in `locus_pipeline.py` reproduces this because the materialiser builds its caches with explicit `order_by='start'` at every call site (see §5.6.7), and the legacy code always passes `order_by='start'` where order matters.

`children_batched` / `parents_batched` (`interface.py:642-685`) are the vectorised analogues that return PyArrow/pandas/polars tables via `_batched_relation` (`interface.py:687-765`); they use the identical cache-vs-dynamic decision (one-time, not per-row) and the recursive CTE seeded by *every* anchor in the batch. Not on LiftOn's hot path; documented for API completeness.

#### 5.6.7 The exact `children()` signatures LiftOn uses, and the proxy contract

The per-locus body (`run_liftoff.process_liftoff` and helpers) issues **only these four `children()` shapes** against the Liftoff DB, plus one against the miniprot DB. The proxy classes in `locus_pipeline.py` hard-code dispatch on exactly these, raising `NotImplementedError` for any other shape so a future new call site is caught immediately. A reimplementer of `FeatureDB.children` must return, for each, the row set described:

| # | Call (against `l_feature_db`) | Cache field (`_FeaturePreFetch`) | Proxy dispatch (`locus_pipeline.py`) |
|---|---|---|---|
| 1 | `children(feature, featuretype="exon", level=1, order_by="start")` | `exon_children_l1` | `featuretype=="exon" and level==1` → `:219-220` |
| 2 | `children(feature, level=1)` (any type) | `children_l1` | `featuretype is None and level==1` → `:225-226` |
| 3 | `children(feature, featuretype="exon", order_by="start")` (no level) | `exon_children_full` | `featuretype=="exon" and level is None` → `:221-222` |
| 4 | `children(feature, featuretype=("CDS","stop_codon"), order_by="start")` | `cds_stop_children` | `featuretype==("CDS","stop_codon")` → `:223-224` |
| 5 (miniprot) | `children(m_entry, featuretype=("CDS","stop_codon"), order_by="start")` | `_MiniprotPreFetch.cds_stop_children` | `_MFeatureDbProxy.children`, `:261-262` |

**Gotcha (ordering byte-identity):** every shape that LiftOn cares about passes `order_by="start"`. The real `FeatureDB` then orders by `f.start ASC`. The proxy returns `iter(cached_list)`, where the cached list was built by the parent thread calling the real DB with the same `order_by="start"` (`_walk_and_cache_features`, `locus_pipeline.py:353-407`). Therefore the proxy's iteration order is byte-identical to the real DB's. If a reimplementer changes `children()`'s tie-breaking for equal `start`, both the real and proxied paths shift together only if the cache is rebuilt — but the *default* ORDER BY (`file_order`) is never used at these call sites, so the relevant ordering is always `start ASC`.

**Gotcha (de-dup optimisation, `locus_pipeline.py:384-396`):** Phase 17c short-circuits shape #3 — when `exon_children_l1` is non-empty, `exon_children_full` is set to a *copy* of it without re-querying, on the proven invariant that for transcript-shaped features all exons are at level 1 (leaf-exon convention), so shape #3 ⊇ shape #1 with equality. The real DB query for #3 is only issued when there are no level-1 exon children. A `FeatureDB` reimplementer must ensure that, for the GFF3 hierarchies LiftOn produces (gene → mRNA → exon leaves), `children(feature, featuretype="exon")` and `children(feature, featuretype="exon", level=1)` return the same rows in the same order; otherwise the proxy path diverges from a real-DB path.

The proxies also reproduce `__getitem__` (returning the cached `Feature` / `_RefFeatureStub`) and raise `NotImplementedError` for un-cached `children()` signatures (`locus_pipeline.py:231-235, 263-267`). `_RefDbProxy` returns a `_RefFeatureStub` exposing only `.id` and `.attributes` — the only members `run_liftoff` reads from `ref_db[id]` (`locus_pipeline.py:179-187`).

#### 5.6.8 `iter_by_parent_childs` — grouped parent+children iteration

```python
def iter_by_parent_childs(self, featuretype="gene", level=None,
                          order_by=None, reverse=False,
                          completely_within=False) -> Iterator[List[Feature]]
```
(`interface.py:1360-1368`). For each `parent` in `features_of_type(featuretype, order_by=order_by, reverse=reverse)`, materialise `kids = list(children(parent, level=level, order_by=order_by, reverse=reverse, completely_within=completely_within))` and **yield `[parent, *kids]`** — a list whose first element is the parent and the rest are its descendants. Default `featuretype="gene"`. The grouping order follows the parent scan order (`file_order ASC` by default).

#### 5.6.9 Counts, distinct iterators, escape hatch, maintenance

| Method | Signature / line | Semantics |
|---|---|---|
| `count_features_of_type` | `(featuretype=None)`, `interface.py:207-212` | `None` → `SELECT COUNT(*) FROM features`; else `… WHERE featuretype = ?`. Returns the scalar `int`. |
| `featuretypes` | `interface.py:214-218` | yields each `SELECT DISTINCT featuretype FROM features ORDER BY featuretype`. |
| `seqids` | `interface.py:220-224` | yields each `SELECT DISTINCT seqid FROM features ORDER BY seqid`. |
| `execute` | `(query)`, `interface.py:1374-1378` | `self.conn.execute(query.rstrip(";"))` — returns DuckDB's relation cursor. The arbitrary-SQL escape hatch. **Gotcha:** strips a single trailing `;` only. |
| `analyze` | `interface.py:1380-1382` | `ANALYZE`; sets `_analyzed_flag = True`. |
| `set_pragmas` | `(pragmas)`, `interface.py:1384-1391` | iterate `pragmas.items()`, run `PRAGMA {k} = {v}`; **silently skip** any that raise `duckdb.Error` (legacy callers pass SQLite-only pragmas like `journal_mode`). |
| `schema` (property) | `interface.py:170-177` | concatenated DDL of all tables + views in the current database. |
| `_analyzed` (property) | `interface.py:179-181` | returns `_analyzed_flag`. |

`_yield_features(sql, params)` (`interface.py:1397-1404`) is the shared streaming primitive: execute, then loop `cur.fetchmany(10_000)` and `yield feature_from_row(row, dialect=self.dialect)` per row until exhausted. **Gotcha:** the batch size is 10000; this bounds peak memory on large scans but does not affect ordering.

`_read_meta` (`interface.py:149-154`): `SELECT key, value FROM meta` → `{k: v}`; returns `{}` on `duckdb.Error` (tolerates a DB without a `meta` table). `_has_rtree_index` (`interface.py:156-164`) and `_parse_dialect` (`interface.py:140-147`) as described in §5.6.2.

#### 5.6.10 Mutation & synthesis (not on the LiftOn read path)

The class also exposes write/synthesis methods that LiftOn does not exercise on the lift-over hot path, but which a faithful reimplementation should provide for drop-in parity. These all mutate `conn`:

| Method | Line | Behaviour summary |
|---|---|---|
| `delete(features, …)` | `1049-1064` | Coerce IDs via `_coerce_ids`; `DELETE FROM features/attributes/edges/closure` for those IDs. Returns `self`. |
| `update(data, …)` | `1066-1122` | Append `Feature`/`ParsedFeature`s via `_ArrowBatchBuilder` seeded with `_seqid_y_map`; refresh closure (`DELETE FROM closure` then re-run `CLOSURE_RECURSIVE_CTE` with `[_max_depth]`). New IDs default to `f"{featuretype}_{order}"`. |
| `add_relation(parent, child, level=1, …)` | `1124-1159` | `INSERT INTO edges`; optional `parent_func`/`child_func` callbacks; rebuild closure from edges via `CLOSURE_RECURSIVE_CTE`. |
| `interfeatures(features, …)` | `1181-1214` | yield gap `Feature`s between consecutive features (`new_start = prev.end+1`, `new_end = cur.start-1`, skip if `new_end < new_start`); optionally merge+dedupe attributes. |
| `merge(features, merge_criteria=None, multiline=False)` | `1216-1240` | sort by `(seqid, start, end)`; accumulate while all `merge_criteria` predicates hold (default `seqid, overlap_end_inclusive, strand, feature_type`); extend `accum.end = max(accum.end, f.end)`. |
| `merge_all`, `create_introns`, `create_splice_sites`, `children_bp`, `bed12` | `1251-1358` | gffutils-compat synthesis helpers (intron/splice-site generation, exonic-bp count, BED12 emission). |

**`CLOSURE_RECURSIVE_CTE`** (`schema.py:261-272`) is the canonical closure-materialisation statement these methods reuse: a recursive walk seeded from `edges` (`SELECT parent, child, 1`), recursing `JOIN edges e ON e.parent = w.descendant WHERE w.depth < ?`, inserting `(ancestor, descendant, depth)` rows; the `?` bound is `_max_depth`.

#### 5.6.11 `Feature` reconstruction & UCSC bin computation

**`feature_from_row(row, dialect=None)`** (`feature.py:477-497`) converts a 12-tuple (the `_DB_ROW_FIELDS` order, `feature.py:471-474`) into a `Feature`:
- `score`/`strand`/`frame` default to `"."` when the DB value is `None`.
- `attributes` is passed the **raw `attributes_blob` bytes** (or `None`) — parsing is *deferred* (`_LazyAttributes`, `feature.py:110-216`); a row that is never inspected never decodes col-9.
- `extra` is the `extra_blob` (or `None`); `dialect` defaults to `{"fmt": "gff3"}`.

`Feature` (`feature.py:218-456`) is 1-based-inclusive; key fields: `seqid, source, featuretype, start, end, score, strand, frame, attributes (_LazyAttributes), extra, id, dialect, file_order`. `__len__ = end - start + 1` (0 if either bound is `None`). `__str__`/`_format_line` (`feature.py:397-406`) reconstruct the 9 GFF3 columns + `extra`, with col-9 emitted **byte-faithfully from the original blob** when the attributes mapping was never materialised or mutated (`feature.py:357-366`). `__eq__`/`__hash__` are defined on `_format_line()` output (`feature.py:330-342`).

**UCSC binning (`_bins.py`).** Used *only* by the legacy SQLite export path (`Feature.astuple`/`calc_bin`); the runtime DuckDB query layer never reads a `bin` column (`_bins.py:4-7`). Constants:
```python
_BINOFFSETS  = (512+64+8+1, 64+8+1, 8+1, 1, 0)   # = (585, 73, 9, 1, 0)
_BINFIRSTSHIFT = 17
_BINNEXTSHIFT  = 3
```
`bin_from_coords(start, end)` (`_bins.py:16-30`) — algorithm (1-based inclusive input):
1. `start0 = start - 1`; `end0 = end` (convert to 0-based half-open).
2. `start_bin = start0 >> 17`; `end_bin = (end0 - 1) >> 17`.
3. For each `offset` in `_BINOFFSETS`: if `start_bin == end_bin`, **return `offset + start_bin`**; else `start_bin >>= 3`, `end_bin >>= 3`.
4. Fallthrough: return `0`.
This returns the smallest UCSC bin fully containing `[start, end]`.

#### 5.6.12 Exception hierarchy (`exceptions.py`)

| Class | Base | Constructor / attributes | Raised by |
|---|---|---|---|
| `GFFFormatError` | `ValueError` | `(message="", *, line_no=0, kind="")`; stores `.message`, `.line_no`, `.kind` | parser on a spec-violating line (Rust `._native` variant when extension loaded) |
| `FeatureNotFoundError` | `Exception` | `(feature_id)`; message `"feature not found: {id}"`; stores `.feature_id` | `FeatureDB.__getitem__` on missing ID (`interface.py:193`) |
| `DuplicateIDError` | `Exception` | — | ingestion with `merge_strategy='error'` |
| `AttributeStringError` | `Exception` | — | malformed col-9 attributes |
| `EmptyInputError` | `Exception` | — | input file/iterable yields no features |

**Gotcha:** `GFFFormatError` deliberately subclasses `ValueError` so legacy `pytest.raises(ValueError)` and `except ValueError` callers keep working; a reimplementation that makes it a bare `Exception` would break those. `FeatureNotFoundError` is *not* a `ValueError` and is *not* a `KeyError` — but consumers that mix the real DB and the proxy must catch `(KeyError, FeatureNotFoundError)` because the proxies raise `KeyError` while the real DB raises `FeatureNotFoundError` (§5.6.3).

#### 5.6.13 Dialect template (`dialect.py`)

`default_dialect()` (`dialect.py:13-25`) returns the canonical GFF3 dialect dict consumed by `Feature._format_attributes` (`feature.py:357-395`):

| Key | Default | Meaning |
|---|---|---|
| `"fmt"` | `"gff3"` | format flavour |
| `"field separator"` | `";"` | between key=value pairs in col-9 |
| `"keyval separator"` | `"="` | between key and value (GTF uses a space) |
| `"multival separator"` | `","` | between repeated values of one key |
| `"leading semicolon"` / `"trailing semicolon"` | `False` | edge punctuation |
| `"quoted GFF2 values"`, `"repeated keys"`, `"semicolon in quotes"` | `False` | GTF/edge-case flags |
| `"order"` | `[]` | first-appearance key order |

`merge_dialects(samples)` (`dialect.py:28-62`) reconciles per-line observations: format is `"gtf"` iff GTF lines strictly outnumber GFF3 lines; `keyval separator` follows the chosen format; `field separator` is the plurality vote (`max(set(seps), key=seps.count)`); the five boolean flags are OR'd across samples; `order` is first-appearance union of each sample's `order`. Empty input ⇒ `default_dialect()`.

#### 5.6.14 Summary of byte-identity-critical invariants for reimplementers

1. **Column projection order** must be exactly `_SELECT_FEATURE` (§5.6.3); `feature_from_row` unpacks positionally.
2. **Default ORDER BY is `file_order ASC`** for scans and relation queries; **`region()` orders by `start ASC`**. Do not unify.
3. **`children()`/`parents()` return generators of `Feature`** (lazy attributes), streamed in 10000-row batches; the proxy returns `iter(list)` — same iteration order, same `Feature` shape.
4. **The four (+1 miniprot) `children()` shapes in §5.6.7** are the entire contract LiftOn depends on; `order_by="start"` at every site is what makes the proxy and the real DB interchangeable.
5. **`FeatureNotFoundError` vs `KeyError`** must both be catchable by the consumer; do not narrow.
6. **`"end"` must be SQL-quoted** everywhere; it is a reserved word in DuckDB.
