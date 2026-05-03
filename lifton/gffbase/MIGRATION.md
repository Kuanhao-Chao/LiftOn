# Migrating from `gffutils` to `gffbase`

GFFBase is a drop-in successor to legacy
[`gffutils`](https://github.com/daler/gffutils). For most users, the
migration is one import change.

> ## ⚠️  READ THIS FIRST — the OLAP/OLTP gotcha
>
> **There is exactly one common code pattern that gets *slower*, not
> faster, when you migrate to gffbase.** It's the per-id Python loop:
>
> ```python
> # ❌ ANTI-PATTERN with gffbase: 50 000 small queries.
> # Pays DuckDB's vectorization startup × 50 000 + per-row Feature
> # construction × 1.6 M. ≥ 10 minutes wall on GENCODE v49.
> for transcript_id in fifty_thousand_transcript_ids:
>     for exon in db.children(transcript_id, featuretype="exon"):
>         starts.append(exon.start)
>         ends.append(exon.end)
> ```
>
> DuckDB is an **OLAP** engine — designed for big set-based queries.
> Iterating it row-by-row pays vectorization startup *per call* and
> never amortizes. SQLite (legacy gffutils) is **OLTP** — its B-tree
> seek on a cache-warm file is microseconds per call.
>
> ### ✅ The fix — one canonical PyArrow snippet
>
> ```python
> # ✅ ONE set-based SQL query for all 50 000 transcripts.
> # Returns a zero-copy pyarrow.Table — no `Feature` object is ever
> # constructed. 1.16 s wall on GENCODE v49 → a 553× speedup.
> exons = db.children_batched(
>     fifty_thousand_transcript_ids,
>     featuretype="exon",
>     format="arrow",         # or "df" / "polars"
> )
>
> # NumPy / PyTorch / JAX / Hugging Face datasets — all native.
> starts = exons.column("start").to_numpy()
> ends   = exons.column("end").to_numpy()
>
> # The "anchor" column carries the input transcript_id for each row,
> # so you can groupby in Python or downstream Arrow tooling without
> # re-issuing N queries:
> import pyarrow.compute as pc
> per_tx_exon_count = pc.value_counts(exons.column("anchor"))
> ```
>
> If your code has a `for x in ids: db.children(x, …)` loop and you
> care about wall time, **convert it now**, before you migrate. This
> is the only common change required of legacy `gffutils` users.

---

## 1. Drop-in compatibility — the easy part

Every public surface from legacy `gffutils` is preserved verbatim:

| `gffutils` symbol | `gffbase` equivalent |
|---|---|
| `gffutils.create_db(path, dbfn, ...)` | `gffbase.create_db(path, dbfn, ...)` |
| `gffutils.FeatureDB(dbfn)` | `gffbase.FeatureDB(dbfn)` |
| `gffutils.Feature(...)` | `gffbase.Feature(...)` |
| `gffutils.DataIterator(...)` | `gffbase.DataIterator(...)` |
| `gffutils.GFFWriter(...)` | `gffbase.GFFWriter(...)` |
| `gffutils.merge_criteria.*` | `gffbase.merge_criteria.*` |
| `gffutils.example_filename(name)` | `gffbase.example_filename(name)` |
| Exceptions (`FeatureNotFoundError`, …) | same names |

```python
# Before
import gffutils
db = gffutils.create_db("annotation.gff3", "annotation.db")

# After
import gffbase as gffutils      # one-line alias migration
db = gffutils.create_db("annotation.gff3", "annotation.duckdb")
```

All `FeatureDB` methods (`children`, `parents`, `region`,
`features_of_type`, `interfeatures`, `merge`, `bed12`, `update`,
`delete`, `add_relation`, `execute`, …) accept the same arguments and
return generators of `Feature` objects — identical to the legacy API.

The **storage backend** changes (DuckDB instead of SQLite). This is
transparent for almost all callers, but raw SQL queries that hit the
legacy schema directly via `db.execute(...)` need rewriting against
the GFFBase schema (or against the SQLite-compat views; see §4). We
also ship `gffbase.export_sqlite(con, path)` to dump a GFFBase
database into a legacy `.sqlite` file when you need the old format.

---

## 2. What you gain immediately, no code changes

The internal Big Four mega-bench, with the v0.1.0 GFF3 ingest optimizations
applied — head-to-head against legacy `gffutils`:

| Corpus                  | Format | gffbase ingest | legacy ingest | speedup       | spatial qps (gffbase R-tree) | batched 5 k anchors |
| ----------------------- | :----: | -------------: | ------------: | ------------: | ---------------------------: | ------------------: |
| GENCODE v49 (basic)     |  GTF   |   4 min 37 s   | ≥ 2 hr 30 min | **🚀 ≥ 32×**  |                    **1,204** | 172 ms / 596 k desc |
| GENCODE v49 (basic)     |  GFF3  |   6 min 7 s    |  11 min 23 s  |    **1.86×**  |                    **1,292** | 422 ms / 1.93 M desc|
| RefSeq GRCh38.p14       |  GFF3  |   4 min 12 s   |     6 min 5 s | **1.45×**     |                    **1,011** | 263 ms / 999 k desc |
| MANE v1.5 (Ensembl)     |  GFF3  |       21.6 s   |        45.1 s | **2.09×**     |                    **1,766** |  78 ms / 156 k desc |
| CHESS 3.1.3             |  GFF3  |       53.6 s   |   2 min 13.1 s| **2.48×**     |                    **1,175** |  91 ms / 161 k desc |

| Single-call workload                                  | Speedup      |
|---|---|
| Spatial overlap (`db.region(...)`, p50 latency)       | **8.35× lower** (0.72 ms vs 6.01 ms) |
| `db.children(id, level=1)` indexed lookup             | comparable   |
| `db.children_batched(ids, format="arrow")` (50 k ids) | **🚀 36.68×** |

Your existing `gffutils` script gets the ingest, spatial, and
attribute-query wins the moment you swap the import. To unlock the
36.68× ML-batched win, see the warning at the top of this page.

---

## 3. ⚠️  Deep-dive: the OLAP vs OLTP tradeoff

DuckDB is an **OLAP** engine. It's optimized for big set-based queries
(JOINs, aggregations, scans of millions of rows). SQLite is an **OLTP**
engine — optimized for tiny indexed point lookups against cache-warm
pages. **For tiny, repeated point queries against a cache-warm DB,
SQLite (and therefore legacy gffutils) is faster.**

The fix is the canonical PyArrow snippet at the top of this page.
Full benchmark: at 50 000 transcripts the row-by-row gffbase loop is
≥ 10 minutes; the batched call is **1.16 s — a 553× speedup over the
gffbase loop, 36.68× over legacy** (`PERFORMANCE_COMPARISON.md` §4b).

### Vectorized methods at a glance

| Vectorized method | Replaces this loop |
|---|---|
| `db.children_batched(ids, level=…, featuretype=…, format='arrow')` | `for x in ids: db.children(x, …)` |
| `db.parents_batched(ids, …, format='arrow')` | `for x in ids: db.parents(x, …)` |
| `db.region_batched(regions, …, format='arrow')` | `for r in regions: db.region(r, …)` |

`format` accepts `"arrow"` (default — `pyarrow.Table`), `"df"`
(`pandas.DataFrame`), or `"polars"` (`polars.DataFrame`). All three
share memory with DuckDB's query buffers — no per-row Python
materialization happens at any layer.

### When you don't need to migrate the pattern

- One-off scripts that ask `db[gene_id]` or `db.children(gene)` for
  fewer than ~100 anchors.
- Small annotations (< 100 k features) where SQL startup overhead is
  not visible.

For everything else — ML feature extraction, BED12 export of every
transcript, "for each peak in this 50 000-row BED file find every
overlapping CDS" — switch to `*_batched`.

---

## 4. SQL-compat views (for raw `execute()` users)

Legacy code that did `db.execute("SELECT * FROM features WHERE …")`
hits the new DuckDB schema (`features`, `attributes`, `edges`,
`closure`). Two compatibility views provide the legacy column shapes:

```sql
-- features_compat: legacy SQLite-style 12-column features table.
SELECT * FROM features_compat WHERE seqid = 'chr1' LIMIT 5;

-- relations_compat: legacy parent/child/level table.
SELECT parent, child, level FROM relations_compat WHERE level = 1;
```

The `attributes` column on `features_compat` is the **raw col-9
bytes** (UTF-8), not legacy-style JSON. If your raw-SQL code parses
JSON out of that column, switch to querying the normalized
`attributes` table directly:

```sql
SELECT a.value FROM attributes a
WHERE a.feature_id = ? AND a.key = 'gene_biotype';
```

This is also faster — `attributes_kv` indexes `(key, value)`, so
attribute filters become indexed seeks.

---

## 5. SQLite export — the safety valve

If a downstream tool only knows how to read legacy
`gffutils`-compatible SQLite files:

```python
from gffbase import export_sqlite
export_sqlite(db.conn, "legacy_compatible.sqlite")
```

Produces a SQLite database with the original gffutils schema,
populated UCSC `bin` column, and the closure flattened back into
`relations(parent, child, level)`. The downstream tool can open this
file with `gffutils.FeatureDB("legacy_compatible.sqlite")`.

---

## 6. Things that changed (small list)

- **Storage backend**: SQLite → DuckDB. Database file extension is
  `.duckdb` by convention. The legacy SQLite layout is reachable via
  `export_sqlite()` (above) or the compat views.
- **Disk size**: GFFBase databases are ~1.5× larger than legacy
  SQLite (the price of materializing the closure + R-tree). Worth it
  for the 5–550× query speedups.
- **Peak ingest RSS**: ~10× higher (~1.6 GB vs ~150 MB on GENCODE
  v49). DuckDB allocates a vectorized ingest buffer pool; reduce with
  `PRAGMA memory_limit='512MB'` if needed.
- **Hierarchy depth**: GFFBase materializes the closure to depth 8 by
  default (vs depth 2 in legacy). Anything past 8 falls through to a
  dynamic recursive CTE — the dispatcher is automatic.
- **Attributes column shape**: in raw SQL, the legacy single-cell
  JSON blob is replaced by a normalized
  `attributes(feature_id, key, value, idx)` long-form table. Filtering
  by attribute is now an indexed query, not a full scan.
- **Duplicate IDs**: NCBI RefSeq emits multiple GFF3 rows that share
  `ID=cds-NP_xxx`. gffbase auto-suffixes repeats with `__N` (mirroring
  `gffutils.merge_strategy="create_unique"`) and records the remap in
  the `duplicates` table. No config change required.

---

## 7. Migration checklist

- [ ] `pip install gffbase`
- [ ] Replace `import gffutils` with `import gffbase as gffutils` (or
      use the new name directly).
- [ ] Re-ingest your annotations (`create_db`) — old `.sqlite` files
      can still be read by legacy gffutils; they're not GFFBase
      databases.
- [ ] **Audit your code for `for x in ids: db.children(x, …)` loops
      and convert them to `db.children_batched(ids, format='arrow')`.**
      This is the only common change that requires user action.
- [ ] If you have raw `db.execute(...)` SQL: use
      `features_compat` / `relations_compat` views, or move attribute
      filters onto the normalized `attributes` table.
- [ ] Run your existing test suite. Everything else should be
      identical.

If anything breaks, please open an issue at
<https://github.com/your-org/gffbase/issues> with a minimal
reproducer.
