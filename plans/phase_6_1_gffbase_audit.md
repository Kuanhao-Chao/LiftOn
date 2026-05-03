# LiftOn — Phase 6.1: `gffbase` Ingestion Audit

> **Status:** vendored package landed at `lifton/gffbase/`. Public API
> imports cleanly (`from lifton.gffbase import FeatureDB, create_db,
> Feature, …`); Rust extension `_native.abi3.so` loads under the new
> namespace; functional smoke test ingests `hierarchy.gff3` end-to-end;
> all 299 LiftOn tests still pass.

---

## 1. What was copied (and what was excluded)

### Copied
| Path under `lifton/gffbase/` | Source (gffbase repo) | Notes |
|---|---|---|
| `__init__.py`, `_bins.py`, `create_db.py`, `dialect.py`, `exceptions.py`, `feature.py`, `gffwriter.py`, `helpers.py`, `ingest.py`, `interface.py`, `iterators.py`, `merge_criteria.py`, `parser.py`, `schema.py`, `sqlite_export.py` | `python/gffbase/` | The 15 Python modules (~3.6 KLOC). `from gffbase.X` imports rewritten to relative `from .X` |
| `_pyfallback/{__init__,attributes,parser}.py` | `python/gffbase/_pyfallback/` | Pure-Python correctness oracle |
| `_native.abi3.so` | `python/gffbase/_native.abi3.so` | Pre-built arm64 macOS Rust extension; loads without rebuild |
| `_rust/Cargo.toml`, `_rust/Cargo.lock`, `_rust/src/*.rs` | `rust/` (sans `target/`) | 6 source files — kept for forward maintenance / rebuild |
| `data/{hierarchy,simple}.gff3`, `data/{simple,synthesize}.gtf` | `tests/data/` | 4 small fixtures used by `helpers.example_filename` |
| `LICENSE` | `LICENSE` | MIT (upstream attribution) |
| `MIGRATION.md` | `MIGRATION.md` | In-tree gffutils → gffbase migration reference |

**Total footprint:** 788 KB.

### Excluded
`.git/`, `.github/`, `.gitignore`, `.pytest_cache/`, `.coverage`,
`coverage.xml`, `htmlcov/`, `site/` (rendered MkDocs output), `dist/`
(wheels), `bench/`, `benchmarks/`, `plans/` (gffbase's own planning
docs), `docs/` (MkDocs source), `mkdocs.yml`, `MANIFEST.in`,
`pyproject.toml` (gffbase's own — replaced by entries in LiftOn's
`setup.py`), `PERFORMANCE_COMPARISON.md`, the `tests/` suite (we'll
write LiftOn-specific tests in Phase 6.2), and `rust/target/` (213 MB
of build artifacts).

---

## 2. Import-rewrite summary

Every internal `from gffbase.<x> import …` was rewritten to a relative
import so the package is fully self-contained under `lifton.gffbase`:

| Pattern (before) | Pattern (after) |
|---|---|
| `from gffbase.feature import Feature` (top-level module) | `from .feature import Feature` |
| `from gffbase._native import GFFFormatError` | `from ._native import GFFFormatError` |
| `from gffbase._pyfallback.attributes import …` | `from ._pyfallback.attributes import …` |
| `from gffbase.exceptions import …` (inside `_pyfallback/`) | `from ..exceptions import …` |
| `from gffbase import merge_criteria as mc` | `from . import merge_criteria as mc` |

Verified: `grep -rEn "(from gffbase|^import gffbase)" lifton/gffbase/` returns 0 hits.

The Rust `module-name = "gffbase._native"` setting in `_rust/Cargo.toml`
is now stale relative to the install path, but the prebuilt
`_native.abi3.so` works fine because PyO3 only uses the module name at
the Python init site (`#[pymodule] fn _native(...)`) and the function
name `_native` matches the file's basename. **Future rebuild caveat:**
when running `maturin develop` from `_rust/`, set
`module-name = "lifton.gffbase._native"` in `_rust/Cargo.toml` (or
`MATURIN_PYTHON_PATH=lifton/gffbase`). The current binary is stable;
flagging this only because the Rust `Cargo.toml` is now an
out-of-context fragment.

---

## 3. Public API surface (drop-in for `gffutils`)

Re-exported from `lifton/gffbase/__init__.py`:

| Symbol | Source module | Purpose |
|---|---|---|
| `create_db(path, dbfn=…, force=…, **kw)` | `create_db.py` | Build a DuckDB-backed feature database from a GFF3/GTF file (drop-in for `gffutils.create_db`) |
| `FeatureDB(dbfn)` | `interface.py` | Query interface — see method table below |
| `Feature`, `ParsedFeature` | `feature.py` | Public feature class + slotted dataclass emitted by the parser |
| `DataIterator(...)` | `iterators.py` | Stream features from a GFF/GTF file without persisting to DB |
| `GFFWriter(...)` | `gffwriter.py` | Write features back to GFF3 |
| `merge_criteria.*` | `merge_criteria.py` | Predicates for `FeatureDB.merge()` |
| `parse_gff`, `parse_bytes`, `detect_dialect`, `native_available` | `parser.py` | Lower-level parser (gffbase extras) |
| `ingest`, `export_sqlite` | `ingest.py`, `sqlite_export.py` | Direct ingest API + legacy-SQLite export |
| `FeatureNotFoundError`, `DuplicateIDError`, `AttributeStringError`, `EmptyInputError`, `GFFFormatError` | `exceptions.py` | Drop-in exception names |
| `example_filename(name)` | `helpers.py` | Returns absolute path to a bundled fixture under `lifton/gffbase/data/` |
| `__version__` | constant | `"0.1.0"` |

### `FeatureDB` method inventory (40 public methods)

Mirroring legacy `gffutils.FeatureDB`:

| Category | Methods |
|---|---|
| **Indexing / metadata** | `__getitem__`, `__contains__`, `count_features_of_type`, `featuretypes`, `seqids`, `schema`, `analyze`, `set_pragmas`, `execute` |
| **Sequential scan** | `all_features`, `features_of_type` |
| **Spatial query** | `region`, `region_batched` (auto-routes between R-tree and B-tree based on `meta.rtree_built`) |
| **Hierarchy traversal** | `children`, `parents`, `children_batched`, `parents_batched` (closure-cache vs recursive-CTE dispatch) |
| **Mutation** | `delete`, `update`, `add_relation` |
| **Derivation** | `interfeatures`, `merge`, `merge_all`, `create_introns`, `create_splice_sites`, `children_bp`, `bed12`, `iter_by_parent_childs` |

---

## 4. Side-by-side: `gffutils` calls in LiftOn vs `gffbase` equivalents

Inventory of every `db_connection.<method>` call in `lifton/`:

| Call site | Legacy `gffutils` call | `gffbase` equivalent | Compatibility |
|---|---|---|---|
| `lifton/annotation.py:20` | `gffutils.FeatureDB(file_name)` | `gffbase.FeatureDB(dbfn)` | ✅ same signature; backend swap is transparent |
| `lifton/annotation.py:37, 60` | `gffutils.create_db(file, dbfn, merge_strategy=, id_spec=, force=, verbose=, disable_infer_transcripts=, disable_infer_genes=, transform=)` | `gffbase.create_db(file, dbfn=…, force=…, …)` | ⚠ `disable_infer_transcripts` / `disable_infer_genes` → check ingest module flag names; `transform` callback may need adapter (see §5) |
| `lifton/annotation.py:83, 94, 123, 146` | `db.features_of_type(featuretype=…)` | `db.features_of_type(featuretype, …)` | ✅ |
| `lifton/annotation.py:85, 96, 130, 131, 136, 137` | `db.children(feature, featuretype=…, level=…)` | `db.children(feature, featuretype=…, level=…)` | ✅ same signature |
| `lifton/annotation.py:171` | `db.parents(feature_name)` | `db.parents(feature_name)` | ✅ |
| `lifton/annotation.py:159, 167, 175, 184, 187` | `db[feature_name]`, `db.children(...)` | `db[feature_name]`, `db.children(...)` | ✅ |
| `lifton/extract_sequence.py:30, 40, 41, 54` | `ref_db.db_connection.features_of_type(...)`, `.children(feature, featuretype=…, level=…)` | identical method names | ✅ |
| `lifton/lifton_utils.py:326, 327, 348, 352, 356, 369` | `ref_db.db_connection.features_of_type(...)`, `.children(...)` | identical | ✅ |
| `lifton/lifton.py:285, 290, 291, 327, 339, 340, 352` | `tgt_feature_db.features_of_type(feature, limit=…)`, `m_feature_db.features_of_type('mRNA')`, `l_feature_db.features_of_type(feature)` | `db.features_of_type(featuretype)` + spatial filter via `db.region(seqid, start, end, featuretype=…)` for the `limit=` kwarg | ⚠ `limit=("chr1", a, b)` keyword needs mapping (see §5) |
| `lifton/run_miniprot.py:90, 91, 120` | `m_feature_db.children(m_feature, featuretype='CDS')` | identical | ✅ |
| `lifton/run_liftoff.py:70, 73` | `l_feature_db.children(locus, featuretype='exon', order_by='start')`, `featuretype=('CDS', 'stop_codon')` | identical (`order_by`, tuple featuretype both supported) | ✅ |
| `lifton/lifton_utils.py:256-300` (`LiftOn_miniprot_alignment`) | `m_feature_db[m_id]`, `m_feature_db.children(m_entry, featuretype=('CDS', 'stop_codon'), order_by='start')` | identical | ✅ |
| `lifton/run_evaluation.py` | various `db_connection.children(...)` | identical | ✅ |

### Method-signature parity quick reference

| Method | Args we pass today | gffbase support |
|---|---|---|
| `features_of_type` | `featuretype` (str or tuple), `limit=(seqid,start,end)` (rare) | `featuretype`: ✅ both forms. `limit`: see §5 — gffbase routes spatial filtering through `region()`, not a `limit=` kwarg on `features_of_type` |
| `children` | `feature_or_id`, `featuretype=…`, `level=…`, `order_by=…` | ✅ all four kwargs |
| `parents` | `feature_or_id`, optionally `level=…` | ✅ |
| `__getitem__` | feature id (str) | ✅ raises `FeatureNotFoundError` on miss (same exception class name) |
| `__contains__` | feature id | ✅ |
| `merge_strategy="create_unique"` (in `create_db`) | string | ✅ — but verify supported values in gffbase ingest |

---

## 5. Behaviour deltas to validate during the swap (Phase 6.2)

These are NOT broken yet — flagging now so Phase 6.2 covers them:

1. **`create_db` keyword fidelity.** Legacy `gffutils.create_db` accepts
   `disable_infer_transcripts`, `disable_infer_genes`, `id_spec`,
   `transform`, `merge_strategy`. gffbase has `_synthesize_transcripts`
   / `_synthesize_genes` helpers (`ingest.py:455, 474`); confirm the
   public `create_db` exposes equivalent flags. If `transform=` (a
   callback for mutating each parsed feature) is missing, our
   `lifton/annotation.py:72-76` `transform_func` will need a
   pre/post-pass via `db.update()` instead.

2. **`features_of_type(limit=("chr1", a, b))`.** Used in the commented
   code path `lifton/lifton.py:290, 339`. If a future test enables
   that kwarg, gffbase steers spatial filtering through `db.region()`.
   Drop-in adapter:
   ```python
   def features_of_type_compat(db, featuretype, limit=None):
       if limit is None:
           yield from db.features_of_type(featuretype)
       else:
           seqid, start, end = limit
           yield from db.region(seqid=seqid, start=start, end=end,
                                featuretype=featuretype)
   ```

3. **`gffutils.FeatureDB(<existing-db-file>)`.** Currently
   `lifton/annotation.py:20` tries to attach an existing
   `<gff>_db` cache. gffbase will accept either a `.duckdb` path or a
   live `duckdb.DuckDBPyConnection`. The legacy SQLite caches will NOT
   be readable directly — gffbase reads only DuckDB-format files. We
   either (a) re-ingest, (b) run `export_sqlite` in reverse via
   gffbase's import path (not supplied; would need conversion script),
   or (c) keep gffutils alongside for one release and migrate
   gradually. **Recommendation:** plan a one-shot migration in Phase
   6.2 — produce `.gff3.duckdb` siblings the first time we see an
   `.gff3_db` legacy cache.

4. **`feature.attributes`.** Legacy `gffutils.Feature.attributes` is a
   dict-like (`{key: [values]}`) that supports both `feat.attributes['ID']`
   and assignment. gffbase's `Feature.attributes` is a lazy proxy
   (`feature.py:165` decodes `attributes_blob` on demand) but exposes
   the same `MutableMapping` interface. Pinning behaviour in tests
   during Phase 6.2 will catch any subtle divergence (e.g., does
   reassignment persist back into `attributes_blob`?).

5. **`feature.id`, `feature.seqid`, `feature.start`, `feature.end`,
   `feature.strand`, `feature.frame`, `feature.featuretype`.** All
   present on `gffbase.Feature`. ✅

6. **String coercion `str(feature)`.** Legacy gffutils renders the
   GFF3 line via `__str__`. gffbase preserves this for round-trip
   fidelity (the `attributes_blob` path) — verify in Phase 6.2 that
   `Lifton_GENE.write_entry` (which does `fw.write(str(self.entry) + "\n")`)
   produces byte-identical GFF3 lines vs the gffutils output.

7. **`Lifton_TRANS.add_cds` mutates `gffutil_entry_cds.attributes['Parent']`.**
   At `lifton/lifton_class.py:254`. gffbase's lazy `attributes`
   property must support assignment (`feat.attributes['Parent'] = ['tx1']`)
   and re-serialise on `str(feat)`. Phase 6.2 must add a regression
   test specifically for this mutation pattern.

---

## 6. Dependency impact

`gffbase` introduces two new runtime dependencies:

| Dependency | Lower bound | Why |
|---|---|---|
| `duckdb>=1.0` | as upstream gffbase requires | columnar storage backend |
| `pyarrow>=14` | as upstream gffbase requires | Rust → Python → DuckDB record-batch hand-off |

Added to `setup.py:install_requires` (LiftOn doesn't have a
`pyproject.toml` yet — Phase 4 roadmap §3.1 deferred that migration).
Also added `package_data` entries so `pip install .` ships the
`_native.abi3.so`, `data/*.gff3`, `data/*.gtf`, the Rust source
fragments under `_rust/`, and the LICENSE.

`gffbase` does NOT pull in `gffutils` itself, so we can decommission
`gffutils>=0.10.1` from `install_requires` once Phase 6.2 finishes
the cutover. Holding off until then to keep the bridge tested.

### Disk footprint

| Item | Size |
|---|---|
| `lifton/gffbase/` total | 788 KB |
| `lifton/gffbase/_native.abi3.so` | ~600 KB (the bulk) |
| `lifton/gffbase/_rust/` | 64 KB (Rust source only — `target/` excluded) |
| `lifton/gffbase/data/` | 12 KB (4 small fixture files) |

---

## 7. Verification (executed at audit time)

```bash
source /opt/anaconda3/etc/profile.d/conda.sh && conda activate lifton-test
pip install --quiet duckdb pyarrow

# Import test
python -c "
from lifton.gffbase import FeatureDB, Feature, create_db, native_available
print('native_available:', native_available())   # True
"

# End-to-end ingest + query
python -c "
from lifton.gffbase import create_db, example_filename
db = create_db(example_filename('hierarchy.gff3'), dbfn='/tmp/h.duckdb', force=True)
print(list(db.seqids()))                # ['chr1']
print(list(db.featuretypes()))          # ['CDS', 'exon', 'gene', 'mRNA']
print(db.count_features_of_type('gene'))  # 1
"

# LiftOn suite remains green
pytest tests/ -q                        # 299 passed
```

All three checks pass at the time this audit was written.

---

## 8. What Phase 6.2 will do (preview, not part of this phase)

1. Add `tests/test_gffbase_smoke.py` to pin the gffbase API surface
   we depend on.
2. Build a thin adapter layer in `lifton/annotation.py` that swaps
   `gffutils` calls for the corresponding `gffbase` calls behind a
   single feature flag (`LIFTON_USE_GFFBASE`).
3. Run the full LiftOn suite with the flag flipped; resolve the
   six items in §5 one by one.
4. Once all 299 tests pass under gffbase, remove `gffutils` from
   `install_requires` and delete the legacy code paths.
5. Output expected: byte-identical lifton.gff3 vs the pre-cutover
   golden output, plus a substantial wall-clock improvement on the
   `--measure_time` step that builds the gffutils SQLite DB
   (`lifton/lifton.py:238, 311, 314`).

**Phase 6.1 stops here. No LiftOn logic has been modified.**
