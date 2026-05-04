# LiftOn — Phase 6.2: `gffutils` → `gffbase` Migration Strategy

## Context

Phase 6.1 vendored `lifton/gffbase/` (788 KB; Rust-accelerated GFF3
parser with DuckDB-backed `FeatureDB`) and the audit confirmed
**99 %** of LiftOn's database calls (`features_of_type`, `children`,
`parents`, `__getitem__`) are signature-compatible with the legacy
`gffutils.FeatureDB`. We must now actually swap backends without
breaking the Phase 5 zero-bug baseline (**297 tests, 0 failures**).

The risk is concentrated in one place: `lifton/annotation.py` is the
sole direct importer of `gffutils`. Every other LiftOn module
consumes a `FeatureDB`-like object that already came out of
`Annotation.db_connection`. That gives us a clean single-cut surface,
which is why the strategy below is **a swap-in backend behind a
feature flag, not a parallel rewrite**.

The intended outcome of Phase 6.2 is: with `LIFTON_USE_GFFBASE=1` set
(or `--backend gffbase` on the CLI) all 297 tests still pass and the
end-to-end golden-path GFF3 output is byte-identical to the gffutils
baseline. Once that is demonstrated, Phase 6.3 will flip the default
and drop `gffutils` from `setup.py:install_requires`.

## Critical files

| File | Role | What changes in Phase 6.2 |
|---|---|---|
| `lifton/annotation.py` | The **only** importer of `gffutils`; wraps a `FeatureDB` and exposes it via the `db_connection` property | Add `_build_database_gffbase()` parallel to the existing 3-strategy `_build_database()`; route via `args.backend` / env-var |
| `lifton/gffbase_adapter.py` (NEW) | Tiny shim that maps the LiftOn-specific `gffutils` calls onto `gffbase` semantics where they diverge (transform callback, DB file extension, merge strategy semantics, id_spec) | Single new file, ~120 LOC, no logic — just translation |
| `lifton/lifton.py` | Adds `--backend {gffutils,gffbase}` CLI flag (default `gffutils`) and threads it into every `annotation.Annotation(...)` constructor call | Three call sites: `lifton.py:336, 362, 375` |
| `lifton/extract_sequence.py`, `lifton/lifton_utils.py`, `lifton/run_miniprot.py`, `lifton/run_liftoff.py`, `lifton/run_evaluation.py`, `lifton/intervals.py` | Consumers of `db_connection.X(...)` | **No source edits** — they receive whichever backend `Annotation` chose |
| `tests/test_gffbase_swap.py` (NEW) | Mirror of every existing LiftOn integration test, parametrised over backend | New test file ~200 LOC |
| `tests/conftest.py` | Already has `hermetic_pipeline` and the small fixtures | Add a `backend` parametrise fixture |

### Reusable functions / utilities (no need to re-implement)

- `lifton.gffbase.create_db(...)` — already accepts the legacy
  `disable_infer_genes`, `disable_infer_transcripts`, `force`,
  `verbose`, `merge_strategy`, `id_spec`, `transform`, `from_string`,
  `gtf_subfeature` keyword arguments. The gffbase docstring explicitly
  notes that `merge_strategy`, `id_spec`, `transform` are
  "accepted-but-no-op for now" — Phase 6.2 plumbs them.
- `lifton.gffbase.FeatureDB(dbfn)` — accepts a path, a live
  `duckdb.DuckDBPyConnection`, or `(con, stats)`. Maps 1:1 to
  `gffutils.FeatureDB(dbfn)`.
- `lifton.gffbase.example_filename(name)` — already returns
  `lifton/gffbase/data/<name>` (verified Phase 6.1).
- `lifton.gffbase.Feature` — exposes `id`, `seqid`, `start`, `end`,
  `strand`, `frame`, `featuretype`, `attributes` (mutable
  `_LazyAttributes`), `__str__()` (verified: `g.attributes['Parent']
  = ['demo']` mutates and `str(g)` round-trips).

### Confirmed compatibility (verified Phase 6.1, no work needed)

The following gffbase behaviours match gffutils as LiftOn calls them:

| LiftOn call | gffbase support |
|---|---|
| `db.features_of_type('gene')` / `db.features_of_type(featuretype='gene')` | ✅ both forms |
| `db.children(feat, featuretype='exon', level=1, order_by='start')` | ✅ all four kwargs |
| `db.children(feat, featuretype=('CDS','stop_codon'), order_by='start')` | ✅ tuple featuretype supported |
| `db.parents(feat, level=1)` | ✅ |
| `db[feature_id]` raising `FeatureNotFoundError` on miss | ✅ same exception name |
| `db.featuretypes()`, `db.seqids()`, `db.count_features_of_type()` | ✅ |
| `feat.attributes['Parent'] = [...]` (mutation) followed by `str(feat)` (round-trip) | ✅ verified |

## Strategy: Adapter / shim with a feature flag

### A) Adapter module — `lifton/gffbase_adapter.py` (new, ~120 LOC, behaviour-only)

```python
"""Phase 6.2 shim — translate the few LiftOn-specific gffutils calls
into gffbase semantics. Pure translation; no LiftOn logic lives here."""

from __future__ import annotations
import os, sys
from typing import Optional, Callable

from lifton import gffbase as _gffbase


def db_path_for(file_name: str) -> str:
    """LiftOn caches gffutils SQLite at <file>_db. gffbase emits a
    DuckDB file; use a different suffix so the two backends can
    coexist on disk without clobbering each other."""
    return file_name + ".duckdb"


def open_existing_db(file_name: str) -> Optional[_gffbase.FeatureDB]:
    """Equivalent of `try: gffutils.FeatureDB(<file>_db)`. Returns
    None if no on-disk DuckDB cache exists yet."""
    p = db_path_for(file_name)
    if not os.path.exists(p):
        return None
    try:
        return _gffbase.FeatureDB(p)
    except Exception:
        return None


def build_database(
    *,
    file_name: str,
    infer_genes: bool,
    infer_transcripts: bool,
    merge_strategy: str,
    id_spec: Optional[str],
    force: bool,
    verbose: bool,
    transform: Optional[Callable] = None,
) -> _gffbase.FeatureDB:
    """The three-strategy retry from lifton/annotation.py:_build_database
    collapsed to a single call. gffbase's ingestor already deduplicates
    IDs internally; the legacy `create_unique` retry is therefore not
    strictly needed, but we still pass merge_strategy through for
    forward compatibility."""
    db = _gffbase.create_db(
        file_name,
        dbfn=db_path_for(file_name),
        force=force,
        verbose=verbose,
        merge_strategy=merge_strategy,
        id_spec=id_spec,
        disable_infer_genes=not infer_genes,
        disable_infer_transcripts=not infer_transcripts,
        transform=transform,
    )
    return db
```

The shim is **deliberately tiny** — three functions, all pure
translation. No LiftOn algorithmic logic moves into it.

### B) Backend selection in `lifton/annotation.py`

Single decision point inside `Annotation.__init__` (around line 64
where `gffutils.FeatureDB` is first called). The selection is driven
by, in priority order:

1. The `backend` kwarg passed to `Annotation(...)` directly (used by
   tests).
2. `args.backend` if `args` is supplied (used by the CLI flow).
3. The `LIFTON_USE_GFFBASE` env var (set to anything truthy).
4. Default: `"gffutils"` (preserves Phase 5 baseline byte-for-byte).

```python
class Annotation:
    def __init__(self, file_name, infer_genes, infer_transcripts,
                 merge_strategy, id_spec, force, verbose,
                 auto_convert_gtf=True, *, backend: Optional[str] = None):
        # ... existing init ...
        self.backend = self._resolve_backend(backend)
        self._get_db_connection()

    @staticmethod
    def _resolve_backend(explicit: Optional[str]) -> str:
        if explicit in ("gffutils", "gffbase"):
            return explicit
        if os.environ.get("LIFTON_USE_GFFBASE"):
            return "gffbase"
        return "gffutils"

    def _get_db_connection(self):
        if self.backend == "gffbase":
            from lifton import gffbase_adapter
            db = gffbase_adapter.open_existing_db(self.file_name)
            if db is None:
                db = gffbase_adapter.build_database(
                    file_name=self.file_name,
                    infer_genes=self.infer_genes,
                    infer_transcripts=self.infer_transcripts,
                    merge_strategy=self.merge_strategy,
                    id_spec=self.id_spec,
                    force=self.force,
                    verbose=self.verbose,
                    transform=self._get_transform_func(),
                )
            self._db_connection = db
            return
        # ── existing gffutils path unchanged ─────────────────────────
        ... existing 3-strategy retry ...
```

The gffutils branch is left **byte-for-byte unchanged** so the
default-on case can never regress. Only the new `gffbase` branch is
added.

### C) CLI flag in `lifton/lifton.py`

```python
parser.add_argument(
    '--backend', choices=['gffutils', 'gffbase'], default='gffutils',
    help='Annotation database backend. gffbase (DuckDB-backed, '
         'Rust-parsed) is faster and stricter; gffutils (SQLite) is '
         'the legacy default.'
)
```

Threaded into the three `annotation.Annotation(...)` calls at
`lifton/lifton.py:336, 362, 375` via a new `backend=args.backend`
kwarg.

## Test-driven swap plan

### Step 1 — Pin the gffbase API surface (no LiftOn changes yet)

New file `tests/test_gffbase_smoke.py` (~30 tests). Asserts the
gffbase methods + return-types we depend on. Each test **must pass
today against vanilla gffbase** (no LiftOn glue) and serve as the
contract our adapter relies on:

- `create_db(path, dbfn, force=True)` returns a `FeatureDB`.
- `db.features_of_type('gene')` yields `Feature`s with `.id`,
  `.seqid`, `.start`, `.end`, `.strand`, `.frame`, `.featuretype`,
  `.attributes`.
- `db.children(feat, featuretype='exon', level=1, order_by='start')`
  yields children in start order.
- `db.children(feat, featuretype=('CDS','stop_codon'), order_by='start')`
  honours tuple featuretypes.
- `db.parents(feat, level=1)` returns the gene.
- `db['gene1']` raises `FeatureNotFoundError` on miss.
- `feat.attributes['Parent'] = ['gA']` followed by `str(feat)` round-trips.
- `feat.attributes['extra_copy_number']` round-trips.
- `merge_strategy='create_unique'` accepted (currently no-op; pin it).
- `disable_infer_genes`, `disable_infer_transcripts` accepted.
- `transform=lambda x: x` accepted (no-op; pin it).
- `example_filename('hierarchy.gff3')` resolves under
  `lifton/gffbase/data/`.

If any of these fail, Phase 6.2 stops immediately and the gap goes
back to the gffbase repo.

### Step 2 — Build the adapter; backend defaults to gffutils

Implement `lifton/gffbase_adapter.py` and the `_resolve_backend` /
`_get_db_connection` branch in `lifton/annotation.py`. No call-site
changes anywhere else. Existing 297 tests must still pass (default
backend unchanged). Verification:

```bash
pytest tests/ -v          # 297 passed; xfails: 0
```

### Step 3 — Parametrise integration tests across both backends

New file `tests/test_gffbase_swap.py`. Every test in it is the
gffbase-flag-on counterpart of an existing LiftOn integration test:

```python
@pytest.fixture(params=["gffutils", "gffbase"])
def backend(request):
    return request.param

def test_run_all_lifton_steps_golden_path(backend, integration_workspace,
                                          hermetic_pipeline):
    # Same body as test_integration_pipeline.py, but threads
    # backend=backend through args / Annotation.
    ...
    assert out_gff_gffbase_text == out_gff_gffutils_text   # byte-identical
```

The fixture-parametrise pattern keeps the legacy assertions intact
and adds a *byte-identical golden output* gate for the gffbase path.

### Step 4 — Run, fix, repeat

For each test that fails under `backend="gffbase"`, classify the
failure and route to one of:

| Failure category | Fix location |
|---|---|
| Missing kwarg / signature mismatch on a gffbase method | Adapter (`lifton/gffbase_adapter.py`) |
| Returned `Feature` lacks an attribute LiftOn reads | Adapter wraps the result with a `_FeatureView` proxy; do NOT modify gffbase |
| `transform=` callback not honoured | Pre-pass the GFF file through the callback line-by-line into a temp file, then ingest; document as a Phase 6.3 follow-up |
| `merge_strategy='create_unique'` no-op produces wrong dedup | Adapter emits its own pre-pass that renames duplicates (mirrors gffutils' `create_unique` strategy) |
| `id_spec` no-op produces wrong ID column | Adapter pre-rewrites attributes |
| Byte-identical-output diff | Source-of-truth comparison: gffutils output is the baseline; if gffbase emits semantically-equivalent but byte-different GFF, normalise both before comparing (e.g., sort attribute keys) |

The adapter is allowed to grow; LiftOn's algorithmic code is not.

### Step 5 — Default flip + cleanup

When the parametrised suite is fully green for **both** backends:

1. Change the default in `_resolve_backend` from `"gffutils"` to
   `"gffbase"`.
2. Re-run the full suite — all 297 (now ~327) tests must still pass.
3. Drop `gffutils>=0.10.1` from `setup.py:install_requires`.
4. Delete the gffutils branch in `Annotation._get_db_connection` and
   the three-strategy retry in `_build_database`.
5. Delete `lifton/gffbase_adapter.py` if everything has graduated to
   first-class gffbase calls; otherwise keep it as the only place
   that imports `lifton.gffbase`.

This step is **gated on user approval** — Phase 6.2 stops at the end
of Step 4 with the default still on gffutils.

## Validation plan

### V1 — Smoke test the gffbase API surface
Run `pytest tests/test_gffbase_smoke.py -v`. Expected: every
contract assertion passes against vanilla gffbase. If not, the
delta gets reported back to the gffbase repo as a bug; we don't
work around it inside LiftOn.

### V2 — Default-off regression check
Run `pytest tests/ -v` with no env var. Expected: **297 passed, 0
failed, 0 xfailed** — identical to the Phase 5 baseline because the
default backend is unchanged.

### V3 — Default-on parametrised swap check
Run `pytest tests/test_gffbase_swap.py -v` and
`LIFTON_USE_GFFBASE=1 pytest tests/test_integration_pipeline.py -v`.
Expected: every parametrised test passes for both backends.

### V4 — Byte-identical golden GFF3
On the chr22 fixture (and on the small synthetic fixture in
`tests/conftest.py:integration_workspace`), assert that the output
GFF3 produced via `--backend gffbase` is byte-identical to the
output produced via `--backend gffutils`. If gffbase emits a
semantically-equivalent but byte-different ordering (e.g., attribute
key order), introduce a normaliser BUT document the diff in
`plans/phase_6_2_diff_log.md`.

### V5 — Memory + wall-clock benchmark (informational, non-gating)
Use the existing `--measure_time` harness (`lifton/lifton.py:486-525`)
to record:

```
backend  step1_db_build  step5_db_build  total_walltime  peak_rss
gffutils ...
gffbase  ...
```

Phase 6.2 expectation per the gffbase MIGRATION.md numbers: the
two `Annotation(...)` calls that build a gffutils SQLite DB
(`lifton/lifton.py:362, 375`) drop from O(seconds) to O(milliseconds)
on real annotations, and total wall time on chr22 falls by ≥30 %.
This is informational, not a gate — Phase 6.2 ships even if the
numbers show only parity, as long as V2-V4 are green.

### V6 — Coverage gate maintained
`coverage report --include="lifton/lifton_class.py,lifton/lifton_utils.py"
--fail-under=90` — both core modules must remain ≥90 %.

## Rollback discipline

Each step lands as a separate commit on `devel`:

1. `Add gffbase smoke contract test (no logic changes)`
2. `Add gffbase backend adapter + opt-in flag (default off)`
3. `Parametrise integration tests over backend`
4. `Fix gffbase divergence: <category> #N` (one commit per category)
5. *(deferred to Phase 6.3)* `Flip default backend to gffbase + drop gffutils`

If any step breaks the V2 default-off check, that commit is reverted
immediately. The adapter is allowed to absorb arbitrary complexity
to keep V2 green; LiftOn's algorithmic code is not.

## Verification commands

```bash
source /opt/anaconda3/etc/profile.d/conda.sh && conda activate lifton-test

# V1 — smoke
pytest tests/test_gffbase_smoke.py -v

# V2 — default-off baseline (must remain green)
pytest tests/ -v                                         # 297 passed

# V3 — parametrised swap
pytest tests/test_gffbase_swap.py -v                     # all green
LIFTON_USE_GFFBASE=1 pytest tests/test_integration_pipeline.py -v

# V4 — byte-identical golden output
lifton --backend gffutils -g chr22.gff3 chr22.fa chr22.fa -o /tmp/utils.gff3
lifton --backend gffbase  -g chr22.gff3 chr22.fa chr22.fa -o /tmp/base.gff3
diff /tmp/utils.gff3 /tmp/base.gff3                      # empty

# V5 — bench
lifton --backend gffutils --measure_time ... && cat time.txt
lifton --backend gffbase  --measure_time ... && cat time.txt

# V6 — coverage gate
coverage run --source=lifton --omit="lifton/liftoff/*,lifton/gffbase/*" \
    -m pytest tests/ -q
coverage report --include="lifton/lifton_class.py,lifton/lifton_utils.py" \
    --fail-under=90
```

---

## Phase 6.2 deliverables (what this phase produces)

1. `lifton/gffbase_adapter.py` — new (~120 LOC).
2. `lifton/annotation.py` — `_resolve_backend` + opt-in `gffbase`
   branch in `_get_db_connection`. **No edits to the gffutils
   branch.**
3. `lifton/lifton.py` — `--backend` CLI flag + thread-through to the
   three `Annotation(...)` call sites.
4. `tests/test_gffbase_smoke.py` — new (~30 tests).
5. `tests/test_gffbase_swap.py` — new (~10 parametrised tests).
6. `plans/phase_6_2_diff_log.md` — itemised log of every behavioural
   delta we had to absorb in the adapter, with NCBI compliance
   notes.

Phase 6.2 STOPS before flipping the default. Phase 6.3 (separate
approval required) does the cleanup: flip default + drop gffutils.
