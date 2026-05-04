# LiftOn — Phase 8: Vendored-Liftoff In-Memory Refactor Execution Report

> **Status:** **365 tests passing, 0 failures, 0 xfailed.** Output GFF3
> byte-identical across the full 2×2 flag matrix (`--stream={off,on}`
> × `--inmemory-liftoff={off,on}`). Phase 7 baseline preserved;
> coverage gates intact on every core module.

---

## 1. What landed

### Production code (4 files)

| File | Change | Purpose |
|---|---|---|
| `lifton/liftoff/liftoff_main.py` | refactored `run_all_liftoff_steps`; new `run_all_liftoff_steps_inmemory` + private `_run_liftoff_pipeline` | Decompose Liftoff orchestration so the legacy entrypoint and the in-memory entrypoint share the same pipeline body. The legacy entrypoint still calls `write_new_gff` after the pipeline returns; the in-memory entrypoint stops at the four-tuple `(lifted_feature_list, feature_db, ref_parent_order, unmapped_features)`. |
| `lifton/liftoff/inmemory_emitter.py` (NEW, 56 LOC) | new module | `lifted_features_to_gff3_bytes(...)` mirrors `write_new_gff.write_new_gff` exactly but emits to a `StringIO` buffer and returns UTF-8 bytes. Reuses `liftoff_utils.get_parent_list`, `write_new_gff.finalize_parent_features`, `build_parent_dict`, `write_feature` so the byte-output is provably identical. |
| `lifton/run_liftoff.py` | +25 LOC | Branches on `args.inmemory_liftoff`: when set, calls `run_all_liftoff_steps_inmemory(...)`, runs the bytes through `inmemory_emitter`, and returns `bytes`. When unset (default), preserves the Phase 5 path that returns a file path. |
| `lifton/lifton.py` | +9 LOC | New `--inmemory-liftoff` CLI flag (default off); flows through `args.inmemory_liftoff` to the call site at `lifton.py:362`. |

### Test code (1 new file, 19 new tests)

| File | Tests | Layer |
|---|---:|---|
| `tests/test_liftoff_inmemory.py` | 19 | 5 unit (emitter/header/entrypoint surface), 4 unit (run_all_liftoff_steps_inmemory shape + legacy-still-writes-disk), 5 integration (each cell of the 2×2 matrix individually + the full 4-cell byte-identity gate), 3 unit (CLI flag plumbing), 2 differential (`run_liftoff` return shape switches between path and bytes under the flag) |

The roadmap target was **~10 unit + ~3 integration + ~5 differential = ~18 new tests**. We landed **19** (target met).

---

## 2. Test counts

```
Phase 5 baseline:           297 passed
Phase 6.1 smoke:             +2          = 299
Phase 7 streaming adapter:   +47         = 346
Phase 8 (this phase):        +19         = 365
                            ─────────────────
Total:                      365 passed, 0 failed, 0 xfailed
```

Wall-clock: **55 s** for the full suite on the `lifton-test` conda env.
Roadmap target: **~356 tests**. Achieved: **365** (+9 over).

---

## 3. Coverage report

```
Name                        Stmts   Miss  Cover
-----------------------------------------------
lifton/gffbase_adapter.py      25      2    92%
lifton/lifton_class.py        671     69    90%
lifton/lifton_utils.py        284     17    94%
lifton/run_liftoff.py          89     18    80%   ← was 71%, +9pp from
                                                    new in-memory branch
TOTAL (these modules)        1069    106    90%
```

Gates required: `lifton_class.py ≥ 90 %`, `lifton_utils.py ≥ 90 %`. **Both met.**

`lifton/run_liftoff.py` jumped to 80 % (from 71 % in Phase 7) because
the new differential tests in `test_liftoff_inmemory.py::TestRunLiftoffReturnShape`
exercise both the path and bytes return-shape branches.

---

## 4. Golden output gate — 2×2 flag matrix

CLI-driven end-to-end against the synthetic chr1 fixture used by
`test_integration_pipeline.py`:

```
GOLDEN OUTPUT — 2x2 flag matrix
=======================================================
  stream=0,inmem=0  →  391 bytes  [OK]
  stream=0,inmem=1  →  391 bytes  [OK]
  stream=1,inmem=0  →  391 bytes  [OK]
  stream=1,inmem=1  →  391 bytes  [OK]
=======================================================
Verdict: BYTE-IDENTICAL ACROSS ALL 4 COMBINATIONS ✅
```

The `tests/test_liftoff_inmemory.py::TestInmemoryLiftoffMatrix::test_all_four_combinations_byte_identical`
test enforces the same gate inside the suite for CI, so any future
divergence trips a hard test failure.

---

## 5. Phase 8 contract upheld

| Contract | Status |
|---|---|
| `--stream` and `--inmemory-liftoff` change I/O, not algorithms | ✅ verified by 2×2 matrix |
| Legacy entrypoint (`run_all_liftoff_steps`) still produces a disk file | ✅ test `test_legacy_entrypoint_still_writes_disk` |
| In-memory entrypoint (`run_all_liftoff_steps_inmemory`) skips the polish-intermediate disk write | ✅ test `test_returns_four_tuple` asserts `polish_intermediate_write=False` |
| Emitter byte-equivalence with `write_new_gff` | ✅ test `test_emitter_matches_disk_write_byte_for_byte` |
| Hermetic pipeline (no real `minimap2`/`miniprot` invoked) | ✅ Phase 7's `hermetic_pipeline` fixture preserved across all 2×2 cells |
| NCBI compliance (`--strict-gff` validator unchanged) | ✅ ingest path is gffbase + the validator runs at Step 0.5 regardless of flags |

---

## 6. What gets eliminated under the streamed in-memory path (`--stream --inmemory-liftoff`)

| Disk artefact | Phase 5 default | Phase 8 streamed |
|---|---|---|
| `<liftoff_outdir>/liftoff.gff3` (Liftoff write) | written, then re-parsed | **never written** |
| `<liftoff_outdir>/liftoff.gff3.duckdb` (gffbase cache) | n/a | n/a (in-memory dbfn) |
| `<liftoff_outdir>/miniprot/miniprot.gff3` (miniprot write) | written, then re-parsed | **never written** |
| `<liftoff_outdir>/miniprot/miniprot.gff3.duckdb` (gffbase cache) | n/a | n/a (in-memory dbfn) |
| `<liftoff_outdir>/unmapped_features.txt` (Liftoff unmapped log) | written | written (kept — diagnostic value) |
| `<final_output>/lifton_output/lifton.gff3` (final result) | written | written (the user's answer) |

So the streamed in-memory path leaves only the unmapped-features
diagnostic log and the final answer on disk. Two intermediate GFF3s
and zero SQLite caches removed.

---

## 7. Verification commands

```bash
source /opt/anaconda3/etc/profile.d/conda.sh && conda activate lifton-test

# Full suite
pytest tests/ -q                                            # 365 passed

# Phase 8 specifically
pytest tests/test_liftoff_inmemory.py -q                     # 19 passed

# Coverage gates
coverage run --source=lifton --omit="lifton/liftoff/*,lifton/gffbase/*" \
    -m pytest tests/ -q
coverage report --include="lifton/lifton_class.py,lifton/lifton_utils.py" \
    --fail-under=90                                          # 90/94 — both pass

# 2x2 golden output gate (CLI-driven)
python <<'PY'
... (driver shown in §4 above) ...
PY
# Output: BYTE-IDENTICAL ACROSS ALL 4 COMBINATIONS ✅
```

---

## 8. Acceptance against the roadmap

| Roadmap gate (Phase 6.4.B) | Required | Achieved |
|---|---|---|
| All 297 + 49 (Phase 7) + new tests green | yes | **365 passed, 0 failed** |
| Coverage on `lifton_class.py` ≥ 90 % | yes | 90 % |
| Coverage on `lifton_utils.py` ≥ 90 % | yes | 94 % |
| Byte-identical output across full 2×2 matrix | yes | ✅ 391 bytes ×4 |
| `lifton/liftoff/inmemory_emitter.py` shipped | yes | ✅ 56 LOC |
| `run_all_liftoff_steps_inmemory` shipped | yes | ✅ legacy `run_all_liftoff_steps` preserved |
| `--inmemory-liftoff` CLI flag (default off) | yes | ✅ |
| Phase 8 report written | yes | this file |

Phase 8 deliverable complete. **Default behaviour is unchanged** — the
new flag is opt-in. The Phase 5 → Phase 7 → Phase 8 progression has
removed both intermediate GFF3 disk materialisations without ever
breaking the byte-identity contract.

---

## 9. Surface for Phase 6.4.C (next sub-phase)

Phase 8 closes the I/O cliff. Phase 6.4.C (locus-major fusion +
parallelism) starts from a clean state where:
- The reference / liftoff / miniprot DBs all live in RAM under `--stream --inmemory-liftoff`.
- `Annotation` polymorphism on bytes is established.
- The `run_all_liftoff_steps_inmemory` four-tuple is the natural
  hand-off to a per-locus pipeline that consumes
  `lifted_feature_list` lazily.

**Awaiting your approval to begin Phase 6.4.C (locus-major fusion +
parallelism).**
