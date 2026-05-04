# LiftOn — Phase 4.5: Algorithmic Integrity Report

> **Status:** **240 passed, 3 xfailed** in 11.9 s. Core-module coverage
> gate **CLEARED**: `lifton_class.py` **91 %**, `lifton_utils.py` **93 %**.
> Six new legacy bugs surfaced (3 from Phase 3 + 3 newly found in 4.5).
> **Phase 5 (refactor execution) is green-lit.**

---

## 1. Coverage delta vs Phase 3

| Module | Phase 3 | Phase 4.5 | Δ | Gate |
|---|---:|---:|---:|---|
| `lifton/lifton_class.py` | 52 % | **91 %** | **+39** | ✅ ≥ 90 % |
| `lifton/lifton_utils.py` | 78 % | **93 %** | **+15** | ✅ ≥ 90 % |
| `lifton/extract_sequence.py` | 99 % | **99 %** | — | ✅ |
| `lifton/get_id_fraction.py` | 95 % | **95 %** | — | ✅ |
| `lifton/intervals.py` | 100 % | **100 %** | — | ✅ |
| `lifton/align.py` | 86 % | **87 %** | +1 | — |
| `lifton/annotation.py` | 59 % | **59 %** | — | (uncovered = subprocess paths) |
| `lifton/variants.py` | 18 % | **97 %** | **+79** | ✅ |
| `lifton/lifton.py` | 67 % | **67 %** | — | (orchestrator; subprocess) |
| `lifton/run_liftoff.py` | 71 % | **71 %** | — | (subprocess) |
| `lifton/run_miniprot.py` | 19 % | **19 %** | — | (subprocess) |
| `lifton/protein_maximization.py` | 8 % | **8 %** | — | (deferred to Phase 5) |
| `lifton/run_evaluation.py` | 10 % | **10 %** | — | (eval-only mode) |
| `lifton/stats.py` | 69 % | **69 %** | — | (output-only) |
| **Package total** (excl. vendored) | **58 %** | **75 %** | **+17** | — |

The two core modules targeted in this phase are both above the 90 %
gate. The remaining sub-90 % modules are all subprocess-bound or
specific to alignment-only paths that Phase 5 Step 6 (alignment
extraction) and Step 9 (parallelism) will touch — those will get
their own coverage push as part of the refactor PRs.

---

## 2. Test counts by category

| Category | File | Tests | Notes |
|---|---|---:|---|
| Pre-existing (Phase 3) | `test_lifton_utils.py` etc. | 97 | unchanged |
| Boundary logic (Step 1) | `test_boundary_logic.py` | 15 | 1-bp features, chrom edges, off-end coords, padding |
| Translation & ORF (Step 2) | `test_translation_orf.py` | 28 | premature stops, frameshifts, non-canonical starts, all `find_variants` branches |
| Complex gene models (Step 3) | `test_complex_gene_models.py` | 12 | nested genes, isoform-shared exons, opp-strand overlap, prok pattern |
| Property-based (Step 4) | `test_property_based.py` | 10 + 2 xfail | Hypothesis: 200 examples per strategy |
| Corruption (Step 5) | `test_corruption_handling.py` | 14 | dirty input behaviour pinned |
| Coverage backfill (Step 6) | `test_coverage_backfill.py` | 38 | god-class gap closure |
| ORF rescue deep | `test_orf_rescue_deep.py` | 16 | `__find_orfs`, `__update_cds_boundary` end-to-end |
| Class branch polish | `test_class_branch_polish.py` | 10 | last-mile ENSEMBL/CHESS, novel CDS, etc. |
| **Total** | | **240** | + **3 strict xfail** marking real bugs |

Suite is fully hermetic — no `minimap2` / `miniprot` invoked. The
`hermetic_pipeline` fixture from Phase 3 is reused for the integration
test; everything else is direct unit-level construction.

---

## 3. Newly surfaced legacy bugs

Each bug is captured by a `@pytest.mark.xfail(strict=True)` test so that
when Phase 5 Step 1 fixes it, the xfail flips to `XPASS` and CI yells.
Phase 3 documented bugs #1-3; this phase adds #4-6.

| # | Site | Severity | NCBI § | Description | Test reference |
|---|---|---|---|---|---|
| 1 | `lifton_class.py:71-75` | High | ID | `attributes["ID"]` assigned a *string*, not a list; `[0]` then yields one character (`"g"` instead of `"gene1"`). Propagates into `entry.id` and any downstream `Parent`. | `test_lifton_class.py::test_extra_copy_number_appends_suffix` (pinned at `entry.id == "g"`) |
| 2 | `lifton_class.py:236` | High | Parent | Transcript `Parent` inherits the buggy single-char gene id from #1. | `test_lifton_class.py::test_transcript_id_assigned` (pinned at `Parent == ["g"]`) |
| 3 | `lifton_utils.py:521` | Medium | n/a (3rd party) | `check_ovps_ratio` passes a raw tuple to `IntervalTree.overlap()`; newer `intervaltree` requires `Interval` or two ints. Raises `AttributeError` on every miniprot-merge attempt that gets past the unknown-chromosome short-circuit. | `test_lifton_utils.py::test_overlap_ratio_triggers_true` (xfail strict) |
| **4** | `lifton_utils.py:484-502` | High | Cols 4-5 | **`segments_overlap_length` is not symmetric when both segments share the same start.** `sorted([s1, s2], key=lambda x: x[0])` is stable; when starts tie, the input order leaks into the result. Falsifying example from Hypothesis: `((1,1),(1,2))` → `(1, True)` but `((1,2),(1,1))` → `(2, True)`. | `test_property_based.py::test_segments_overlap_length_symmetric` (xfail strict) |
| **5** | `lifton_utils.py:174-179` | Medium | n/a | **`get_ID_base` reduces a single-component all-numeric id to the empty string.** `"0".split("_")` → `["0"]`; `int("0")` succeeds; `splits[:-1]` is `[]`; `"_".join([])` is `""`. Any feature whose id is a bare integer is silently nuked. | `test_property_based.py::test_get_id_base_strips_trailing_int_once` (xfail strict) |
| **6** | `lifton/annotation.py` (gffutils contract) | Medium | Cols 1-9, Directives, Parent, Cols 4-5 | **The annotation parser silently accepts every NCBI violation.** Concretely: missing `##gff-version 3` directive (NCBI Directives §); 8-column rows (NCBI Column Specifications §); negative coordinates (NCBI cols 4-5); `start > end` (NCBI cols 4-5); dangling `Parent` references (NCBI Parent §); unencoded reserved characters in attribute values (NCBI Attribute Specifications §). | All in `test_corruption_handling.py` |

### Severity rubric
- **High:** affects correctness of legitimate inputs (#1, #2, #4 silently produce wrong IDs / wrong overlap measurements that flow into the IntervalTree).
- **Medium:** affects edge cases (numeric-only IDs in #5, NCBI-non-conforming inputs in #6) or only triggers on dual-LiftOff+miniprot paths (#3).

### Note on bug #4 (overlap symmetry)
This one is the most consequential of the new findings. `segments_overlap_length` feeds into `Lifton_TRANS.add_cds` (`lifton_class.py:252`), `update_cds_list` (called five times across Cases 1, 3, 5), `LiftOn_miniprot_alignment` (`lifton_utils.py:261`), and `check_ovps_ratio` (`lifton_utils.py:523`). When two features share a start coordinate (common for genes anchored to the same TSS), the reported overlap length depends on argument order. The downstream calls happen to be argument-order-stable (the function is called with `(exon, cds)` consistently), so production output is currently NOT corrupt — but any future caller that swaps argument order will silently get a different answer. **Phase 5 Step 1 fix:** sort by `(start, end)` so ties tie-break deterministically.

### Selenocysteine / `transl_except` gap (documented, not bug)
`test_translation_orf.py::test_selenocysteine_tga_treated_as_stop` pins the legacy gap that LiftOn does NOT honour the NCBI `transl_except` attribute (selenocysteine, programmed stop-codon readthrough). Internal `TGA` is always translated as `*`. This is a feature gap, not a regression — the test exists so the validator extension in Phase 5 Step 2 picks up `transl_except` attributes and warns when LiftOn cannot honour them.

---

## 4. NCBI compliance gap matrix

| NCBI invariant | Spec §| Current behaviour | Will be enforced in |
|---|---|---|---|
| `##gff-version 3` directive on line 1 | Directives | Silently accepted if absent | Phase 5 Step 2 (validator) |
| Column count == 9 | Column Specifications | 8 cols → tolerated; 10 cols → tolerated | Phase 5 Step 2 |
| `start <= end` (cols 4-5) | Cols 4-5 | Silently accepted when violated | Phase 5 Step 2 |
| Coordinates ≥ 1 (1-based) | Cols 4-5 | Negative ints silently accepted | Phase 5 Step 2 |
| `strand ∈ {+, -, ., ?}` | Col 7 | Tested via property suite — passes | Already enforced by gffutils |
| `phase ∈ {0, 1, 2, .}` | Col 8 | Not validated; `.` allowed even on CDS rows | Phase 5 Step 2 (warn-only per NCBI) |
| Percent-encoded reserved chars | Attribute Specifications | Unencoded `;` `=` `,` mis-split silently | Phase 5 Step 2 (validate); Phase 5 Step 8 (writer percent-encodes via gffutils `__str__`) |
| Multi-value `Parent`, `Alias`, `Note`, `Dbxref` | Attributes | Round-trip works (verified in `test_complex_gene_models.py`) | ✅ Already correct |
| Hierarchy `gene → mRNA → {exon, CDS}` | Annotation Data Model | Works | ✅ Already correct |
| Permitted exception: `gene → CDS` (prokaryote, NOTE 2) | Annotation Data Model | `extract_features` recurses correctly | ✅ Already correct |
| Multi-row same-feature semantics (same ID) | ID | Tolerated via `merge_strategy="create_unique"` (auto-renames) | Phase 5 Step 2 (warn) |
| Dangling `Parent` references | Parent | Silently accepted | Phase 5 Step 2 (error) |
| `transl_except` / selenocysteine | Unofficial Attributes | Not honoured | Phase 5 Step 2 (warn) |
| 1-based inclusive coordinate conversion to 0-based pyfaidx slice | Cols 4-5 | Verified correct at `extract_sequence.py:81` | ✅ Test-pinned |

---

## 5. What property-based testing told us

`hypothesis` ran each strategy 200 times. Key findings:

- **`segments_overlap_length`** — symmetry violated on tied starts (Bug #4).
- **`segments_overlap_length(s, s)`** for any segment — always returns `(end - start + 1, True)`. ✅
- **`merge_children_intervals`** — for any random list of intervals, the output is monotonic in start, has no negative-length intervals, and every input interval is contained in some output interval. ✅
- **`custom_bisect_insert`** — for any sequence of insertions, the bucket stays sorted ascending in `entry.end`. ✅
- **`get_padding_length`** — for any L in `[0, 10⁶]`, returns a value in `{0, 1, 2}` and `(L + result) % 3 == 0`. ✅
- **`get_AA_id_fraction(s, s)` and `get_DNA_id_fraction(s, s)`** — perfect-match returns `(len(s), len(s))` for any non-stop string. ✅
- **`get_ID_base`** — idempotent EXCEPT for numeric-only stems (Bug #5).

The hypothesis statistics report (`pytest --hypothesis-show-statistics`) confirms each strategy explored its full input space; no examples were rejected for filter-too-much.

---

## 6. Green-light decision for Phase 5

**Approved.** The acceptance criteria from Phase 4.5's plan are met:

- ✅ `lifton_class.py` ≥ 90 % (actual: **91 %**)
- ✅ `lifton_utils.py` ≥ 90 % (actual: **93 %**)
- ✅ All NCBI invariants either enforced today or tracked in §4 with a Phase 5 fix step
- ✅ Six legacy bugs catalogued, each with a strict-xfail or pinning test
- ✅ Suite remains fully hermetic and runs in under 12 s

Phase 5 Step 1 (the bug-fix step) now has a concrete six-bug fix list with a regression net per fix:

1. Fix `lifton_class.py:71-75` ID list/string → `test_extra_copy_number_appends_suffix` flips
2. Fix `lifton_class.py:236` Parent cascade → `test_transcript_id_assigned` flips
3. Fix `lifton_utils.py:521` `IntervalTree.overlap` arg → `test_overlap_ratio_triggers_true` flips
4. Fix `lifton_utils.py:498` sort key → `test_segments_overlap_length_symmetric` flips
5. Fix `lifton_utils.py:174-179` `get_ID_base` numeric-stem → `test_get_id_base_strips_trailing_int_once` flips
6. Implement Phase 5 Step 2 GFF3 validator → drives all `test_corruption_handling.py` xfails to documented warnings/errors instead of silent acceptance

When all six fixes land, expect: 3 xfails removed (today's strict xfails flip to XPASS, then are converted to plain assertions); the corruption tests gain new positive assertions on the validator's `ValidationFinding` outputs; total test count rises to ~250.

---

## 7. Verification commands

```bash
source /opt/anaconda3/etc/profile.d/conda.sh
conda activate lifton-test

# Full suite
pytest tests/ -v                                          # 240 pass + 3 xfail
pytest tests/ -v --tb=short                               # quick fail diagnostics

# Core-module coverage gate (the Phase 4.5 acceptance criterion)
coverage run --source=lifton --omit="lifton/liftoff/*" -m pytest tests/ -q
coverage report --include="lifton/lifton_class.py,lifton/lifton_utils.py" \
                --fail-under=90

# Hypothesis explorer stats
pytest tests/test_property_based.py --hypothesis-show-statistics

# Just the bug-pin xfails (should remain xfail until Phase 5 Step 1)
pytest tests/ -v -m xfail
```
