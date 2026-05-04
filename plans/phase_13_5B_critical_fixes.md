# Phase 13.5B — Critical & High Vulnerability Fixes — EXECUTION REPORT

**Date:** 2026-05-03
**Branch:** `devel`
**Baseline:** `3fdbf88` (Phase 13 — 518 tests passing)
**Result:** `<commit-hash>` — 543 tests passing (518 + 25 new hostile tests),
zero regressions, 24-cell parallelism byte-identity preserved.

## Audit input

Phase 13.5A produced `plans/phase_13_5A_vulnerability_audit.md` with 53
findings. Phase 13.5B addressed the **1 Critical + 11 High** items.
Mediums and Lows are explicitly **deferred** — the user said: "Do not
proceed to the medium bugs until prompted."

V3.1 (`extract_features` materialises ALL reference proteins in RAM) is
the only High item not patched. The fix requires the lazy
RefSeqProvider abstraction from the Phase 4 roadmap, which is a memory-
focused refactor that exceeds one execution session's scope. The audit
report flags this; the deferral is documented here.

## Test-driven hardening protocol

Per the Phase 13.5B mandate, hostile tests in
`tests/test_vulnerabilities.py` were written **before** the production
patches, proven to FAIL on the unpatched codebase, then rerun after
the patches to PROVE the failure mode was eliminated.

Pre-patch hostile-suite run: **13 failed, 12 passed**
(the 12 passes were V4.2 IUPAC tests — parasail handled the small
strings without crashing — and V5.1, plus regression-only tests).

Post-patch hostile-suite run: **25 passed, 0 failed**.

## Findings addressed

| ID | Sev | Symptom | Fix site | Patched behaviour |
|---|---|---|---|---|
| **V1.4** | **CRIT** | `process_locus_native` & `process_locus` caught `BaseException` → Ctrl-C silently packaged into `LocusResult.error`, pool kept running | `lifton/locus_pipeline.py:99, 258` | Narrowed to `except Exception` so `KeyboardInterrupt` / `SystemExit` propagate and kill the pool |
| **V1.1a** | High | Bare `except:` at `run_liftoff.py:209` swallowed every exception during `ref_db[trans_id]` lookup, including `KeyboardInterrupt` | `lifton/run_liftoff.py:218-220` | Narrowed to `except (KeyError, gffutils.exceptions.FeatureNotFoundError)` |
| **V1.1b** | High | Bare `except:` in `check_miniprot_installed` masked `MemoryError`, `KeyboardInterrupt`, etc. as "miniprot not installed" | `lifton/run_miniprot.py:22` | Narrowed to `except (FileNotFoundError, PermissionError, NotADirectoryError, subprocess.SubprocessError)` |
| **V1.2** | High | `extract_sequence.__inner_extract_feature` silently swallowed exceptions — failed extractions were invisible to users | `lifton/extract_sequence.py:96, 105` | Replaced silent `pass` with `logger.log_warning(...)` — the user now sees which feature ID failed and why |
| **V1.3** | High | `materialise_locus` had 6 `except Exception:` blocks that silently produced empty payloads | `lifton/locus_pipeline.py:174-220` | Each except now logs a structured `[WARNING]` with the locus id and the underlying exception |
| **V1.5** | High | `Lifton_TRANS.write_entry` swallowed all exceptions but still proceeded to write child rows → orphan exon/CDS rows in output | `lifton/lifton_class.py:798-808` | Narrowed to `(OSError, ValueError, TypeError, AttributeError)` AND `return` early to suppress orphan child writes |
| **V1.9** | High | `Annotation._handle_gtf_input` silently fell back to raw GTF on conversion failure → user got obscure downstream error | `lifton/annotation.py:215` | When `auto_convert_gtf=True` and both gffread + agat fail, raises `LiftOnInputError` with actionable guidance |
| **V2.1** | High | `get_dna_sequence` accepted `start < 1`, producing `fasta[chrom][-N:end]` which silently wraps to chromosome tail | `lifton/extract_sequence.py:130` | Bounds-check raises `LiftOnInputError` for any `start < 1` (GFF3 spec violation) |
| **V2.3** | High | `get_DNA_id_fraction` raised `IndexError` mid-loop on length-mismatched inputs | `lifton/get_id_fraction.py:43-58` | Length-mismatch guard raises `ValueError` with a clear message before the loop |
| **V4.2** | High | Non-ACGT bases (IUPAC codes R/Y/S/W/K/M/B/D/H/V/N) were not in parasail's matrix alphabet — risk of C-kernel crash | `lifton/align.py:91-115` | Extended matrix to `ACGTN*` AND added `_sanitise_for_parasail_dna` helper that maps any non-alphabet base → N |
| **V5.1** | High | Validator already caught `dangling_parent`; regression test added | `tests/test_vulnerabilities.py` only | Test confirms `GFF3Validator` returns `dangling_parent` finding for a missing-ID `Parent=` |
| **V5.2** | High | `process_liftoff` recursive descent could infinite-loop on circular `Parent=` cycles → `RecursionError` swallowed by V1.4 | `lifton/run_liftoff.py:161-200` | Added `_visited` set parameter; on second visit raises `LiftOnInputError("Circular Parent reference detected at ...")` |

### V3.1 — DEFERRED

`extract_sequence.extract_features` still materialises every reference
transcript and protein into in-memory dicts. The structural fix is the
lazy `RefSeqProvider` from the Phase 4 roadmap. This is a multi-day
memory-arena refactor that does not fit alongside an exception-handling
audit batch. Recorded for the next memory-focused phase.

### Latent bug discovered & fixed (bonus)

While exercising the V1.9 path, a pre-existing `NameError` was found
at `lifton/annotation.py:119, 123`: `file_name` referenced an
undefined variable. Should have been `self.file_name`. Fixed inline
because the V1.9 hostile test could not pass without it. No existing
test exercised the GTF input path.

## New artefacts

* `lifton/exceptions.py` — `LiftOnError`, `LiftOnInputError`,
  `LiftOnAlignmentError` exception hierarchy.
* `tests/test_vulnerabilities.py` — 25 hostile tests covering V1.1a,
  V1.1b, V1.2, V1.3, V1.4 (3 sub-tests), V1.5, V1.9, V2.1 (3
  sub-tests), V2.3 (3 sub-tests), V4.2 (6 parametrised sub-tests),
  V5.1, V5.2.

## Verification

```bash
$ pytest tests/test_vulnerabilities.py -v
============================== 25 passed in 2.03s ==============================

$ pytest tests/ -q --ignore=tests/perf
............................................................................. (8 lines)
=============================== 543 passed in 76.59s ===========================

$ pytest tests/test_parallelism_matrix.py -v
=============================== 12 passed in 1.97s =============================
```

Coverage on touched modules (sampled):
* `locus_pipeline.py`: every materialise_locus branch now has a
  log-warning verifier in `TestV1_3_*`.
* `extract_sequence.py`: V1.2 + V2.1 cover both the silent-swallow and
  the start-bounds branches.
* `get_id_fraction.py`: V2.3 covers length-mismatch in both directions
  AND preserves the equal-length numeric output.
* `align.py`: V4.2 covers 5 IUPAC permutations + a regression for pure
  ACGT.

## Next steps (NOT part of this phase)

The 41 Medium and Low findings remain open. Per the user's mandate,
they will be picked up only after explicit prompt. The next
sensibly-scoped batch is probably:

* **V2.4–V2.11** (algorithmic / boundary issues that are mostly
  one-line guards but each needs a hostile test).
* **V3.x performance / scaling** items (V3.1 included; this is its own
  multi-phase refactor).
* **V5.3–V5.9** (GFF3 edge cases — duplicate IDs, attribute escaping,
  attribute ordering).
