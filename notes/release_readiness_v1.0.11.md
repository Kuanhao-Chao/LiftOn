# Release readiness — v1.0.11

Written 2026-08-17. v1.0.11 was published on 2026-08-01; this is a verification
pass over the released artifact, not a preparation for a new release. Nothing
here moves a tag or uploads to PyPI.

## The four release surfaces

Checking a tag alone is what let v1.0.10 sit unreleased for a month, so all
four are recorded:

| Surface | State |
|---|---|
| Git tag | `v1.0.11` → `c623f0b` |
| GitHub Release | "LiftOn v1.0.11", published 2026-08-01T09:20:07Z |
| PyPI | `lifton 1.0.11`, wheel and sdist, uploaded 2026-08-01T09:21 |
| `main` ancestry | reachable from `origin/main` and `origin/devel` |

`lifton/__init__.py` reads `v1.0.11`; `CITATION.cff`, `docs/source/conf.py`,
`README.md` and `docs/source/index.rst` all agree. `CHANGELOG.md` and
`docs/source/content/changelog.rst` both carry a v1.0.11 section.

## What was verified

| Gate | Result |
|---|---|
| `Run tests` on `main` | **green on 3.10, 3.11 and 3.12** |
| `Build and Deploy Sphinx` on `main` | **green**; khchao.com/LiftOn rebuilt and serves v1.0.11 |
| chr22 example, real minimap2 2.28 + miniprot | exit 0; gene **878/890 (98.7 %)**, pseudogene **358/370 (96.8 %)** |
| `gff3-validate` on that output | **0 errors**, 59 warnings, 77,788 data rows |
| Released wheel in a clean venv | `pip install --no-cache-dir lifton==1.0.11` succeeds; reports `1.0.11`; `gff3-validate` present |

The chr22 counts and the validator totals reproduce the v1.0.10 readiness note
exactly, which is the strongest available statement that the released build
behaves as recorded.

`--no-cache-dir` is not decoration. A cached wheel hides a broken sdist, which
is how `cigar` shipped broken for years while every local install looked fine.

## Defects this pass found and fixed

Three were in the released repository's automation, two in the benchmark
harness. None was a defect in the lifted output.

1. **Python 3.10 was never actually tested.** `tests/test_packaging_metadata.py`
   imported `tomllib` unconditionally, which is standard-library only from
   3.11, so the 3.10 job died at collection while `setup.py` continued to
   declare `python_requires>=3.10`. The import now falls back to `tomli`,
   which the `[test]` extra installs below 3.11.
2. **The documentation deploy had been failing since 2026-08-10**, so the site
   served a stale build for a week. Two dependabot bumps were responsible, and
   the first masked the second: `docutils==0.23` (PR #68) against
   `Sphinx==9.1.0` and `sphinx_rtd_theme==3.1.0`, both of which require
   `<0.23`; and then `sphinx-panels==0.6.0` (PR #66), which declares
   `sphinx<5`. Both are pinned back to the versions the last successful deploy
   used, and dependabot is told to hold them.
3. **Two test modules could not be collected without matplotlib**, which the
   figure generators import at module scope. It is now in the `[test]` extra.
   LiftOn itself never imports it.
4. **The benchmark scheduler admitted seven whole-genome cells against
   `max_active: 2`.** A launched worker counted as neither active nor pending
   while its tmux session was not yet listable and its identity not yet
   recorded. Two repetitions of the biology study were lost to
   `Cannot allocate memory` because of it. Admission and orphan detection now
   share one liveness rule and grace window.
5. **One malformed row discarded an entire annotation.** The target-truth
   parser raised on a structurally unusable record, which threw away four
   fully executed paired cells because LiftOn v1.0.8 emits a single
   inverted-coordinate mRNA on the dog-to-cat transfer. Such a row was never
   scored under the old behaviour either — the difference is that the parser
   used to take the other 90,000 transcripts down with it. It is now recorded
   and skipped, and the count is reported as evidence.

## Open items

- **`c8daa44` is released nowhere.** The trans-spliced-copies fix is on `main`
  and `devel` but is not in v1.0.11, and neither changelog mentions it. It
  needs its own version entry whenever it ships. It repairs the one invalid
  candidate output in the sealed qualification campaign.
- **`sphinx-panels` is abandoned** (last release 2021) and the documentation
  uses its `dropdown` directive in fifteen places. Migrating to
  `sphinx-design` is the real fix; until then the pin is held.
- **Three dependabot security alerts** on the default branch (1 high, 1
  moderate, 1 low), untouched by this pass.
- The benchmark evidence for the v1.0.11 technical report is a separate
  workstream; see `benchmarks/manifests/lifton_v111_biology_study.json`.
