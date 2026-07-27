# Repository Guidelines

## Project Structure & Module Organization

`lifton/` contains the Python package and CLI pipeline; `lifton/lifton.py` is the main entry point. Treat `lifton/liftoff/` and `lifton/gffbase/` as vendored subsystems and avoid broad refactors there without targeted regression coverage. Current unit, integration, property, and native-path tests live in `tests/`, with performance cases in `tests/perf/`. The singular `test/` directory holds legacy tests, shell examples, and FASTA/GFF fixtures. Documentation is under `docs/`, benchmark tooling under `benchmarks/`, and project artwork under `graphics/`.

## Build, Test, and Development Commands

- `python -m pip install -e ".[test]"` installs LiftOn plus pytest, Hypothesis, and coverage tools.
- `pytest tests/ -q` runs the complete Python suite.
- `make test-fast LIFTON_PY=python` runs the byte-identity matrix and integration tests. Override `LIFTON_PY` because the Makefile default is maintainer-specific.
- `make benchmark-gate LIFTON_PY=python` checks regressions and the fast benchmark before output- or performance-sensitive changes.
- `python -m build` creates wheel and source distributions.
- `lifton -g reference.gff3 -o out.gff3 target.fa reference.fa` runs the CLI; real runs require `minimap2` and `miniprot` on `PATH`.

## Coding Style & Naming Conventions

Use four-space indentation and established Python naming: `snake_case` for modules and functions, `PascalCase` for classes, and `UPPER_SNAKE_CASE` for constants. Keep changes focused and lines near the CI advisory limit of 127 characters. No automatic formatter is configured. Before submitting, run the blocking lint check: `flake8 . --select=E9,F63,F7,F82 --show-source`.

## Testing Guidelines

Name files and functions `test_*.py` and `test_*`; use `Test*` for classes and shared fixtures from `tests/conftest.py`. Add a focused regression test for each bug fix. Changes to alignment, GFF3 writing, IDs, or per-locus results must keep `tests/test_native_matrix.py` byte-identical across all 24 configurations. Coverage is available via `coverage run --source=lifton --omit="lifton/liftoff/*" -m pytest tests/ -q`; no repository-wide percentage threshold is enforced.

## Commit & Pull Request Guidelines

Follow the established `type(scope): imperative summary` pattern, for example `fix(gff3): preserve parent IDs (GH #32)`. PRs should explain what changed and why, link relevant issues, list validation commands/results, and call out output, compatibility, or performance effects. Include screenshots only for visual documentation or benchmark-figure changes.

## Security & Configuration Tips

Do not commit secrets, local environments, generated outputs, or large datasets. Before adding a runtime dependency, confirm it installs in a clean environment with `pip install --no-cache-dir`: LiftOn shipped for years depending on `cigar` 0.1.3, whose sdist cannot build at all on a machine without a cached wheel.
