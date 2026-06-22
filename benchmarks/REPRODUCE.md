# Reproducing the LiftOn v1.0.9 comparison benchmarks

This guide describes how to reproduce the published comparison between LiftOn
(devel / v1.0.9), the previous stable release (v1.0.8), and the two single-method
baselines it builds on, **Liftoff** (DNA-DNA) and **miniprot** (protein-DNA).

All scripts referenced below live under `benchmarks/compare/`. The numbers in the
committed `*.json` / `*.md` artifacts were produced with the deterministic
settings described in **Determinism** — reproduce those first, then run the
scripts.

## Code version

Check out the `devel` branch at the commit the v1.0.9 artifacts were generated
from:

```bash
git clone https://github.com/Kuanhao-Chao/LiftOn
cd LiftOn
git checkout devel        # v1.0.9 release line
git rev-parse HEAD        # expected: d11b30b (or the v1.0.9 release tag)
pip install -e .          # see CLAUDE.md / lifton.yml for the conda env
```

External tools `minimap2` and `miniprot` must be on `PATH`.

## Determinism (run all comparison cells with these)

The published comparison is **byte-deterministic**. Pin every run with:

- `PYTHONHASHSEED=0` — stabilizes any dict/set ordering.
- `-t1` (single thread) — `--threads N` is byte-identical to `-t1`, but a single
  thread removes any scheduling-related ambiguity when diffing outputs.
- **No `-copies`** — Liftoff's multi-copy (`-copies`) alignment is non-deterministic
  across fresh runs (tRNA repeat copies differ run to run), which would inject
  noise into a fresh A/B diff. Leave `-copies` off for reproducible comparisons.

```bash
export PYTHONHASHSEED=0
# pass -t1 and do NOT pass -copies to every lifton / Liftoff invocation
```

## Apples-to-apples baseline: rescue is OFF for the 4-way cells

As of v1.0.9 the **miniprot-only rescue pass is ON by default** (see
`--no-miniprot-rescue` in `lifton --help`). To keep the 4-way comparison
baselines apples-to-apples against v1.0.8 (which has no rescue), the devel /
v1.0.9 benchmark cells are run with **`--no-miniprot-rescue`**. The harness does
this automatically: `version_compare._build_argv` (reused by `fourway_compare.py`)
appends `--no-miniprot-rescue` to the devel / devel-legacy cells, so the committed
`fourway_results.json` stays comparable across versions.

The recall gain from the rescue pass is reported as a **separate arm**, not folded
into the 4-way devel column — see `miniprot_rescue_ab.py` below.

## Scripts → artifacts

| Script | Produces | What it does |
|---|---|---|
| `benchmarks/compare/fourway_compare.py` | `fourway_results.json` | The headline 4-way comparison (Liftoff vs miniprot vs LiftOn v1.0.8 vs LiftOn devel/v1.0.9) across the subset and full-genome cells, using a neutral re-scorer. Devel cells run with `--no-miniprot-rescue`. |
| `benchmarks/compare/version_compare.py` | `version_compare.results.json` | The v1.0.8-vs-devel version comparison (the layer `fourway_compare.py` builds on; `_build_argv` is the single place that pins `--no-miniprot-rescue` for devel cells). |
| `benchmarks/compare/joint_metrics.py` | apples-to-apples joint recall/identity metrics | Computes common-set protein identity, `covPI` (= recall x accuracy), recall@PI, and a paired sign-test — the metric that separates the per-transcript accuracy edge from the recall (completeness) axis. |
| `benchmarks/compare/enrich_joint_metrics.py` | joint-metric fields enriched into `fourway_results.json` | Post-processes the 4-way cells to add the `joint_metrics.py` numbers per cell. |
| `benchmarks/compare/miniprot_rescue_ab.py` | `miniprot_rescue_ab.json` / `.md` | The **separate** rescue-on vs rescue-off A/B arm: `--no-miniprot-rescue` (off) vs default (on). Reports n_added / n_lost / n_redundant / mean-PI-of-added / completeness delta per dataset (the strict gate the rescue passes 8/8). |
| `benchmarks/compare/fourway_report.py` | `fourway_report.json` / `.md` | Renders the 4-way results into a human-readable report. |
| `benchmarks/compare/subset_builder.py` | per-dataset subset inputs | Builds the chromosome/gene subsets and reference proteins used by the comparison cells. |
| `benchmarks/compare/evaluator.py` | per-tool scores | The neutral re-scorer (protein identity, completeness/recall) used by `fourway_compare.py`. |

## Typical run order

```bash
export PYTHONHASHSEED=0

# 1. Build subset inputs (if not already cached)
python benchmarks/compare/subset_builder.py

# 2. Run the 4-way comparison (devel cells pinned --no-miniprot-rescue, -t1, no -copies)
python benchmarks/compare/fourway_compare.py        # -> fourway_results.json

# 3. Enrich with the apples-to-apples joint metrics
python benchmarks/compare/enrich_joint_metrics.py   # uses joint_metrics.py

# 4. Render the report
python benchmarks/compare/fourway_report.py         # -> fourway_report.{json,md}

# 5. Separately, the rescue recall arm (rescue-on vs rescue-off)
python benchmarks/compare/miniprot_rescue_ab.py     # -> miniprot_rescue_ab.{json,md}
```

Run each script with `--help` for its exact flags and dataset selectors. The
committed `*.json` / `*.md` files under `benchmarks/compare/` are the reference
artifacts to diff against.
