# LiftOn — Phase 16 Biological Benchmark Harness

End-to-end validation of the optimised LiftOn pipeline against the
five published JHU CCB same-species liftover datasets:

| ID | Species | Approx. RAM needed | Approx. download |
|---|---|---|---|
| `human` | Human (GRCh38 → T2T-CHM13) | 16–24 GB | ~6 GB |
| `mouse` | Mouse (GRCm39 → NOD/SCID) | 16–24 GB | ~5 GB |
| `bee` | Honey bee (HAv3.1 → ASM1932182v1) | 4–8 GB | ~0.5 GB |
| `arabidopsis` | Thale cress (TAIR10 → ASM2311539v1) | 2–4 GB | ~0.3 GB |
| `rice` | Rice (IRGSP → ASM3414082v1) | 4–8 GB | ~0.5 GB |

Sizes assume the Phase 15 lazy `pyfaidx` + streaming-Popen + bounded
heap configuration is active.

## What it does

`run_benchmarks.py` is a Python driver (no Snakemake / Nextflow
dependency) that:

1. **Downloads** each dataset's reference + target FASTA + reference
   GFF and, where available, a published comparison GFF from
   <ftp://ftp.ccb.jhu.edu/pub/data/LiftOn/>. Resumable; idempotent
   (skips files already on disk that meet a minimum-size sentinel).
2. **Lifts** with the Phase 11/15 fastest configuration:
   ```
   lifton --stream --inmemory-liftoff --locus-pipeline -t 8 \
          --native -copies -g <ref.gff> -o <out.gff3> \
          <target.fa> <ref.fa>
   ```
3. **Profiles** every invocation under `/usr/bin/time -v` (Linux) or
   `/usr/bin/time -l` (macOS) — captures peak RSS + wall-clock +
   user/sys CPU.
4. **Evaluates** an isolated copy of the newly lifted candidate by
   re-invoking LiftOn with `-E` against the reference annotation and
   reference/target sequences. The downloaded comparison GFF is preserved
   as an input artifact but is not used as truth. The harness parses LiftOn's
   `lifton_output/{score.txt,eval.txt,stats/...}` for mapped / lost /
   extra-copy / mean-identity metrics.
5. **Reports** as a stdout summary table plus a `summary_<UTC>.json`
   roll-up under `benchmarks/results/`.

## Files

```
benchmarks/
├── README.md                ← this file
├── datasets.json            ← dataset registry (URLs + flags)
├── run_benchmarks.py        ← Python driver (the work)
├── run_benchmarks.sh        ← Bash wrapper (env + log tee)
├── data/                    ← inputs (created on first run)
└── results/                 ← outputs + logs + roll-up JSON
```

## Prerequisites on the host

The benchmark uses LiftOn's optional `--native` path which depends
on `mappy` (PyO3 binding to minimap2). For the full pipeline you
also need:

* `lifton` (the CLI; `pip install -e .` from the repo root)
* `minimap2` on `PATH`
* `miniprot` on `PATH`
* `parasail-python`, `pysam`, `pyfaidx`, `gffutils` (already in
  `lifton.yml`)
* No extra Python deps — the registry is plain JSON
  (`benchmarks/datasets.json`).

The harness's first action is a runtime probe that prints whether
each binary was found.

## Quick start

```bash
# From the repo root, after installing LiftOn:
./benchmarks/run_benchmarks.sh --datasets bee
```

Run all 5 datasets:

```bash
./benchmarks/run_benchmarks.sh
```

Just download (e.g. on a Slurm login node):

```bash
./benchmarks/run_benchmarks.sh --download-only
```

Lift only, skip the evaluation pass:

```bash
./benchmarks/run_benchmarks.sh --no-evaluation
```

Re-run a dataset whose `.lifton.done` flag exists:

```bash
./benchmarks/run_benchmarks.sh --datasets human --force
```

The sentinel is a versioned JSON cache record, not a presence-only flag. It
stores the exact LiftOn argv, SHA-256 fingerprints for all three inputs, LiftOn
source and executable fingerprints, required dependency versions/metadata, the
output GFF3 hash, and the measured nonzero profile. Any mismatch or legacy
sentinel is stale and reruns the lift.

## Output layout

For dataset `<id>` the harness creates:

```
benchmarks/data/<id>/
    <reference.fna>
    <target.fa>
    <reference.gff>
    <published-comparison.gff3>  (optional; not used as evaluation truth)

benchmarks/results/<id>/
    lifton.gff3            ← the lifted annotation
    lifton_output/         ← LiftOn's own output dir
        score.txt
        chain.txt          (when --write_chains is set)
        stats/
            unmapped_features.txt
            extra_copy_features.txt
            mapped_feature.txt
            mapped_transcript.txt
    logs/
        lift.stdout.log
        lift.stderr.log
        lift.time.log      ← the parsed /usr/bin/time output
        evaluation.{stdout,stderr,time}.log
    evaluation/
        candidate.gff3      ← isolated copy evaluated with -E
        lifton_output/
            eval.txt
            run_manifest.json
    .lifton.done           ← idempotency sentinel
```

The roll-up JSON is `benchmarks/results/summary_<UTC>.json` and
contains every dataset's profile + parsed eval summary in one
machine-readable file.

## Publishing a paired release report

Release publication is stricter than ad-hoc aggregation. Run complete
`paired-subset`, `paired-full`, and `paired-e2e` controller stages with the
same source/toolchain provenance and a positive even repetition count (four is
the normal final-run setting). Canary stages are diagnostic and cannot produce
a release PASS.

Final roots must use the approved shared-host protocol exactly: 8 threads per
cell, 4 active cells, at most 2 full/E2E cells, 32 total worker threads, load-1
limit 32, minimum 256 GiB available memory, 15-second launch staggering, and a
30-second controller poll. Canary or custom-policy runs remain diagnostic
because their concurrency, memory envelope, and host-load exposure are not
comparable to the reviewed release campaign.

Generate and review a `campaign.json` that declares every canonical ID and
repetition:

```bash
LIFTON_PY=/path/to/release-environment/bin/python
"$LIFTON_PY" -c 'import json; from benchmarks.compare.release_report import \
canonical_campaign_spec; print(json.dumps(canonical_campaign_spec(), indent=2))' \
> campaign.json
```

The reporter independently re-derives the canonical 34 subset, 17 full, and
five E2E IDs. The reviewed lists must match those canonical sets in order and
must also exactly match each controller `plan.json`; a non-canary run created
with partial `--ids` remains diagnostic. Publish only from the three completed
controller roots:

```bash
"$LIFTON_PY" -m benchmarks.compare.release_report \
  --runs-root benchmarks/compare/_runs/<paired-subset-run> \
  --runs-root benchmarks/compare/_runs/<paired-full-run> \
  --runs-root benchmarks/compare/_runs/<paired-e2e-run> \
  --campaign-spec campaign.json \
  --candidate-sha <40-character-SHA> \
  --reference-sha <40-character-SHA> \
  --output-dir benchmarks/compare/release-publication
```

The reporter recomputes plan/cell fingerprints and checks the frozen source,
registry, mode, input, executable, and policy contracts. It live-rehashes every
successful pair result, arm manifest, GFF3, validation report, evaluator TSV,
and controller-recorded tooling source. Correctness gates require positive
coding recovery and protein-identity coverage, nonempty common-set PI, stable
subset/full baseline floors, and minimum E2E completeness/identity. Publication
uses an atomic directory replacement; a failed rerun removes the canonical
path and quarantines any prior PASS rather than leaving stale approval. For
partial exploration, pass `--diagnostic` without `--campaign-spec`; diagnostic
reports are explicitly labeled and can never claim a release PASS.

## Slurm / SGE submission template

```bash
#!/bin/bash
#SBATCH --job-name=lifton-bench
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=lifton-bench-%j.out

cd /path/to/LiftOn
./benchmarks/run_benchmarks.sh
```

For a single quick smoke-test on a login node (uses ~5 GB RAM and
takes a few minutes):

```bash
./benchmarks/run_benchmarks.sh --datasets bee --no-evaluation
```

## Submitting results back

After a run finishes, paste the contents of
`benchmarks/results/summary_<UTC>.json` plus the printed summary
table back into the Phase 17 chat. The JSON has every metric needed
for the manuscript table.
