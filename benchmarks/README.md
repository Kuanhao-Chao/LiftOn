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
   GFF + ground-truth GFF from
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
4. **Evaluates** by re-invoking LiftOn with `-E` against the
   pre-computed truth GFF, then parses LiftOn's own
   `lifton_output/{score.txt,eval.txt,stats/...}` for mapped / lost /
   extra-copy / mean-identity metrics.
5. **Reports** as a stdout summary table plus a `summary_<UTC>.json`
   roll-up under `benchmarks/results/`.

## Files

```
benchmarks/
├── README.md                ← this file
├── datasets.yaml            ← dataset registry (URLs + flags)
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

## Output layout

For dataset `<id>` the harness creates:

```
benchmarks/data/<id>/
    <reference.fna>
    <target.fa>
    <reference.gff>
    <truth.gff3>           (optional)

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
    .lifton.done           ← idempotency sentinel
```

The roll-up JSON is `benchmarks/results/summary_<UTC>.json` and
contains every dataset's profile + parsed eval summary in one
machine-readable file.

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
