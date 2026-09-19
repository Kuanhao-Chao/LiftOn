# Qualification controller protocol

The `qualification` stage uses the existing build controller's tmux scheduler,
watchdog, terminal status, resource admission and retry handling. It does not
replace the historical paired-release stages or their readers.

Create a JSON configuration with `schema_version: 1`, an absolute `python`
interpreter, `candidate` and `reference` objects containing absolute `root` and
exact lowercase 40-character `sha`, and an explicit `cells` object. Each cell
has a unique path-safe key and these fields:

```json
{
  "inputs": {
    "ref_gff": "/absolute/inputs/reference.gff3",
    "ref_fa": "/absolute/inputs/reference.fa",
    "tgt_fa": "/absolute/inputs/target.fa"
  },
  "species": "human",
  "cross_species": false,
  "annotation_database": "RefSeq",
  "copies": false,
  "candidate_mode": "safe",
  "reference_mode": "safe",
  "full_job": true
}
```

Modes are `safe`, `stream`, `inmemory`, or `stream-inmemory`. Copy search applies
to both arms. Use separate cell IDs for default and copy comparisons. Modes and
all three input paths are explicit; relocation never consults a checkout-relative
data directory. Target annotations/truth are not passed to the LiftOn CLI.

```bash
python -m benchmarks.compare.build_controller start \
  --stage qualification --qualification-config /absolute/campaign.json \
  --run-id v1013-correctness-UNIQUE --dry-run
python -m benchmarks.compare.build_controller start --run-id v1013-correctness-UNIQUE
python -m benchmarks.compare.build_controller status v1013-correctness-UNIQUE --json
python -m benchmarks.compare.build_controller reconcile v1013-correctness-UNIQUE --deep --json
python -m benchmarks.compare.build_controller retry v1013-correctness-UNIQUE
```

The plan pins the selected cells, sources, interpreter/dependencies/native tools,
inputs and evaluator. The configuration file can move after planning because its
normalized content is embedded in the plan; changed supplied configuration is
rejected on resume. Source and content evidence are rechecked before execution,
after execution, at success reuse and during qualification reconciliation.

Each cell uses `qualification_attempts/0001`, `0002`, etc. A failed or cancelled
attempt remains in place. A retry starts fresh; the controller's successful
result must point to the current attempt and its sealed report, summary,
campaign and history. It cannot reuse an earlier successful result after a
newer attempt fails. Exit zero plus `partial_success` does not qualify.

The default is 8 LiftOn threads with a conservative 16-thread scheduler charge
per cell, including native/helper reserve. At most two whole-genome cells,
32 aggregate scheduling threads, and at least 256 GiB available memory are
enforced; smaller thread counts are permitted for canaries. Native miniprot is
explicitly bounded and BLAS/OpenMP helper threads are set to one. The profiling
wrapper receives the exact recorded environment without inherited LiftOn
overrides. Full paired cells receive twice the full-cell watchdog allowance.

This stage establishes correctness evidence. Its timings are descriptive;
performance claims still require exclusive alternating paired replicates and
independent RSS/PSS measurements. Final release qualification also needs the
separate independent biological truth checks in the release plan.
