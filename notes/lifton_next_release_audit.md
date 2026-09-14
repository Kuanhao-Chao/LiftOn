# LiftOn next-release audit and implementation record

Started 2026-09-14 against `devel` at
`e30101dc314d89a6a07e21038fd2e0d19d2d3b1e`. The released comparison baseline is
`v1.0.11`, `c623f0bddc5b5051a3670a3b0064c52cc61bb719`. This document distinguishes
observations from release qualification; an unfinished gate is not a pass.

## Objective and constraints

Prepare v1.0.12 with correct annotation behavior, reproducible evidence, and
measured improvements to rescue runtime and working memory. Correctness takes
precedence over resource gains. Qualify Linux on Python 3.10–3.12. Complete the
advertised flat-CDS and GTF input paths. Additional-copy recovery is an optional,
truth-gated experiment, not a default change. Publishing and contacting reporters
are separate release actions.

Long runs use detached tmux and immutable source snapshots. Qualification uses
8 threads per cell, at most two simultaneous whole-genome cells, an aggregate
32-thread scheduling budget including concurrent aligners, and a 256-GiB host
memory reserve. Timed comparisons run exclusively, three alternating paired
replicates. Retain old results rather than overwriting release evidence.

## How the current algorithm works

```mermaid
flowchart TD
  I[Reference annotation + reference and target FASTA] --> V[Scan, validate, build/reuse annotation DB]
  V --> E[Select roots; extract reference transcripts and proteins]
  E --> L[Liftoff: minimap2 DNA alignment and coordinate lifting]
  E --> M[miniprot: protein to target-genome alignment]
  L --> D[Feature DBs, reference-ID maps, target interval indexes]
  M --> D
  D --> S7[Step 7: DNA lift, chained CDS, native miniprot candidate; ORF rescue]
  S7 --> S8[Step 8: miniprot loci outside emitted genes]
  S8 --> R[Rescue A: span gate; rescue B: protein-coverage gate]
  R --> ISO[Attach co-located isoforms after placement is final]
  ISO --> W[Normalize hierarchy; allocate IDs; serialize and validate]
  W --> P[Atomic GFF3 publication + scores, statistics, run manifest]
```

The reference database defaults to gffutils/SQLite. Optional streaming and
in-memory paths use vendored gffbase/DuckDB. Sequence access uses pyfaidx.
Parasail performs alignment; Python code performs chaining, model assembly,
variant classification, and ORF rescue. Protein identity against the reference
is an optimization objective, not an independent biological truth label.

Step 7 threads coordinate copy allocation and cross-locus interval reads through
`Step7StateCoordinator`. Detached feature materialization prevents workers from
using another thread's SQLite connection. Results are committed in submission
order. Step 8 similarly evaluates candidates and rechecks acceptance serially.
Rescue placement updates the suppression tree and emitted-reference-gene set;
these decisions must remain ordered even when scoring is parallel. Isoform
scoring already uses forked workers; attachment and publication remain serial.

The output transaction stages and validates the GFF3 before replacing the final
path. Worker processes must not execute the parent's signal handler. FASTA file
offsets must not be shared across forked workers: workers reopen the files.
Input and annotation caches have separate provenance and locking requirements.

## Claude history: sources and evidence boundary

The available local LiftOn history includes sessions beginning 2026-07-02,
2026-08-16, and 2026-09-11, ending as recently as 2026-09-14. The two latest
substantial sessions are `e1194ba0-8508-473d-8be2-194ad4f7741b` and
`25545d1d-1c84-417b-b380-c5b5473fd1db`. Local project memories summarize earlier
iterations. Raw conversations are private working records and are not copied
into the repository. Git history and reproducible artifacts take precedence
over an assistant's claims. Missing/truncated history is not reconstructed as
fact.

| Work | Current evidence and interpretation |
|---|---|
| v1.0.10 alignment, hierarchy, IDs, and installation fixes | Historical audit in `lifton_correctness_audit.md`; `cigar` clean-install failure establishes why empty-cache packaging tests matter. |
| v1.0.11 duplicate-ID collision repair | Exact tag `c623f0b`; keep this release's sealed scientific evaluation unchanged. |
| Large-target native scheduling and diagnostics | `307abc6`; public 22-Gb indexing surrogate establishes large independent peaks, not the cause of the private issue #71 failure. |
| Coverage and isoform rescue | `79a212a`; 13-cell feature A/B reports support additions without losses under their measured protocol. |
| Default locus pipeline and parallel Liftoff | `d4c1f02`, `8816ef3`, `804f01e`; measured byte-equivalent paths. Timings must retain cache/thread/protocol context. |
| Terminal-stop completion | `d3fa4da`, `cff6ec0`, `8b91f3a`; 13-cell A/B restricts changes to three-base terminal extensions; reference stop convention and post-ORF timing are essential. |
| Feature clone and streaming validator | Already implemented before this audit. Older notes calling them unfinished are superseded by current code and tests. |
| Cross-locus replacement | Still opt-in: revised experiment loses 161 transcripts net. No default promotion. |
| Flat annotation fallback | Root selection accepts CDS, but extraction only visits children. The full input-to-output path remains unresolved at this baseline. |
| `--stream`/`--inmemory-liftoff` at scale | Historical under-recovery/hangs are recorded. Current synthetic parity does not settle whole-genome behavior. |
| Latest 17-cell release campaign | Outputs and rescored records exist locally; harness is initially untracked and Markdown summary stale. Requires provenance and stronger gates. |

The latest campaign uses `-t 16 -copies`, fresh native alignment, and reusable
reference DBs. It is not a no-option default or fully cold-cache campaign.
Direct set comparison of the 17 transcript tables found zero lost recovered
coding IDs and no duplicate reference rows. Zero identity regressions applies
only to the numerically scored common set. For chicken, 1,258 recovered models
in each arm have no protein identity and status `map_failed`; their recovery
cannot be described as independently confirmed accuracy.

Follow-up inspection resolved the chicken class: all 1,258 rows have
`n_cds_lifted=0` and `lifted_prot_len=0`; an inspected emitted model has exons
but no CDS and is explicitly tagged `mutation=no_protein`. These are mapped
reference IDs without recovered coding sequence, not unexplained failed
alignments. The revised comparison records this class separately and gates on
loss of previously scored coding models. A blank score without corresponding
CDS/sequence evidence remains unresolved.

The evaluator uses LiftOn's extraction/alignment code. It provides a consistent
yardstick across tools but can share bugs with the implementation. Independent
sequence extraction, target-coordinate concordance, and constructed truth must
supplement it. An annotation's own invalidity is context, never blanket
permission to introduce new output defects.

## Prioritized work and acceptance

| Milestone | Required result | Status |
|---|---|---|
| 1. History and architecture | Reconciled evidence ledger and current data/state flow | In progress: this report |
| 2. Release evidence | Fail-closed gates, explicit lost/unscored models, pinned snapshots and resumable provenance | In progress |
| 3. Inputs | Flat CDS and actual GTF work end to end; preserve source semantics | Pending |
| 4. Parity and reliability | Real fast-path equivalence, reduced regressions, injected-failure coverage | Pending |
| 5. Rescue efficiency | Bounded isoform work, bulk reads, ordered parallel placement scoring | Pending |
| 6. Independent verification | Sequence, coordinate and known-truth evidence with explicit exclusions | Pending |
| 7. Additional copies | Optional experiment retained only if independent gates pass | Pending |
| 8. Qualification | Full suite, Linux CI matrix, whole genomes, wheel/sdist execution, documentation | Pending |

Performance targets are 20% lower rescue wall time and 25% lower rescue working
memory on affected workloads. These are experimental targets, not promised
results. Report parent/process-tree RSS and PSS separately: summed RSS counts
fork-shared pages more than once. Overlapping manifest phases must not be summed
as end-to-end wall time. Small improvements must exceed paired-run variability.

No broad vendored refactors, process-based Step 7, concurrent gffutils builds,
or speculative rescue-threshold sweeps are planned. Prior negative experiments
stay rejected unless new evidence changes their premises. Final-model start/stop
quality will be measured independently of historical rescue-trigger labels.

## Validation log

- Planning baseline: 74 tests passed in 246.99 seconds: native matrix,
  integration pipeline, sequence extraction, streaming validator. This was a
  focused check, not a new full-suite qualification.
- Current-head GitHub tests were successful at planning time (run
  `34767497257`). New implementation commits require their own qualification.
- Release-shell fault tests: **8 failed / 1 passed before the repair; 9 passed
  after**. Injected failures cover lifting, validation, build, installation,
  wheel execution, serial/thread parity, wheel parity, and preservation of
  previous evidence. The script now fails on the first error, compares actual
  bytes, and runs the wheel outside the source tree without PYTHONPATH leakage.
- `make test` now includes the three previously ignored property-test modules.
- Release-comparison regressions: **14 failed / 2 passed before the repair**;
  the expanded comparison/provenance suite now has **26 passing tests**.
  Coverage includes ID loss despite net gain, nonfinite/missing identity,
  missing CDS versus scoring failure, per-model validity regressions, changed
  binaries/source/commands, altered outputs, and incomplete receipts. Legacy
  rescoring is written to a new campaign and cannot manufacture provenance.

## Release state

Not release-qualified. Milestones above must be resolved with evidence before
building the final release packet. No new tag, package publication, website
deployment, or reporter message is part of this implementation checkpoint.
