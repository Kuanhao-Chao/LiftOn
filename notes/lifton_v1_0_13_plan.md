# LiftOn v1.0.13 and subsequent efficiency work

User-approved implementation plan, refreshed 2026-09-14 after v1.0.12 release.
Binding specification: the user's approved correctness-first release plan in this
session, captured below. Baseline is v1.0.12, peeled commit
`6c86d1bf8dc0d53744a4410dcb2d6d2deb17f5e3`. Worktree:
`/tmp/lifton-v1.0.13-build`, branch `next-v1.0.13`. Earlier working trees and
historical evidence are preserved. v1.0.11 comparisons remain historical.

## Global constraints

- Correctness precedes efficiency. Ordinary supported inputs retain the 24-way
  byte-identity contract. Explicitly corrected inputs require independently
  specified expected differences, not blanket baseline equality.
- Linux Python 3.10–3.12. Audit and pin dependencies; upgrade only for a demonstrated
  need, separately tested. No new runtime dependency without no-cache installation.
- Focused regression tests before behavior fixes. No broad vendored refactors.
  Full suite on frozen final source; focused checks during iteration.
- Long runs use tmux and frozen source. Eight threads/cell, at most two full
  genomes, 32 aggregate scheduling threads including native aligners, and at
  least 256 GiB host memory reserve. Timed experiments exclusive, three
  alternating paired replicates. Record RSS and PSS separately.
- Preserve prior results. Unique campaign directories, exact source/import path,
  commands, hashes of inputs/tools/dependencies/evaluator, verified resume only.
  Independent target truth is evaluation-only; reference PI alone is insufficient.
- Keep default CLI partial-output compatibility; qualification rejects
  partial_success even with exit code zero. Experimental copy rescue stays off.
- No publication, tagging, push, deployment, external messages, or destruction of
  historical work. Endpoint is a reviewed, locally qualified candidate.

## Execution status (2026-09-15)

| Step | Status | Evidence / next action |
|---|---|---|
| Release synchronization and dependency audit | Complete | Fresh origin fetch and release API confirm v1.0.12 at 6c86d1b; dependencies audited without upgrades. |
| Task 1: release evidence | Review fixes committed; independent re-review pending | b7a9033 + c0fbe4a; 138 evidence/evaluator/inventory tests pass. Required campaign/attempt history, guarded imports and annotated coding inventory implemented. |
| Task 2: controller integration | Implemented; native canary passed | 0620294; 189 controller/evidence tests pass. Frozen candidate stream/inmemory equals v1.0.12 control; deep reconciliation verifies sealed first-attempt evidence. |
| Task 3: sparse coding models | Implemented; independent review pending | 54df064; 31 focused tests, supplied-FASTA combinations, both strands, mixed inputs and 30 fast/native tests pass. |
| Task 4: GTF hierarchy | Review fix implemented; re-review pending | f0d63fa + 9fe9f7e; real panel and fast benchmark pass. Direct GTF backend fallback: 44 focused tests pass in isolated Python 3.11. |
| Task 5: translation | Initial phase implemented; genetic codes next | Eight independent phase tests, including native runs on both strands, pass; 160 translation/alignment regressions and all 24 byte configurations pass. |
| Task 6: reproducibility | Pending | Fresh copy variation and complete copy hierarchies remain explicit release gates. |
| Task 7: correctness release qualification | Pending | No v1.0.13 release-readiness claim yet. |
| Tasks 8–11: efficiency, optional copies and final review | Pending | Follow correctness checkpoint. |

Independent Task 1/4 reviews completed and identified the gaps listed above.
The review/implementation agents subsequently reached their service usage limit;
saved work is preserved and implementation continues locally. No independent
approval is claimed for tasks whose review did not finish. Task 4 preceded
Tasks 2/3 because its confirmed input bug was
independent of scheduling and sparse normalization. The execution workspace
retains the dependency metadata, draft hashes, task reports and review packages.

The frozen f0d63fa real-GTF check retained all 8 genes, 26 transcripts (17 coding,
9 noncoding), 116 exons and 69 CDS. Raw GTF and independently retained GFF3 control,
each in ordinary and stream mode, all finished with manifest status `success`
and identical output SHA256
`ba2596d2f351f42a5d2bceadad50f7f4f77265dd49ac103355b68f803f875ff0`.
Commands, source/input hashes and manifests are in
`/tmp/lifton-v1013-gtf-evidence/evidence.json`. The human-MANE fast benchmark gate
also passed; its single timing observation is not a performance improvement
claim. These are intermediate checks, not qualification of the final branch.

Test environment correction: the shared editable installation can import a
missing LiftOn submodule from another checkout even when the top-level package
comes from this worktree. Focused sparse-model tests disable that finder; final
qualification must use isolated installs and verify actual module origins.

Task 2's real controller campaign is
`/tmp/lifton-controller-qualification-runs/native-controller-20260915`.
Frozen candidate 0620294 and reference 6c86d1b retained all 26 transcripts,
including all 17 coding transcripts with protein and DNA identity 1.0 and no
lost IDs. Candidate stream/inmemory and reference output have the same
`ba2596d2...f875ff0` SHA256 as the GTF control above. The first attempt and deep
reconciliation passed. This single canary qualifies the controller path, not
whole-genome correctness or performance of later candidate commits.

Task 5 initial-phase validation uses explicitly constructed coding sequences,
including split codons, phase 1/2, both strands, repeated scoring and stop
extension. Extraction and scoring trim the initial phase once; downstream
frames retain junction-spanning codons. Ordinary phase-zero extraction retains
its existing padding convention. Alternative tables and mixed-code native
alignment remain unfinished and are required before Task 5 is complete.

Isolated Python 3.10.21, 3.11.15 and 3.12.14 dependency environments are installed,
with no editable finders and no broken requirements. Their resolved locks and
download reports are retained in the execution workspace. This is environment
preparation, not cross-version candidate qualification. Python 3.10 uses the
compatible NumPy 2.2.6, NetworkX 3.4.2 and Matplotlib 3.10.8; the other interpreters
retain the audited direct dependency versions. The initial Python 3.10 install
used pip 23.0.1's legacy interlap setup path and emitted a deprecation; repeat
clean package qualification with a current isolated build frontend.

## Task 1: Complete portable release evidence and version-neutral roles

Review and finish the four unfinished draft files from the old build tree:
release_validation.py, release_provenance.py, test_release_provenance.py,
test_release_evidence_boundaries.py. Legacy inputs must come from consistent
recorded manifests and verify hashes, independent of worktree location. Never
invent missing historical provenance. Evaluator indexes must be private; source
sidecars unchanged. Retain every ERROR identity, bounded nonerror examples and
exact counts. Fail closed on changed inputs/tools/dependencies/evaluator, stale,
malformed or unsealed reports, unsuccessful latest retries, absent/empty,
duplicate/unexpected/missing expected cells. Distinguish ID presence, CDS
recovery, scored recovery and unresolved evaluation. Generalize hardcoded
v1.0.11/v1.0.12 arms to candidate/reference roles with version and commit as
separate metadata, retaining historical readers. Add focused behavior tests and
register tooling in benchmark inventory. Controller adapter belongs to Task 2.

## Task 2: Integrate qualification into the existing controller

Add a narrow qualification stage to benchmarks/compare/build_controller.py,
reusing existing selection, resource policy, scheduler, watchdog, attempts and
retry records. Do not build another scheduler. Freeze role-based source and
inputs, explicit expected set and campaign configuration. Preserve historical
paired-schema readers. Test dispatch/resume/failure boundaries and resource
accounting with lightweight cells before launching long runs.

## Task 3: Normalize sparse coding references

At annotation intake normalize parentless CDS and gene-to-CDS into deterministic
gene/transcript/exon/CDS models. Preserve CDS IDs, coordinates, strand, phase and
attributes. Group by logical source ID, allocate collision-safe generated IDs,
emit versioned original-to-normalized map and manifest evidence. The same
normalized reference drives extraction, Liftoff and lookups. Resolve supplied
protein/transcript aliases via the map and reject ambiguity. Reject dangling
parents, mixed strand/seqid groups and malformed mixed hierarchies. Ordinary
models remain a no-op. Tests include multisegment/both-strand/partial inputs,
independent sequence and coordinate truth, collisions and native CLI execution.

## Task 4: Preserve GTF hierarchy and attributes

Use verified gffread conversion options to retain/infer genes, transcripts and
attributes in a run-private converted file with conversion provenance. Cover
explicit/inferred hierarchies, coding/noncoding, exon/CDS attributes, stop codons,
Ensembl biotypes and direct-GTF opt-out separately. Qualify the retained real
Ensembl panel (8 genes, 26 transcripts, 17 coding, 69 CDS) in ordinary and stream
modes with no processing_error. Ordinary GFF behavior remains unchanged.

## Task 5: Correct phase and translation-table semantics

Trim initial CDS phase once in transcript orientation, never each segment;
preserve partial meaning and recompute downstream frames. Honor explicit tables
through extraction/scoring/ORF/stop logic. Keep native miniprot code consistent;
run mixed-code queries in sequential groups with separate indexes/provenance.
Reject conflicting or unsupported codes clearly. Independent expected proteins
for both strands, split codons, phases 1/2, partial models and alternative codes.
Ordinary table-1 behavior remains byte-identical outside explicit corrections.

## Task 6: Completeness and reproducibility

Retain released exact-ID-first and gene-membership copy guards. Measure complete
copy hierarchies and target loci, not only reference IDs. Qualification rejects
partial success. Diagnose fresh -copies variation with repeated native runs and
shared -L/-M replay. Reduce ordering defects before targeted vendored changes;
stable ordering must preserve biological ranking. Add a regression for every
proved nondeterminism; record remaining limits without claiming determinism.

## Task 7: Qualify and assemble the correctness release candidate

Freeze code, run full suite including properties, all 24 native configurations,
copy/new-input checks, fatal flake8, make test-fast and benchmark-gate. Qualify
Linux Python 3.10/3.11/3.12 on the same candidate. Native controls plus real GTF,
independent gffread and constructed sequence/coordinate truth. Default and
-copies comparisons versus v1.0.12 for Drosophila, Arabidopsis, rice, dog-to-cat,
human-to-zebrafish; real stream/inmemory/both on the three small genomes.
Failure/cancellation/cache/atomic publication coverage. Build wheel and sdist,
install each without cache outside source and execute lifts with import-path and
output parity checks. Update both changelogs, CLI docs, limits, architecture and
exact evidence packet. Review candidate; no publication. Performance tasks follow
this correctness checkpoint rather than delaying it with speculative changes.

## Task 8: Bound isoform work and bulk materialization

Add --rescue-max-inflight N (positive, default twice effective workers), bound
jobs/results, ordered attachment, worker-private reopened FASTA. Record queue
high-water and submitted/completed metrics. Bound beneficial bulk DB reads.
Adversarial worker ordering/failure and many-isoform tests, byte parity, independent
process-tree memory measurement. Preserve all suppression, IDs and placement.

## Task 9: Parallelize detached placement scoring

Only expensive pure placement scoring may run in workers. Candidate ranking,
acceptance, suppression, copy allocation, IDs and serialization remain serial in
original order. No shared mutable coordinator/SQLite connections. Reuse inflight
bounds; test competing/colliding candidates and delayed/errors. Three paired
performance replicates: experimental targets 20% lower rescue wall and 25% lower
working memory; retain measurements even when targets are not met.

## Task 10: Truth-gated optional copy experiment

Prototype --coortholog-rescue requiring -copies and default off. Deterministic
free-locus competition under existing identity/coverage/overlap/collision rules.
Human-to-zebrafish development, Arabidopsis-to-rice and Drosophila-to-bee holdouts,
synthetic positives/negatives. Retain runtime only with independently supported
extra loci and no precision/regression loss or processing errors; otherwise
remove runtime and keep diagnostics. Prior negative experiments stay rejected.

## Task 11: Final review and evidence refresh

Review the complete branch and integrate reviewed fixes. Repeat qualification
only where post-checkpoint code changed or evidence is unresolved. Update exact
source and performance/evaluation evidence, compatibility and limitations. Keep
a clear boundary between qualified v1.0.13 and any subsequent experimental work.
