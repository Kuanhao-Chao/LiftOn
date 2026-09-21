# Release readiness — v1.0.14 (cycle 3)

Written 2026-09-21 on branch `v1014-integration`. Cycles 1 and 2 landed 23
commits; this cycle adds three. Two of them are correctness fixes for defects
that shipped in every release; the third is a measured speed win on the
windowed aligner.

Nothing here has been pushed, tagged, released or published.

## What changed

| | what | class |
|---|---|---|
| P1 | `lifton_add_trans_exon_cds` asks for **level-1** exons, so `--threads 1` matches `--threads N` | output-corrective |
| P2 | no transcript emits overlapping exons: miniprot's redundant `stop_codon` is no longer ingested as a second exon, and `update_cds_list` reconciles a rebuilt exon against the one it ran into | output-corrective |
| P3 | `windowed_align._unique_anchors` indexes the reference only over the query's k-mers | byte-neutral |

Details and the reasoning behind each are in
`notes/threading_exon_divergence_2026-09.md`,
`notes/overlapping_exons_2026-09.md` and
`notes/windowed_anchor_construction_2026-09.md`.

## What was wrong the first time

Two claims in this cycle's own working notes had to be withdrawn after
measurement, and both are worth keeping visible.

**The P1 blast radius was counted in the wrong database.** The audit located
the divergence correctly and then counted how often it could fire by scanning
the *reference* annotation (1,915 human loci, 46 on chr22). The query reads
Liftoff's *output*, which contains none of that shape, because Liftoff does not
lift the nested miRNA. The real trigger is a different one entirely — RefSeq's
organellar convention — and it appears on rice (17) and arabidopsis (7), not on
human at all. A whole-genome human → CHM13 A/B was run on the strength of the
wrong count and correctly showed no change.

**The P2 gate was specified too strictly.** "Protein identity unchanged on
every already-valid transcript" failed on the first run — two rice transcripts
moved. Both had moved *up*, because repairing the terminal CDS let the miniprot
candidate be scored on its real sequence. The invariant that matters is that no
already-valid transcript gets worse, and the gate was corrected to that rather
than the result being explained away.

## Verification

Every A/B arm ran in a detached tmux session from a build pinned to an explicit
worktree, and asserted on load which `lifton/__init__.py` it had imported.
Paired arms shared one cached `-L`/`-M` so the build is the only difference.
Where an arm was started before a pinned worktree existed, it was re-run from
the pinned build and the two outputs compared byte-for-byte (identical).

### P1 — does `--threads 1` equal `--threads N`?

| | before: `-t 1` vs `-t 8` | after: `-t 1` vs `-t 8` |
|---|---|---|
| rice | **DIFFER** — 187,083,978 vs 187,068,200 bytes | identical |
| human → CHM13 | identical (1,237,342,491 both) | identical |

On rice the pre-fix **serial** arm is the only one of the four that differs:
`before -t 8`, `after -t 1` and `after -t 8` are byte-identical to each other.
Seventeen transcripts lose 35 exon rows and 35 CDS rows; all 17 carried
duplicate-coordinate exons before and none do after. `-t 1` is the default, so
the default path was the wrong arm.

Human → CHM13 is inert, which the Liftoff-output scan predicted (0 features
carrying the shape) and three measurements confirm: pre-fix `-t 1` equals
pre-fix `-t 8`; the pinned P1 build at `-t 8` is byte-identical to the pinned
pre-P1 build at `-t 8`; and the P1 build at `-t 1` closes the pair.

### P2 — overlapping exons, five whole genomes

Both arms pinned, one shared cached `-L`/`-M` per pair.

| | human → zebrafish | drosophila | CHM13 | rice | bee |
|---|---:|---:|---:|---:|---:|
| overlapping-exon transcripts | 96 → **0** | 47 → **0** | 27 → **0** | 13 → **0** | 13 → **0** |
| overlapping-CDS transcripts | 95 → **0** | 16 → **0** | 13 → **0** | 4 → **0** | 8 → **0** |
| genes / transcripts | unchanged | unchanged | unchanged | unchanged | unchanged |
| already-valid transcripts worse | 0 | 0 | 0 | 0 | 0 |
| already-valid transcripts better | 4 | 28 | 4 | 2 | 3 |
| pairs crossing a strand/seqid | 0 | 0 | 0 | 0 | 0 |
| "spans 2 exons" warnings | — | 2,632 → 0 | 8,546 → 0 | 4,468 → 0 | 3,614 → 0 |

(rice's 13 is what remains after P1 removed the 17 it was responsible for.)

### P3 — is the aligner change output-safe?

Whole-genome dog → cat, both arms pinned to frozen worktrees differing only by
this change: **523,466,820 bytes, byte-identical.**

### Suite

2,399 passed, 2 skipped, 0 failed (2,358 at the start of the cycle). 24-cell
matrix green with no golden edit. Fatal flake8 clean.

### P5 — is the drop ledger visible on real data?

The counter reaches `run_manifest.json` on every real run, with all seven
classes recorded including the new `hierarchy_depth_exceeded`, and the
end-of-run summary correctly stays silent when nothing was dropped.

It has still not been seen firing outside a test, because on these corpora
nothing is dropped. The 550 `Skipping … was not found` lines the plan expected
from rice were the `-copies` resolution bug, fixed in v1.0.12. That is the
right answer for these inputs, not a gap in the instrument — but it does mean
the classes are exercised only by unit tests.

## Known, not fixed

* **`nad5` in rice** — `rna-OrsajM_p05` is written on `CP132246.1` while its
  gene `gene-OrsajM_p05` stays on `CP132245.1`. A gene and its transcript on
  different sequences is invalid GFF3. It is a trans-spliced mitochondrial
  model (`exception=trans-splicing`), present identically before and after this
  cycle's changes. Surfaced by the new validator check; not caused by it and
  not fixed by it.
* **Two disjoint coding blocks under one overlapping exon pair** — an exon
  holds one CDS, so `reconcile_overlapping_exons` refuses rather than inventing
  coding sequence. No such pair occurs in the five genomes measured.

## A process failure worth recording

Mid-cycle I rewrote `p2_ab/arm.sh` in place to add a fifth genome. Two
`human → zebrafish` arms were still running, and bash reads a script
incrementally from an open file descriptor: `open(path, 'w')` truncates and
rewrites the **same inode**, so every byte offset after the insertion shifted
and both shells resumed at a misaligned position, re-executing the tail of the
script and starting a second lift on top of a finished one.

The artifacts survived, and were verified rather than assumed:

* `out.gff3` was byte-unchanged from the copy taken the moment the problem was
  noticed, so the second run never reached the publish step;
* there is no `*.partial.gff3` and `run_manifest.json` records
  `status: success` for both arms — `OutputTransaction` publishes only on
  success;
* the `before` arm is **byte-identical to `c1_ab/on`**, an independent,
  earlier, complete run of the same configuration.

That third check also independently confirms the P1 analysis: `c1_ab/on`
predates P1 and `before` includes it, and they are the same file, which is what
the Liftoff-output scan predicted for a genome carrying none of the affected
shape.

The rule this cycle already had — *pin every A/B arm to an explicit build* —
did not cover the driver script itself. A running script is as much live state
as a running tree. The arm scripts are now read-only, and a new variant goes in
a new file.

**A second one, same family.** The human → CHM13 A/B began as a four-arm serial
driver. To parallelise it I wrote empty placeholder `out.gff3` files so the
driver's `[ -f out.gff3 ]` guard would skip the arms I was moving, and gave the
parallel arms an `rm -f out.gff3` so they would not inherit a placeholder. The
two interact: the parallel arm deleted the placeholder, the driver reached that
arm before the parallel one had published, saw no file, and ran it again — from
the **live** tree, which by then carried P2 and P3. Its `after_t1` output
therefore had 0 overlapping-exon transcripts where its `after_t8` sibling had
27, and the difference read at first like a second threading divergence.

Diagnosis came from the driver's own log (`[after_t1] build: .../src/...`,
where the pinned arms log a worktree path) and the output timestamp, two hours
after the parallel arm had finished.

**The claim was then rebuilt on arms that are actually pinned, and it holds.**
Commit `3a92400`'s human → CHM13 row rests on these three measurements, none of
which involve the discarded arm:

| | |
|---|---|
| pre-fix `-t 1` vs pre-fix `-t 8` | identical, 1,237,342,491 bytes (both pinned to the pre-P1 worktree) |
| P1 build `-t 8` vs pre-P1 build `-t 8` | byte-identical — the fix changes nothing here |
| P1 build `-t 1` vs P1 build `-t 8` | byte-identical (`p1_chm13_clean/after_t1` vs `p2_ab/chm13/before`) |

All four arms are the same bytes, which is what the Liftoff-output scan
predicted for a genome carrying none of the affected shape. Anyone reproducing
this from the artifacts on disk should use `p1_chm13_clean/after_t1`, not
`p1_ab/after_t1` — the latter is the discarded arm and is left in place only so
this note can point at it.

The general shape, for the third time this cycle: a guard is only a guard if
nothing else is allowed to change what it tests.

## The regenerated CHM13 annotation

`/ccb/salz3/kh.chao/lifton_chm13_regen2/` — 1,237,332,649 bytes, 1 h 06 m,
peak RSS 30.3 GiB. Same recipe as the cycle-2 regeneration so the two compare
directly. **Staged, not published.**

| | staged (cycle 2) | regenerated (cycle 3) |
|---|---:|---:|
| genes / transcripts | 42,689 / 184,596 | unchanged |
| transcripts with overlapping exons | 27 | **0** |
| transcripts with overlapping CDS | 13 | **0** |
| `"reference model is malformed"` warnings | 8,537 | **0** |
| `gff3-validate` | `False`, 40 errors | **`True`, 0 errors** |

The 40 errors in the previous file were *all* of the kind this cycle fixed (27
`exon_overlap`, 13 `cds_overlap`), so the regenerated annotation is completely
clean. The cycle-1 genetic-code fix still holds: all 13 mitochondrial CDS match
reference length exactly.

This is what lets issue #26 be answered as fixed on both halves rather than
half-fixed — with the caveat, stated in every draft reply, that the **posted**
annotation has not been replaced.
