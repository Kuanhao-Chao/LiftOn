# Refreshing the public T2T-CHM13 annotation (2026-09-16)

Steven asked whether the CHM13 annotation posted on the LiftOn page should be
updated. It should. Two reasons, the second stronger than the first.

## What was posted

`ftp://ftp.ccb.jhu.edu/pub/data/LiftOn/JHU_LiftOn_v1.0_chm13v2.0.gff3`,
25 April 2024, LiftOn v1.0 — twelve releases old.

**It carried a known defect.** 1,094 of its 3,362 extra gene copies were emitted
as a bare `gene` line with no transcript, exon or CDS: the `-copies` identifier
bug fixed in v1.0.12.

**It was also invalid GFF3 and could not be loaded.** 29,502 duplicated exon IDs
across 447,932 rows. Every gffutils strategy fails on it with
`UNIQUE constraint failed: features.id`, so the per-transcript comparison in this
work could not be run against it at all, and any user pulling it into a standard
toolchain hits the same wall. `gff3-validate` reports 154 errors.

## What replaces it

`JHU_LiftOn_v1.0.12_chm13v2.0.gff3` — md5 `8ee2604760b622e1169f5f8fbc524181`,
1,237,312,969 bytes, also published in place as `human_refseq/lifton.gff3`.

Lifted with **LiftOn v1.0.12** from **GCF_000001405.40-RS_2025_08** (GRCh38.p14,
primary assembly, rRNA removed) onto T2T-CHM13v2.0, using the command the docs
publish:

```
lifton -g <reference> -o lifton.gff3 -copies chm13v2.0.fa GCF_000001405.40_GRCh38.p14_genomic.fna
```

| | posted v1.0 | **new v1.0.12** |
|---|---:|---:|
| rows | 4,037,192 | **4,047,289** |
| mRNA | 130,780 | **131,809** |
| exon | 2,113,964 | **2,119,501** |
| CDS | 1,680,900 | **1,688,385** |
| duplicated non-CDS IDs | 29,502 (447,932 rows) | **0** |
| `gff3-validate` errors | 154 | **0** |
| childless gene copies from LiftOn | 1,094 | **0** |
| coding transcript recall | not measurable | **0.99754** |
| mean protein identity | not measurable | **0.99789** |

"Not measurable" is literal: the old file cannot be loaded, which is itself part
of the answer to the question.

162 childless gene copies remain. All 162 are faithful passthroughs — Liftoff
itself emitted them without children — and none is attributable to LiftOn. The
run logs zero `Skipping … reference transcript not found` warnings.

## Two corrections made along the way

**The posted file used a primary-assembly-only reference.** The first run used
the full GRCh38 reference, which adds ~8,400 genes on alt/fix contigs. Those are
duplicate copies of primary genes and CHM13 has no alt contigs, so they compete
for a single locus: `gene-AATF-2` (chr17_KI270857v1_alt) displaced `gene-AATF` at
identical coordinates. Read literally that looked like "v1.0.12 loses 1,660
genes" — the opposite of the truth. Both references are now primary-only. The
full-reference run is kept at `lifton_chm13_2026/a_fullref/` as the evidence.

**Even primary-only, the reference is not identical to 2024's.** RS_2023_03 has
58,647 gene-like features against the 2024 run's 59,115, and 576 of that run's
genes do not exist in RS_2023_03. The posted file came from an older release, so
old-vs-new is *near*, not exactly, like-for-like.

## One open observation

Extra gene copies fell from 4,215 to 1,343, concentrated in small repeat
families (`MIR663A` 52→1, `LOC100419985` 49, `MIR10396A` 42, `DUX4L6` 25).

Traced `MIR663A`: **Liftoff's own `liftoff.gff3` contains 2 copies and LiftOn
reproduces it row for row.** The 2024 file placed 53 copies of this 92 bp miRNA
across chr13 satellite regions. So the change originates in the alignment /
copy-mapping stage, not in LiftOn's merge or ORF logic. Plausible drivers, not
separated: the minimap2 version (2.28-r1209) and the reference release. Fewer
repeat-region copies is plausibly higher precision, but it is a real change in
what the file contains and is recorded here rather than assumed benign.

## Still derived from v1.0 — follow-up

`human_refseq/UCSC_genome_browser/lifton.bb`, `mutations/`, `ref_chm13_cmp/` and
`visualization/` were not regenerated; the BigBed needs UCSC tooling that could
not be validated here. The CHM13 tutorial page says so explicitly.

## Host constraint worth recording

`vm.overcommit_memory=2`, `overcommit_ratio=97`, no swap → 977 GB CommitLimit.
`fork()` reserves the parent's address space per child with no copy-on-write
credit, so a whole-human parent (~35 GB) forking 32 workers reserves ~1.1 TB and
fails with ENOMEM **while ~930 GB of RAM is physically free**. Cap both fork
sites: `LIFTON_PARALLEL_LIFT_WORKERS` and `LIFTON_RESCUE_ISOFORM_WORKERS`. `-t`
may stay high — Steps 7/8 use threads, which cost no reservation.

**Robustness gap this exposed (not yet fixed):** `parallel_lift` already has a
serial fallback (`workers <= 1 → return False`), but an ENOMEM raised by `Pool()`
construction is not caught, so it aborts the whole genome instead of falling
back. A `try/except OSError` → serial would make whole-genome runs robust on
strict-overcommit hosts.

## Provenance

Runs, logs, per-arm manifests and the results JSON:
`/ccb/salz3/kh.chao/lifton_chm13_2026/`. The replaced files are backed up with
verified md5s at `ftp_backup_20260915/`, and the pre-change site at
`/ccb/salz3/kh.chao/lifton_docs_live_backup_20260916.tar.gz`.
