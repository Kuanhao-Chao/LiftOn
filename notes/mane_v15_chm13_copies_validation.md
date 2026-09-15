# MANE v1.5 → CHM13v2.0 with `-copies`: confirming the childless-gene fix

The `-copies` childless-gene bug was reported against MANE. Every gate up to this point
used the benchmark corpus, so this run reproduces the reported command on the reported
input and confirms the fix there.

## The command

Exactly as reported, plus `-t 16` — which the pinned 24-cell byte-identity matrix makes
neutral to the output, and which is what lets this finish in 21 minutes:

```
lifton -t 16 -g MANEv1.5.gff -chroms chrom_mapping.txt -copies -sc 0.9 \
       chm13v2.0.fa GCF_000001405.40_GRCh38.p14_genomic.fna
```

| Input | Provenance |
|---|---|
| `MANEv1.5.gff` | `MANE.GRCh38.v1.5.refseq_genomic.gff` from NCBI, sha256 `040f0d4056de2e9a416cd52bc20ff07ef403baadf5f5968faf188907893f6002` — 19,363 genes / 19,367 mRNA |
| reference genome | `benchmarks/data/human/GCF_000001405.40_GRCh38.p14_genomic.fna` (GRCh38.p14, `chr`-named) |
| target genome | `benchmarks/data/human/chm13v2.0.fa` |
| `chrom_mapping.txt` | 24 primary pairs, derived from the intersection of the two `.fai` files — covers 19,299 of 19,363 genes; the 64 excluded sit on GRCh38 alt/fix contigs, which CHM13 has no counterpart for |

Two arms, identical inputs: **off** = the pre-fix tree (detached worktree at `f1f16e0`)
pinned on `PYTHONPATH`; **on** = the fixed tree.

## Result

| | off (pre-fix) | on (fixed) |
|---|---:|---:|
| Exit / wall / peak RSS | 0 / 1,248 s / 30.1 GB | 0 / 1,294 s / 30.1 GB |
| Rows emitted | 262,927 | 266,278 |
| Extra gene copies | 509 | 509 |
| **…emitted as a bare gene line** | **472** | **0** |
| **`Skipping … not found` warnings** | **472** | **0** |
| Coding gene recall | 19,244 / 19,293 = 0.99746 | 0.99746 |
| Coding transcript recall | 19,318 / 19,367 = 0.99747 | 0.99747 |
| **Common-set mean protein identity** | 0.99569 | **0.99670 (+0.00101)** |
| Common set improved / regressed | — | **87 / 0** |
| `gff3-validate` errors | 0 | 0 |
| `gff3-validate` issues | 101 | 51 |
| Rows lost | — | **0** |

Recall is unchanged by design: the fix restores transcripts of genes that were already
counted, so it cannot move a gene from missing to found. The identity gain is real —
restoring a copy's model gives the evaluator a better representative for 87 reference
transcripts, and none get worse. The validator issue count halves because a childless gene
raises a `gene_has_transcripts` warning.

## The reported example

`gene-GAGE12J_1`, chrX:48,781,522-48,788,849 — the gene and coordinates in the report.

Before, the gene line alone. After:

```
chrX  LiftOn  gene  48781522  48788849  .  +  .  ID=gene-GAGE12J_1
chrX  LiftOn  mRNA  48781522  48788849  .  +  .  ID=rna-NM_001098406.4_1
chrX  LiftOn  exon  48781522  48781629  .  +  .
chrX  LiftOn  exon  48782670  48782761  .  +  .
chrX  LiftOn  exon  48783227  48783347  .  +  .
chrX  LiftOn  exon  48786170  48786295  .  +  .
chrX  LiftOn  exon  48788736  48788849  .  +  .
chrX  LiftOn  CDS   48782678  48782761  .  +  0
chrX  LiftOn  CDS   48783227  48783347  .  +  0
chrX  LiftOn  CDS   48786170  48786295  .  +  2
chrX  LiftOn  CDS   48788736  48788758  .  +  2
```

Five exons and four CDS rows — the hierarchy the report said Liftoff produces. The
recovered exon/CDS coordinate set is **identical** to this run's own `liftoff.gff3`
(9 of 9 features), so the models are Liftoff's, not reconstructed.

## Why this run needs a second, controlled pair

Running both arms fresh is faithful to the reported command, but it means the two arms do
**not** share an input: **Liftoff's `-copies` placement is not deterministic across fresh
runs**, a property of Liftoff independent of any LiftOn change (already noted in
`CLAUDE.md`). The two arms' `liftoff.gff3` differ; their `miniprot.gff3` are byte-identical.

That shows up as 492 changed rows rather than the expected 472:

| Changed rows | Count | Explanation |
|---|---:|---|
| Parents that gained children, `source=Liftoff` added | 472 | the fix — see below |
| Columns 1-8 moved | 36 | **36 of 36 had a different Liftoff input between the arms** |
| Not a parent that gained children | 20 | 5 complete `USP17L` hierarchies (gene+mRNA+exon+CDS), the chr4/chr8 tandem array |
| `dna_identity` changed | 1 | `rna-NM_001256872.1` — Liftoff placed it 14 kb apart in the two arms (chr4:9,369,413 vs 9,383,651) |

Every discrepancy traces to Liftoff's placement; none to LiftOn. The `source=Liftoff`
delta on the 472 is the provenance attribute a gene gains once a transcript is actually
processed under it — in the exact branch that used to return early.

To isolate the fix, both trees were re-run against **one shared** `liftoff.gff3` and
`miniprot.gff3` (`-L`/`-M`), the same protocol the benchmark A/B used.

<!-- CONTROL RESULT APPENDED BELOW -->
## Controlled pair — the fix in isolation

Same two trees, same MANE input, both given **one shared** `liftoff.gff3` and
`miniprot.gff3` via `-L`/`-M`. 525 s and 515 s, both exit 0.

| | off_ctl | on_ctl |
|---|---:|---:|
| Extra gene copies | 509 | 509 |
| …emitted as a bare gene line | **472** | **0** |
| `Skipping … not found` warnings | **472** | **0** |
| Duplicate non-CDS ids | 0 | 0 |

| Diff | |
|---|---:|
| Rows lost | **0** |
| Rows added | **3,351** |
| Genes that gained children | **472** |
| Child rows recovered / matching the shared Liftoff | **472 / 472** |
| Rows changed | **472** — every one a parent that gained children |
| Rows changed that are *not* such a parent | **0** |
| Rows with columns 1-8 moved | **0** |
| Attribute deltas | `source` × 472, nothing else |

With Liftoff held constant the changed-row count falls from 492 to exactly 472 and the
coordinate moves vanish, which is what the fresh-pair analysis predicted. The fix is
**purely additive**: it adds back 3,351 rows under genes that were already being emitted,
moves nothing, and loses nothing. Columns 1-8 byte-identical on every changed row means
seqid, coordinates, strand, score and phase provably did not move — so no protein can have
changed.

## Verdict

The reported bug reproduces on the reported input, the fix resolves it there, and the
change is additive under a controlled input. Recall is unchanged, common-set protein
identity improves slightly with nothing regressing, and validator errors stay at zero.

## Artifacts

Under `/ccb/salz3/kh.chao/lifton_mane_chm13/`:

| Path | What |
|---|---|
| `on/on.gff3` | the lifted annotation, MANE v1.5 → CHM13v2.0 (266,278 rows) |
| `mane_chm13_results.json` | fresh-pair arms, recall, identity, validity, gates |
| `mane_controlled_results.json` | controlled-pair diff and gates |
| `off/`, `on/`, `off_ctl/`, `on_ctl/` | per-arm `argv.json`, logs, `completion.json`, intermediates |
| `inputs/` | MANE GFF (with sha256), `chrom_mapping.txt` |
