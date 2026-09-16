# JHU_LiftOn_v1.0.12_chm13v2.0.gff3 — how it was made, and what is in it

A LiftOn annotation of **T2T-CHM13v2.0**, lifted from the GRCh38.p14 RefSeq annotation.
It replaces `JHU_LiftOn_v1.0_chm13v2.0.gff3` (April 2024). Every figure below is derived
from the published file by `make_stats_report.py`.

## How it was generated

| | |
|---|---|
| LiftOn version | **v1.0.12** |
| Reference annotation | NCBI RefSeq `GCF_000001405.40-RS_2025_08` (GRCh38.p14), primary assembly only, rRNA removed |
| …gene-like features in it | 59,015 |
| Reference genome | `GCF_000001405.40_GRCh38.p14_genomic.fna` |
| Target genome | `chm13v2.0.fa` (T2T-CHM13v2.0) |
| Runtime / peak memory | 1.68 h / 34.5 GB |
| Size | 1,237,312,969 bytes |
| md5 | `8ee2604760b622e1169f5f8fbc524181` |

```bash
lifton -t 32 -g RefSeq_RS_2025_08_primary.gff -o lifton.gff3 -copies \
       chm13v2.0.fa GCF_000001405.40_GRCh38.p14_genomic.fna
```

The reference is restricted to the primary assembly (chr1–22, X, Y, M) on purpose: the
alt/fix contigs of GRCh38 carry duplicate copies of primary genes, and T2T-CHM13 has no
alt contigs, so lifting them only creates competing models at a single locus.

## Contents

**4,047,289 feature rows** across **25 sequences**, **22 feature types**.

| feature type | count |
|---|---:|
| `exon` | 2,119,501 |
| `CDS` | 1,688,385 |
| `mRNA` | 131,809 |
| `gene` | 42,675 |
| `lnc_RNA` | 30,195 |
| `pseudogene` | 17,015 |
| `transcript` | 13,571 |
| `primary_transcript` | 2,138 |
| `snoRNA` | 1,186 |
| `tRNA` | 529 |
| `snRNA` | 173 |
| `scaRNA` | 49 |
| `antisense_RNA` | 41 |

Gene-like total: **59,690** (`gene` 42,675 + `pseudogene` 17,015).

### Genes by biotype

| gene_biotype | count |
|---|---:|
| protein_coding | 20,032 |
| lncRNA | 18,053 |
| pseudogene | 15,581 |
| miRNA | 2,138 |
| transcribed_pseudogene | 1,230 |
| snoRNA | 1,186 |
| tRNA | 529 |
| V_segment_pseudogene | 279 |
| V_segment | 244 |
| snRNA | 173 |
| J_segment | 81 |
| ncRNA | 49 |
| misc_RNA | 29 |
| C_region | 21 |
| antisense_RNA | 19 |
| other | 12 |
| J_segment_pseudogene | 9 |
| C_region_pseudogene | 6 |
| scRNA | 4 |
| vault_RNA | 4 |
| Y_RNA | 4 |
| ncRNA_pseudogene | 4 |
| RNase_P_RNA | 1 |
| telomerase_RNA | 1 |
| RNase_MRP_RNA | 1 |

## Accuracy

Scored by LiftOn's neutral evaluator against the same reference annotation.

| measure | value |
|---|---:|
| Coding **gene** recall | **0.99598** (19,818 / 19,898) |
| Coding **transcript** recall | **0.99754** (130,943 / 131,266) |
| Mean **protein identity** | **0.99789** |
| Mean DNA identity | 0.99832 (median 0.99949) |
| ORF valid / starts with M / ends in stop | 0.99507 / 0.99816 / 0.99796 |
| Intron chain exact | 0.98881 |
| Exon sensitivity / specificity | 0.99296 / 0.99342 |
| Extra gene copies | 239 |

### Recovery by feature type

| feature type | in reference | recovered | fraction |
|---|---:|---:|---:|
| `exon` | 2,125,424 | 2,117,200 | 0.9961 |
| `CDS` | 1,691,722 | 1,687,064 | 0.9972 |
| `mRNA` | 131,266 | 130,943 | 0.9975 |
| `gene` | 42,141 | 41,919 | 0.9947 |
| `lnc_RNA` | 30,198 | 30,116 | 0.9973 |
| `pseudogene` | 16,874 | 16,733 | 0.9916 |
| `transcript` | 13,605 | 13,560 | 0.9967 |
| `primary_transcript` | 1,915 | 1,886 | 0.9849 |
| `snoRNA` | 1,193 | 1,184 | 0.9925 |
| `tRNA` | 451 | 446 | 0.9889 |

## Validity

| | |
|---|---:|
| `gff3-validate` errors | **0** |
| advisory warnings | 103 |
| duplicated non-CDS IDs | **0** |

## Three things to know when reading these numbers

**1. Some reference feature types are deliberately not lifted.** They are alignment
evidence and regulatory annotations, not gene models:

| type | in reference |
|---|---:|
| `match` | 149,645 |
| `biological_region` | 128,237 |
| `silencer` | 34,232 |
| `cDNA_match` | 3,277 |
| `miRNA` | 2,875 |
| `transcriptional_cis_regulatory_region` | 2,033 |

Because of these, an "all feature types" completeness figure reads
**0.35417** and is *not* a meaningful measure of this
file. The gene-model numbers above are.

**2. Mature miRNA features are absent.** They sit at the third level of the hierarchy
(gene → primary_transcript → miRNA). LiftOn emits the gene and the primary_transcript
but not the mature miRNA. The 2024 annotation has none either — long-standing
behaviour, not a regression introduced here.

**3. 162 extra gene copies carry no child features.** Every one is a
faithful passthrough — Liftoff itself emitted them without children — and **none is
attributable to LiftOn**. The run logs zero unresolved-transcript warnings.

## Versus the annotation it replaces

| | posted v1.0 (Apr 2024) | this file |
|---|---:|---:|
| feature rows | 1,908,328 | 4,047,289 |
| mRNA | 130,780 | 131,809 |
| exon | 1,666,000 | 2,119,501 |
| CDS | 0 | 1,688,385 |
| duplicated non-CDS IDs | 29,502 | **0** |
| `gff3-validate` errors | 154 | **0** |
| childless gene copies from LiftOn | 1,094 | **0** |

The 2024 file could not be scored for recall or identity: it cannot be loaded at all
(29,502 duplicated exon IDs over 447,932 rows; every gffutils strategy fails with
`UNIQUE constraint failed: features.id`). That is part of why it was replaced.

Full account of the refresh, including two corrections made along the way:
`notes/chm13_annotation_refresh_2026.md`.
