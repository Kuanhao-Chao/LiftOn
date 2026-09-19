# Preferring the primary-assembly gene: measured NO-GO

14.6 % of GRCh38 RefSeq coding genes (3,392 of 23,307) sit on alt or fix
contigs as copies of a gene on the primary assembly. They compete for the same
target locus, and the improvement plan proposed making LiftOn prefer the
primary member deterministically, on the reading that "the model is right but
carries the copy's ID".

The premise is half true, and the half that is false is the half that mattered.

## What actually happens

On a full-GRCh38 → CHM13 lift, of the 1,611 GeneID groups that have both a
primary and an alt/fix member:

| | groups |
|---|---:|
| represented only by the alt/fix copy | **1,058** |
| represented only by the primary | 415 |
| represented by both | 84 |
| represented by neither | 54 |

So two thirds of contested groups are indeed carried by the copy. Liftoff makes
that choice, not LiftOn: `fix_overlapping_features.find_feature_to_remap`
resolves an overlap by alignment score, and an alt contig is frequently a
different haplotype that matches CHM13 better than GRCh38's primary does. The
copy usually wins because it genuinely aligns better.

## What the copy actually carries

RefSeq gives the alt copy the **same gene symbol, the same Dbxrefs, and
unsuffixed transcript accessions**:

```
chr1                  ID=gene-PLCH2     Name=PLCH2  gene=PLCH2  Dbxref=GeneID:9651,HGNC:…
chr1_KI270762v1_alt   ID=gene-PLCH2-2   Name=PLCH2  gene=PLCH2  Dbxref=GeneID:9651,HGNC:…
                      ID=rna-NM_001303012.2-2   Name=NM_001303012.2   transcript_id=NM_001303012.2
```

Only the GFF3 `ID` string carries NCBI's `-2`, and it carries it **in the
reference**. LiftOn is reproducing its input faithfully. Every identifier a
reader or a downstream tool matches on — symbol, GeneID, HGNC, MIM,
`transcript_id`, `Name` — is already correct.

## Whether it costs annotation

The one non-cosmetic risk is the transcript catalogue: an alt contig carries
its own isoforms. Comparing catalogues by accession across all 1,611 contested
groups:

| | groups |
|---|---:|
| copy has **fewer** transcripts | 231 |
| copy has **more** transcripts | 259 |
| same count | 1,121 |

4,555 accessions appear only on the primary; 5,627 only on the copy. It is a
wash, tilted slightly toward the copy.

## Decision

**NO-GO.** Forcing a primary preference would change which model wins on 1,058
groups, against a better-scoring alignment, to fix an ID suffix that the
reference itself assigns and that no biologically meaningful identifier shares
— with no measurable gain in the transcript catalogue.

## What ships instead

The real problem this investigation surfaced is that a user cannot tell their
reference has alt/fix contigs until they measure it: the CHM13 refresh had to
discover it by hand, after a first run whose results read backwards. So LiftOn
should **say so** — count alt/fix contigs in the reference at startup and record
them in the run manifest, naming the option to restrict to the primary
assembly. Additive, no algorithm change, and it addresses what actually went
wrong.

Reproduce: `a2/baseline.py` and `a2/isoform_cost.py` (evaluation-only).
