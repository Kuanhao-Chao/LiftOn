# LiftOn emits overlapping exons within one transcript

*2026-09-21. Cycle 3, item P2. The open half of issue #26.*

Two exons of one transcript cannot overlap: an overlap says there is no intron
between them, so they are one exon. Reference annotations have none. LiftOn
introduces them.

| annotation | with an overlapping exon pair | after |
|---|---:|---:|
| human → zebrafish | 96 of 69,123 (0.144 %) | **0** |
| drosophila | 47 of 33,722 | **0** |
| rice | 30 (17 of them from P1's separate bug) | **0** |
| bee | 13 of 28,112 | **0** |
| CHM13 (regenerated) | 27 of 145,394 (0.019 %) | pending |
| every reference annotation checked | **0** | — |

Not one overlapping pair, on any of those five, crosses a strand or a seqid —
so none of them is a trans-spliced model where an overlap would be legitimate.

## Two causes, not one

Counting the 27 CHM13 pairs by the size of the smaller side separates them
cleanly. All 27 carry `status=LiftOn_chaining_algorithm`.

### 1. miniprot's `stop_codon` was ingested twice (11 of 27)

miniprot writes the terminal CDS of a hit **including** its stop codon, then
repeats those same three bases as a nested `stop_codon` row:

```
CDS        1135793 1135923  ... StopCodon=1;Target=rna-XM_047432838.1 291 333
stop_codon 1135921 1135923  ... Parent=MP000228
```

Of the 3,963 `stop_codon` rows in a whole-genome human miniprot run,
**3,963 are fully contained in a sibling CDS**. Not one carries a base the CDS
list does not already have.

`lifton_utils.LiftOn_miniprot_alignment` queried
`children(featuretype=('CDS','stop_codon'))` and added **every** row as an exon
*and* as a CDS. So each redundant row became a 3 bp exon nested inside the
terminal exon, and `add_cds` overwrote the real terminal CDS with the 3 bp one.
It also fired the "spans 2 exons … each will carry its own copy" warning 8,537
times in one CHM13 run.

This also accounts for 11 of the 13 overlapping **CDS** pairs — five of which
were an exactly duplicated row.

`coreutils.drop_redundant_stop_codons` drops a `stop_codon` a sibling CDS
already covers. A stop codon genuinely outside every CDS (the GTF convention,
where the CDS stops short of it) carries sequence and is kept, so the filter
tests containment rather than deleting the featuretype.

On human → zebrafish, 95 of the 96 are the stop-codon class and exactly one is
a rebuild overlap — the distant transfer leans almost entirely on miniprot, so
almost every overlap there comes from the duplicated stop codon.

### The warning that blamed the reference

Adding the 3 bp stop codon as an exon *inside* the terminal exon meant the real
terminal CDS then overlapped two exons, so `Lifton_TRANS.add_cds` logged

> CDS 'CDS_242553' (…) spans 2 exons of MP045063; each will carry its own copy.
> **The reference model is malformed here.**

The reference was not malformed. The second exon was one LiftOn had just
created. Each redundant `stop_codon` produced two of these — one for itself and
one for the CDS it duplicated — and the whole class disappears with the fix:

| | before | after |
|---|---:|---:|
| rice | 4,468 | **0** |
| bee | 3,614 | **0** |
| drosophila | 2,632 | **0** |
| CHM13 (the regeneration run) | 8,537 | — |

Zero, not halved: there were no genuinely intron-spanning reference CDS behind
any of them.

### 2. A rebuilt exon ran into the next Liftoff exon (16 of 27)

`Lifton_TRANS.update_cds_list` Case 3 takes a rebuilt exon's end from a chained
CDS end, which can lie inside the *next* Liftoff exon; the "append any
remaining 3′ UTR exons" step then appends that exon verbatim. Nothing compared
the two.

The canonical case is POLR2A:

```
reference   CDS 7513021-7513776 (756 bp)  +  CDS 7513779-7514179 (401 bp)
lifted      exon 7417085-7418241 with CDS 7417085-7418241 (1,157 bp = 756+401)
            exon 7417841-7418678                      <- Liftoff's, appended
```

The two reference CDS blocks chain into one contiguous block in the target; the
second Liftoff exon is then appended on top of it.

`reconcile_overlapping_exons` collapses an overlapping pair into the contiguous
exon it describes. It returns the original list object untouched when nothing
overlaps, so a valid transcript keeps its exact bytes — including its exon
order, which is what lets `LIFTON_NO_CONTAINMENT_NORMALIZE=1` still reproduce
the older output. `LIFTON_NO_EXON_OVERLAP_RECONCILE=1` opts out.

It runs at the birth site rather than at the write funnel because the
best-of-outcome compare and the ORF rescue both read `self.exons` straight out
of `update_cds_list`; a pair fixed only at write time would still have been
scored and searched.

An exon holds at most one CDS. When both sides of an overlap carry one, the two
coding blocks are merged if they overlap or abut — which is what they are, and
the overlap means the protein was counting those bases twice (FAM118A
`rna-XM_024452255.2`: CDS 45813076-45813324 and 45813298-45813384, 27 bp
double-counted). If two coding blocks are genuinely disjoint, merging would
have to invent coding sequence, so the pair is left intact and the validator
reports it.

## Why the mandatory gate passed anyway

`normalize_containment` does not merge overlapping exons — and its collision
renumbering actively hides them. The pair starts out sharing an exon ID, which
the **mandatory** `duplicate_id` check would have caught as a fatal error;
renumbering gives them distinct IDs first, so the file validates clean.

And `gff3_validator` had no sibling-overlap check at any severity, for exons or
for CDS — while its module docstring had claimed "No overlapping CDS within one
transcript" since the file was written. `_check_sibling_overlap` makes that
claim true and adds the exon case.

## What the new validator check says about a whole genome

`gff3-validate` on both arms, uncapped:

| | errors before | errors after |
|---|---|---|
| rice | 19 (14 `exon_overlap`, 4 `cds_overlap`, 1 `seqid_consistency`) | **1** |
| bee | 21 (13 `exon_overlap`, 8 `cds_overlap`) | **0 — `is_valid: True`** |

Bee's lifted annotation now validates clean. Rice keeps one error, and it is a
different defect this check happened to surface:

```
rna-OrsajM_p05: seqid 'CP132246.1' differs from parent 'gene-OrsajM_p05'
                seqid 'CP132245.1'
```

`nad5` is a trans-spliced mitochondrial gene (`exception=trans-splicing`). The
gene row stayed on one sequence while its mRNA, exons and CDS were written on
another — a gene and its transcript on different sequences, which is invalid
GFF3 and impossible as written. It is present identically in both arms, so it
is not something this change introduced and not something it fixes. Recorded
here as a real, separate finding.

It is also why `reconcile_overlapping_exons` refuses to merge a pair that
disagrees on strand *or* seqid: rice genuinely contains models of that shape.
No overlapping pair in any of the five genomes crosses either, so the refusal
never fires on this corpus — it is there so that the one case where merging
would be meaningless is refused rather than performed silently.
