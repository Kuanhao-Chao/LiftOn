# v1.0.14 second pre-release scan (2026-09-24)

*Candidate before the scan: `e688ff5` (qualified 2026-09-23, CI green). After:
`3ab51cf`. Every fix below has a test that fails on `e688ff5` (N2: on
`ac7e863`) and passes on `3ab51cf`; every one is inert where its defect is
absent (the 24-cell byte-identity matrix is green with no golden edit).
`3ab51cf` differs from `ac7e863` only in the validator's phase check (N2),
which cannot change lift output: a subset lifted with each is byte-identical.
The performance campaign that found N1 and N2 is in
`notes/v1.0.14_performance_evaluation.md`.*

## Why a second scan found anything

The same blind spot as the 2026-07 audit: the test fixtures are gene-only,
coding-only, 600-bp, two-level, RefSeq-shaped. Every defect below needs an
input shape the suite did not have -- a GenBank (GCA) annotation, `-E` on an
older output, a GENCODE/Ensembl annotation, a mitochondrial protein that only
miniprot found, a malformed `transl_except`, a trans-spliced transcript.
Three read-only surveys (code risk, benchmark infrastructure, verification
data) ran first; each suspicion was then confirmed in code and reproduced
before it was fixed.

## Fixed

| id | defect | shipped since | fix (commit) | proof |
|---|---|---|---|---|
| C1 | reference normalization raised `LiftOnInputError` -- aborting the whole run -- on any sparse coding model it could not rebuild unambiguously; stock NCBI GenBank yeast (the repo's own `GCA_000146045.2_R64`, 47 Ty gag-pol frameshifts written as two CDS under the gene beside an mRNA) | v1.0.14 candidate only | ambiguous models are lifted as written (v1.0.13's treatment), warned once, counted in the manifest; `--strict-gff` keeps it fatal (`af68f4a`) | yeast R64 self-lift: `e688ff5` exit 2, no output (`Gene 'gene-YBL005W-B' mixes direct CDS with ['mRNA'] children`); `af68f4a` exit 0, 25,957 rows |
| H1 | an unsplittable CDS raised `ValueError`: in Step 7 it escaped to the per-locus handler and the whole gene was left out; in `-E` it aborted the run | v1.0.14 candidate only | the transcript is removed and skipped (still counted as `cds_spanning_exons`); `-E` also drops a stop_codon row a sibling CDS covers, which used to replace the CDS (`246df9d`) | 4 tests, serial and threaded |
| M2 | Ensembl/GENCODE selenocysteine (`stop_codon_redefined_as_selenocysteine` rows in GFF3, `Selenocysteine` in GTF -- which gffread drops, 0 of 130 surviving) not read at all: GENCODE users still got truncated selenoproteins | v1.0.14 fix was RefSeq-only | rows read as `transl_except` declarations (one query; GTF pre-scanned before conversion); output carries `transl_except` in target coordinates (`73e2762`) | GENCODE v49 Sec genes (+447 chr22 controls) GRCh38 → CHM13, same aligner output: 71 Sec transcripts 0.680 → 0.995 (GFF3) and 0.661 → 0.997 (GTF), below 0.9 48 → 0, 0 worse; 103/107 written codons all TGA; every non-Sec gene unchanged |
| M3 | the ORF search stopped every ORF at a Sec UGA and scored without read-through, so a selenoprotein with a frameshift or lost start could still be replaced by a truncated ORF | v1.0.14 candidate | scan skips declared codons; candidates scored with read-through (`09c0c78`) | unit test on the SEPHS2 shape |
| M4 | one malformed `transl_except` value (even `, aa:Sec`) discarded all of a transcript's values; only the first CDS row's values were read | v1.0.14 candidate | per-value parsing; union over CDS rows (`6bfb7b5`) | tests |
| L1 | a declaration that is neither a read-through nor a start was written wherever it aligned (`aa:TERM` onto a sense codon, dog → cat `NM_001003050.1`) | v1.0.14 candidate | written only onto an identical codon; the reference codon is recorded at Step 3 (`6bfb7b5`) | tests; dog → cat: two models drop the value their target codon no longer carries |
| L10 | a malformed `transl_except` range was materialised in memory | v1.0.14 candidate | capped at 1 Mb (`6bfb7b5`) | test (`1..6910` still parses) |
| L11 | written/not-applicable counts double-counted re-renders | v1.0.14 candidate | per-model tally (`6bfb7b5`) | test |
| L3 | a failure scanning declarations on `-P/-T` runs aborted the lift | v1.0.14 candidate | warn and continue (`6bfb7b5`) | test |
| M1 | every miniprot-derived model (Step 8, candidate scaffold, rescue, isoforms, cross-locus) used the standard genetic code even when its reference declares 2/5 -- TGA read as a stop, stop completion appending TGA | v1.0.14 candidate (declared codes are new in v1.0.14) | `Lifton_TRANS.transl_table()` falls back to the reference transcript's code from the Step-3 sidecar (`014850c`) | tests; CHM13 unchanged (its mitochondrial genes are Liftoff models); **dog → cat and human → zebrafish mitochondrial genes placed from miniprot hits**: COX2 0.259 → 0.965, CYTB 0.063 → 0.889, ND5 0.100 → 0.842 (dog → cat); COX1 0.053 → 0.848, COX3 0.260 → 0.805 (human → zebrafish) -- the false `stop_codon_gain` calls are gone |
| M5 | the second-locus sub-pass queried every hit's CDS before checking its locus: 40–55 s per whole genome | v1.0.14 candidate | occupied loci refused first (exact: the tree only grows) (`52cb07c`) | test; `miniprot_rescue_ms_subpass_c_second_locus` at `-t 1`, same input: drosophila 24.3 s → 0.49 s, rice 40.1 s → 0.60 s, dog → cat 45.8 s → 0.90 s (whole genomes; output byte-identical) |
| L4 | drop counts: rescue passes counted the same feature repeatedly; forked isoform workers' records lost (manifest differed between `-t 1` and `-t N`); report ran before cross-locus | v1.0.14 candidate | distinct (class, id); worker journals merged; report last (`aabc2fa`) | tests |
| L5 | second locus did not defer to cross-locus enabled by env | v1.0.14 candidate | (`c056644`) | test |
| L9 | parent overlay missed gffutils `FeatureNotFoundError` | v1.0.14 candidate | (`c056644`) | test |
| N2 | `gff3-validate`'s `cds_phase_consistency` assumed every transcript's first CDS starts in phase 0, so a 5'-partial model (start lost) had every later segment flagged. v1.0.14 keeps the initial phase (`9bdf8a1`), so the check reported v1.0.14's own output as far worse than v1.0.13's. Found by the performance campaign's validator comparison | v1.0.14 candidate (made visible by `9bdf8a1`) | expected phase reads the first row's phase (`3ab51cf`); warning-only | tests fail on `ac7e863`; re-validated: dog → cat 5,648 → 0, gorilla 1,465 → 0, chicken 597 → 0, drosophila 202 → 0; v1.0.13's genuine 26/0/4/1 unchanged; nothing else in any tally moves. An independent check (`phase_check.py`) finds 0 phase-inconsistent transcripts in v1.0.14's dog and gorilla outputs |
| L7 | validator overlap check ignored seqid/strand; a declared −1 ribosomal-slippage CDS overlap (RefSeq PEG10; 49 CDS rows in human RefSeq) was an ERROR | v1.0.14 candidate | per (seqid, strand); declared slippage ≤ 3 bp → WARNING (`ac7e863`) | tests |

## Verified, left as designed

- **N1** a transcript whose two Liftoff CDS rows share a base is refused
  and counted (`cds_spanning_exons`) -- 1 in the 34 subsets (C. elegans →
  C. briggsae W07G4.3), 0 on the 8 genomes, CHM13 and MANE. Found by the
  campaign's recall accounting. A fix was written (`1403433`: attach a CDS to
  the one exon that contains it) and **reverted** (`1081a2f`) when the full
  suite failed it: `test_miniprot_candidate_cannot_drop_gene` pins that
  overlapping CDS rows must be rejected, because attaching both counts the
  shared bases twice -- and W07G4.3's rows do exactly that at one base. The
  probe's recovered model was valid only because the merge won. v1.0.13
  emitted the transcript as an invalid model (overlapping exons, a duplicated
  CDS row). Keeping it needs the double-counted base resolved, which changes
  the protein: next version. Lesson: only the targeted test file was run
  before committing; the full suite caught it.
- **L8** `except OSError` around the isoform pool and the parallel lift also
  catches an error from a worker. Every worker wraps each job in
  `except Exception`, so only pool-infrastructure failures (a worker dying)
  reach it, and the fallback reruns in-process with identical output.
  Narrowing it would turn a graceful fallback into an abort.

## Deferred (next version), with reason

- **L2** a second full miniprot run for `transl_table=11` (plastid), which
  translates like table 1: performance only; merging changes miniprot's input.
- **L6** stats biotype fallback for GENCODE/ENSEMBL/CHESS: report text only.
- **L12** initial phase kept when a merged CDS's 5′ end moves: bounded by the
  best-of-outcome check; no observed case.
- Pre-existing: 6 mitochondrial cosRNAs of arabidopsis (15–46 bp `ncRNA`
  children of coding genes such as nad7 and rps3) are absent from both
  v1.0.13's and v1.0.14's arabidopsis → rice output and are not counted;
  v1.0.13 also logged them as failures, v1.0.14 does not.
- Pre-existing: two CDS rows in one exon (ribosomal slippage) -- the exon keeps
  one; organellar genes with exons listed directly under the gene get no
  reference protein (`mutation=no_protein`, e.g. rice `rpl2`).

## Reach on the qualification corpus (`e688ff5` → `ac7e863`, same aligner input)

| arm | changed transcripts | declared exception | declared genetic code | anything else |
|---|---:|---:|---:|---:|
| rice `-t 1` | 1 (0.480 → 0.509) | 1 | 0 | 0 |
| drosophila `-t 1` | 0 | 0 | 0 | 0 |
| dog → cat `-t 1` | 13, none worse | 8 | 5 | 0 |
| human → zebrafish, rescue off | 9, none worse | 4 | 5 | 0 |
| bee / rice / drosophila / dog → cat / human → zebrafish, `-t 8`, defaults | 4 / 1 / 0 / 13 / 9 | 1 / 1 / 0 / 8 / 4 | 3 / 0 / 0 / 5 / 5 | 0 |
| MANE v1.5 → CHM13 (the reporter's command) | 0 -- byte-identical | 0 | 0 | 0 |
| GRCh38 RefSeq → CHM13 (whole genome, staged file) | 0 -- byte-identical | 0 | 0 | 0 |
| chr22 installed-package lifts | 0 -- byte-identical | 0 | 0 | 0 |

Nothing was added or removed anywhere, and `-t 8` with every default on (second
locus included) attributes exactly as `-t 1` does. The "declared exception" rows on the
mitochondrial genes (ND2/ND3/ND4, COX3, ND1) are the genetic-code fix too:
those transcripts also declare a partial `TERM` stop. The dog → cat
LCE6A model recorded as a −0.007 limit yesterday is now 0.467 (it was 0.407
before the transl_except fix).

## What v1.0.13 recovered and v1.0.14 does not

The campaign compared every transcript ID in the paired outputs (same aligner
input). v1.0.14 recovers far more than it gives up, and each transcript it
gives up was traced to a cause. None of them is a defect; N1 is a refusal by design.

| where | gained | lost | why the lost ones are gone |
|---|---:|---:|---|
| 34 subsets | 246 | 13 | 12: a neighbour's model is now complete and its span covers them (below); 1: N1 (refused by design, counted) |
| human → gorilla | 21 | 0 | -- |
| human → zebrafish | 707 | 2 | 1 behind a completed neighbour, 1 behind a newly placed model |
| human → chicken | 148 | 12 | 4 SBF1 isoforms behind PPP6R2, which v1.0.13 lost whole to an empty-protein error (`ab4cbcf`); 1 behind a newly rescued selenoprotein; 3 at identity 0.390, just under the adaptive rescue floor, which v1.0.14's extra recovered genes nudged up; 4 isoforms whose widening would now reach a neighbour |
| arabidopsis → rice | 168 | 1 | behind a completed neighbour |
| bee, rice, drosophila, dog → cat | 2, 5, 2, 15 | 0 | -- |
| **8 genomes** | **1,068** | **15** | every one traced; none is a defect |

**A completed neighbour.** Bisected on the D. melanogaster → D. pseudoobscura
subset (7 of 7 present at `3a92400`, 0 at `1963d5b`): `1963d5b` stopped
miniprot's `stop_codon` row being read twice, which had corrupted the terminal
CDS of every miniprot model and so every merge. With it fixed, merges that used
to fail succeed: CG8031 0.557 → 0.922, CG1455 0.772 → 0.907, CG3953 0.763 →
0.898, CG7917 0.627 → 0.954, pug 0.857 → 0.881 -- each now with the
reference's exon structure (CG8031: 6 CDS of 49/82/69/113/112/463 bp, exactly
the reference's). A complete model spans its introns, and a gene nested in an
intron (here CG14683/CG46459 inside pug, as in the D. melanogaster reference;
five of the seven on the opposite strand) is then suppressed by Step 8's
span-overlap rule, which predates v1.0.14. Suppression by span rather than by
exon is the lever; changing it is a recall feature with its own duplicate risk
(Iteration 13/15), deferred.
