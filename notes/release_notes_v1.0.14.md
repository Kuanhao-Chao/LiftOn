# LiftOn v1.0.14

v1.0.14 is a correctness release. Every lift it produces on the seven
whole-genome runs we qualify against — including T2T-CHM13 — passes
`gff3-validate` with **zero errors**, and the validator itself now checks
several defects it used to miss. It also recovers genes a whole-genome
duplication leaves at a second locus, and honours the genetic code an
annotation declares.

## Output changes you should know about

Two defaults change the annotation relative to v1.0.13. Each has an opt-out.

- **Second-locus rescue (default on).** A reference gene can now be placed at a
  second target locus when miniprot finds it there and no emitted model
  reaches it — what a whole-genome duplication produces. Measured against
  zebrafish's own GRCz11 annotation, refusing it hid 690 real target genes on
  human → zebrafish; on rice → sorghum, 51 of 53 new placements match distinct
  protein-coding sorghum genes (shifted-locus null: maximum 2 of 1,000,
  p = 0.001). `--no-rescue-second-locus` restores v1.0.13 behaviour.
- **Declared genetic codes (`transl_table`)** are honoured when extracting,
  translating, ORF-searching and stop-completing. Translating a vertebrate
  mitochondrial CDS with the standard code reads its TGA tryptophans as stops:
  the published CHM13 annotation truncated four of the 13 human mitochondrial
  genes. Standard-code output is unchanged.

## Fixed — defects that shipped in earlier releases

- **Selenoproteins and other declared recodings (`transl_except`) are lifted
  whole.** LiftOn read a selenocysteine's UGA as a premature stop: on
  GRCh38 → CHM13 all 25 human selenoprotein genes were mis-scored or
  truncated (SEPHS2 lost its first 118 residues). Declared recoded stops are
  now read through, and `transl_except` is written in the lifted model's own
  coordinates instead of the reference's. The same fix lets a distant lift
  place selenoprotein genes it used to miss — human → zebrafish gains seven
  (GPX4, DIO1, DIO2, DIO3, SEPHS2, SELENOT, SELENOM), each on its zebrafish
  ortholog — and scores stop-readthrough isoforms correctly (drosophila: 482
  of 486 improve, none worse). Ensembl and GENCODE mark selenocysteine as
  separate rows (`Selenocysteine` in GTF, which gffread drops;
  `stop_codon_redefined_as_selenocysteine` in GFF3); those are read the same
  way — on GENCODE v49 → CHM13 the 71 selenoprotein transcripts go from
  0.66–0.68 to 0.995–0.997. The ORF search reads through declared codons
  too. Thanks to the user who reported it.
- **Models built from miniprot hits use their reference's genetic code.**
  Step 8, the rescue and its isoforms translated mitochondrial genes with the
  standard code: dog → cat COX2 0.259 → 0.965, CYTB 0.063 → 0.889.
- **A reference model LiftOn cannot rebuild no longer aborts the run.** NCBI
  GenBank annotations write a frameshifted gene as CDS rows beside an mRNA
  (yeast R64: 47 Ty genes); such a model is lifted as written and counted, and
  `--strict-gff` keeps it fatal. A CDS that cannot be split at exon
  boundaries costs its transcript, not its gene, and no longer stops `-E`.
- **`--threads 1` (the default) and `--threads N` produce the same
  annotation.** On RefSeq organellar genes the single-threaded path emitted
  duplicated exons and a doubled CDS (17 genes on rice).
- **No transcript has overlapping exons or overlapping CDS** (the open half of
  #26). On CHM13, human → zebrafish, drosophila, rice and bee, LiftOn emitted
  27, 96, 47, 30 and 13 such transcripts; now 0 on all five.
- **A CDS crossing an intron is split into exonic segments** with the correct
  per-segment phase, instead of being duplicated onto several exons or
  stretched across the intron.
- **A gene is no longer lost over one untranslatable transcript.** A
  transcript whose lifted CDS is shorter than one codon took its whole gene
  out of the annotation (dog → cat: a 37-transcript gene), in every release
  since v1.0.9.
- **Trans-spliced genes stay on their own sequence.** A transcript of a
  duplicate-ID trans-spliced gene was bound to a gene fragment that did not
  contain it — on another sequence (rice mitochondrial `nad5`) or on the same
  one (drosophila `mod(mdg4)`) — and that fragment's gene row was written at
  its child's coordinates.
- **Run status reports what actually failed.** A gene's tRNA/rRNA children
  were revisited as loci of their own and recorded as pipeline failures: 395
  false failures and a `partial_success` status on dog → cat.
- Setting `LIFTON_RESCUE_ISOFORM_WORKERS` (as LiftOn's own fork-failure
  warning suggests) no longer aborts a run whose rescued genes have no other
  isoform (v1.0.12, v1.0.13).
- The miniprot-only rescue **counts every candidate it abandons**; a gene whose
  exons RefSeq lists twice (organellar convention) is now indexed, so its
  rescue candidates are no longer silently dropped.
- A worker pool that cannot fork (strict overcommit) falls back to in-process
  work instead of aborting; a CDS naming an undeclared `Parent` no longer
  aborts a lift that previously worked.

## Validator (`gff3-validate`, `--validate-output`)

New ERROR checks: overlapping exons, overlapping CDS, and a CDS lying outside
every exon of its transcript. **A file that validated before can now report
errors.** On the reference annotations we measured (six RefSeq assemblies) the
CDS-in-exon rule reports none, and on every v1.0.14 output none of the three
fires. Overlaps are compared within one sequence and strand, and a −1
ribosomal-frameshift overlap the annotation declares (`exception=ribosomal
slippage`) is a warning. The CDS phase check now reads a 5′-partial model's
own first phase; it used to flag every later segment of such a model.

## Faster and leaner

- Step 7 roughly a tenth faster on same-species and mammalian runs
  (translation 2.1×, attribute encoding 3.2×, attribute cloning 1.7× — all
  byte-identical); the windowed aligner's anchor construction 1.3–1.9×.
- The second-locus sub-pass takes under a second per genome (it took 24–46 s).
- `--rescue-max-inflight` bounds isoform-rescue memory: 24.7 % lower peak on
  human → zebrafish for 4.4 % more wall time.

## Compared with v1.0.13

Eight whole genomes, from same-species to human → chicken, both versions on
the same aligner output and scored by one evaluator:

- mean protein identity higher on all 8 (paired: 3,484 transcripts better,
  607 worse) and on 32 of 34 single-chromosome subsets (none lower);
- 1,068 more transcripts and 15 fewer, each of the 15 traced (mostly genes
  nested in a neighbour's intron once that neighbour's model is complete);
  962 of the gains are second-locus placements, 655 of 689 of them on human →
  zebrafish confirmed by zebrafish's own annotation;
- `gff3-validate` errors from 21–1,693 per genome to 0 on all 8;
- 1.06× faster by geometric mean, with 1.1–1.4 GiB less peak memory on the
  largest distant lifts.

## T2T-CHM13 annotation

The posted GRCh38 → T2T-CHM13 annotation is regenerated with v1.0.14:
[`JHU_LiftOn_v1.0.14_chm13v2.0.gff3`](ftp://ftp.ccb.jhu.edu/pub/data/LiftOn/JHU_LiftOn_v1.0.14_chm13v2.0.gff3)
(selenoprotein transcripts 0.906 → 0.998 mean identity with 11 → 0 below 0.9,
scored by one evaluator; mitochondrial proteins at reference length 8 → 13 of
13; `gff3-validate` errors 42 → 0). A statistics sheet sits beside the file.

## Installation

`pip install lifton==1.0.14` needs no compiler: mappy is an optional extra
(`lifton[mappy]`), used only by the experimental in-process Liftoff path.
Install `minimap2` and `miniprot` separately, or use Bioconda.

Full details: [CHANGELOG.md](https://github.com/Kuanhao-Chao/LiftOn/blob/main/CHANGELOG.md).
