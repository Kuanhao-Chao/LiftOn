# What LiftOn still misses after v1.0.12

Produced by `benchmarks/compare/recall_gap_diagnosis.py` pointed at the v1.0.12
default whole-genome outputs (the `on` arms of
`rescue_extension_ab.py --experiment isoforms`), with miniprot's own output and
the neutral evaluator's per-transcript tables as the other inputs. Read-only: no
lift was run to produce it. Results: `recall_gap_after_v1012.json`.

## Remaining gap

For every coding gene LiftOn does not recover that miniprot recovers at protein
identity >= 0.5, the first gate LiftOn's shipped rescue fails. Counted on
primary-assembly genes.

| transfer | primary gene recall | still missed | (a) a gene already holds the locus | (b) pseudogene filter | (c) coverage or length bound | (d) placed, lost to the identity floor or ORF search |
|---|---:|---:|---:|---:|---:|---:|
| human -> zebrafish | 0.593 | 1,162 | 1,018 (88 %) | 31 | 57 | 17 |
| human -> chicken | 0.615 | 1,840 | 1,669 (91 %) | 27 | 27 | 37 |
| human -> xenopus | 0.638 | 1,461 | 1,259 (86 %) | 44 | 48 | 39 |
| arabidopsis -> rice | 0.375 | 3,940 | 3,598 (91 %) | 180 | 32 | 10 |
| drosophila -> honey bee | 0.296 | 530 | 421 (79 %) | 58 | 23 | 3 |

Before v1.0.12 the dominant class was the genomic-span band, at 76-84 % of the
misses on the vertebrate transfers. That class is gone. What is left is 79-91 %
one class: **miniprot places the gene where LiftOn has already put a different
gene.** A reference paralog family has no separate locus in the target, so
emitting both would duplicate a locus rather than recover a gene; Iteration 15
produced exactly those duplicates and synteny rescue had too few anchors.

Two consequences for planning:

- **Lowering the protein-coverage gate is not worth it.** The whole (c) class is
  23-57 genes per transfer, at most +0.003 primary gene recall on zebrafish, and
  those are by definition the partial hits. Recorded as a measured NO-GO; no
  sweep was run.
- The identity floor and ORF search lose only 3-39 genes per transfer, so they
  are not mis-tuned either.

## Rescued-model ORF validity

Of the emitted miniprot-only models, the fraction that begin with M, end in a
stop, and are fully ORF-valid:

| transfer | models | starts with M | ends in a stop | ORF-valid |
|---|---:|---:|---:|---:|
| human -> zebrafish | 59,147 | 0.456 | 0.502 | 0.241 |
| human -> chicken | 47,335 | 0.553 | 0.585 | 0.329 |
| human -> xenopus | 62,751 | 0.517 | 0.556 | 0.291 |
| arabidopsis -> rice | 6,122 | 0.327 | 0.392 | 0.136 |
| drosophila -> honey bee | 4,989 | 0.308 | 0.421 | 0.148 |

Internal stops account for about 1 % — the failures are the termini. A miniprot
CDS ends at the last aligned codon, and such a model has no UTR, so
`Lifton_TRANS.__find_orfs` scans a sequence that cannot reach the stop codon
immediately downstream. This is what `lifton/orf_completion.py` addresses.

## Co-ortholog opportunity (measured, deferred)

Miniprot hits at loci no emitted gene occupies, at identity >= 0.5 and protein
coverage >= 0.8, whose reference gene LiftOn already emitted elsewhere. The
rescue deduplicates on the reference gene id, so it can never place a second
copy — on a whole-genome-duplication target these are the co-orthologs it cannot
reach.

| transfer | hits | reference genes | non-overlapping loci |
|---|---:|---:|---:|
| human -> zebrafish | 15,144 | 2,135 | 2,460 |
| arabidopsis -> rice | 1,982 | 767 | 318 |
| human -> xenopus | 1,166 | 432 | 202 |
| human -> chicken | 1,073 | 400 | 193 |
| drosophila -> honey bee | 61 | 34 | 15 |

Emitting these would change what LiftOn annotates rather than how well it
recovers the reference, so it is recorded as an option and deliberately not
built: gene recall cannot move, and the claim would need target-annotation truth
(`benchmarks/compare/target_truth.py`) before any promotion.
