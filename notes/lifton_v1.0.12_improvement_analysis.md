# LiftOn v1.0.12 — algorithm, performance, limitations, and the improvement program

Written 2026-09-11. This note records what LiftOn does, where it is accurate and where it is not, where the time goes, which edge cases it misses, and which changes v1.0.12 makes in response. Every number is re-derivable. The accuracy evidence comes from `benchmarks/compare/recall_gap_diagnosis.py` run on the shipped-default whole-genome outputs, and the per-transfer results are in `benchmarks/compare/recall_gap_diagnosis.json`. The timing evidence is the v1.0.11 biology study's `run_manifest.json` phase clocks.

## 1. The algorithm as shipped

LiftOn lifts a reference annotation onto a target genome by combining two independent homology signals and choosing, locus by locus, the one that best reproduces the reference protein.

| Step | What happens |
|---|---|
| 0–2 | Read both genomes; validate the reference GFF3; build (or reuse) its gffutils database; select gene and gene-like parents to lift. |
| 3 | Extract reference transcripts and proteins to FASTA. |
| 4 | **Liftoff** (vendored: minimap2 whole-gene DNA alignment, a DAG over aligned blocks, coordinate conversion) and **miniprot** (protein-to-genome, splice-aware) run concurrently. |
| 5–6 | Load both outputs as feature databases; map miniprot hits to reference transcripts; build per-sequence interval trees. |
| 7 | For every Liftoff locus: parasail protein and DNA alignment against the reference; *protein-maximization chaining* of Liftoff and miniprot CDS chunks; ORF rescue; a best-of-outcome choice among merge+ORF, Liftoff+ORF, and miniprot's native CDS model (candidate 3). |
| 8 | Emit miniprot-only mRNAs that overlap no lifted gene (≤ 10 %) and whose genomic span is 0.9–1.5× the reference CDS span. There is no identity floor. |
| Rescue | A separate pass after Step 8. It emits one miniprot-only gene per reference gene the DNA lift missed entirely: span ratio in (0.5, 2.0), protein identity above a floor that adapts from 0.5 down to 0.30 as DNA-lift recall falls. |
| Write | Containment normalization, the GFF3 writer, and validation. |

Liftoff's parameters are fixed: `-a 0.5 -s 0.5`, minimap2 `-a --end-bonus 5 --eqx -N 50 -p 0.5` with no `-x` preset, and `-sc 1.0` for extra copies. Nothing adapts to divergence except the rescue floor.

## 2. Accuracy

### 2.1 Close species: at the ceiling for distinct loci

Human → macaque (primary-assembly coding genes, n = 19,901):

| | LiftOn | Liftoff | miniprot |
|---|---|---|---|
| gene recall, all genes | 0.841 | 0.846 | 0.997 |
| gene recall, primary assembly | **0.933** | 0.935 | 0.998 |
| gene recall, GeneID-collapsed | **0.967** | 0.971 | 0.998 |

Of the 1,209 primary genes that miniprot recovers at protein identity ≥ 0.5 and LiftOn misses, 96 % lie at a locus already occupied by a *different* lifted gene. These are paralog families (OR, PRAMEF, CDK11A/B, HNRNPCL) with no separate locus in the target. The rescue gates reject only 45 genes for any other reason. miniprot's near-perfect recall counts overlapping paralog models separately, so it is not a distinct-locus recall.

### 2.2 Distant and very distant species: recall lost at one gate

| transfer | LiftOn primary gene recall | Liftoff | miniprot | missed genes miniprot finds (PI ≥ 0.5) | span gate | overlap gate | pseudogene filter | span ratio median | coverage of span-rejected hits (median / ≥ 0.8) | ceiling |
|---|---|---|---|---|---|---|---|---|---|---|
| human → zebrafish | 0.315 | 0.033 | 0.797 | 5,755 | **4,819 (84 %)** | 866 (15 %) | 50 | 0.30 | 1.00 / 98.5 % | +0.242 |
| human → chicken | 0.364 | 0.201 | 0.835 | 6,468 | **4,901 (76 %)** | 1,455 (22 %) | 50 | 0.30 | 1.00 / 99.2 % | +0.246 |
| human → xenopus | 0.361 | 0.057 | 0.833 | 6,398 | **5,216 (82 %)** | 1,050 (16 %) | 77 | 0.33 | 1.00 / 98.9 % | +0.262 |
| arabidopsis → rice | 0.330 | 0.016 | 0.725 | 4,953 | 1,392 (28 %) | **3,357 (68 %)** | 193 | 1.60 | 0.99 / 97.3 % | +0.051 |
| drosophila → bee | 0.240 | 0.003 | 0.566 | 1,087 | **639 (59 %)** | 384 (35 %) | 60 | 1.18 | 0.98 / 95.8 % | +0.046 |

"Ceiling" is the span-gate class as a fraction of primary coding genes: the most that fixing that gate alone can add, before rescued candidates compete for loci.

**Root cause 1: the length gate compares intron-inclusive spans across species.** The reference length is the span from the first CDS start to the end of the last-starting CDS, taken over all isoforms (`lifton_utils.get_ref_liffover_features`, `ref_features_len_dict`). It is divided into one miniprot mRNA's genomic span in Step 8 (`run_miniprot.process_miniprot`, band 0.9–1.5) and in the rescue (`miniprot_rescue`, band 0.5–2.0). Intron length scales with genome size:
- fish, bird, and frog genomes are compact relative to human, so the median ratio is 0.30–0.33;
- rice introns are longer than arabidopsis introns, so 38 % of ratios exceed 2.0.

The rejected hits cover the whole reference protein. The gate discards complete, high-identity alignments because the target's introns are a different length. It is valid only when intron lengths are conserved, that is, within a species.

**Root cause 2: the rescue collapses isoforms.** The rescue deduplicates by reference gene, so it emits one transcript per gene. For the genes it already rescues, the other reference transcripts with a co-located miniprot hit at PI ≥ 0.30 number:

| | zebrafish | chicken | xenopus | rice | bee | macaque |
|---|---|---|---|---|---|---|
| rescued genes | 3,199 | 2,110 | 3,472 | 2,309 | 889 | 36 |
| co-located isoforms left out | 16,847 | 11,635 | 20,187 | 1,948 | 1,840 | 318 |
| per rescued gene | 5.3 | 5.5 | 5.8 | 0.8 | 2.1 | 8.8 |

The same collapse explains the opt-in cross-locus replacement (`lifton/cross_locus_rescue.py`). On human → zebrafish it raised mean protein identity from 0.597 to 0.632 but lost 401 transcripts net, because it swapped multi-isoform genes for single-isoform models.

**Minor: the rescue walks candidates in coordinate order.** When two hits compete for one gene, the earlier coordinate wins, not the better hit. On human → zebrafish, 471 of 3,199 rescues (14.7 %) used a secondary (`Rank` > 1) hit. The median identity cost is zero, because these are mostly teleost ohnologs, so this is an ordering defect rather than an accuracy loss.

**Collisions between paralogs.** These are 15–22 % of vertebrate misses and 68 % on arabidopsis → rice. Recovering them needs orthology evidence. Two earlier attempts did not survive: Iteration 15 produced duplicate models, and synteny rescue had too few anchors.

### 2.3 A measurement caveat for every human-source number

GRCh38 RefSeq annotates alternate-locus and fix-patch copies as separate genes: 3,391 of 23,292 coding genes (14.6 %).
- On human → macaque, 64 % of LiftOn's "missed" genes are such copies.
- In 771 GeneID groups the target locus went to the alt or fix copy instead of the primary gene. The model is right but carries the copy's ID.
- Only 3 duplicate overlapping models were found.

Transcript and all-gene recall therefore understate human-source performance. v1.0.12 results are reported on primary-assembly and GeneID-collapsed denominators (`benchmarks/compare/gene_level.py`).

## 3. Performance

v1.0.11 biology study, `--stream --inmemory-liftoff --locus-pipeline -t 8 --native -copies`:

| transfer | total | Step 7 | aligners | reference DB | alignment DBs |
|---|---|---|---|---|---|
| bee | 914 s | 63 % | 21 % | 5 % | 4 % |
| drosophila | 1,173 s | 60 % | 28 % | 3 % | 3 % |
| dog → cat | 5,812 s | 38 % | 52 % | 3 % | 2 % |
| mouse → caroli | 7,066 s | 52 % | 35 % | 4 % | 3 % |
| human → gorilla / macaque / marmoset | 10.3–11.6 ks | 43–49 % | 36–45 % | 4–6 % | 2–3 % |

- **Default runs do not parallelize Steps 7 and 8.** `--locus-pipeline` is opt-in (`lifton.py`, `default=False`). Every benchmark and the study pass it, so the published speed is not what `lifton -t 8` delivers.
- **Step 7 is capped at about 4× on 8 threads by the GIL.** The top profiled cost is parasail's `nw_trace_scan_sat`, which releases the GIL. After it come GIL-bound Python costs: SQLite reads, `__find_orfs`, feature copying, and variant classification.
- **About two-thirds of the aligner phase is Liftoff's serial Python** (SAM parsing, DAG mapping, coordinate conversion), not minimap2. lifton2 has a byte-exact multiprocess version of that loop (3.6× on the mammalian lift loop).
- **Memory.** The 157 GiB peak was a transient spike from 12 forked Liftoff workers against a 13–16 GiB steady state. Summed RSS double-counts copy-on-write pages. The #71 change (`307abc6`) caps Liftoff workers at one per alignment task.

## 4. Limitations and edge cases

1. The length gate is biased across species, in Step 8 and in the rescue (§2.2).
2. The rescue emits one isoform per gene.
3. The rescue walks candidates in coordinate order.
4. Whole-genome-duplication targets (teleosts, polyploid plants) get one co-ortholog. Liftoff's `-copies` needs identity 1.0 by default, so it never finds diverged copies.
5. Paralog families collide at shared loci.
6. Alt or fix copies can displace the primary gene's ID.
7. miniprot-derived models are CDS-only: no UTRs, and exons equal CDS.
8. Liftoff's thresholds and minimap2 seeding do not adapt to divergence.
9. At distance, ncRNAs have DNA evidence only.
10. `--locus-pipeline` is off by default.
11. Known and deferred:
    - `start_lost` labelling (NO-GO to change)
    - the three-source strand ambiguity in `Lifton_TRANS`
    - case-1 merged exons inheriting the downstream exon ID
    - 5.64 GB validator peak on dog → cat
    - the second half of the Step-7 SQL collapse

## 5. The v1.0.12 program and its results

Promotion rule, unchanged from earlier iterations. An output-changing idea ships behind a flag and becomes the default only if a strict A/B passes on the eight-cell ladder and on whole genomes:
- 0 lost transcripts
- 0 redundant models
- 0 regressions
- validity no worse

Byte-neutral changes need byte identity on real outputs.

| # | Change | Status |
|---|---|---|
| R1 | #71 resource-aware scheduling and native-failure diagnostics | `307abc6`: full suite, real chr22 byte identity, clean-wheel smoke pass |
| E1 | gene-level, primary-assembly, and GeneID-collapsed recall (`gene_level.py`, opt-in evaluator hook) | `e282103` |
| A1 | protein-coverage rescue sub-pass | `79a212a`, default on; 13/13 A/B cells pass |
| A2 | isoform-aware rescue, scored in forked workers | `79a212a`, default on; 13/13 A/B cells pass |
| S1 | `--locus-pipeline` by default when `-t > 1` | `d4c1f02`; byte-identical on the drosophila and dog → cat whole genomes |
| S2 | Liftoff: bisect overlap test in SAM parsing, index in the GFF writer (`8816ef3`); forked parallel lift loop, default with `-t > 1` (`804f01e`) | byte-identical on drosophila and dog → cat, old code vs new, serial vs parallel |
| S4 | second half of the Step-7 SQL collapse | deferred (see 5.5) |
| — | forked workers could abort the parent's output transaction | `16f67d2` |
| — | intermediate files written to `lifton_outputliftoff/` since v1.0.10 | `764e1eb` |

Standing NO-GOs, not retried:
- Step-8 threshold relaxation (Iteration 13)
- miniprot-only models at weak Liftoff loci (Iteration 15)
- ORF best-match (Iteration 9)
- concurrent Step-5 DB builds (Iteration 11)
- Strategy B (P4)
- the `start_lost` fix
- the codon-table rewrite

### 5.1 A1 — protein-coverage rescue sub-pass

Design:
- **Placement.** Sub-pass B runs after the existing rescue (sub-pass A) and considers only candidates the length band rejected. Sub-pass A is unchanged, and B fills only loci still free in the final suppression tree, so B's output is appended and the OFF output is a byte prefix of the ON output.
- **Gates B replaces.** The band becomes two checks:
  - the fraction of the reference protein the hit aligns (`Target=<id> <start> <end>`, terminal stop excluded) must be ≥ 0.8 (`LIFTON_RESCUE_COVERAGE_MIN`);
  - the model's CDS may be at most 1.5× the reference coding length, because the adaptive identity floor (as low as 0.30) would not reject a mostly inserted model on its own.
- **Gates B keeps.** Dedup, overlap, the processed-pseudogene filter, and the identity floor are those of sub-pass A.
- **Ordering.** Candidates are tried best first: miniprot `Rank`, then `Identity`, then coverage, then position.
- **Switches.** `--coverage-rescue-gate` / `--no-coverage-rescue-gate`, overridden by `LIFTON_RESCUE_COVERAGE_GATE`. The default is off until the A/B passes.
- **Output tags.** Rescued mRNAs carry `rescue_gate=protein_coverage` and `miniprot_protein_coverage`.

Tests (`tests/test_rescue_coverage_gate.py`, 26) cover:
- coverage arithmetic and candidate order;
- flag and environment resolution;
- end to end on a fixture whose target intron is 150 bp against a 1,300 bp reference intron (span ratio 0.18, coverage 1.0):
  - rescued with the gate on, absent with it off, OFF output a byte prefix of ON;
  - a partial hit refused;
  - the CDS-length bound deciding both ways.

A/B: `benchmarks/compare/rescue_extension_ab.py --experiment coverage_gate`, 13 cells: the eight-cell ladder at `-t 1` plus five distant whole genomes at `-t 8`, all on cached aligner inputs. Every cell passes:
- 0 lost transcripts
- 0 duplicate models
- 0 identity regressions
- validator errors unchanged
- OFF output a byte prefix of ON

| cell | added transcripts | gene recall | primary gene recall | wall |
|---|---|---|---|---|
| human → zebrafish (whole genome) | 5,638 | 0.274 → 0.516 | 0.315 → **0.593** | +10.5 % |
| human → chicken (whole genome) | 5,022 | 0.316 → 0.531 | 0.364 → **0.615** | +12.8 % |
| human → xenopus (whole genome) | 5,549 | 0.311 → 0.550 | 0.361 → **0.638** | +7.4 % |
| arabidopsis → rice (whole genome) | 1,253 | 0.330 → 0.375 | 0.330 → 0.375 | +2.5 % |
| drosophila → bee (whole genome) | 790 | 0.240 → 0.296 | 0.240 → 0.296 | +8.9 % |
| human chr20 → chicken / xenopus (ladder) | 164 / 171 | 0.386 → 0.692 / 0.382 → 0.701 | — | < +10 % |
| zebrafish → medaka, drosophila → anopheles, *C. elegans* → *briggsae* | 136 / 217 / 178 | +0.074 / +0.063 / +0.031 | — | < +11 % |
| rice → sorghum, human → mouse | 31 / 9 | +0.008 / +0.017 | — | ~0 |
| drosophila (same species) | 0 | inert | — | 0 |

- **Quality of the added models.** On zebrafish their mean protein identity is 0.646, and 81 % have PI ≥ 0.5. The earlier rescues score 0.624 and 76 %.
- **Agreement with the released target annotation.** 99.7 % of added models overlap an annotated CDS on the same strand (zebrafish, chicken, rice, bee). The earlier rescues score 99.7–99.8 % and the DNA lift 95.8–99.4 %.
- **Unscored additions.** The evaluator cannot score a handful of added models per human transfer (10 on zebrafish). They are immunoglobulin and T-cell-receptor V segments, which RefSeq models as `V_gene_segment` rather than `mRNA`. The first version of the harness miscounted them as redundant. The corrected gate checks directly for duplicate IDs and reports these separately.
- **Measured recall versus the diagnostic.** The gains land at or above the diagnostic ceilings. Those were counted at PI ≥ 0.5, while the adaptive floor admits hits down to 0.30.

### 5.2 A2 — isoform-aware rescue

Design:
- **Timing.** The pass runs after all rescue placement, so it cannot change which genes are placed or where.
- **Candidates.** For each rescued gene it takes the other transcripts of that reference gene whose miniprot hits share the placed hit's sequence and strand and overlap it, using the best hit per transcript.
- **Admission.** Each candidate must clear the same identity floor. A gene is widened to cover an isoform only if the extension reaches no other gene's interval. Widened intervals join the tree, so neighbours cannot widen into each other.
- **Scoring.** Scoring is mostly GIL-bound ORF search: 59 % of it is the full-transcript DNA alignment. So it runs in forked worker processes when `-t > 1`, with the jobs shared copy-on-write. On a 1,494-isoform subset the pass takes 5.4 s instead of 14.3 s. Attachment stays serial, and the output is byte-identical at `-t 1`, `-t 8`, and to the original serial implementation on the ladder.

A/B, `--experiment isoforms`, measured with A1 on in both arms. 13/13 cells pass, and the gene set is identical in both arms of every cell.

| whole genome | added transcripts | coding transcript recall (v1.0.11 → A1 → A1 + A2) | primary-assembly transcript recall |
|---|---|---|---|
| human → zebrafish | 50,310 | 0.067 → 0.106 → **0.455** | 0.496 |
| human → chicken | 40,203 | 0.258 → 0.293 → **0.572** | 0.626 |
| human → xenopus | 53,739 | 0.112 → 0.151 → **0.523** | 0.576 |
| arabidopsis → rice | 2,560 | 0.196 → 0.222 → 0.275 | — |
| drosophila → bee | 3,310 | 0.113 → 0.139 → 0.247 | — |

On the ladder, human → mouse transcript recall rises from 0.864 to 0.908, and the same-species control is inert. The added isoforms match the earlier rescues on identity (for example 0.628 vs 0.638 on zebrafish) and on ORF validity. 99.9 % overlap an annotated CDS of the target annotation.

Cost:
The serial implementation added 24–48 % wall time on the three human-source genomes and about
2 GB of peak RSS. With the worker pool and the prefetch fix (§5.6), measured from the run
manifests of the whole-genome A/B arms:

| whole genome | rescue phase | sub-pass A | sub-pass B | isoform prefetch | isoform scoring | isoform attach | isoform total |
|---|---:|---:|---:|---:|---:|---:|---:|
| human → zebrafish | 1,016 s | 148 s | 177 s | 125 s | 120 s | 29 s | 274 s (27 %) |
| human → xenopus | 774 s | 112 s | 134 s | 107 s | 101 s | 33 s | 241 s (31 %) |
| human → chicken | 671 s | 93 s | 158 s | 97 s | 71 s | 27 s | 196 s (29 %) |
| drosophila → honey bee | 121 s | 22 s | 17 s | 6 s | 6 s | 1 s | 13 s (11 %) |
| arabidopsis → rice | 156 s | 30 s | 18 s | 4 s | 3 s | 1 s | 8 s (5 %) |

So the isoform pass now costs 27–31 % of the rescue phase on the human-source transfers and
5–11 % elsewhere — and within it, the serial **prefetch** has overtaken the pooled scoring as
the larger half. The largest serial block in the rescue is not the isoform pass at all but
**candidate placement**: sub-passes A and B together are 325 s of the zebrafish phase against
the isoform pass's 274 s.

### 5.3 S1 — `--locus-pipeline` by default

Drosophila whole genome, cached aligner inputs, `-t 8`:
- 606 s instead of 899 s, with Step 7 taking 488 s instead of 789 s.
- The output is byte-identical (md5 `6a8acadb`).
- Peak RSS rises from 0.7 to 1.9 GiB for the in-flight loci.

Dog → cat, same setup:
- 1,868 s instead of 2,792 s (−33 %), with Step 7 taking 1,462 s instead of 2,377 s.
- The output is byte-identical (md5 `3fa3597d`).
- Peak RSS rises from 1.1 to 2.8 GiB.

### 5.4 S2 — Liftoff speedups

A cProfile of a fresh drosophila run put Liftoff at 710 s of 2,244 s:
- The lift loop took 260 s, in `find_best_mapping`.
- SAM parsing took 183 s. `find_overlapping_children` recomputed every child's parent-relative span for every aligned block: 3.55 M calls and 42.6 M overlap checks.
- Writing the GFF3 took 110 s. A per-child-root scan over every parent, added with the trans-spliced-copy fix, cost 64 s.

Three changes, each exact by construction:
1. **Parsing.** Spans are computed once per alignment and overlap is tested with a bisect, which is identical because merged child intervals are disjoint. Empty blocks, from two consecutive gap operations, never overlap, as in the per-child formula; a randomized comparison caught this edge case before any real run. An aligned segment is built only for kept blocks. Equivalence: 120,000 random alignments on both strands, with and without hard clips.
2. **Writing.** An `(id, seqid)` index replaces the scan, and a property test shows it returns the same candidates in the same order.
3. **Parallel lift.** `--parallel-lift` lifts each reference chromosome in its own forked worker. This is exact because the only link between features in the primary and unmapped passes is the upstream-neighbour hint, which never crosses a reference chromosome. Workers also see earlier passes' lifted features, and unmapped parents come back as the original objects.

Fresh Liftoff on the drosophila whole genome (`-t 8`, cached miniprot):
- All four arms (old or new code, serial or parallel lift) give the same final output (md5 `a23d9edc`) and identical Liftoff GFF3 bodies.
- The aligner phase takes 261 s with old code, 201 s with the new parsing and writing (−23 %), and **134 s** adding the parallel lift (−49 %).
- Total wall drops from 928 s to 738 s, and peak RSS is unchanged at 4.3 GiB.

Dog → cat, fresh Liftoff at `-t 8`:
- All four arms produce identical intermediate Liftoff GFF3 bodies (`82f9ce19`, 1.78 M lines).
- The new-code arms produce identical final outputs (`36846eee`).
- With the new code, parallel lift takes the aligner phase from 1,765 s to 1,206 s (−32 %) and total wall from 3,710 s to 3,092 s. Peak RSS is unchanged at 13.7 GiB.
- Against v1.0.11 code, the aligner phase goes from 3,020 s to 1,765 s with the new parsing and writing alone (−42 %), and to **1,206 s** with parallel lift as well (−60 %). Total wall goes from 4,878 s to 3,092 s (−37 %).
- All four arms produce the same final output (`36846eee`) at the same peak RSS.
- Old code with parallel lift was not faster than old code serial (3,435 s), because the old parsing and writing dominated that run. That arm ran on a heavily loaded host, so its number is not a clean comparison.

Parallel lift is therefore the default whenever `-t > 1`; `--no-parallel-lift` opts out.

### 5.5 Second round — what v1.0.12 still misses, and what was done about it

The first round's diagnostic explained the gap the coverage sub-pass then closed. Pointed at
a v1.0.12 output it had nothing to say, because it replays the gates as they were *before*
that sub-pass existed. It now also replays the gates LiftOn actually ships, and reports the
ORF validity of the emitted miniprot-only models and the co-ortholog count
(`benchmarks/compare/recall_gap_after_v1012.{json,md}`, no new lift — it reads the A/B arms).

**The recall well is nearly dry.** Of the coding genes still missed that miniprot finds at
identity ≥ 0.5:

| transfer | primary gene recall | still missed | a gene already holds the locus | pseudogene filter | coverage or length bound | placed, lost to the floor or ORF search |
|---|---:|---:|---:|---:|---:|---:|
| human → zebrafish | 0.593 | 1,162 | 1,018 (88 %) | 31 | 57 | 17 |
| human → chicken | 0.615 | 1,840 | 1,669 (91 %) | 27 | 27 | 37 |
| human → xenopus | 0.638 | 1,461 | 1,259 (86 %) | 44 | 48 | 39 |
| arabidopsis → rice | 0.375 | 3,940 | 3,598 (91 %) | 180 | 32 | 10 |
| drosophila → honey bee | 0.296 | 530 | 421 (79 %) | 58 | 23 | 3 |

Before v1.0.12 the dominant class was the genomic-span band, 76–84 % of the misses on the
vertebrate transfers. That class is gone. What is left is 79–91 % one class: **miniprot places
the gene where LiftOn has already put a different gene** — a reference paralog family with no
separate locus in the target. Emitting both would duplicate a locus rather than recover a
gene, which is exactly what Iteration 15 did.

Two consequences. **Lowering the protein-coverage gate is a measured NO-GO**: the whole class
is 23–57 genes per transfer, at most +0.003 primary gene recall on zebrafish, and those are by
definition the partial hits — so no sweep was run. And the identity floor and ORF search lose
only 3–39 genes per transfer, so they are not mis-tuned either.

**Model quality was the real gap.** miniprot reports a coding alignment, so its CDS ends at
the last aligned codon and excludes the stop, while the reference convention — and every other
model LiftOn emits — includes it. `Lifton_TRANS.__find_orfs` scans the spliced *transcript*,
and a miniprot model has no UTR, so the stop sitting immediately downstream in the genome was
outside the sequence being searched. Measured on the emitted output, rescued models ended in a
stop only 39–59 % of the time and internal stops were about 1 % of the failures: it was the
termini.

`lifton/orf_completion.py` closes that. The A/B forced two design corrections, both kept:

- **After, not before, the ORF search.** Completing the model first suppressed the
  `stop_missing` mutation that had been triggering the search, so two drosophila → anopheles
  transcripts took a different ORF path and one lost identity. Running after it means the
  search sees exactly the sequence it always did, and the per-transcript shape check went from
  2 exceptions to 0.
- **Only when the reference protein ends in a stop**, and only when re-scoring shows the model
  did not get worse. Appending a residue can shift a global alignment; 1 of 122 completed
  transcripts on that cell lost 0.0013 identity, so the extension is now reverted in that case.

A/B, **13/13 PASS** — the eight-cell ladder and five whole genomes. On every cell: 0 lost,
0 duplicates, 0 regressions, validity unchanged, and **0 transcripts changed by anything other
than a three-base terminal CDS+exon extension**. The start-codon fraction is unchanged
everywhere, as it must be; only the 3′ end moves.

| whole genome | transcripts | extended | ends in a stop | ORF-valid |
|---|---:|---:|---|---|
| human → zebrafish | 66,299 | 2,943 | 0.502 → 0.550 | 0.241 → 0.262 |
| human → chicken | 82,638 | 2,362 | 0.585 → 0.633 | 0.329 → 0.356 |
| human → xenopus | 75,763 | 3,552 | 0.556 → 0.611 | 0.291 → 0.319 |
| arabidopsis → rice | 13,982 | 771 | 0.392 → 0.454 | 0.136 → 0.156 |
| drosophila → honey bee | 7,612 | 359 | 0.421 → 0.468 | 0.148 → 0.159 |

The ladder shows the same pattern (*C. elegans* → *briggsae* 0.653 → 0.700, rice → sorghum
0.438 → 0.527, human → xenopus 0.574 → 0.639), and the same-species control moves one
transcript. It is free: human → zebrafish takes 2,435 s against 2,431 s, and peak memory is
unchanged.

**Cross-locus replacement keeps its isoforms.** The opt-in pass replaces a weakly lifted gene
with a better miniprot model on another chromosome, dropping every block the weak gene had. It
emitted one transcript, so on human → zebrafish it raised mean protein identity 0.597 → 0.632
and still lost **401 transcripts net** — the reason it never became a default. It now also
receives the gene's other transcripts whose miniprot hits sit at the replacement locus, scored
by the same detached scorer the rescue's isoform pass uses, and may widen only if the wider
span reaches no other emitted model. Two further corrections: a replacement must clear the
coverage gate, so a high-identity *partial* hit cannot displace a weak but full-length lift;
and the transcript ids lose the copy suffix the gene id already loses.

**Cross-locus stays opt-in, and its motivation has shrunk.** The A/B on human →
zebrafish, both arms on the v1.0.12 default:

| | before the isoform repair | after |
|---|---:|---:|
| genes replaced | 146 | 144 |
| isoforms attached | — | 244 |
| cross-locus transcripts emitted | 146 | 388 |
| net transcripts | **−401** | **−161** |
| mean protein identity | 0.597 → 0.632 | 0.623 → 0.628 |
| apples-to-apples deficit vs miniprot | −0.067 → −0.031 | −0.0088 → −0.0038 |

The population it acts on is unchanged — 144 genes against 146, of which 140
were DNA lifts and 4 were miniprot-only rescues — and what changed is that each
replacement now brings its isoforms, so the same genes emit 388 transcripts
instead of 146. The repair cuts the transcript cost by 60 %, but the pass still
ends with fewer transcripts than it started, so it **fails the promotion gate and
stays opt-in**. Duplicate-safe, no regression on the common set, validity
unchanged.

The more interesting number is the baseline. The deficit against miniprot that
cross-locus was built to close was −0.067 when it was written; on the v1.0.12
default it is **−0.0088**, because the coverage sub-pass and the isoform rescue
already closed 87 % of it. Cross-locus now moves it to −0.0038 at a cost of 161
transcripts. A second cell was not run: the gate is decided by zebrafish, and
the motivation is no longer there to justify the compute.

**Two reported failures fixed.** A flat annotation — a prokaryotic bakta GFF, a miniprot GFF —
has top-level `CDS` rows and no `gene`, so the gene-like auto-detection found nothing, fell
back to `gene`, selected nothing, and the run died several steps later inside vendored Liftoff
with a bare "Use -f …" (GH #37). Detection now falls back to the top-level types the
annotation actually has, and an empty selection stops the run at once naming them.
`-dir/--intermediate-dir` gives a run its own artifact directory, so concurrent jobs sharing an
output directory stop sharing one `lifton_output/` (GH #14).

### 5.6 Where the rescue's time goes

The rescue is the largest phase of a distant-species run — `process_miniprot_loci` is 37–79 %
of wall on the five distant whole genomes — and the split across its passes was not visible at
all. The run manifest now records it. On the human → xenopus subset at `-t 8`: sub-pass A
2.4 s, the coverage sub-pass 2.8 s, isoform prefetch 2.0 s, isoform scoring 2.6 s (11.3 s at
`-t 1`, so the worker pool gives 4.4×), attachment 0.5 s.

That measurement redirected the speed work, and corrected two assumptions:

- The `copy.deepcopy` hot spot CLAUDE.md names as "the strongest remaining target" is **already
  fixed** — `Lifton_EXON` and `Lifton_CDS` have custom `__deepcopy__` built on
  `coreutils.clone_feature`, from the 2026-07-25 batch.
- Input fingerprinting looked like 4–8 % of wall, but it runs on a background thread and the
  measured `join_wait_seconds` is **23 µs**. It costs the run nothing; summing manifest phases
  double-counts it.
- Step 8 is already thread-parallel and smaller than the rescue (≈6.5 s against 11.2 s on that
  subset), so converting it to processes is deferred with numbers rather than built.

What the measurement did find: the isoform prefetch re-fetched each miniprot mRNA by its own
ID although the enumeration already held that row. Removing the query cut prefetch by 22 %
(xenopus) and 28 % (chicken) at `-t 8`, with byte-identical output at `-t 1` and `-t 8`.

The remaining serial cost is candidate *placement* scoring in sub-passes A and B. Those
decisions cascade through the suppression tree, but the scores do not depend on it, so scoring
could be precomputed in the existing pool and the decision loop left untouched. Not attempted
here.

### 5.7 Deferred with reasons

- **S4, the Step-7 SQL collapse.** The materialisation walkers already derive the no-level exon list from the level-1 query and split containers from terminals. The default multi-threaded path, after S1, therefore issues two queries per transcript, and the change would help only `-t 1`, by roughly 5 % of Step 7.
- **Parallel placement scoring in the rescue — now the top remaining speed lever.** Sub-passes
  A and B score candidates serially, and on whole genomes that is the single largest block in
  the rescue: 325 s of the 1,016 s zebrafish phase, ahead of the isoform pass's 274 s (§5.2).
  The *decisions* cascade through the suppression tree, but the *scores* do not depend on it,
  so scoring could be precomputed in the existing fork pool and the decision loop left
  untouched — the shape `_score_isoform` already uses. Deliberately not attempted in this
  cycle: it touches a default-on path whose byte-identity is the release's load-bearing claim,
  and it wants its own gate rather than one taken on the way to a tag.
- **Serial isoform prefetch.** After the redundant-lookup fix it is 97–125 s on the
  human-source genomes, now larger than the pooled scoring beside it. The two remaining
  queries per candidate are a `children()` call and a reference-attribute read; batching needs
  an API the default gffutils backend does not have.
- **Process-based Step 7.** `Step7StateCoordinator` gates copy-number allocation on a
  `threading.Condition` and serves live cross-locus interval reads through `_JournalTreeDict`;
  a forked child can satisfy neither. Converting it means porting Step 8's
  evaluate-then-rebase shape, which is not byte-neutral by construction.
- **Process-based Step 8.** Already thread-parallel, and smaller than the rescue (§5.6).
- **Co-ortholog models.** Miniprot hits at free loci whose reference gene LiftOn already
  emitted elsewhere: 2,460 non-overlapping loci on human → zebrafish, 318 on arabidopsis →
  rice, 15 on drosophila → honey bee. The rescue deduplicates on the reference gene id, so it
  can never place a second copy. Emitting them would change what LiftOn annotates rather than
  how well it recovers the reference — gene recall cannot move — so the claim needs
  target-annotation truth (`benchmarks/compare/target_truth.py`) before any promotion.

### 5.8 Defects found along the way

- **Forked workers could abort the parent's output transaction** (`16f67d2`). `Pool.terminate` sends SIGTERM, and the inherited handler renamed the staged GFF3 to `*.partial.gff3` while the parent was still writing. No release was exposed. The first post-Step-8 pool hit it: a run's staged output vanished and the matching md5 was a stale file. The regression test forks a child that signals itself.
- **Intermediate files in `lifton_outputliftoff/` and `lifton_outputminiprot/`** (v1.0.10 and v1.0.11, `764e1eb`).
- **The quadratic writer scan** introduced by the post-v1.0.11 trans-spliced-copy fix (S2).
