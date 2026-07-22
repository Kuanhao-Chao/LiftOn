# Canonical-v2 Benchmark Coverage

Canonical-v2 contains **16 campaigns**, **411 cells**, and **751 arms**. The release matrix has 40 subset, 22 full-genome, and 11 end-to-end dataset IDs, plus 2 exact synthetic cases.

| Scenario | Kind | Panels | Format | Campaigns |
|---|---|---|---|---|
| `v2_truth_human_grch38_chm13` | biological | full, subset, truth | GFF3 | full-truth, subset-truth |
| `v2_dialect_ensembl116_gtf` | biological | full, subset, truth | GTF | full-truth, subset-truth |
| `v2_dialect_flybase_dmel_dere` | biological | cross_e2e, full, subset, truth | GTF | e2e-cross-truth, full-truth, subset-truth |
| `v2_dialect_wormbase_celegans_cbriggsae` | biological | cross_e2e, full, subset, truth | GFF3 | e2e-cross-truth, full-truth, subset-truth |
| `v2_truth_soybean_w82_lee` | biological | full, subset, truth | GFF3 | full-truth, subset-truth |
| `v2_deep_zebrafish_xenopus` | biological | full, subset, truth | GFF3 | full-deep-diagnostic, subset-deep-diagnostic |
| `v2_deep_tomato_rice` | biological | full, subset, truth | GFF3 | full-deep-diagnostic, subset-deep-diagnostic |
| `v2_truth_rat_mouse_e2e` | biological | cross_e2e, subset, truth | GFF3 | e2e-cross-truth, subset-truth |
| `v2_synth_chr22_fragmented` | synthetic | synthetic | n/a | synthetic-truth |
| `v2_synth_chr22_sv` | synthetic | synthetic | n/a | synthetic-truth |
| `v2_protocol_thread_scaling_bee` | protocol | protocol | n/a | e2e-canary, e2e-same-species, full-legacy, subset-legacy, v2_protocol_thread_scaling_bee |
| `v2_protocol_io_modes_arabidopsis` | protocol | protocol | n/a | e2e-same-species, full-canary, full-legacy, subset-legacy, v2_protocol_io_modes_arabidopsis |

## Prioritized expansion backlog

- **P1 — genome-topology / `polyploid-homeologs`:** Polyploid transfer with homeolog-aware truth and copy-number scoring.
- **P1 — genome-topology / `alternate-haplotype-pangenome`:** Primary-to-alternate-haplotype or graph-derived assembly transfer.
- **P1 — difficult-biology / `segmental-duplication-paralogy`:** Segmental duplications, recent paralogs, and explicit copy truth.
- **P1 — difficult-biology / `gene-rearrangement-truth`:** Controlled translocations, fusion, splitting, inversion, and pseudogene cases.
- **P2 — annotation-content / `noncoding-and-partial-annotations`:** Noncoding-only, partial, missing-parent, and mixed gene-model inputs.
- **P2 — genome-topology / `organelle-and-circular-sequences`:** Mitochondrial, chloroplast, and circular-coordinate annotations.
- **P2 — input-robustness / `compressed-streaming-malformed-inputs`:** Compressed, reordered, very large-attribute, and intentionally malformed inputs.
- **P2 — performance-protocol / `cold-warm-storage-protocol`:** Cold and warm cache measurements across local storage classes.
- **P3 — software-compatibility / `runtime-compatibility-matrix`:** Pinned Python, DuckDB, spatial-extension, and native-backend combinations.
- **P3 — comparative-methods / `external-tool-baselines`:** Reproducible comparisons with actively maintained annotation-transfer tools.

Backlog cases remain diagnostic until every promotion requirement in `canonical_v3_backlog.json` is satisfied.
