# Figure-4 outlier deep-dive (LiftOn v1.0.9 vs miniprot, common recovered coding set)

`fix_headroom` = mean over common of max(0, miniprot_pi − lifton_pi) — the additive upside a miniprot-only merge candidate would deliver.

| pair | n_common | mean Δ | regressed | improved | headroom | catastrophic (LiftOn<0.3, mini>0.7) | of which frameshift (dna≥0.6) |
|---|---|---|---|---|---|---|---|
| human_to_zebrafish | 6391 | -0.09505 | 1573 | 284 | 0.10227 | 571 (36% of reg) | 0 |
| t4_human_to_xenopus | 12175 | -0.0141 | 930 | 663 | 0.02439 | 228 (24% of reg) | 0 |
| t4_human_to_chicken | 33874 | -0.01086 | 2117 | 2243 | 0.02369 | 690 (33% of reg) | 1 |
| arabidopsis_to_rice | 7138 | -0.00848 | 510 | 173 | 0.01153 | 52 (10% of reg) | 1 |
| drosophila_to_bee | — (no eval TSV) | | | | | | |
| drosophila | 29885 | 0.00442 | 272 | 4610 | 0.0018 | 17 (6% of reg) | 8 |
