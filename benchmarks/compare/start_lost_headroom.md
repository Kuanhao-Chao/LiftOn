## `start_lost` headroom (byte-neutral diagnostic)

`variants.find_variants` decides `start_lost` from four AND-ed clauses whose first two read `align_dna.query_aln[0:3]` — the first three columns of the FULL transcript alignment, which for an mRNA with a 5'UTR is UTR sequence. When that UTR matches, the chain short-circuits and a genuine start loss is never flagged, and `start_lost` gates ORF rescue.

`LIFTON_START_LOST_DIAG` re-evaluates the same four clauses against the three columns AT the CDS start (`cds_span`, from GH #46) and records both verdicts. **newly rescued** — a start loss the shipped test misses AND that no other ORF trigger would have caught anyway — is the honest figure; **missed** overstates it.

| Dataset | divergence | scored | current | scoped | missed | **newly rescued** | % | false pos |
|---|---|---|---|---|---|---|---|---|
| celegans_to_briggsae | distant_cross_species | 5767 | 1888 | 1925 | 46 | **0** | 0.0000% | 9 |
| drosophila_to_anopheles | very_distant_cross_species | 2435 | 1534 | 1539 | 6 | **0** | 0.0000% | 1 |
| zebrafish_to_medaka | distant_cross_species | 1065 | 620 | 625 | 5 | **0** | 0.0000% | 0 |
| rice_to_sorghum | distant_cross_species | 10704 | 4654 | 4724 | 157 | **0** | 0.0000% | 87 |
| t4_human_to_xenopus | very_distant_cross_species | 637 | 404 | 401 | 0 | **0** | 0.0000% | 3 |
| t4_human_to_chicken | very_distant_cross_species | 2007 | 1043 | 1033 | 2 | **0** | 0.0000% | 12 |
| human_to_mouse | distant_cross_species | 6451 | 2154 | 2144 | 74 | **10** | 0.1550% | 84 |
| drosophila | D. melanogaster -> D. erecta | 13756 | 1437 | 1516 | 102 | **32** | 0.2326% | 23 |

**TOTAL newly-rescued: 42 of 42822 scored (0.0981%).**
