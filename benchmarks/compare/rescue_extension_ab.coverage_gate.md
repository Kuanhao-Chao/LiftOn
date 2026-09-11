## Rescue-extension A/B: protein-coverage rescue sub-pass (A1)

Two arms per cell on cached Liftoff/miniprot inputs; only the experiment switch differs (`{'LIFTON_RESCUE_COVERAGE_GATE': '0'}` vs `{'LIFTON_RESCUE_COVERAGE_GATE': '1'}`). Scored by the neutral evaluator. Gate: 0 lost, 0 redundant, 0 regressed, validity not worse, OFF output a byte prefix of ON.

| cell | mode | new tx | added mean PI | added ORF-valid | lost | regr | redundant | gene recall off→on | primary gene off→on | tx recall off→on | val off→on | gate |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| celegans_to_briggsae | ladder | 178 | 0.73584 | 0.4831 | 0 | 0 | 0 | 0.41986→0.4507 | 0.41986→0.4507 | 0.4151→0.43868 | 0→0 | **PASS** |
| drosophila_to_anopheles | ladder | 217 | 0.62596 | 0.1751 | 0 | 0 | 0 | 0.33821→0.40164 | 0.33821→0.40164 | 0.20852→0.23845 | 0→0 | **PASS** |
| zebrafish_to_medaka | ladder | 136 | 0.65345 | 0.3456 | 0 | 0 | 0 | 0.29781→0.37213 | 0.29781→0.37213 | 0.24159→0.28735 | 0→0 | **PASS** |
| rice_to_sorghum | ladder | 31 | 0.63754 | 0.2581 | 0 | 0 | 0 | 0.81255→0.82036 | 0.81255→0.82036 | 0.7788→0.7841 | 0→0 | **PASS** |
| t4_human_to_xenopus | ladder | 171 | 0.69674 | 0.3275 | 0 | 0 | 0 | 0.38246→0.70149 | 0.38246→0.70149 | 0.11277→0.16818 | 0→0 | **PASS** |
| t4_human_to_chicken | ladder | 164 | 0.75327 | 0.4207 | 0 | 0 | 0 | 0.38619→0.69216 | 0.38619→0.69216 | 0.26701→0.32016 | 0→0 | **PASS** |
| human_to_mouse | ladder | 9 | 0.63265 | 0.1111 | 0 | 0 | 0 | 0.86754→0.88433 | 0.86754→0.88433 | 0.86099→0.8639 | 0→0 | **PASS** |
| drosophila | ladder | 0 | None | None | 0 | 0 | 0 | 0.97545→0.97545 | 0.97545→0.97545 | 0.98359→0.98359 | 0→0 | **PASS** |
| human_to_zebrafish | full | 5638 | 0.64615 | 0.227 | 0 | 0 | 0 | 0.27383→0.51589 | 0.31536→0.59299 | 0.06725→0.10631 | 0→0 | **PASS** |
| t4_human_to_chicken | full | 5022 | 0.71516 | 0.313 | 0 | 0 | 0 | 0.3156→0.53121 | 0.3639→0.61499 | 0.25833→0.29313 | 0→0 | **PASS** |
| t4_human_to_xenopus | full | 5549 | 0.69009 | 0.2822 | 0 | 0 | 0 | 0.31139→0.54963 | 0.36053→0.63761 | 0.11224→0.15068 | 0→0 | **PASS** |
| arabidopsis_to_rice | full | 1253 | 0.6093 | 0.1484 | 0 | 0 | 0 | 0.32962→0.37508 | 0.32962→0.37508 | 0.196→0.22196 | 1→1 | **PASS** |
| t4_drosophila_to_bee | full | 790 | 0.58603 | 0.138 | 0 | 0 | 0 | 0.23988→0.29637 | 0.23988→0.29637 | 0.11341→0.13906 | 0→0 | **PASS** |

### Added models against the earlier rescues

Earlier rescues are the miniprot-only models already in the OFF output, the natural comparison for what the experiment adds.

| cell | added n | added mean PI | added PI≥0.5 | added ORF-valid | earlier n | earlier mean PI | earlier PI≥0.5 | earlier ORF-valid |
|---|---|---|---|---|---|---|---|---|
| celegans_to_briggsae (ladder) | 178 | 0.73584 | 0.9831 | 0.4831 | 581 | 0.71748 | 0.9656 | 0.4664 |
| drosophila_to_anopheles (ladder) | 217 | 0.62596 | 0.8111 | 0.1751 | 287 | 0.59148 | 0.7317 | 0.1324 |
| zebrafish_to_medaka (ladder) | 136 | 0.65345 | 0.8015 | 0.3456 | 171 | 0.65999 | 0.8596 | 0.3216 |
| rice_to_sorghum (ladder) | 31 | 0.63754 | 1.0 | 0.2581 | 103 | 0.66278 | 1.0 | 0.3107 |
| t4_human_to_xenopus (ladder) | 171 | 0.69674 | 0.8889 | 0.3275 | 105 | 0.67115 | 0.8476 | 0.4571 |
| t4_human_to_chicken (ladder) | 164 | 0.75327 | 0.9573 | 0.4207 | 56 | 0.68941 | 0.8929 | 0.2857 |
| human_to_mouse (ladder) | 9 | 0.63265 | 1.0 | 0.1111 | 31 | 0.70414 | 1.0 | 0.2903 |
| drosophila (ladder) | 0 | None | None | None | 1 | 0.55952 | 1.0 | 0.0 |
| human_to_zebrafish (full) | 5638 | 0.64615 | 0.8131 | 0.227 | 3199 | 0.62438 | 0.7643 | 0.2663 |
| t4_human_to_chicken (full) | 5022 | 0.71516 | 0.9156 | 0.313 | 2110 | 0.69844 | 0.8877 | 0.382 |
| t4_human_to_xenopus (full) | 5549 | 0.69009 | 0.885 | 0.2822 | 3463 | 0.66079 | 0.8244 | 0.3237 |
| arabidopsis_to_rice (full) | 1253 | 0.6093 | 0.8021 | 0.1484 | 2309 | 0.58787 | 0.7337 | 0.1252 |
| t4_drosophila_to_bee (full) | 790 | 0.58603 | 0.6937 | 0.138 | 889 | 0.58121 | 0.6412 | 0.1327 |

### Agreement with the released target annotation

Fraction of coding models whose CDS overlaps an annotated CDS on the same strand, by model class in the ON output.

| cell | experiment-added | earlier rescues | DNA lift |
|---|---|---|---|
| human_to_zebrafish | 0.9973 (n=5648) | 0.9969 (n=3199) | 0.9711 (n=7023) |
| t4_human_to_chicken | 0.9974 (n=5024) | 0.9981 (n=2110) | 0.9939 (n=34028) |
| arabidopsis_to_rice | 0.9968 (n=1253) | 0.9957 (n=2309) | 0.9579 (n=7813) |
| t4_drosophila_to_bee | 0.9937 (n=790) | 0.9955 (n=889) | 0.992 (n=2621) |
