## Rescue-extension A/B: isoform-aware rescue (A2), on top of the coverage sub-pass

Two arms per cell on cached Liftoff/miniprot inputs; only the experiment switch differs (`{'LIFTON_RESCUE_COVERAGE_GATE': '1', 'LIFTON_RESCUE_ISOFORMS': '0'}` vs `{'LIFTON_RESCUE_COVERAGE_GATE': '1', 'LIFTON_RESCUE_ISOFORMS': '1'}`). Scored by the neutral evaluator. Gate: 0 lost, 0 redundant, 0 regressed, validity not worse.

| cell | mode | new tx | added mean PI | added ORF-valid | lost | regr | redundant | gene recall off→on | primary gene off→on | tx recall off→on | val off→on | gate |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| celegans_to_briggsae | ladder | 254 | 0.73117 | 0.4252 | 0 | 0 | 0 | 0.4507→0.4507 | 0.4507→0.4507 | 0.43868→0.47232 | 0→0 | **PASS** |
| drosophila_to_anopheles | ladder | 793 | 0.63004 | 0.1677 | 0 | 0 | 0 | 0.40164→0.40164 | 0.40164→0.40164 | 0.23845→0.34781 | 0→0 | **PASS** |
| zebrafish_to_medaka | ladder | 319 | 0.67807 | 0.3511 | 0 | 0 | 0 | 0.37213→0.37213 | 0.37213→0.37213 | 0.28735→0.39468 | 0→0 | **PASS** |
| rice_to_sorghum | ladder | 69 | 0.64264 | 0.2464 | 0 | 0 | 0 | 0.82036→0.82036 | 0.82036→0.82036 | 0.7841→0.7959 | 0→0 | **PASS** |
| t4_human_to_xenopus | ladder | 1389 | 0.67423 | 0.3384 | 0 | 0 | 0 | 0.70149→0.70149 | 0.70149→0.70149 | 0.16818→0.61828 | 0→0 | **PASS** |
| t4_human_to_chicken | ladder | 1004 | 0.72039 | 0.3954 | 0 | 0 | 0 | 0.69216→0.69216 | 0.69216→0.69216 | 0.32016→0.6455 | 0→0 | **PASS** |
| human_to_mouse | ladder | 136 | 0.71227 | 0.3162 | 0 | 0 | 0 | 0.88433→0.88433 | 0.88433→0.88433 | 0.8639→0.90797 | 0→0 | **PASS** |
| drosophila | ladder | 0 | None | None | 0 | 0 | 0 | 0.97545→0.97545 | 0.97545→0.97545 | 0.98359→0.98359 | 0→0 | **PASS** |
| human_to_zebrafish | full | 50310 | 0.62791 | 0.2411 | 0 | 0 | 0 | 0.51589→0.51589 | 0.59299→0.59299 | 0.10631→0.45489 | 0→0 | **PASS** |
| t4_human_to_chicken | full | 40203 | 0.7022 | 0.3284 | 0 | 0 | 0 | 0.53121→0.53121 | 0.61499→0.61499 | 0.29313→0.57168 | 0→0 | **PASS** |
| t4_human_to_xenopus | full | 53739 | 0.67141 | 0.2894 | 0 | 0 | 0 | 0.54963→0.54963 | 0.63761→0.63761 | 0.15068→0.52302 | 0→0 | **PASS** |
| arabidopsis_to_rice | full | 2560 | 0.59498 | 0.1387 | 0 | 0 | 0 | 0.37508→0.37508 | 0.37508→0.37508 | 0.22196→0.275 | 1→1 | **PASS** |
| t4_drosophila_to_bee | full | 3310 | 0.60254 | 0.1541 | 0 | 0 | 0 | 0.29637→0.29637 | 0.29637→0.29637 | 0.13906→0.24653 | 0→0 | **PASS** |

### Added models against the earlier rescues

Earlier rescues are the miniprot-only models already in the OFF output, the natural comparison for what the experiment adds.

| cell | added n | added mean PI | added PI≥0.5 | added ORF-valid | earlier n | earlier mean PI | earlier PI≥0.5 | earlier ORF-valid |
|---|---|---|---|---|---|---|---|---|
| celegans_to_briggsae (ladder) | 254 | 0.73117 | 0.9803 | 0.4252 | 759 | 0.72179 | 0.9697 | 0.4704 |
| drosophila_to_anopheles (ladder) | 793 | 0.63004 | 0.8537 | 0.1677 | 504 | 0.60633 | 0.7659 | 0.1508 |
| zebrafish_to_medaka (ladder) | 319 | 0.67807 | 0.8871 | 0.3511 | 307 | 0.6571 | 0.8339 | 0.3322 |
| rice_to_sorghum (ladder) | 69 | 0.64264 | 1.0 | 0.2464 | 134 | 0.65694 | 1.0 | 0.2985 |
| t4_human_to_xenopus (ladder) | 1389 | 0.67423 | 0.8805 | 0.3384 | 276 | 0.68701 | 0.8732 | 0.3768 |
| t4_human_to_chicken (ladder) | 1004 | 0.72039 | 0.9273 | 0.3954 | 220 | 0.73701 | 0.9409 | 0.3864 |
| human_to_mouse (ladder) | 136 | 0.71227 | 1.0 | 0.3162 | 40 | 0.68806 | 1.0 | 0.25 |
| drosophila (ladder) | 0 | None | None | None | 1 | 0.55952 | 1.0 | 0.0 |
| human_to_zebrafish (full) | 50310 | 0.62791 | 0.7812 | 0.2411 | 8837 | 0.63827 | 0.7954 | 0.2413 |
| t4_human_to_chicken (full) | 40203 | 0.7022 | 0.9028 | 0.3284 | 7132 | 0.71022 | 0.9073 | 0.3334 |
| t4_human_to_xenopus (full) | 53739 | 0.67141 | 0.8539 | 0.2894 | 9012 | 0.67883 | 0.8617 | 0.2982 |
| arabidopsis_to_rice (full) | 2560 | 0.59498 | 0.7844 | 0.1387 | 3562 | 0.59541 | 0.7577 | 0.1334 |
| t4_drosophila_to_bee (full) | 3310 | 0.60254 | 0.7196 | 0.1541 | 1679 | 0.58348 | 0.6659 | 0.1352 |

### Agreement with the released target annotation

Fraction of coding models whose CDS overlaps an annotated CDS on the same strand, by model class in the ON output.

| cell | experiment-added | earlier rescues | DNA lift |
|---|---|---|---|
| human_to_zebrafish | 0.9992 (n=50310) | 0.9972 (n=8847) | 0.9711 (n=7023) |
| t4_human_to_chicken | 0.9988 (n=40203) | 0.9976 (n=7134) | 0.9939 (n=34028) |
| arabidopsis_to_rice | 0.9992 (n=2560) | 0.9961 (n=3562) | 0.9579 (n=7813) |
| t4_drosophila_to_bee | 0.997 (n=3310) | 0.9946 (n=1679) | 0.992 (n=2621) |
