## Rescue-extension A/B: terminal-stop completion of miniprot-derived models

Two arms per cell on cached Liftoff/miniprot inputs; only the experiment switch differs (`{'LIFTON_ORF_STOP_COMPLETION': '0'}` vs `{'LIFTON_ORF_STOP_COMPLETION': '1'}`). Scored by the neutral evaluator. Gate: 0 lost, 0 redundant, 0 regressed, validity not worse.

| cell | mode | new tx | added mean PI | added ORF-valid | lost | regr | redundant | gene recall off→on | primary gene off→on | tx recall off→on | val off→on | gate |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| celegans_to_briggsae | ladder | 0 | None | None | 0 | 0 | 0 | 0.4507→0.4507 | 0.4507→0.4507 | 0.47232→0.47232 | 0→0 | **PASS** |
| drosophila_to_anopheles | ladder | 0 | None | None | 0 | 0 | 0 | 0.40164→0.40164 | 0.40164→0.40164 | 0.34781→0.34781 | 0→0 | **PASS** |
| zebrafish_to_medaka | ladder | 0 | None | None | 0 | 0 | 0 | 0.37213→0.37213 | 0.37213→0.37213 | 0.39468→0.39468 | 0→0 | **PASS** |
| rice_to_sorghum | ladder | 0 | None | None | 0 | 0 | 0 | 0.82036→0.82036 | 0.82036→0.82036 | 0.7959→0.7959 | 0→0 | **PASS** |
| t4_human_to_xenopus | ladder | 0 | None | None | 0 | 0 | 0 | 0.70149→0.70149 | 0.70149→0.70149 | 0.61828→0.61828 | 0→0 | **PASS** |
| t4_human_to_chicken | ladder | 0 | None | None | 0 | 0 | 0 | 0.69216→0.69216 | 0.69216→0.69216 | 0.6455→0.6455 | 0→0 | **PASS** |
| human_to_mouse | ladder | 0 | None | None | 0 | 0 | 0 | 0.88433→0.88433 | 0.88433→0.88433 | 0.90797→0.90797 | 0→0 | **PASS** |
| drosophila | ladder | 0 | None | None | 0 | 0 | 0 | 0.97545→0.97545 | 0.97545→0.97545 | 0.98359→0.98359 | 0→0 | **PASS** |
| human_to_zebrafish | full | 0 | None | None | 0 | 0 | 0 | 0.51589→0.51589 | 0.59299→0.59299 | 0.45489→0.45489 | 0→0 | **PASS** |
| t4_human_to_chicken | full | 0 | None | None | 0 | 0 | 0 | 0.53121→0.53121 | 0.61499→0.61499 | 0.57168→0.57168 | 0→0 | **PASS** |
| t4_human_to_xenopus | full | 0 | None | None | 0 | 0 | 0 | 0.54963→0.54963 | 0.63761→0.63761 | 0.52302→0.52302 | 0→0 | **PASS** |
| arabidopsis_to_rice | full | 0 | None | None | 0 | 0 | 0 | 0.37508→0.37508 | 0.37508→0.37508 | 0.275→0.275 | 1→1 | **PASS** |
| t4_drosophila_to_bee | full | 0 | None | None | 0 | 0 | 0 | 0.29637→0.29637 | 0.29637→0.29637 | 0.24653→0.24653 | 0→0 | **PASS** |

### What changed, and the quality of the models it changed

`changed otherwise` must be 0: every transcript that moved did so only by a three-base terminal CDS+exon extension. The ORF columns cover the miniprot-only rescued models, the population this affects.

| cell | transcripts | extended | changed otherwise | ends in a stop off→on | starts with M off→on | ORF-valid off→on |
|---|---|---|---|---|---|---|
| celegans_to_briggsae (ladder) | 3600 | 76 | 0 | 0.6525→0.6999 | 0.6861→0.6861 | 0.459→0.4946 |
| drosophila_to_anopheles (ladder) | 2536 | 121 | 0 | 0.4017→0.4611 | 0.3362→0.3362 | 0.1611→0.1773 |
| zebrafish_to_medaka (ladder) | 1230 | 26 | 0 | 0.5895→0.615 | 0.5048→0.5048 | 0.3419→0.3578 |
| rice_to_sorghum (ladder) | 4834 | 44 | 0 | 0.4384→0.5271 | 0.5714→0.5714 | 0.2808→0.3202 |
| t4_human_to_xenopus (ladder) | 1912 | 108 | 0 | 0.5742→0.639 | 0.5345→0.5345 | 0.3447→0.3724 |
| t4_human_to_chicken (ladder) | 1992 | 57 | 0 | 0.6087→0.6544 | 0.6005→0.6005 | 0.3938→0.4216 |
| human_to_mouse (ladder) | 2804 | 3 | 0 | 0.5227→0.5398 | 0.6648→0.6648 | 0.3011→0.3125 |
| drosophila (ladder) | 7141 | 1 | 0 | 0.0→0.0 | 1.0→1.0 | 0.0→0.0 |
| human_to_zebrafish (full) | 66299 | 2943 | 0 | 0.5024→0.5495 | 0.4564→0.4564 | 0.2411→0.2623 |
| t4_human_to_chicken (full) | 82638 | 2362 | 0 | 0.585→0.6333 | 0.5526→0.5526 | 0.3291→0.3557 |
| t4_human_to_xenopus (full) | 75763 | 3552 | 0 | 0.5559→0.6106 | 0.5174→0.5174 | 0.2906→0.3189 |
| arabidopsis_to_rice (full) | 13982 | 771 | 0 | 0.3915→0.4543 | 0.327→0.327 | 0.1356→0.1562 |
| t4_drosophila_to_bee (full) | 7612 | 359 | 0 | 0.4205→0.4676 | 0.3081→0.3081 | 0.1477→0.1585 |

### Added models against the earlier rescues

Earlier rescues are the miniprot-only models already in the OFF output, the natural comparison for what the experiment adds.

| cell | added n | added mean PI | added PI≥0.5 | added ORF-valid | earlier n | earlier mean PI | earlier PI≥0.5 | earlier ORF-valid |
|---|---|---|---|---|---|---|---|---|
| celegans_to_briggsae (ladder) | 0 | None | None | None | 1013 | 0.72414 | 0.9724 | 0.459 |
| drosophila_to_anopheles (ladder) | 0 | None | None | None | 1297 | 0.62082 | 0.8196 | 0.1611 |
| zebrafish_to_medaka (ladder) | 0 | None | None | None | 626 | 0.66778 | 0.861 | 0.3419 |
| rice_to_sorghum (ladder) | 0 | None | None | None | 203 | 0.65208 | 1.0 | 0.2808 |
| t4_human_to_xenopus (ladder) | 0 | None | None | None | 1665 | 0.67635 | 0.8793 | 0.3447 |
| t4_human_to_chicken (ladder) | 0 | None | None | None | 1224 | 0.72337 | 0.9297 | 0.3938 |
| human_to_mouse (ladder) | 0 | None | None | None | 176 | 0.70677 | 1.0 | 0.3011 |
| drosophila (ladder) | 0 | None | None | None | 1 | 0.55952 | 1.0 | 0.0 |
| human_to_zebrafish (full) | 0 | None | None | None | 59147 | 0.62946 | 0.7833 | 0.2411 |
| t4_human_to_chicken (full) | 0 | None | None | None | 47335 | 0.70341 | 0.9035 | 0.3291 |
| t4_human_to_xenopus (full) | 0 | None | None | None | 62751 | 0.67248 | 0.855 | 0.2906 |
| arabidopsis_to_rice (full) | 0 | None | None | None | 6122 | 0.59523 | 0.7689 | 0.1356 |
| t4_drosophila_to_bee (full) | 0 | None | None | None | 4989 | 0.59613 | 0.7015 | 0.1477 |

### Agreement with the released target annotation

Fraction of coding models whose CDS overlaps an annotated CDS on the same strand, by model class in the ON output.

| cell | experiment-added | earlier rescues | DNA lift |
|---|---|---|---|
| human_to_zebrafish | 1.0 (n=2943) | 0.9989 (n=56368) | 0.9704 (n=6869) |
| t4_human_to_chicken | 0.9992 (n=2362) | 0.9986 (n=45050) | 0.9939 (n=33953) |
| arabidopsis_to_rice | 0.9922 (n=771) | 0.9972 (n=5738) | 0.9565 (n=7426) |
| t4_drosophila_to_bee | 0.9833 (n=359) | 0.9968 (n=4754) | 0.9924 (n=2497) |
