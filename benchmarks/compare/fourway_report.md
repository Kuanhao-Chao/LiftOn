# 4-way tool comparison — Liftoff vs miniprot vs LiftOn v1.0.8 vs LiftOn devel

Each tool's predicted annotation is scored by the **same version-agnostic neutral evaluator** (`benchmarks.compare.evaluator`): every predicted protein is re-aligned to the reference protein with LiftOn's own parasail kernel, and completeness is measured by reference-id recovery — identically for all four tools. The two LiftOn columns run on the **same** standalone Liftoff (`-L`) + miniprot (`-M`) outputs that the Liftoff/miniprot columns are scored on, with the **neutral common flag set only** (`-t 1 -copies -ad -g -L -M -o`; no devel-only flags, so the v1.0.8 binary is valid). So LiftOn is measured against exactly the two inputs it combines.

Provenance: LiftOn v1.0.8 = the real `v1.0.8` tag binary (`e503643d`, no `--native`) in an isolated `lifton_stable` env; LiftOn devel = current `devel` HEAD. miniprot's `completeness_feature_total` is **n/a** by construction (its `MP*` ids never match reference feature ids; its protein identity and coding-completeness are still scored).

_Records: 16 subset benchmark(s), 4 full-genome headline(s)._

## Executive summary — the two headline deltas

| benchmark | mode | class | devel−v1.0.8 meanPI | devel−v1.0.8 compl | devel−best(LO,MP) meanPI | devel−best(LO,MP) compl |
|---|---|---|---|---|---|---|
| arabidopsis | subset | same | +0.00011 | +0.00000 | +0.00026 | -0.00119 |
| bee | subset | same | +0.00000 | +0.00000 | +0.00924 | -0.00064 |
| human_ensembl_gtf | subset | same | +0.00015 | +0.00000 | +0.00665 | +0.00000 |
| human_mane | subset | same | +0.00000 | +0.00000 | +0.00117 | -0.00236 |
| human_pseudogene_stress | subset | same | +0.00066 | +0.00000 | +0.00316 | +0.00000 |
| mouse | subset | same | +0.00000 | +0.00000 | +0.00382 | -0.01750 |
| rice | subset | same | +0.00000 | +0.00000 | +0.00041 | +0.00000 |
| candida_albicans_to_dubliniensis | subset | close_cross_species | +0.00137 | +0.00000 | +0.01426 | -0.01035 |
| celegans_to_briggsae | subset | distant_cross_species | +0.00083 | +0.00517 | +0.07210 | -0.40040 |
| chicken_to_quail | subset | distant_cross_species | +0.00755 | +0.00000 | +0.04987 | -0.02393 |
| drosophila | subset | cross_species | +0.00368 | +0.00000 | +0.00499 | -0.00676 |
| fly_mel_to_pseudoobscura | subset | distant_cross_species | +0.01362 | +0.00000 | +0.00623 | -0.12646 |
| human_to_chimp | subset | cross_species | +0.00000 | +0.00000 | +0.01529 | -0.00904 |
| human_to_mouse | subset | distant_cross_species | +0.01306 | +0.00000 | -0.00931 | -0.10921 |
| mouse_to_rat | subset | cross_species | +0.00399 | +0.00000 | +0.05609 | -0.04183 |
| yeast_cerevisiae_to_paradoxus | subset | close_cross_species | +0.00009 | +0.00131 | +0.00873 | -0.01955 |
| arabidopsis | full | same | +0.00020 | +0.71665 | +0.00068 | -0.00037 |
| bee | full | same | +0.00009 | +0.00000 | +0.00962 | -0.00307 |
| rice | full | same | -0.00030 | +0.22532 | +0.00175 | -0.00082 |
| drosophila | full | cross_species | +0.00134 | -0.00319 | +0.00925 | -0.02030 |

## Subset 4-way matrix

### arabidopsis — Arabidopsis TAIR10 -> ASM2311539v1  (same-species; RefSeq; n_coding=12653)

| tool | completeness_coding | feature_total | mean PI | median PI | %identical | n recovered | wall (s) | peak RSS (MB) | errors/warn |
|---|---|---|---|---|---|---|---|---|---|
| Liftoff | 0.99834 | 0.91601 | 0.99871 | 1 | 0.99533 | 12632 | 68.8 | 2281 | 100/100 |
| miniprot | 0.99953 | n/a | 0.99585 | 1 | 0.92283 | 12647 | 8 | 402 | 50/200 |
| LiftOn v1.0.8 | 0.99834 | 0.91601 | 0.99886 | 1 | 0.99533 | 12632 | 208.5 | 1335 | 50/18 |
| LiftOn devel | 0.99834 | 0.91601 | 0.99897 | 1 | 0.99533 | 12632 | 190.5 | 538 | 50/18 |

- **LiftOn devel − v1.0.8:** mean PI +0.00011, completeness +0.00000, n_recovered +0
- **LiftOn devel − best(Liftoff, miniprot):** mean PI +0.00026, completeness -0.00119
- **LiftOn devel speedup vs v1.0.8:** 1.09x

Feature-type recovery where LiftOn devel ≠ v1.0.8 (gene-like lift):

| feature type | n_reference | v1.0.8 recovered | devel recovered | |
|---|---|---|---|---|
| CDS | 77644 | 12575 | 12580 | |
| exon | 87301 | 85056 | 85168 | |

### bee — Honey bee HAv3.1 -> ASM1932182v1  (same-species; RefSeq; n_coding=3130)

| tool | completeness_coding | feature_total | mean PI | median PI | %identical | n recovered | wall (s) | peak RSS (MB) | errors/warn |
|---|---|---|---|---|---|---|---|---|---|
| Liftoff | 0.99744 | 0.96665 | 0.98236 | 1 | 0.72998 | 3122 | 32.8 | 1267 | 100/100 |
| miniprot | 0.99808 | n/a | 0.98063 | 1 | 0.63284 | 3124 | 4.3 | 1532 | 50/200 |
| LiftOn v1.0.8 | 0.99744 | 0.96665 | 0.9916 | 1 | 0.73062 | 3122 | 215.2 | 10892 | 50/0 |
| LiftOn devel | 0.99744 | 0.96665 | 0.9916 | 1 | 0.73062 | 3122 | 95.6 | 464 | 50/0 |

- **LiftOn devel − v1.0.8:** mean PI +0.00000, completeness +0.00000, n_recovered +0
- **LiftOn devel − best(Liftoff, miniprot):** mean PI +0.00924, completeness -0.00064
- **LiftOn devel speedup vs v1.0.8:** 2.25x

### human_ensembl_gtf — Human GRCh38 chr22 (Ensembl GTF) -> T2T-CHM13  (same-species · _annotation_source_ensembl_gtf_; ENSEMBL; n_coding=2371)

| tool | completeness_coding | feature_total | mean PI | median PI | %identical | n recovered | wall (s) | peak RSS (MB) | errors/warn |
|---|---|---|---|---|---|---|---|---|---|
| Liftoff | 0.99452 | 0.99029 | 0.94419 | 1 | 0.76039 | 2358 | 39.4 | 1118 | 50/100 |
| miniprot | 0.9941 | n/a | 0.9407 | 1 | 0.55537 | 2357 | 4 | 722 | 50/158 |
| LiftOn v1.0.8 | 0.99452 | 0.99058 | 0.95069 | 1 | 0.69042 | 2358 | 124.3 | 9466 | 93/4 |
| LiftOn devel | 0.99452 | 0.99058 | 0.95084 | 1 | 0.69042 | 2358 | 77 | 388 | 86/4 |

- **LiftOn devel − v1.0.8:** mean PI +0.00015, completeness +0.00000, n_recovered +0
- **LiftOn devel − best(Liftoff, miniprot):** mean PI +0.00665, completeness +0.00000
- **LiftOn devel speedup vs v1.0.8:** 1.61x

Feature-type recovery where LiftOn devel ≠ v1.0.8 (gene-like lift):

| feature type | n_reference | v1.0.8 recovered | devel recovered | |
|---|---|---|---|---|
| CDS | 18381 | 14803 | 14826 | |
| exon | 34428 | 30268 | 30301 | |

### human_mane — Human GRCh38 -> CHM13 (MANE)  (same-species; RefSeq; n_coding=423)

| tool | completeness_coding | feature_total | mean PI | median PI | %identical | n recovered | wall (s) | peak RSS (MB) | errors/warn |
|---|---|---|---|---|---|---|---|---|---|
| Liftoff | 0.99764 | 0.99764 | 0.99308 | 1 | 0.72275 | 422 | 15.3 | 553 | 100/100 |
| miniprot | 1 | n/a | 0.98974 | 1 | 0.60993 | 423 | 2.3 | 722 | 50/157 |
| LiftOn v1.0.8 | 0.99764 | 0.99764 | 0.99425 | 1 | 0.72275 | 422 | 22.4 | 1447 | 50/10 |
| LiftOn devel | 0.99764 | 0.99764 | 0.99425 | 1 | 0.72275 | 422 | 14.7 | 366 | 50/10 |

- **LiftOn devel − v1.0.8:** mean PI +0.00000, completeness +0.00000, n_recovered +0
- **LiftOn devel − best(Liftoff, miniprot):** mean PI +0.00117, completeness -0.00236
- **LiftOn devel speedup vs v1.0.8:** 1.53x

### human_pseudogene_stress — Human GRCh38 chr19 (pseudogene/ncRNA-rich) -> T2T-CHM13  (same-species · _feature_type_stress_; RefSeq; n_coding=7518)

| tool | completeness_coding | feature_total | mean PI | median PI | %identical | n recovered | wall (s) | peak RSS (MB) | errors/warn |
|---|---|---|---|---|---|---|---|---|---|
| Liftoff | 1 | 0.48537 | 0.99318 | 1 | 0.73916 | 7518 | 61.1 | 176 | 100/103 |
| miniprot | 1 | n/a | 0.98983 | 1 | 0.56371 | 7518 | 12.2 | 1131 | 50/200 |
| LiftOn v1.0.8 | 1 | 0.48537 | 0.99568 | 1 | 0.73889 | 7518 | 413.8 | 783 | 50/39 |
| LiftOn devel | 1 | 0.48537 | 0.99634 | 1 | 0.73889 | 7518 | 249.3 | 509 | 50/34 |

- **LiftOn devel − v1.0.8:** mean PI +0.00066, completeness +0.00000, n_recovered +0
- **LiftOn devel − best(Liftoff, miniprot):** mean PI +0.00316, completeness +0.00000
- **LiftOn devel speedup vs v1.0.8:** 1.66x

Feature-type recovery where LiftOn devel ≠ v1.0.8 (gene-like lift):

| feature type | n_reference | v1.0.8 recovered | devel recovered | |
|---|---|---|---|---|
| CDS | 68317 | 48454 | 48859 | |
| exon | 88920 | 71793 | 72162 | |

### mouse — Mouse GRCm39 -> NOD_SCID  (same-species; RefSeq; n_coding=8170)

| tool | completeness_coding | feature_total | mean PI | median PI | %identical | n recovered | wall (s) | peak RSS (MB) | errors/warn |
|---|---|---|---|---|---|---|---|---|---|
| Liftoff | 0.93966 | 0.78821 | 0.99026 | 1 | 0.74439 | 7677 | 152.3 | 3308 | 100/113 |
| miniprot | 0.95716 | n/a | 0.97958 | 1 | 0.60806 | 7820 | 33.8 | 2082 | 50/200 |
| LiftOn v1.0.8 | 0.93966 | 0.78821 | 0.99408 | 1 | 0.74439 | 7677 | 1831.7 | 44612 | 50/16 |
| LiftOn devel | 0.93966 | 0.78821 | 0.99408 | 1 | 0.74439 | 7677 | 352.3 | 758 | 50/16 |

- **LiftOn devel − v1.0.8:** mean PI +0.00000, completeness +0.00000, n_recovered +0
- **LiftOn devel − best(Liftoff, miniprot):** mean PI +0.00382, completeness -0.01750
- **LiftOn devel speedup vs v1.0.8:** 5.20x

### rice — Rice IRGSP -> ASM3414082v1  (same-species; RefSeq; n_coding=5850)

| tool | completeness_coding | feature_total | mean PI | median PI | %identical | n recovered | wall (s) | peak RSS (MB) | errors/warn |
|---|---|---|---|---|---|---|---|---|---|
| Liftoff | 1 | 0.97572 | 0.99864 | 1 | 0.99538 | 5850 | 27.9 | 927 | 100/100 |
| miniprot | 1 | n/a | 0.9951 | 1 | 0.91573 | 5850 | 4.3 | 606 | 50/186 |
| LiftOn v1.0.8 | 1 | 0.97572 | 0.99905 | 1 | 0.99538 | 5850 | 120.9 | 1445 | 50/10 |
| LiftOn devel | 1 | 0.97572 | 0.99905 | 1 | 0.99538 | 5850 | 106.4 | 615 | 50/10 |

- **LiftOn devel − v1.0.8:** mean PI +0.00000, completeness +0.00000, n_recovered +0
- **LiftOn devel − best(Liftoff, miniprot):** mean PI +0.00041, completeness +0.00000
- **LiftOn devel speedup vs v1.0.8:** 1.14x

### candida_albicans_to_dubliniensis — C. albicans -> C. dubliniensis  (close_cross_species · _new_taxa_fungi_; RefSeq; n_coding=1353)

| tool | completeness_coding | feature_total | mean PI | median PI | %identical | n recovered | wall (s) | peak RSS (MB) | errors/warn |
|---|---|---|---|---|---|---|---|---|---|
| Liftoff | 0.95787 | 0.93944 | 0.87204 | 0.92011 | 0.02006 | 1296 | 11.9 | 211 | 100/100 |
| miniprot | 0.97635 | n/a | 0.88266 | 0.91971 | 0.01514 | 1321 | 1.1 | 237 | 50/152 |
| LiftOn v1.0.8 | 0.966 | 0.94718 | 0.89555 | 0.91979 | 0.01989 | 1307 | 17.3 | 27 | 43/16 |
| LiftOn devel | 0.966 | 0.94718 | 0.89692 | 0.91979 | 0.01989 | 1307 | 17.6 | 34 | 35/16 |

- **LiftOn devel − v1.0.8:** mean PI +0.00137, completeness +0.00000, n_recovered +0
- **LiftOn devel − best(Liftoff, miniprot):** mean PI +0.01426, completeness -0.01035
- **LiftOn devel speedup vs v1.0.8:** 0.98x

Feature-type recovery where LiftOn devel ≠ v1.0.8 (gene-like lift):

| feature type | n_reference | v1.0.8 recovered | devel recovered | |
|---|---|---|---|---|
| CDS | 1463 | 42 | 62 | |

### celegans_to_briggsae — C. elegans -> C. briggsae  (distant_cross_species · _new_taxa_nematode_; RefSeq; n_coding=7550)

| tool | completeness_coding | feature_total | mean PI | median PI | %identical | n recovered | wall (s) | peak RSS (MB) | errors/warn |
|---|---|---|---|---|---|---|---|---|---|
| Liftoff | 0.22278 | 0.20275 | 0.41024 | 0.32435 | 0.00716 | 1682 | 79.7 | 3038 | 100/106 |
| miniprot | 0.73921 | n/a | 0.64726 | 0.70505 | 0.00466 | 5581 | 25.5 | 3126 | 50/200 |
| LiftOn v1.0.8 | 0.33364 | 0.29655 | 0.71853 | 0.76376 | 0.00517 | 2519 | 176.7 | 7551 | 56/119 |
| LiftOn devel | 0.33881 | 0.30093 | 0.71936 | 0.76263 | 0.00509 | 2558 | 68.5 | 5 | 54/99 |

- **LiftOn devel − v1.0.8:** mean PI +0.00083, completeness +0.00517, n_recovered +39
- **LiftOn devel − best(Liftoff, miniprot):** mean PI +0.07210, completeness -0.40040
- **LiftOn devel speedup vs v1.0.8:** 2.58x

Feature-type recovery where LiftOn devel ≠ v1.0.8 (gene-like lift):

| feature type | n_reference | v1.0.8 recovered | devel recovered | |
|---|---|---|---|---|
| CDS | 43432 | 408 | 1040 | |
| exon | 47983 | 4251 | 4820 | |
| gene | 7004 | 2266 | 2258 | |
| mRNA | 7550 | 2519 | 2558 | |
| pseudogene | 867 | 0 | 47 | ⬅ gene-like |

### chicken_to_quail — Chicken GRCg7b -> Japanese quail  (distant_cross_species · _new_taxa_bird_; RefSeq; n_coding=9068)

| tool | completeness_coding | feature_total | mean PI | median PI | %identical | n recovered | wall (s) | peak RSS (MB) | errors/warn |
|---|---|---|---|---|---|---|---|---|---|
| Liftoff | 0.96967 | 0.91032 | 0.82107 | 0.97116 | 0.06074 | 8793 | 422.9 | 15373 | 100/134 |
| miniprot | 0.99382 | n/a | 0.86142 | 0.97115 | 0.05282 | 9012 | 26.9 | 2016 | 50/200 |
| LiftOn v1.0.8 | 0.96989 | 0.91058 | 0.90374 | 0.97247 | 0.06173 | 8795 | 1634.7 | 10417 | 54/102 |
| LiftOn devel | 0.96989 | 0.91058 | 0.91129 | 0.9726 | 0.06173 | 8795 | 466.5 | 667 | 54/77 |

- **LiftOn devel − v1.0.8:** mean PI +0.00755, completeness +0.00000, n_recovered +0
- **LiftOn devel − best(Liftoff, miniprot):** mean PI +0.04987, completeness -0.02393
- **LiftOn devel speedup vs v1.0.8:** 3.50x

Feature-type recovery where LiftOn devel ≠ v1.0.8 (gene-like lift):

| feature type | n_reference | v1.0.8 recovered | devel recovered | |
|---|---|---|---|---|
| CDS | 121208 | 4918 | 10984 | |
| exon | 144266 | 40849 | 47055 | |

### drosophila — D. melanogaster -> D. erecta  (cross_species; RefSeq; n_coding=7251)

| tool | completeness_coding | feature_total | mean PI | median PI | %identical | n recovered | wall (s) | peak RSS (MB) | errors/warn |
|---|---|---|---|---|---|---|---|---|---|
| Liftoff | 0.98166 | 0.90811 | 0.90845 | 0.95789 | 0.04767 | 7118 | 118 | 3127 | 131/107 |
| miniprot | 0.99021 | n/a | 0.92079 | 0.95789 | 0.04095 | 7180 | 6.6 | 431 | 50/200 |
| LiftOn v1.0.8 | 0.98345 | 0.91009 | 0.9221 | 0.95793 | 0.04758 | 7131 | 252.6 | 2449 | 52/78 |
| LiftOn devel | 0.98345 | 0.91009 | 0.92578 | 0.958 | 0.04758 | 7131 | 178.1 | 646 | 51/63 |

- **LiftOn devel − v1.0.8:** mean PI +0.00368, completeness +0.00000, n_recovered +0
- **LiftOn devel − best(Liftoff, miniprot):** mean PI +0.00499, completeness -0.00676
- **LiftOn devel speedup vs v1.0.8:** 1.42x

Feature-type recovery where LiftOn devel ≠ v1.0.8 (gene-like lift):

| feature type | n_reference | v1.0.8 recovered | devel recovered | |
|---|---|---|---|---|
| CDS | 37963 | 1761 | 2182 | |
| exon | 44245 | 19968 | 20233 | |

### fly_mel_to_pseudoobscura — D. melanogaster -> D. pseudoobscura  (distant_cross_species · _divergence_ladder_distant_insect_; RefSeq; n_coding=7251)

| tool | completeness_coding | feature_total | mean PI | median PI | %identical | n recovered | wall (s) | peak RSS (MB) | errors/warn |
|---|---|---|---|---|---|---|---|---|---|
| Liftoff | 0.71438 | 0.61709 | 0.63161 | 0.80933 | 0.00904 | 5180 | 237.9 | 3696 | 131/150 |
| miniprot | 0.93297 | n/a | 0.78948 | 0.83179 | 0.00813 | 6765 | 13.4 | 1263 | 50/200 |
| LiftOn v1.0.8 | 0.80651 | 0.71863 | 0.78209 | 0.82167 | 0.0085 | 5848 | 237.2 | 2448 | 74/153 |
| LiftOn devel | 0.80651 | 0.71863 | 0.79571 | 0.82737 | 0.0085 | 5848 | 199.6 | 583 | 68/137 |

- **LiftOn devel − v1.0.8:** mean PI +0.01362, completeness +0.00000, n_recovered +0
- **LiftOn devel − best(Liftoff, miniprot):** mean PI +0.00623, completeness -0.12646
- **LiftOn devel speedup vs v1.0.8:** 1.19x

Feature-type recovery where LiftOn devel ≠ v1.0.8 (gene-like lift):

| feature type | n_reference | v1.0.8 recovered | devel recovered | |
|---|---|---|---|---|
| CDS | 37963 | 656 | 2658 | |
| exon | 44245 | 12695 | 13985 | |

### human_to_chimp — Human GRCh38 -> Chimp NHGRI_mPanTro3-v1.1  (cross_species; RefSeq; n_coding=12842)

| tool | completeness_coding | feature_total | mean PI | median PI | %identical | n recovered | wall (s) | peak RSS (MB) | errors/warn |
|---|---|---|---|---|---|---|---|---|---|
| Liftoff | 0.99034 | 0.32282 | 0.96489 | 0.99516 | 0.2065 | 12718 | 301.4 | 4745 | 100/116 |
| miniprot | 0.99938 | n/a | 0.96528 | 0.99392 | 0.15763 | 12834 | 60.5 | 6075 | 50/200 |
| LiftOn v1.0.8 | 0.99034 | 0.32282 | 0.98057 | 0.99517 | 0.2069 | 12718 | 1121.8 | 11058 | 50/66 |
| LiftOn devel | 0.99034 | 0.32282 | 0.98057 | 0.99517 | 0.2069 | 12718 | 464.4 | 675 | 50/66 |

- **LiftOn devel − v1.0.8:** mean PI +0.00000, completeness +0.00000, n_recovered +0
- **LiftOn devel − best(Liftoff, miniprot):** mean PI +0.01529, completeness -0.00904
- **LiftOn devel speedup vs v1.0.8:** 2.42x

Feature-type recovery where LiftOn devel ≠ v1.0.8 (gene-like lift):

| feature type | n_reference | v1.0.8 recovered | devel recovered | |
|---|---|---|---|---|
| CDS | 166690 | 163377 | 163361 | |

### human_to_mouse — Human GRCh38 chr20 -> Mouse GRCm39  (distant_cross_species · _divergence_ladder_distant_mammal_; RefSeq; n_coding=3086)

| tool | completeness_coding | feature_total | mean PI | median PI | %identical | n recovered | wall (s) | peak RSS (MB) | errors/warn |
|---|---|---|---|---|---|---|---|---|---|
| Liftoff | 0.84413 | 0.24948 | 0.49835 | 0.5981 | 0.01197 | 2605 | 282.2 | 40603 | 100/115 |
| miniprot | 0.96047 | n/a | 0.79104 | 0.86492 | 0.01215 | 2964 | 17.4 | 3172 | 50/200 |
| LiftOn v1.0.8 | 0.85126 | 0.25257 | 0.76867 | 0.84229 | 0.01225 | 2627 | 154.7 | 6911 | 51/104 |
| LiftOn devel | 0.85126 | 0.25257 | 0.78173 | 0.84536 | 0.01225 | 2627 | 107.1 | 604 | 51/78 |

- **LiftOn devel − v1.0.8:** mean PI +0.01306, completeness +0.00000, n_recovered +0
- **LiftOn devel − best(Liftoff, miniprot):** mean PI -0.00931, completeness -0.10921
- **LiftOn devel speedup vs v1.0.8:** 1.44x

Feature-type recovery where LiftOn devel ≠ v1.0.8 (gene-like lift):

| feature type | n_reference | v1.0.8 recovered | devel recovered | |
|---|---|---|---|---|
| CDS | 32174 | 454 | 1756 | |
| exon | 42253 | 8554 | 9660 | |

### mouse_to_rat — Mouse GRCm39 chr18 -> Rat mRatBN7.2  (cross_species; RefSeq; n_coding=2343)

| tool | completeness_coding | feature_total | mean PI | median PI | %identical | n recovered | wall (s) | peak RSS (MB) | errors/warn |
|---|---|---|---|---|---|---|---|---|---|
| Liftoff | 0.90738 | 0.74875 | 0.85311 | 0.95814 | 0.05705 | 2126 | 161.2 | 3766 | 100/105 |
| miniprot | 0.94921 | n/a | 0.85239 | 0.95472 | 0.05621 | 2224 | 10.7 | 1039 | 50/200 |
| LiftOn v1.0.8 | 0.90738 | 0.74875 | 0.90521 | 0.95946 | 0.05705 | 2126 | 243.2 | 1062 | 50/8 |
| LiftOn devel | 0.90738 | 0.74875 | 0.9092 | 0.9596 | 0.05705 | 2126 | 119.9 | 455 | 50/5 |

- **LiftOn devel − v1.0.8:** mean PI +0.00399, completeness +0.00000, n_recovered +0
- **LiftOn devel − best(Liftoff, miniprot):** mean PI +0.05609, completeness -0.04183
- **LiftOn devel speedup vs v1.0.8:** 2.03x

Feature-type recovery where LiftOn devel ≠ v1.0.8 (gene-like lift):

| feature type | n_reference | v1.0.8 recovered | devel recovered | |
|---|---|---|---|---|
| CDS | 26986 | 1114 | 1877 | |
| exon | 36459 | 13314 | 14008 | |

### yeast_cerevisiae_to_paradoxus — S. cerevisiae -> S. paradoxus  (close_cross_species · _new_taxa_fungi_; RefSeq; n_coding=767)

| tool | completeness_coding | feature_total | mean PI | median PI | %identical | n recovered | wall (s) | peak RSS (MB) | errors/warn |
|---|---|---|---|---|---|---|---|---|---|
| Liftoff | 0.95437 | 0.90224 | 0.89639 | 0.92824 | 0.03415 | 732 | 6.2 | 292 | 85/100 |
| miniprot | 0.98044 | n/a | 0.89896 | 0.92704 | 0.03059 | 752 | 0.6 | 197 | 50/151 |
| LiftOn v1.0.8 | 0.95958 | 0.90695 | 0.9076 | 0.9279 | 0.03397 | 736 | 10.8 | 301 | 10/9 |
| LiftOn devel | 0.96089 | 0.90813 | 0.90769 | 0.92812 | 0.03392 | 737 | 9.6 | 266 | 11/10 |

- **LiftOn devel − v1.0.8:** mean PI +0.00009, completeness +0.00131, n_recovered +1
- **LiftOn devel − best(Liftoff, miniprot):** mean PI +0.00873, completeness -0.01955
- **LiftOn devel speedup vs v1.0.8:** 1.13x

Feature-type recovery where LiftOn devel ≠ v1.0.8 (gene-like lift):

| feature type | n_reference | v1.0.8 recovered | devel recovered | |
|---|---|---|---|---|
| CDS | 808 | 39 | 42 | |
| mRNA | 767 | 736 | 737 | |
| pseudogene | 1 | 0 | 1 | ⬅ gene-like |

## Full-genome headlines

### arabidopsis — Arabidopsis TAIR10 -> ASM2311539v1  (same-species; RefSeq; n_coding=48265)

| tool | completeness_coding | feature_total | mean PI | median PI | %identical | n recovered | wall (s) | peak RSS (MB) | errors/warn |
|---|---|---|---|---|---|---|---|---|---|
| Liftoff | 0.99896 | 0.8813 | 0.99836 | 1 | 0.99355 | 48215 | n/a | n/a | 137/150 |
| miniprot | 0.99917 | n/a | 0.99619 | 1 | 0.86559 | 48225 | n/a | n/a | 50/200 |
| LiftOn v1.0.8 | 0.28215 | 0.25154 | 0.99884 | 1 | 0.9942 | 13618 | n/a | n/a | 149/55 |
| LiftOn devel | 0.9988 | 0.97951 | 0.99904 | 1 | 0.99369 | 48207 | n/a | n/a | 107/114 |

- **LiftOn devel − v1.0.8:** mean PI +0.00020, completeness +0.71665, n_recovered +34589
- **LiftOn devel − best(Liftoff, miniprot):** mean PI +0.00068, completeness -0.00037

Feature-type recovery where LiftOn devel ≠ v1.0.8 (gene-like lift):

| feature type | n_reference | v1.0.8 recovered | devel recovered | |
|---|---|---|---|---|
| CDS | 286264 | 13545 | 284355 | |
| antisense_RNA | 92 | 20 | 92 | |
| exon | 324848 | 90468 | 322744 | |
| gene | 33468 | 9476 | 33175 | |
| lnc_RNA | 3878 | 1197 | 3878 | ⬅ gene-like |
| mRNA | 52177 | 13618 | 52086 | |
| miRNA | 428 | 0 | 1 | ⬅ gene-like |
| ncRNA | 540 | 96 | 314 | ⬅ gene-like |
| primary_transcript | 326 | 93 | 325 | |
| pseudogene | 4851 | 0 | 4816 | ⬅ gene-like |
| rRNA | 14 | 4 | 14 | ⬅ gene-like |
| snRNA | 82 | 23 | 77 | ⬅ gene-like |
| snoRNA | 287 | 91 | 283 | ⬅ gene-like |
| tRNA | 684 | 235 | 663 | ⬅ gene-like |
| three_prime_UTR | 23 | 0 | 1 | |
| transcript | 1826 | 197 | 1822 | |

### bee — Honey bee HAv3.1 -> ASM1932182v1  (same-species; RefSeq; n_coding=23471)

| tool | completeness_coding | feature_total | mean PI | median PI | %identical | n recovered | wall (s) | peak RSS (MB) | errors/warn |
|---|---|---|---|---|---|---|---|---|---|
| Liftoff | 0.99544 | 0.94493 | 0.98071 | 1 | 0.70189 | 23364 | n/a | n/a | 100/150 |
| miniprot | 0.99851 | n/a | 0.98046 | 1 | 0.56584 | 23436 | n/a | n/a | 50/200 |
| LiftOn v1.0.8 | 0.99544 | 0.94493 | 0.99024 | 1 | 0.70283 | 23364 | n/a | n/a | 103/96 |
| LiftOn devel | 0.99544 | 0.94565 | 0.99033 | 1 | 0.70275 | 23364 | n/a | n/a | 53/99 |

- **LiftOn devel − v1.0.8:** mean PI +0.00009, completeness +0.00000, n_recovered +0
- **LiftOn devel − best(Liftoff, miniprot):** mean PI +0.00962, completeness -0.00307

Feature-type recovery where LiftOn devel ≠ v1.0.8 (gene-like lift):

| feature type | n_reference | v1.0.8 recovered | devel recovered | |
|---|---|---|---|---|
| CDS | 226465 | 16432 | 138297 | |
| exon | 267977 | 189742 | 194744 | |
| gene | 12356 | 12296 | 12295 | |
| pseudogene | 42 | 0 | 33 | ⬅ gene-like |
| rRNA | 57 | 56 | 55 | ⬅ gene-like |

### rice — Rice IRGSP -> ASM3414082v1  (same-species; RefSeq; n_coding=42580)

| tool | completeness_coding | feature_total | mean PI | median PI | %identical | n recovered | wall (s) | peak RSS (MB) | errors/warn |
|---|---|---|---|---|---|---|---|---|---|
| Liftoff | 0.99899 | 0.95623 | 0.99643 | 1 | 0.98747 | 42537 | n/a | n/a | 137/154 |
| miniprot | 0.99958 | n/a | 0.99352 | 1 | 0.85666 | 42562 | n/a | n/a | 50/200 |
| LiftOn v1.0.8 | 0.77344 | 0.72901 | 0.99848 | 1 | 0.99022 | 32933 | n/a | n/a | 174/56 |
| LiftOn devel | 0.99876 | 0.97587 | 0.99818 | 1 | 0.98791 | 42527 | 10174 | 8172 | 129/102 |

- **LiftOn devel − v1.0.8:** mean PI -0.00030, completeness +0.22532, n_recovered +9594
- **LiftOn devel − best(Liftoff, miniprot):** mean PI +0.00175, completeness -0.00082

Feature-type recovery where LiftOn devel ≠ v1.0.8 (gene-like lift):

| feature type | n_reference | v1.0.8 recovered | devel recovered | |
|---|---|---|---|---|
| CDS | 238915 | 32617 | 42032 | |
| exon | 334017 | 254006 | 331837 | |
| gene | 33624 | 25583 | 33565 | |
| lnc_RNA | 6277 | 4487 | 6271 | ⬅ gene-like |
| mRNA | 42580 | 32933 | 42527 | |
| pseudogene | 1605 | 0 | 1593 | ⬅ gene-like |
| pseudogenic_tRNA | 3 | 0 | 2 | ⬅ gene-like |
| rRNA | 200 | 9 | 188 | ⬅ gene-like |
| snRNA | 71 | 59 | 71 | ⬅ gene-like |
| snoRNA | 610 | 482 | 610 | ⬅ gene-like |
| tRNA | 672 | 476 | 666 | ⬅ gene-like |
| transcript | 3840 | 2703 | 3836 | |

### drosophila — D. melanogaster -> D. erecta  (cross_species; RefSeq; n_coding=30799)

| tool | completeness_coding | feature_total | mean PI | median PI | %identical | n recovered | wall (s) | peak RSS (MB) | errors/warn |
|---|---|---|---|---|---|---|---|---|---|
| Liftoff | 0.97075 | 0.81686 | 0.88695 | 0.95397 | 0.04477 | 29898 | n/a | n/a | 131/181 |
| miniprot | 0.99114 | n/a | 0.91315 | 0.95441 | 0.03882 | 30526 | n/a | n/a | 50/200 |
| LiftOn v1.0.8 | 0.97403 | 0.82015 | 0.92106 | 0.95588 | 0.04591 | 29999 | 2755.1 | 19659 | 150/154 |
| LiftOn devel | 0.97084 | 0.81939 | 0.9224 | 0.95605 | 0.04587 | 29901 | 5945 | 7105 | 100/135 |

- **LiftOn devel − v1.0.8:** mean PI +0.00134, completeness -0.00319, n_recovered -98
- **LiftOn devel − best(Liftoff, miniprot):** mean PI +0.00925, completeness -0.02030
- **LiftOn devel speedup vs v1.0.8:** 0.46x

Feature-type recovery where LiftOn devel ≠ v1.0.8 (gene-like lift):

| feature type | n_reference | v1.0.8 recovered | devel recovered | |
|---|---|---|---|---|
| CDS | 163306 | 6413 | 7221 | |
| exon | 190719 | 82817 | 83552 | |
| gene | 17559 | 16458 | 16372 | |
| mRNA | 30799 | 29999 | 29901 | |
| pseudogene | 339 | 0 | 120 | ⬅ gene-like |
| rRNA | 134 | 66 | 83 | ⬅ gene-like |
| tRNA | 317 | 281 | 282 | ⬅ gene-like |

