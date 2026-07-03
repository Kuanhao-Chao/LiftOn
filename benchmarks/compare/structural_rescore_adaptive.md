## Phase-3 structural re-score — adaptive (fixed vs adaptive)

Re-evaluates the existing A/B output GFF3s with the coordinate-independent structural metrics (ORF-validity, intron-chain exactness, exon Sn·Sp). **COMMON Δ** (changed − base on ref transcripts recovered by BOTH arms) is the real regression check — it must be >= −0.005. **ADDED** = the mean structure of the newly-recovered genes (informational: expected weaker — that is the honest recall/precision trade, NOT a regression).

| Dataset | common | added | COMMON Δ orf/intron/sn/sp | ADDED orf/intron/sn/sp | verdict |
|---|---|---|---|---|---|
| drosophila_to_anopheles | 1431 | 73 | 0.0/0.0/0.0/0.0 | 0.09589/0.09589/0.16712/0.15183 | **OK** |
| zebrafish_to_medaka | 690 | 23 | 0.0/0.0/0.0/0.0 | 0.04348/0.04348/0.0942/0.1087 | **OK** |
| t4_human_to_xenopus | 331 | 15 | 0.0/0.0/0.0/0.0 | 0.2/0.13333/0.18556/0.18794 | **OK** |
| celegans_to_briggsae | 3113 | 18 | 0.0/0.0/0.0/0.0 | 0.16667/0.05556/0.21667/0.28148 | **OK** |
