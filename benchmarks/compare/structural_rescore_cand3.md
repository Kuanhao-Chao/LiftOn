## Phase-3 structural re-score — cand3 (off vs on)

Re-evaluates the existing A/B output GFF3s with the coordinate-independent structural metrics (ORF-validity, intron-chain exactness, exon Sn·Sp). **COMMON Δ** (changed − base on ref transcripts recovered by BOTH arms) is the real regression check — it must be >= −0.005. **ADDED** = the mean structure of the newly-recovered genes (informational: expected weaker — that is the honest recall/precision trade, NOT a regression).

| Dataset | common | added | COMMON Δ orf/intron/sn/sp | ADDED orf/intron/sn/sp | verdict |
|---|---|---|---|---|---|
| t4_human_to_chicken | 822 | 0 | -0.00122/0.0/0.00021/9e-05 | None/None/None/None | **OK** |
| human_to_mouse | 2657 | 0 | 0.00753/0.00489/0.00562/0.00025 | None/None/None/None | **OK** |
| rice_to_sorghum | 4556 | 0 | 0.00307/0.00022/0.0007/-0.00064 | None/None/None/None | **OK** |
| celegans_to_briggsae | 3113 | 0 | -0.00032/0.00064/0.00051/-0.00011 | None/None/None/None | **OK** |
