## Step-8 rescue: default vs relaxed thresholds A/B — Iteration 13

`relaxed` = `-overlap 0.2 -min_miniprot 0.8 -max_miniprot 1.6` (values from the headroom near-miss buckets); cached `-L`/`-M`, scored by the independent re-alignment evaluator. **n_added** = transcripts the relaxed filters newly emit; **mean PI added** = their protein identity (legit rescues high, garbage low — the false-rescue discriminator). Gate: n_added>0 AND mean PI added ≥ 0.75 AND 0 regressed AND no new validity errors, on BOTH datasets.

| Dataset | n_added | mean PI added | frac ≥0.7 | added PIs | common impr/regr | validity def/relax | gate |
|---|---|---|---|---|---|---|---|
| drosophila | 2 | 0.34599 | 0.0 | [0.132, 0.56] | 1/0 | 3/3 | NO-GO |
| mouse_to_rat | 0 | None | None | [] | 0/0 | 2/2 | NO-GO |
