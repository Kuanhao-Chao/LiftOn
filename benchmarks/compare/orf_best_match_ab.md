## --orf-best-match A/B (ORF-rescue best-matching ORF) — Iteration 9

Both states run on the same cached `-L`/`-M` (Step 4 skipped → only the ORF-rescue selection differs), scored by the independent re-alignment evaluator. **default** = legacy longest-per-frame; **orf_best** = `LIFTON_ORF_BEST_MATCH=1` top-K best-match. Gate: mean-PI Δ ≥ +0.003 on BOTH datasets AND completeness identical AND regressed ≤ 0.1×improved.

| Dataset | default mean PI | orf_best mean PI | mean Δ | improved/regressed | complete | wall d/o | gate |
|---|---|---|---|---|---|---|---|
| drosophila | 0.92578 | 0.92575 | -0.000030 | 4/59 | **NO** | 146.62/141.12 | NO-GO |
| mouse_to_rat | 0.9092 | 0.90928 | +0.000080 | 3/23 | **NO** | 99.12/93.49 | NO-GO |
