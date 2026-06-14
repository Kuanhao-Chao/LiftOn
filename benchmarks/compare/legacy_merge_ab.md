## Promoted merge: default best-of-outcome vs `--legacy-merge` A/B (Step 7)

Self-contained — both states run fresh on the same cached `-L`/`-M` (genes mode, transcript-space miniprot), scored by the independent re-alignment evaluator. **default** = the promoted best-of-outcome merge (no flag); **legacy** = the pre-promotion unconditional merge (`--legacy-merge`). Acceptance: default >= legacy mean protein identity with a favorable improved/regressed ratio; completeness not reduced.

| Dataset | default mean PI | legacy mean PI | default vs legacy (impr/regr, net) | merge-kept default/legacy | complete default/legacy | wall default/legacy |
|---|---|---|---|---|---|---|
| drosophila | 0.92578 | 0.92208 | 115/1, +0.003741 | 6594/6699 | 339/339 | 145.69/116.4 |
| mouse_to_rat | 0.9092 | 0.90398 | 91/0, +0.005221 | 1894/1983 | 121/121 | 126.47/110.25 |
