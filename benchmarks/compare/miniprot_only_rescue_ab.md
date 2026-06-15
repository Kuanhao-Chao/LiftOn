## Weak-Liftoff miniprot-rescue A/B — Iteration 15 (Gap #2)

`weak_rescue` = `--miniprot-rescue-weak-liftoff`; cached `-L`/`-M`, scored by the independent re-alignment evaluator (keyed by ref_mrna_id). **n_added** = coding ref transcripts weak_rescue emits but default does not (genuine completeness); **mean PI added** = their identity. **n_lost** = the off ⊆ on gate (must be 0). **extra gene models** = raw gene-line delta (>> n_added ⇒ duplicate/overlapping models, a quality concern). Gate: n_added>0 AND mean PI added ≥ 0.75 AND n_lost=0 AND 0 regressed AND no new validity errors.

| Dataset | n_added | mean PI added | frac ≥0.7 | n_lost | common impr/regr | extra gene models | validity def/weak | gate |
|---|---|---|---|---|---|---|---|---|
| drosophila | 5 | 0.90508 | 1.0 | 0 | 36/0 | 44 | 3/3 | **PASS** |
| mouse_to_rat | 0 | None | None | 0 | 8/0 | 8 | 2/2 | NO-GO |
