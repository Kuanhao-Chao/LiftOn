## Regime-gated miniprot-only rescue A/B — Iteration 22

`rescue_<floor>` = `--miniprot-rescue` with PI floor `<floor>`; cached `-L`/`-M`, scored by the independent re-alignment evaluator (keyed by ref_mrna_id). **n_added** = coding ref transcripts rescue emits that the DNA lift missed (genuine completeness); **n_redundant** = rescued models that did NOT add a new ref id (the Iter-15 duplicate failure mode — MUST be 0); **n_lost** = off ⊆ on gate (MUST be 0). Gate: n_added>0 AND n_lost=0 AND 0 regressed AND n_redundant=0 AND mean PI added ≥ floor AND validity not worse.

| Dataset | floor | n_added | mean PI added | Δcompl | n_lost | regr | tagged | redundant | val d→r | gate |
|---|---|---|---|---|---|---|---|---|---|---|
| drosophila | 0.5 | 1 | 0.55952 | 0.00014 | 0 | 0 | 1 | 0 | 3→3 | **PASS** |
| human_to_mouse | 0.5 | 30 | 0.71014 | 0.00972 | 0 | 0 | 30 | 0 | 3→3 | **PASS** |
| celegans_to_briggsae | 0.5 | 555 | 0.72775 | 0.07351 | 0 | 0 | 555 | 0 | 3→3 | **PASS** |
| rice_to_sorghum | 0.5 | 100 | 0.66174 | 0.01709 | 0 | 0 | 100 | 0 | 3→3 | **PASS** |
| drosophila_to_anopheles | 0.5 | 201 | 0.64026 | 0.02772 | 0 | 0 | 201 | 0 | 3→3 | **PASS** |
| zebrafish_to_medaka | 0.5 | 143 | 0.69431 | 0.04812 | 0 | 0 | 143 | 0 | 2→2 | **PASS** |
| t4_human_to_chicken | 0.5 | 51 | 0.71979 | 0.01653 | 0 | 0 | 51 | 0 | 3→3 | **PASS** |
| t4_human_to_xenopus | 0.5 | 87 | 0.71272 | 0.02819 | 0 | 0 | 87 | 0 | 2→2 | **PASS** |
