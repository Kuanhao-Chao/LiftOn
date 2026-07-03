## candidate-3 divergence-ladder A/B (Phase 1)

`off` = `LIFTON_MINIPROT_CANDIDATE=0` (2-way merge); `on` = default (3-way, miniprot-only candidate). Cached `-L`/`-M`, `-t1` no-`copies`, scored by the independent re-alignment evaluator (keyed by ref_mrna_id). candidate-3 SWAPS a transcript's model (never adds/drops), so **n_lost** and **common_regressed** MUST be 0 (the additive guarantee, re-verified neutrally); **common_improved** is the win; **a2a vs miniprot** is the mean(LiftOn − miniprot) deficit, which should shrink. **validity** is the structured `gff3_validator` error count (the `Errors : N` summary line). Gate: n_lost=0 AND 0 regressed AND validity not worse AND deficit not worse AND (improved>0 OR the cell is inert — a legitimate no-op).

| Dataset | divergence | improved | regressed | n_lost | tagged | mean PI off→on | a2a off→on | val off→on | gate |
|---|---|---|---|---|---|---|---|---|---|
| celegans_to_briggsae | distant_cross_species | 25 | 0 | 0 | 30 | 0.72086→0.7215 | -0.02614→-0.02549 | 0→0 | **PASS** |
| drosophila_to_anopheles | very_distant_cross_species | 2 | 0 | 0 | 2 | 0.56333→0.56335 | -0.01982→-0.01979 | 0→0 | **PASS** |
| zebrafish_to_medaka | distant_cross_species | 1 | 0 | 0 | 1 | 0.58714→0.58716 | -0.07522→-0.0752 | 0→0 | **PASS** |
| rice_to_sorghum | distant_cross_species | 39 | 0 | 0 | 42 | 0.67517→0.67745 | -0.04364→-0.04135 | 0→0 | **PASS** |
| t4_human_to_xenopus | very_distant_cross_species | 0 | 0 | 0 | 0 | 0.63453→0.63453 | -0.05439→-0.05439 | 0→0 | **PASS** |
| t4_human_to_chicken | very_distant_cross_species | 20 | 0 | 0 | 20 | 0.72046→0.72146 | -0.02743→-0.02643 | 0→0 | **PASS** |
| human_to_mouse | distant_cross_species | 36 | 0 | 0 | 36 | 0.78091→0.78266 | -0.03209→-0.03034 | 0→0 | **PASS** |
| drosophila | D. melanogaster -> D. erecta | 2 | 0 | 0 | 2 | 0.92573→0.92576 | 0.00267→0.00271 | 0→0 | **PASS** |
