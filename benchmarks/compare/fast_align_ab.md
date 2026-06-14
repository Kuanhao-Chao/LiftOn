## Fast-align (band everything) vs default A/B — Iteration 3

Same cached `-L`/`-M`, single-threaded, differing only by `LIFTON_FAST_ALIGN`. Promotion gate: **speedup ≥ 5% wall AND mean identity Δ ≥ −1e-3 AND completeness unchanged**.

| Dataset | wall default→fast | Δwall | speedup | meanΔ id | improved | regressed | worst id Δ | complete | PROMOTE |
|---|---|---|---|---|---|---|---|---|---|
| mouse_to_rat | 234.5s→109.6s | +53.3% | 2.14× | +0.000000 | 0 | 0 | ('rna-NM_001001181.4', 0.0) | yes | yes |
| drosophila | 273.8s→191.8s | +30.0% | 1.428× | -0.000003 | 3 | 3 | ('rna-NM_141762.6', -0.012) | yes | no |
| mouse | 855.4s→330.1s | +61.4% | 2.591× | +0.000000 | 0 | 0 | ('rna-NM_001001178.2', 0.0) | yes | yes |
