## CDS attribute parity divergence-ladder A/B

`off` = `LIFTON_NO_CDS_ATTR_CARRY=1` (the pre-change `{Parent}`-only CDS rows); `on` = default (the reference CDS's descriptive attributes carried through the chaining merge and the ORF-rescue rebuild). Cached `-L`/`-M`, `-t1` no-`copies`, identical argv — only the env differs.

The change touches column 9 only, so the gate is stronger than a mean-PI comparison: **columns 1-8 must be byte-identical** on every row, which is exactly a CDS's coordinates, strand and phase — proving the translated protein cannot have moved. On top of that the ON attribute set must be a strict superset of OFF (nothing lost, no shared value rewritten), LiftOn's own `protein_identity` must be unchanged on every transcript, and `gff3-validate` must not report more errors.

| Dataset | divergence | rows | cols1-8 Δ | attrs lost | attrs changed | CDS coverage off→on | PI changed | validity off→on | size Δ | gate |
|---|---|---|---|---|---|---|---|---|---|---|
| celegans_to_briggsae | distant_cross_species | 48313 | 0 | 0 | 0 | 0.41132→1.0 | 0 | 0→0 | 30.56% | PASS |
| drosophila_to_anopheles | very_distant_cross_species | 18539 | 0 | 0 | 0 | 0.54315→1.0 | 0 | 0→0 | 19.52% | PASS |
| zebrafish_to_medaka | distant_cross_species | 15294 | 0 | 0 | 0 | 0.62728→1.0 | 0 | 0→0 | 12.47% | PASS |
| rice_to_sorghum | distant_cross_species | 67458 | 0 | 0 | 0 | 0.14728→1.0 | 0 | 0→0 | 27.78% | PASS |
| t4_human_to_xenopus | very_distant_cross_species | 5773 | 0 | 0 | 0 | 0.61456→1.0 | 0 | 0→0 | 19.48% | PASS |
| t4_human_to_chicken | very_distant_cross_species | 18041 | 0 | 0 | 0 | 0.13451→1.0 | 0 | 0→0 | 43.02% | PASS |
| human_to_mouse | distant_cross_species | 63228 | 0 | 0 | 0 | 0.08602→1.0 | 0 | 0→0 | 39.9% | PASS |
| drosophila | D. melanogaster -> D. erecta | 92679 | 0 | 0 | 0 | 0.05996→1.0 | 0 | 0→0 | 30.21% | PASS |

**GATE: 8/8 cells pass.**
