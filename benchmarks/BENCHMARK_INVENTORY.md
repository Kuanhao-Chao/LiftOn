# Curated Benchmark Inventory

This inventory is generated deterministically from Git-tracked benchmark files and explicitly managed provenance files. Local run trees are excluded from release evidence.

- Curated files: **294** (13,543,609 bytes)
- Registered benchmark IDs: **35**
- Frozen artifacts verified: **1**
- Inventory digest: `1933fbad4859c4e80a471d1e47898ff81c8d695a46ea6fceb27e3b8ed77dfd60`

## Classifications

| Family | Class | State | Claim eligible | Files | Superseded by |
| --- | --- | --- | --- | ---: | --- |
| `canonical_release_baseline` | canonical | complete | yes | 12 | — |
| `canonical_registries` | canonical | supporting | no | 5 | — |
| `canonical_release_tooling` | canonical | supporting | no | 35 | — |
| `canonical_v2_provenance` | canonical | supporting | no | 16 | — |
| `canonical_tier_full_results` | canonical | complete | yes | 11 | — |
| `diagnostic_headroom` | diagnostic | complete | no | 33 | — |
| `historical_ab_experiments` | historical | superseded | no | 75 | canonical_release_baseline |
| `historical_v1_review` | historical | superseded | no | 79 | canonical_release_baseline |
| `obsolete_legacy_runners` | obsolete | retired | no | 17 | canonical_release_tooling |
| `ineligible_incomplete_soybean` | ineligible | incomplete | no | 1 | — |
| `promotion_ab_v1012` | diagnostic | complete | no | 9 | — |
| `promotion_ab_v1014` | diagnostic | complete | no | 1 | — |

## Canonical-v2 Expansion

The expansion fixes **18** source packages (15 remote and 3 existing-registry packages) across **12** approved scenario families.

| Scenario | Kind | Panels |
| --- | --- | --- |
| `v2_truth_human_grch38_chm13` | biological | subset, full, truth |
| `v2_dialect_ensembl116_gtf` | biological | subset, full, truth |
| `v2_dialect_flybase_dmel_dere` | biological | subset, full, cross_e2e, truth |
| `v2_dialect_wormbase_celegans_cbriggsae` | biological | subset, full, cross_e2e, truth |
| `v2_truth_soybean_w82_lee` | biological | subset, full, truth |
| `v2_deep_zebrafish_xenopus` | biological | subset, full, truth |
| `v2_deep_tomato_rice` | biological | subset, full, truth |
| `v2_truth_rat_mouse_e2e` | biological | subset, cross_e2e, truth |
| `v2_synth_chr22_fragmented` | synthetic | synthetic |
| `v2_synth_chr22_sv` | synthetic | synthetic |
| `v2_protocol_thread_scaling_bee` | protocol | protocol |
| `v2_protocol_io_modes_arabidopsis` | protocol | protocol |

## Evidence Policy

- `canonical` files may support a claim only when the family is also marked claim-eligible.
- `diagnostic` and `historical` results explain design decisions but must not be pooled into the v1.0.10 release campaign.
- `obsolete` files are retained solely for reproducibility of older reviews.
- `ineligible` files are incomplete, mutable, local, or otherwise unsuitable as publication evidence.
- `benchmarks/compare/fourway_results.json` is byte-frozen; the scanner fails if its SHA-256 changes.

## Regeneration

```bash
python -m benchmarks.inventory --check
python -m benchmarks.inventory --stdout json > /tmp/inventory.json
python -m benchmarks.inventory --stdout markdown > /tmp/inventory.md
```

Review and apply regenerated files deliberately; do not mix local `_runs/`, `work/`, rerun figures, or untracked outputs into this index.
