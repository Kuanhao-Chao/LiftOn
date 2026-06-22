## Frame-consistency gate: 3-state A/B (Step 7 merge)

All three states on the same cached `-L`/`-M` (genes mode), scored by the independent re-alignment evaluator. **state 1** = default (unconditional merge); **state 2** = `--optimize` best-of-outcome alone (`LIFTON_MERGE_FRAME_GATE=0`); **state 3** = best-of-outcome + frame gate (`=1`). Acceptance: state3 ≥ state2 ≥ state1 mean protein identity; more merges kept (fewer fallbacks) under state 3.

| Dataset | s1 mean PI | s2 mean PI | s3 mean PI | s3 vs s2 (impr/regr, net) | merge-kept s2→s3 | fallback s2→s3 |
|---|---|---|---|---|---|---|
| drosophila | 0.92457 | 0.92578 | 0.92578 | 0/0, +0.00000 | 6594→6593 | 1330→1331 |
