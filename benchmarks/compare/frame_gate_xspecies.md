## Frame-consistency gate: mammalian cross-species 3-state A/B (Step 7 merge)

Self-contained — all three states run fresh on the same cached `-L`/`-M` (genes mode, transcript-space miniprot), scored by the independent re-alignment evaluator. **state 1** = default (unconditional merge); **state 2** = `--optimize` best-of-outcome alone (`LIFTON_MERGE_FRAME_GATE=0`); **state 3** = best-of-outcome + frame gate (`=1`). Acceptance: state3 ≥ state2 ≥ state1 mean protein identity; completeness not reduced.

| Dataset | s1 mean PI | s2 mean PI | s3 mean PI | s2 vs s1 (impr/regr, net) | s3 vs s2 (impr/regr, net) | merge-kept s1→s2→s3 | complete s1→s2→s3 | wall s1/s2/s3 |
|---|---|---|---|---|---|---|---|---|
| mouse_to_rat | 0.91109 | 0.91778 | 0.91719 | 294/2, +0.00679 | 0/23, -0.00059 | 6848→6550→6545 | 359→359→359 | 1306.74/1272.63/1289.97 |
