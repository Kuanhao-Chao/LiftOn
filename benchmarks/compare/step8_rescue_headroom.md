## Step-8 miniprot-rescue headroom diagnostic (Iteration 13, read-only)

Default run on cached `-L`/`-M`, `LIFTON_STEP8_DIAG` set (byte-frozen output). Per candidate miniprot mRNA: drop reason + tunable-filter metrics. **Near-miss** = drops just past the threshold (overlap in (0.1,0.2]; length in [0.8,0.9) ∪ [1.5,1.6)) — the UPPER BOUND on what relaxing could recover. GO if near-miss/emitted ≥ 0.05 on any dataset; else NO-GO (defaults well-tuned).

| Dataset | candidates | emitted | overlap-drop | length-drop | overlap near-miss | length near-miss | near-miss/emitted | verdict |
|---|---|---|---|---|---|---|---|---|
| drosophila | 7253 | 22 | 7222 | 5 | 3 | 2 | 0.2273 | GO (worth an A/B) |
| mouse_to_rat | 2225 | 0 | 2165 | 58 | 0 | 0 | 0.0 | NO-GO / negligible headroom |
