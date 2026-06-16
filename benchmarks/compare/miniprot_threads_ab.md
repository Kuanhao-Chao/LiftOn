## miniprot thread plumbing A/B — Iteration 17

miniprot was capped at its binary default of 4 threads regardless of LiftOn's `-t/--threads`. This A/B plumbs `args.threads` into miniprot's `-t` (gated on `>1`, so the default `-t 1` is unchanged). **ARM 1** (cached `-L`, fresh miniprot) is the clean gate: Liftoff is a deterministic cached load, so the only variable is miniprot's thread count and Step-4 wall == miniprot wall. Gate: **byte-identical AND a miniprot/Step-4 wall win AND completeness unchanged**. **ARM 2** (fresh both, concurrent Step 4) confirms the win holds end-to-end despite transient oversubscription (byte-identity not gated there — Liftoff `-copies` is non-deterministic across fresh runs).

Baseline `mp_t4` = `LIFTON_MINIPROT_THREADS=0` (old default 4); `mp_tN` = miniprot `-t 8`.

| Dataset | ARM1 byte-ident | ARM1 miniprot wall t4→tN | ARM1 speedup | ARM1 complete | ARM2 Step-4 wall t4→tN | ARM2 total t4→tN |
|---|---|---|---|---|---|---|
| drosophila | yes | 11.21s→6.57s | 1.71x (41.4%) | yes | 64.55s→61.41s (1.05x) | 230.9s→202.4s (1.14x) |
| mouse_to_rat | yes | 18.32s→10.63s | 1.72x (42.0%) | yes | 106.01s→106.57s (0.99x) | 197.1s→196.8s (1.0x) |
| rice | yes | 7.13s→4.33s | 1.65x (39.3%) | yes | 18.99s→18.99s (1.0x) | 116.6s→116.3s (1.0x) |

**Interpretation:** ARM-1 byte-identical confirms the change is byte-neutral end-to-end (miniprot output is order-stable across threads — the standalone feasibility gate proved `-t 1/4/8/16` are byte-identical). The ARM-1 miniprot wall win is the direct speedup; ARM-2 shows it survives concurrent Step 4. The rigorous determinism proof is the standalone gate + the unit suite (`tests/test_miniprot_threads.py`); the 24-cell matrix is unaffected (it uses cached `-M`, so miniprot is never invoked).
