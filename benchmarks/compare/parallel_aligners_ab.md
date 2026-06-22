## --parallel-aligners (concurrent Step 4) vs serial A/B — Iteration 6

Each state runs LiftOn **fresh** (no cached `-L`/`-M`, so Step 4 actually runs both aligners), differing only by `--parallel-aligners`. The flag overlaps Liftoff and miniprot → wall = max() instead of sum. Output is byte-identical (speed-only change). Gate: **byte-identical AND a wall-clock win**.

Threads: `-t 8`. **Step-4 wall** is the direct measure of the overlap (probe via `LIFTON_PERF_STEP4`); **total wall** is end-to-end and is diluted by serial Step 7 on these subsets (no `--locus-pipeline`).

| Dataset | byte-identical | Step-4 wall serial→parallel | Step-4 saved | total wall serial→parallel | total saved | peak RSS Δ | complete |
|---|---|---|---|---|---|---|---|
| bee | yes | 26.7s→20.88s | 5.8s (21.8%) | 136.7s→131.9s | 4.8s (+3.5%) | -0.0% | yes |
| drosophila | **NO** | 73.71s→64.12s | 9.6s (13.0%) | 276.5s→269.6s | 6.9s (+2.5%) | +0.0% | NO |
| mouse_to_rat | yes | 125.58s→106.17s | 19.4s (15.5%) | 256.2s→235.2s | 21.0s (+8.2%) | -0.4% | yes |

**Caveat on `byte-identical=NO`:** the flag itself changes no output (Liftoff runs on the main thread identically in both states; miniprot output is byte-identical). But Liftoff's `-copies` multi-copy alignment is **non-deterministic across fresh runs** (confirmed: two identical-config standalone Liftoff runs on drosophila diverge in tRNA copy lines) — so a fresh serial-vs-parallel diff on a repeat-rich genome is Liftoff noise, not the flag. The rigorous byte-identity proof is the unit test on fixed `-L`/`-M` inputs (`tests/test_parallel_aligners.py`).
