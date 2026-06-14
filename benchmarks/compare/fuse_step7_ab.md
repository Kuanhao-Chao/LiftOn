## Fused Step-7 pool vs two-phase A/B — Iteration 10

Both states run on cached `-L`/`-M` (Step 4 skipped → Step-7 dispatch isolated, byte-identity clean) with the SAME argv (`-t N --locus-pipeline`), differing only by the `LIFTON_FUSE_STEP7` env gate. **two_phase** = `LIFTON_FUSE_STEP7=0` (pre-Iteration-10 prefetcher pool + barrier); **fused** = default (materialise + process in one pool). The change is pure scheduling. Gate: **byte-identical (hard) AND a Step-7 wall win**.

Threads: `-t 8`. **Step-7 wall** (probe `LIFTON_PERF_STEP7`) is the direct measure; **total wall** is end-to-end.

| Dataset | byte-ident | Step-7 wall two_phase→fused | Step-7 speedup | total wall two_phase→fused | total speedup | peak RSS Δ | complete |
|---|---|---|---|---|---|---|---|
| drosophila | yes | 107.86s→102.49s | 1.05× (5.0%) | 144.0s→138.4s | 1.04× (3.9%) | -43.8% | yes |
| mouse_to_rat | yes | 79.23s→67.09s | 1.18× (15.3%) | 99.6s→87.2s | 1.14× (12.5%) | -25.5% | yes |
| rice | yes | 71.79s→64.03s | 1.12× (10.8%) | 105.6s→98.3s | 1.07× (6.9%) | -37.8% | yes |

**Interpretation:** byte-identical is the hard gate (the fusion is pure scheduling). A positive Step-7 speedup confirms the SQLite-bound materialise now overlaps the GIL-released parasail instead of running in a separate barrier'd phase. The rigorous byte-identity proof is the unit suite (`tests/test_fuse_step7.py` + the 24-cell matrix `tests/test_native_matrix.py`, whose on-disk cells take the fused path).
