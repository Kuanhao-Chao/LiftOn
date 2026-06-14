## Fresh parallel Step 7 without `--native` — Iteration 8

Each state reuses cached `-L`/`-M` (Step 4 skipped → Step-7 dispatch isolated, byte-identity clean) and runs at the same `-t`, differing only by `--locus-pipeline` / `--native`. The change routes parallel Step 7 through the materialise + proxy-DB path on **any** backend, so `--locus-pipeline` now fans out WITHOUT `--native`. Gate: **byte-identical (hard) AND a Step-7 wall win**.

Threads: `-t 8`. **Step-7 wall** (probe `LIFTON_PERF_STEP7`) is the direct measure; **total wall** is end-to-end.

| Dataset | byte-ident serial=parallel | native parity | Step-7 wall serial→parallel | Step-7 speedup | total wall serial→parallel | total speedup | complete |
|---|---|---|---|---|---|---|---|
| drosophila | yes | yes | 174.33s→113.61s | 1.53× (34.8%) | 211.6s→150.7s | 1.4× (28.8%) | yes |
| mouse_to_rat | yes | yes | 106.56s→82.65s | 1.29× (22.4%) | 127.6s→103.9s | 1.23× (18.6%) | yes |
| rice | yes | yes | 98.6s→73.24s | 1.35× (25.7%) | 133.4s→107.4s | 1.24× (19.5%) | yes |

**Interpretation:** byte-identical across all three states is the hard gate (the change is pure scheduling). `parallel ≈ parallel_native` Step-7 wall confirms the speedup no longer requires `--native`. The rigorous byte-identity proof is the unit suite (`tests/test_fresh_parallel_step7.py` + the now-non-native-exercising 24-cell matrix `tests/test_native_matrix.py`).
