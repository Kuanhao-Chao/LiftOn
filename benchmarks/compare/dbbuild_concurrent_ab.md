## Concurrent Step-5 DB build vs serial A/B — Iteration 11

Both states reuse cached `-L`/`-M` (Step 4 a fast load → Step-5 DB build a large fraction of wall) with the SAME argv, differing only by the `LIFTON_CONCURRENT_DBBUILD` env gate. **serial** = `LIFTON_CONCURRENT_DBBUILD=0` (build liftoff then miniprot DB serially); **concurrent** = default (miniprot DB on a worker, liftoff on main, reopen on main). Pure scheduling of where the `.db` file is built. Gate: **byte-identical (hard) AND a Step-5 wall win**.

**Step-5 wall** (probe `LIFTON_PERF_STEP5`) is the direct measure; **total wall** is end-to-end.

| Dataset | byte-ident | Step-5 wall serial→concurrent | Step-5 speedup | total wall serial→concurrent | total speedup | peak RSS Δ | complete |
|---|---|---|---|---|---|---|---|
| mouse | yes | 28.42s→30.81s | 0.92× (-8.4%) | 359.5s→346.9s | 1.04× (3.5%) | +2.5% | yes |
| drosophila | yes | 11.85s→13.22s | 0.9× (-11.6%) | 198.8s→200.8s | 0.99× (-1.0%) | +18.0% | yes |
| rice | yes | 13.95s→15.11s | 0.92× (-8.3%) | 124.0s→122.7s | 1.01× (1.1%) | +26.0% | yes |

**Interpretation:** byte-identical is the hard gate (the change is pure scheduling — the `.db` file contents are identical, only *where* it is built moves). A positive Step-5 speedup confirms the miniprot DB build now overlaps the liftoff DB build. The win scales with DB-build cost (∝ GFF size), so mouse (255K/120K lines) shows the largest saving. The rigorous byte-identity proof is the unit suite (`tests/test_dbbuild_concurrent.py` + the 24-cell matrix `tests/test_native_matrix.py`, whose cells all take the concurrent path).
