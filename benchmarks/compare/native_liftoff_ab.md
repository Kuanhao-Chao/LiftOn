## --native vs subprocess Liftoff A/B — Iteration 7 (fresh-Liftoff fix)

Each state runs LiftOn **fresh** (no cached `-L`/`-M`). **subprocess** = Liftoff drives minimap2 as a subprocess (baseline); **native** = `--native` drives minimap2 via mappy in-process. Before the `mm2_options`/`--eqx` fix, native mapped **0** genes; the gate is that it now maps a comparable gene set at ≈equal protein identity. Byte-identity is NOT expected (mappy ≠ subprocess minimap2; Liftoff `-copies` is non-deterministic across runs).

Threads: `-t 8`.

| Dataset | genes sub→native (ratio) | mRNA sub→native | CDS sub→native | shared mRNA | mean protein-id Δ | improved/regressed | wall sub→native |
|---|---|---|---|---|---|---|---|
| bee | 1544→1544 (1.0) | 3680→3680 | 27495→27420 | 3680 | -0.001963 | 5/27 | 138.1s→160.0s |
| drosophila | 4049→4058 (1.0022) | 7947→7943 | 37305→36434 | 7925 | -0.015951 | 44/477 | 279.9s→358.8s |
| mouse_to_rat | 945→950 (1.0053) | 3114→3119 | 23633→23632 | 3098 | -0.000165 | 0/2 | 242.6s→543.1s |
