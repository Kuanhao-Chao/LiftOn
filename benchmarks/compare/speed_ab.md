## Speed: serial vs threaded (byte-identical)

LiftOn run twice on the same cached `-L`/`-M` inputs (genes mode, non-`--optimize`): **serial** `-t 1` (no fast path) vs **threaded** `-t 8 --native --locus-pipeline` (Phase-17b/17c materialised). `Byte-identical` MUST be ✓ — the harness-side echo of the 24-cell matrix gate.

| Dataset | Serial wall (s) | Threaded wall (s) | Speedup | Serial RSS (MB) | Threaded RSS (MB) | Byte-identical |
|---|---|---|---|---|---|---|
| human_mane | 24.35 | 13.77 | 1.77× | 1516.0 | 2280.0 | ✓ |
| drosophila | 262.04 | 137.69 | 1.9× | 2520.0 | 4798.0 | ✓ |
| arabidopsis | 244.05 | 231.29 | 1.06× | 1558.0 | 5057.0 | ✓ |
| rice | 135.85 | 111.44 | 1.22× | 1474.0 | 3162.0 | ✓ |
