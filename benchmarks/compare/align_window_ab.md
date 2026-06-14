## Giant-gene windowed alignment: windowed (default) vs forced full-DP A/B

Same cached `-L`/`-M`, single-threaded, differing only by `LIFTON_ALIGN_WINDOW_{AA,NT}`. Acceptance: large speed + RSS reduction; **0 non-giant transcripts differ** (normal genes byte-identical); per-giant protein_identity delta ~0.

| Dataset | wall windowed/fulldp (×) | peak RSS windowed/fulldp (×) | giants | non-giant diffs | max giant PI Δ |
|---|---|---|---|---|---|
| mouse | 799.0s / 1783.0s (2.2×) | 2649.0MB / 44860.0MB (16.9×) | 26 | 0 | ('rna-XM_036160536.1', 0.0) |
