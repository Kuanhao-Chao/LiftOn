## Chaining headroom diagnostic (Iteration 4, read-only)

Default best-of-outcome run on cached `-L`/`-M`; per-transcript chain log parsed. **#2 empty-drop** = transcripts with a both-zero chunk whose CDS is currently dropped (headroom for the Liftoff-fallback refinement). **#1 coarseness** = single-chunk fraction (headroom for relaxing the exact-genomic-end sync). Decision: material #2 -> implement #2; mostly coarse -> consider #1 (higher risk); neither -> chaining near ceiling.

| Dataset | merge-fired | kept/liftoff (valve) | #2 empty-drop txpts (%) | dropped chunks | #1 single-chunk (%) | mean chunks |
|---|---|---|---|---|---|---|
| drosophila | 6699 | 6594/1330 | 5 (0.1%) | 5 | 2215 (33.1%) | 3.98 |
| mouse_to_rat | 1983 | 1894/1220 | 19 (1.0%) | 30 | 338 (17.0%) | 9.05 |
