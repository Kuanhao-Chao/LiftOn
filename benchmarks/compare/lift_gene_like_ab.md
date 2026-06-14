## --lift-gene-like vs default A/B (Iteration 5)

Both states share ONE gene-like Liftoff `-L` (built via the standalone liftoff binary) + cached `-M`. **off** = default (LiftOn processes `gene` only); **on** = `--lift-gene-like` (processes auto-detected gene-like types). Same `-L` → the gene lift is identical; `on` only ADDS gene-like loci. Promote iff: more real gene-like features lifted, validates clean, off-lines ⊆ on-lines.

| Dataset | gene-like types | pseudogene off/on | gene off/on | off⊆on (lift preserved) | on-only lines | validate-on (exit/errs) | PROMOTE |
|---|---|---|---|---|---|---|---|
| drosophila | gene,pseudogene | 0/21 | 4048/4048 | yes | +71 | 1/3 | yes |
| mouse_to_rat | gene,pseudogene | 0/124 | 945/945 | yes | +175 | 1/2 | yes |
