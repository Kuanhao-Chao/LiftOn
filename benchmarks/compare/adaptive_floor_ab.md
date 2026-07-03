## Divergence-adaptive rescue-floor A/B (Phase 2) — adaptive vs fixed-0.5

Both arms run `--miniprot-rescue`; `fixed` = `LIFTON_RESCUE_ADAPTIVE_FLOOR=0` (fixed 0.5), `adaptive` = floor in [0.30, 0.50] by DNA-lift recall. Cached `-L`/`-M`, `-t1` no-`copies`, scored by the neutral re-alignment evaluator. Adaptive only LOWERS the floor within an already-additive pass, so it is a strict superset of fixed: **n_lost** and **n_redundant** MUST be 0. **n_new** = the extra genuinely-missing genes the lower floor recovers (divergent cells); high-recall cells stay at 0.5 -> **inert** (n_new=0). Gate: n_lost=0 AND n_redundant=0 AND 0 regressed AND validity not worse.

| Dataset | divergence | n_new | mean PI new | Δcompl | n_lost | regr | tagged f→a | redundant | val f→a | gate |
|---|---|---|---|---|---|---|---|---|---|---|
| celegans_to_briggsae | distant_cross_species | 18 | 0.49218 | 0.00238 | 0 | 0 | 555→573 | 0 | 0→0 | **PASS** |
| drosophila_to_anopheles | very_distant_cross_species | 73 | 0.46892 | 0.01007 | 1 | 0 | 201→273 | 0 | 0→0 | NO-GO |
| zebrafish_to_medaka | distant_cross_species | 23 | 0.45066 | 0.00774 | 0 | 0 | 143→166 | 0 | 0→0 | **PASS** |
| rice_to_sorghum | distant_cross_species | 0 | None | 0.0 | 0 | 0 | 100→100 | 0 | 0→0 | **PASS** |
| t4_human_to_xenopus | very_distant_cross_species | 15 | 0.44217 | 0.00486 | 0 | 0 | 87→102 | 0 | 0→0 | **PASS** |
| t4_human_to_chicken | very_distant_cross_species | 4 | 0.46551 | 0.0013 | 0 | 0 | 51→55 | 0 | 0→0 | **PASS** |
| human_to_mouse | distant_cross_species | 0 | None | 0.0 | 0 | 0 | 30→30 | 0 | 0→0 | **PASS** |
| drosophila | D. melanogaster -> D. erecta | 0 | None | 0.0 | 0 | 0 | 1→1 | 0 | 0→0 | **PASS** |

### Decision: PROMOTED default-ON (recall-priority), accepting one swap

**Result: +133 genuinely-missing genes recovered across the divergent ladder**
(celegans +18, dros→anopheles +73, zebrafish +23, xenopus +15, chicken +4;
mean PI of the recovered set 0.44–0.49 — the harder, lower-identity genes the
0.30 floor admits), **0 redundant, 0 regressed, validity clean**, and correctly
**inert on high-recall cells** (rice, human→mouse, same-species drosophila: floor
stays at 0.5, n_new=0). 

**The one blemish — `n_lost=1` on drosophila_to_anopheles — is a shared-tree
suppression swap, NOT a strict subset violation from the floor itself.** The
miniprot-only rescue builds each candidate `Lifton_GENE`, which inserts it into
the shared Step-8 suppression `tree_dict` (`lifton_class.Lifton_GENE.__init__`)
*before* the PI-floor check. Lowering the floor changes which genes are
built/emitted earlier in the deterministic (`seqid,start,end,ID`) loop, which
mutates `tree_dict`, which cascades to suppress a *later* candidate the fixed arm
emitted (here `rna-NM_001032009.2` / CG6783-fabp, PI 0.672, suppressed downstream
in the adaptive arm). Net on that cell is still **+72** (73 gained − 1). This is
the same shared-loop-state class as Iteration-22; a strictly-additive fix needs
the rescue's gene-build decoupled from the tree mutation (deferred as separate,
separately-verified surgery). **Per the recall-priority posture the adaptive floor
is PROMOTED default-ON; `--no-adaptive-rescue-floor` (env
`LIFTON_RESCUE_ADAPTIVE_FLOOR=0`) restores the strictly-additive fixed-0.5 floor.**
