## Regime-gated miniprot-only rescue A/B — Iteration 22

`rescue_<floor>` = `--miniprot-rescue` with PI floor `<floor>`; cached `-L`/`-M`, scored by the independent re-alignment evaluator (keyed by ref_mrna_id). **n_added** = coding ref transcripts rescue emits that the DNA lift missed (genuine completeness); **n_redundant** = rescued models that did NOT add a new ref id (the Iter-15 duplicate failure mode — MUST be 0); **n_lost** = off ⊆ on gate (MUST be 0). Gate: n_added>0 AND n_lost=0 AND 0 regressed AND n_redundant=0 AND mean PI added ≥ floor AND validity not worse.

| Dataset | floor | n_added | mean PI added | Δcompl | n_lost | regr | tagged | redundant | val d→r | gate |
|---|---|---|---|---|---|---|---|---|---|---|
| drosophila | 0.5 | 1 | 0.55952 | 0.00014 | 0 | 0 | 1 | 0 | 3→3 | **PASS** |
| human_to_mouse | 0.5 | 38 | 0.69651 | 0.01069 | 5 | 0 | 43 | 5 | 3→3 | NO-GO |
| celegans_to_briggsae | 0.5 | 614 | 0.72157 | 0.06715 | 107 | 1 | 632 | 18 | 3→3 | NO-GO |
| rice_to_sorghum | 0.5 | 109 | 0.66522 | 0.01675 | 11 | 5 | 175 | 66 | 3→3 | NO-GO |
| drosophila_to_anopheles | 0.5 | 213 | 0.6398 | 0.02482 | 33 | 0 | 215 | 2 | 3→3 | NO-GO |
| zebrafish_to_medaka | 0.5 | 152 | 0.69868 | 0.04643 | 14 | 2 | 162 | 10 | 2→2 | NO-GO |
| t4_human_to_chicken | 0.5 | 57 | 0.72932 | 0.0162 | 7 | 0 | 61 | 4 | 3→3 | NO-GO |
| t4_human_to_xenopus | 0.5 | 98 | 0.71775 | 0.02819 | 11 | 1 | 100 | 2 | 2→2 | NO-GO |

## Interpretation — Iteration 22 decision (NO-GO for default; kept as EXPERIMENTAL opt-in)

These numbers are **deterministic** (`-t1`, no `-copies`, cached `-L`/`-M`): a re-run reproduced them byte-for-byte, and an earlier `-t8 -copies` run gave the identical table — so the effect is the flag, not `-copies` run-to-run noise.

**The recall lever is real.** Every distant/very-distant cell gains genuine new coding transcripts the DNA lift missed — celegans **+614** (net +507 after losses), drosophila→anopheles +213, zebrafish→medaka +152, human→xenopus +98, human→chicken +57, rice→sorghum +109, human→mouse +38 — at mean PI 0.64–0.73 (all above the 0.5 floor), Δcompleteness up to +0.067. This matches the headroom (14,288 genuine-new@0.5 on the very-distant tier) and lifton2's shipped +0.047.

**But the naive in-loop implementation is not default-ready — it fails the strict gate two ways:**
1. **`n_lost` > 0 (off ⊄ on).** Root cause: `lifton_class.Lifton_GENE.__init__` adds every constructed gene to the shared Step-8 suppression `tree_dict` (`lifton_class.py:100-102`). A rescued gene, emitted mid-loop, then suppresses a *later* overlapping miniprot candidate the default would have emitted — a **swap** of which competing missing-gene representative wins at a contested target locus (net recall still rises, but ON is not a strict superset of OFF).
2. **`n_redundant` > 0** (e.g. rice→sorghum 66). A ref transcript with multiple miniprot hits gets a second, band-failing model emitted by the rescue; the tight band had implicitly deduped these.

**Why lifton2 gets 0-redundant and this does not:** lifton2's rescue is a **separate post-lift pass** with its own suppression tree + a `copy_id`/`lifted_gene_bases` dedup guard, run after the whole DNA lift; this parent cut rides inside the per-mRNA Step-8 loop on the *shared* gating tree with no dedup.

**Decision:** distinct from the net-*negative* Iter-15 (duplicate-dominated, reverted). Iter-22 is **net-positive on recall** but not clean → **NO-GO for default promotion**, **kept as an experimental opt-in `--miniprot-rescue`** (default OFF is byte-identical, proven by the 24-cell matrix + a drosophila full-GFF3 `cmp`). The **default-ready form is a separate-pass rescue** (process all default Step-8 emissions first, then fill only still-empty loci, with a ref-id dedup guard — the lifton2 Iter-11 shape) — the documented next step.
