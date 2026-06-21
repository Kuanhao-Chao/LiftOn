# LiftOn deep audit — software, benchmark, evaluation directions + design/analysis/conclusion verification

**2026-06-20. Read-only deep scan.** Method: an 8-agent parallel audit (engine / benchmark / evaluation scans + two adversarial design-and-claims auditors) → adversarial verification of each finding (filter dead wells + ungrounded) → synthesis, cross-checked against independent data re-derivation by the main loop. The multi-agent run was partially rate/session-limited; the engine-dimension agents and the auto-synthesis were finished by the main loop from `improvement_headroom.{md,py}`, the CLAUDE.md iteration history, the lifton2 shipped wins, and direct re-derivation from the per-transcript eval TSVs. Every load-bearing claim below was re-derived from `fourway_results.json` / the eval TSVs, not taken on an agent's word.

---

## Executive summary (the four things that matter)

1. **The headline accuracy claim is overstated on the very-distant tier — a real flaw, confirmed two independent ways.** "LiftOn leads the best single baseline on all 17, the lead growing to +0.129" is computed as *mean-PI-over-each-tool's-own-recovered-set*. On the very-distant tier devel recovers only the easy DNA-anchored 4–24%, so its set-mean is compared against miniprot's mean over a far larger, harder set. **Apples-to-apples on the transcripts both recover, miniprot is as accurate or *more* accurate than devel on 4 of 5 very-distant cells** (chicken −0.011, xenopus −0.014, At→rice −0.009, zebrafish **−0.095**). devel's per-transcript accuracy advantage is **genuine and strongest on the *distant* tier** (apples-to-apples +0.023 to +0.032), not the very-distant tier. **Fix the conclusion**, and add the metric that would have caught it (Part 2 §A).

2. **The single highest-value *evaluation* upgrade is a joint recall-vs-identity metric** (apples-to-apples common-set PI + a recall@PI curve / F-like score). It is the direct remedy for finding #1 and converts the benchmark's biggest weakness (two independent axes that can be traded against each other and reported separately) into an honest combined view.

3. **The single highest-value *software* direction is regime-gated miniprot-only rescue** (lifton2 already shipped it: +0.047 mean completeness, 0 redundant on distant pairs). It is the only lever that addresses the real very-distant deficit, which finding #1 shows is **recall, not per-transcript accuracy**. The close-pair accuracy/completeness levers remain dead wells (Iter-4/9/13/15) — do **not** re-test them on close pairs.

4. **The benchmark's biggest credibility gaps are coverage-shape, not size**: the very-distant tier is dominated by human-as-reference (3/5 fulls), there is one annotation source (33/34 RefSeq), and there is no same-species mammalian *full* genome and no curated-target ground truth. These don't change today's conclusions but they bound how far the conclusions generalize.

---

## Part 1 — Verification of the design, analysis, and conclusions

### 1A. FLAW (confirmed) — very-distant accuracy "lead" is a denominator artifact
**Status: material flaw in the *conclusion* (not the data). Confirmed by (i) the design-audit agent reading `evaluator.py:386-426` + `fourway_compare.py:262`, and (ii) main-loop apples-to-apples re-derivation from the per-transcript TSVs.**

`devel_vs_best_baseline.meanpi = mean_PI(devel set) − mean_PI(best-baseline set)`, each over its own recovered set. When completeness differs wildly the two means are over different populations. Re-derivation on the **common** recovered transcripts (the fair per-transcript comparison):

| Tier | reported lead | apples-to-apples (devel − stronger baseline, common set) | verdict |
|---|---|---|---|
| same-species (5) | +0.0007 … +0.0124 | +0.0007 … +0.0119 (positive) | genuine ✓ |
| cross (drosophila) | +0.0092 | +0.0044 | genuine ✓ |
| close (3) | +0.0169 … +0.0264 | +0.013, +0.017, **−0.003 (tomato→potato)** | mostly genuine |
| **distant (3)** | +0.0305 … +0.0373 | **+0.023 … +0.032 (positive — the strongest real win)** | genuine ✓ |
| **very-distant (5)** | +0.016 … **+0.129** | **−0.095, −0.014, −0.011, −0.008, +0.001** | **artifact ✗** |

**The report claim** (§3.1, "where LiftOn places a transcript at extreme distance, it places it markedly more accurately than either single-method baseline") **is false vs miniprot on the very-distant tier.** Honest reframing: *LiftOn's per-transcript accuracy advantage is real and peaks on the distant cross-species tier; at very-distant divergence it trades recall heavily and its higher aggregate mean reflects only the easy DNA-anchored subset it recovers — on the shared transcripts miniprot is at least as accurate.*

**Recommended report edits:** soften Table 3 / Fig 1B framing at the very-distant tier; replace the §3.1 "more accurately than either baseline" sentence with the apples-to-apples truth; report the common-set metric (Part 2 §A) alongside the set-means.

### 1B. Minor concern — "grows monotonically with divergence"
Tier-**means** are monotonic (same +0.0064 → cross +0.0093 → close +0.0205 → distant +0.0331 → very-distant +0.0744), so the tier-level statement holds. But **per-cell it is not strictly monotonic** — `human_to_zebrafish` (+0.0164) sits below every distant-tier cell (+0.030–0.037). Soften "monotonically" to "the tier-average lead grows with divergence," or annotate zebrafish as the exception. (And note: under the apples-to-apples metric of 1A the very-distant "growth" reverses, which is the deeper reason to reframe.)

### 1C. Secondary concern — very-distant tier is reference-confounded
3 of the 5 very-distant fulls use **full human** as the reference (zebrafish, chicken, xenopus; n_coding 144329 each — *not* chr20; the zebrafish cell's `species` string "Human GRCh38 chr20" is a stale label to fix). So "lead grows with divergence" partly measures "human genes onto phylogenetically distant targets." Mitigated (not removed) by At→rice and dros→bee being non-human very-distant cells that show the same denominator pattern. Recommend diversifying very-distant references (Part 3).

### 1D. What is SOUND (held up under adversarial check — use with confidence)
- **`best(LO,MP)` baseline = the stricter comparison, not flattering.** liftoff is the better baseline on 3 cells, miniprot on 14; `max()` forces devel to beat the *higher* of the two. Confirmed sound.
- **Cross-tool transcript matching is fair** (copy-suffix-aware id for Liftoff/LiftOn, `Target` for miniprot; transcript-space proteins handled). Sound.
- **v1.0.8 partial-crash scoring is fair and correctly caveated** — the partial mean PI (arabidopsis 0.99884) is *above* devel's complete-run mean, i.e. survivorship works *against* devel, and the report flags the partials as not-comparable. Sound.
- **The crash-fix attributions are correct**: all 5 v1.0.8 crashes are the same `__str__` defect (verified in the stable stderr logs); the dog→cat devel crash + recovery is real (re-lift 11,048 → 61,276, byte-identical on valid annotations). Sound.
- **The completeness/PI numbers are internally consistent** — all 17 cells re-derive 0-mismatch from the per-transcript TSVs.

### 1E. Lower-severity hygiene
- `t1_soybean` is a registered-but-empty cell (`tools:[]`) — correctly excluded from the merge/figures, but it is a silent coverage hole; either score it (its inputs are fetched) or drop it from the registry.
- subset chromosome pinning (AUTO_LARGEST_CODING ref / dominant-syntenic target / WHOLE) is a reasonable, documented heuristic (minor selection-bias caveat only; the subsets are the controlled-engine arm, not the headline).

---

## Part 2 — Evaluation comprehensiveness roadmap (ranked)

The evaluator currently measures protein identity (mean/median/%identical over recovered coding), coding/feature completeness, dna identity, and a capped GFF3-error count. It is **alignment-tolerant and recovery-denominated** — which is exactly what hid finding #1. Ranked additions:

1. **[HIGH — fixes 1A] Joint recall-vs-identity view.** (a) **common-set (apples-to-apples) mean PI**: for each tool pair, mean PI over transcripts both recover — the honest head-to-head. (b) a **recall@PI-threshold** curve / a single **F-like score** combining completeness and identity so a tool cannot win accuracy by sacrificing recall. Cheap: computable now from the existing per-transcript TSVs (already prototyped in this audit). This is the most important evaluation change.
2. **[HIGH] Paired statistical significance for tool deltas.** Same-species leads are +0.0007–0.0018; with no per-transcript paired test (sign-test / Wilcoxon / bootstrap CI) "leads on all 17" can't distinguish a real edge from noise. Add a paired test + an `n_improved : n_regressed` count per cell (the engine A/Bs already use this convention). Cheap, byte-neutral, high defensibility.
3. **[HIGH] Structural accuracy: intron-chain exactness + exon-level Sn/Sp.** Protein/DNA identity is alignment-tolerant — a shifted/merged exon boundary that still translates near-identically scores ~1.0. Splice-junction (intron-chain) exactness is the gold-standard structural metric for a lift-over and is exactly where a homology lift can silently differ. Not computed anywhere today. Most discriminating on close/distant. (Also the natural metric to justify the deferred Iter-21 `update_cds_list` sort-repair.)
4. **[HIGH] Start/stop-codon correctness + ORF-validity rate of the *emitted* model.** LiftOn's ORF rescue is a core algorithm but the evaluator never checks whether the output is a valid ORF (correct start, single in-frame stop, no internal stop) — only identity to reference. A direct correctness axis.
5. **[MED] Uncap the validity metric.** The GFF3-error count is hard-capped at 50/check (`gff3_validator.py:168`), so several cells report saturated round numbers — expose per-check counts / uncap for the benchmark so the validity comparison is real.
6. **[MED] Conserved/single-copy (BUSCO-style) completeness slice.** Overall completeness at very-distant (e.g. 0.148) is dominated by genes with no true ortholog in the target; a conserved-gene slice separates "missed a real ortholog" from "no ortholog exists" and makes the recall story interpretable.
7. **[MED] Multi-copy gene handling.** Evaluation dedups to the best copy per ref id; surface copy-level precision/recall (quantifies the duplicate-model failure mode behind the Iter-15 NO-GO).
8. **[REJECTED] UTR coordinate-concordance** — verified near-zero headroom on this corpus (the references largely lack explicit UTR features to score against); revisit only with a UTR-rich annotation source.

---

## Part 3 — Benchmark completeness roadmap (ranked)

Current: 17 full + 34 subset across a 4-tier divergence ladder; strong on divergence depth, weak on coverage *shape*.

1. **[HIGH] Diversify the very-distant references off human.** 3/5 very-distant fulls are human-referenced (finding 1C). Add non-human very-distant fulls (e.g. plant→plant monocot↔eudicot beyond At→rice, a non-human vertebrate→vertebrate deep split) so the very-distant conclusions aren't a human-as-reference statement.
2. **[HIGH] Add a second annotation source/dialect.** 33/34 cells are RefSeq. Ecosystem dialect handling (Ensembl/GENCODE, FlyBase, WormBase, TAIR/Phytozome — `transcript` vs `mRNA`, `gene_biotype`, Parent/ID conventions) is divergence-independent and historically where product bugs hide (Iter-20 was exactly an mRNA-vs-transcript dialect bug). Empirically ~1 latent dialect/scale bug surfaced per expansion — high bug-finding yield.
3. **[MED] Add a same-species *full* mammalian cell** (strain/haplotype/pangenome, e.g. an HPRC human haplotype or mouse strain). The same-species full tier is all plant/insect today; a mammalian same-species full anchors the easy end for the clade the very-distant tier leans on.
4. **[MED] One curated-target ground-truth cell** (a target with its *own* independent manual annotation), enabling an *absolute* accuracy read instead of self-recovery-vs-reference. Most impactful on same/close; modest for the relative claims (verified marginal), but uniquely able to validate the evaluator itself.
5. **[LOW/MED] A reference-quality contrast** (draft-ref vs T2T-ref of the same pair) to substantiate the report's "T2T/pangenome era" framing, which the current all-RefSeq-grade corpus never tests.
6. **[REJECTED] gene-density/genome-size stratification** beyond clade — verified false premise (size co-varies with clade here; already inspectable for free).

---

## Part 4 — Software / engine improvement roadmap (ranked)

The engine is the **most-explored** dimension: close-pair accuracy and completeness levers are mapped dry (Iter-4 chaining sync, Iter-9 ORF best-match, Iter-13 Step-8 relax, Iter-15 miniprot-only fallback — all NO-GO *on close pairs*). The new very-distant tier + finding #1 reopen exactly the **recall** lever, and the sibling lifton2 has already shipped the relevant wins — so the highest-value engine work here is **regime-aware backports under the byte-identity contract**, not novel invention.

1. **[HIGH] Regime-gated miniprot-only rescue (completeness).** Finding #1 proves the very-distant deficit is *recall*, and `improvement_headroom` quantifies the clean additive headroom (genuine-new @PI≥0.5: ~14,288 transcripts on the very-distant tier with 0 overlap to existing devel features). **lifton2 already shipped this** (Iter-11, default-on floor-0.5: +0.047 mean completeness, **0 redundant on all distant datasets**) — the parent's Iter-15 NO-GO was close-pair-specific and does not apply here. Backport recipe: flag-OFF → byte-neutral feasibility diagnostic (already have the headroom) → A/B on the distant+very-distant tier with the `off⊆on` + `0-redundant` gate → promote. Highest expected impact.
2. **[MED-HIGH] 3-way / best-of-structure merge (per-transcript accuracy at distance).** Add a standalone miniprot model as a merge candidate so a broken Liftoff stub loses to miniprot's full model. The apples-to-apples data shows real per-transcript headroom precisely where the current merge underperforms (very-distant common set, miniprot ahead by up to 0.095). lifton2's **Strategy-B splice-junction rebuild** (Iter-9, +0.0055–0.0147 final mean-PI end-to-end) is the deeper, shipped version. Regime-gate to distant/very-distant; close-pair gain is dry (Iter-15 Gap#1 +0.00005).
3. **[MED] map_failed recovery (regime-agnostic, clean).** ~the cleanest small win in `improvement_headroom` idea C (e.g. maize 513 map_failed, ~80% recoverable): a region lifts but the protein didn't map → recover via the existing ORF path. Orthogonal to 1–2, smallest blast radius.
4. **[MED] Proper `update_cds_list` non-monotonic root-cause fix (robustness).** Iter-21 shipped only the *write-guard* (drop-and-log the inverted-coordinate transcript); the deeper cause — Case-2 rebuilding `self.exons` from a non-monotonic chained `cds_list` without re-sorting — is documented but unfixed. A re-sort would *keep* the transcript correctly instead of dropping it, but it is **not byte-neutral** (changes exon order on some passing transcripts) → its own flag + golden re-verification + the new intron-exactness metric (Part 2 §3) to prove it. Low frequency, real correctness.
5. **[LOW] miRNA-under-primary_transcript** nested-product handling (report 1/428) — a known structural gap; small, isolated.
- **Dead wells — do NOT re-open on close pairs:** Iter-4 chaining-sync, Iter-9 ORF best-match (net-negative), Iter-13 Step-8 threshold relax, Iter-15 weak-Liftoff/miniprot-only fallback. Only the *very-distant/distant* regime legitimately reopens #1–#2.

---

## The coherent thesis

Finding #1 (very-distant "accuracy lead" is a recall sacrifice, not a per-transcript win) ties the three roadmaps together: the **evaluation** fix (joint recall-vs-identity, Part 2 §A) is what *measures* it honestly; the **software** fix (regime-gated miniprot-only rescue, Part 4 §1, lifton2-proven) is what *closes* the real deficit (recall at distance); and the **benchmark** fix (non-human very-distant references + a second dialect, Part 3 §1–2) is what makes the conclusion *generalize*. None of these is a dead well, and the largest single action is correcting the report's very-distant accuracy framing (Part 1 §1A) — which the data already supports.
