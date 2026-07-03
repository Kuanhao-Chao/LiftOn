## Phase-4 Strategy-B feasibility probe — incremental headroom OVER candidate-3

Byte-neutral (reuses the candidate-3 ladder A/B `on` + standalone `miniprot` TSVs, neutral re-scorer). **strategy_b_incremental_pi** = mean(pi_miniprot − pi_lifton_on) on the transcripts candidate-3 FIRED (tagged LiftOn_miniprot): Strategy B's hybrid uses the SAME miniprot CDS, so this is the protein it would add beyond candidate-3's emitted model. **tagged |gap|** = how close candidate-3's tagged emissions already are to miniprot (≈0 => candidate-3 emits miniprot's model). **oracle headroom** = mean(max(0, pi_miniprot − pi_lifton_on)) over the common set — the upper bound of a perfect best-of selection (a firing-criterion lever, NOT Strategy B). Gate: strategy_b_incremental_pi >= +0.003.

| Dataset | tagged∩common | Strat-B incr PI | tagged \|gap\| vs mp | oracle headroom | verdict |
|---|---|---|---|---|---|
| celegans_to_briggsae | 30 | -0.26671 | 0.3347 | 0.04258 | **NO-GO** |
| drosophila_to_anopheles | 2 | -0.07356 | 0.07356 | 0.024 | **NO-GO** |
| zebrafish_to_medaka | 1 | -0.14054 | 0.14054 | 0.07831 | **NO-GO** |
| rice_to_sorghum | 41 | -0.35539 | 0.40824 | 0.05802 | **NO-GO** |
| t4_human_to_chicken | 20 | -0.51259 | 0.5578 | 0.04868 | **NO-GO** |

### Verdict: NO-GO — candidate-3 subsumes (and exceeds) Strategy B on the protein axis

**Strategy B is not worth building as an accuracy lever.** Two independent
arguments, both giving incremental protein Δ < +0.003:

1. **By construction (Δ = 0):** Strategy B's hybrid and candidate-3 both derive
   the coding structure from *miniprot's own CDS junctions*; the only difference
   (keep DNA-lift UTR) is untranslated, so they emit the **same protein**. And a
   LiftOn-native Strategy B would apply the same ORF-rescue candidate-3 does →
   identical CDS → identical protein identity. Zero incremental.

2. **Empirically (Δ < 0):** on the transcripts candidate-3 fired, its emitted
   model already scores **0.07–0.56 ABOVE raw miniprot** (`Strat-B incr PI`
   strongly negative) — because candidate-3 = miniprot CDS **+ LiftOn ORF-rescue
   + containment-normalize**, which beats miniprot's raw protein-to-genome model.
   Strategy B's building block (miniprot CDS) is therefore already fully captured,
   and then improved, by candidate-3. A hybrid rebuild cannot recover protein
   candidate-3 has already surpassed.

Candidate-3's decision bar (beat the 2-way merge `max(Liftoff+ORF, chain+ORF)`) is
also **stronger** than Strategy B's (beat a pure DNA chain), so wherever Strategy
B would fire, candidate-3 has already considered and matched-or-beaten it.

**Residual `oracle headroom` (0.024–0.078) is a candidate-3 FIRING-criterion
lever, NOT Strategy B.** It is the gap between LiftOn's emitted model and a
perfect best-of-{lifton, miniprot} *selection* (which no method can reach at lift
time — the neutral re-alignment score is unknown while lifting). Closing it needs
a better firing criterion (which miniprot-only transcripts to emit), not miniprot
*structure* (candidate-3 already gives that). Prior iterations already probed this
firing lever and found it net-negative: Iter-9 (ORF best-match, net −), Iter-15
(miniprot-only rescue, duplicate-dominated), Iter-22 (in-loop rescue, off⊄on
swaps). So the residual headroom is known-hard and out of Strategy B's scope.

**Decision: NO-GO. No `lifton/` source written** (matches the Iter-4/9/13/15
feasibility NO-GO precedent); the probe + this note are the audit trail. Strategy
B's only genuine residual value — UTR/transcript-structure completeness on
candidate-3's CDS-only emissions — is a separate STRUCTURAL-completeness lever
(now measurable via the Phase-3 exon Sn·Sp metrics), deferred.
