## miniprot-only merge-candidate headroom (Iteration 15, Gap #1, read-only)

Default run on cached `-L`/`-M`, feature OFF, `LIFTON_MINIPROT_ONLY_DIAG` set (byte-frozen output). Per merge-fired transcript: the APPLIED best-of-{merge, Liftoff} identity vs what a **miniprot-only** CDS would score after ORF-rescue. **would-change (same-locus)** = miniprot-only strictly beats the applied result AND its CDS span overlaps the Liftoff transcript (coordinate-compatible with an in-place swap). **projected corpus mean-PI delta** = sum of same-locus gains / all transcripts (the number the A/B must beat the +0.003 bar on). Decision: would-change ~0 or projected < +0.002 on both -> NO-GO (document); else build the A/B.

| Dataset | corpus | merge-fired | same-locus | would-change | mean Δ over changed | projected corpus PI Δ | diff-locus would-change |
|---|---|---|---|---|---|---|---|
| drosophila | 7946 | 6699 | 6699 | 26 | 0.01479 | +0.00005 | 0 |
| mouse_to_rat | 3114 | 1983 | 1983 | 20 | 0.06325 | +0.00041 | 0 |
