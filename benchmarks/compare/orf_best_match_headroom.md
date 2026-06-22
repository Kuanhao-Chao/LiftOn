## ORF best-match headroom diagnostic (Iteration 9, read-only)

Default run on cached `-L`/`-M`, flag OFF, `LIFTON_ORF_DIAG` set (byte-frozen output). Per rescued transcript: applied longest-per-frame identity vs the top-K(=3)/unbounded best-match identity. **gate-crossed** = the broadened pick differs from the longest AND clears the 1 % rescue gate (would change the emitted CDS). **projected corpus mean-PI delta** = sum of gains over gate-crossed transcripts / all transcripts (the number the A/B must beat the +0.003 bar on). Decision: gate-crossed ~0 or projected < ~+0.002 on both -> NO-GO (keep opt-in, document); otherwise build the A/B.

| Dataset | corpus | rescued | multi-ORF | gate-crossed | mean Δ over crossed | projected corpus PI Δ | K=3 insufficient |
|---|---|---|---|---|---|---|---|
| drosophila | 7946 | 6979 | 6958 | 5 | 0.16888 | +0.00011 | 28 |
| mouse_to_rat | 3114 | 2305 | 2304 | 12 | 0.07894 | +0.00030 | 6 |
