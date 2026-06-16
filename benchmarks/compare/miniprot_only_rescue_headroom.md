## Weak-Liftoff miniprot-rescue headroom (Iteration 15, Gap #2, read-only)

Default run on cached `-L`/`-M`, `LIFTON_MINIPROT_RESCUE_DIAG` set (byte-frozen output). Per SUPPRESSED (overlapped) miniprot mRNA: the drop reason + what it WOULD score (`mp_aa`), coordinate-joined against `score.txt` for the **local Liftoff PI**. **recoverable** = a clean rescue at a weak locus: mp_aa >= 0.7 AND mp_aa >= local Liftoff PI + 0.05. Distinct from the Iter-13 threshold lever (this is Liftoff-quality-aware). GO if recoverable >= 10 & mean mp_aa >= 0.75 on any dataset; else NO-GO.

| Dataset | suppressed | pass-filters | mp_aa buckets | recoverable | mean mp_aa (rec) | mean local Liftoff PI (rec) | verdict |
|---|---|---|---|---|---|---|---|
| drosophila | 7222 | 5474 | {'<0.5': 33, '[0.5,0.7)': 111, '[0.7,0.85)': 542, '>=0.85': 4788} | 44 | 0.92496 | 0.68205 | GO (worth an A/B) |
| mouse_to_rat | 2165 | 1040 | {'<0.5': 26, '[0.5,0.7)': 40, '[0.7,0.85)': 67, '>=0.85': 907} | 8 | 0.987 | 0.91707 | NO-GO / negligible |
