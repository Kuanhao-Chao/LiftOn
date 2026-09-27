# Ranking rescue candidates by quality: measured marginal

The rescue walks miniprot's transcripts in file order and takes the first one
that clears every gate. miniprot ranks its own hits (`Rank=1` is its best), so
walking by position can accept a `Rank > 1` hit while a better one for the same
reference transcript was available and would also have passed. The improvement
plan proposed ordering candidates by score before the placement loop.

Measured on the shipped human → zebrafish whole-genome output — the cell with
by far the most rescue activity:

| | |
|---|---:|
| rescued transcripts in the output | 59,157 |
| matched back to their miniprot hit | 56,531 |
| **built from a `Rank > 1` hit** | **4,754 (8.4 %)** |
| **median identity cost** | **0.0000** |
| mean identity cost (on that subset) | 0.0031 |
| worst single case | 0.2558 |

So the behaviour is real and more common than the plan estimated — 8.4 %, not
the ~15 % of a smaller sample, but across far more transcripts. What is not
real is the cost. **The median secondary hit is exactly as good as the best one
available.** miniprot's rank orders by its own alignment score, and where two
hits for the same protein both clear LiftOn's floor and its length band, they
are usually the same model found twice.

Weighting the 0.0031 mean by the 8.4 % that are affected puts the corpus-wide
mean-identity effect at **≈ 0.0003** — two orders of magnitude below the
+0.003 promotion bar this repo uses, and in the same range as Iteration 9's
ORF-best-match projection (+0.00011 / +0.00030), which the A/B then showed to
be net-negative once downstream effects were included.

## Decision

**Deprioritised, not refuted.** Identity is not where the value is. There is a
second, untested channel: a better-ranked hit might clear a gate the
worse-ranked one fails, which would change *which* genes get placed rather than
how well. That is a recall effect, and this measurement — taken from emitted
output — cannot see it.

It is worth revisiting only after the free-locus work
(`notes/coortholog_recall_measurement_2026-09.md`), which addresses the same
axis and is worth 2,051 genes on this cell rather than 0.0003 mean identity.

Reproduce: `a3/baseline.py` (evaluation-only; does not run or change LiftOn).
