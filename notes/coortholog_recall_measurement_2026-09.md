# What the target's own annotation says LiftOn is missing

Recall has always been measured against the **source** annotation: of the
reference's genes, how many reached the target. That denominator cannot see the
thing that matters most on a distant transfer — a target genome that has *more*
genes than the source, because a whole-genome duplication gave it two copies of
what the source has one of. Source recall scores "one source gene, one target
model" as a perfect result whether or not the target actually has two genes
there.

Human → zebrafish is the clean case: the lift targets GRCz11, and GRCz11's own
RefSeq annotation is available as an independent answer key. Every number below
is against that annotation, not against human.

## The measurement

| | genes |
|---|---:|
| GRCz11 protein-coding genes | 33,564 |
| covered by some LiftOn gene | 13,783 |
| **not covered at all** | **19,781** |
| — of those, miniprot finds nothing either | 17,028 |
| — of those, **miniprot finds them** | **2,753** |
| ⟶ at a locus LiftOn already occupies | 326 |
| ⟶ **at a locus nothing occupies** | **2,427** |

The 17,028 are zebrafish genes with no human counterpart to lift from; no
homology method reaches them and they are not a defect. The interesting number
is the 2,427: miniprot found them, nothing was in the way, and LiftOn emitted
nothing.

A gene counts as "covered" if *any* LiftOn gene overlaps it by a single base.
That is deliberately generous, and generosity can only shrink the gap, which is
the safe direction for a number used to justify building something.

## Why the rescue does not reach them

Attributing each of the 2,427 to the gate that excludes it:

| | genes | |
|---|---:|---|
| reference gene already emitted elsewhere (dedup) | **2,051** | 84.5 % |
| outside the length-ratio sanity band | 317 | 13.1 % |
| unexplained | 59 | 2.4 % |

The separate-pass rescue dedups by **reference gene**: a reference gene already
emitted anywhere is never rescued again. Iteration 23 added that rule, and it is
exactly right at an *occupied* locus — it is what makes the pass strictly
additive and what kept `n_redundant` at 0 on all eight ladder cells.

In a duplicated genome it is the wrong key. The second zebrafish copy is not a
redundant model of the gene already emitted; it is a different gene, at a
different locus, that the target's own annotation lists separately.

## Why this is safe to change, and how it differs from the failures

Iteration 15 and Iteration 22 both failed by adding models that competed with
models already emitted — duplicates in one case, a suppression swap in the
other. Neither failure mode applies to a locus **nothing occupies**: no emitted
model is displaced, and `off ⊆ on` still holds by construction, which is the
property that let Iteration 23 pass a gate Iteration 22 could not.

So the change is narrow: allow a reference gene to be rescued a second time
**when the hit lands on a locus no emitted model reaches**, keeping the dedup
everywhere else.

## What still has to be proven

Headroom is upside-only, and the last three times it was trusted alone it was
wrong (Iterations 9, 15, 22 — each time the A/B found what the projection could
not). The decision belongs to the strict A/B on emitted output: `n_lost = 0`,
`n_redundant = 0`, no common-set identity regression, validity unchanged — and
this time also re-scored against the target's own annotation, since source
recall is precisely the instrument that cannot see the gain.

Reproduce: `c1/coortholog_truth.py` and `c1/why_free_loci_are_skipped.py`
(evaluation-only; neither runs or changes LiftOn).
