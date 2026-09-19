# Second-locus rescue: three attempts at one gate

The measurement said a reference gene's second target copy is worth 2,051 genes
on human → zebrafish (`notes/coortholog_recall_measurement_2026-09.md`). Getting
it past the strict A/B took three designs, and each round's failure was
invisible to the round before it.

Every number below is human → zebrafish, whole genome, shared cached `-L`/`-M`,
`-t 8`, both arms on the same build. **Control:** the off arm run twice is
byte-identical, so `-t 8` is deterministic on this cell and every difference
below is signal.

## Round 1 — relax the dedup inside sub-pass A

| | |
|---|---:|
| genes added | 830 |
| **genes lost** | **119** |
| added models overlapping an existing one | 165 |
| transcripts regressed / improved | 33 / 7 |

The commit message claimed `off ⊆ on` held "by construction" because a second
copy was only placed at a free locus. That is false, and the A/B said so: each
acceptance commits its interval into the shared suppression tree, so an extra
gene placed early in the walk takes a locus a *later default* candidate wanted.
The Iteration-22 swap, in a new place.

## Round 2 — its own sub-pass, after sub-passes A and B

| | |
|---|---:|
| genes added | 694 |
| genes lost | 0 |
| genes moved | 6 |
| transcripts regressed | 0 |

The gene-level gate reads clean. The six moved genes are what matters: chasing
them showed `gene-CCND3` keeping its ID and start while **four of its isoforms,
identity 0.584–0.628, silently disappeared**, replaced by a second-locus model
scoring 0.356.

Sub-pass C ran before the isoform pass, so it occupied ground the default gene
needed to widen into, and `_extension_collides` then refused the isoform.
Counted at the transcript level: **21 lost, mean identity 0.662**.

**The gate was measuring the wrong unit.** A gene that survives while its
transcripts vanish scores as `n_lost = 0`. The scorer now counts transcripts.

## Round 3 — sub-pass C runs last, after the isoform pass

| | |
|---|---:|
| off arm vs canonical default | **byte-identical** (`e320cca9…`) |
| genes added, all tagged | **689**, mean identity 0.651 |
| genes lost / **transcripts lost** | 0 / **0** |
| default genes moved | **0** |
| transcripts regressed | **0** |
| **target genes newly covered** | **+690** (13,783 → 14,473 of 33,564) |
| added models overlapping an existing one | 8, all below the pipeline's own 0.10 gate |

Running last costs five genes (694 → 689): a widened default gene now holds
ground a second-locus model wanted. That is the right trade. It buys the
property the design claimed twice and only now has — the default output, every
gene, every isoform, every span, is provably untouched.

Second-locus genes get no isoforms, since the isoform pass has already run.
Worth revisiting; a single correct gene beats a wrong one with isoforms.

## Still owed before promotion

This is one cell. The divergence ladder decides whether the gain generalises or
is a zebrafish-shaped artefact of the teleost duplication — which is the honest
prior, since that duplication is exactly why this cell was chosen. The flag
stays off until then.

Reproduce: `c1_ab/run.sh` (arms), `c1_ab/replicate.sh` (the control),
`c1_ab/score.py` (the gate).
