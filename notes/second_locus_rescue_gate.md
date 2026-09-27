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

## The ladder: does it generalise?

Human → zebrafish is the most favourable cell that exists for this idea — the
teleost whole-genome duplication is the reason the idea exists — so one cell
proves nothing about anything else. Eight cells, deterministic `-t 1`, cached
`-L`/`-M`, both arms pinned to one build.

**What the ladder can and cannot say.** The neutral evaluator keys on
`ref_mrna_id` and keeps the best model per reference transcript, so a second
target gene for a reference gene that is *already scored* is invisible to it.
Reference-keyed recall is exactly the instrument that cannot see this feature.
The ladder therefore runs a **safety** gate — nothing lost, nothing regressed,
nothing placed on top of an emitted model, validity not worse — and reports how
many genes were placed and how good they look. Whether they are real needs the
target's own annotation, which is on disk for zebrafish only.

| cell | placed | per 1,000 ref proteins | mean PI |
|---|---:|---:|---:|
| rice → sorghum | 53 | **9.06** | 0.773 |
| human → zebrafish *(whole genome)* | 689 | 4.74 | 0.651 |
| zebrafish → medaka | 6 | 2.02 | 0.685 |
| C. elegans → briggsae | 13 | 1.72 | 0.673 |
| human → chicken | 5 | 1.62 | 0.780 |
| human → xenopus | 4 | 1.30 | 0.707 |
| drosophila → anopheles | 5 | 0.69 | 0.548 |
| human → mouse | 2 | 0.65 | 0.715 |
| **drosophila (same species)** | **0** | **0.00** | — |

**Safety gate: 8/8.** Zero lost, zero regressed, zero placed above the
pipeline's own `-overlap 0.10` gate, validity unchanged, on every cell.

Three things make this more than a count:

* **The same-species control places exactly nothing.** A transfer between two
  *D. melanogaster* assemblies has no second copies to find, and the feature
  finds none. Had it placed anything there, the rest of the table would mean
  nothing.
* **It is not a zebrafish artefact.** The strongest per-gene rate is rice →
  sorghum, not the teleost cell.
* **The ordering tracks duplication history.** Grasses (ancient WGD plus heavy
  segmental duplication) above teleosts, above nematode/bird/amphibian, above
  fly and mouse, above zero for same-species. That is the ranking the mechanism
  predicts, arrived at without being told.

Every mean identity sits above the 0.5 rescue floor.

## Where this leaves promotion

The evidence supports the feature doing what it claims, safely, across
divergence regimes. It does not yet establish that the placed genes are correct
anywhere except zebrafish, where the target's own annotation says +690 of them
are real. Whether that is enough to turn it on by default is a judgement about
how much weight one target annotation carries, and that belongs to the
maintainer, not to this note.

Reproduce: `python -m benchmarks.compare.second_locus_ab`
(`LIFTON_AB_PYTHONPATH` pins the build, `LIFTON_AB_WORK` points at the subset
trees).
