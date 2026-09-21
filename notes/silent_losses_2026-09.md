# The things LiftOn discarded without counting

*2026-09-21. Cycle 4, parts 1 and 2.*

Four losses that were invisible, and one reported so badly it could not be
seen. The unifying fact: every instrument built to expose a silent loss was
itself silent, or drowned in its own noise.

## The rescue reported no losses because it counted none

`miniprot_rescue.py` is default-ON, abandons candidates in roughly forty
places, and contained **no `drop_ledger` call at all**. Five whole genomes
reported `dropped_features_total: 0`, which reads as "nothing was lost" and
meant "nothing was counted".

Classifying those forty sites:

* **~25 are deliberate filters** — identity floor, overlap suppression, length
  band, dedup by reference gene, the processed-pseudogene rule. Those are
  decisions, and they stay uncounted. Counting them would recreate the problem
  in the other direction.
* **~15 are lookups that failed** — an id that did not resolve, a reference
  protein that was never extracted, a length index with no entry, an exception
  swallowed into a `continue`.

Only the second group got classes: `miniprot_hit_unmapped`,
`miniprot_gene_unresolved`, `reference_protein_sequence` (reused),
`reference_feature_length_missing`, `rescue_candidate_error`.

Two distinctions that already existed and were being thrown away:

* `get_ref_ids_miniprot` returns `(None, None)` when a hit is not in the id map
  and `(None, ref_trans_id)` when the transcript is known but no gene could be
  found for it. All five callers collapsed both into one `continue`.
* `ref_features_len_dict.get(x)` returning `None` (not in the index) and
  returning `0` (in the index, genuinely no CDS rows) were both consumed by
  `if not ref_len: continue`. The first is a failure; the second is a correct
  decision about a non-coding gene.

Two swallows had no log at all. One shrank the recall universe, which silently
makes the adaptive identity floor **stricter genome-wide**; the other dropped a
hit out of the isoform index so its gene never saw that isoform. And a single
database error at the top of the pass made the whole rescue a no-op
indistinguishable, in every report, from "there was nothing to rescue" — on a
distant transfer that is most of the recovered genes.

### What the instrument then said

After wiring, the answer on this corpus is **still zero**:

| | rescued genes | drop total | output |
|---|---:|---:|---|
| human → zebrafish | 8,846 → 8,846 | 0 | byte-identical |
| drosophila | 18 → 18 | 0 | byte-identical |
| rice | 0 → 0 | 0 | byte-identical |

151,967 miniprot candidates passed through those guards on human → zebrafish
and not one was abandoned for a failed lookup. Every rejection really is a
decision.

That is worth stating plainly rather than dressing up: **nothing was fixed
here, because nothing was broken.** What changed is that the zero is now
evidence instead of silence, and a regression would be visible.
`tests/test_rescue_drop_accounting.py` exists to keep that true — it checks the
counters fire when there *is* something to count, which is the failure mode
this whole change is about.

## The childless-gene counter was 99.98 % false positives

`genes_emitted_without_children` exists because the Liftoff `-copies`
resolution bug emitted ~4,400 bare gene lines across the benchmark corpus and
nothing added them up. It counted **every** emitted gene with no children —
including the 10,626 single-row pseudogenes RefSeq itself declares in the CHM13
reference.

Compared against the reference that produced each file:

| | reported | real |
|---|---:|---:|
| human → CHM13 | 10,720 | **0** |
| dog → cat | 272 | **0** |
| human → zebrafish | 113 | **0** |
| rice | 25 | **2** |

11,130 reported, 2 real. A counter at that ratio cannot show a recurrence of
the bug it was built for, which is its only purpose. It now reports the
actionable case — the reference gave this gene children and we emitted none —
and keeps the raw tally as `bare_gene_lines_emitted`.

## A CDS spanning two exons was emitted twice

`Lifton_TRANS.add_cds` attached a CDS to **every** exon it overlapped, cloning
the feature, so the transcript emitted the same coding block more than once and
the protein counted those bases twice. An earlier fix had made each exon get
its own *copy* — correct for the aliasing bug of the day, but it left the
doubling in place.

It also blamed the input: *"The reference model is malformed here"*. On real
data that was false. The extra exon was one LiftOn had just created by
ingesting miniprot's redundant `stop_codon`. With that fixed the branch stops
firing entirely — CHM13 8,546 → 0, rice 4,468 → 0, bee 3,614 → 0, drosophila
2,632 → 0 — so what remains is a reference CDS that genuinely spans an intron.

It is now kept on the exon it overlaps most (ties broken toward the earlier
exon, deterministically) and the rest is counted as `cds_spanning_exons`.
Emitting it once loses the overhanging bases; emitting it twice is a wrong
protein, silently. The tests are constructed, because the corpus no longer
reaches this path.

## A reference index built for only one annotation shape

`get_ref_liffover_features` populated `ref_features_reverse_dict` and
`ref_trans_exon_num_dict` only in its 3-level branch. RefSeq's organellar
convention writes a plastid gene's exons twice — once directly under the gene
and once under its mRNA — so such a gene took the other branch and its mRNA was
never indexed. **12 genes in the rice reference, 0 in human RefSeq.**

A miniprot hit on one could not be mapped back to a gene, so every rescue
candidate for it was abandoned in silence.

### A correction I had to make to my own plan

I wrote that this was issue #37's bacterial `gene → CDS` shape and that "every
miniprot hit on such a genome is unrescuable". That was wrong twice over. A
bacterial annotation has no exon rows, so it takes the *other* branch and
indexes correctly — verified by running the real bacterial fixture, not by
reading the code. And the affected population is 12 genes, not a genome class.

This is the same error as the threading work earlier in the programme: the
mechanism was located correctly and its reach was counted in the wrong place.

### The two halves were one bug

Before this fix the new childless counter reported **1** real loss on rice
where the reference plainly showed 2. The missing one, `gene-OrsajCp001`, is an
organellar gene taking exactly the unindexed branch: `Lifton_feature.children`
stayed empty, so the counter read "the reference had no children either". With
the index fixed it reports 2.

## Verification

rice, drosophila and human → zebrafish are **byte-identical** to the previous
build, with rescued-gene counts unchanged. So the index fix closes a blind spot
without changing this corpus — no recall win is claimed. What changed is the
accounting: rice's real childless count 1 → 2, human → zebrafish 113 → 0,
against 113 and 25 bare gene lines respectively.
