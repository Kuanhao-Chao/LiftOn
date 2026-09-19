# Where Step 7 actually spends its time

`LIFTON_PROFILE_STEP7`, whole genomes, cached `-L`/`-M` so the aligners do not
dominate, `-t 1` (cProfile serialises Python; the goal is *where*, not *how
fast*). Same build for both cells: `devel` at the v1.0.13 release.

| cell | Step-7 dispatch | features |
|---|---:|---|
| rice (same species) | 730 s | 130,893 translations |
| dog → cat (mammal) | 3,299 s | 539,740 translations |

CLAUDE.md warns that this profiler has contradicted a guessed ranking twice.
It did again: the improvement plan expected "parasail on long mammalian
transcripts". Parasail *is* first — but it is a C kernel reached through a thin
binding, and the only way to spend less time in it is to align less, which
changes output. The addressable cost is elsewhere, and it is substantial.

## Ranked, with the addressable share marked

| # | site | rice | dog → cat | addressable |
|---|---|---:|---:|---|
| 1 | `parasail.nw_trace_scan_sat` | 210 s (29 %) | 1,003 s (30 %) | no — C kernel |
| 2 | `sqlite3.Cursor.execute` | 151 s (21 %) | 298 s (9 %) | yes — query **count** |
| 3 | `Bio.Seq._translate_str` + `CodonTable.__getitem__` | **66 s (9 %)** | **269 s (8 %)** | **yes** |
| 4 | `windowed_align` (`_unique_anchors`, `_chain_colinear`, `_cigar_from_aln`) | — | **320 s (10 %)** | yes — mammalian only |
| 5 | `gff3_writer.encode_attribute_value` | 38 s (5 %) | 111 s (3 %) | yes |
| 6 | `coreutils.clone_attributes` | 16 s (2 %) | 148 s (4 %) | yes |
| 7 | `protein_maximization.get_protein_reference_length_single` | — | 58 s (2 %) | yes — 1.31 M calls |

## The finding worth acting on first

`Bio.Data.CodonTable.__getitem__` is called **62 million times on rice and 257
million times on dog → cat** — one Python-level `__getitem__` per codon per
translation, 15 s and 60 s of pure interpreter overhead. Together with
`_translate_str` around it, translation is **8–9 % of Step 7 on both regimes**,
which is more than the entire GFF3 writer.

That is the same code path the `transl_table` work just centralised into
`coding.translate`, so the optimisation has exactly one seam, and its
correctness is provable rather than argued: a flat 64-entry `{codon: amino
acid}` dict per table gives the same answer as Biopython for every
unambiguous codon, and anything else (an `N`, a gap, a lowercase base) falls
back to Biopython untouched. Exhaustive comparison over all 64 codons × every
NCBI table, plus fuzzing with ambiguity codes, is a complete equivalence proof
for a function with a 64-element domain.

Second is `encode_attribute_value` at 17.3 M calls on dog → cat, and third
(mammalian only) is the windowed aligner's anchor/chain phase at ~10 %.

## What this rules out

Nothing here supports a Step-7 *algorithmic* rewrite. The dispatch is already
fused and parallel; the remaining cost is per-feature constant factors in four
identifiable functions. Each is separately measurable and separately provable,
so S1 should be a sequence of narrow byte-identical changes, not one large one.
