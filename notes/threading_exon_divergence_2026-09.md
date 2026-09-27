# `--threads N` did not equal `--threads 1`

*2026-09-21. Cycle 3, item P1.*

## What was wrong

`run_liftoff.lifton_add_trans_exon_cds` asked the feature database for

```python
l_feature_db.children(locus, featuretype='exon', order_by='start')
```

No `level`, so gffutils answers recursively: every exon **descendant**, at any
depth. The Step-7 proxy (`locus_pipeline._LFeatureDbProxy`) has always served
that signature from a cache both walkers populate with the **level-1** result.

So a serial run read the real database and a threaded run read the proxy, and
on any locus where those two answers differ the emitted GFF3 differed. That is
the one contract this repository is organised around.

## Which input actually carries the shape — and the first answer was wrong

A feature must have a level-1 exon **and** a deeper exon. The obvious candidate
is RefSeq's microRNA precursor (`primary_transcript` → its own exon, plus a
nested `miRNA` → that miRNA's exon), and the reference annotations are full of
them: 1,915 in the human RefSeq reference, 46 in `test/GRCh38_chr22.gff3`.

**That is the wrong database.** `lifton_add_trans_exon_cds` reads
`l_feature_db` — Liftoff's *output*, not the reference. Liftoff does not lift
the nested miRNA transcript, so the shape never reaches the query:

| Liftoff output | features whose two exon queries disagree |
|---|---:|
| human → CHM13 (2,138 `primary_transcript` rows) | **0** |
| the chr22 example | **0** |
| bee, human → zebrafish, dog → cat | **0** |
| **rice** | **17** (`gene`) |
| **arabidopsis** | **7** (`gene`) |

A whole-genome human → CHM13 four-arm A/B confirms it: `-t 1` and `-t 8` agree
on the pre-fix build, and the fix changes nothing there. The count that matters
was measured against the reference, and the reference is not what is queried.

## What does carry it: RefSeq's organellar convention

A chloroplast gene lists the same exons **twice** — once directly under the
gene and once under its mRNA, with different IDs:

```
gene  gene-OrsajCp058
  exon id-OrsajCp058-1    34105781-34105788   <- the gene's own copy
  exon id-OrsajCp058-2    34106537-34107011
  mRNA rna-OrsajCp058
    exon exon-OrsajCp058-1  34105781-34105788 <- the same two coordinates
    exon exon-OrsajCp058-2  34106537-34107011
```

The gene has direct level-1 exons, so `process_liftoff` classifies it as
transcript-shaped and calls `lifton_add_trans_exon_cds` on it. The recursive
query then returned **4** exons for a 2-exon gene — each coordinate pair twice.

So the level-1 answer is the correct one and the **serial** path was the buggy
arm (the opposite of what the audit note assumed), and its extra exons were a
duplicate-coordinate pair: the same defect cycle 3's P2 is about, produced by
this bug. Rice's 30 overlapping-exon transcripts split 17 from here and 13 from
P2's own two causes.

## Why nothing caught it

`test_native_matrix`, `test_parallelism_matrix` and `test_fresh_parallel_step7`
all import the same `integration_workspace`: one gene, one transcript, leaf
exons. In that shape a level-1 query and a recursive one are equal by
construction, so all 24 cells agree no matter which is used.

The one gate that could have caught it, `benchmarks/release_gates.sh`, does a
real chr22 `-t 1` vs `-t 8` `cmp` — but it is a manual script, CI's chr22
example passes no `-t` at all, and (per the table above) the chr22 example does
not carry the shape anyway.

`tests/test_materialise_query_budget.py` looked like it pinned the query shape.
Its `_CountingDB.children` ignored `level` and answered every exon query with
the same list, so the assertion was vacuous. It now asserts `level == 1`.

## What is pinned now

`tests/test_parallelism_matrix.py::TestNestedHierarchyThreadAgreement` drives
the whole pipeline over a 3-level fixture at `-t 1` and `-t 4` and requires the
bytes to match, and separately requires the parent to emit only its own exon —
so agreeing on the *wrong* answer cannot pass either.
`tests/test_locus_materialise.py::TestNestedExonAgreement` pins the coupling at
the unit level by spying on the query `lifton_add_trans_exon_cds` actually
sends, rather than restating it.

All four fail on the pre-fix build and pass after.

## Also here: the depth guard dropped subtrees silently

`_walk_and_cache_features` warned that exceeding `max_depth=8` "will surface as
KeyError on the proxy". It does not. `_LFeatureDbProxy.children` answers an
un-cached id with an **empty iterator**, so the runtime reads the feature as
childless: its own row is still emitted and its entire subtree — transcripts,
exons, CDS — disappears with no error, no failure record and no count. The
batched twin exits on its loop condition and logged nothing at all.

That is the same shape as the `-copies` loss that shipped in every release up
to v1.0.11. Both walkers now record a `drop_ledger` class
(`hierarchy_depth_exceeded`) and the warning says what actually happens.

## Lesson

The audit located the divergence in the code correctly and then counted its
reach in the wrong database. A query's blast radius is whatever the *queried*
object contains, which here is a pipeline artifact two steps removed from the
input the count was taken from. Measuring the reference felt like measuring the
input; it was measuring something the query never sees.
