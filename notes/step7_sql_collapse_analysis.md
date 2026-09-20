# The reverted Step-7 SQL collapse: what a second attempt has to handle

`notes/lifton_correctness_audit.md` records a collapse that was implemented,
verified byte-identical on both drosophila anchors, measured at −20.5 % of
Step-7 SQL statements, and then reverted because the 24-cell matrix failed:
`children(locus, order_by='start')` is a signature `_LFeatureDbProxy` does not
cache. The note says landing it needs an `all_children_full` entry on the
prefetch, populated in both walkers and served by the proxy.

That prescription is right as far as it goes. Two things it does not say, found
by reading the code before starting:

## 1. The threaded path does not gain what the note implies

The saving is in `run_liftoff.lifton_add_trans_exon_cds`, which issues two
filtered `children()` calls per transcript. Under `--threads N` the runtime
never touches the database: it reads `_LFeatureDbProxy`, and the queries are
issued instead by `_walk_and_cache_features` on the parent thread.

That walker's terminal branch already issues exactly two queries per transcript
— one to classify (`featuretype='exon', level=1`) and one for CDS/stop-codon —
and the classification query cannot be collapsed away, because container-vs-
terminal is a level-1 question and the collapsed query is recursive.

So the collapse removes one query per transcript on the **serial** path, and
none on the threaded path unless the walker is restructured too. Since
`--threads N > 1` implies the locus pipeline, that is most real runs. The
−20.5 % figure was measured serially and should not be quoted for a threaded
run.

## 2. The proxy currently answers a recursive question with a level-1 answer

`_FeaturePreFetch.exon_children_full` is documented as
`children(featuretype='exon', order_by='start')` — recursive, no level — but
the walker populates it as `list(entry.exon_children_l1)`, the level-1 result,
justified by the leaf-exon convention (`locus_pipeline.py:953-955`).

That is fine while nothing exposes the difference. Building `all_children_full`
from a real unfiltered recursive query does expose it: partitioning that result
yields the **recursive** exon set, where the proxy today yields the level-1
set. On a standard annotation they are equal. On a hierarchy with a non-leaf
exon they are not, and the threaded path would change behaviour — in the
direction of matching the serial path, which is arguably a latent divergence
being fixed, but it is still an output change arriving through a commit whose
stated purpose is "fewer queries".

A second attempt therefore has to decide explicitly which of the two answers is
correct, and verify that decision on a hierarchy that distinguishes them. The
24-cell fixtures are all two-level and cannot.

## What this is worth

On rice, `sqlite3.Cursor.execute` is 151 s of a 730 s Step 7. The second
collapse alone moved statements 65,901 → 57,977 on the drosophila subset,
about −12 %. Applied to the serial-path share, that is a small single-digit
percentage of Step 7 — worth having, not worth landing without the hierarchy
fixture above.

Recommended order for the next attempt:

1. Add a three-level fixture with a non-leaf exon and assert serial and
   threaded agree on it *before* changing anything. It may already fail.
2. Decide and document which exon set is correct.
3. Then collapse, with the drosophila anchor `cmp` and the 24-cell matrix.
