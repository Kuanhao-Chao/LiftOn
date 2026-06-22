## Step-3 gffutils query-collapse micro-bench — Iteration 18

Step 3 issued up to **3** `gffutils.children()` queries per feature (exon; CDS-tuple; recursive-descent fallback). Iteration 18 collapses them to **one** `order_by='start'` query + an in-Python featuretype partition (byte-neutral: `merge_children_intervals` re-sorts the exon/CDS lists before concatenation, so SQL order is discarded for those branches; the recursive-descent branch keeps `order_by='start'`).

This times **only Step 3** on a real reference DB built ONCE and reused. `baseline` = the pre-Iter-18 3-query extractor loaded verbatim from `git show HEAD:lifton/extract_sequence.py`; `branch` = the working-tree 1-query extractor. The `children()` call count is the deterministic (noise-free) mechanism proof; the wall-clock median±IQR is the GO/NO-GO measurement.

**Gate:** byte-identical FASTAs AND median wall reduction > 2× the run-to-run noise band.

| Dataset | byte-ident | children() calls | reduction | Step-3 wall median | speedup | wall-GO |
|---|---|---|---|---|---|---|
| drosophila | yes | 28816→12315 | 2.34× | 6077.931ms→4219.745ms (Δ1858.186ms, ±73.667 noise) | 1.44× (30.57%) | **yes** |
| mouse_to_rat | yes | 10072→4504 | 2.24× | 2804.449ms→2099.811ms (Δ704.638ms, ±14.325 noise) | 1.336× (25.13%) | **yes** |
| rice | yes | 27704→11616 | 2.38× | 5106.079ms→3379.346ms (Δ1726.733ms, ±52.204 noise) | 1.511× (33.82%) | **yes** |

**Interpretation:** the `children()` call-count reduction is exact and deterministic (the optimization's mechanism). The wall-clock gate decides GO/NO-GO: if the median Step-3 reduction does not clearly exceed the noise band on any dataset, the change is byte-correct but unmeasurable → NO-GO + revert, keeping this script + results as the audit trail (the Iter-9/11/13/15 precedent).
