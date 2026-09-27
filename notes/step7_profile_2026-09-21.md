# Where Step 7 spends its time now

*2026-09-21. Re-profile of `notes/step7_profile_2026-09.md` on the current
build. Same recipe so the tables compare: `LIFTON_PROFILE_STEP7`, whole genome,
cached `-L`/`-M` so the aligners cannot dominate, `-t 1` (cProfile serialises
Python; the question is **where**, not how fast).*

The previous profile named seven targets. Four have since been addressed —
translation (`transl_table` centralisation), `clone_attributes`, the Step-3
query collapse, and the windowed aligner's anchor construction. This is what
the same measurement says afterwards.

## Rice, whole genome

Step-7 dispatch 667 s (was 730 s), 130,900 translations.

| # | site | now | previous | note |
|---|---|---:|---:|---|
| 1 | `parasail.nw_trace_scan_sat` | 210.1 s (31 %) | 210 s (29 %) | not addressable — C kernel |
| 2 | **`sqlite3.Cursor.execute`** | **165.5 s (25 %)** | 151 s (21 %) | **grew, and is now the top addressable cost** |
| 3 | `sqlite3.Cursor.fetchone` | 17.4 s | — | |
| 4 | `gffutils.interface._relation` | 13.2 s (cum 101.5 s) | — | 1.48 M generator steps |
| 5 | `coding.py` translate listcomp | **10.5 s** | **66 s** | the `transl_table` flattening: 6× |
| 6 | simplejson `raw_decode` | 9.6 s (cum 18.4 s) | — | 2.86 M attribute decodes |
| 7 | `pyfaidx.__getitem__` | 7.8 s (cum 31.7 s) | — | 1.74 M sequence fetches |
| 8 | `gff3_writer.format_attributes` | 6.6 s (cum 25.3 s) | 38 s | |
| — | `encode_attribute_value` | **2.6 s** | **38 s** | 5.84 M calls, now nearly free |

The two targets the previous note recommended are gone from the top of the
profile. Translation fell 66 s → 10.5 s and attribute encoding 38 s → 2.6 s.

What they were hiding is the database.

## Dog → cat, whole genome

Step-7 dispatch **2,848.9 s, down from 3,299 s — −13.7 %** on the same cached
inputs. That is the accumulated effect of the programme's byte-neutral work,
measured rather than projected.

| # | site | now | previous | note |
|---|---|---:|---:|---|
| 1 | `parasail.nw_trace_scan_sat` | 995.5 s | 1,003 s | flat, as expected — C kernel |
| 2 | **`sqlite3.Cursor.execute`** | **274.8 s** | 298 s | same rank as rice |
| 3 | **`lifton_class.__find_orfs`** | **73.7 s** (cum 221.8 s) | *not in the top 7* | newly visible |
| 4 | `windowed_align._unique_anchors` | **61.2 s** | **116 s** | **1.9×** — this cycle's change, on real data |
| 5 | `variants._coding_subalignment` | 58.9 s | — | |
| 6 | `get_protein_reference_length_single` | 57.5 s | 58 s | unchanged; 1.31 M calls |
| 7 | `windowed_align._cigar_from_aln` | 52.9 s | — | |
| 8 | `windowed_align._chain_colinear` | 51.7 s | — | |
| 9 | `coreutils.segments_overlap_length` | 42.4 s | — | 43.6 M calls |
| 10 | `coding.py` translate listcomp | **42.2 s** | **269 s** | **6.4×** — the `transl_table` flattening |
| — | `encode_attribute_value` | out of the top 12 | 111 s | |
| — | `clone_attributes` | out of the top 12 | 148 s | |

Two of the previous note's recommendations are confirmed on real mammalian
data, not just in a microbenchmark: translation is 6.4× cheaper and the
windowed aligner's anchor construction 1.9×, the latter matching the top of the
range its A/B predicted (1.3–1.9×, widest at high divergence — and dog → cat is
the divergent case).

The windowed group as a whole is 165.8 s, down from 320 s; `_cigar_from_aln`
and `_chain_colinear` are the untouched half of it.

**`__find_orfs` is new to the top of the profile** — 73.7 s of self time and
221.8 s cumulative, 8 % of dispatch. It was always there; it became visible
when the costs above it fell. It is absent from rice, so like the windowed
aligner it is a mammalian and divergent-transfer cost. An earlier cycle
measured an ORF-selection change here as **net-negative** and reverted it
(`notes/`, Iteration 9), so any return to this code needs the A/B before the
change, not after.

## The finding: it is query volume, and the rows each query builds

It is the top addressable cost on **both** regimes, which the previous
profile's ranking did not show — there it sat behind translation, attribute
encoding and `clone_attributes`, all of which have since been fixed.

Grouping the gffutils layer on rice:

| | |
|---|---:|
| `sqlite3.Cursor.execute` | 165.5 s |
| `fetchone` | 17.4 s |
| `_relation` (generator) | 13.2 s |
| `Feature.__init__` + `_feature_returner` | 9.0 s |
| attribute construction + `_unjsonify` + simplejson | 23.7 s |
| **total** | **≈ 229 s, 34 % of Step 7** |

476,000 SQL executions build 1.43 M `Feature` objects, and every one of them
JSON-decodes its attribute column.

Counted directly on the repository's own chr22 example, by signature:

| featuretype | level | order_by | calls | rows | rows/call |
|---|---|---|---:|---:|---:|
| `exon` | 1 | start | 10,042 | 81,167 | 8.08 |
| `None` | 1 | — | 9,751 | 7,660 | **0.79** |
| `('CDS','stop_codon')` | — | start | 6,609 | 64,320 | 9.73 |
| `CDS` | — | start | 4,952 | 64,950 | 13.12 |
| `None` | 1 | start | 4,952 | 76,688 | 15.49 |
| `CDS` | — | — | 2,332 | 28,158 | 12.07 |

38,638 calls in total. Two things stand out:

* **A quarter of all queries return less than one row.** `children(level=1)`
  with 0.79 rows per call is the container probe, asked of features that turn
  out to be leaves. It costs a round-trip to learn nothing.
* **A transcript's CDS children are fetched under three different
  signatures** — `('CDS','stop_codon')` ordered, `CDS` ordered, and `CDS`
  unordered — 13,893 calls returning 157,428 rows between them.

## What this rules in, and what it rules out

It rules **out** another pass at translation, attribute encoding or the
aligner's anchor phase: those are done, and the profile says so.

It rules **in** one thing, and it is the same mechanism found independently in
two other places this cycle:

> **gffutils has no batched children API.** `children_batched_features` exists
> only on the gffbase backend. `HierarchyBatchLoader` checks for it
> (`locus_pipeline.py:772`) and falls back to row-at-a-time on gffutils, which
> is the default. Step 8 reports `step8_child_batch_calls: 0` and
> `step8_child_scalar_materializations: 132,781` on human → zebrafish, and the
> rescue's isoform prefetch spends 117 s on ~112,000 scalar round-trips that
> it could batch.

So Step 7, Step 8 and the rescue all issue one query per feature for the same
reason, and a batched `children()` for the gffutils backend would be a single
change benefiting all three.

## Decision for this cycle: NO-GO, and why

Not attempted here. It is a new query path on the default backend, under the
24-cell byte-identity contract, touching the three hottest loops in the
program. That deserves its own cycle with its own gate — a batched result must
return the same rows in the same order as the scalar one, per signature, which
is provable but needs to be proved rather than assumed.

Recorded as the top addressable Step-7 target, with the mechanism measured
rather than guessed. The previous note's warning applies to this one too: this
profiler has contradicted a guessed ranking three times, so the next cycle
should re-measure before acting, not act on this table alone.
