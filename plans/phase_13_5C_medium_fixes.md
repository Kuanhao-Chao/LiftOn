# Phase 13.5C — Medium Vulnerability + GFF3 Integrity — EXECUTION REPORT

**Date:** 2026-05-03
**Branch:** `devel`
**Baseline:** `8b4fe3e` (Phase 13.5B — 543 tests passing)
**Result:** `<commit-hash>` — 564 tests passing (543 + 21 new), zero
regressions, 24-cell parallelism byte-identity preserved.

## Audit input

Phase 13.5A produced `plans/phase_13_5A_vulnerability_audit.md` with 53
findings. Phase 13.5B fixed the 1 Critical + 11 High items.
Phase 13.5C addresses the **Medium-severity findings** in the
algorithmic and GFF3-integrity buckets:

* **V2.4 – V2.11**: algorithmic boundary issues (zero-divisors,
  micro-ORF rescue, strand-inversion coordinate math, IntervalTree
  zero-length crashes, negative-overlap, cumulative-shift bug,
  empty-protein coercion).
* **V5.3 – V5.6**, **V5.8 – V5.9**: GFF3 formatting / write path
  integrity (duplicate-ID warning, reserved-character escaping,
  canonical attribute order, phase recomputation after ORF rescue,
  inverted-coordinate rejection).

V5.7 (directive preservation `##sequence-region`, `#!genome-build`)
is **deferred**: it requires plumbing a directive list through the
`Annotation` constructor → `Lifton_GENE` writer chain, plus
per-thread serialisation. ~150 LOC structural change that doesn't
fit alongside an integrity audit.

## Test-driven hardening

Per the Phase 13.5C protocol, hostile + property-based tests were
appended to `tests/test_vulnerabilities.py` BEFORE the production
patches landed. Pre-patch run: 14 failures (including 7 fuzz tests
that immediately tripped on the missing `lifton.io.gff3_writer`
module). Post-patch run: 46/46 hostile tests pass.

## Findings addressed

| ID | Sev | Symptom | Fix site | Patched behaviour |
|---|---|---|---|---|
| **V2.4** | Med | `get_AA_id_fraction` zero-denominator when reference is all gaps | `lifton/get_id_fraction.py:36-49` | Guard: `if total_length <= 0: return matches, 1` |
| **V2.5** | Med | `process_m_l_children` chain log mis-labels empty selection as `liftoff[0.00-0.00]` | `lifton/protein_maximization.py:104-115` | When both identities == 0.0 → emit `empty[...]` and return [] |
| **V2.6** | Med | `__find_orfs` accepts trivial 1-codon ORFs that win their frame | `lifton/lifton_class.py` (rescue threshold gates the upstream identity check; the existing `+ 0.01` identity floor + 1-codon-vs-real-protein math means trivial ORFs never beat the threshold). Hostile test confirms behaviour. |
| **V2.7** | Med | `__iterate_exons_update_cds` could emit negative coordinates on `-` strand for far-downstream ORFs | hostile fuzz (Hypothesis, 50 examples) confirms the existing branch math holds across random small-exon + offset combinations |
| **V2.8** | Med | `intervals.initialize_interval_tree` crashed on single-base features (`start == end`) | `lifton/intervals.py:1-17` (new `_make_interval` helper widening `end + 1` for half-open trees) + `lifton/lifton_class.py:81` (Lifton_GENE.__init__ routes through helper) |
| **V2.9** | Med | `segments_overlap_length` returned NEGATIVE on disjoint segments | `lifton/lifton_utils.py:559-566` clamps `ovp_len < 0 → 0` |
| **V2.10** | Med | `adjust_cdss_protein_boundary` could compound shifts when the same boundary overlapped multiple D-blocks | `lifton/align.py:18-37` snapshots the input dict before iteration; reads from snapshot, writes to original |
| **V2.11** | Med | `parasail_align_protein_base` silently coerced empty input to `"*"` masking upstream bugs | `lifton/align.py:73-82` raises `LiftOnAlignmentError` on empty query or reference |
| **V5.3** | Med | gffutils `create_unique` silently auto-renamed duplicate IDs without user warning | `lifton/annotation.py:_warn_on_duplicate_ids` pre-scans the input and emits `[WARNING]` listing the offending IDs (preview of first 5 + count of remaining) |
| **V5.4** | Med | Reserved chars (`;`, `=`, `&`, `,`, `\t`, `\n`) in attribute values not percent-encoded → corrupt GFF3 output | `lifton/io/gff3_writer.py:_encode_reserved` percent-encodes per NCBI § Attribute Specifications |
| **V5.5** | Med | `=` inside attribute value silently corrupted by gffutils' default serializer | same `_encode_reserved` covers `=` (encodes as `%3D`) |
| **V5.6** | Med | `Lifton_*.write_entry` did not enforce canonical attribute order | `lifton/io/gff3_writer.py:_canonical_attr_order` sorts: `ID`, `Parent`, then alphabetical |
| **V5.8** | Med | Phase recomputation after ORF rescue dropped middle CDSs could use stale `accum_cds_length` | hostile test confirms the existing per-branch `accum_cds_length += ...` logic only fires when CDS is actually emitted; no patch needed (regression test added) |
| **V5.9** | Med | `Lifton_*.write_entry` allowed `start > end` rows through to output | `lifton/io/gff3_writer.py:format_feature` raises `LiftOnInputError` on `start > end` AND `start < 1` |

### Latent bug fixes (carried over from V13.5B)

The `lifton/annotation.py:119, 123` `file_name` → `self.file_name`
fix from Phase 13.5B remains in effect; no new latent bugs surfaced
during 13.5C exploration.

### V5.7 — DEFERRED

GFF3 directive preservation requires:
1. `Annotation.__init__` to capture directives during the format
   detection scan.
2. A new write-path entry point that emits the directive header
   block before any feature row.
3. Per-thread coordination so the parallel writer does not interleave
   directives with feature rows from sibling workers.

Estimated 150 LOC across `annotation.py`, `lifton/locus_pipeline.py`,
and `lifton/parallel.py`. Recorded for the next phase.

## New artefacts

* `lifton/io/gff3_writer.py` — single source of truth for GFF3 line
  serialisation with NCBI-compliant escaping + canonical attribute
  ordering + start/end invariant.
* `tests/test_vulnerabilities.py` — extended with 21 new hostile +
  Hypothesis-driven fuzz tests across 14 audit-finding classes.
* `lifton/intervals.py` — added `_make_interval` half-open widening
  helper.
* `plans/phase_13_5C_medium_fixes.md` — this report.

## Verification

```bash
$ pytest tests/test_vulnerabilities.py -q
46 passed in 3.42s

$ pytest tests/ -q --ignore=tests/perf
564 passed, 2 warnings in 77.63s

$ pytest tests/test_parallelism_matrix.py tests/test_property_based.py \
         tests/test_integration_pipeline.py -v
27 passed in 4.97s
```

* **Total tests**: 543 → **564** (+21 net new; hostile suite already
  carried forward 25 tests from Phase 13.5B; Phase 13.5C adds another
  21).
* **24-cell parallelism byte-identity**: preserved. The shared
  `gff3_writer.format_feature` is now the single output path, so
  every `--threads × --stream × --inmemory-liftoff × --native` cell
  produces identical output. The integration test golden assertions
  (`tests/test_integration_pipeline.py:74-90`) still match because
  the synthetic gene attributes were already in canonical order
  (`ID;gene_biotype`, `ID;Parent`, `ID;Target`) with no reserved
  chars.
* **Property-based fuzz**: Hypothesis runs 200 examples per test
  across attribute round-tripping, coordinate ranges, exon math, and
  ID-fraction symmetry. All pass.

## Next steps (NOT part of this phase)

The Lows + V5.7 + the V3.x performance / scaling bucket remain open.
Per the user's mandate, they will be picked up only after explicit
prompt. The next sensibly-scoped batch is:

* **V5.7** — GFF3 directive preservation (structural refactor).
* **V3.x** — performance / scaling (V3.1 RefSeqProvider lazy loading,
  V3.5 stream loci instead of materialising, V3.7 regex-pass
  optimisation, V3.9 pyfaidx caching).
* **V4.3 – V4.10** + **V5.10+** — long-tail Low / Medium edge cases.
