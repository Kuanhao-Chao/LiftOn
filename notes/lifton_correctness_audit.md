# LiftOn full-repository correctness audit (2026-07, `devel` @ v1.0.10.1)

Three parallel code-verified audits over the whole repository — the core algorithmic
path; the ten newest subsystems; the test/perf/memory landscape — against `devel` at
`e7a2aeb`. Baseline at the time: **1382 passed, 2 skipped** (1384 collected). Every
finding below was read in code; the ones marked *reproduced* were executed.

The recurring theme: **the suite is green but structurally blind.** Every synthetic
fixture is gene-only / coding-only / 600 bp and 2-level, and the newest parallel Step-8
dispatcher was absent from the load-bearing 24-cell matrix. That is the same shape as the
Iter-20/21/24 bugs, which all shipped past hundreds of green tests and surfaced only on
full genomes.

---

## FIXED in this pass

| # | Sev | Issue | Commit |
|---|-----|-------|--------|
| 1 | CRITICAL | `lifton_utils.LiftOn_miniprot_alignment` dereferenced a `None` parasail alignment (miniprot mRNA with no CDS children) → `AttributeError` → **the entire Liftoff gene locus was dropped** despite a fine DNA lift. Landed with the paired ordering fix: the chaining whole-side fallbacks returned CDS ascending while `update_cds_list`'s `reverse()` expects the descending order the main path emits, so fixing the crash alone would have activated a latent corruption. | `c8ef69b` |
| 2 | CRITICAL | Any per-locus failure blocked publication → **exit 2, only `*.partial.gff3`**; one bad gene in a 60k-gene genome produced no output. Now publishes + reports `partial_success`; only failures that mean the output cannot be trusted block. `--strict-completeness` restores the old behaviour. | `46b2fda` |
| 3 | CRITICAL | The **mandatory** structural validator used a closed `TRANSCRIPT_TYPES` allowlist while the lift auto-detects types and preserves the reference's middle level verbatim → LiftOn lifted an Ensembl `pseudogene → pseudogenic_transcript → exon` correctly, then **rejected its own output** and exited 2. Recognition is now shape-based (`is_transcript_type`), applied to both validators. | `1a4378a` |
| 4 | HIGH | `exec_miniprot` preflighted miniprot *before* the `-M` short-circuit, so the documented cached `-L`/`-M` workflow exited 1 without miniprot installed. Also hardened `check_miniprot_installed` (suppress output — the version banner corrupted `-o stdout` GFF3 — and check the return code). | `43785a9` |
| 5 | HIGH | `update_cds_list` Case 1 re-emitted the CDS-bearing exon for **every** downstream exon when `optimize=False` (`--legacy-merge` and every direct caller). Measured on the pre-fix baseline: a 7-exon transcript emitted **12 exons**, all but one a duplicate. | `e375c93` |
| ~~6~~ | — | **REVERTED — see the NO-GO section below.** Restoring an ambiguous chunk's CDS was structurally right but net-negative on real data. | `f63324b` → reverted in `e1470c8` |
| 7 | HIGH | Input fingerprinting (a multi-GB SHA-256 on a non-daemon executor the interpreter joins at exit) started **before** the aligner preflight, so a typo'd path or missing binary appeared to hang. | `36cbbec` |
| 8 | HIGH | `CdsIdAllocator` claims were never released, so a **rejected** rescue candidate permanently burned stable CDS IDs and the rightful owner was pushed to a `-1` variant — output depending on how many candidates happened to be rejected. New `tentative()` savepoint, wired into all three stage-then-discard paths. | `eb265b2` |
| 9 | HIGH | 3-level hierarchies (`gene → primary_transcript → miRNA → exon`) silently lost the transcript **and all its exons/CDS**: ids were resolved from the intermediate feature, and `LiftOn_FEATURE` lacked `add_transcript`. Added the delegates, fixed id resolution, made the skip log, and fixed `_render_feature`, which recursed blindly into `.features`. | `96a3e9f` |
| P1 | PERF | Scalar materialise walked exon/CDS **leaves**, ~4 redundant SQLite queries each (millions on human RefSeq). Mirrored the batched sibling's container/terminal split. | `da8d610` |
| P2 | PERF | `get_id_fraction` per-character Python loops ran ~6–12× per transcript in the Step-7 hot loop, holding the GIL. Vectorized: **2.4–2.8× protein, 1.6–1.8× DNA**, 126 fuzz cases against the literal originals. | `b1009a9` |
| P3 | PERF | `Step7StateCoordinator._logical_intervals` was never pruned, so every `overlap()` (two per transcript/candidate pair) rescanned a list that grew all run, under one lock — O(n²) capping `--threads N`. Committed intervals are provable duplicates of the real tree, so pruning them is a no-op for results. | `c130c59` |
| V1 | GAP | Parallel Step 8 had **no** serial-vs-threaded byte gate (the 24-cell fixtures yield ~0 Step-8 candidates). Added one with a fixture that genuinely fires, plus a guard so it cannot go vacuous. | `6126700` |

### Measured effect

Drosophila, `-t8 --locus-pipeline`, cached `-L`/`-M`, 2 repeats, arms alternated,
sibling DBs rebuilt per arm (base = `e7a2aeb` in a detached worktree):

| arm | Step-7 wall | peak RSS |
|---|---:|---:|
| base r1 / r2 | 158.9 s / 164.3 s (mean **161.6 s**) | 1.06 GB / 1.10 GB |
| head r1 / r2 | 132.4 s / 122.5 s (mean **127.4 s**) | 1.10 GB / 1.07 GB |

**Step-7 wall −21.2 % (1.27×)**, consistent across repeats; peak RSS flat (within noise).
Total wall on r1: 3:17 → 2:51 (−13.4 %).

The 12-line output delta observed in this A/B came entirely from the ambiguous-chunk
change, which was subsequently **reverted** (below). **With the final tree, the batch is
byte-identical on drosophila and exactly identity-neutral on mammalian data.**

---

## NO-GO — reverted after full-genome validation

**Restoring an ambiguous chunk's CDS** (`f63324b`, reverted in `e1470c8`).

A chunk where neither aligner scores a matching column returned `[]`, dropping its CDS
blocks and leaving a hole in the merged list — structurally wrong, so the fix emitted
Liftoff's blocks instead (the documented tie-break). The mammalian dog→cat validation
(`-t8`, cached `-L`/`-M`) showed the hole was frequently scoring *higher* protein identity
than the completed model, because restoring the block reintroduces its mismatches:

| arm | improved | regressed | mean Δ protein identity |
|---|---:|---:|---:|
| base → HEAD, all fixes | 0 | **8** | **−0.00038106** |
| base → HEAD minus this commit | 0 | **0** | **+0.00000000** |

Individual losses were large (0.844→0.615, 0.849→0.616, 0.840→0.609). Only 1 of the 8 was
a status flip, so this was not about which candidate wins: the merge **candidate itself**
changed, got worse, and still beat pure Liftoff, so it was kept.

Two lessons worth keeping:

1. **The gate would not have caught it.** −0.00038 is inside
   `scripts/benchmark_gate.py`'s ≤0.005 tolerance, and the drosophila subset showed only
   3 benign transcripts with *identical* coordinates and *identical* protein identity.
   Only the mammalian run exposed the cost — the Iteration-20/21/24 lesson again, now
   also true for accuracy and not just crashes.
2. **Structural correctness and the scoring metric can disagree.** Reverting matches the
   established discipline for net-negative accuracy changes (Iterations 4, 9 and 15;
   Iter-9 was "4 improved / 59 regressed", this was 0/8).

The paired CDS-ordering fix (`_cds_in_chain_order`, from `c8ef69b`) is independent and
stays. To revisit, the real question — is a hole better than a mismatched block? — needs a
divergence-ladder A/B judged on ORF-validity and the structural metrics
(`benchmarks/compare/structural_rescore.py`), not on protein identity alone.

---

## DEFERRED — verified findings not addressed in this pass

Kept here so nothing is lost. Each has a rationale for deferring.

### Correctness / robustness
- **Proxy `NotImplementedError` and depth-8 truncation are per-locus fatal**
  (`locus_pipeline.py` `_LFeatureDbProxy.children`, `_walk_and_cache_features`). The
  user-facing half — a whole-run abort — is **already fixed** by #2 above. The residual is
  a developer guardrail that only fires if a future refactor adds a new `children()`
  signature; wiring a real-DB fallback risks using a parent SQLite connection from a
  worker, the exact hazard the proxy architecture exists to prevent. Fix by giving the
  worker's *thread-local* DB to the proxy, only on the `factory.viable` path.
- **`validate_gff3_file` is not memory-bounded** (`--validate-output` materialises a
  `GFF3Record` per row; ~4 M rows on a human-scale output). Rewrite the hierarchy /
  containment / phase checks as a streaming pass over the contiguous top-level blocks
  LiftOn already guarantees.
- **`validate_gff3_structure` under-reports counts**: it never populates
  `severity_totals`/`issue_totals`, so `manifest.record_validation(..., errors=...)` is
  capped at the per-check limit and `is_valid`'s first clause is inert.
- **gffbase root-scan ordering is not total**: `_order_clause(None)` emits
  `file_order ASC` with no tiebreak, so rows with NULL `file_order` (inferred roots) can
  order differently between the serial and parallel scans → different submission indices.
  One-line fix (`, id ASC`), plus a matrix cell with inferred roots.
- **Read-only reference directories hard-fail**: `annotation_cache.cache_lock` opens a
  `.lock` next to the DB, so a shared/read-only reference mount (a common HPC layout)
  raises an uncaught `PermissionError` where the DB used to be reused read-only.
- **Durability gap**: `stats.print_report`'s `finally` does a bare `fw.close()` on the
  transaction's staging stream, so `OutputTransaction.close()` skips its `fsync` and
  `commit()` publishes data that was never flushed to disk.
- **Miniprot child proxies ignore `order_by=None`** and always serve start-sorted rows,
  while the real backends return file order — currently masked because miniprot emits
  ascending.
- **Three-source strand ambiguity in `Lifton_TRANS`** (`get_coding_trans_seq` keys on
  `exons[0].entry.strand`, `__update_cds_boundary` on `self.entry.strand`,
  `get_coding_seq` per-CDS). Enforce one source of truth at `add_exon`/`add_cds` time.
- **`start_lost` is effectively untestable as written** (`variants.py`): it compares the
  first three bases of the full spliced transcript, which is 5′UTR for most transcripts,
  so a genuine start loss is not flagged — and it gates ORF rescue. The CDS span is now
  computed in `orf_search_protein` (added for GH #46) and can be reused.
- **`Lifton_TRANS.add_cds` attaches one CDS to every overlapping exon** and mutates the
  caller's `Parent`; "a CDS lies in exactly one exon" is assumed, never checked.
- **Case-1 merged exon inherits the DOWNSTREAM exon's ID** (masked today because
  `normalize_containment` renumbers on collision — see `test_cds_ids.py`).
- **`get_id_fraction` guard drift**: `get_partial_id_fraction` tests `total_length == 0`
  where its sibling uses `<= 0`. Unreachable today; left alone deliberately because
  changing it would change output in the negative case.

### Output quality
- **Merged CDS lines lose their descriptive attributes.** Measured on drosophila: of the
  CDS beneath `status=LiftOn_chaining_algorithm` transcripts, **35,145 carry only
  `ID`/`Parent` while 1 keeps the rich set** — `Dbxref`, `Name`, `gene`, `gbkey`,
  `locus_tag`, `product`, `protein_id` are dropped whenever the merge fires, because
  `create_lifton_entries` assigns a shared `parent_attrs` and `add_lifton_cds` then resets
  attributes to `{Parent}`. Transcripts that stay on the pure-Liftoff path keep everything,
  so a single output is inconsistent and most CDS rows are impoverished. This is
  long-standing behaviour (not introduced by this batch — it was measured on the *base*
  commit) but it is a real usability regression against the reference annotation and a good
  candidate for the next pass: carry the source CDS attributes through the merge the way
  the CDS **ID** is now carried (`Lifton_CDS._source_id`).

### Performance / memory
- **`scan_annotation` runs full ID bookkeeping + NCBI per-line validation on the derived
  Liftoff/miniprot intermediates**, not just the reference — two extra full-file passes
  and multi-hundred-MB ID dicts per run. Add a `fingerprint_only` scan level.
- **`AnnotationScanResult` is retained for the whole run**; only `cds_namespace_ids` and
  `copy_suffix_ids` (reference only) are needed after the allocator is seeded.
- **`trunc_ref_proteins` builds a whole dict only to take `len()`** (~5.6 s at human
  scale). Replace with a counter.
- **The gffbase adapter re-hashes the reference up to 4×** because `_expected_manifest`
  omits `source=` (the gffutils path passes it correctly).
- **`align.py` computes the coding sequence twice per transcript**
  (`get_coding_seq` then `get_coding_trans_seq`), paying two full sets of pyfaidx fetches
  in the Step-7 hot loop.
- **`orf_search_protein` re-reads the reference protein/transcript FASTA records on every
  call** — up to 3× per merge-firing transcript.
- **`__find_orfs` accumulates a write-only `orf_seq` string** (4–9 % of the scan) and its
  outer loop is worst-case quadratic when an ATG has no downstream stop.
- **Step 5 (2 gffutils DB builds, ~15 % of a cached-aligner run) is NOT a threading
  target** — that is the documented Iteration-11 NO-GO (GIL-bound; measured 8–12 % slower).

### Measurement
- **`--measure_time` uses `time.process_time()`** for all steps, so it cannot see the
  miniprot subprocess and double-counts parallel threads. The correct per-phase wall clock
  already exists in `run_manifest.json`; `time.txt` should be sourced from it.

### Documentation
- **`CLAUDE.md` is stale**: it describes 524–805 tests and a ~649-line `lifton.py`
  (actually 1384+ tests and 1850 lines) and does not mention nine shipped modules
  (`run_manifest`, `output_transaction`, `annotation_cache`, `annotation_validator`,
  `cds_id_allocator`, `cross_locus_rescue`, `miniprot_pipeline`, `variants`, `coreutils`).

---

## Verified sound (checked, no action needed)

Recorded so they are not re-audited: the windowed aligner's parasail-convention CIGAR and
`adjust_cdss_protein_boundary`; CDS phase arithmetic and `__iterate_exons_update_cds` on
both strands; `__find_orfs` frame handling; the `find_variants` early-return ladder and the
GH #46 coding-scoped frameshift window; the best-of-outcome snapshot/restore protocol
(including `_optimize_fast` skip byte-neutrality); `chaining_algorithm` termination;
`_ordered_bounded_map` progress guarantees; the ordered copy-number gate's release-on-
exception; SQLite thread affinity in `_ThreadLocalCtxFactory`; worker/parent state
separation in Step 7; deferring delta commits until after a successful write; the Step-8
monotonic pre-overlap elision; `_has_unencoded_reserved`'s `isdisjoint` optimisation;
`source_fingerprint` change detection; `atomic_write_json`; the `OutputTransaction` state
machine; `CdsIdAllocator._copy_bases` termination; and the DB publish-then-manifest
ordering.
