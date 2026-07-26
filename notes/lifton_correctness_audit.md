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

## Follow-up pass (2026-07-25) — what was taken off this list

The next batch closed the top output-quality item and four robustness items, and
replaced the guessed performance ranking with a measured one. Resolved entries
are struck from the lists below; the reasoning that produced them is kept.

- **Merged CDS attribute loss — FIXED.** The measurement grew when re-taken on a
  full output rather than a subset: **151,277 CDS rows under
  `status=LiftOn_chaining_algorithm`, of which 21 carried descriptive
  attributes**, against 7,127 of 7,142 on the pure-Liftoff path — and all 28,051
  merged mRNA rows kept theirs, confirming an internal inconsistency rather than
  a policy. A transcript-level template captured in `Lifton_TRANS.add_cds` now
  feeds both rebuild paths. Proven on the eight-cell divergence ladder
  (`benchmarks/compare/cds_attr_parity_ab.{py,json,md}`): CDS coverage 0.06–0.63
  → 1.0, columns 1–8 byte-identical on every row, 0 attributes lost or
  rewritten, `protein_identity` unchanged, validity flat, +12–43% file size.
  *Known caveat, accepted for parity:* a handful of reference CDS attributes
  describe the **reference** feature rather than the lifted one (`partial`,
  `end_range`, `exception`, `transl_except` — 900/699/494/59 rows on the
  drosophila anchor). Carrying them can therefore be slightly stale after a
  merge that moved a boundary. The pure-Liftoff path has always carried them, so
  emitting them is the *consistent* choice; treating them specially would mean
  diverging from the reference in a second, undocumented way. Revisit only with
  a rule that applies to both paths.
- **`validate_gff3_structure` under-reporting — FIXED**, including the inert
  `is_valid` clause. The deep `validate_gff3_file` / `gff3-validate` CLI keeps
  its documented per-check cap on purpose: changing that number would silently
  re-baseline every committed benchmark validity figure.
- **`stats.print_report` fsync bypass — FIXED** (close through the transaction).
- **gffbase `_order_clause(None)` ordering — FIXED** (`, id ASC`).
- **`annotation_cache.cache_lock` on read-only references — FIXED**
  (temp-dir fallback, then unlocked with a warning).
- **`trunc_ref_proteins` dict built for its `len()` — FIXED**
  (`count_truncated_proteins`).

### The Step-7 profile, and what it changed

A `LIFTON_PROFILE_STEP7` cProfile gate (output-neutral) was added and run at
`-t1` on drosophila and dog→cat. It contradicted parts of the guessed ranking
above, which is the point:

| rank | drosophila `tottime` | note |
|---|---|---|
| 1 | `parasail nw_trace_scan_sat` 132.1 s | the kernel; irreducible |
| 2 | `copy.deepcopy` 36.0 s / 73.6 s cum, **36.2 M calls** | **not on the list** |
| 3 | `sqlite3.Cursor.execute` 24.0 s | |
| 4 | `__find_orfs` 18.9 s | on the list |
| 5 | `variants._coding_subalignment` 13.5 s | **not on the list** |
| 8 | `align.py:113` alphabet genexpr 8.0 s + `all` 6.0 s, **71.3 M calls** | **not on the list** |

Fixed from this: the alphabet check (set membership + a C translation table),
`_coding_subalignment` (the selected columns are contiguous, so only the two
boundaries need locating — plus an ungapped fast path), `__find_orfs` (a
write-only accumulator, and a provably equivalent break that removes the
quadratic re-scan), and the discarded coding sequence in `LiftOn_translate`.

**`copy.deepcopy` was addressed in the next pass (2026-07-25, second batch).**
`coreutils.clone_feature` copies a `gffutils.Feature`'s mutable state and shares
its immutable `dialect`: **33.06 µs → 5.00 µs** on a rich CDS. `Lifton_EXON` and
`Lifton_CDS` gained `__deepcopy__` so the hottest caller sped up without
changing its call site. One trap worth remembering: the first cut used
`feature.__dict__`, which is right for `gffutils.Feature` but wrong for the
vendored `lifton.gffbase.feature.Feature` — that one declares `__slots__` and has
no `__dict__`, and it is the class `--inmemory-liftoff` uses. It did not raise;
it silently produced a different lift. The 24-cell matrix caught it, no unit test
did. The original reasoning for deferring is kept below.

**`copy.deepcopy` was deferred in the first pass** even though it ranks second.
Its callers are `run_liftoff._snapshot_merge_state` (20,690 calls, **32.0 s
cumulative** — the best-of-outcome snapshot/restore protocol),
`Lifton_TRANS.get_coding_seq` (15.5 s) and `update_cds_list` (15.3 s). Two
reasons to leave it: the profiled figure is inflated because per-call profiler
overhead falls hardest on 36 M tiny calls, so the real share is smaller than it
looks; and a cheaper clone must reproduce gffutils `Feature` semantics exactly
while `update_cds_list` mixes copied and original exon objects in one list, and
`normalize_containment` later mutates an attribute dict in place. Getting that
subtly wrong changes output. Note the interaction with the CDS-attribute fix:
richer attribute dicts make this snapshot *more* expensive, so it is the
strongest remaining performance target — it just needs its own careful pass.

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
- **`validate_gff3_file` is not memory-bounded — now MEASURED, with the rewrite designed
  and its baselines captured** (2026-07-26). On one 92,675-row output: the deep validator
  peaks at **406 MB** where the streaming `validate_gff3_structure` uses **36 MB** (11.2×,
  and 2.8× slower). On the full dog→cat output (1,815,199 rows) it peaks at
  **5.64 GB**.

  *Design.* `_parse_gff3` materialises every row into `records`, then `_check_hierarchy`,
  `_check_containment`, `_check_cds_phase`, `_check_lifton_attrs` and `_compute_stats` each
  re-walk that list. The five check functions can be kept **unchanged** and fed one
  contiguous top-level block at a time — LiftOn guarantees blocks are contiguous — with
  three things made global so the report does not move: (a) per-check issue counters, so
  `max_issues_per_check` still caps across the file rather than per block; (b) the issue
  lists accumulated per check-type and concatenated in the original check order at the
  end, since today all hierarchy issues precede all containment issues; (c) a first pass
  collecting the **ID set only** (strings, not records), because duplicate-ID detection and
  `_check_hierarchy` rule 1 are inherently cross-block. That last item makes memory
  O(distinct ids) rather than truly O(1) — roughly 200 MB at dog→cat scale, still ~28×
  better than 5.64 GB.

  *Baselines are already captured* for the equivalence gate: the full `ValidationResult`
  (every issue in order, with its message) on the drosophila anchor (92,675 rows, 66
  issues, sha `856609d89b079a22…`) and the full dog→cat output (1,815,199 rows, 176 issues,
  sha `7b290aafe82a0e83…`). These reports feed the committed benchmark validity figures and
  must come back byte-identical.
- ~~**`validate_gff3_structure` under-reports counts**~~ — **FIXED** (2026-07-25).
- ~~**gffbase root-scan ordering is not total**~~ — **FIXED** (`, id ASC`).
- ~~**Read-only reference directories hard-fail**~~ — **FIXED** (temp-dir lock fallback).
- ~~**Durability gap** (bare `fw.close()` skips the transaction `fsync`)~~ — **FIXED**.
- ~~**Miniprot child proxies ignore `order_by=None`**~~ — **FIXED**: the proxy re-sorts
  the cached rows in memory (no database round-trip from a worker thread) and fails
  loudly for an ordering it cannot serve.
- **Three-source strand ambiguity in `Lifton_TRANS`** (`get_coding_trans_seq` keys on
  `exons[0].entry.strand`, `__update_cds_boundary` on `self.entry.strand`,
  `get_coding_seq` per-CDS). Enforce one source of truth at `add_exon`/`add_cds` time.
- **`start_lost` reads the wrong three columns — real, quantified, and NO-GO to change**
  (2026-07-26). `variants.find_variants` ANDs four clauses, and the first two read
  `align_dna.query_aln[0:3]`: the first three columns of the FULL transcript alignment,
  which for any mRNA with a 5′UTR is UTR sequence. When that UTR matches, the chain
  short-circuits and a genuine start loss is never flagged — and `start_lost` gates ORF
  rescue. `cds_span` (GH #46) locates the real start codon.

  A byte-neutral `LIFTON_START_LOST_DIAG` probe evaluated the same four clauses at the
  CDS start across the eight-cell ladder
  (`benchmarks/compare/start_lost_headroom.{py,json,md}`). **42,822 transcripts scored:**

  | | count |
  |---|---|
  | start losses the shipped test MISSES | 392 |
  | of those, that would NEWLY enter ORF rescue | **42 (0.098 %)** |
  | FALSE POSITIVES (shipped says yes, scoped says no) | **219** |

  The mechanism prediction held — the miss appears only on close pairs (drosophila 32,
  human→mouse 10; **0 on all five distant/very-distant cells**), which is exactly where a
  matching 5′UTR is likely. But the effect is negligible, because a lost start codon
  almost always comes with a frameshift or stop anomaly, so those transcripts are already
  re-searching. And the probe caught what an upside-only count could not: **the 219 false
  positives outnumber the 42 genuine gains 5:1**, so correcting the condition is a SWAP,
  not an addition — and ORF-rescue swaps have gone net-negative before (Iteration 9: 4
  improved, 59 regressed, reverted).

  **Decision: NO-GO** on changing the condition. This is a *labelling* defect — the
  `mutation=start_lost` annotation in the output and `score.txt` is wrong for ~1.4 % of
  transcripts in both directions — with essentially no effect on the emitted models. The
  probe was reverted (Iteration-9/13/15 precedent); the harness and results are retained,
  and re-running needs the `LIFTON_START_LOST_DIAG` instrumentation restored in
  `variants.find_variants` plus the `diag_id` pass-through from
  `Lifton_TRANS.orf_search_protein`.
- ~~**`Lifton_TRANS.add_cds` attaches one CDS to every overlapping exon**~~ — **FIXED**:
  each overlapping exon now gets its own copy and the malformed model is reported.
- **Case-1 merged exon inherits the DOWNSTREAM exon's ID** (masked today because
  `normalize_containment` renumbers on collision — see `test_cds_ids.py`).
- **`get_id_fraction` guard drift**: `get_partial_id_fraction` tests `total_length == 0`
  where its sibling uses `<= 0`. Unreachable today; left alone deliberately because
  changing it would change output in the negative case.

### Output quality
- ~~**Merged CDS lines lose their descriptive attributes.**~~ — **FIXED** (2026-07-25).
  See the follow-up section above; the full-output re-measurement was far larger than
  this subset figure (151,277 rows affected, not 35,145), and the fix is proven on the
  eight-cell divergence ladder.

### Performance / memory
- **`scan_annotation` runs full ID bookkeeping + NCBI per-line validation on the derived
  Liftoff/miniprot intermediates**, not just the reference — two extra full-file passes
  and multi-hundred-MB ID dicts per run. Add a `fingerprint_only` scan level.
- **`AnnotationScanResult` is retained for the whole run**; only `cds_namespace_ids` and
  `copy_suffix_ids` (reference only) are needed after the allocator is seeded.
- ~~**`trunc_ref_proteins` builds a whole dict only to take `len()`**~~ — **FIXED**
  (`lifton_utils.count_truncated_proteins`).
- **The gffbase adapter re-hashes the reference up to 4×** because `_expected_manifest`
  omits `source=` (the gffutils path passes it correctly).
- ~~**`align.py` computes the coding sequence twice per transcript**~~ — **FIXED**
  (`get_coding_seq(..., include_sequence=False)`).
- **`orf_search_protein` re-reads the reference protein/transcript FASTA records on every
  call** — up to 3× per merge-firing transcript. Still deferred: it did **not** appear
  anywhere in the top 45 of either Step-7 profile, so the predicted cost was wrong. Worth
  revisiting only if a profile puts it there.
- ~~**`__find_orfs` accumulates a write-only `orf_seq` string** and its outer loop is
  worst-case quadratic~~ — **FIXED**.
- **Step 5 (2 gffutils DB builds, ~15 % of a cached-aligner run) is NOT a threading
  target** — that is the documented Iteration-11 NO-GO (GIL-bound; measured 8–12 % slower).

### Measurement
- **`--measure_time` uses `time.process_time()`** for all steps, so it cannot see the
  miniprot subprocess and double-counts parallel threads. The correct per-phase wall clock
  already exists in `run_manifest.json`; `time.txt` should be sourced from it. Still
  deferred — but `LIFTON_PROFILE_STEP7` (added 2026-07-25) now gives a trustworthy
  *per-function* view of the hot phase, which is what the optimisation work actually
  needed; `--measure_time` remains a reporting wart rather than a blocker.

### Documentation
- ~~**`CLAUDE.md` is stale**~~ — **FIXED** in the first batch (the nine undocumented
  modules and the true test/line counts are now recorded there).

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
