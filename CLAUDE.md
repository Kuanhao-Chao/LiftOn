# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What LiftOn is

LiftOn is a homology-based genome-annotation lift-over CLI. Given a reference genome `R`, its annotation `R_A`, and a target genome `T`, it produces an annotation `T_A` for `T`. It combines DNA-level alignment (a vendored fork of Liftoff, which can shell out to `minimap2` *or* drive it in-process via the `mappy` PyO3 binding) with protein-level alignment (`miniprot` invoked as a subprocess, with a `pyminiprot`-shaped facade ready to drop in a real PyO3 binding when one ships) and merges the two via an in-house "protein-maximization" chaining algorithm + ORF-rescue pass. CLI entry: `lifton.lifton:main` (`setup.py:34` → `lifton/lifton.py:649`, which calls `run_all_lifton_steps` at `lifton/lifton.py:283`). A second console script `gff3-validate = lifton.gff3_validator:_main` ships from the same package.

## Environment & commands

The project requires native deps (`parasail`, `pysam`, `pyfaidx`, `gffutils`, `duckdb`, `pyarrow`) that ship as bioconda/conda-forge wheels. On macOS/ARM, `pip install parasail` will fail to build from source — use conda. The vendored `lifton/gffbase/` ships a pre-built Rust extension (`_native*.so`); a missing extension falls back to the pure-Python parser at `lifton/gffbase/_pyfallback/`.

```bash
# one-time env (mirrors lifton.yml — Python 3.11)
conda create -n lifton-test -y python=3.11
conda activate lifton-test
conda install -y -c bioconda -c conda-forge \
    parasail-python pysam pyfaidx gffutils intervaltree \
    biopython networkx ujson cigar pytest coverage duckdb pyarrow
pip install mappy   # Phase 16 Tier 5: real --native minimap2 path; also unblocks test_native_bindings.py
pip install -e .

# Run the test suite (524 tests collect; fully hermetic re: minimap2/miniprot —
# they're monkey-patched to raise). 3 files (test_property_based, test_streaming_property,
# test_vulnerabilities) ERROR on collection unless `hypothesis` is installed; without `mappy`
# in the env, ~5 test_native_bindings cases fail (they assert mappy is present).
pytest tests/ -v

# Single test / single class / single case
pytest tests/test_lifton_class.py -v
pytest tests/test_lifton_class.py::TestUpdateCdsListSingleCds -v
pytest tests/test_lifton_class.py::TestWriteEntry::test_tmp_gene_skips_gene_line -v

# The 24-cell byte-identity gate (the load-bearing contract for refactors)
pytest tests/test_native_matrix.py::TestFullNativeMatrix::test_all_24_combinations_byte_identical -v

# Coverage (excludes vendored Liftoff)
coverage run --source=lifton --omit="lifton/liftoff/*" -m pytest tests/ -q
coverage report --skip-covered

# Run LiftOn against real data (requires minimap2 + miniprot on PATH)
lifton -g <ref.gff3> <ref.fa> <target.fa> -o lifton.gff3
# Skip the external runners by reusing prior outputs:
lifton -g <ref.gff3> -L <prior_liftoff.gff3> -M <prior_miniprot.gff3> <ref.fa> <target.fa> -o out.gff3
# Built-in timing harness:
lifton ... --measure_time   # writes time.txt
# Validate a standalone GFF3 file:
gff3-validate path/to/out.gff3
```

CI (`.github/workflows/tests.yml`, Python 3.12) installs `minimap2` + `miniprot` and runs `flake8` (E9/F63/F7/F82 hard-fail, then a soft-warn lint pass) followed by the `test/lifton_chr22_example.sh` shell example. **CI does not run pytest** — pytest is the right local gate to run before any change lands. Note the **Python-version skew** between `setup.py` (`python_requires='>=3.6'`), `lifton.yml` (3.11), and CI (3.12); unifying these is on the cleanup list (`plans/phase_2_bottlenecks.md`).

### Optional fast-path flags

Several flags toggle alternate paths. The four I/O / scheduling fast-paths (`--stream`, `--inmemory-liftoff`, `--threads`/`--locus-pipeline`, `--native`) are **byte-identical to the default output** — the 24-cell matrix test pins this; use them to optimise wall-clock or memory, never to change algorithms. Two (`--strict-gff`, `--validate-output`) are validation gates that don't alter output bytes but can change the exit code. Two **change output**: **`--legacy-merge`** restores the pre-promotion *unconditional* Liftoff↔miniprot merge (the default now runs the verified **best-of-outcome** merge), and **`--full-dp-align`** restores the pre-Iteration-3 *exact giant-only* alignment (the default now "bands everything" — anchor-windowed above ~2500 aa / 8000 nt; see *The byte-identity contract*). Two are kept **no-op aliases**: `--optimize` (best-of-outcome is now default) and `--fast-align` (band-everything is now default — use `--full-dp-align` to opt out).

| Flag | Phase | What it changes |
|---|---|---|
| `--stream` | 7 | Pipe `miniprot` stdout into an in-memory gffbase FeatureDB; skip the `miniprot.gff3` disk round-trip + SQLite re-ingest. |
| `--inmemory-liftoff` | 8 | Serialise vendored Liftoff's `lifted_feature_list` to bytes in-process and feed gffbase directly; skip `liftoff.gff3` disk write + SQLite re-ingest. |
| `--threads N` + `--locus-pipeline` | 9 | Dispatch Step 7 per-locus work through a `ThreadPoolExecutor`; emit in submission order via a heap-backed ordered writer (Phase 15d adds bounded spill-to-disk). Both flags required. |
| `--native` | 10/11 | Drive `minimap2` through the `mappy` PyO3 binding in-process; route miniprot through the `pyminiprot`-shaped facade. **As of Phase 17b**, `--native` unlocks Step-7 threading on *any* backend (including SQLite-backed gffutils) — workers no longer touch the DB connection directly; they read from materialised proxy DBs. The old "gffutils serialises under `--native`" behaviour is now opt-in via `LIFTON_PARALLEL_BLOCK_GFFUTILS=1`. `LIFTON_PARALLEL_FORCE=1` still force-enables threading even without `--native`. |
| `--strict-gff` | 5 | Run the NCBI GFF3 input-side validator on the reference annotation; exit non-zero on any spec violation. |
| `--validate-output` (+ `--validate-verbose`) | 13.5C | After writing, re-validate the output GFF3 with the in-tree `gff3_validator`; print a structured report. |
| `--legacy-merge` | 7 | **Changes output bytes.** Restore the pre-promotion *unconditional* Liftoff↔miniprot merge. The default now keeps, per transcript, the better of {merge+ORF, Liftoff+ORF} (best-of-outcome, `run_liftoff.process_liftoff_with_protein`); `--legacy-merge` reverts to applying the chained CDS unconditionally (manuscript reproduction / A-B baseline). `--optimize` is a kept no-op alias of the default. |
| `--full-dp-align` | Align (Iter 3) | **Changes output bytes.** Restore the pre-Iteration-3 *exact giant-only* alignment: full-DP `nw_trace_scan_sat` for every non-giant gene (gate 8000 aa / 25000 nt), giants still memory-bounded-windowed. The default now "bands everything" (anchor-windowed above 2500 aa / 8000 nt, cap 1500; `lifton/align.py:configure_alignment`). Env equivalents: `LIFTON_FULL_DP_ALIGN=1` (giant-only, for subprocess A/Bs); `LIFTON_ALIGN_WINDOW_{AA,NT}=<huge>` (pure full DP incl. giants — manuscript repro, OOM-prone). `--fast-align` is a kept no-op alias of the default. |

## Big-picture architecture

### The pipeline (`lifton/lifton.py:283-647`, `run_all_lifton_steps`)

Eleven steps + a strict-GFF gate, all in one strictly sequential function. Step numbers follow the in-source comments. (Line numbers drift with every phase — `grep -n "# Step"` to re-anchor.)

1. **Step 0 (~286):** Open target + reference FASTAs via `pyfaidx`, set up output dirs.
2. **Strict-GFF gate / input validator (~327-365):** Always runs `lifton.io.gff3_validator.GFF3Validator` on the reference annotation. Under `--strict-gff` it dumps findings per-row and **exits 2** on any spec error; otherwise (Phase 16 Tier 4) it writes findings to `lifton_output/stats/gff3_input_validation.txt` and emits a single summary line — this silenced the ~500K-line Dbxref stderr flood on RefSeq inputs.
3. **Step 1 (~368):** Build the reference `Annotation` (`lifton/annotation.py`) — gffutils SQLite by default, gffbase-DuckDB when `LIFTON_USE_GFFBASE=1`. GTF is auto-converted to GFF3 unless `--no-auto-convert-gtf`.
4. **Step 2 (~376):** Partition features into protein-coding vs non-coding via `lifton_utils.get_parent_features_to_lift` + `get_ref_liffover_features`.
5. **Step 3 (~383):** **Streaming extractor by default** (Phase 15b). `extract_sequence.extract_features_to_fasta` writes `transcripts.fa` + `proteins.fa` directly to `intermediate_files/`, then re-opens via `pyfaidx.Fasta` for lazy mmap-style access. The user's `-T` / `-P` branch takes the same `pyfaidx` path. The legacy `extract_features` (full in-memory dict) is still in the file for reference but not on the default path. **The historical "dominant memory hot-spot" caveat no longer applies on the default path.**
6. **Step 4 (~441):** Run vendored Liftoff in-process and miniprot as a subprocess. Both can be **skipped** by passing pre-built GFFs via `-L` / `-M` (`lifton_utils.exec_liftoff` / `exec_miniprot` short-circuit on `os.path.exists`). Phase 16 wrapped the vendored-Liftoff call in `run_liftoff.py:57-83` with a `sys.setrecursionlimit(max(orig, 10000))` save/restore guard and a `traceback.format_exc()` dump on failure (the old code logged only `str(e)`, hiding which vendored frame recursed).
7. **Step 5 (~450):** Re-load the Liftoff and miniprot GFFs as gffutils DBs; open output writers. The "Creating … annotation database" log lines go through `_describe_annotation_source` (`lifton.py:9`) so an in-memory `bytes` blob (under `--inmemory-liftoff`/`--stream`) renders as `<in-memory bytes, N bytes>` instead of dumping a 30 MB GFF3 to stderr. Phase 15a emits the GFF3 directive prologue (`##gff-version 3` + every other `##` carried from the input) **before** any feature row, via `lifton.io.gff3_writer.format_directives`.
8. **Step 6 (~500):** Build the ref-id ↔ miniprot-id map (`lifton_utils.miniprot_id_mapping`) and seed per-chromosome `IntervalTree`s (`lifton/intervals.py`).
9. **Step 7 (~512):** Walk every Liftoff gene → build a `Lifton_GENE` → align with parasail → optionally chain with miniprot via `protein_maximization.chaining_algorithm` → ORF-rescue → write the entry. **Wall-clock hot-spot.** Dispatched through `lifton.parallel.parallel_step7` (`StepContext` is the immutable bundle; `LocusResult` is the submission-indexed return). With `--locus-pipeline` + `--threads N > 1`, work runs on a `ThreadPoolExecutor` (threads, not processes — gffutils/gffbase connections aren't picklable; parasail releases the GIL). **Phase 17b** materialises each locus on the parent thread into proxy DBs (`_RefDbProxy`/`_LFeatureDbProxy`/`_MFeatureDbProxy` in `locus_pipeline.py`) so worker threads never touch the SQLite connection — this is what made gffutils + `--native` thread-safe. **Phase 17c** added `_ThreadLocalCtxFactory`: when DBs are backed by on-disk files, the parent-thread materialise loop itself fans out across a 4-thread prefetcher pool (each thread re-opens `gffutils.FeatureDB(dbfn)`); falls back to a serial loop for in-memory/blob DBs. Output is emitted in submission order so `--threads N` is byte-identical to `--threads 1`.
10. **Step 8 (~548):** Walk miniprot-only mRNAs and emit any that pass the overlap/length filters.
11. **Step 9 (~566):** `stats.print_report`; close handles.
12. **Step 10 (~583):** When `--validate-output` is set, re-validate the just-written GFF3 with `lifton.gff3_validator` (hierarchy / CDS phase / containment / LiftOn-attr checks).

An orthogonal **evaluation mode** (`-E` / `-EL`, `lifton.py:413`) re-runs the per-locus loop through `run_evaluation.evaluation` to score an existing annotation instead of producing one; it writes `score.txt` and is the path the benchmark harness's eval pass drives.

### Vendored subtrees

Two subtrees of `lifton/` are external code in spirit and should be treated like frozen dependencies during refactors of the rest of `lifton/`:

- **`lifton/liftoff/`** — complete in-tree fork of upstream Liftoff (~3 KLOC, 18 files, no provenance marker). Invoked **as a library** (`run_liftoff.py:67` calls `liftoff_main.run_all_liftoff_steps(liftoff_args, ref_db)`, or `run_all_liftoff_steps_inmemory` at line 62 under `--inmemory-liftoff`). Liftoff itself shells out to `minimap2` (`liftoff/align_features.py`). Phase 11 added `lifton/liftoff/native_align.py`: when `args.native` is set, alignment is routed through `mappy` in-process via a `_PysamShim` that adapts mappy hits to the existing pysam-shaped consumer.
- **`lifton/gffbase/`** — first-party DuckDB-backed successor to gffutils (~3.6 KLOC + Rust extension `_native*.so`; MIT-licensed; same author). Imported under the `lifton.gffbase` namespace. Optional opt-in for the reference Annotation backend; fully drop-in for FeatureDB / Feature / DataIterator / GFFWriter. The `gffbase_adapter.py` shim translates between LiftOn's `Annotation` shape and gffbase's API. Phases 6.1–6.4 are the relevant plans.

### The `Lifton_TRANS` god module

`lifton/lifton_class.py` (896 LOC) is the central object model and the single most important file to understand. The class hierarchy:

```
Lifton_GENE  →  Lifton_TRANS  →  Lifton_EXON  →  Lifton_CDS
(57-217)        (267-820)        (826-873)       (876-895)
                god class —      thin gffutil    thin
                mixes 5          wrapper, ≤1     wrapper
                concerns         cds child
```

`Lifton_TRANS` (lines 267-820) currently mixes five orthogonal concerns:
- **Construction & exon/CDS bookkeeping** (~268-313)
- **Coordinate reconciliation** — the 5-case `update_cds_list` god method (315-)
- **Sequence assembly + alignment dispatch**
- **ORF rescue + CDS boundary patching** — `__find_orfs` (645-), `__update_cds_boundary` (718-)
- **GFF serialisation** — `write_entry` (794-)

The four-level write chain is `Lifton_GENE.write_entry → Lifton_TRANS.write_entry → Lifton_EXON.write_entry → Lifton_CDS.write_entry`, each delegating to gffutils' `__str__()` for the actual GFF3 line. The Phase 9 dispatcher splits parent-side write from worker-side processing: workers return a `LocusResult` and the parent thread is the only one that calls `write_entry`. This separation is what makes parallel output deterministic.

### Coupling shape

There is a star-shaped dependency around `lifton_utils` and `lifton_class`. `lifton_class` imports `lifton_utils`; `lifton_utils` imports `lifton_class`. This is a benign-but-fragile package-level cycle, currently kept working only by deferred attribute access. Newer modules (`parallel`, `locus_pipeline`, `io/gff3_writer`) sit cleanly outside the cycle and lazy-import where needed (e.g., `lifton.py:523` does `from lifton import parallel as _parallel`). Be cautious about reorganising imports — particularly anything that adds an eager top-level import from `lifton_class` to a previously-cycle-free module.

## The byte-identity contract

The defining test of the whole codebase is the **24-cell matrix** at `tests/test_native_matrix.py::TestFullNativeMatrix::test_all_24_combinations_byte_identical`: every combination of `--stream × --inmemory-liftoff × --threads ∈ {1,2,4} × --native` (= 2·2·3·2 = 24 cells) must produce the same 391-byte output GFF3 as the default path. Smaller matrix gates exist for individual subsets:

- `tests/test_pipeline_streaming.py` — `--stream` byte-identity
- `tests/test_liftoff_inmemory.py` — `--inmemory-liftoff` 4-cell matrix
- `tests/test_parallelism_matrix.py` — `--threads` matrix
- `tests/test_native_matrix.py` — full 24-cell

Any change that touches the alignment kernel, the writer, the directive carrier, the ID-formation logic, or the per-locus result shape must keep this gate green. Phase 14's banded-alignment / mappy-seeded extension proposals were explicitly **deferred** because they would mutate alignment output by ~0.1% and break the gate.

**Scope clarification (best-of-outcome merge promotion, Iteration 1).** The 24-cell matrix pins *fast-path equivalence* (the `--stream × --inmemory-liftoff × --threads × --native` axes all reproduce the default), which is **orthogonal to the merge algorithm**. As of the best-of-outcome promotion, the default *merge* output is **no longer byte-frozen against the pre-promotion default**: when miniprot supplies a higher-identity chunk, the default now emits whichever of {merge+ORF-rescue, Liftoff+ORF-rescue} yields the higher emitted protein identity per transcript, instead of applying the chained CDS unconditionally. The 24-cell matrix stays green because the synthetic fixture has a perfect ORF (`liftoff_aln.identity == 1`) so the merge branch never fires there; `--legacy-merge` reproduces the pre-promotion bytes, and `tests/test_native_matrix.py::TestMergePromotion` (on the `merge_firing_workspace` fixture) pins the promoted behaviour: the merge fires (`status=LiftOn_chaining_algorithm`), `--optimize` is a byte-for-byte no-op alias, and the default never scores below `--legacy-merge`. The aggregate accuracy win is proven on real divergent data by `benchmarks/compare/legacy_merge_ab.py` (drosophila +0.0037 mean protein identity, 115 improved / 1 regressed; mouse_to_rat chr18 +0.0052, 91 improved / 0 regressed — both completeness-preserving, by reverting corrupting merges; the larger mouse chr2 run earlier confirmed +0.0067, 294/2). `mouse_to_rat` is pinned to chr18 (`NC_000084.7`, 2343 mRNA) for a ~3.5×-faster mammalian A/B; `benchmarks/compare/build_inputs.py` rebuilds one benchmark's `-L`/`-M` after a `ref_chrom` change. Promoting any further output-changing accuracy win follows the same recipe: keep the 24-cell matrix green (or deliberately re-baseline with a manuscript erratum), add a merge-firing/behaviour test, and prove the aggregate gain in the benchmark loop.

**Scope clarification (giant-gene windowed alignment, Iteration 2).** The two parasail base aligners (`lifton/align.py:parasail_align_protein_base`/`parasail_align_DNA_base`) call full-DP `nw_trace_scan_sat`, which is O(L²) in **time and memory** — a single titin alignment (~106 kb transcript) is tens of GB and OOM-crashes on mammalian genomes (the mouse benchmark peaked ~48 GB / 1264 s, almost all in ~14 giant genes; see memory [[lifton-giant-gene-align-blowup]]). Above a length gate the aligners route to `lifton/windowed_align.py`: unique-k-mer anchors → co-linear LIS chain → per-window full parasail → concatenated `.traceback` + a parasail-convention `.cigar` shim. The Iteration-2 gate was giant-only (8000 aa / 25000 nt); Iteration 3 (below) lowered it. `tests/test_windowed_align.py` pins losslessness, identity==full-DP, cigar compatibility, and gate routing. The giant-only gate (now reached via `--full-dp-align` / `LIFTON_FULL_DP_ALIGN=1`) changed output **only for above-gate giants**, in score attributes only (coding-giant `protein_identity` Δ = 0.0; long non-coding `dna_identity` ±0.001). Proven by `benchmarks/compare/align_window_ab.py` (mouse, giant-only-windowed vs forced full-DP, `-t 1`): **2.2× faster, 16.9× less peak RSS (44.9 GB → 2.6 GB), 0 non-giant transcripts differ.**

**Scope clarification (band-everything alignment is the default, Iteration 3).** Iteration 3 makes the anchor-windowed aligner the **default for all sizes**, not just giants: the gate is lowered to **2500 aa / 8000 nt** with a fine band (`WINDOW_CAP` 1500) and the exact-DP fallback raised to the giant boundary (`lifton/align.py:configure_alignment(band=True)`, applied at import). This **mutates normal-gene output** and is therefore *not* byte-frozen against the pre-Iteration-3 default — but the 24-cell matrix **stays green with no golden edit** because the synthetic fixtures (~33 aa / 600 nt) are far below the 2500 gate, so they still take exact full DP. Windowing is exact (== full DP) wherever anchors exist (the common homologous case) or a divergent region is below the giant boundary (the raised `max_fulldp` exact fallback); only true-giant anchor-less regions take the bounded approximate split. Escape hatches reproduce the older behaviour: `--full-dp-align` (exact giant-only path) and `LIFTON_ALIGN_WINDOW_{AA,NT}=<huge>` (pure full DP incl. giants, manuscript repro); `--fast-align` is a kept no-op alias. The aggregate win is proven on real data by `benchmarks/compare/fast_align_ab.py` (default vs `--full-dp-align`, `-t 1`): **mouse 2.59×, mouse_to_rat 2.14×, drosophila 1.43× wall; 3.5–4.5× less peak RSS**. Accuracy splits by divergence: **same-species mouse + mouse_to_rat = identity-exact (0 of 13,780 transcripts changed)**; cross-species drosophila is **mean-neutral** (16/7946 changed: 3 improved, 3 regressed by ≤1.2% on already-imperfect alignments; mean Δ ≈ 0) — completeness unchanged on all three. Promotion followed the loop's decision gate (speedup ≥5% wall AND mean identity Δ ≥ −1e-3 AND completeness unchanged). See memory [[lifton-giant-gene-align-blowup]] and [[lifton-optimization-flag-gating]].

## Test-suite guarantees

- The integration tests in `tests/test_integration_pipeline.py` use the `hermetic_pipeline` fixture, which monkey-patches `run_liftoff.run_liftoff` and `run_miniprot.run_miniprot` to **raise** if invoked. CI machines do not need `minimap2` or `miniprot` installed for pytest to pass, and any future regression that accidentally falls through to the external runners trips a hard failure.
- Golden-path testing strategy: provide `-L` and `-M` with pre-baked GFFs so `exec_liftoff` / `exec_miniprot` short-circuit on `os.path.exists`. The same trick is the right way to add new integration scenarios.
- `tests/conftest.py` builds 600-bp synthetic chromosomes whose positions 101-199 (`ATG…`) and 301-399 (`…TAA`) form a clean ORF. Reuse those fixtures rather than authoring new FASTA from scratch.
- Opt-in performance / memory profiling lives under `tests/perf/` and is gated by the `perf` pytest marker (registered in the root `conftest.py`); nothing under `tests/perf/` runs by default.
- A small handful of property-based tests (`test_property_based.py`, `test_streaming_property.py`, `test_vulnerabilities.py`) `import hypothesis` at module scope, so they **ERROR on collection** (not skip) when `hypothesis` is absent — `pytest tests/` reports "524 collected, 3 errors". Install `hypothesis` or deselect those three files for a clean run. The rest of the suite has no extra dependency beyond what conda installs.
- `tests/test_native_bindings.py` deliberately asserts `mappy` is installed (the test env ships it as of Phase 16); without `mappy` those cases fail rather than skip.

### Biological validation benchmarks (`benchmarks/`)

Phase 16 added a real-data harness (`benchmarks/run_benchmarks.py`, registry `datasets.json`) that downloads, lifts, and evaluates the five published JHU CCB same-species datasets (human, mouse, bee, arabidopsis, rice). It shells out to a real `lifton` install (not the hermetic test fixtures), parses GNU/BSD `/usr/bin/time -v` into a `summary_*.json`, and is **not** part of `pytest` — its own unit tests live in `tests/test_benchmark_harness.py`. Phase-17 wall-clock runs use `benchmarks/phase17_rerun.sh` (parameterised by `$DATASET`) and `benchmarks/phase17_human_subset.sh`. Cached inputs/outputs live under `benchmarks/data/` and `benchmarks/results/`. Known open issue: the optional `-E` evaluation pass exits 1 with `mapped: 0` (a harness-side score-parsing mismatch, not a lift-step bug).

## Refactor history & what's in flight

The staged refactor that the existing audit anticipated has largely landed. Plans in `plans/` document each phase end-to-end:

- `plans/phase_1_audit.md` — read-only architectural audit (historical baseline)
- `plans/phase_2_bottlenecks.md` — the original prioritised-targets list (memory, god-module split, parallelism, packaging)
- `plans/phase_3_test_plan.md` — the test-suite shape and the original three pinned bugs
- `plans/phase_4_5_algorithmic_trust.md`, `plans/phase_5_bug_elimination.md` — **landed** (commit `3e2db85`): all three of the original "pinned legacy bugs" are now fixed (`Lifton_GENE` ID built as a list `[gene_id]`, transcript Parent inherits the corrected gene id, `check_ovps_ratio` no longer passes a raw tuple). The `xfail` markers were removed.
- `plans/phase_6_*.md` — vendored gffbase landing (audit → migration → integration)
- `plans/phase_7_streaming_execution.md` — `--stream` (Phase 7)
- `plans/phase_8_liftoff_inmemory_execution.md` — `--inmemory-liftoff` (Phase 8)
- `plans/phase_9_locus_fusion_execution.md` — `parallel.py` + `locus_pipeline.py` + `--threads` / `--locus-pipeline` (Phase 9)
- `plans/phase_10_native_bindings_execution.md` — `lifton/native_bindings/` mappy + pyminiprot facade (Phase 10)
- `plans/phase_11_inloop_rewiring_execution.md` — `--native` in-loop rewiring + thread-safety unlock (Phase 11)
- `plans/phase_12_manuscript_audit.md` — manuscript-vs-code parameter audit
- `plans/phase_13_edge_cases.md`, `plans/phase_13_5{A,B,C}_*.md` — algorithmic hardening + GFF3 formatting fixes (Phase 13)
- `plans/phase_14_optimization_proposal.md` — proposal-only HPC redesign; banded alignment / mappy-seeded extension excluded as byte-output-mutating
- `plans/phase_15_optimization_execution.md` — byte-safe Phase 14 subset shipped: directive carrier (15a), streaming RefSeqProvider (15b), bounded miniprot drain (15c), bounded ordered-writer with on-disk spill (15d)
- *Latest commit (`7152cca`):* `phase-16` — biological validation benchmark harness under `benchmarks/`.
- `plans/phase_16_diagnostic_and_fixes.md` — **uncommitted in the working tree.** Fixes that made the bee benchmark exit 0: GNU-time parser (Tier 1), traceback + 10K recursion-limit guard in `run_liftoff.py` (Tier 2/3), gated the GFF3-validator stderr flood (Tier 4), declared `mappy` in `setup.py`/`lifton.yml` (Tier 5), and the `_describe_annotation_source` bytes-blob log helper.
- `plans/phase_17b_materialisation_refactor.md` — **uncommitted.** The materialisation refactor: per-locus workers read from proxy DBs built from a pre-fetched payload, so `lifton/parallel.py:_backend_supports_threads` now returns True for any backend under `--native` (opt out with `LIFTON_PARALLEL_BLOCK_GFFUTILS=1`). Byte-identical by construction; bee wall-clock was flat (Step 4 subprocess, not Step 7, dominates at bee scale).
- `plans/phase_17c_followon.md` — **uncommitted.** Verified the unlock pays off at scale (rice: 2.6× per-feature speedup vs bee's flat 1.0×) and shrank the parent-thread materialise overhead (Variant-3 query dedup + the `_ThreadLocalCtxFactory` 4-thread prefetcher pool). Items 3 (Step-4 attack) and 4 (DAG-fusion algorithmic rewrite) are explicitly deferred — both would mutate byte output.

**Working-tree state:** the CLAUDE.md you're reading and the Phase 16/17 changes above are **uncommitted** (`git status` shows modified `lifton/{lifton,parallel,locus_pipeline,run_liftoff,annotation}.py`, `setup.py`, `lifton.yml`, benchmark/test additions). 524 tests collect and pass under the `lifton_devel` env.

Read the phase plans for any structural change. The core invariant for any future work is the **24-cell byte-identity matrix**: refactors below that line are free; refactors above it require a deliberate test edit and a manuscript erratum.
