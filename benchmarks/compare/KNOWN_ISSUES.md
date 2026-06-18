# Known issues surfaced by the distant-pair benchmark expansion (2026-06-18)

The full very-distant flagships (`arabidopsis_to_rice` eudicot→monocot, `human_to_zebrafish` mammal→fish) initially showed LiftOn **devel** apparently "crashing + under-recovering ~⅓ of v1.0.8." Investigation resolved this into **three distinct issues**. Net result after the two fixes below: **devel's full-genome output is at parity with v1.0.8** (At→rice devel == stable 0.148 completeness; human→zebrafish 0.045 vs 0.046, within fresh-`-copies` noise) and beats the best single baseline on protein identity on both.

---

## 1. Eval stale-gffutils-DB bug — *the dominant cause of the bogus number* — **FIXED**

`evaluator._build_db(gff_path)` reused a **stale `<gff>_db` sidecar** instead of rebuilding from the current GFF. When a full flagship's `devel.gff3` was re-run (overwritten with the corrected 7,856-mRNA output), the eval still scored the **leftover 2,730-mRNA DB** from the earlier fast-path run (`devel.gff3` mtime 12:47 vs `devel.gff3_db` mtime 03:44, 2,730 mRNA, missing `rna-ArthCp003`). That fabricated "devel recovers 2,540 / completeness 0.053." `force=True` is passed to the LiftOn `Annotation` loader but did **not** rebuild the on-disk DB.

**Proof devel was never deficient:** all 4,615 "missed" transcripts were present in the new `devel.gff3`; devel output = 7,856 mRNA vs stable 7,817, identical featuretype + dna-identity distributions; devel-recovered ⊂ stable-recovered, devel-only = 0.

**Fix:** `evaluator._build_db` now unlinks `<gff>_db` and `<gff>.eval_db` before building, guaranteeing a fresh ingest of the current GFF on every eval (mirrors `tool_runners._clean_input_dbs`). After the fix, re-scoring jumped devel-full from 0.053 → **0.148 (At→rice, exact parity)** and 0.029 → **0.045 (human, ≈parity)**. Harness-only change; no engine impact.

---

## 2. `--stream` / `--inmemory-liftoff` under-recover (and can segfault) at genome scale — **OPEN (engine), harness-mitigated**

Independent of #1: on a full very-distant lift these in-memory fast paths genuinely write fewer mRNA than the default/disk path, and a related combo segfaults:

| devel config (full At→rice) | mRNA written |
|---|---|
| default / parallel `-t8 --locus-pipeline` (disk I/O) | **7,856** (correct) |
| `-t8 --stream --inmemory-liftoff --locus-pipeline` | **2,730** |
| `-t8 --stream --locus-pipeline -L <cached>` | **SEGFAULT** (exit 139) |

The shortfall is in **Step 8 (miniprot-only rescue)** — likely the in-memory miniprot drain / gffbase native path dropping records on a large stream (At→rice miniprot = 42,890 records). The synthetic 391-byte 24-cell matrix can't see this (a byte-identity-at-scale gap). A genome-scale equivalence test (default vs `--stream`/`--inmemory-liftoff` on a real ref, asserting equal feature counts) would close it.

**Mitigation:** `version_compare._build_argv` emits **only `--locus-pipeline`** for `devel_fast` (parallel Step 7 + default disk I/O — the path pinned byte-identical at the matrix's stream=off/inmemory=off/-t4 cell), so full-genome benchmark runs stay fast and correct. Restore `--stream --inmemory-liftoff` there once this is fixed.

**Reproduce (fast, reuses cached `-L`/`-M`):** see `work/arabidopsis_to_rice/_repro_fixed/` (default path → 7,856) vs `_fastpath_fixed/` (fast path → 2,730).

---

## 3. Step-7 `consume()` crash on a broken-`__str__` exception — **FIXED** (engine)

devel's gene-like lift (Iter 12) processes Liftoff `-copies` organellar gene-like features (`rna-DA397_mg{t,r,p}NN-2`) that v1.0.8 (gene-only) never touches. One hits a gffutils `FeatureNotFoundError` (a copy whose `ref_gene_id` isn't in the reference DB; `run_liftoff.py:112`) whose **own `__str__` is broken** (`<exception str() failed>`), so the old `consume()` log f-string raised `TypeError: __str__ returned non-string` **inside the error handler** and crashed the whole run (default path, 743 mRNA). Fixed in `lifton/locus_pipeline.py` (defensive formatting in `consume()` + `error_tb` capture in the workers) and `lifton/parallel.py`. Byte-neutral (off the happy path); 24-cell matrix green + full suite 684 passed. The 27 skipped organellar copies are benign (un-liftable).
