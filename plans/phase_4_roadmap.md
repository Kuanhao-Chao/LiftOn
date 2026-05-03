# LiftOn — Phase 4: Strategic Roadmap for Robustness & Modernization

> **Reference standard:** [NCBI GFF3 Format Guidelines](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/file-formats/annotation-files/about-ncbi-gff3/) (fetched 2026-05-03; key invariants extracted into `lifton/io/ncbi_gff3_spec.py` per Step 2.1).
>
> **Pre-condition:** Phase 3 test suite (97 passing, 1 xfail) is the regression net. Every step below ends with "all 97 tests still pass" — any deliberate test edit is itself a reviewed change.

This roadmap is sequential. Steps 1–4 are foundational fixes that must land before the architectural work in Steps 5–8. Step 9 (parallelism) is intentionally last because it depends on the picklability won by the decomposition in Step 7.

---

## NCBI GFF3 invariants this roadmap enforces

Distilled from the NCBI spec page; each invariant is wired into the validator (Step 2) and into the per-module validation hooks (Steps 5–8).

**Column-level (1-9):**
| # | Column | Invariant |
|---|---|---|
| 1 | `seqid` | Must be present; NCBI recommends accession.version but our validator only requires non-empty + present in target FASTA index |
| 2 | `source` | Free text; LiftOn always rewrites to `LiftOn` |
| 3 | `type` | SOFA term; we accept the SO/INSDC superset gffutils tolerates and warn on unknown types |
| 4-5 | `start`, `end` | 1-based inclusive integers; **`start <= end`** is mandatory (this is currently violated silently — see `gff_malformed_coords` fixture) |
| 6 | `score` | `.` or float |
| 7 | `strand` | `+`, `-`, `.`, or `?` |
| 8 | `phase` | `0`, `1`, `2` for CDS; `.` otherwise. NCBI explicitly notes phase may be wrong for pseudogenes/internal frameshifts — validator must therefore *warn*, not reject |
| 9 | `attributes` | Semicolon-delimited `tag=value`; reserved characters (`;` `=` `&` `,` `\t` `\n` `\r` `%`) MUST be percent-encoded; tags are case-sensitive; uppercase initial = official, lowercase = unofficial |

**Hierarchy (Annotation Data Model):**
- Central dogma `gene → mRNA → exon, CDS` for protein-coding eukaryotes; `gene → {ncRNA, rRNA, tRNA, …} → exon` for non-coding.
- A single feature may legally span **multiple rows with the same ID** (multi-line CDS, alignments). LiftOn currently treats each row as a distinct entry — must be tolerated on input, preserved on output.
- Permitted exceptions: prokaryotes/organelles can have `gene → CDS` with no mRNA; gene segments use `gene → C/V/D/J_gene_segment → CDS`; some pseudogenes have `CDS` with no `mRNA` parent; some tRNAs have no `gene` parent.
- `Parent` value MUST be the `ID` of an entry that appears in the same file (or in `##sequence-region` scope).

**Official attributes (capitalised):** `ID`, `Parent`, `Name`, `Alias`, `Target`, `Gap`, `Derives_from`, `Note`, `Dbxref`, `Ontology_term`, `Is_circular`. **Multi-value attributes** (`Parent`, `Alias`, `Note`, `Dbxref`, `Ontology_term`) are comma-separated lists; the parser must split them, the writer must join them back.

**Directives we must preserve when writing output:** `##gff-version 3` (mandatory first line), `##sequence-region`, NCBI-specific `#!gff-spec-version`, `#!processor`, `#!genome-build`, `#!genome-build-accession`, `#!annotation-date`, `#!annotation-source`. LiftOn currently emits none of these — Step 8.3 fixes this.

---

## Step 1 — Fix the three pinned legacy bugs

**Goal:** Eliminate the bugs documented in `plans/phase_3_test_plan.md §5` so downstream refactors stand on correct foundations.

### 1.1 `Lifton_GENE.entry.id` is one character (`lifton_class.py:71-75`)
**Current (wrong):**
```python
self.entry.attributes["ID"] = self.ref_gene_id + "_" + str(self.copy_num) if self.copy_num > 0 else self.ref_gene_id
self.entry.id = self.entry.attributes["ID"][0]   # takes first char of a string
```
**Fix:**
```python
gene_id = (f"{self.ref_gene_id}_{self.copy_num}"
           if self.copy_num > 0 else self.ref_gene_id)
self.entry.attributes["ID"] = [gene_id]          # list, per gffutils contract
self.entry.id = gene_id
```
**Files:** `lifton/lifton_class.py:57-91` (`Lifton_GENE.__init__`).
**NCBI compliance:** GFF3 §"ID" — IDs are strings, multi-value not permitted. gffutils stores them as `[str]`; this fix restores invariant.
**Test impact:** flip these assertions in `tests/test_lifton_class.py`:
- `test_extra_copy_number_appends_suffix`: `gene.entry.id == "g"` → `"gene1_2"`
- `test_constructor_seeds_tree_dict`: `iv.data == "g"` → `"gene1"`
- `test_transcript_id_assigned`: `Parent == ["g"]` → `["gene1"]`

### 1.2 Transcript `Parent` cascades the bug (`lifton_class.py:236`)
Already fixed transitively by 1.1. Add a positive assertion to `tests/test_lifton_class.py` that `Parent` is a non-empty list of length 1 with full gene id, then drop the legacy-bug comment.

### 1.3 `check_ovps_ratio` passes a tuple to `IntervalTree.overlap` (`lifton_utils.py:521`)
**Fix:**
```python
ovps = tree_dict[mtrans.seqid].overlap(mtrans_interval[0], mtrans_interval[1])
```
**Files:** `lifton/lifton_utils.py:505-530`.
**Test impact:** remove `xfail` mark on `test_overlap_ratio_triggers_true`; the strict xfail will fail loudly if missed.

### 1.4 Verification
```bash
pytest tests/ -v          # 97 still pass; xfail removed
```

---

## Step 2 — Build the GFF3 sanity-check / validator layer

**Goal:** Detect malformed input *before* any heavy work begins. Reject hard errors; warn-and-continue on soft issues per NCBI's own pragmatism (e.g. phase may legitimately be wrong for pseudogenes).

### 2.1 New module: `lifton/io/ncbi_gff3_spec.py`
Single source of truth for the spec invariants used at runtime. Constants only — no logic.
```python
RESERVED_CHARS = {";", "=", "&", ",", "\t", "\n", "\r", "%"}
OFFICIAL_ATTRS = frozenset({"ID", "Parent", "Name", "Alias", "Target",
                            "Gap", "Derives_from", "Note", "Dbxref",
                            "Ontology_term", "Is_circular"})
MULTI_VALUE_ATTRS = frozenset({"Parent", "Alias", "Note", "Dbxref",
                               "Ontology_term"})
VALID_STRANDS = frozenset({"+", "-", ".", "?"})
VALID_PHASES = frozenset({"0", "1", "2", "."})
DIRECTIVE_PREFIX = "##"
NCBI_DIRECTIVE_PREFIX = "#!"
```

### 2.2 New module: `lifton/io/gff3_validator.py`
Streaming, single-pass validator. Emits structured findings; never mutates input.

```python
@dataclass(frozen=True)
class ValidationFinding:
    severity: Literal["error", "warning"]
    line_no: int
    rule: str          # e.g. "start_gt_end", "missing_parent"
    message: str

class GFF3Validator:
    def __init__(self, *, target_seqids: set[str] | None = None,
                 strict: bool = False): ...
    def validate(self, path: Path) -> list[ValidationFinding]: ...
```

**Hard errors (raise after validate completes; default strict=False = log only):**
- Missing `##gff-version 3` directive on line 1.
- Column count != 9 on a non-comment line.
- `start` or `end` not parseable as int, or `start > end` (NCBI cols 4-5).
- `strand` not in `{+ - . ?}`.
- `phase` not in `{0 1 2 .}` for CDS rows.
- `Parent=ID` references an `ID` that never appears in the file (closed after EOF).
- Duplicate `ID` on rows whose `(seqid, type)` differ — true ID collision (multi-row same-feature is allowed when `(seqid, type)` matches per NCBI "ID" §).
- Reserved characters in attribute values that are not percent-encoded (NCBI §Attribute Specifications).

**Warnings (log; don't reject):**
- `seqid` not present in target FASTA index (when `target_seqids` provided).
- Unknown SOFA `type` (warn, not block — SO evolves).
- `phase` = `.` on a CDS row that is not partial (NCBI explicitly says phase may be wrong for pseudogenes).
- Attribute tags starting lowercase that match an official tag spelled wrong (e.g. `parent` vs `Parent`).
- mRNA without an `exon` child (allowed for prokaryote/organelle by NCBI "NOTE 2" but worth surfacing).
- Multi-row feature where rows are not contiguous (NCBI permits but parsers vary).

### 2.3 Wire into the pipeline
**File:** `lifton/lifton.py:208-238` (Step 0–1 of `run_all_lifton_steps`).
After loading FASTAs, before building the gffutils DB, call:
```python
findings = GFF3Validator(
    target_seqids=set(tgt_fai.keys()) | set(ref_fai.keys()),
    strict=args.strict_gff,
).validate(args.reference_annotation)
for f in findings:
    logger.log(f"[GFF3:{f.severity}] line {f.line_no}: {f.rule} — {f.message}",
               debug=True)
if args.strict_gff and any(f.severity == "error" for f in findings):
    sys.exit(2)
```
**New CLI flag** `--strict-gff` (default off for backward compat) added to `lifton/lifton.py:85-133` (`args_optional`).

### 2.4 Verification
- New tests under `tests/test_gff3_validator.py`:
  - Each NCBI invariant has at least one positive (clean fixture) and one negative (corrupt fixture) case.
  - Reuse `gff_malformed_coords` fixture; add `gff_dangling_parent`, `gff_bad_phase`, `gff_unencoded_semicolon`, `gff_missing_directive`, `gff_duplicate_id_collision`, `gff_strand_q`.
  - Assert `validate()` returns the expected findings; assert `--strict-gff` exits non-zero.
- All 97 existing tests still pass (validator is opt-in; the integration test runs without it).

---

## Step 3 — Environment & packaging modernization

**Goal:** Single source of truth for Python, deps, and metadata. Unblocks `uv` / `pdm` / modern resolvers and removes the four-way Python version drift documented in `plans/phase_1_audit.md §3`.

### 3.1 Create `pyproject.toml` (PEP 621 + setuptools backend)
```toml
[build-system]
requires = ["setuptools>=68", "wheel"]
build-backend = "setuptools.build_meta"

[project]
name = "lifton"
version = "1.0.5"
description = "Combining DNA and protein alignments to improve genome annotation"
readme = "README.md"
requires-python = ">=3.10"
license = {text = "GPL-3.0-or-later"}
authors = [{name = "Kuan-Hao Chao", email = "kh.chao@cs.jhu.edu"}]
dependencies = [
    "biopython>=1.81,<2",
    "cigar>=0.1.3,<1",
    "gffutils>=0.13,<0.14",
    "intervaltree>=3.1.0,<4",
    "networkx>=3.3,<4",
    "numpy>=1.22,<3",
    "parasail>=1.3.4,<2",
    "pyfaidx>=0.8,<1",
    "pysam>=0.22,<1",
]

[project.optional-dependencies]
test = ["pytest>=7", "coverage>=7"]

[project.scripts]
lifton = "lifton.lifton:main"

[tool.setuptools.packages.find]
include = ["lifton*"]

[tool.pytest.ini_options]
testpaths = ["tests"]
```

### 3.2 Drop dead deps
- Remove `interlap` from declared deps (unused — `grep -rn "interlap" lifton/` returns 0 hits).
- Move `pytest` from runtime to `[project.optional-dependencies] test`.
- Verify `ujson` import via `grep`; remove if unused.

### 3.3 Unify Python everywhere on 3.10+
| File | Action |
|---|---|
| `setup.py` | **Delete** (replaced by `pyproject.toml`) |
| `lifton.yml` | Pin `python=3.11` (3.10 minimum, but conda env is reproducible at 3.11) |
| `Dockerfile:2` | `FROM python:3.11-slim` (was 3.8-slim) |
| `.github/workflows/tests.yml` | Matrix `python-version: ["3.10", "3.11", "3.12"]` |

### 3.4 Add the missing `minimap2` startup check
**Files:** `lifton/lifton.py:417-435` (`main()`), `lifton/run_liftoff.py` (new helper).
```python
def check_minimap2_installation():
    if shutil.which("minimap2") is None:
        sys.exit("ERROR: minimap2 not found on PATH. "
                 "Install via `conda install -c bioconda minimap2`.")
```
Wire into `main()` immediately after `check_miniprot_installation()`.

### 3.5 Verification
```bash
pip install -e .                 # works on 3.10, 3.11, 3.12
pytest tests/ -v                 # 97 still pass
docker build .                   # uses python:3.11-slim
conda env create -f lifton.yml   # 3.11 env builds
```
**NCBI compliance:** N/A for this step (purely toolchain).

---

## Step 4 — Lazy reference-sequence loading (`RefSeqProvider`)

**Goal:** Replace the eager `extract_features` dict materialisation (`extract_sequence.py:25-55`) with on-demand lookup. Per `plans/phase_2_bottlenecks.md §1a`, this drops peak RSS by hundreds of MB on mammalian genomes.

### 4.1 New module: `lifton/seq/ref_provider.py`
```python
class RefSeqProvider:
    """Lazy, LRU-bounded resolver for reference transcript / protein
    sequences. Backed by either a gffutils FeatureDB (compute on demand)
    or a pre-built FASTA (pyfaidx). Same API in both modes."""

    def __init__(self, *, ref_db=None, ref_fai=None,
                 trans_fasta_path: str | None = None,
                 protein_fasta_path: str | None = None,
                 cache_size: int = 4096): ...

    @lru_cache(maxsize=None)
    def trans(self, ref_id: str) -> str | None: ...
    @lru_cache(maxsize=None)
    def protein(self, ref_id: str) -> str | None: ...
    def keys(self) -> Iterable[str]: ...    # for backward compat with .keys() callers
    def __contains__(self, key: str) -> bool: ...
```

### 4.2 Drop-in replacement for the dict consumers
The only call sites that look at `ref_trans` / `ref_proteins` are:
- `lifton/lifton.py:256, 263, 264, 265` (construction + `len()` logging)
- `lifton/lifton_utils.py:24` (`get_truncated_protein`) — iterates `.keys()`
- `lifton/lifton_class.py:140-144` (`orf_search_protein`) — does `id in ref_proteins.keys()` then `ref_proteins[id]`
- `lifton/run_liftoff.py`, `run_miniprot.py` — pass them through

The `RefSeqProvider` exposes `.keys()` and `__contains__` so all sites continue to work without further edits. The dict-style `provider[id]` is **not** supported — call sites must change `ref_proteins[id]` → `ref_proteins.protein(id)` and `ref_trans[id]` → `ref_trans.trans(id)`.

### 4.3 Secondary fixes in the same pass
- `extract_sequence.py:73-85` — replace `+=` string growth with `''.join(parts)` (O(n²) → O(n)).
- `extract_sequence.py:78` — precompute `chrom_set = set(fasta.keys())`; don't call `fasta.keys()` per transcript.
- `annotation.py:80-88` — `len(CDS_list) > 0` → `any(c.featuretype=='CDS' for c in db.children(...))` to avoid materialising a list.

### 4.4 Verification
- `tests/test_ref_provider.py` (new): cold-call → trans/protein returned; warm-call → LRU hit (mocked clock); `keys()` exhaustive iteration matches the gffutils pass; `__contains__` semantics.
- `tests/test_extract_sequence.py::TestExtractFeatures` updated to assert `provider.protein("tx1")` instead of dict lookup.
- All 97 existing tests still pass after wiring `RefSeqProvider` into `run_all_lifton_steps`.
- New benchmark in `tests/perf/` (skipped by default, `pytest -m perf`): `mprof run` shows >50% RSS reduction on a chr1 GFF.

**NCBI compliance:** N/A (memory layer only).

---

## Step 5 — Extract `lifton/seq/assembly.py` (sequence assembly)

**Goal:** First god-module split. Lift the pure functions (`get_coding_seq`, `get_coding_trans_seq`, `translate_coding_seq`) out of `Lifton_TRANS` into stateless free functions.

### 5.1 New file: `lifton/seq/assembly.py`
```python
def coding_seq(exons: list[Lifton_EXON], fai) -> tuple[str, list, list[int]]:
    """Concatenate CDS sequence across exons. Replaces
    Lifton_TRANS.get_coding_seq (lifton_class.py:453-468)."""

def coding_and_trans_seq(exons: list[Lifton_EXON], fai,
                         strand: str) -> tuple[str, str]:
    """Replaces Lifton_TRANS.get_coding_trans_seq (lifton_class.py:470-494).
    Strand becomes an explicit argument."""

def translate(coding_seq: str) -> str | None:
    """Pure translation. Replaces Lifton_TRANS.translate_coding_seq."""
```

### 5.2 Validation hook
Each function asserts:
- `start <= end` for every exon/CDS (NCBI cols 4-5 invariant).
- `strand in {"+", "-"}` (no `.` or `?` for translatable features).
- Concatenated CDS length is a multiple of 3 *or* the function pads with `N` (current behaviour preserved; NCBI permits partial features which `partial=true` declares).
Validation failures emit a `ValidationFinding` via the same channel as Step 2; do not raise (refactor must be byte-identical).

### 5.3 Backward-compat shim on `Lifton_TRANS`
Methods on the class become one-liners delegating to the new functions, preserving the existing API for one release cycle:
```python
def get_coding_seq(self, fai):
    return seq.assembly.coding_seq(self.exons, fai)
```

### 5.4 Verification
- `tests/test_seq_assembly.py` (new): identity tests against the existing `test_lifton_class.py::TestSequenceAssembly` cases.
- All 97 existing tests still pass (the shim guarantees this).

---

## Step 6 — Extract `lifton/align/transcript_align.py` (alignment dispatch)

**Goal:** Move the stateless alignment wrappers out of the god class so they become safe to call from worker processes (Step 9).

### 6.1 New file: `lifton/align/transcript_align.py`
```python
def align_protein(protein_seq: str, ref_protein_seq: str
                  ) -> tuple[Lifton_Alignment | None, list[str] | None]:
    """Replaces Lifton_TRANS.align_coding_seq (lifton_class.py:502-514).
    Returns (alignment, peptides) — both may be None on empty input."""

def align_transcript(trans_seq: str, ref_trans_seq: str
                     ) -> Lifton_Alignment | None:
    """Replaces Lifton_TRANS.align_trans_seq (lifton_class.py:516-524)."""
```

Both functions return their result instead of mutating a `Lifton_Status` argument; the orchestrator (which still lives on `Lifton_TRANS` for now) merges identity values into the status. This decouples the worker code path from shared state.

### 6.2 Backward-compat shim on `Lifton_TRANS`
Mirrors Step 5.3.

### 6.3 Verification
- `tests/test_transcript_align.py` (new) + identity assertions against `test_lifton_class.py::TestSequenceAssembly`.
- All 97 tests still pass.

---

## Step 7 — Extract `lifton/reconcile/cds_exon.py` (the 5-case god method)

**Goal:** The single highest-value decomposition. `update_cds_list` (`lifton_class.py:266-451`, 185 LOC, 5 branches) becomes five named pure functions in a sibling module. This is the method the Phase 4 parallelism work depends on, because the current implementation mutates `self.exons` in place.

### 7.1 New file: `lifton/reconcile/cds_exon.py`
```python
def reconcile_cds_with_exons(
    exons: list[Lifton_EXON],
    cds_list: list[Lifton_CDS],
    strand: str,
) -> list[Lifton_EXON]:
    """Top-level dispatch. Replaces Lifton_TRANS.update_cds_list."""

def _reconcile_single_cds(exons, only_cds): ...        # Case 1 (lines 275-307)
def _reconcile_single_exon(exons, cds_list): ...       # Case 2 (lines 311-328)
def _reconcile_init_head(exons, cds_list): ...         # Case 3 step 1 (337-381)
def _reconcile_inner(exons, cds_list, ...): ...        # Case 3 step 2 (387-395)
def _reconcile_tail(exons, cds_list, ...): ...         # Case 3 steps 3-5 (399-449)
```

### 7.2 Validation hook
Per-function assertions, all wrapped to emit `ValidationFinding` not raise:
- Each input exon has `start <= end` and `strand` matches the transcript strand (NCBI cols 4-5, 7).
- For strand `-`, `cds_list` is reversed in a *copy*, never in place (kills the current side-effect at `lifton_class.py:271`).
- Output exons are sorted by `(start, end)` and have no zero-length intervals.
- Every output CDS is contained in its parent exon (`exon.start <= cds.start <= cds.end <= exon.end`).

### 7.3 Backward-compat shim on `Lifton_TRANS`
```python
def update_cds_list(self, cds_list):
    self.exons = reconcile.cds_exon.reconcile_cds_with_exons(
        self.exons, cds_list, self.entry.strand,
    )
    self.update_boundaries()
```

### 7.4 Verification
- `tests/test_reconcile_cds_exon.py` (new): unit tests for each of the five cases. The two cases already covered by `test_lifton_class.py::TestUpdateCdsListSingleCds` and `TestUpdateCdsListSingleExonMultipleCds` are mirrored, plus three new fixtures for Cases 3-5 (multi-CDS × multi-exon overlap permutations).
- All 97 existing tests still pass via the shim.
- Coverage of `reconcile/cds_exon.py` ≥ 90 % (was ~30 % effective coverage of the inlined god method).

---

## Step 8 — Extract `lifton/orf/rescue.py` and `lifton/io/gff_writer.py`

**Goal:** Complete the god-module decomposition by moving the ORF rescue logic and the GFF3 serialisation chain into single-purpose modules. Step 8.3 also fixes the missing NCBI directives in our output.

### 8.1 New file: `lifton/orf/rescue.py`
```python
@dataclass
class ORFCandidate:
    start: int
    end: int
    identity: float

def find_best_orf(trans_seq: str,
                  ref_protein_seq: str) -> ORFCandidate | None:
    """Replaces Lifton_TRANS.__find_orfs (lifton_class.py:553-599)."""

def patch_cds_boundary(exons: list[Lifton_EXON],
                       orf: ORFCandidate, strand: str) -> None:
    """Replaces __update_cds_boundary + __iterate_exons_update_cds
    (lifton_class.py:600-657)."""

def cds_frame(accum_cds_length: int) -> int:
    """Replaces __get_cds_frame. Returns 0/1/2 per NCBI col 8."""
```

**Validation hook:** `cds_frame` output is asserted to be in `{0, 1, 2}` (NCBI col 8). `patch_cds_boundary` validates that resulting CDS coordinates remain inside their parent exon.

### 8.2 New file: `lifton/io/gff_writer.py`
```python
def write_gene(fw, gene: Lifton_GENE,
               transcripts_stats_dict: dict) -> None: ...
def write_transcript(fw, trans: Lifton_TRANS) -> None: ...
def write_exon(fw, exon: Lifton_EXON) -> None: ...
def write_cds(fw, cds: Lifton_CDS) -> None: ...
```

Replaces the four `write_entry` methods (`lifton_class.py:158, 662, 718, 739`).

### 8.3 NCBI-compliant directive header (new)
**File:** `lifton/io/gff_writer.py` exposes:
```python
def write_header(fw, *, target_assembly: str | None = None,
                 source_annotation: str | None = None) -> None:
    fw.write("##gff-version 3\n")
    if target_assembly:
        fw.write(f"#!genome-build {target_assembly}\n")
    fw.write("#!processor LiftOn v1.0.5\n")
    fw.write(f"#!annotation-date {datetime.date.today().isoformat()}\n")
    if source_annotation:
        fw.write(f"#!annotation-source {source_annotation}\n")
```
Wired into `lifton/lifton.py:315` immediately after opening `fw`. This brings LiftOn output into compliance with the NCBI Directives §.

### 8.4 NCBI-compliant attribute writer
The current code relies on gffutils' `__str__()` for serialisation. gffutils does percent-encode reserved characters per the spec. We add a defensive post-write check (debug mode only) that re-parses each emitted line through `GFF3Validator.validate_line()` and surfaces any drift. This catches future regressions where someone bypasses gffutils and writes raw strings.

### 8.5 Backward-compat shims
`Lifton_GENE.write_entry` / `Lifton_TRANS.write_entry` / `Lifton_EXON.write_entry` / `Lifton_CDS.write_entry` become one-line delegators to the `gff_writer` functions.

### 8.6 Verification
- `tests/test_orf_rescue.py` (new): direct unit tests for `find_best_orf`, `patch_cds_boundary`, `cds_frame`. Use synthetic transcripts with frameshifts, premature stops, missing-start-codon to drive each ORF-rescue branch.
- `tests/test_gff_writer.py` (new): asserts (a) header includes `##gff-version 3` + `#!processor LiftOn`, (b) every emitted line round-trips through `GFF3Validator` with zero errors, (c) reserved characters in attributes (`;`, `=`, `,`) are percent-encoded on output.
- `tests/test_lifton_class.py::TestWriteEntry` updated: header now precedes the gene line; emitted GFF passes the validator.
- All 97 (now ~140) tests pass.

---

## Step 9 — Gene-level parallelism (`--threads`)

**Goal:** Last because it depends on Steps 5–7 producing picklable, side-effect-free units. Per `plans/phase_2_bottlenecks.md §3`, we add `concurrent.futures.ProcessPoolExecutor` around the Step 7 loop with a deterministic ordered-writer.

### 9.1 New file: `lifton/parallel/step7_runner.py`
```python
def parallel_step7(features, l_feature_db_path, ref_db_path,
                   ref_genome_path, target_genome_path,
                   provider_state, tree_dict, args, *,
                   n_workers: int) -> Iterable[Lifton_GENE]:
    """Producer/ordered-consumer over the Step 7 loop. Yields completed
    Lifton_GENE objects in *submission* order so the GFF writer sees a
    deterministic stream."""

def _init_worker(...): ...     # opens its own pyfaidx + gffutils handles
def _process_one(locus_id, feature_type, ...): ...
```

The worker rehydrates each `gffutils.Feature` from `db_connection[locus_id]` (cheap with SQLite cache); only IDs cross the process boundary.

### 9.2 CLI wiring
`lifton/lifton.py:85-133` (`args_optional`) — `-t/--threads` already exists (default 1) but is forwarded only to minimap2. Extend `run_all_lifton_steps` Step 7 to dispatch to `parallel_step7` when `args.threads > 1`. Leave `--threads 1` on the existing serial path so regression diffs are trivial.

### 9.3 What stays serial (unchanged)
- Step 4 (`exec_liftoff` / `exec_miniprot`) — they already manage their own threads.
- Step 5 (gffutils DB build) — single-writer SQLite.
- Step 8 (miniprot-only loop) — depends on Step 7's IntervalTree; revisit in a follow-up.
- Final writer — single file handle.

### 9.4 Verification
- `tests/test_parallel_step7.py` (new): runs the integration fixture with `--threads 1` and `--threads 4`, asserts the output GFFs are byte-identical (this is the determinism contract).
- All existing tests still pass with `--threads 1` (default).
- Hand-verified speedup on a chr22 fixture: target ≥3× wall-clock reduction on a 4-core box.

**NCBI compliance:** Step 9 is purely scheduling — output is byte-identical to the serial path, so all NCBI invariants enforced in Steps 5-8 carry through.

---

## Cross-cutting test budget

| Step | New tests added | Cumulative |
|---|---:|---:|
| 1 | 0 (edits 3 existing) | 97 |
| 2 | ~25 | ~122 |
| 3 | 0 | 122 |
| 4 | ~10 | ~132 |
| 5 | ~6 | ~138 |
| 6 | ~5 | ~143 |
| 7 | ~12 | ~155 |
| 8 | ~15 | ~170 |
| 9 | ~3 | ~173 |

Coverage gate raised from `--fail-under=55` (Phase 3) to `--fail-under=75` after Step 8 lands.

---

## Rollback discipline

Each step is a separate branch / PR. Merge order is the order above; a failed step blocks all later steps, never just the next one. If Step 7's reconciliation rewrite produces any GFF3 diff vs the Phase 3 golden output, the PR is rejected — no "we'll fix it forward" exceptions, because every later step assumes byte-identical Step 7 output.

---

## Verification matrix (run after every step)

| Check | Command | Expected |
|---|---|---|
| Unit + integration suite | `pytest tests/ -v` | all pass |
| Coverage | `coverage run --source=lifton --omit="lifton/liftoff/*" -m pytest tests/ -q && coverage report --fail-under=<step-target>` | meets gate |
| Golden GFF3 | `lifton -g chr22_ref.gff3 chr22_ref.fa chr22_tgt.fa -o /tmp/out.gff3 && diff /tmp/out.gff3 tests/golden/chr22_lifton.gff3` | empty diff (or reviewed change) |
| GFF3 spec compliance | `lifton --strict-gff -g chr22_ref.gff3 chr22_ref.fa chr22_tgt.fa -o /tmp/out.gff3 && python -m lifton.io.gff3_validator /tmp/out.gff3` | zero errors |
| Memory (Step 4+) | `mprof run lifton ... ; mprof peak` | < 60 % of pre-Step-4 peak |
| Wall-clock (Step 9+) | `lifton --measure_time --threads 8 ...` vs `--threads 1` | Step 7 ≥ 3× faster |
