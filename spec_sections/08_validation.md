## 8. Validation Subsystems

LiftOn ships **three independent validation modules**, each owning a distinct point in the lifecycle:

| § | Module | File | When it runs | Entry trigger |
|---|---|---|---|---|
| 8.1 | Input-side NCBI GFF3 validator | `lifton/io/gff3_validator.py` + `lifton/io/ncbi_gff3_spec.py` | Pipeline Step 0→1 boundary, **always** (gate behaviour differs by `--strict-gff`) | `lifton.py:327-366` |
| 8.2 | Output-side GFF3 validator | `lifton/gff3_validator.py` | (a) standalone `gff3-validate` console tool; (b) pipeline Step 10 when `--validate-output` set | `setup.py:35`; `lifton.py:585-603` |
| 8.3 | Pre-flight annotation validator | `lifton/annotation_validator.py` | Inside `Annotation.__init__` before any gffutils DB build, **always**; also inside vendored Liftoff's `extract_features` | `annotation.py:139,267`; `liftoff/extract_features.py:83` |

These three modules are **completely separate implementations** that do not share code. Each defines its own `ValidationResult`/`ValidationFinding` dataclass, its own severity vocabulary, and its own report formatter. A reimplementer must build all three independently. The two `print_validation_report` functions (in §8.2 and §8.3) and the two `_wrap` helpers are distinct functions with the same name in different modules.

---

### 8.1 Input-side NCBI GFF3 validator (`--strict-gff`)

A **streaming, single-pass** validator that reads the reference annotation line-by-line and emits structured findings for violations of the NCBI GFF3 specification. It never mutates the file.

#### 8.1.1 Spec constants (`lifton/io/ncbi_gff3_spec.py`)

This file is **constants only — no logic**. A reimplementer must hard-code these exact values:

| Constant | Type | Value |
|---|---|---|
| `RESERVED_CHARS` (`ncbi_gff3_spec.py:12`) | `frozenset[str]` | `{ ";", "=", "&", ",", "\t", "\n", "\r" }` — characters that MUST be percent-encoded inside an attribute value |
| `OFFICIAL_ATTRS` (`ncbi_gff3_spec.py:16`) | `frozenset[str]` | `{ "ID", "Parent", "Name", "Alias", "Target", "Gap", "Derives_from", "Note", "Dbxref", "Ontology_term", "Is_circular" }` |
| `MULTI_VALUE_ATTRS` (`ncbi_gff3_spec.py:22`) | `frozenset[str]` | `{ "Parent", "Alias", "Note", "Dbxref", "Ontology_term" }` (defined but **not referenced** by the validator logic; imported but unused) |
| `VALID_STRANDS` (`ncbi_gff3_spec.py:27`) | `frozenset[str]` | `{ "+", "-", ".", "?" }` — note `"?"` is allowed here (differs from §8.2 which forbids it) |
| `VALID_PHASES` (`ncbi_gff3_spec.py:32`) | `frozenset[str]` | `{ "0", "1", "2", "." }` — string values, `"."` permitted (NCBI tolerates it for pseudogenes / frameshifts) |
| `DIRECTIVE_PREFIX` (`ncbi_gff3_spec.py:35`) | `str` | `"##"` |
| `NCBI_DIRECTIVE_PREFIX` (`ncbi_gff3_spec.py:36`) | `str` | `"#!"` |
| `GFF_VERSION_DIRECTIVE` (`ncbi_gff3_spec.py:39`) | `str` | `"##gff-version 3"` |

#### 8.1.2 Finding data model (`ValidationFinding`, `io/gff3_validator.py:36-45`)

A `@dataclass(frozen=True)`:

| Field | Type | Meaning |
|---|---|---|
| `severity` | `str` | `"error"` or `"warning"` (lowercase) |
| `line_no` | `int` | 1-based line number; `0` means "file-level" |
| `rule` | `str` | short stable identifier (table in §8.1.4) |
| `message` | `str` | human-readable explanation |

`__str__` (`io/gff3_validator.py:43-45`) formats as: `[GFF3:{severity}] {loc}: {rule} — {message}` where `loc` is `f"line {line_no}"` when `line_no` is truthy, else the literal `"file"`. The em-dash `—` is U+2014.

#### 8.1.3 `GFF3Validator` class (`io/gff3_validator.py:63-252`)

Constructor `__init__(self, *, target_seqids: set[str] | None = None, strict: bool = False)` (`io/gff3_validator.py:77`). Keyword-only args. Instance state:

| Field | Type | Purpose |
|---|---|---|
| `self.target_seqids` | `set[str] \| None` | seqids known present in target+reference FASTA; `None` disables the seqid-membership check |
| `self.strict` | `bool` | stored but **not used by the validator itself** — the caller reads it (the class always returns the full finding list regardless) |
| `self._findings` | `list[ValidationFinding]` | accumulator |
| `self._declared_ids` | `set[str]` | every `ID=` value seen |
| `self._referenced_parents` | `list[tuple[int, str]]` | `(line_no, parent_id)` for every Parent reference, for the deferred dangling-parent check |

Public methods:
- `validate(path) -> list[ValidationFinding]` (`:87`) — resets all three accumulators, opens the file with `encoding="utf-8", newline=""`, calls `_validate_line(line_no, raw)` for each 1-based line, then runs the deferred parent-resolution pass, then returns `list(self._findings)`.
- `validate_line(line_no, raw) -> list[ValidationFinding]` (`:104`) — public per-line entry that does **not** reset `_declared_ids`; returns only findings appended by this call (slices `self._findings[before:]`). Intended for post-write sanity checks.
- `has_errors() -> bool` (`:112`) — `any(f.severity == "error" for f in self._findings)`.
- `findings` property (`:115`) — copy of `_findings`.

Module-level convenience: `validate_path(path, *, target_seqids=None, strict=False)` (`:255`) constructs a `GFF3Validator` and calls `.validate(path)`.

Helpers:
- `_PCT_RE = re.compile(r"%[0-9A-Fa-f]{2}")` (`:49`) — matches a valid percent-encoded byte.
- `_has_unencoded_reserved(value) -> bool` (`:52`): substitute every `_PCT_RE` match out of `value`, then return `True` iff any remaining character is in `RESERVED_CHARS`. **Gotcha:** because `=` is in `RESERVED_CHARS`, but the attribute parser splits on the *first* `=` via `partition("=")`, a second `=` inside the value triggers this rule. This is the source of the high-volume `unencoded_reserved_char` findings on RefSeq `Dbxref` values (the `,` separator inside `Dbxref=GeneID:1,Genbank:NM_…` is reserved, and individual `DBTAG:ID` tokens can carry reserved chars).

#### 8.1.4 `_validate_line` algorithm (`io/gff3_validator.py:127-252`), numbered

1. **BOM strip (line 1 only):** if `line_no == 1` and `raw` starts with the UTF-8 BOM `"﻿"` (written literally as `"﻿"` in source), emit `warning / utf8_bom`, then `raw = raw.lstrip("﻿")`.
2. `line = raw.rstrip("\r\n")`. If `line` is empty (falsy), **return immediately** (blank lines are skipped, contribute nothing).
3. **Directive lines:** if `line.startswith("##")` OR `line.startswith("#!")`:
   - if `line_no == 1` and `line` does **not** start with `"##gff-version 3"`, emit `error / missing_gff_version` ("First line must be '##gff-version 3' (NCBI § Directives).").
   - **return** (directives are otherwise unvalidated).
4. **Comment lines:** if `line.startswith("#")` (single hash, not a directive), **return** (tolerated).
5. **First-line-is-a-feature guard:** if `line_no == 1` (and we reached here, i.e. it is a feature row), emit `error / missing_gff_version` ("First line must be '##gff-version 3' before any feature rows…"). Execution continues (does not return).
6. **Column count:** `cols = line.split("\t")`. If `len(cols) != 9`, emit `error / bad_column_count` and **return**. Otherwise unpack `seqid, _src, _type, start_s, end_s, _score, strand, phase, attrs = cols`.
7. **Col 1 (seqid):**
   - if `not seqid` (empty) → `error / empty_seqid`.
   - elif `self.target_seqids is not None and seqid not in self.target_seqids` → `warning / unknown_seqid` (message includes `seqid!r`). **Gotcha:** when `target_seqids is None`, this check is fully disabled (no warning ever).
8. **Cols 4-5 (start/end):** `try: start=int(start_s); end=int(end_s)`. On `ValueError` → `error / bad_coordinate` and **return**. Then:
   - `start < 1` → `error / negative_start` (rule name is `negative_start` despite checking `< 1`).
   - `end < 1` → `error / negative_end`.
   - `start > end` → `error / start_gt_end`.
   - These three are independent (all can fire on one row).
9. **Col 7 (strand):** if `strand not in VALID_STRANDS` (i.e. not one of `+ - . ?`) → `error / bad_strand` (message lists `sorted(VALID_STRANDS)`).
10. **Col 8 (phase):**
    - if `phase not in VALID_PHASES` (not one of `"0" "1" "2" "."`) → `error / bad_phase`.
    - **elif** `_type == "CDS" and phase == "."` → `warning / cds_missing_phase`. **Gotcha:** the `elif` means a CDS with valid numeric phase `"0"/"1"/"2"` is fine; only `"."` triggers the warning, and a totally invalid phase (e.g. `"3"`) takes the `if` branch and never reaches this elif.
11. **Col 9 (attributes):** build `attr_dict: dict[str,str]`. Split `attrs` on `";"`; for each `piece`:
    - `piece = piece.strip()`; skip if empty.
    - if `"=" not in piece` → `error / bad_attribute` ("Attribute … is not 'key=value'.") and `continue`.
    - `key, _, value = piece.partition("=")` (splits on **first** `=` only). Store `attr_dict[key] = value`. **Gotcha:** `key` is not stripped, so a leading space (already removed by the `piece.strip()`) cannot occur, but trailing whitespace before `=` is preserved into `key`.
    - **Miscapitalisation check:** if `key not in OFFICIAL_ATTRS` AND `key[:1].isupper()` AND `key.lower() in {a.lower() for a in OFFICIAL_ATTRS}` → `warning / miscapitalised_attribute`. I.e. an attribute that starts uppercase and case-folds onto an official name but is not exactly official (e.g. `Parent` mis-typed `PArent`). Note `parent` (lowercase initial) fails the `key[:1].isupper()` guard and is **not** flagged.
    - **Reserved-char check:** if `_has_unencoded_reserved(value)` → `error / unencoded_reserved_char`.
12. **ID / Parent tracking (for deferred check):**
    - if `"ID" in attr_dict` → `self._declared_ids.add(attr_dict["ID"])`.
    - if `"Parent" in attr_dict` → split the value on `","`, strip each, and for every non-empty `parent_id` append `(line_no, parent_id)` to `_referenced_parents`.

#### 8.1.5 Deferred parent-resolution pass (`io/gff3_validator.py:96-101`)

After the whole file is consumed, iterate `self._referenced_parents`; for each `(line_no, parent_id)` where `parent_id not in self._declared_ids`, emit `error / dangling_parent` ("Parent={parent_id!r} references an ID that never appears in the file (NCBI § Parent)."). **Gotcha:** this is order-independent (declared IDs from anywhere in the file count), but a Parent that points forward to an ID declared later still resolves correctly because the check runs only after the full pass.

#### 8.1.6 Complete rule catalogue (8.1)

| Rule id | Severity | Trigger |
|---|---|---|
| `utf8_bom` | warning | line 1 starts with UTF-8 BOM |
| `missing_gff_version` | error | line 1 is a `##`/`#!` directive other than `##gff-version 3`, OR line 1 is a feature row |
| `bad_column_count` | error | feature row split on tab ≠ 9 columns |
| `empty_seqid` | error | column 1 empty |
| `unknown_seqid` | warning | `target_seqids` provided and seqid not in it |
| `bad_coordinate` | error | start or end not an integer |
| `negative_start` | error | `start < 1` |
| `negative_end` | error | `end < 1` |
| `start_gt_end` | error | `start > end` |
| `bad_strand` | error | strand ∉ `{+ - . ?}` |
| `bad_phase` | error | phase ∉ `{0 1 2 .}` |
| `cds_missing_phase` | warning | type == `CDS` and phase == `.` |
| `bad_attribute` | error | attribute piece has no `=` |
| `miscapitalised_attribute` | warning | uppercase-initial key case-folds onto an official attr but is not exactly official |
| `unencoded_reserved_char` | error | attribute value (after stripping `%XX`) contains any of `; = & , \t \n \r` |
| `dangling_parent` | error | a `Parent=` value never appears as any `ID=` |

#### 8.1.7 Strict-mode vs non-strict gate (`lifton.py:323-366`)

The pipeline always runs the validator (`lifton.py:329-332`) with `target_seqids = set(tgt_fai.keys()) | set(ref_fai.keys())` (union of target and reference FASTA seqids) and `strict=getattr(args, "strict_gff", False)`, against `args.reference_annotation`. Then `_strict_gff = getattr(args, "strict_gff", False)` selects:

1. **Strict mode (`--strict-gff` set):** every finding is logged per-row to stderr via `logger.log(str(f), debug=True)` (`lifton.py:342-343`).
2. **Non-strict (default) with findings present:** write **all** findings to a side-car file `gff3_input_validation.txt` under the stats directory (`os.path.join(stats_dir, "gff3_input_validation.txt")`, `lifton.py:345`), one `str(finding)` per line. Then count `n_err = #error`, `n_warn = len(findings) - n_err` and emit one summary line via `logger.log_info`: `">> GFF3 input validator: {N} finding(s) ({n_err} error, {n_warn} warning) written to {path}; pass --strict-gff to also dump per-row to stderr."`. **Fallback:** if the side-car write raises `OSError`, emit a `log_warning` and fall back to the per-row stderr dump.
3. **Non-strict with no findings:** nothing printed.
4. **Exit gate:** `if _strict_gff and any(f.severity == "error" for f in findings): sys.exit(2)` (`lifton.py:365-366`). **Gotcha:** exit code is **2** (distinct from the output validator's exit 1). In non-strict mode the pipeline **never aborts** on input findings, no matter how many errors — they are advisory only. Warnings alone never cause exit even in strict mode.

The rationale comment (`lifton.py:333-339`) documents that real RefSeq inputs trigger hundreds of thousands of `unencoded_reserved_char` findings on `Dbxref` values, which is why the default path was changed to a side-car file rather than a 100+ MB stderr dump.

---

### 8.2 Output-side GFF3 validator (`gff3-validate` / `--validate-output`)

A **full in-memory** validator (parses the entire file into records, then runs cross-record checks). It validates a LiftOn-produced GFF3 against the NCBI spec **and** LiftOn conventions. Used two ways:
- standalone console tool `gff3-validate = lifton.gff3_validator:_main` (`setup.py:35`); also `python -m lifton.gff3_validator output.gff3`.
- pipeline Step 10 when `--validate-output` is set (`lifton.py:585-603`).

#### 8.2.1 Module constants (`gff3_validator.py:37-76`)

| Constant | Value |
|---|---|
| `GENE_TYPES` (`:37`) | `{ "gene", "pseudogene", "transposable_element", "LiftOn-gene" }` |
| `TRANSCRIPT_TYPES` (`:41`) | `{ "mRNA", "ncRNA", "transcript", "lncRNA", "lnc_RNA", "nc_RNA", "rRNA", "tRNA", "miRNA", "snoRNA", "snRNA", "scRNA", "primary_transcript", "processed_pseudogene", "three_prime_overlapping_ncrna", "antisense_RNA", "antisense", "guide_RNA", "RNase_MRP_RNA", "RNase_P_RNA", "SRP_RNA", "vault_RNA", "Y_RNA", "telomerase_RNA", "C_gene_segment", "V_gene_segment", "D_gene_segment", "J_gene_segment" }` |
| `EXON_TYPES` (`:51`) | `{ "exon" }` |
| `CDS_TYPES` (`:52`) | `{ "CDS" }` |
| `REGION_TYPES` (`:53`) | `{ "region" }` (defined, not used in checks) |
| `VALID_STRANDS` (`:55`) | `{ "+", "-", "." }` — note **no `"?"`** (differs from §8.1) |
| `VALID_PHASES` (`:56`) | `{ 0, 1, 2 }` — **integers** (differs from §8.1's strings) |
| `OFFICIAL_ATTRS` (`:59`) | same 11 names as §8.1 (defined, not referenced by checks) |
| `LIFTON_TRANS_ATTRS` (`:64`) | `{ "protein_identity", "dna_identity", "extra_copy_number", "annotation_source" }` (defined, not directly referenced) |
| Column index constants (`:68-76`) | `COL_SEQID=0, COL_SOURCE=1, COL_TYPE=2, COL_START=3, COL_END=4, COL_SCORE=5, COL_STRAND=6, COL_PHASE=7, COL_ATTRS=8` |

`class Severity` (`:83`) defines string constants `ERROR="ERROR"`, `WARNING="WARNING"`, `INFO="INFO"` (uppercase, differs from §8.1's lowercase).

#### 8.2.2 Data model

`GFF3Issue` (`@dataclass`, `:93-104`):

| Field | Type | Meaning |
|---|---|---|
| `severity` | `str` | one of `ERROR`/`WARNING`/`INFO` |
| `lineno` | `int` | 1-based; `0` → "global" |
| `feature_id` | `str` | feature ID (may be `""`) |
| `check` | `str` | check name (grouping key in the report) |
| `message` | `str` | detail |

`__str__` (`:101`) → `[{severity}] {loc}{fid} — {message}` where `loc = "line N"` if `lineno>0` else `"global"`, `fid = " [feature_id]"` if non-empty.

`GFF3Record` (`@dataclass`, `:107-128`) — one parsed data line:

| Field | Type | Meaning |
|---|---|---|
| `lineno` | `int` | source line |
| `seqid, source, ftype` | `str` | cols 1-3 |
| `start, end` | `int` | 1-based coords (cols 4-5) |
| `score, strand, phase` | `str` | cols 6-8 (kept as strings) |
| `attrs` | `Dict[str, List[str]]` | parsed col 9, comma-split values |
| `raw` | `str` | the original line |
| `feat_id` (property) | `str` | `attrs.get("ID", [""])[0]` |
| `parent_id` (property) | `str` | `attrs.get("Parent", [""])[0]` (**first** Parent only) |

`ValidationResult` (`@dataclass`, `:131-150`):

| Field | Type | Meaning |
|---|---|---|
| `file_path` | `str` | input path |
| `total_lines`, `data_lines`, `comment_lines` | `int` | counts |
| `issues` | `List[GFF3Issue]` | all issues |
| `stats` | `Dict[str,int]` | feature-type histogram |
| `errors` (property) | list | issues with severity `ERROR` |
| `warnings` (property) | list | issues with severity `WARNING` |
| `is_valid` (property) | `bool` | `len(errors) == 0` (warnings/INFO do **not** invalidate) |

#### 8.2.3 Public API `validate_gff3_file` (`gff3_validator.py:157-277`)

Signature:
```python
def validate_gff3_file(gff3_path, check_hierarchy=True, check_cds_phase=True,
                       check_containment=True, check_lifton_attrs=True,
                       max_issues_per_check=50) -> ValidationResult
```

Algorithm:
1. **File pre-checks:** if `not os.path.exists(path)` → append `ERROR / file_exists` and return early. If `os.path.getsize(path) == 0` → append `ERROR / file_not_empty` and return early.
2. **Parse pass:** call `_parse_gff3(path, max_issues_per_check)` inside `try/except`; any exception → append `ERROR / parse_error` and return. Store `total/data/comment_lines` from `meta`, extend `issues` with parse issues.
3. **Build feature index:** iterate `records` once:
   - track `seen_ids: dict[id→first_lineno]`. If a record's `feat_id` already seen → increment `duplicate_id` counter and (if ≤ `max_issues_per_check`) append `ERROR / duplicate_id` ("Duplicate ID '…' (first seen on line N)"). Else register `seen_ids[fid]=lineno` and `id_to_record[fid]=rec`. **Gotcha:** a duplicate is *not* added to `id_to_record`, so only the first occurrence is the canonical parent for later containment/hierarchy lookups.
   - `parent_to_children[pid].append(rec)` for every record with a parent.
4. Run, gated by flags: `_check_hierarchy`, `_check_containment`, `_check_cds_phase`, `_check_lifton_attrs` (each extends `issues`).
5. `result.stats = _compute_stats(...)` (feature-type histogram).
6. Return `result`.

Every check applies a **per-check-type cap**: each module keeps a local `issue_counts` defaultdict; a finding is appended only while `issue_counts[key] <= max_issues_per_check` (default **50**). The counter increments unconditionally, so the cap suppresses display but not counting.

#### 8.2.4 `_parse_gff3` column-level checks (`:351-531`)

Reads with `open(path, "r", errors="replace")`, enumerating from line 1. Per line:
1. `line = raw.rstrip("\n\r")`.
2. **Comment/blank** (`line.startswith("#")` or `line.strip()==""`): increment `comment_lines`; **special case** — if `lineno==1` and line does not start with `"##gff-version"` → `WARNING / gff3_header` ("First line should be '##gff-version 3' directive"). `continue`.
3. `cols = line.split("\t")`. If `len(cols) != 9` → `ERROR / column_count` and `continue`.
4. increment `data_lines`. Strip cols into `seqid, source, ftype, score, strand, phase, attrs_str`.
5. **seqid:** if empty or `"."` → `ERROR / seqid_empty`.
6. **start/end:** `int(cols[3])`, `int(cols[4])`; on `ValueError` → `ERROR / coord_not_int` and `continue`. Then:
   - `start < 1` → `ERROR / coord_1based`.
   - `end < start` → `ERROR / coord_order`.
7. **score:** if `score != "."` and `float(score)` raises → `WARNING / score_format`.
8. **strand:** if `strand not in {+ - .}` → `ERROR / strand_valid`.
9. **phase:** if `ftype in CDS_TYPES`:
   - `phase == "."` → `ERROR / cds_phase_required`.
   - else `int(phase)`; if not in `{0,1,2}` (or non-int) → `ERROR / cds_phase_value`.
   - else (non-CDS) if `phase != "."` → `WARNING / non_cds_phase`.
10. **attributes:** `_parse_attributes(...)` (below).
11. **required ID:** if `ftype in GENE_TYPES or ftype in TRANSCRIPT_TYPES` and `"ID" not in attrs` → `ERROR / missing_id`.
12. Build `GFF3Record` and append.

`_parse_attributes` (`:534-574`): if `attrs_str` empty or `"."` → empty dict. Split on `";"`; per part: strip; skip empty; if no `"="` → `ERROR / attr_format`; `key, _, value = part.partition("=")`, strip both; if empty key → `WARNING / attr_empty_key` and skip; else `attrs[key] = [v.strip() for v in value.split(",") if v.strip()]` (comma-split, empties dropped). **Note** `_ATTR_RE = re.compile(r'^[a-zA-Z_][a-zA-Z0-9_]*=')` is defined at `:369` but is **unused dead code**.

#### 8.2.5 `_check_hierarchy` (`:581-731`)

Per-record (`:604-685`):
- **Rule 1 — orphan_parent (ERROR):** `pid` non-empty and `pid not in all_ids` → `orphan_parent`.
- **Rule 2 — gene_has_parent (ERROR):** `ftype in GENE_TYPES` and `pid` truthy.
- **Rule 3 — transcript parent (ERROR):** `ftype in TRANSCRIPT_TYPES`: if no `pid` → `transcript_no_parent`; elif `pid in id_to_record` and the parent's `ftype not in GENE_TYPES` → `transcript_parent_type`.
- **Rule 4 — exon parent (ERROR):** `ftype in EXON_TYPES`: no parent → `exon_no_parent`; else if parent type ∉ TRANSCRIPT_TYPES → `exon_parent_type`.
- **Rule 5 — CDS parent (ERROR):** `ftype in CDS_TYPES`: no parent → `cds_no_parent`; else if parent type ∉ TRANSCRIPT_TYPES → `cds_parent_type`.

Post-loop over `id_to_record` (`:688-729`):
- **Rules 6/7 (WARNING):** for each transcript-type record: gather `exon_children`/`cds_children` from `parent_to_children[trans_id]`. If neither exons nor CDS → `transcript_has_exons` warning. If `ftype == "mRNA"` (exact string) and no CDS children → `mrna_has_cds` warning.
- **Rule 8 (WARNING):** for each gene-type record: if it has **no children at all** (`not trans_children and not other_children`, where `other_children` is the full child list) → `gene_has_transcripts` warning. **Gotcha:** the condition is logically `not other_children` alone (since `other_children` is a superset of `trans_children`), so a gene with any child of any type passes.

#### 8.2.6 `_check_containment` (`:738-790`)

Per record with a resolvable parent (`pid in id_to_record`):
- **seqid_consistency (ERROR):** `rec.seqid != parent.seqid`.
- **strand_consistency (ERROR):** both strands ∉ `{".", ""}` and `rec.strand != parent.strand`.
- **coord_containment (ERROR):** `rec.start < parent.start or rec.end > parent.end`.

**Gotcha:** containment is checked only against the **first Parent** (`parent_id` property) and only if that parent's ID was registered (non-duplicate). Records whose parent is missing/duplicate are silently skipped.

#### 8.2.7 `_check_cds_phase` (`:797-851`)

For each `(trans_id, children)` in `parent_to_children`:
1. `cds_list = [c for c in children if c.ftype in CDS_TYPES]`. If `len(cds_list) < 2` → skip (single-CDS transcripts are not validated).
2. Look up `trans_rec = id_to_record.get(trans_id)`; skip if absent. `strand = trans_rec.strand`.
3. **Sort 5'→3':** if `strand == "-"` → `sorted(cds_list, key=lambda r: r.end, reverse=True)`; else `sorted(key=lambda r: r.start)`.
4. `accum_len = 0`. For each `i, cds` in the sorted list:
   - `expected_phase = (3 - accum_len % 3) % 3`.
   - `actual_phase = int(cds.phase)`; on `ValueError`, add `cds.end - cds.start + 1` to `accum_len` and `continue`.
   - **if `i > 0` and `actual_phase != expected_phase`** → `WARNING / cds_phase_consistency` ("CDS #{i+1} of transcript … expected phase {expected_phase} (accum_len=…), got {actual_phase}"). **Gotcha:** the very first CDS (`i==0`) is never phase-checked — its phase is assumed correct and only seeds `accum_len`. The formula is `(3 - (cumulative_length_of_preceding_CDS % 3)) % 3`.
   - `accum_len += cds.end - cds.start + 1`.

#### 8.2.8 `_check_lifton_attrs` (`:858-938`)

Per record:
- **lifton_source (INFO):** if `ftype in GENE_TYPES | TRANSCRIPT_TYPES | EXON_TYPES | CDS_TYPES` and `source not in {"LiftOn","miniprot","Liftoff","."}` → INFO finding.
- **lifton_protein_identity (ERROR):** if `"protein_identity"` present, take `attrs["protein_identity"][0]`; `float(val)`; must satisfy `0.0 <= val <= 1.0` else (or `ValueError`/`IndexError`) → ERROR.
- **lifton_dna_identity (ERROR):** identical logic for `dna_identity`.
- **Transcript identity presence (WARNING):** for transcript-type records, `has_cds_children = any(child.ftype in CDS_TYPES ...)`:
  - if `has_cds_children` and no `protein_identity` → `lifton_attrs_present` warning ("Coding transcript is missing 'protein_identity' attribute").
  - if no `dna_identity` → `lifton_attrs_present` warning ("Transcript is missing 'dna_identity' attribute"). **Gotcha:** the `dna_identity` warning fires for **every** transcript lacking it, coding or not.

#### 8.2.9 Complete check catalogue (8.2)

| Check name | Severity | Source | Trigger |
|---|---|---|---|
| `file_exists` | ERROR | validate_gff3_file | path missing |
| `file_not_empty` | ERROR | validate_gff3_file | zero-byte file |
| `parse_error` | ERROR | validate_gff3_file | parse raised |
| `gff3_header` | WARNING | parse | line 1 not `##gff-version` |
| `column_count` | ERROR | parse | not 9 columns |
| `seqid_empty` | ERROR | parse | seqid empty or `.` |
| `coord_not_int` | ERROR | parse | start/end not int |
| `coord_1based` | ERROR | parse | `start < 1` |
| `coord_order` | ERROR | parse | `end < start` |
| `score_format` | WARNING | parse | score not `.`/float |
| `strand_valid` | ERROR | parse | strand ∉ `{+ - .}` |
| `cds_phase_required` | ERROR | parse | CDS phase == `.` |
| `cds_phase_value` | ERROR | parse | CDS phase ∉ `{0,1,2}` |
| `non_cds_phase` | WARNING | parse | non-CDS phase ≠ `.` |
| `attr_format` | ERROR | parse | attribute has no `=` |
| `attr_empty_key` | WARNING | parse | empty attribute key |
| `missing_id` | ERROR | parse | gene/transcript lacks ID |
| `duplicate_id` | ERROR | validate_gff3_file | ID seen twice |
| `orphan_parent` | ERROR | hierarchy | Parent not a known ID |
| `gene_has_parent` | ERROR | hierarchy | gene has a Parent |
| `transcript_no_parent` | ERROR | hierarchy | transcript lacks Parent |
| `transcript_parent_type` | ERROR | hierarchy | transcript parent not a gene |
| `exon_no_parent` | ERROR | hierarchy | exon lacks Parent |
| `exon_parent_type` | ERROR | hierarchy | exon parent not transcript |
| `cds_no_parent` | ERROR | hierarchy | CDS lacks Parent |
| `cds_parent_type` | ERROR | hierarchy | CDS parent not transcript |
| `transcript_has_exons` | WARNING | hierarchy | transcript has no exon/CDS children |
| `mrna_has_cds` | WARNING | hierarchy | mRNA has no CDS children |
| `gene_has_transcripts` | WARNING | hierarchy | gene has no children |
| `seqid_consistency` | ERROR | containment | child seqid ≠ parent seqid |
| `strand_consistency` | ERROR | containment | child strand ≠ parent strand (both definite) |
| `coord_containment` | ERROR | containment | child extends outside parent |
| `cds_phase_consistency` | WARNING | cds_phase | computed phase mismatch (non-first CDS) |
| `lifton_source` | INFO | lifton_attrs | unexpected source column |
| `lifton_protein_identity` | ERROR | lifton_attrs | protein_identity not float in [0,1] |
| `lifton_dna_identity` | ERROR | lifton_attrs | dna_identity not float in [0,1] |
| `lifton_attrs_present` | WARNING | lifton_attrs | coding transcript missing protein_identity / any transcript missing dna_identity |

#### 8.2.10 Report formatting (`print_validation_report`, `:280-332`)

Prints a box-drawing report to **stderr**, width `w=72`. Header line shows `GFF3 VALIDATION REPORT — ✅  VALID` if `result.is_valid` else `❌  INVALID`. Then: file path, line counts, a `FEATURE COUNTS:` block (`sorted(result.stats.items())`), then `Errors : N` / `Warnings : N`. The `ISSUES:` block shows `result.errors + (result.warnings if verbose else [])`, grouped by `check` name (a blank line + `[severity] Check: <name>` header whenever `issue.check` changes), each issue wrapped to `w-6` via `_wrap` (`:335`, a hard character chop, not word-aware). If not verbose and warnings exist, a trailing hint line is printed. `_wrap` splits text into fixed-width chunks.

#### 8.2.11 Exit codes & pipeline integration

- **Standalone `_main` (`:959-1030`):** argparse flags `gff3_file` (positional), `-v/--verbose`, `--no-hierarchy`, `--no-phase`, `--no-containment`, `--no-lifton-attrs`, `--max-issues N` (default 50), `--json`. The four `--no-*` flags map to `check_*=not args.no_*`. If `--json`, dumps `{file, is_valid, stats, issues:[…]}` to **stdout** with `indent=2`; else prints the box report. **Exit code: `sys.exit(0 if result.is_valid else 1)`** — i.e. **1** on any ERROR (output-side uses 1; input-side strict uses 2).
- **Pipeline Step 10 (`lifton.py:585-603`):** runs only when `--validate-output` set, against `args.output` with all four checks enabled. Prints a banner, the box report (`verbose=args.validate_verbose`), and if `not is_valid` an extra `[LiftOn] Output GFF3 has N error(s).` line. **Gotcha:** in the pipeline path the result is reported but **does not change LiftOn's exit code** — `sys.exit` is not called here, so an invalid output GFF3 is advisory only when run via `--validate-output`.

---

### 8.3 Pre-flight annotation validator (`annotation_validator.py`)

A **GFF3/GTF input sanity checker** that runs **before any gffutils database build**, giving an actionable error instead of a cryptic SQLite `UNIQUE` failure. It is a third, separate implementation.

#### 8.3.1 When it runs

1. **`Annotation.__init__` path:** `Annotation.__init__` calls `self._run_preflight_validation()` (`annotation.py:139`) before `self._get_db_connection()` (`:142`). `_run_preflight_validation` (`annotation.py:262-295`) calls `validate_annotation_file(self.file_name, max_duplicate_examples=20, check_orphan_parents=True)`.
2. **Vendored Liftoff path:** `liftoff/extract_features.py:83` calls `validate_annotation_file(annotation_file, max_duplicate_examples=20)` (default `check_orphan_parents=True`) before its own gffutils build.

Both run unconditionally (no flag toggles them).

#### 8.3.2 `ValidationResult` (`annotation_validator.py:20-33`)

| Field | Type | Default | Meaning |
|---|---|---|---|
| `file_path` | `str` | — | input |
| `is_valid` | `bool` | `True` | `False` → not 100% clean; **not** always fatal |
| `file_format` | `str` | `"unknown"` | `"GFF3"` / `"GTF"` / `"unknown"` |
| `total_lines`, `data_lines`, `comment_lines` | `int` | `0` | counts |
| `duplicate_ids` | `Dict[str, List[int]]` | `{}` | id → line numbers (only ids seen >1×) |
| `orphan_parents` | `List[str]` | `[]` | up to 50 unresolved Parent ids |
| `warnings`, `errors`, `fix_suggestions` | `List[str]` | `[]` | human-readable strings (not structured findings) |

#### 8.3.3 `validate_annotation_file` algorithm (`:40-245`), numbered

1. **File checks (each sets `is_valid=False` and returns early):** not `os.path.exists` → error "File not found"; not `os.path.isfile` → "Path is not a regular file"; `os.path.getsize==0` → "File is empty"; not `os.access(..., os.R_OK)` → "File is not readable".
2. **Single streaming pass** (`open(..., errors="replace")`), per non-comment/non-blank line:
   - skip lines starting `#` or blank (counted as comment).
   - `cols = line.split("\t")`; if `!= 9` → record lineno (first 5 only) in `bad_column_count_lines` and skip (not counted as data).
   - else `data_lines += 1`; `attrs = cols[8]`.
   - **Format vote:** if `any(k in attrs for k in {"ID=","Parent=","Name="})` and `"=" in attrs` → `gff_count += 1`; elif regex `(gene_id|transcript_id)\s+"` matches → `gtf_count += 1`.
   - **ID extraction:** regex `(?:^|;)ID=([^;]+)` → append lineno to `id_linenos[feat_id]`, add to `all_ids`.
   - **Parent extraction:** regex `(?:^|;)Parent=([^;]+)`; split on `","`, strip, add each to `all_parent_refs`.
3. **Format determination:** if `gff_count>0 and gff_count>=gtf_count` → `"GFF3"`. Elif `gtf_count>0` → `"GTF"` (+ warning about GTF auto-conversion + a `gffread`/`agat` fix suggestion). Else `"unknown"` (+ warning, will try as GFF3).
4. **Structural sanity:** if `data_lines == 0` → `is_valid=False`, error "File contains no valid 9-column data lines", fix suggestion, return. If `bad_column_count_lines` non-empty → warning listing the first 5.
5. **Duplicate IDs:** `duplicates = {id: linenos for ... if len(linenos)>1}`. If any: warning with up to `max_duplicate_examples=20` examples (each showing first 3 line numbers, `…` if more), a fix suggestion mentioning RefSeq duplicate CDS IDs + the auto-retry unique-ID transform + `grep` commands for the first 5, and **set `is_valid=False`**. **Gotcha:** duplicate IDs are a WARNING that nevertheless flips `is_valid` to `False` — this signals the caller to use the fallback DB-build strategy, **not** to abort.
6. **Orphan parents** (only if `check_orphan_parents` and `file_format=="GFF3"`): `orphans = [p for p in all_parent_refs if p and p not in all_ids]`, capped at 50 into `orphan_parents`. If any → warning (first 5) + fix suggestion with `grep` for first 3. **Gotcha:** orphan parents do **not** set `is_valid=False`.

#### 8.3.4 Fatal vs recoverable in the caller

`_run_preflight_validation` (`annotation.py:262-295`): prints the report if any errors/warnings; then defines `fatal_keywords = ["not found", "not readable", "empty", "no valid 9-column", "not a regular file"]` and sets `truly_fatal = any(kw in err.lower() for err in result.errors for kw in fatal_keywords)`. If `truly_fatal` → `logger.log_error(...)` + `sys.exit(1)`. **Gotcha:** the decision keys on **error-string substrings**, not on `is_valid`. So `is_valid=False` caused solely by duplicate IDs is **not** fatal (the duplicate-ID message contains none of the keywords); only the four file-level errors abort. The Liftoff path (`extract_features.py:88-99`) uses the identical keyword list and `sys.exit(1)`.

#### 8.3.5 Reporting helpers (`annotation_validator.py:252-367`)

- `print_validation_report(result, always_show=False)` (`:252`): returns early if `is_valid and not always_show and not warnings`. Prints a box-drawing report to stderr (width 70) with a red `❌  ANNOTATION FILE ERROR` header if `result.errors` else a yellow `⚠️  ANNOTATION FILE WARNING` header, then File/Format/Lines, then `ERRORS:`, `WARNINGS:`, `SUGGESTED FIXES:` blocks, each line wrapped by the module-local **word-aware** `_wrap` (`:323`, splits on spaces, preserves embedded `\n`). Output lines are truncated to `width+4`.
- `print_db_build_error(file_path, strategy, exception)` (`:344`): boxed `❌  gffutils DB BUILD FAILED (strategy=…)` with the exception text wrapped.
- `print_db_build_success(file_path, strategy)` (`:359`): a single `✓ gffutils database built successfully [file=…, strategy=…]` line to stderr.

These three helpers are invoked around each of the 3 progressive gffutils build strategies in `annotation._build_database` (`annotation.py:377,381,407,415,434,441`) and the 2-strategy build in `liftoff/extract_features.py:116,120,142,150`, narrating each attempt.
