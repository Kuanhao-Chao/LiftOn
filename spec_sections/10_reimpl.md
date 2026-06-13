## 10. Reimplementation Guidance

This section is the operational closeout of the specification. It tells a developer **in what order** to build the component, **which exact behaviours must be preserved byte-for-byte**, and **how to prove the reimplementation is correct**. The four appendices are the consolidated reference tables (CLI flags, numeric constants, the status/mutation vocabularies, and the file-tree map) that every preceding section refers to.

The single non-negotiable correctness contract is the **24-cell byte-identity matrix** (`tests/test_native_matrix.py::TestFullNativeMatrix::test_all_24_combinations_byte_identical`): every combination of `--stream × --inmemory-liftoff × --threads ∈ {1,2,4} × --native` (2·2·3·2 = 24 cells) must emit the same bytes as the default path. Every "Gotcha:" below is a way to break that gate.

---

### 10.1 Recommended build order

Build bottom-up. Each stage has a **testable milestone** that can be validated in isolation before the next stage is started. The arrows show the hard data dependencies.

```
data model → annotation backend → sequence extraction → alignment kernel
   → variant classification → chaining → ORF rescue → writers
   → orchestration → parallelism
```

#### Stage 1 — Object model (`Lifton_GENE → Lifton_TRANS → Lifton_EXON → Lifton_CDS`, `Lifton_Alignment`, `Lifton_Status`, `Lifton_ORF`)
Implement the class hierarchy in `lifton/lifton_class.py` (`Lifton_GENE` 57-217, `Lifton_TRANS` 267-820, `Lifton_EXON` 826-873, `Lifton_CDS` 876-895). These are thin wrappers around a gffutils-`Feature`-shaped record (`seqid, source, featuretype, start, end, score, strand, frame, attributes` where `attributes` is a `dict[str, list[str]]`). Implement `Lifton_Status` first (`lifton_class.py:12-21`): float fields `liftoff`, `miniprot`, `lifton_dna`, `lifton_aa`, `eval_dna`, `eval_aa` (all default `0`); `annotation: str|None` (default `None`); `status: list[str]` (default `[]`). There is no `parent` field.
- **Milestone:** Construct a single-CDS `Lifton_GENE`, call `write_entry`, and emit a syntactically valid 1-gene/1-mRNA/1-exon/1-CDS GFF3 block.

#### Stage 2 — Annotation backend (`lifton/annotation.py` + `gffbase_adapter.py`)
A read-only `Annotation` over a feature DB. Two interchangeable backends: gffutils-SQLite (default) and gffbase-DuckDB (`LIFTON_USE_GFFBASE=1`). The only methods the rest of the pipeline calls are `db[feature_id]` (returns a Feature), `db.children(feature, featuretype=..., level=..., order_by=...)`, and iteration. **The `children()` signatures are part of the byte-identity contract** (see §10.2, Gotcha PROXY).
- **Milestone:** Round-trip an input GFF3 → DB → feature lookup → children iteration, matching gffutils' ordering.

#### Stage 3 — Sequence extraction (`lifton/extract_sequence.py`)
Stream-extract `transcripts.fa` and `proteins.fa` to `intermediate_files/`, then re-open via `pyfaidx.Fasta` for lazy mmap access (`lifton.py` Step 3, ~382-408). Translation uses the standard genetic code via `Bio.Seq.Seq.translate()`. Coordinates are **1-based inclusive** throughout (see §10.2, Gotcha COORD).
- **Milestone:** For the synthetic fixture chromosome (positions 101-199 `ATG…`, 301-399 `…TAA`), extract a transcript whose translation is a clean ORF ending in `*`.

#### Stage 4 — Alignment kernel (`lifton/align.py`, `lifton/get_id_fraction.py`)
This is the most byte-sensitive stage. Implement exactly:
1. `parasail_align_protein_base` (`align.py:84-110`): `parasail.Matrix("blosum62")`, `gap_open=11`, `gap_extend=1`, kernel `parasail.nw_trace_scan_sat`. Raise `LiftOnAlignmentError` on empty query or ref.
2. `parasail_align_DNA_base` (`align.py:134-155`): `parasail.matrix_create("ACGTN*", 1, -3)`, `gap_open=5`, `gap_extend=2`, same kernel. **Sanitise both inputs** through `_sanitise_for_parasail_dna` (`align.py:16-31`): any character not in `frozenset("ACGTN*")` (after `.upper()`) → `N`.
3. Identity = `matches/length` from `get_id_fraction.get_AA_id_fraction` / `get_DNA_id_fraction` over the traceback `ref`/`query` strings.
- **Milestone:** Two identical protein sequences → `identity == 1.0`; a single substitution → `< 1.0`. An IUPAC `R` in the DNA input does not crash the kernel.

#### Stage 5 — Variant classification (`lifton/variants.py`)
Implement `find_variants` (`variants.py:45-120`) exactly as the 9-branch decision tree in §10.2 / Appendix C. Helpers `has_stop_codon` (1-15) and `is_frameshift` (18-42).
- **Milestone:** Identical DNA alignment → `['identical']`; identical protein but non-identical DNA → `['synonymous']`; a `'-'` run of length not divisible by 3 → `['frameshift']`.

#### Stage 6 — Protein-maximization chaining (`lifton/protein_maximization.py`)
`chaining_algorithm` (`protein_maximization.py:176`) walks Liftoff and miniprot CDS blocks synchronized by accumulated reference-AA counts, calling `push_cds_idx` per chunk to emit the higher-identity tool's CDS list. Honour the empty-input guards, the both-zero `empty[...]` label, and **Liftoff-wins-on-tie** (see §10.2, Gotchas TIE and BOTHZERO).
- **Milestone:** When miniprot is absent → return Liftoff's CDS list unchanged. When both tools score 0 on a chunk → that chunk contributes no CDS and logs `empty[...]`.

#### Stage 7 — ORF rescue + CDS-boundary patching (`Lifton_TRANS.__find_orfs` 645-, `__update_cds_boundary` 718-, `__get_cds_frame` 791-792)
Scan the spliced transcript in all 3 frames; keep one best (longest) ORF per frame; align each candidate's translation to the reference protein; apply the winner **only if** `max_identity > lifton_status.lifton_aa + 0.01` (`lifton_class.py:713-714`). Recompute CDS `frame` as `(3 - accum_cds_length % 3) % 3` (see §10.2, Gotchas ORFAPPLY and FRAME).
- **Milestone:** An ORF improving identity by 2% is applied; one improving by 0.5% is not.

#### Stage 8 — Writers (`lifton/io/gff3_writer.py`, `Lifton_*.write_entry`)
Implement `format_directives` (`gff3_writer.py:84-111`), `format_attributes` (114-128), and `_canonical_attr_order` (70-81). Attribute order is **ID, Parent, then alphabetical** (see §10.2, Gotcha ATTRORDER). Percent-encode reserved characters. The four-level write chain is `Lifton_GENE.write_entry → Lifton_TRANS.write_entry → Lifton_EXON.write_entry → Lifton_CDS.write_entry`.
- **Milestone:** Re-serialising a feature with attributes `{Name, ID, Parent, x}` yields `ID=…;Parent=…;Name=…;x=…`.

#### Stage 9 — Orchestration (`lifton/lifton.py:run_all_lifton_steps` 283-599)
Wire the 11 steps + strict-GFF gate sequentially: open FASTAs → optional strict gate → build ref `Annotation` → partition coding/non-coding → extract sequences → run Liftoff (in-process) + miniprot (subprocess) → reload GFFs and open writers → emit directive prologue → walk Liftoff genes (Step 7) → walk miniprot-only mRNAs (Step 8) → stats → optional output validation.
- **Milestone:** End-to-end run with `-L`/`-M` pre-baked GFFs produces a complete GFF3 without touching minimap2/miniprot.

#### Stage 10 — Parallelism (`lifton/parallel.py`, `lifton/locus_pipeline.py`)
Dispatch Step 7 through a `ThreadPoolExecutor` (`--locus-pipeline` + `--threads N`). Each worker returns a submission-indexed `LocusResult`; a heap-backed `OrderedWriter` (`parallel.py:43+`) emits in submission order. Build the proxy DBs (`_RefDbProxy`, `_LFeatureDbProxy`, `_MFeatureDbProxy`) with the exact `children()` signatures (see §10.2, Gotcha PROXY). This stage MUST NOT change output bytes versus Stage 9.
- **Milestone:** `--threads 4` output is `diff`-identical to `--threads 1`; the full 24-cell matrix passes.

---

### 10.2 Determinism-critical invariants ("gotchas" catalog)

Every item is a behaviour that, if implemented differently, silently changes the output bytes and trips the 24-cell gate. The reimplementer MUST honour each one.

| Tag | Invariant | Why it matters | Source |
|---|---|---|---|
| **SUBIDX** | Every locus carries an integer **submission index** assigned in walk order; results are buffered and emitted strictly by ascending index, never by completion order. | Threads finish out of order; emitting by completion order scrambles the GFF3. | `lifton/parallel.py:62-64` (`pending` min-heap keyed by `result.index`, `next_to_emit`) |
| **HEAP** | The ordered writer emits only the **contiguous prefix**: pop `pending[0]` iff `pending[0][0] == next_to_emit`, increment, repeat; otherwise hold. Spilled future indices are restored from disk in the same order. | A gap at index *k* must block *k+1*; emitting eagerly reorders genes. | `lifton/parallel.py:121-139` (`offer`/`_drain_contiguous`) |
| **PROXY** | Worker DB proxies serve only the **exact pre-cached `children()` signatures**; any other signature raises `NotImplementedError`. The four cached `_LFeatureDbProxy` shapes: `(featuretype="exon", level=1)`, `(featuretype="exon", level=None)`, `(featuretype=("CDS","stop_codon"))`, `(featuretype=None, level=1)`. `_MFeatureDbProxy` serves only `(featuretype=("CDS","stop_codon"))`. | The proxy must return children in the **identical order** the real DB would; a missing signature would silently fall through to `iter([])` and drop features if not guarded. | `lifton/locus_pipeline.py:209-235` (L proxy), `255-267` (M proxy) |
| **FIDCONS** | Proxy lookups key on `getattr(feature, "id", feature)` — the feature's `id` must equal the cache key used at materialisation time. | A mismatch between the id used to cache and the id used to look up yields `iter([])` (empty children) and drops the gene's exons. | `lifton/locus_pipeline.py:213`, `257` |
| **EXCEPT** | Worker bodies catch `except Exception`, **never** `except BaseException`. | `BaseException` would swallow `KeyboardInterrupt`/`SystemExit`/`GeneratorExit`; the narrowed catch lets the orchestrator shut down cleanly and matches the serial path's error surface. | `lifton/locus_pipeline.py:100` ("narrowed from BaseException"), `359`, `368`, `393`, `404` |
| **TIE** | In per-chunk chaining, **Liftoff wins ties**: take miniprot's CDS only when `m_identity > l_identity`; the `else` branch (including equality) keeps Liftoff. | Equal-identity chunks must deterministically pick one tool; choosing miniprot on tie would mutate CDS boundaries. | `lifton/protein_maximization.py:110-122` |
| **BOTHZERO** | When `m_identity == 0.0 and l_identity == 0.0` on a chunk, emit **no CDS** (`return []`) and append the label `empty[start-end]`, **not** `liftoff[0.00-0.00]`. | The both-zero case is provenance-ambiguous; the `empty[...]` label and empty CDS list are pinned. | `lifton/protein_maximization.py:104-108` |
| **CHAINGUARD** | If either alignment has 0 CDS children, return Liftoff's CDS list as-is (or empty if Liftoff is also empty); never divide by zero on a zero-length protein window — fall back to Liftoff. | Degenerate inputs must not raise; the fallback target is always Liftoff. | `lifton/protein_maximization.py:75-80`, `189-206` |
| **ORFAPPLY** | An ORF-rescue candidate is applied **only if** `final_orf is not None and max_identity > (lifton_status.lifton_aa + threshold_orf)` where `threshold_orf = 0.01`. A `>` (strict) comparison, not `>=`. | Applying an ORF that improves by exactly 0.01 or less would mutate CDS boundaries; the +1% floor is load-bearing. | `lifton/lifton_class.py:713-714` |
| **FRAME** | CDS `frame` (GFF3 phase) is computed as `(3 - accum_cds_length % 3) % 3`, where `accum_cds_length` is the running sum of preceding CDS lengths in transcription order. The outer `% 3` maps the `accum%3 == 0` case to `0` (not `3`). | An off-by-one in the phase column changes every multi-CDS feature line. | `lifton/lifton_class.py:791-792`, applied at `577`, `758`, `784` |
| **ATTRORDER** | Column-9 attribute order is **`ID` first, `Parent` second, then all remaining keys sorted alphabetically (case-sensitive Python `sorted`)**. Values that are `None` or empty lists are skipped; list values are joined with `,`. | Any other ordering changes every feature line. | `lifton/io/gff3_writer.py:70-81` (`_canonical_attr_order`), `114-128` (`format_attributes`) |
| **COORD** | All coordinates are **1-based inclusive**. CDS/exon length is `end - start + 1`. The miniprot trans-ratio uses `(mtrans.end - mtrans.start + 1) / ref_features_len_dict[ref_gene_id]`. | A 0-based/half-open mix-up shifts every coordinate and every derived length. | `lifton/run_miniprot.py:306`, `lifton/lifton_class.py:728` (`exon.entry.end - exon.entry.start + 1`) |
| **STRAND** | Strand reversal points: ORF/CDS boundary updates iterate exons in **forward order on `+`** and **reversed order (`exons[::-1]`) on `-`**. Sequence assembly reverse-complements on the minus strand before translation. | Iterating the wrong direction on `-` strand assigns CDS pieces to the wrong exons. | `lifton/lifton_class.py:718-722` (`__update_cds_boundary`) |
| **DIRECTIVES** | The output begins with the directive block from `format_directives`: always `##gff-version 3` first, then every other input `##`/`#!` directive in **input order**, de-duplicated; the whole block ends with exactly one `\n`, written **before any feature row**. | A missing or reordered directive prologue changes the first bytes of the file. | `lifton/io/gff3_writer.py:84-111` |
| **DNASANITISE** | DNA inputs to parasail are upper-cased and every non-`ACGTN*` byte → `N` **before** alignment. | Changing the sanitisation set or skipping `.upper()` changes alignment identity on any sequence with ambiguity codes or lowercase. | `lifton/align.py:13`, `16-31`, `151-152` |
| **EMPTYPROT** | `parasail_align_protein_base` **raises** `LiftOnAlignmentError` on an empty query or reference; it does not silently coerce to `"*"`. | Silent coercion would hide a zero-CDS bug behind a near-zero identity and change downstream chaining decisions. | `lifton/align.py:102-107` |
| **MINIPROTRATIO** | A miniprot-only transcript is accepted only if its CDS count is not 1-while-ref->1 (processed-pseudogene reject) **and** `args.min_miniprot < ratio < args.max_miniprot` (strict both sides; defaults 0.9 and 1.5). The accepted ratio is written as `miniprot_annotation_ratio` formatted `f"{ratio:.3f}"`. | The strict `<`/`>` bounds and the 3-decimal format are pinned in output attributes. | `lifton/run_miniprot.py:304-309` |
| **MP2DEFAULTS** | `parse_args` post-processes `mm2_options`: each of `-a`, `--eqx`, `-N 50`, `-p 0.5`, `--end-bonus 5` is appended **only if its flag substring is absent**. | Different default mm2 options change Liftoff's DNA alignments and the final lift-over. | `lifton/lifton.py:266-275` |
| **RECURSION** | Before running Liftoff, raise the recursion limit to `max(original, 10000)` and **restore it afterward**. | Deep gene hierarchies overflow the default limit; not restoring leaks state across runs. | `lifton/run_liftoff.py:56-57`, `83` |
| **VARORDER** | `find_variants` returns early on the first matching of: non-coding → full_transcript_loss → no_protein → identical (`dna.identity == 1.0`) → synonymous (`protein.identity == 1.0`); only then accumulates the multi-label set. The status list ordering is **fixed by evaluation order**. | Re-ordering the checks changes the emitted `status=` attribute. | `lifton/variants.py:69-120` |

---

### 10.3 Test strategy

Mirror the existing pytest harness exactly; it is the only evidence the reimplementation is correct.

1. **The 24-cell byte-identity matrix (primary gate).** Re-run the full pipeline under all 24 combinations of `--stream × --inmemory-liftoff × --threads ∈ {1,2,4} × --native` and assert every output is byte-for-byte equal to the default-path output. The reference file in the existing suite is a **391-byte** synthetic GFF3. Mirror at `tests/test_native_matrix.py::TestFullNativeMatrix::test_all_24_combinations_byte_identical`. Subset gates (build these first, then the full matrix): `tests/test_pipeline_streaming.py` (`--stream`), `tests/test_liftoff_inmemory.py` (4-cell `--inmemory-liftoff`), `tests/test_parallelism_matrix.py` (`--threads`).

2. **Hermetic `-L`/`-M` golden-path.** Supply pre-baked Liftoff and miniprot GFFs via `-L`/`-M`; `exec_liftoff`/`exec_miniprot` short-circuit on `os.path.exists`, so no minimap2/miniprot binary is needed. The integration fixture (`tests/test_integration_pipeline.py`, fixture `hermetic_pipeline`) monkey-patches `run_liftoff.run_liftoff` and `run_miniprot.run_miniprot` to **raise if invoked** — any accidental fall-through to the external runners is a hard test failure. Use this same trick to author new integration scenarios.

3. **Synthetic fixtures.** `tests/conftest.py` builds 600-bp synthetic chromosomes whose positions 101-199 (`ATG…`) and 301-399 (`…TAA`) form a clean ORF; reuse these rather than authoring FASTA from scratch. The minimal end-to-end output is the **391-byte** GFF3 referenced by the matrix test.

4. **Pytest layout.** Unit tests per module (`test_align.py`, `test_variants.py`, `test_protein_maximization.py`, `test_lifton_class.py`, …); integration tests in `test_integration_pipeline.py`; matrix gates in the four files above; the byte-identity gate in `test_native_matrix.py`. Opt-in perf/memory profiling under `tests/perf/` is gated by the `perf` marker and never runs by default. Property-based tests (`test_property_based.py`, `test_streaming_property.py`) require `hypothesis`; the rest of the suite needs only the conda deps. Run the whole suite with `pytest tests/ -v`; run the load-bearing gate alone with `pytest tests/test_native_matrix.py::TestFullNativeMatrix::test_all_24_combinations_byte_identical -v`. CI does **not** run pytest — it is the correct local pre-merge gate.

---

### Appendix A — Complete CLI reference

Console entry points (`setup.py`): `lifton = lifton.lifton:main`, `gff3-validate = lifton.gff3_validator:_main`. Defined in `parse_args` (`lifton/lifton.py:209-280`) and its helper groups. Positionals first: `target` (target FASTA), `reference` (reference FASTA).

| Flag(s) | Default | Type / choices | Effect |
|---|---|---|---|
| `target` (positional) | — | path | Target genome FASTA to lift genes **to**. |
| `reference` (positional) | — | path | Reference genome FASTA to lift genes **from**. |
| `-g`, `--reference-annotation` | — (**required**) | GFF/GTF path or DB name | Reference annotation to lift over. GTF auto-detected and converted to GFF3 unless `--no-auto-convert-gtf`. |
| `-o`, `--output` | `lifton.gff3` | FILE | Output GFF3. `stdout` is a recognised special value (outdir = `.`). |
| `-u` | `unmapped_features.txt` | FILE | Write unmapped features here. |
| `-exclude_partial` | off | store_true | Send sub-threshold partial mappings to the unmapped file instead of annotating them. |
| `-mm2_options` | `-a --end-bonus 5 --eqx -N 50 -p 0.5` | str | minimap2 params; missing defaults re-appended post-parse (Gotcha MP2DEFAULTS). |
| `-mp_options` | `''` (empty) | str | miniprot params, space-delimited; appended after `--gff-only`. |
| `-a` | `0.5` | float | Min alignment coverage for a feature to be mapped. |
| `-s` | `0.5` | float | Min child-feature sequence identity for a feature to be mapped. |
| `-min_miniprot` | `0.9` | float | Min length ratio for a miniprot-only transcript (strict lower bound). |
| `-max_miniprot` | `1.5` | float | Max length ratio for a miniprot-only transcript (strict upper bound). |
| `-d` | `2.0` | float | Distance scaling factor for the alignment graph (Liftoff). |
| `-flank` | `0` | float | Flanking fraction [0-1] of gene length to align. |
| `-V`, `--version` | — | action=version | Print version and exit. |
| `-D`, `--debug` | off | store_true | Debug mode. |
| `-t`, `--threads` | `1` | int | Worker thread count (Step 7 fan-out with `--locus-pipeline`). |
| `-m` | — | path | minimap2 path override. |
| `-f`, `--features` | — | TYPES | List of feature types to lift over. |
| `-infer-genes` | off | store_true | Infer genes when input has only transcripts/exon/CDS (auto for GTF). |
| `-infer_transcripts` | off | store_true | Infer transcripts when input has only exon/CDS (auto for GTF). |
| `-chroms` | — | TXT | Comma-separated reference,target chromosome pairs. |
| `-unplaced` | — | TXT | Unplaced-sequence names to map after chroms; **requires** `-chroms`. |
| `-copies` | off | store_true | Look for extra gene copies in target. |
| `-sc` | `1.0` | float | With `-copies`, min identity for a copy; must be ≥ `-s` (validated). |
| `-overlap` | `0.1` | float | Max overlap fraction [0-1] between two features. |
| `-mismatch` | `2` | int | Exon mismatch penalty (best-mapping search). |
| `-gap_open` | `2` | int | Exon gap-open penalty. |
| `-gap_extend` | `1` | int | Exon gap-extend penalty. |
| `-subcommand` | — | str (SUPPRESS) | Hidden internal subcommand. |
| `-polish` | off | store_true | Polishing pass. |
| `-cds` | `True` | store_true | Annotate per-CDS status (partial/missing start/missing stop/inframe stop). |
| `-time`, `--measure_time` | off | store_true | Per-step timing → `time.txt`. |
| `--validate-output` | off | store_true | Re-validate output GFF3 after writing (Phase 13.5C). |
| `--validate-verbose` | off | store_true | With `--validate-output`, also print warnings. |
| `--strict-gff` | off | store_true | Run NCBI GFF3 input validator on the reference; exit non-zero on any violation. |
| `--stream` | off | store_true | Pipe miniprot stdout into in-memory gffbase DB; skip `miniprot.gff3` disk round-trip. |
| `--inmemory-liftoff` | off | store_true | Serialise Liftoff `lifted_feature_list` in-process; skip `liftoff.gff3` disk write. |
| `--locus-pipeline` | off | store_true | Dispatch Step 7 through a `ThreadPoolExecutor` (needs `--threads N>1` for parallelism). |
| `--native` | off | store_true | Drive minimap2 via `mappy` and miniprot via the pyminiprot facade in-process. |
| `-E`, `--evaluation` | off | store_true | Evaluation mode. |
| `-EL`, `--evaluation-liftoff-chm13` | off | store_true | Evaluation mode (Liftoff CHM13). |
| `-c`, `--write_chains` | `True` | store_true | Write chaining files. |
| `--no-orf-search` | off | store_true | Skip ORF-rescue pass. |
| `-P`, `--proteins` | `None` | FASTA | Reference protein FASTA (else extracted from annotation). |
| `-T`, `--transcripts` | `None` | FASTA | Reference transcript FASTA (else extracted). |
| `-L`, `--liftoff` | `None` | gff/db | Pre-built Liftoff annotation (short-circuits the Liftoff run). |
| `-M`, `--miniprot` | `None` | gff/db | Pre-built miniprot annotation (short-circuits the miniprot run). |
| `--merge-strategy` | `create_unique` | choice | gffutils DB merge strategy (`create_unique\|merge\|error\|warning\|replace`). |
| `--id-spec` | `None` | str | Attribute used as feature ID (default `ID`). |
| `--force` | off | store_true | Overwrite existing DB. |
| `--verbose` | off | store_true | Verbose gffutils output. |
| `-ad`, `--annotation-database` | `RefSeq` | str | Reference annotation source (`RefSeq`/`GENCODE`/`Ensembl`/`CHESS`/other); affects ID extraction. |
| `--no-auto-convert-gtf` | off | store_true | Disable automatic GTF→GFF3 conversion. |

**Post-parse validation** (`lifton.py:266-279`): mm2 default re-append (Gotcha MP2DEFAULTS); `parser.error` if `float(args.s) > float(args.sc)`; `parser.error` if `chroms is None and unplaced is not None`.

---

### Appendix B — Constants & default-parameter table

| Constant | Value | Where | Role |
|---|---|---|---|
| Protein parasail matrix | `parasail.Matrix("blosum62")` | `align.py:95` (protein base; called via `protein_align` at `:124`) | Substitution scores for protein NW alignment. |
| Protein gap open | `11` | `align.py:96` | Affine gap-open for protein alignment. |
| Protein gap extend | `1` | `align.py:97` | Affine gap-extend for protein alignment. |
| DNA parasail matrix | `parasail.matrix_create("ACGTN*", 1, -3)` | `align.py:148` | Match `+1`, mismatch `-3` over the 6-char DNA alphabet. |
| DNA gap open | `5` | `align.py:149` | Affine gap-open for DNA alignment. |
| DNA gap extend | `2` | `align.py:150` | Affine gap-extend for DNA alignment. |
| parasail kernel | `parasail.nw_trace_scan_sat` | `align.py:109`, `154` | Global (Needleman-Wunsch) traceback kernel, saturated, for both DNA and protein. |
| DNA alphabet (sanitise set) | `frozenset("ACGTN*")` | `align.py:13` | Chars passed through; everything else → `N`. |
| `threshold_orf` | `0.01` | `lifton_class.py:713` | ORF rescue applied only if identity gain `> 1%` (strict). |
| Recursion limit | `max(orig, 10000)` | `run_liftoff.py:56-57` | Raised before Liftoff; restored after (line 83). |
| Prefetcher pool cap | `min(4, threads, …)` | `parallel.py:351` | Max concurrent locus-prefetch threads. |
| Stream drain chunk | `1 << 16` (65536 B) | `run_miniprot.py:6,171` | Bounded per-stream read buffer when draining miniprot stdout/stderr. |
| Popen bufsize (stream) | `1 << 20` | `run_miniprot.py:168` | Pipe buffer for streamed miniprot. |
| mm2 default options | `-a --end-bonus 5 --eqx -N 50 -p 0.5` | `lifton.py:66-69`, re-append `266-275` | Default minimap2 parameters for Liftoff DNA alignment. |
| `-a` (coverage) | `0.5` | `lifton.py:73` | Min alignment coverage threshold. |
| `-s` (identity) | `0.5` | `lifton.py:77` | Min child sequence identity threshold. |
| `min_miniprot` | `0.9` | `lifton.py:82` | Lower length-ratio bound, strict. |
| `max_miniprot` | `1.5` | `lifton.py:86` | Upper length-ratio bound, strict. |
| `-d` | `2.0` | `lifton.py:90` | Distance scaling factor. |
| `-overlap` | `0.1` | `lifton.py:136` | Max feature-overlap fraction. |
| `-mismatch`/`-gap_open`/`-gap_extend` | `2`/`2`/`1` | `lifton.py:138-143` | Exon best-mapping penalties (Liftoff). |
| `-sc` | `1.0` | `lifton.py:132` | Copy-detection identity floor. |
| CDS frame formula | `(3 - accum_cds_length % 3) % 3` | `lifton_class.py:792` | GFF3 phase column for each CDS. |
| miniprot ratio format | `f"{ratio:.3f}"` | `run_miniprot.py:309` | 3-decimal `miniprot_annotation_ratio` attribute. |

---

### Appendix C — Status vocabularies

#### C.1 `lifton_status.annotation` (provenance, single string)
Set per-transcript to record which engine produced the final annotation.

| Value | Trigger | Source |
|---|---|---|
| `Liftoff` | Liftoff transcript processed normally. | `run_liftoff.py:246` |
| `no_ref_protein` | Liftoff alignment is `None` (no reference protein available). | `run_liftoff.py:172` |
| `LiftOn_chaining_algorithm` | Liftoff identity `< 1` and a valid miniprot alignment exists → chaining runs. | `run_liftoff.py:179` |
| `miniprot` | Transcript came from the miniprot-only Step 8 path. | `run_miniprot.py:281` |

Note: when Liftoff identity `== 1` exactly, `lifton_status.lifton_aa = 1` is set and `annotation` retains its prior value (`Liftoff`); the chaining branch is skipped (`run_liftoff.py:173-175`).

#### C.2 `lifton_status.status` — the mutation/variant vocabulary
`find_variants` (`variants.py:45-120`) sets `lifton_status.status` to a **list** of labels. Early-return labels are singletons; the late block can accumulate several. `align_dna` = the DNA pairwise alignment, `align_protein` = the protein pairwise alignment, `peps` = the protein split on `*`.

| # | Label | Triggering condition | Source line |
|---|---|---|---|
| pre | `non_coding` | `is_non_coding` is true. | `variants.py:69-72` |
| pre | `full_transcript_loss` | `align_dna is None`. | `73-76` |
| pre | `no_protein` | `align_protein is None`. | `77-80` |
| 1 | `identical` | `align_dna.identity == 1.0`. | `82-85` |
| 2 | `synonymous` | `align_protein.identity == 1.0` (DNA not identical). | `88-91` |
| 3 | `frameshift` | `is_frameshift(align_dna.query_aln)` OR (`is_frameshift(align_dna.ref_aln)` and not already set). A `'-'`-run length not divisible by 3. Sets `frameshift = True`. | `93-98` |
| 4 | `start_lost` | First codon differs from ref AND `!= 'ATG'` AND first protein residue differs from ref AND `!= 'M'` (full 4-part `and`). | `99-103` |
| 5 | `inframe_insertion` | `'-' in align_dna.ref_aln` and not `frameshift`. | `104-105` |
| 6 | `inframe_deletion` | `'-' in align_dna.query_aln` and not `frameshift`. | `106-107` |
| 7 | `nonsynonymous` | `len(peps) == 2 and peps[1] == ""` (valid protein ending in `*`, alignment not 100%) **and** no other label yet accumulated. | `108-112` |
| 8 | `stop_missing` | `len(peps) == 1` (no stop codon). | `113-115` |
| 9 | `stop_codon_gain` | else branch (premature internal stop): an internal `*` before the last residue. | `116-117` |

`has_stop_codon` (`variants.py:1-15`): returns True if any target-align position is `*` while the corresponding ref-align position is not `*`. `is_frameshift` (`18-42`): scans for any maximal run of `'-'` whose length `% 3 != 0`.

---

### Appendix D — File-tree map of `lifton/`

Top-level package (`lifton/`), with one-line roles. Subtrees marked **[vendored]** should be treated as frozen dependencies during refactors of the rest.

| Path | Role |
|---|---|
| `lifton/lifton.py` | CLI entry (`main`), `parse_args`, and `run_all_lifton_steps` — the 11-step sequential pipeline orchestrator. |
| `lifton/lifton_class.py` | Object model: `Lifton_GENE`/`Lifton_TRANS`/`Lifton_EXON`/`Lifton_CDS`, `Lifton_Alignment`, `Lifton_Status`, `Lifton_ORF`; ORF search, CDS-boundary patching, `write_entry` chain. |
| `lifton/lifton_utils.py` | Shared helpers: feature partitioning, id mapping, overlap-ratio checks, alignment wrappers, status writers. (In a package-level cycle with `lifton_class`.) |
| `lifton/align.py` | parasail alignment kernels (protein + DNA), DNA sanitisation, CDS↔protein boundary mapping, translation dispatch. |
| `lifton/get_id_fraction.py` | `get_AA_id_fraction` / `get_DNA_id_fraction`: identity = matches/length over traceback strings. |
| `lifton/variants.py` | `find_variants` 9-type mutation classifier + `has_stop_codon` / `is_frameshift` helpers. |
| `lifton/protein_maximization.py` | CDS chaining algorithm: `chaining_algorithm`, `push_cds_idx`, `create_lifton_entries` — merge Liftoff/miniprot CDS by per-chunk identity. |
| `lifton/annotation.py` | `Annotation` reference-DB abstraction (gffutils-SQLite default, gffbase-DuckDB opt-in). |
| `lifton/gffbase_adapter.py` | Shim translating LiftOn's `Annotation` shape ↔ gffbase API. |
| `lifton/annotation_validator.py` | Reference-annotation structural validation helpers. |
| `lifton/extract_sequence.py` | Stream-extract transcript/protein FASTAs to `intermediate_files/`; re-open via pyfaidx. |
| `lifton/intervals.py` | Per-chromosome `IntervalTree` construction/seeding for overlap queries. |
| `lifton/run_liftoff.py` | Drive vendored Liftoff in-process; `process_liftoff` per-gene walk; recursion-limit raise/restore; cycle detection. |
| `lifton/run_miniprot.py` | Run miniprot subprocess (`--gff-only`); bounded stdout drain; failure-mode handling; miniprot-only Step 8 processing. |
| `lifton/run_evaluation.py` | Evaluation-mode runner (`-E`/`-EL`). |
| `lifton/parallel.py` | `parallel_step7`, `StepContext`, `LocusResult`, heap-backed `OrderedWriter` with spill-to-disk; thread-pool dispatch. |
| `lifton/locus_pipeline.py` | `materialise_locus`, `MaterialisedLocus`, the read-only worker DB proxies (`_RefDbProxy`/`_LFeatureDbProxy`/`_MFeatureDbProxy`). |
| `lifton/stats.py` | `print_report` — final run statistics. |
| `lifton/logger.py` | Logging helpers (`log`, `log_error`). |
| `lifton/exceptions.py` | Custom exceptions (`LiftOnAlignmentError`, `LiftOnInputError`, …). |
| `lifton/gff3_validator.py` | Output-side GFF3 validator + `_main` (`gff3-validate` console script). |
| `lifton/__init__.py` | Package init, `__version__`. |
| `lifton/io/gff3_writer.py` | Canonical GFF3 serialisation: directives, attribute ordering/encoding, `format_feature`. |
| `lifton/io/gff3_validator.py` | Input-side strict-GFF (NCBI) validator used by `--strict-gff`. |
| `lifton/io/ncbi_gff3_spec.py` | NCBI GFF3 spec rules/tables for the strict validator. |
| `lifton/native_bindings/minimap_facade.py` | `mappy`-backed in-process minimap2 facade (`--native`). |
| `lifton/native_bindings/miniprot_facade.py` | `pyminiprot`-shaped miniprot facade (subprocess today; PyO3-swappable). |
| `lifton/native_bindings/types.py` | Shared dataclasses/types for the native bindings. |
| `lifton/liftoff/` | **[vendored]** Complete in-tree fork of upstream Liftoff (~3 KLOC, 18 files). Invoked as a library; shells out to minimap2 unless `--native` routes through `native_align.py`. |
| `lifton/gffbase/` | **[vendored]** First-party DuckDB-backed gffutils successor (~3.6 KLOC + Rust `_native*.so`; pure-Python fallback under `_pyfallback/`). Drop-in FeatureDB/Feature/DataIterator/GFFWriter. |
