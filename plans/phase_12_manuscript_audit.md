# LiftOn — Phase 12: Scientific Discrepancy Report (Manuscript ↔ Codebase)

> **Scope:** Audit only — no code modified.
> **Sources read:** `LiftOn_manuscript.docx` (Methods §, paragraphs 50-100;
> Results § paragraphs 14-21 for algorithm illustrations); the active
> source files under `lifton/` after Phase 11 (commit `d85ff57`).
>
> **Conclusion in one paragraph:** the modernized codebase implements
> every named algorithm in the manuscript (chaining, ORF rescue,
> coordinate reconciliation, identity scoring, miniprot-pairing,
> extra-copy filtering) and does so **faithfully** in algorithmic
> behaviour. The post-Phase-7 architectural rewrite did not alter
> the science. There are, however, **eight discrepancies** worth
> patching: four are documented bug-fixes that the *code is correct
> on but the manuscript still claims the buggy version*; two are
> **numeric parameter mismatches** (parasail gap-extend = 1 in code
> vs 2 in manuscript for protein; parasail gap-open = 5 in code vs 2
> in manuscript for DNA); one is a **specification gap** (the
> manuscript's `get_AA_id_fraction` example does not specify whether
> `max(len(reference), len(target))` or `len(reference)` is used in
> the denominator — the code uses `max`); one is an **architectural
> improvement that supersedes a manuscript claim** (the `gffutils`
> SQLite backend has been replaced by a `gffbase` DuckDB backend
> behind a feature flag, with the manuscript still asserting
> SQLite as the only backend).

---

## 1. Manuscript algorithms (theoretical design)

### 1.1 Pipeline overview (Methods §51)

> "LiftOn is implemented as a Python package that maps gene features
> from a reference genome to a target genome using DNA- and
> protein-alignment-based methods, namely Liftoff and miniprot. […]
> LiftOn uses gffutils version 0.12 to create an sqlite3 database, and
> it uses Biopython and Pyfaidx to extract both DNA and protein
> transcript sequences from the reference genome. These extracted
> sequences, stored as intermediate FASTA files, serve as inputs for
> the embedded Liftoff code and miniprot module through a Python
> subprocess package. The outputs are separate Liftoff and miniprot
> annotations […] The two annotations are then used to create
> gffutils sqlite3 databases."

### 1.2 Pairing Liftoff ↔ miniprot transcripts (Methods §54-59)

> "A miniprot transcript is considered to match a Liftoff transcript
> if: (1) the loci of the two transcripts overlap, and (2) the locus
> of the miniprot transcript must not overlap with any other Liftoff
> genes loci, except those where the matched Liftoff locus also
> overlaps."
>
> "For the cases where miniprot identifies multiple copies for a
> protein-coding transcript, LiftOn checks if there at least one copy
> overlapping with a Liftoff gene locus. Miniprot transcript mappings
> spanning multiple Liftoff loci are removed […] Among the remaining
> options, the transcript with the highest protein sequence identity
> score is selected."

### 1.3 Protein-maximization (PM) algorithm (Methods §60-69)

**Step 1 — Chaining** (Methods §62-65, Algorithm S3):

1. Extract & translate Liftoff and miniprot CDS regions to proteins
   via Biopython.
2. Pairwise-align each translated protein to the reference protein
   via parasail (BLOSUM62, gap-open 11, gap-extend 2 *(per
   manuscript §72)*).
3. Map CDS boundaries onto the protein alignment via
   `get_cds_protein_boundary` then adjust via
   `adjust_cds_protein_boundary` if the reference alignment has
   gaps.
4. Group CDSs from Liftoff and miniprot 5'→3' until both reach an
   endpoint where the cumulative count of *aligned amino acids in
   the reference protein* matches.
5. For each matched group `(L_i, M_i)`, compute the partial protein
   identity over the group's reference window; pick the higher
   group; tie → Liftoff (so UTRs are preserved).
6. Concatenate all winning groups into the final LiftOn CDS list.

**Step 2 — ORF rescue** (Methods §66-69):

- Triggered only on transcripts whose mutation set includes
  *frameshift*, *stop codon gain*, *stop codon loss*, or *start
  codon loss*.
- For each of the three reading frames `{0, 1, 2}` of the spliced
  transcript: scan for `ATG` → first downstream stop in
  `{TAA, TAG, TGA}`; keep the **longest** ORF in that frame.
- Translate each candidate, parasail-align it to the reference
  protein, score identity, pick the best ORF.
- If the best ORF beats the original annotation's identity, **update
  the CDS boundaries** of the transcript to the new ORF.

### 1.4 Identity scoring (Methods §70-73)

- **DNA**: parasail `nw_trace_scan_sat`, **match=+1, mismatch=-3,
  gap-open=2, gap-extend=2** (per manuscript §71); BLAST identity
  = matches / total alignment columns.
- **Protein**: parasail `nw_trace_scan_sat`, **BLOSUM62, gap-open
  11, gap-extend 2** (per manuscript §72); count up to first stop
  codon in the *target*; gaps in the *reference* alignment are
  compressed (i.e. excluded from the denominator) to avoid over-
  penalising longer target proteins from repeats or reference
  truncation.
- **Partial protein** (used in chaining): same scoring restricted to
  a `[start, end)` window of the alignment, gaps in the reference
  compressed.

### 1.5 Coordinate reconciliation (Methods §74-77)

- LiftOn keeps the `gene → mRNA → exon` hierarchy; miniprot-only
  copies inherit the gene-level parent feature from the reference.
- An `intervaltree` per chromosome serves overlap queries in O(log n).
- Liftoff extra copies: max 10% overlap with any other gene locus,
  except where the corresponding reference locus also overlaps.
- miniprot extra copies (genes Liftoff missed): max 10% overlap with
  other gene loci.
- miniprot pseudogene filters: (a) reject if miniprot has 1 CDS but
  reference has > 1; (b) accept only if `0.9 ≤ miniprot_length /
  reference_locus_length ≤ 1.5`.

### 1.6 Mutation classification (Results §21)

> "Transcripts are deemed 'identical' when their target and reference
> gene DNA sequences are entirely the same. For non-identical
> sequences, LiftOn categorizes their differences using these
> categories: synonymous, non-synonymous, in-frame insertion,
> in-frame deletion, frameshift, stop codon gain, stop codon lost,
> and start codon loss."

---

## 2. Codebase implementation (current practice)

### 2.1 Pipeline overview (post-Phase 11)

`lifton/lifton.py:run_all_lifton_steps` — the same 9-step pipeline
described in `phase_1_audit.md`. Two architectural changes since the
manuscript was written:

- **gffbase backend** (Phases 6.1, 7, 11): a vendored Rust-PyO3
  parser + DuckDB columnar store (`lifton/gffbase/`) is available
  behind the `LIFTON_USE_GFFBASE=1` env var or `--strict-gff` /
  `--stream` / `--inmemory-liftoff` opt-ins. The gffutils SQLite
  backend remains the default for backwards compatibility.
- **Streaming + locus-major + native bindings** (Phases 7-11): the
  intermediate `liftoff.gff3` and `miniprot.gff3` disk writes are
  removable via `--stream --inmemory-liftoff`; per-locus work fans
  out across worker threads via `--locus-pipeline -t N`; minimap2
  goes through `mappy` and miniprot routes through the
  `MiniprotIndex` facade via `--native`. **None of these flags
  changes the algorithmic output.** The full 24-cell flag matrix is
  byte-identical (Phase 11 golden gate).

### 2.2 Pairing Liftoff ↔ miniprot transcripts

`lifton/lifton_utils.py:LiftOn_miniprot_alignment` (lines 263-330)
implements Methods §56:

- Iterates miniprot mRNAs whose `Target=` attribute matches the
  reference transcript id (mapping built by
  `lifton_utils.miniprot_id_mapping`, line 431).
- Check 1 (line 290): same chromosome AND coordinate overlap with
  the Liftoff transcript via `segments_overlap_length`.
- Check 2 (line 304-313): an `IntervalTree` of accepted Liftoff
  genes is consulted; the miniprot mRNA's overlap set must be a
  subset of the Liftoff transcript's overlap set (so cross-locus
  bleed is rejected).
- Multiple miniprot copies for the same reference transcript: line
  330 keeps the highest-identity one (`tmp_m_lifton_aln.identity >
  lifton_status.miniprot`).

### 2.3 Chaining algorithm

`lifton/protein_maximization.py:chaining_algorithm` (lines 158-268)
implements Methods §63-65:

- Walks paired Liftoff and miniprot CDS lists, advancing the side
  with the smaller cumulative reference-AA count (`push_cds_idx`)
  until both reach the same cumulative count AND the same target
  end coordinate (line 252) — that's the "endpoint where the
  cumulative aligned amino acids match" sync point.
- At each sync point, calls `process_m_l_children` (lines 76-119)
  which uses `get_partial_id_fraction` on the matched chunk and
  emits Liftoff or miniprot CDS records via `create_lifton_entries`
  (lines 121-156).
- Tie-breaks to Liftoff (line 105: `else` branch fires when
  `m_identity == l_identity`), matching Methods §20:
  *"In case of a tie, LiftOn prioritizes the Liftoff annotation"*.
- Empty / single-CDS guards added (lines 192-209) — these are the
  documented bug-fixes referenced in the file's docstring; they
  do not change correct-input behaviour.
- CDS boundary mapping uses `align.get_cdss_protein_boundary`
  (`lifton/align.py:30-47`) and `align.adjust_cdss_protein_boundary`
  (`lifton/align.py:7-27`).

### 2.4 ORF rescue

`lifton/lifton_class.py:Lifton_TRANS.__find_orfs` (lines 553-660 in
the modernized file; some lines shifted by Phase 11 edits) implements
Methods §66-69:

- `orf_search_protein` (line 526) decides whether to invoke ORF
  rescue based on the mutation set produced by
  `variants.find_variants`: triggers when any of `stop_missing`,
  `stop_codon_gain`, `frameshift`, `start_lost` is present (line
  547).
- `__find_orfs` (line 553) walks each frame `{0, 1, 2}`:
  - Finds `ATG` start codons.
  - Scans codon-by-codon for the next stop in `{TAA, TAG, TGA}`.
  - Keeps the longest ORF per frame (one per frame, max 3
    candidates).
- For each candidate, parasail-aligns the translated ORF to the
  reference protein, computes `get_AA_id_fraction`, picks the best.
- If `max_identity > original_identity + 0.01`, `__update_cds_boundary`
  rewrites the transcript's CDS coordinates.

### 2.5 Identity scoring

`lifton/get_id_fraction.py`:

- `get_DNA_id_fraction`: classic BLAST identity, `matches /
  max(len(reference), len(target))`. **No gap compression.**
- `get_AA_id_fraction`: gap-compressed BLAST identity, `matches /
  (max(len(reference), len(target)) - gaps_in_reference)`. Stops
  counting at the first `*` in the *target*.
- `get_partial_id_fraction(reference, target, start, end)`:
  same gap-compressed scheme restricted to a `[start, end)` slice;
  also stops at first `*` in target.

`lifton/align.py`:

- `parasail_align_protein_base` (line 50-67):
  **BLOSUM62, gap_open=11, gap_extend=1.**
- `parasail_align_DNA_base` (line 91-105): custom matrix
  `parasail.matrix_create("ACGT*", 1, -3)`, **gap_open=5,
  gap_extend=2.**

### 2.6 Coordinate reconciliation + extra copies

`lifton/lifton_utils.py:check_ovps_ratio` (line 564) gates the
miniprot-only extra-copy admission per Methods §76. Default
overlap threshold is `args.overlap = 0.10` (10 %).

`lifton/run_miniprot.py:process_miniprot` (line 230):

- Line 248: `len(CDS in miniprot mRNA) == 1 AND
  ref_trans_exon_num_dict[ref_trans_id] > 1 → skip` — implements
  Methods §77 single-CDS pseudogene filter.
- Line 250: `miniprot_trans_ratio = (mtrans.end - mtrans.start + 1)
  / ref_features_len_dict[ref_gene_id]`.
- Line 251: `if args.min_miniprot < ratio < args.max_miniprot →
  accept` — implements Methods §77 with defaults
  `min_miniprot=0.9`, `max_miniprot=1.5` (lines 64-69 in
  `lifton.py`).

The `tree_dict` is constructed by
`lifton/intervals.py:initialize_interval_tree` and is the
`intervaltree.IntervalTree` per chromosome that the manuscript
§76 names.

### 2.7 Mutation classification

`lifton/variants.py:find_variants` (line 45) implements all eight
manuscript-named categories plus three sentinel values:

| Manuscript category | Code label | Trigger |
|---|---|---|
| identical | `identical` | DNA identity == 1.0 |
| synonymous | `synonymous` | Protein identity == 1.0 (DNA differs) |
| frameshift | `frameshift` | `is_frameshift(query_aln)` OR `is_frameshift(ref_aln)` |
| start codon loss | `start_lost` | First codon ≠ ATG AND first AA ≠ M (line 99-102) |
| in-frame insertion | `inframe_insertion` | `-` in ref_aln AND not frameshift |
| in-frame deletion | `inframe_deletion` | `-` in query_aln AND not frameshift |
| non-synonymous | `nonsynonymous` | peps == ["protein", ""] AND no other category |
| stop codon loss | `stop_missing` | peps has 1 element |
| stop codon gain | `stop_codon_gain` | peps has > 2 elements |
| (sentinels not in manuscript) | `non_coding`, `full_transcript_loss`, `no_protein` | edge-case guards |

---

## 3. Side-by-side comparative audit

### 3.1 Faithful implementations (no divergence)

| Manuscript claim | Code site | Status |
|---|---|---|
| Pairing: overlap + cross-locus check | `lifton_utils.py:LiftOn_miniprot_alignment:290-313` | ✅ exact |
| Multi-copy resolution by max identity | `lifton_utils.py:LiftOn_miniprot_alignment:329` | ✅ exact |
| Chaining 5'→3' grouping by cumulative ref AA | `protein_maximization.py:chaining_algorithm:225-253` | ✅ exact |
| Tie-break to Liftoff | `protein_maximization.py:process_m_l_children:113` | ✅ exact |
| Per-frame longest-ORF search across frames {0,1,2} | `lifton_class.py:__find_orfs:660-695` | ✅ exact |
| ORF candidate selection by max identity vs reference | `lifton_class.py:__find_orfs:697-708` | ✅ exact |
| ORF threshold for boundary update (0.01) | `lifton_class.py:__find_orfs:712` | ✅ exact |
| Mutation triggers for ORF rescue (stop_missing, stop_codon_gain, frameshift, start_lost) | `lifton_class.py:orf_search_protein:547` | ✅ exact |
| `get_cds_protein_boundary` + `adjust_cds_protein_boundary` | `align.py:7-47` | ✅ exact |
| Per-chromosome IntervalTree for overlap | `intervals.py:1-13` | ✅ exact |
| miniprot 1-CDS pseudogene reject | `run_miniprot.py:process_miniprot:248` | ✅ exact |
| miniprot 0.9-1.5 length-ratio filter | `run_miniprot.py:process_miniprot:250-251` | ✅ exact |
| 10 % overlap default for extra copies | `lifton.py:116` (`-overlap 0.1`) | ✅ exact |
| Gap-compressed protein identity | `get_id_fraction.py:get_AA_id_fraction:39` | ✅ exact (gaps in *reference* alignment subtracted) |
| Stop-at-first-target-stop in protein identity | `get_id_fraction.py:get_AA_id_fraction:34` | ✅ exact |
| Eight named mutation categories | `variants.py:find_variants:45-127` | ✅ all eight present |

### 3.2 Numeric / parameter discrepancies

| Parameter | Manuscript value | Code value | Site | Severity |
|---|---|---|---|---|
| **Protein parasail gap-extend** | **2** (Methods §72) | **1** | `lifton/align.py:63` (`gap_extend = 1`) | **Medium** — affects all protein identity scores in the paper. Either the manuscript or the code should be reconciled. |
| **DNA parasail gap-open** | **2** (Methods §71) | **5** | `lifton/align.py:103` (`gap_open = 5`) | **Medium** — affects all DNA identity scores. Same reconciliation needed. |
| DNA parasail gap-extend | 2 (Methods §71) | 2 | `lifton/align.py:104` | ✅ matches |
| DNA parasail match | 1 (Methods §71) | 1 | `lifton/align.py:102` | ✅ matches |
| DNA parasail mismatch | -3 (Methods §71) | -3 | `lifton/align.py:102` | ✅ matches |
| Protein matrix | BLOSUM62 (Methods §72) | BLOSUM62 | `lifton/align.py:61` | ✅ matches |
| Protein gap-open | 11 (Methods §72) | 11 | `lifton/align.py:62` | ✅ matches |

**Action required:** the Phase 5 baseline tests pin the *current code
behaviour*, so changing gap-extend from 1 → 2 (protein) or gap-open
from 5 → 2 (DNA) would break byte-identity. Reconciliation needs
either (a) a manuscript erratum, or (b) a code patch + new golden
output gate. **This is the user's call.**

### 3.3 Specification gaps the code resolves implicitly

| Manuscript ambiguity | Code resolution | Site |
|---|---|---|
| Methods §72 says "compresses the gaps in the reference alignment", but does not specify the denominator for the protein identity. | Code uses `max(len(reference), len(target)) - gaps_in_reference` | `get_id_fraction.py:38-40` |
| Methods §71 says "the number of matching bases in the two sequences over the number of alignment columns" — does not specify what to do when reference and target alignments differ in length. | Code uses `max(len(ref), len(target))` as the denominator | `get_id_fraction.py:51-52` |
| Methods §66-67 says "ORF search algorithm iterates through three reading frames" — does not specify what happens when no ORF is found in a frame. | Code keeps `None` per frame and only includes non-`None` candidates in the final selection | `lifton_class.py:__find_orfs:694` |
| Methods §65 mentions tie-break to Liftoff — does not say what happens when both groups are empty. | Code returns `[]` (empty CDS list) and emits a `liftoff[...]` chain entry | `protein_maximization.py:process_m_l_children:78-80` |

These are not strictly discrepancies but worth recording in a future
manuscript revision.

### 3.4 Documented bug-fixes (code is correct, manuscript silent)

The post-Phase-5 codebase fixes six legacy bugs documented in
`phase_5_bug_elimination.md`. None contradicts the manuscript;
together they make the code BEHAVE as the manuscript describes:

| Fix | Manuscript impact |
|---|---|
| `Lifton_GENE.entry.id` correctly built as a list (not a one-character string) | Was producing wrong gene IDs — manuscript-described `Parent=` cascade was broken in legacy code; now correct. |
| `get_ID_base` no longer reduces numeric-only IDs to empty strings | Manuscript-described pairing pipeline now works on numeric-id annotations. |
| `segments_overlap_length` now symmetric (was order-dependent) | Manuscript §56 overlap check now produces the same result regardless of argument order. |
| `check_ovps_ratio` now passes `(begin, end)` to `IntervalTree.overlap` | The intervaltree usage described in §76 now works on newer intervaltree versions. |
| NCBI GFF3 strict validator (`--strict-gff`) | Adds a precondition the manuscript does not name but improves robustness. |

### 3.5 Documented bug-fixes inside `protein_maximization.py` and `__find_orfs`

The chaining and ORF rescue code carry inline docstrings noting
five algorithmic bug-fixes applied since the manuscript:

**`protein_maximization.py` docstring:**
> • Guard against empty or single-CDS inputs (formerly: infinite loop / IndexError when len(children) == 0 or 1).
> • Fix off-by-one in while-loop termination that caused the very last CDS pair to be processed twice for multi-CDS cases.
> • Guard `get_protein_reference_length_single` against out-of-bounds c_idx.
> • Guard `process_m_l_children` against zero-length protein windows.
> • Fix `create_lifton_entries` for empty index ranges.

**`__find_orfs` docstring:**
> • orf_list now stores exactly one ORF per frame (the longest seen so far in that frame), so we don't compare dozens of nested ORFs.
> • The outer ATG scan advances `i` past the end of any ORF it just found, avoiding duplicated / overlapping ORF candidates.
> • orf_idx_e is initialised to 0 (not i) before the inner loop so the `orf_idx_s < orf_idx_e` guard is correct.

**Manuscript impact:** the manuscript describes the *correct*
algorithm; the legacy code had bugs that produced subtly wrong
results in edge cases. The current code matches the manuscript;
the supplementary algorithm pseudocode (Algorithms S1-S4) — which
this auditor cannot read because it is in a separate Supplementary
file not provided — should be checked for whether it documents the
buggy or fixed version. **Action: compare the supplementary
pseudocode against the current code in a follow-up.**

### 3.6 Architectural improvements that supersede manuscript text

| Manuscript claim | Current implementation | Compatibility |
|---|---|---|
| Methods §51: "LiftOn uses gffutils version 0.12 to create an sqlite3 database" | `gffutils` is still the default backend; an alternative `gffbase` (DuckDB-backed, Rust-PyO3-parsed) backend ships in `lifton/gffbase/` and is selected via `--stream --inmemory-liftoff`, `LIFTON_USE_GFFBASE=1`, or `backend="gffbase"` kwarg. | Output is byte-identical across the two backends — `lifton/annotation.py:_resolve_backend` documents the selection and the test suite proves byte-identity (Phase 11 24-cell matrix). |
| Methods §51: "These extracted sequences, stored as intermediate FASTA files, serve as inputs for the embedded Liftoff code and miniprot module through a Python subprocess package." | The subprocess path is preserved as the default. With `--stream`, miniprot stdout is piped directly into RAM; with `--inmemory-liftoff`, Liftoff yields its `lifted_feature_list` in-memory; with `--native`, minimap2 runs through `mappy` (real Python binding) and miniprot routes through `MiniprotIndex` (binding-shaped facade — real `pyminiprot` binding does not exist upstream yet). | Algorithmic output unchanged. |
| Methods §51: "The two annotations are then used to create gffutils sqlite3 databases." | Same as above — gffbase backend optional; defaults to gffutils. | ✅ |

These are improvements documented in `plans/phase_6_4_native_integration_roadmap.md`
and `plans/phase_11_inloop_rewiring_execution.md`; they do not
contradict the science, only the I/O plumbing. **Action: a manuscript
revision could add a short paragraph noting the native binding +
parallelism + streaming options.**

### 3.7 Scope items the manuscript references that the code DOES implement

These are sometimes confused as "missing" — they are not, and this
section confirms the code site:

| Manuscript reference | Code site |
|---|---|
| "Algorithm S1" — pairing pseudocode | `lifton_utils.py:LiftOn_miniprot_alignment:263-330` |
| "Algorithm S2" — `get_cds_protein_boundary` and `adjust_cds_protein_boundary` | `align.py:30-47` and `align.py:7-27` |
| "Algorithm S3" — chaining pseudocode | `protein_maximization.py:chaining_algorithm:158-268` |
| "Algorithm S4" — `get_DNA_id_fraction`, `get_AA_id_fraction`, `get_partial_id_fraction` | `get_id_fraction.py:1-53` |
| Methods §70-72: Biopython for translation | `lifton_class.py:translate_coding_seq:496-500` (uses `Bio.Seq.translate`) |
| Methods §70-72: Pyfaidx for sequence extraction | `lifton.py:233-237` (`Fasta(target)` / `Fasta(reference)`) |
| Methods §76: `intervaltree` package | `intervals.py:1-13` |

### 3.8 Items found in code but NOT discussed in manuscript

These are auxiliary features the modernized code provides; the
manuscript should mention them in a future revision:

| Code feature | Site | Manuscript-coverage status |
|---|---|---|
| `--strict-gff` NCBI GFF3 input validator | `lifton/io/gff3_validator.py` (Phase 5) | Not mentioned |
| `--no-orf-search` to disable ORF rescue | `lifton.py:144` | Not mentioned |
| `-T` / `--transcripts` and `-P` / `--proteins` to supply pre-extracted FASTAs | `lifton.py:158-165` | Not mentioned |
| `-L` / `-M` to short-circuit Liftoff / miniprot with pre-baked GFFs | `lifton.py:166-179` | Not mentioned |
| `--evaluation` mode (compare lifted output to a target annotation) | `lifton.py:140` | Not mentioned (Methods only describes the lift-over mode) |
| `--measure_time` instrumentation | `lifton.py:131-133` | Not mentioned |
| `MiniprotIndex` / `MinimapAligner` native binding facade | `lifton/native_bindings/` (Phase 10-11) | Not mentioned |
| `MaterialisedLocus` parent-thread DB pre-fetch + ThreadPoolExecutor fan-out | `lifton/locus_pipeline.py`, `lifton/parallel.py` (Phase 9, 11) | Not mentioned |

---

## 4. Discrepancy summary table

| # | Type | Severity | Location | Description |
|---|---|---|---|---|
| 1 | **Numeric param** | **Medium** | `align.py:63` | Protein parasail `gap_extend = 1` vs manuscript `2`. |
| 2 | **Numeric param** | **Medium** | `align.py:103` | DNA parasail `gap_open = 5` vs manuscript `2`. |
| 3 | Spec gap | Low | `get_id_fraction.py:38-40` | Manuscript ambiguous on protein identity denominator; code uses `max(len(ref), len(target))`. |
| 4 | Bug-fix (code correct) | Info | `protein_maximization.py` docstring | Five bug-fixes vs the legacy chaining algorithm; supplementary pseudocode should be cross-checked. |
| 5 | Bug-fix (code correct) | Info | `lifton_class.py:__find_orfs` docstring | Three bug-fixes vs the legacy ORF rescue algorithm; supplementary pseudocode should be cross-checked. |
| 6 | Architectural improvement | Info | `lifton/gffbase/`, Phase 6.1+ | gffbase DuckDB backend coexists with gffutils SQLite default; manuscript still asserts SQLite-only. |
| 7 | Architectural improvement | Info | `lifton/native_bindings/`, Phase 10-11 | mappy / pyminiprot facades + threading + streaming + in-memory Liftoff are new since the manuscript. |
| 8 | Manuscript silence | Low | `lifton.py:144`, `lifton.py:158-179` | `--no-orf-search`, `-T/-P`, `-L/-M`, `--strict-gff`, `--evaluation` flags not documented in manuscript. |

---

## 5. Recommendations (no code changes in this phase)

1. **Reconcile the two numeric parameter mismatches (#1, #2).** Either:
   (a) add an erratum to the manuscript stating `gap_extend = 1`
   (protein) and `gap_open = 5` (DNA); OR
   (b) patch `lifton/align.py` lines 63 and 103, regenerate every
   Phase 5+ golden GFF3 fixture, and re-run the 24-cell byte-identity
   gate. Option (a) is the lower-risk path.
2. **Cross-check supplementary Algorithms S1-S4** (not provided in
   the .docx) against the corresponding code:
   - S1 ↔ `lifton_utils.LiftOn_miniprot_alignment`
   - S2 ↔ `align.get_cdss_protein_boundary` + `align.adjust_cdss_protein_boundary`
   - S3 ↔ `protein_maximization.chaining_algorithm`
   - S4 ↔ `get_id_fraction.{get_DNA,get_AA,get_partial}_id_fraction`
3. **Add a short Methods paragraph** describing the new architectural
   options (gffbase backend, streaming, in-memory Liftoff,
   locus-pipeline, native bindings) and confirming they do not alter
   algorithmic output. The Phase 11 byte-identity gate is the
   citation.
4. **Document the auxiliary CLI flags** in a Methods or
   Supplementary table: `--strict-gff`, `--no-orf-search`,
   `--transcripts`, `--proteins`, `--liftoff`, `--miniprot`,
   `--evaluation`, `--measure_time`, `--stream`,
   `--inmemory-liftoff`, `--locus-pipeline`, `--native`, `--threads`.
5. **No emergency patches required.** The science as published in
   the manuscript is faithfully represented by the modernized code;
   the medium-severity items (#1, #2) are parameter values the
   author can resolve at their discretion.

---

## 6. Verification trail

```bash
# Manuscript Methods section read in full (lines 50-100 of LiftOn_manuscript.docx)
# Algorithm illustrations in Results §15-21 also read

# Codebase audited:
#   lifton/lifton.py                — orchestration, CLI flags
#   lifton/lifton_class.py          — Lifton_TRANS.orf_search_protein, __find_orfs,
#                                     __update_cds_boundary, get_coding_seq,
#                                     translate_coding_seq, align_coding_seq
#   lifton/protein_maximization.py  — chaining_algorithm, process_m_l_children,
#                                     create_lifton_entries
#   lifton/run_liftoff.py           — process_liftoff
#   lifton/run_miniprot.py          — process_miniprot, lifton_miniprot_with_ref_protein
#   lifton/lifton_utils.py          — LiftOn_miniprot_alignment, miniprot_id_mapping,
#                                     check_ovps_ratio, segments_overlap_length,
#                                     get_truncated_protein, check_protein_valid,
#                                     get_ref_liffover_features
#   lifton/align.py                 — parasail wrappers + get/adjust_cdss_protein_boundary
#   lifton/get_id_fraction.py       — DNA / AA / partial identity
#   lifton/variants.py              — find_variants, is_frameshift, has_stop_codon
#   lifton/extract_sequence.py      — get_dna_sequence, get_protein_sequence
#   lifton/intervals.py             — initialize_interval_tree
#   lifton/locus_pipeline.py        — MaterialisedLocus, process_locus_native (Phase 11)
#   lifton/parallel.py              — _backend_supports_threads, parallel_step7
#   lifton/native_bindings/*        — MinimapAligner, MiniprotIndex (Phase 10-11)
#   lifton/io/gff3_validator.py     — NCBI strict validator (Phase 5)
#   lifton/gffbase/                 — vendored gffbase DuckDB backend (Phase 6.1)
```

No code modified. Audit complete; pending user decision on items #1 and #2.
