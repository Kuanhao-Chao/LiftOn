## 3. Core Algorithms & Logic

This section specifies the lowest layer of LiftOn's per-transcript machinery: how a spliced coding/transcript nucleotide sequence is assembled from per-exon CDS pieces and translated to protein (§3.1), how protein and DNA sequences are aligned with parasail and how the resulting alignment is post-processed into a `Lifton_Alignment` with per-CDS protein boundaries (§3.2), and the exact arithmetic of the three sequence-identity functions that score those alignments (§3.3). Everything here is byte-identity-relevant: the alignment kernel, the matrix parameters, the boundary math, and the identity denominators are all pinned by the 24-cell matrix test, so the constants and formulas below must be reproduced exactly.

### 3.1 Sequence extraction, assembly, and translation

The transcript object `Lifton_TRANS` (`lifton/lifton_class.py:267`) holds an ordered list `self.exons` of `Lifton_EXON` objects, each of which may carry exactly one `Lifton_CDS` child at `exon.cds` (or `None`). Each `Lifton_EXON` and `Lifton_CDS` wraps a gffutils/gffbase feature `.entry` whose `.sequence(fai)` method extracts the genomic substring (already reverse-complemented by the underlying library when `.strand == '-'`). Three methods turn these per-exon pieces into the strings that feed alignment.

#### 3.1.1 `get_coding_seq(self, fai)` — `lifton/lifton_class.py:545`

Builds the **spliced CDS nucleotide sequence** (concatenation of CDS pieces only, exons without a CDS contribute nothing), plus a parallel list of per-CDS lengths and a deep-copied list of CDS feature entries.

Signature: `get_coding_seq(self, fai) -> (coding_seq: str, cds_children: list, cdss_lens: list[int])`

Algorithm (numbered):

1. Initialise `coding_seq = ""`, `cdss_lens = []`, `cds_children = []`.
2. Iterate `exon` over `self.exons` **in stored list order**.
3. If `exon.cds is None`, skip this exon entirely (no contribution to any of the three outputs).
4. Otherwise append `copy.deepcopy(exon.cds.entry)` to `cds_children`. (Deep copy — the returned entries are independent of the live feature objects.)
5. Extract the CDS nucleotide piece `p_seq = exon.cds.entry.sequence(fai)`. This is already strand-correct (reverse-complemented for `-` strand) by the FASTA library.
6. **Strand-dependent assembly** (the key gotcha):
   - If `exon.cds.entry.strand == '-'`: **prepend** — `coding_seq = p_seq + coding_seq`, and **insert the length at the front** — `cdss_lens.insert(0, exon.cds.entry.end - exon.cds.entry.start + 1)`.
   - If `exon.cds.entry.strand == '+'`: **append** — `coding_seq = coding_seq + p_seq`, and append the length — `cdss_lens.append(exon.cds.entry.end - exon.cds.entry.start + 1)`.
   - Any other strand value (`.`, `?`, empty) matches **neither** branch: the piece is silently dropped from `coding_seq` and `cdss_lens`, even though its entry was already appended to `cds_children` in step 4. (Subtle asymmetry — `cds_children` and `cdss_lens` can differ in length for malformed strand values.)
7. Return `(coding_seq, cds_children, cdss_lens)`.

**CDS length formula:** `end - start + 1` (1-based inclusive genomic coordinates → nucleotide count). For a `-`-strand transcript the `.insert(0, …)` reverses iteration order so `cdss_lens` ends up in **5'→3' transcript order** matching the assembled `coding_seq`. For `+`-strand the natural append order is already 5'→3'.

Gotcha: `self.exons` is assumed to be in genomic 5'→3' order for `+` and the *reverse* genomic order is corrected by the prepend/insert(0) trick for `-`; the stored list iteration order is therefore load-bearing. `cdss_lens` produced here is the exact input to `get_cdss_protein_boundary` (§3.2.4), so any reordering changes the per-CDS protein boundary map.

#### 3.1.2 `get_coding_trans_seq(self, fai)` — `lifton/lifton_class.py:562`

Builds **both** the full spliced transcript sequence (all exons, including UTR) and the spliced CDS sequence, and as a side effect **rewrites the CDS `frame` attribute** of each CDS entry.

Signature: `get_coding_trans_seq(self, fai) -> (coding_seq: str, trans_seq: str)` (note: returns `(coding_seq, trans_seq)`, in that order).

Algorithm:

1. Initialise `trans_seq = ""`, `coding_seq = ""`, `accum_cds_length = 0`.
2. Set `lcl_exons = self.exons`. **Reverse for negative strand:** if `len(self.exons) > 0 and self.exons[0].entry.strand == '-'`, set `lcl_exons = self.exons[::-1]` (a reversed *copy*, leaving `self.exons` untouched). This makes iteration proceed 5'→3' along the transcript for `-`-strand genes.
3. For each `exon` in `lcl_exons`, in order:
   1. `p_trans_seq = exon.entry.sequence(fai)` — the exon nucleotide piece (strand-correct).
   2. `p_trans_seq = Seq(p_trans_seq).upper()` — wrap in a BioPython `Seq` and uppercase.
   3. `trans_seq = trans_seq + p_trans_seq` — concatenate (BioPython `Seq` concatenation; converted to `str` at the end).
   4. If `exon.cds is not None`:
      - **Set frame:** `exon.cds.entry.frame = str(self.__get_cds_frame(accum_cds_length))` (mutates the live CDS entry — see §3.1.4).
      - **Accumulate CDS length:** `accum_cds_length += (exon.cds.entry.end - exon.cds.entry.start + 1)` (same `end-start+1` formula).
      - `p_seq = exon.cds.entry.sequence(fai)`; `coding_seq = coding_seq + p_seq`.
4. After the loop: if `trans_seq != None`, set `trans_seq = str(trans_seq).upper()`; if `coding_seq != None`, set `coding_seq = str(coding_seq).upper()`. (Both guards are always true here since both start as `""`; the upper() + str() coercion turns the `Seq` object back into a plain uppercase string.)
5. Return `(coding_seq, trans_seq)`.

Gotcha (frame mutation side effect): this method **writes** `exon.cds.entry.frame` for every CDS, in 5'→3' transcript order, computed from the running `accum_cds_length`. The CDS pieces here are assembled by plain forward append (no prepend trick) because `lcl_exons` was already reversed for `-` strand in step 2. The frames written here are what later get serialised into the output GFF3 phase column.

Gotcha (order divergence from `get_coding_seq`): `get_coding_seq` (§3.1.1) handles `-` strand by *prepending strings and insert(0)-ing lengths while iterating the stored order*, whereas `get_coding_trans_seq` handles `-` strand by *reversing the exon list then appending*. Both produce the same 5'→3' `coding_seq` for a well-formed transcript, but they reach it by different routes; only `get_coding_trans_seq` reverses and only it touches `frame`.

#### 3.1.3 `translate_coding_seq(self, coding_seq)` — `lifton/lifton_class.py:588`

Signature: `translate_coding_seq(self, coding_seq) -> str | None`.

1. Initialise `protein_seq = None`.
2. If `coding_seq != ""`: `protein_seq = str(Seq(coding_seq).translate())`. This uses BioPython `Seq.translate()` with the **standard genetic table (table 1), default behaviour**: translates codon-by-codon, internal and terminal stop codons become `*`, partial trailing codons (length not a multiple of 3) raise/are handled by BioPython's default (the input is expected to be a clean multiple of 3 from CDS assembly). No `to_stop`, no `cds=True` flag is passed — so internal stop codons are **retained as `*`** rather than truncating, and the trailing stop codon (if present) is rendered as `*`.
3. Return `protein_seq` (which is `None` only when `coding_seq == ""`).

The returned protein string is later split on `*` by callers (`align_coding_seq`, §3.2.5) to recover individual peptides (`peps`).

#### 3.1.4 `__get_cds_frame(self, accum_cds_length)` — `lifton/lifton_class.py:791`

Computes the GFF3 **phase** (called `frame` in the code) for a CDS given the number of coding nucleotides that precede it in the transcript.

Exact formula: `return (3 - accum_cds_length % 3) % 3`

This yields the standard GFF3 phase: the number of bases to skip at the start of this CDS to reach the next codon boundary. Worked values: `accum=0 → 0`; `accum=1 → 2`; `accum=2 → 1`; `accum=3 → 0`; `accum=4 → 2`; etc. The double-modulo ensures `accum % 3 == 0` maps to `0` (not `3`).

### 3.2 The alignment kernel

All pairwise alignment is global (Needleman–Wunsch) via parasail's `nw_trace_scan_sat` vectorised kernel, which returns a traceback containing three equal-length strings (`query`, `comp`, `ref`) plus an encoded CIGAR. There are two flavours — protein (BLOSUM62) and DNA (custom match/mismatch matrix) — plus a CDS→protein boundary remapping step used only on the protein path.

#### 3.2.1 Constants (must be reproduced exactly)

| Path | Function | Matrix | gap_open | gap_extend | Kernel |
|---|---|---|---|---|---|
| Protein | `parasail_align_protein_base` (`align.py:84`) | `parasail.Matrix("blosum62")` | 11 | 1 | `parasail.nw_trace_scan_sat` |
| DNA | `parasail_align_DNA_base` (`align.py:134`) | `parasail.matrix_create("ACGTN*", 1, -3)` | 5 | 2 | `parasail.nw_trace_scan_sat` |

The DNA matrix is built fresh on each call with **match score `+1`, mismatch score `-3`** over the 6-character alphabet `ACGTN*` (the `*` column models the terminal-stop / sentinel character). Argument order to both kernels is `nw_trace_scan_sat(query, target, gap_open, gap_extend, matrix)` and the code comments label the return as `(Query, Target)` where **query = the LiftOn/annotated sequence** and **target = the reference sequence** (i.e. `protein_seq` is query, `ref_protein_seq` is target).

#### 3.2.2 DNA input sanitisation — `_sanitise_for_parasail_dna(seq)` (`align.py:16`)

Module constant `_PARASAIL_DNA_ALPHABET = frozenset("ACGTN*")` (`align.py:13`). Both DNA inputs are passed through this normaliser before the kernel (`align.py:151-152`), because `parasail.matrix_create("ACGTN*", …)` has no score column for IUPAC ambiguity codes (`R Y S W K M B D H V`), gaps, spaces, or lowercase, and feeding such a byte can crash parasail's C kernel.

Algorithm:
1. If `seq` is empty/falsy, return it unchanged.
2. `upper = seq.upper()`.
3. **Fast path:** if every character of `upper` is in `_PARASAIL_DNA_ALPHABET`, return `upper` directly.
4. **Slow path:** return `"".join(ch if ch in _PARASAIL_DNA_ALPHABET else "N" for ch in upper)` — every out-of-alphabet character is coerced to `N`.

Gotcha: this is **DNA-only**. The protein path does no sanitisation — BLOSUM62 is a full 24×24 matrix that already covers `B Z X *`, so protein inputs are passed raw to the kernel. Empty protein/DNA inputs are *not* sanitised into something valid; the protein base function instead raises (next subsection).

#### 3.2.3 Base alignment functions

`parasail_align_protein_base(protein_seq, ref_protein_seq)` (`align.py:84`):
1. Build `matrix = parasail.Matrix("blosum62")`, `gap_open = 11`, `gap_extend = 1`.
2. **Empty-input refusal:** if `protein_seq == "" or ref_protein_seq == ""`, raise `LiftOnAlignmentError("parasail_align_protein_base: refusing to align empty protein sequence — caller must provide non-empty query and reference.")`. (V2.11 fix: an empty CDS-derived protein is an upstream bug, not something to silently coerce to `"*"`.)
3. Return `parasail.nw_trace_scan_sat(protein_seq, ref_protein_seq, 11, 1, matrix)`.

`parasail_align_DNA_base(trans_seq, ref_trans_seq)` (`align.py:134`):
1. Build `matrix = parasail.matrix_create("ACGTN*", 1, -3)`, `gap_open = 5`, `gap_extend = 2`.
2. `trans_seq = _sanitise_for_parasail_dna(trans_seq)`; `ref_trans_seq = _sanitise_for_parasail_dna(ref_trans_seq)`.
3. Return `parasail.nw_trace_scan_sat(trans_seq, ref_trans_seq, 5, 2, matrix)`. (No empty-input refusal on the DNA path — empties propagate to parasail.)

#### 3.2.4 `get_cdss_protein_boundary(cdss_lens)` (`align.py:64`)

Maps each CDS's nucleotide span onto **protein-coordinate (amino-acid) boundaries** using the cumulative nucleotide sum divided by 3. Input `cdss_lens` is the 5'→3'-ordered list from `get_coding_seq` (§3.1.1).

Algorithm:
1. `cdss_cumulative = [sum(cdss_lens[:i+1]) for i in range(len(cdss_lens))]` — prefix sums (running total of nucleotides through CDS `i`). E.g. `[90, 60, 30] → [90, 150, 180]`.
2. `cdss_cumulative_div = [x / 3 for x in cdss_cumulative]` — true (float) division by 3 to convert nucleotide count → amino-acid count. E.g. `[30.0, 50.0, 60.0]`.
3. Build dict `cdss_protein_boundary`: for each `idx` in `range(len(cdss_cumulative_div))`, set `cdss_protein_boundary[idx] = (start, end)` where `start = cdss_cumulative_div[idx-1] if idx > 0 else 0` and `end = cdss_cumulative_div[idx]`.
4. Return the dict. The result is half-open-ish in AA space: CDS 0 spans `[0, len0/3)`, CDS 1 spans `[len0/3, (len0+len1)/3)`, etc. Values are floats (a CDS whose nucleotide length is not divisible by 3 yields a fractional boundary, which is intentional — frame is shared across the splice junction).

#### 3.2.5 `adjust_cdss_protein_boundary(cdss_protein_aln_boundary, cigar_accum_len, length)` (`align.py:34`)

Shifts CDS protein boundaries rightward to account for a deletion (`D`) block in the protein alignment CIGAR (gaps in the *query* relative to *reference*, which add columns to the alignment coordinate). Called once per `D` block during CIGAR replay (§3.2.6).

Parameters: `cdss_protein_aln_boundary` is the dict being mutated (`{idx: (start, end)}`); `cigar_accum_len` is the alignment-column offset where the current `D` block begins; `length` is the `D` block's length.

Algorithm (V2.10 snapshot fix is the load-bearing detail):
1. **Snapshot the input dict first:** `snapshot = {i: cdss_protein_aln_boundary[i] for i in range(len(cdss_protein_aln_boundary))}`. All reads in the loop use this *pre-call* state; only the live dict is written. This prevents a single boundary's shift from compounding when the same boundary straddles multiple `D` blocks across successive calls.
2. `cds_boundary_shift = 0`.
3. For each `i` in `range(len(snapshot))`:
   - `cdss_start = snapshot[i][0] + cds_boundary_shift`
   - `cdss_end = snapshot[i][1] + cds_boundary_shift`
   - **Straddle test:** if `(cdss_start <= cigar_accum_len) and (cdss_end >= cigar_accum_len)` — i.e. the `D` block falls inside this CDS's protein span — then `cds_boundary_shift += length` and `cdss_end += length` (the CDS's end is pushed out by the deletion length; subsequent CDSs inherit the accumulated shift via `cds_boundary_shift`).
   - Write `cdss_protein_aln_boundary[i] = (cdss_start, cdss_end)`.
4. Return `cdss_protein_aln_boundary`.

Gotcha: the cumulative `cds_boundary_shift` is applied to **every** CDS (each gets `+= cds_boundary_shift` on both start and end), but the `+= length` extension of `cdss_end` and the bump of `cds_boundary_shift` only fire for the CDS that contains `cigar_accum_len`. Because the snapshot is read-only, processing CDSs in index order is safe; mutating-while-reading (the pre-V2.10 bug) would double-count overlapping shifts.

#### 3.2.6 `lifton_parasail_align(...)` — the protein-path orchestrator (`align.py:202`)

Signature: `lifton_parasail_align(lifton_trans, db_entry, fai, ref_proteins, ref_trans_id) -> Lifton_Alignment | None`.

Algorithm:
1. `aln = None`. If `ref_trans_id not in ref_proteins.keys()`, return `None` (no reference protein to align against).
2. **Translate the annotation** via `LiftOn_translate(lifton_trans, fai, ref_proteins, ref_trans_id)` (`align.py:182`), which:
   - calls `lifton_trans.get_coding_seq(fai)` → `(coding_seq, cds_children, cdss_lens)`,
   - then *re-calls* `lifton_trans.get_coding_trans_seq(fai)` → `(coding_seq, _)` (overwriting `coding_seq` with the trans-method's coding sequence; this second call is what sets the CDS `frame` side effect — §3.1.2),
   - then `protein_seq = lifton_trans.translate_coding_seq(coding_seq)`,
   - returns `(str(ref_proteins[ref_trans_id]), protein_seq, cdss_lens, cds_children)`. (Note `cdss_lens`/`cds_children` come from the *first* `get_coding_seq` call; the coding sequence used for translation comes from the *second* `get_coding_trans_seq` call.)
3. `cdss_protein_boundary = get_cdss_protein_boundary(cdss_lens)` (§3.2.4).
4. If `protein_seq == None`, return `None` (empty coding sequence — nothing to align).
5. `extracted_parasail_res = parasail_align_protein_base(protein_seq, ref_protein_seq)` (§3.2.3; raises if either is empty).
6. Pull traceback strings: `alignment_query = res.traceback.query`, `alignment_comp = res.traceback.comp`, `alignment_ref = res.traceback.ref`.
7. **CIGAR replay to remap boundaries:** decode `cigar = res.cigar`, `decoded_cigar = cigar.decode.decode()`, `cigar_ls = list(Cigar(decoded_cigar).items())` (a list of `(length, symbol)` tuples via the `cigar.Cigar` library). Then:
   - `cigar_accum_len = 0`; `cdss_protein_aln_boundary = cdss_protein_boundary.copy()` (shallow copy of the dict; tuple values are immutable so this is safe).
   - For each `(length, symbol)` in `cigar_ls`: if `symbol == "D"`, call `cdss_protein_aln_boundary = adjust_cdss_protein_boundary(cdss_protein_aln_boundary, cigar_accum_len, length)` (§3.2.5). Then **always** `cigar_accum_len += length` (the accumulator advances over every block, `=`, `X`, `I`, `D` alike, before/after the adjust).
   - Gotcha (order): the `D` adjustment is applied **using the accumulator value as it stands at the start of that block** (`cigar_accum_len` is incremented *after* the adjust call), so a `D` block's boundary shift keys off the alignment column where the deletion *begins*.
8. **Identity:** `extracted_matches, extracted_length = get_id_fraction.get_AA_id_fraction(res.traceback.ref, res.traceback.query)` (§3.3.1; reference first, then query/target); `extracted_identity = extracted_matches / extracted_length`.
9. Build and return `Lifton_Alignment(extracted_identity, cds_children, alignment_query, alignment_comp, alignment_ref, cdss_protein_boundary, cdss_protein_aln_boundary, protein_seq, ref_protein_seq, db_entry)`.

#### 3.2.7 `protein_align(protein_seq, ref_protein_seq)` (`align.py:113`)

A lighter protein aligner used when no CDS-boundary remapping is needed (called from `align_coding_seq`, §3.2.9).
1. `extracted_parasail_res = parasail_align_protein_base(protein_seq, ref_protein_seq)`.
2. Pull `alignment_query/comp/ref` from `res.traceback`.
3. `extracted_matches, extracted_length = get_id_fraction.get_AA_id_fraction(res.traceback.ref, res.traceback.query)`; `extracted_identity = extracted_matches/extracted_length`.
4. Return `Lifton_Alignment(extracted_identity, None, alignment_query, alignment_comp, alignment_ref, None, None, protein_seq, ref_protein_seq, None)` — note `cds_children`, both boundary dicts, and `db_entry` are `None`.

#### 3.2.8 `trans_align(trans_seq, ref_trans_seq)` (`align.py:158`)

The DNA analogue.
1. `extracted_parasail_res = parasail_align_DNA_base(trans_seq, ref_trans_seq)` (sanitises both inputs; §3.2.3).
2. Pull `alignment_query/comp/ref` from `res.traceback`.
3. `extracted_matches, extracted_length = get_id_fraction.get_DNA_id_fraction(res.traceback.ref, res.traceback.query)` (§3.3.2; equal-length requirement); `extracted_identity = extracted_matches/extracted_length`.
4. Return `Lifton_Alignment(extracted_identity, None, alignment_query, alignment_comp, alignment_ref, None, None, trans_seq, ref_trans_seq, None)`.

#### 3.2.9 Transcript-level dispatch wrappers in `Lifton_TRANS`

`align_coding_seq(self, protein_seq, ref_protein_seq, lifton_status)` (`lifton_class.py:594`):
1. If `ref_protein_seq == "" or ref_protein_seq == None`: return `(None, None)` (`lifton_aa_aln`, `peps` both `None`).
2. Elif `protein_seq == "" or protein_seq == None`: return `(None, None)`.
3. Else: `peps = protein_seq.split("*")` (split translated protein into peptides at stop codons); `lifton_aa_aln = align.protein_align(protein_seq, ref_protein_seq)`; update running best identity `lifton_status.lifton_aa = max(lifton_status.lifton_aa, lifton_aa_aln.identity)`; return `(lifton_aa_aln, peps)`.

Gotcha: `lifton_aa` accumulates the **maximum** identity seen across multiple alignment attempts on the same transcript (e.g. across ORF-rescue iterations), whereas the DNA identity below is **overwritten** each call.

`align_trans_seq(self, trans_seq, ref_trans_seq, lifton_status)` (`lifton_class.py:608`):
1. If `ref_trans_seq == "" or ref_trans_seq == None`: return `None`.
2. Elif `trans_seq == "" or trans_seq == None`: return `None`.
3. Else: `lifton_tran_aln = align.trans_align(trans_seq, ref_trans_seq)`; set `lifton_status.lifton_dna = lifton_tran_aln.identity` (direct assignment, **not** `max`); return `lifton_tran_aln`.

#### 3.2.10 `Lifton_Alignment` data structure (`lifton_class.py:24`)

Constructed by all four aligners above. Positional constructor order is fixed; fields:

| Field (attr) | Constructor arg | Type | Meaning |
|---|---|---|---|
| `identity` | `extracted_identity` | float in [0,1] | matches / length from the identity function |
| `cds_children` | `cds_children` | list of feature entries or `None` | deep-copied CDS entries (only set by `lifton_parasail_align`) |
| `query_aln` | `alignment_query` | str | aligned query (annotated) string from traceback |
| `comp` | `alignment_comp` | str | parasail comparison/midline string |
| `ref_aln` | `alignment_ref` | str | aligned reference string from traceback |
| `cdss_protein_boundaries` | `cdss_protein_boundary` | dict `{idx:(start,end)}` or `None` | pre-CIGAR AA boundaries (protein path only) |
| `cdss_protein_aln_boundaries` | `cdss_protein_aln_boundary` | dict `{idx:(start,end)}` or `None` | post-CIGAR-D-shift AA boundaries |
| `query_seq` | `extracted_seq` | str | raw (un-gapped) query sequence aligned |
| `ref_seq` | `reference_seq` | str | raw (un-gapped) reference sequence aligned |
| `db_entry` | `db_entry` | feature or `None` | originating DB feature (protein path only) |

Method `write_alignment(self, outdir, tool_name, mutation, trans_id)` (`lifton_class.py:37`) dumps `ref_aln` and `query_aln` as a 2-record FASTA at `outdir/tool_name/mutation/trans_id.fa` (reference labelled `> Reference`, query labelled `> Target`).

### 3.3 Sequence-identity math (`lifton/get_id_fraction.py`)

Three functions each return a `(matches, denominator)` pair; the caller computes `identity = matches / denominator`. All three uppercase both inputs first. Each returns a **non-zero denominator** in degenerate cases so the caller's division never raises `ZeroDivisionError`. **Argument convention throughout: `reference` first, `target` second** — and callers pass `traceback.ref` then `traceback.query`, so within these functions `target` is the **query/annotated** aligned string.

#### 3.3.1 `get_AA_id_fraction(reference, target)` (`get_id_fraction.py:23`) — gap-collapsed BLAST protein identity

1. Uppercase both. `matches = 0`, `gaps_in_ref = 0`.
2. For `i, letter in enumerate(reference)`:
   - If `letter == '-'`: `gaps_in_ref += 1`.
   - If `letter == target[i]`: `matches += 1`. (A `-`-vs-`-` column counts as a match, since both equal `'-'`; this is the gap-collapse behaviour.)
   - **Early stop:** if `target[i] == "*"`: `break` (stop counting at the first stop codon in the target/annotated protein — downstream sequence after a premature stop is ignored).
3. **Zero guard:** if `max(len(reference), len(target)) == 0`, return `(matches, 1)`.
4. `total_length = max(len(reference), len(target)) - gaps_in_ref` (denominator is the longer of the two aligned strings **minus** the number of reference gaps — collapsing reference-side gaps out of the length).
5. **Second zero guard (V2.4):** if `total_length <= 0` (reference all gaps, or a stop-truncation left nothing), return `(matches, 1)`.
6. Return `(matches, total_length)`.

Formula: `identity = matches / (max(len(ref), len(target)) - gaps_in_ref)`, where `matches` counts equal aligned columns up to the first target `*`, and the denominator floors at 1.

#### 3.3.2 `get_DNA_id_fraction(reference, target)` (`get_id_fraction.py:49`) — equal-length BLAST DNA identity

1. Uppercase both.
2. **Length guard (V2.3, raises before the loop):** if `len(reference) != len(target)`, raise `ValueError("get_DNA_id_fraction: reference length (...) does not match target length (...). The two sequences must be aligned to equal length before identity is computed.")`. (parasail tracebacks are always equal length, so this should never fire in practice; it guards against a future refactor producing unequal strings.)
3. `matches = 0`; for `i, letter in enumerate(reference)`: if `letter == target[i]`: `matches += 1`. (No gap collapse, no `*` early-stop — straight column-by-column match count.)
4. If `max(len(reference), len(target)) == 0`, return `(matches, 1)`.
5. Return `(matches, max(len(reference), len(target)))`.

Formula: `identity = matches / len(aligned_string)` (both strings equal length, so the denominator is simply the alignment length, floored at 1).

#### 3.3.3 `get_partial_id_fraction(reference, target, start, end)` (`get_id_fraction.py:1`) — windowed protein identity

Scores a `[start, end)` slice of the aligned strings (used for per-CDS / per-region identity, keyed off the boundary dicts from §3.2). Note `start`/`end` are alignment-column indices.

1. Uppercase both. `matches = 0`, `gaps_in_ref = 0`.
2. For `i, letter in enumerate(reference[start:end])` (so the absolute index into `target` is `i + start`):
   - If `letter == '-'`: `gaps_in_ref += 1`.
   - If `letter == target[i+start]`: `matches += 1`.
   - **Early stop:** if `target[i+start] == "*"`: `break`.
3. `total_length = (end - start) - gaps_in_ref` (window width minus reference gaps in the window).
4. **Zero guard:** if `total_length == 0`, return `(matches, 1)`.
5. Return `(matches, total_length)`.

Formula: `identity = matches / ((end - start) - gaps_in_ref)`, matches counted over the window up to the first target `*`, denominator floored at 1.

Gotcha (`*` early-stop asymmetry): `get_AA_id_fraction` and `get_partial_id_fraction` both `break` at the first `target == "*"` (so a premature stop in the annotated protein truncates both the match count *and* — because the loop ends early — leaves `gaps_in_ref` reflecting only the prefix), but they do **not** break on `reference == "*"`. `get_DNA_id_fraction` never breaks on `*`. The three denominators are also structurally different: AA = `max(len) - gaps_in_ref`, partial = `(end-start) - gaps_in_ref`, DNA = `len` (no gap subtraction). These differences are byte-identity-relevant for the identity values that gate ORF rescue and miniprot chaining.
