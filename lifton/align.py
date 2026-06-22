import os
import parasail
from Bio.Seq import Seq
from cigar import Cigar
from lifton import get_id_fraction, lifton_class, windowed_align
from lifton.exceptions import LiftOnAlignmentError


# Alignment kernel: `parasail.nw_trace_scan_sat` is full O(L²) Needleman–Wunsch
# in time AND memory, so a titin-scale alignment (protein ~35 k aa, transcript
# ~106 kb) costs tens of GB + tens of seconds and OOM-crashes typical machines
# (memory `lifton-giant-gene-align-blowup`). Above a length gate we route to the
# anchor-windowed aligner (`windowed_align`) — bounded in time/memory, and exact
# (== full DP) wherever anchors exist or a divergent region is below the giant
# boundary; only true giants take an approximate memory-bounded split.
#
# Iteration 3 makes "band everything" the DEFAULT: the gate is low (2500 aa /
# 8000 nt) so the whole O(L²) mid-tail is windowed, with a fine band (cap 1500)
# and the exact-DP fallback raised to the giant boundary. End-to-end 1.4–2.6×
# faster, 3.5–4.5× less peak RSS; identity-exact on same-species lifts,
# mean-neutral on cross-species (benchmarks/compare/fast_align_ab.py). Escape
# hatches restore the pre-Iteration-3 behaviour:
#   --full-dp-align / LIFTON_FULL_DP_ALIGN  -> exact giant-only path: full DP for
#       every non-giant gene (gate 8000 aa / 25000 nt); giants still windowed so
#       it can't OOM. configure_alignment(band=False).
#   LIFTON_ALIGN_WINDOW_AA/NT=<huge>        -> pure full DP incl. giants
#       (manuscript reproduction; OOM-prone on titin-scale genes).
# Tiny genes (fixtures, median ~476 aa) stay below the 2500 gate → exact full DP
# → the 24-cell byte-identity matrix stays green with no golden edit.

# Band-everything (default) vs giant-only (--full-dp-align) tuned constants.
_BAND_AA, _BAND_NT, _BAND_CAP = 2500, 8000, 1500
_GIANT_AA, _GIANT_NT, _GIANT_CAP = 8000, 25000, 6000
# Exact-DP fallback bound for anchor-less regions (band mode only): up to the
# giant boundary stays exact; only beyond it falls back to the bounded split.
_BAND_MAX_FULLDP_AA, _BAND_MAX_FULLDP_NT = _GIANT_AA, _GIANT_NT

# Explicit env pins win over the flag (force pure full DP for reproduction).
_ENV_PINNED_AA = "LIFTON_ALIGN_WINDOW_AA" in os.environ
_ENV_PINNED_NT = "LIFTON_ALIGN_WINDOW_NT" in os.environ
_ENV_PINNED_CAP = "LIFTON_ALIGN_WINDOW_CAP" in os.environ

# Effective thresholds. Default = band-everything; explicit env vars win.
_ALIGN_WINDOW_AA = int(os.environ.get("LIFTON_ALIGN_WINDOW_AA", str(_BAND_AA)))
_ALIGN_WINDOW_NT = int(os.environ.get("LIFTON_ALIGN_WINDOW_NT", str(_BAND_NT)))
# Per-type exact-DP fallback bound passed into windowed_traceback (None → the
# aligner defaults it to WINDOW_CAP, i.e. fallback off / giant-only behaviour).
_ALIGN_MAX_FULLDP_AA = _BAND_MAX_FULLDP_AA
_ALIGN_MAX_FULLDP_NT = _BAND_MAX_FULLDP_NT
_FAST_ALIGN_ACTIVE = True


def configure_alignment(band=True):
    """Switch the alignment kernel between "band everything" (default, Iteration
    3) and the exact giant-only path (--full-dp-align, pre-Iteration-3). Called
    from `lifton.run_all_lifton_steps` only when --full-dp-align flips it off;
    the import-time module state is already band-everything. Symmetric, so it is
    safe to toggle in tests. Explicit `LIFTON_ALIGN_WINDOW_{AA,NT,CAP}` env vars
    are honoured over the flag (pure-full-DP reproduction escape hatch)."""
    global _ALIGN_WINDOW_AA, _ALIGN_WINDOW_NT
    global _ALIGN_MAX_FULLDP_AA, _ALIGN_MAX_FULLDP_NT, _FAST_ALIGN_ACTIVE
    _FAST_ALIGN_ACTIVE = bool(band)
    aa, nt, cap = (_BAND_AA, _BAND_NT, _BAND_CAP) if band else \
                  (_GIANT_AA, _GIANT_NT, _GIANT_CAP)
    if not _ENV_PINNED_AA:
        _ALIGN_WINDOW_AA = aa
    if not _ENV_PINNED_NT:
        _ALIGN_WINDOW_NT = nt
    if not _ENV_PINNED_CAP:
        windowed_align.set_window_cap(cap)
    if band:
        _ALIGN_MAX_FULLDP_AA = _BAND_MAX_FULLDP_AA
        _ALIGN_MAX_FULLDP_NT = _BAND_MAX_FULLDP_NT
    else:
        _ALIGN_MAX_FULLDP_AA = None
        _ALIGN_MAX_FULLDP_NT = None


# Backward-compatible alias (the flag was --fast-align during development).
def configure_fast_align(enabled=True):
    """Deprecated alias of configure_alignment(band=...)."""
    configure_alignment(band=enabled)


# Apply the default alignment mode at import (single source of truth, incl. the
# window cap). Band-everything by default; LIFTON_FULL_DP_ALIGN selects the exact
# giant-only baseline for subprocess benchmarks (the CLI --full-dp-align flag
# does the same in-process). Explicit LIFTON_ALIGN_WINDOW_* env still wins.
_FULL_DP_ALIGN_ENV = os.environ.get(
    "LIFTON_FULL_DP_ALIGN", "") not in ("", "0", "false", "False")
configure_alignment(band=not _FULL_DP_ALIGN_ENV)


# V4.2 fix: bases acceptable to parasail's ACGTN matrix below. Anything
# outside this set (IUPAC ambiguity codes R/Y/S/W/K/M/B/D/H/V, gap '-',
# space, lowercase, ...) is normalised to N before alignment. Crashing
# parasail's C kernel on an IUPAC code is the worst-case behaviour the
# Phase 13.5A audit flagged as High.
_PARASAIL_DNA_ALPHABET = frozenset("ACGTN*")


def _sanitise_for_parasail_dna(seq: str) -> str:
    """Normalise a nucleotide sequence to the {A,C,G,T,N,*} alphabet.

    parasail.matrix_create("ACGTN*", ...) only knows scores for these
    characters; any other byte (e.g. IUPAC ambiguity codes R/Y/W/S/K/M)
    can either crash the C kernel or score garbage. Coerce them to N
    so the kernel handles them as 'unknown nucleotide'.
    """
    if not seq:
        return seq
    upper = seq.upper()
    if all(ch in _PARASAIL_DNA_ALPHABET for ch in upper):
        return upper
    return "".join(
        ch if ch in _PARASAIL_DNA_ALPHABET else "N" for ch in upper
    )


def adjust_cdss_protein_boundary(cdss_protein_aln_boundary, cigar_accum_len, length):
    """
        This function adjust CDS protein boundaries based on CIGAR string on protein alignment.

        Parameters:
        - cdss_protein_aln_boundary: CDS protein alignment boundary
        - cigar_accum_len: accumulated length of CIGAR string
        - length: length of CIGAR string

        Returns:
        cdss_protein_aln_boundary: adjusted CDS protein alignment boundary
    """
    # V2.10 fix: snapshot the input dict so each boundary is read from
    # the PRE-call state. Previously the loop both read and wrote to
    # `cdss_protein_aln_boundary[i]`, and `cds_boundary_shift` could
    # compound when the same boundary overlapped multiple D-blocks
    # processed in successive calls referencing the mutated dict.
    snapshot = {i: cdss_protein_aln_boundary[i]
                for i in range(len(cdss_protein_aln_boundary))}
    cds_boundary_shift = 0
    for i in range(len(snapshot)):
        cdss_start = snapshot[i][0] + cds_boundary_shift
        cdss_end = snapshot[i][1] + cds_boundary_shift
        if (cdss_start <= cigar_accum_len) and (cdss_end >= cigar_accum_len):
            cds_boundary_shift += length
            cdss_end += length
        cdss_protein_aln_boundary[i] = (cdss_start, cdss_end)
    return cdss_protein_aln_boundary


def get_cdss_protein_boundary(cdss_lens):
    """
        This function maps the CDSs boundaries on to the protein.

        Parameters:
        - cdss_lens: list of CDS lengths

        Returns:
        cdss_protein_boundary: CDS protein boundary
    """
    cdss_cumulative = [sum(cdss_lens[:i+1]) for i in range(len(cdss_lens))]
    cdss_cumulative_div = [x / 3 for x in cdss_cumulative]
    cdss_protein_boundary = {}
    for idx in range(len(cdss_cumulative_div)):
        start = cdss_cumulative_div[idx-1] if idx > 0 else 0
        end = cdss_cumulative_div[idx]
        cdss_protein_boundary[idx] = (start, end)
    return cdss_protein_boundary


def parasail_align_protein_base(protein_seq, ref_protein_seq):
    """
        This is the base function that uses parasail to align the protein sequence to the reference protein sequence.

        Parameters:
        - protein_seq: protein sequence in string format
        - ref_protein_seq: reference protein sequence in string format

        Returns:
        extracted_parasail_res: parasail alignment result
    """
    matrix = parasail.Matrix("blosum62")
    gap_open = 11
    gap_extend = 1
    # V2.11 fix: empty inputs are upstream programming errors (the
    # CDS produced no codons). Silently coercing to "*" hid the bug
    # behind a near-zero alignment identity. Raise explicitly so the
    # caller can attribute the cause.
    if protein_seq == "" or ref_protein_seq == "":
        raise LiftOnAlignmentError(
            "parasail_align_protein_base: refusing to align empty "
            "protein sequence — caller must provide non-empty query "
            "and reference."
        )
    # Return: (Query, Target)
    if max(len(protein_seq), len(ref_protein_seq)) > _ALIGN_WINDOW_AA:
        extracted_parasail_res = windowed_align.windowed_traceback(
            protein_seq, ref_protein_seq, gap_open, gap_extend, matrix,
            is_dna=False, max_fulldp=_ALIGN_MAX_FULLDP_AA)
    else:
        extracted_parasail_res = parasail.nw_trace_scan_sat(protein_seq, ref_protein_seq, gap_open, gap_extend, matrix)
    return extracted_parasail_res


def protein_align(protein_seq, ref_protein_seq):
    """
        This function aligns the protein sequence to the reference protein sequence and extracts the alignment information.

        Parameters:
        - protein_seq: protein sequence in string format
        - ref_protein_seq: reference protein sequence in string format

        Returns:
        lifton_aln: Lifton_Alignment object
    """
    extracted_parasail_res = parasail_align_protein_base(protein_seq, ref_protein_seq)
    alignment_query = extracted_parasail_res.traceback.query
    alignment_comp = extracted_parasail_res.traceback.comp
    alignment_ref = extracted_parasail_res.traceback.ref
    extracted_matches, extracted_length = get_id_fraction.get_AA_id_fraction(extracted_parasail_res.traceback.ref, extracted_parasail_res.traceback.query)
    extracted_identity = extracted_matches/extracted_length
    lifton_aln = lifton_class.Lifton_Alignment(extracted_identity, None, alignment_query, alignment_comp, alignment_ref, None, None, protein_seq, ref_protein_seq, None)
    return lifton_aln


def parasail_align_DNA_base(trans_seq, ref_trans_seq):
    """
        This is the base function that uses parasail to align the transcript sequence to the reference transcript sequence.

        Parameters:
        - trans_seq: transcript sequence in string format
        - ref_trans_seq: reference transcript sequence in string format

        Returns:
        extracted_parasail_res: parasail alignment result
    """
    # V4.2 fix: extend the matrix alphabet to include 'N' AND sanitise
    # both inputs so any IUPAC ambiguity code present in real-world
    # reference data is converted to 'N' before reaching the C kernel.
    matrix = parasail.matrix_create("ACGTN*", 1, -3)
    gap_open = 5
    gap_extend = 2
    trans_seq = _sanitise_for_parasail_dna(trans_seq)
    ref_trans_seq = _sanitise_for_parasail_dna(ref_trans_seq)
    # Return: (Query, Target)
    if max(len(trans_seq), len(ref_trans_seq)) > _ALIGN_WINDOW_NT:
        extracted_parasail_res = windowed_align.windowed_traceback(
            trans_seq, ref_trans_seq, gap_open, gap_extend, matrix,
            is_dna=True, max_fulldp=_ALIGN_MAX_FULLDP_NT)
    else:
        extracted_parasail_res = parasail.nw_trace_scan_sat(trans_seq, ref_trans_seq, gap_open, gap_extend, matrix)
    return extracted_parasail_res


def trans_align(trans_seq, ref_trans_seq):
    """
        This function aligns the transcript sequence to the reference transcript sequence and extracts the alignment information.

        Parameters:
        - trans_seq: transcript sequence in string format
        - ref_trans_seq: reference transcript sequence in string format

        Returns:
        lifton_aln: Lifton_Alignment object
    """
    extracted_parasail_res = parasail_align_DNA_base(trans_seq, ref_trans_seq)
    alignment_query = extracted_parasail_res.traceback.query
    alignment_comp = extracted_parasail_res.traceback.comp
    alignment_ref = extracted_parasail_res.traceback.ref
    # print("alignment_query: ", alignment_query) 
    # print("alignment_comp: ", alignment_comp) 
    # print("alignment_ref: ", alignment_ref) 
    extracted_matches, extracted_length = get_id_fraction.get_DNA_id_fraction(extracted_parasail_res.traceback.ref, extracted_parasail_res.traceback.query)
    extracted_identity = extracted_matches/extracted_length
    lifton_aln = lifton_class.Lifton_Alignment(extracted_identity, None, alignment_query, alignment_comp, alignment_ref, None, None, trans_seq, ref_trans_seq, None)
    return lifton_aln


def LiftOn_translate(lifton_trans, fai, ref_proteins, ref_trans_id):
    """
        This function translates the given annotation.
        
        Parameters:
        - lifton_trans: LiftOn transcript instance
        - fai: fasta index
        - ref_proteins: reference protein dictionary
        - ref_trans_id: reference transcript ID

        Returns:
        Referece protein sequence (in string), translated protein sequence (in string), CDS lengths, CDS children
    """
    coding_seq, cds_children, cdss_lens = lifton_trans.get_coding_seq(fai)
    coding_seq, _ = lifton_trans.get_coding_trans_seq(fai)
    protein_seq = lifton_trans.translate_coding_seq(coding_seq)
    # print("protein_seq: ", protein_seq)
    return str(ref_proteins[ref_trans_id]), protein_seq, cdss_lens, cds_children


def lifton_parasail_align(lifton_trans, db_entry, fai, ref_proteins, ref_trans_id):
    """
        This function aligns the annotated protein sequence to the reference protein sequence

        Parameters:
        - lifton_trans: LiftOn transcript instance
        - db_entry: database entry
        - fai: fasta index
        - ref_proteins: reference protein dictionary
        - ref_trans_id: reference transcript ID
        - lifton_status: LiftOn status

        Returns:
        aln: Lifton_Alignment object
    """
    # Step 1: Check if the ref_trans_id is in the reference protein dictionary
    aln = None
    if ref_trans_id not in ref_proteins.keys():
        return aln
    # Step 2: Translate the given annotation
    ref_protein_seq, protein_seq, cdss_lens, cds_children = LiftOn_translate(lifton_trans, fai, ref_proteins, ref_trans_id)
    # Step 3: Get the corresponding protein boundaries for each CDS
    cdss_protein_boundary = get_cdss_protein_boundary(cdss_lens)
    # Step 4: Protein alignment
    if protein_seq == None:
        return aln
    extracted_parasail_res = parasail_align_protein_base(protein_seq, ref_protein_seq)
    # Step 5: Extract the alignment information
    alignment_query = extracted_parasail_res.traceback.query
    alignment_comp = extracted_parasail_res.traceback.comp
    alignment_ref = extracted_parasail_res.traceback.ref
    cigar = extracted_parasail_res.cigar
    decoded_cigar = cigar.decode.decode()
    cigar_ls = list(Cigar(decoded_cigar).items())
    # Step 6: Change the CDS protein boundaries based on CIGAR string.
    cigar_accum_len = 0
    cdss_protein_aln_boundary = cdss_protein_boundary.copy()
    for length, symbol in cigar_ls:
        if symbol == "D":
            # print(length, symbol)
            cdss_protein_aln_boundary = adjust_cdss_protein_boundary(cdss_protein_aln_boundary, cigar_accum_len, length)
        cigar_accum_len += length
        # print("cigar_accum_len: ", cigar_accum_len)
    extracted_matches, extracted_length = get_id_fraction.get_AA_id_fraction(extracted_parasail_res.traceback.ref, extracted_parasail_res.traceback.query)
    extracted_identity = extracted_matches/extracted_length
    aln = lifton_class.Lifton_Alignment(extracted_identity, cds_children, alignment_query, alignment_comp, alignment_ref, cdss_protein_boundary, cdss_protein_aln_boundary, protein_seq, ref_protein_seq, db_entry)
    return aln
