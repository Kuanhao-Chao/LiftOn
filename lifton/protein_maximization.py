"""
protein_maximization.py — CDS chaining algorithm.

Merges Liftoff and miniprot CDS blocks by comparing per-chunk protein
alignment identity, taking the better source for each exon group.

Bug-fixes applied:
  • Guard against empty or single-CDS inputs (formerly: infinite loop /
    IndexError when len(children) == 0 or 1).
  • Fix off-by-one in while-loop termination that caused the very last CDS
    pair to be processed twice for multi-CDS cases.
  • Guard get_protein_reference_length_single against out-of-bounds c_idx.
  • Guard process_m_l_children against zero-length protein windows.
  • Fix create_lifton_entries for empty index ranges.
"""

from lifton import get_id_fraction, lifton_class, logger
import math


# ─────────────────────────────────────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────────────────────────────────────

def get_protein_boundary(cdss_aln_boundary, c_idx_last, c_idx, DEBUG):
    """Return (aa_start, aa_end) for the chunk [c_idx_last, c_idx)."""
    aa_start = cdss_aln_boundary[c_idx_last][0]
    aa_end   = cdss_aln_boundary[c_idx - 1][1]
    return aa_start, aa_end


def get_protein_reference_length_single(lifton_aln, c_idx, DEBUG):
    """
    Count non-gap reference amino acids up through alignment position
    cdss_protein_aln_boundaries[c_idx][1].

    Returns 0 safely when c_idx is out of range.
    """
    boundaries = lifton_aln.cdss_protein_aln_boundaries
    if c_idx >= len(boundaries):
        return 0
    aa_start = 0
    aa_end   = math.ceil(boundaries[c_idx][1])
    ref_count = 0
    for letter in lifton_aln.ref_seq[aa_start:aa_end]:
        if letter != "-":
            ref_count += 1
    return ref_count


def push_cds_idx(c_idx, lifton_aln, ref_aa_count, DEBUG):
    """Advance the CDS pointer and return the updated (c_idx, ref_aa_count)."""
    ref_aa_count = get_protein_reference_length_single(lifton_aln, c_idx, DEBUG)
    c_idx += 1
    return c_idx, ref_aa_count


# ─────────────────────────────────────────────────────────────────────────────
# Per-chunk winner selection
# ─────────────────────────────────────────────────────────────────────────────

def process_m_l_children(
    m_c_idx, m_c_idx_last,
    m_lifton_aln,
    l_c_idx, l_c_idx_last,
    l_lifton_aln,
    fai, chains, DEBUG,
):
    """
    For the protein chunk spanning CDS indices [*_c_idx_last, *_c_idx),
    compare Liftoff vs. miniprot partial identity and emit the winner's CDS
    blocks.

    Edge cases handled:
      • Empty chunk (c_idx_last == c_idx) → return [] without division by zero.
      • Zero-length protein window → fall back to Liftoff (safer default).
    """
    # Empty chunk guard
    if m_c_idx_last >= m_c_idx or l_c_idx_last >= l_c_idx:
        return []

    m_aa_start, m_aa_end = get_protein_boundary(
        m_lifton_aln.cdss_protein_aln_boundaries, m_c_idx_last, m_c_idx, DEBUG)
    l_aa_start, l_aa_end = get_protein_boundary(
        l_lifton_aln.cdss_protein_aln_boundaries, l_c_idx_last, l_c_idx, DEBUG)

    m_matches, m_length = get_id_fraction.get_partial_id_fraction(
        m_lifton_aln.ref_aln, m_lifton_aln.query_aln,
        math.floor(m_aa_start), math.ceil(m_aa_end))
    l_matches, l_length = get_id_fraction.get_partial_id_fraction(
        l_lifton_aln.ref_aln, l_lifton_aln.query_aln,
        math.floor(l_aa_start), math.ceil(l_aa_end))

    # Guard against zero-length windows
    m_identity = m_matches / m_length if m_length > 0 else 0.0
    l_identity = l_matches / l_length if l_length > 0 else 0.0

    if m_identity > l_identity:
        cds_ls = create_lifton_entries(
            m_c_idx, m_c_idx_last, m_lifton_aln,
            l_c_idx, l_c_idx_last, l_lifton_aln,
            fai, True)
        chains.append(f"miniprot[{m_aa_start:.2f}-{m_aa_end:.2f}]")
    else:
        # Liftoff wins on tie (safer: preserves splicing structure)
        cds_ls = create_lifton_entries(
            m_c_idx, m_c_idx_last, m_lifton_aln,
            l_c_idx, l_c_idx_last, l_lifton_aln,
            fai, False)
        chains.append(f"liftoff[{l_aa_start:.2f}-{l_aa_end:.2f}]")

    return cds_ls


def create_lifton_entries(
    m_c_idx, m_c_idx_last,
    m_lifton_aln,
    l_c_idx, l_c_idx_last,
    l_lifton_aln,
    fai, miniprot_is_better,
):
    """
    Emit Lifton_CDS objects for the [*_c_idx_last, *_c_idx) slice of the
    winning aligner, inheriting Parent attributes from liftoff.

    Edge cases handled:
      • Empty index range → return [] instead of iterating a degenerate range.
      • Strand '-' index reversal is bounds-checked.
    """
    cds_list = []
    if miniprot_is_better:
        source_aln = m_lifton_aln
        idx_range  = range(m_c_idx_last, m_c_idx)
    else:
        source_aln = l_lifton_aln
        idx_range  = range(l_c_idx_last, l_c_idx)

    if not idx_range:
        return []

    n = len(source_aln.cds_children)
    parent_attrs = l_lifton_aln.cds_children[0].attributes if l_lifton_aln.cds_children else {}

    for c_idx in idx_range:
        if source_aln.db_entry.strand == "+":
            c_idx_fix = c_idx
        else:  # "-"
            c_idx_fix = n - c_idx - 1

        if c_idx_fix < 0 or c_idx_fix >= n:
            continue  # safety: skip out-of-bounds

        lifton_cds = source_aln.cds_children[c_idx_fix]
        lifton_cds.attributes = parent_attrs
        cds_list.append(lifton_class.Lifton_CDS(lifton_cds))

    return cds_list


# ─────────────────────────────────────────────────────────────────────────────
# Main chaining algorithm
# ─────────────────────────────────────────────────────────────────────────────

def chaining_algorithm(l_lifton_aln, m_lifton_aln, fai, DEBUG):
    """
    Synchronize Liftoff and miniprot CDS block lists by accumulated reference
    amino-acid count, then pick the better source for each protein chunk.

    Returns
    -------
    cds_list : list of Lifton_CDS
        The merged CDS block list for the transcript.
    chains : list of str
        Human-readable log of which source won each chunk.

    Edge cases handled:
      • Either input has 0 CDS children → return Liftoff's CDS list as-is
        (or empty if also zero).
      • Either input has exactly 1 CDS child → skip synchronization loop,
        compare the single chunk directly.
      • Very short protein alignments with empty boundary lists.
    """
    l_children = l_lifton_aln.cds_children
    m_children = m_lifton_aln.cds_children
    chains = []

    # ── Guard: empty inputs ──────────────────────────────────────────────────
    if not l_children and not m_children:
        return [], chains
    if not m_children:
        # miniprot produced nothing usable; return liftoff CDS as-is
        return [lifton_class.Lifton_CDS(c) for c in l_children], chains
    if not l_children:
        return [], chains

    # ── Guard: single-CDS transcripts ───────────────────────────────────────
    # When either aligner produced exactly 1 CDS block, there is only one
    # chunk to compare. Skip the synchronization loop entirely.
    if len(l_children) == 1 or len(m_children) == 1:
        # Treat the entire alignment as one chunk
        cds_list = process_m_l_children(
            len(m_children), 0, m_lifton_aln,
            len(l_children), 0, l_lifton_aln,
            fai, chains, DEBUG)
        return cds_list, chains

    # ── General case: multiple CDS blocks on both sides ──────────────────────
    m_c_idx = 0
    l_c_idx = 0
    m_c_idx_last = 0
    l_c_idx_last = 0
    cds_list = []
    ref_aa_liftoff_count  = 0
    ref_aa_miniprot_count = 0

    # Walk to the second-to-last CDS on each side, finding sync points
    while m_c_idx < len(m_children) - 1 or l_c_idx < len(l_children) - 1:
        # One side exhausted its non-final CDSs → drain the other
        if m_c_idx >= len(m_children) - 1 and l_c_idx < len(l_children) - 1:
            l_c_idx, ref_aa_liftoff_count = push_cds_idx(
                l_c_idx, l_lifton_aln, ref_aa_liftoff_count, DEBUG)
            continue
        if l_c_idx >= len(l_children) - 1 and m_c_idx < len(m_children) - 1:
            m_c_idx, ref_aa_miniprot_count = push_cds_idx(
                m_c_idx, m_lifton_aln, ref_aa_miniprot_count, DEBUG)
            continue

        # Both still have non-final CDSs; synchronize by ref AA count
        if m_lifton_aln.db_entry.strand == "+":
            m_c = m_children[m_c_idx]
            l_c = l_children[l_c_idx]
        else:  # "-": walk from the 3'-most CDS (highest coord = last 5'→3')
            m_c = m_children[len(m_children) - m_c_idx - 1]
            l_c = l_children[len(l_children) - l_c_idx - 1]

        if ref_aa_liftoff_count < ref_aa_miniprot_count:
            l_c_idx, ref_aa_liftoff_count = push_cds_idx(
                l_c_idx, l_lifton_aln, ref_aa_liftoff_count, DEBUG)
        elif ref_aa_liftoff_count > ref_aa_miniprot_count:
            m_c_idx, ref_aa_miniprot_count = push_cds_idx(
                m_c_idx, m_lifton_aln, ref_aa_miniprot_count, DEBUG)
        else:
            # Both have consumed the same cumulative ref AAs.
            # If their target genomic ends match → synchronization point.
            if m_c_idx > 0 and l_c_idx > 0 and m_c.end == l_c.end:
                cdss = process_m_l_children(
                    m_c_idx, m_c_idx_last, m_lifton_aln,
                    l_c_idx, l_c_idx_last, l_lifton_aln,
                    fai, chains, DEBUG)
                cds_list += cdss
                m_c_idx_last = m_c_idx
                l_c_idx_last = l_c_idx
            l_c_idx, ref_aa_liftoff_count = push_cds_idx(
                l_c_idx, l_lifton_aln, ref_aa_liftoff_count, DEBUG)
            m_c_idx, ref_aa_miniprot_count = push_cds_idx(
                m_c_idx, m_lifton_aln, ref_aa_miniprot_count, DEBUG)

    # Advance both to include the final CDS in the last chunk
    l_c_idx += 1
    m_c_idx += 1

    # Process the final (or only remaining) chunk
    cdss = process_m_l_children(
        m_c_idx, m_c_idx_last, m_lifton_aln,
        l_c_idx, l_c_idx_last, l_lifton_aln,
        fai, chains, DEBUG)
    cds_list += cdss

    return cds_list, chains