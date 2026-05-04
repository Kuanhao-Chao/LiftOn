"""Phase 13 — protein-maximization (chaining algorithm) edge-case suite.

Anchored to the manuscript Methods §60-65 + Results §17-20 + Algorithm
S3 description. Each test exercises a single edge of the chaining
algorithm that the legacy Phase 5-12 suite did not cover. The
in-suite docstrings reference the manuscript paragraph that motivates
the test.

Manuscript ground truth (paraphrased):
  • Group CDSs from Liftoff and miniprot 5'→3' until both reach an
    endpoint where cumulative *aligned-amino-acid count* in the
    reference protein matches.
  • Per group, compute partial protein identity over the group's
    reference window; pick the higher-identity source.
  • Tie-break to Liftoff (so UTRs are preserved from the DNA-based
    map).
  • Concatenate all winning groups into the final LiftOn CDS list.
  • Empty / single-CDS inputs must NOT crash the algorithm.

The tests use synthetic, hand-crafted Lifton_Alignment objects so
the chaining algorithm is exercised in isolation without any
parasail / FASTA / IntervalTree dependencies.
"""

from __future__ import annotations

from types import SimpleNamespace

import pytest

from lifton import lifton_class, protein_maximization


# ---------------------------------------------------------------------------
# Helpers — construct synthetic Lifton_Alignment objects
# ---------------------------------------------------------------------------

def _fake_cds_child(seqid="chr1", start=1, end=99, strand="+"):
    """Build a gffutils-shaped CDS feature stub. Carries `attributes`
    so create_lifton_entries' parent-attribute inheritance works."""
    from unittest.mock import MagicMock
    f = MagicMock()
    f.seqid = seqid
    f.start = start
    f.end = end
    f.strand = strand
    f.attributes = {"Parent": ["tx1"], "ID": [f"cds_{start}"]}
    f.featuretype = "CDS"
    f.frame = "0"
    return f


def _fake_alignment(cds_children, cdss_protein_aln_boundaries,
                    ref_aln, query_aln, strand="+"):
    """Build a synthetic Lifton_Alignment with the minimum surface
    `chaining_algorithm` accesses."""
    aln = lifton_class.Lifton_Alignment(
        extracted_identity=0.0,
        cds_children=cds_children,
        alignment_query=query_aln,
        alignment_comp="|" * len(ref_aln),
        alignment_ref=ref_aln,
        cdss_protein_boundary={i: cdss_protein_aln_boundaries[i]
                                for i in range(len(cdss_protein_aln_boundaries))},
        cdss_protein_aln_boundary=cdss_protein_aln_boundaries,
        extracted_seq=query_aln.replace("-", ""),
        reference_seq=ref_aln.replace("-", ""),
        db_entry=SimpleNamespace(strand=strand),
    )
    return aln


# ---------------------------------------------------------------------------
# § 1 — Empty / minimal-input guards (manuscript §65 algorithmic safety)
# ---------------------------------------------------------------------------

class TestChainingEmptyAndSingleCDS:
    """Methods §65 silently assumes both Liftoff and miniprot have
    multiple CDS groups. The Phase 5-12 codebase added explicit
    guards (see protein_maximization.py docstring); these tests pin
    each guard branch."""

    def test_both_empty_returns_empty(self):
        l_aln = _fake_alignment([], [], "", "")
        m_aln = _fake_alignment([], [], "", "")
        result, chains = protein_maximization.chaining_algorithm(
            l_aln, m_aln, fai=None, DEBUG=False,
        )
        assert result == []
        assert chains == []

    def test_miniprot_empty_returns_liftoff_cds_unchanged(self):
        cds = [_fake_cds_child(start=1, end=99)]
        l_aln = _fake_alignment(cds, [(0, 33)], "M" * 33, "M" * 33)
        m_aln = _fake_alignment([], [], "", "")
        result, chains = protein_maximization.chaining_algorithm(
            l_aln, m_aln, fai=None, DEBUG=False,
        )
        # All Liftoff CDSs returned as-is, wrapped in Lifton_CDS
        assert len(result) == 1
        assert isinstance(result[0], lifton_class.Lifton_CDS)
        assert chains == []  # no comparison happened

    def test_liftoff_empty_returns_empty(self):
        cds = [_fake_cds_child(start=1, end=99)]
        l_aln = _fake_alignment([], [], "", "")
        m_aln = _fake_alignment(cds, [(0, 33)], "M" * 33, "M" * 33)
        result, chains = protein_maximization.chaining_algorithm(
            l_aln, m_aln, fai=None, DEBUG=False,
        )
        # Manuscript §20: tie-break to Liftoff; no Liftoff = no
        # output (miniprot-only goes through process_miniprot, not
        # the chaining algorithm).
        assert result == []

    def test_single_cds_each_side_compares_one_chunk(self):
        l_cds = [_fake_cds_child(start=1, end=99)]
        m_cds = [_fake_cds_child(start=1, end=99)]
        # Both have the same alignment quality
        l_aln = _fake_alignment(l_cds, [(0, 33)],
                                ref_aln="M" * 33, query_aln="M" * 33)
        m_aln = _fake_alignment(m_cds, [(0, 33)],
                                ref_aln="M" * 33, query_aln="M" * 33)
        result, chains = protein_maximization.chaining_algorithm(
            l_aln, m_aln, fai=None, DEBUG=False,
        )
        assert len(result) == 1
        assert len(chains) == 1
        # Tie-break to Liftoff → chain entry should mention "liftoff"
        assert "liftoff" in chains[0].lower()


# ---------------------------------------------------------------------------
# § 2 — Liftoff-wins / miniprot-wins selection (manuscript §65)
# ---------------------------------------------------------------------------

class TestChainingSelectionByIdentity:
    def test_miniprot_wins_when_higher_identity(self):
        """Methods §65: 'CDS groups with higher identity scores
        relative to the reference protein are selected for the final
        LiftOn annotation.'"""
        l_cds = [_fake_cds_child(start=1, end=99)]
        m_cds = [_fake_cds_child(start=1, end=99)]
        # Liftoff alignment has a mismatch in the middle
        l_aln = _fake_alignment(l_cds, [(0, 33)],
                                ref_aln="M" * 33,
                                query_aln="M" * 16 + "X" + "M" * 16)
        # Miniprot alignment is a perfect match
        m_aln = _fake_alignment(m_cds, [(0, 33)],
                                ref_aln="M" * 33,
                                query_aln="M" * 33)
        result, chains = protein_maximization.chaining_algorithm(
            l_aln, m_aln, fai=None, DEBUG=False,
        )
        assert len(result) == 1
        # Chain log mentions the winner
        assert "miniprot" in chains[0].lower()

    def test_liftoff_wins_when_tied(self):
        """Methods §20: 'In case of a tie, LiftOn prioritizes the
        Liftoff annotation, so that it will include UTRs in its
        output.'"""
        l_cds = [_fake_cds_child(start=1, end=99)]
        m_cds = [_fake_cds_child(start=1, end=99)]
        l_aln = _fake_alignment(l_cds, [(0, 33)],
                                ref_aln="M" * 33, query_aln="M" * 33)
        m_aln = _fake_alignment(m_cds, [(0, 33)],
                                ref_aln="M" * 33, query_aln="M" * 33)
        result, chains = protein_maximization.chaining_algorithm(
            l_aln, m_aln, fai=None, DEBUG=False,
        )
        assert "liftoff" in chains[0].lower()


# ---------------------------------------------------------------------------
# § 3 — Multi-CDS chunking (manuscript §65 grouping behaviour)
# ---------------------------------------------------------------------------

class TestChainingMultiCDS:
    """Methods §65: 'It starts by grouping the initial CDSs from both
    annotations until it encounters a boundary where the aligned
    amino acid counts in the reference protein up to that boundary
    match for both Liftoff and miniprot.'"""

    def test_two_cds_each_side_synced_at_matching_endpoint(self):
        # Both sides: 30 aa first CDS, 30 aa second CDS, identical
        # boundaries → one sync point at the first CDS boundary
        l_cds = [_fake_cds_child(start=1, end=90),
                 _fake_cds_child(start=200, end=289)]
        m_cds = [_fake_cds_child(start=1, end=90),
                 _fake_cds_child(start=200, end=289)]
        # Each CDS contributes 30 amino acids (90 bp / 3)
        # Boundaries (in aa coords): (0,30), (30,60)
        boundaries = [(0, 30), (30, 60)]
        l_aln = _fake_alignment(l_cds, boundaries,
                                ref_aln="M" * 60, query_aln="M" * 60)
        m_aln = _fake_alignment(m_cds, boundaries,
                                ref_aln="M" * 60, query_aln="M" * 60)
        result, chains = protein_maximization.chaining_algorithm(
            l_aln, m_aln, fai=None, DEBUG=False,
        )
        # Both CDSs from the winning side present
        assert len(result) == 2

    def test_three_cds_with_middle_winner_swap(self):
        """Build three matched groups where Liftoff wins the first
        and last, miniprot wins the middle. The chaining algorithm
        should still produce a coherent CDS list."""
        l_cds = [_fake_cds_child(start=1, end=90),
                 _fake_cds_child(start=200, end=289),
                 _fake_cds_child(start=400, end=489)]
        m_cds = [_fake_cds_child(start=1, end=90),
                 _fake_cds_child(start=200, end=289),
                 _fake_cds_child(start=400, end=489)]
        boundaries = [(0, 30), (30, 60), (60, 90)]
        # Liftoff perfect except middle window
        l_aln = _fake_alignment(
            l_cds, boundaries,
            ref_aln="M" * 90,
            query_aln="M" * 30 + "X" * 30 + "M" * 30,
        )
        # Miniprot perfect everywhere
        m_aln = _fake_alignment(
            m_cds, boundaries,
            ref_aln="M" * 90, query_aln="M" * 90,
        )
        result, chains = protein_maximization.chaining_algorithm(
            l_aln, m_aln, fai=None, DEBUG=False,
        )
        # 3 CDSs concatenated; order preserved
        assert len(result) == 3


# ---------------------------------------------------------------------------
# § 4 — Boundary helper functions (manuscript Algorithm S2)
# ---------------------------------------------------------------------------

class TestBoundaryHelpers:
    def test_get_protein_boundary_returns_aa_window(self):
        boundaries = [(0, 30), (30, 60), (60, 90)]
        aa_start, aa_end = protein_maximization.get_protein_boundary(
            boundaries, c_idx_last=0, c_idx=2, DEBUG=False,
        )
        # aa_start = boundaries[0][0]; aa_end = boundaries[1][1]
        assert (aa_start, aa_end) == (0, 60)

    def test_get_protein_reference_length_single_counts_non_gap(self):
        """The function reads `ref_seq` (the un-aligned reference
        protein, not the alignment), so it returns the count of
        non-gap characters in `ref_seq[0:ceil(boundary_end)]`. By
        construction `ref_seq` has gaps stripped, so this equals
        ceil(boundary_end) for typical inputs."""
        boundaries = [(0, 5), (5, 10), (10, 15)]
        aln = _fake_alignment(
            [], boundaries,
            ref_aln="MAGTA AAGTA AAGTA".replace(" ", ""),  # 15 chars
            query_aln="MAGTQAAGTQAAGTQ",
        )
        n = protein_maximization.get_protein_reference_length_single(
            aln, c_idx=0, DEBUG=False,
        )
        # ref_seq[0:5] = "MAGTA" → 5 non-gap chars
        assert n == 5

    def test_get_protein_reference_length_single_out_of_bounds(self):
        boundaries = [(0, 5)]
        aln = _fake_alignment([], boundaries,
                              ref_aln="MAGTA", query_aln="MAGTA")
        # Ask for an index past the boundaries list
        n = protein_maximization.get_protein_reference_length_single(
            aln, c_idx=99, DEBUG=False,
        )
        assert n == 0  # graceful fallback


# ---------------------------------------------------------------------------
# § 5 — process_m_l_children edge cases
# ---------------------------------------------------------------------------

class TestProcessMLChildrenEdges:
    def test_empty_chunk_returns_empty(self):
        l_aln = _fake_alignment([_fake_cds_child()], [(0, 33)],
                                ref_aln="M" * 33, query_aln="M" * 33)
        m_aln = _fake_alignment([_fake_cds_child()], [(0, 33)],
                                ref_aln="M" * 33, query_aln="M" * 33)
        result = protein_maximization.process_m_l_children(
            m_c_idx=0, m_c_idx_last=0, m_lifton_aln=m_aln,
            l_c_idx=0, l_c_idx_last=0, l_lifton_aln=l_aln,
            fai=None, chains=[], DEBUG=False,
        )
        # Empty chunk (last == idx) → no output
        assert result == []


# ---------------------------------------------------------------------------
# § 6 — Negative-strand chaining (manuscript Figure 1 implies strand-aware
# walk; create_lifton_entries reverses indices for "-" strand)
# ---------------------------------------------------------------------------

class TestChainingNegativeStrand:
    def test_negative_strand_does_not_crash(self):
        """Both annotations on the - strand. The chaining algorithm
        walks 5'→3' which on the - strand means high-coord → low-
        coord; the index reversal in create_lifton_entries handles
        this. The test pins that the algorithm completes without
        raising."""
        l_cds = [_fake_cds_child(start=1, end=99, strand="-")]
        m_cds = [_fake_cds_child(start=1, end=99, strand="-")]
        l_aln = _fake_alignment(l_cds, [(0, 33)],
                                ref_aln="M" * 33, query_aln="M" * 33,
                                strand="-")
        m_aln = _fake_alignment(m_cds, [(0, 33)],
                                ref_aln="M" * 33, query_aln="M" * 33,
                                strand="-")
        result, chains = protein_maximization.chaining_algorithm(
            l_aln, m_aln, fai=None, DEBUG=False,
        )
        assert len(result) == 1
