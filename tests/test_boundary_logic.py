"""Phase 4.5 Step 1 — coordinate boundary tests.

Exercises the 1-based / 0-based conversion at extract_sequence.py:81,
chromosome edges, single-base-pair features, and off-end coordinates.
"""

from __future__ import annotations

from types import SimpleNamespace

import pytest
from pyfaidx import Fasta

from lifton import annotation, extract_sequence


class _IntervalStub:
    def __init__(self, start, end):
        self.start = start
        self.end = end


# ---------------------------------------------------------------------------
# 1-based / 0-based conversion
# ---------------------------------------------------------------------------

class TestOneBasedConversion:
    def test_pyfaidx_zero_based_slice_matches_first_base(self, fasta_standard):
        fa = Fasta(str(fasta_standard))
        # 1-based positions 101..103 == zero-based slice [100:103]
        seq_via_lifton = extract_sequence.get_dna_sequence(
            SimpleNamespace(seqid="chr1", strand="+"), fa,
            [_IntervalStub(101, 103)],
        )
        # The first three bases at position 101 form the start codon ATG
        # (per fixtures._build_chrom_with_gene). Strip trailing N padding.
        assert seq_via_lifton.rstrip("N").startswith("ATG")

    def test_position_one_is_first_base(self, fasta_short_chrom):
        fa = Fasta(str(fasta_short_chrom))
        seq = extract_sequence.get_dna_sequence(
            SimpleNamespace(seqid="chr1", strand="+"), fa,
            [_IntervalStub(1, 3)],
        )
        # tiny.fa starts with "ATG"
        assert seq.rstrip("N") == "ATG"


# ---------------------------------------------------------------------------
# Chromosome edges
# ---------------------------------------------------------------------------

class TestChromosomeEdges:
    def test_feature_at_start_position_1(self, fasta_short_chrom,
                                         gff_at_chrom_edge):
        ann = annotation.Annotation(
            str(gff_at_chrom_edge), False, False, "create_unique", None,
            True, False,
        )
        gene = ann.db_connection["gE"]
        assert gene.start == 1
        # extract DNA across the same range
        fa = Fasta(str(fasta_short_chrom))
        seq = extract_sequence.get_dna_sequence(
            SimpleNamespace(seqid="chr1", strand="+"), fa,
            [_IntervalStub(1, 30)],
        )
        assert len(seq.rstrip("N")) == 30

    def test_feature_at_end_of_chromosome(self, fasta_short_chrom):
        fa = Fasta(str(fasta_short_chrom))
        # Chromosome is 50 bp; ask for last 10 bp (positions 41..50)
        seq = extract_sequence.get_dna_sequence(
            SimpleNamespace(seqid="chr1", strand="+"), fa,
            [_IntervalStub(41, 50)],
        )
        # Strip padding; should be exactly 10 bp of sequence
        assert len(seq.rstrip("N")) == 10

    def test_feature_off_end_returns_truncated_string(self, fasta_short_chrom):
        """Feature whose end exceeds chromosome length. pyfaidx truncates
        silently. We pin that downstream code does not crash."""
        fa = Fasta(str(fasta_short_chrom))
        seq = extract_sequence.get_dna_sequence(
            SimpleNamespace(seqid="chr1", strand="+"), fa,
            [_IntervalStub(40, 1000)],
        )
        # Result is at most chrom_len-39 bp of real sequence + N padding.
        # Crucially: no exception, and length is a multiple of 3 (padded).
        assert len(seq) % 3 == 0


# ---------------------------------------------------------------------------
# 1-bp / single-base features
# ---------------------------------------------------------------------------

class TestSingleBasePairFeatures:
    def test_one_bp_feature_extract(self, fasta_short_chrom):
        fa = Fasta(str(fasta_short_chrom))
        seq = extract_sequence.get_dna_sequence(
            SimpleNamespace(seqid="chr1", strand="+"), fa,
            [_IntervalStub(1, 1)],
        )
        # 1 bp + 2 bp of N padding to reach codon multiple
        assert len(seq) == 3
        assert seq[0] == "A"  # tiny.fa starts with 'A'TG -> position 1 is 'A'

    def test_one_bp_gff_loads(self, gff_one_bp_feature):
        ann = annotation.Annotation(
            str(gff_one_bp_feature), False, False, "create_unique", None,
            True, False,
        )
        feat = ann.db_connection["g1bp"]
        assert feat.start == 100 and feat.end == 100


# ---------------------------------------------------------------------------
# Padding length boundary cases
# ---------------------------------------------------------------------------

class TestPaddingBoundary:
    @pytest.mark.parametrize("L,expected_total_mod3", [
        (0, 0), (1, 0), (2, 0), (3, 0), (100, 0), (101, 0), (102, 0),
    ])
    def test_padding_always_yields_codon_multiple(self, L,
                                                  expected_total_mod3):
        pad = extract_sequence.get_padding_length(L)
        assert (L + pad) % 3 == expected_total_mod3
        assert pad in {0, 1, 2}


# ---------------------------------------------------------------------------
# Missing chromosome path (already covered, expanded)
# ---------------------------------------------------------------------------

class TestMissingChromosome:
    def test_unknown_chrom_yields_empty_string(self, fasta_short_chrom):
        fa = Fasta(str(fasta_short_chrom))
        seq = extract_sequence.get_dna_sequence(
            SimpleNamespace(seqid="chrUNKNOWN", strand="+"), fa,
            [_IntervalStub(1, 10)],
        )
        # Per legacy contract — empty string, no exception
        assert seq == ""
