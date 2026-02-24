import pytest
from lifton.extract_sequence import merge_children_intervals, get_padding_length, get_dna_sequence
from Bio.Seq import Seq

class MockFeature:
    def __init__(self, seqid, start, end, strand="+"):
        self.seqid = seqid
        self.start = start
        self.end = end
        self.strand = strand

def test_merge_children_intervals_empty():
    assert merge_children_intervals([]) == []

def test_merge_children_intervals_non_overlapping():
    features = [ MockFeature("chr1", 10, 20), MockFeature("chr1", 30, 40) ]
    merged = merge_children_intervals(features)
    assert merged == [[10, 20], [30, 40]]

def test_merge_children_intervals_overlapping():
    features = [ MockFeature("chr1", 10, 25), MockFeature("chr1", 20, 40) ]
    merged = merge_children_intervals(features)
    assert merged == [[10, 40]]

def test_merge_children_intervals_contained():
    features = [ MockFeature("chr1", 10, 50), MockFeature("chr1", 20, 30) ]
    merged = merge_children_intervals(features)
    assert merged == [[10, 50]]

def test_get_padding_length():
    assert get_padding_length(0) == 0
    assert get_padding_length(1) == 2
    assert get_padding_length(2) == 1
    assert get_padding_length(3) == 0
    assert get_padding_length(4) == 2

def test_get_dna_sequence_positive_strand():
    fasta = {"chr1": "ACGTACGTACGTACGTACGT"} # 1-indexed: 1=A, 2=C, 3=G, 4=T...
    # Let's extract 2-5: CGTA
    parent = MockFeature("chr1", 1, 20, "+")
    features = [MockFeature("chr1", 2, 5, "+")]
    # Note get_dna_sequence effectively does: start-1 to end
    seq = get_dna_sequence(parent, fasta, features)
    assert seq == "CGTANN", f"Expected CGTANN (length 6), got {seq}"

def test_get_dna_sequence_negative_strand():
    fasta = {"chr1": "ACGTACGTACGTACGTACGT"}
    # extract 2-5: CGTA
    # negative strand reverse complement: TACG
    # pad to multiple of 3: TACGNN
    parent = MockFeature("chr1", 1, 20, "-")
    features = [MockFeature("chr1", 2, 5, "-")]
    seq = get_dna_sequence(parent, fasta, features)
    assert seq == "TACGNN", f"Expected TACGNN, got {seq}"

def test_get_dna_sequence_missing_chrom():
    fasta = {"chr2": "ACGT"}
    parent = MockFeature("chr1", 1, 20, "+")
    features = [MockFeature("chr1", 2, 5, "+")]
    seq = get_dna_sequence(parent, fasta, features)
    assert seq == "", f"Expected empty string for missing chrom, got {seq}"
