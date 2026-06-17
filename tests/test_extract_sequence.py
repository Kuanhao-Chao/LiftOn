"""Unit tests for lifton.extract_sequence."""

from __future__ import annotations

from types import SimpleNamespace

import pytest
from pyfaidx import Fasta

from lifton import annotation, extract_sequence


# ---------------------------------------------------------------------------
# determine_file_format
# ---------------------------------------------------------------------------

class TestDetermineFileFormat:
    def test_detects_gff3(self, gff_standard):
        assert extract_sequence.determine_file_format(str(gff_standard)) == "GFF format"

    def test_detects_gtf(self, tmp_path):
        fp = tmp_path / "f.gtf"
        fp.write_text(
            '#comment\n'
            'chr1\ttest\texon\t100\t200\t.\t+\t.\tgene_id "g1"; transcript_id "t1";\n'
        )
        assert extract_sequence.determine_file_format(str(fp)) == "GTF format"

    def test_unrecognised_file_defaults_to_gff(self, tmp_path):
        # Post-merge: determine_file_format defaults to "GFF format"
        # when no clear indicators are found (instead of returning
        # "Unknown format").
        fp = tmp_path / "f.txt"
        fp.write_text("not\ta\tgff\n")
        assert extract_sequence.determine_file_format(str(fp)) == "GFF format"


# ---------------------------------------------------------------------------
# merge_children_intervals
# ---------------------------------------------------------------------------

class _IntervalStub:
    def __init__(self, start, end):
        self.start = start
        self.end = end


class TestMergeChildrenIntervals:
    def test_disjoint_kept_separate(self):
        merged = extract_sequence.merge_children_intervals(
            [_IntervalStub(1, 10), _IntervalStub(20, 30)]
        )
        assert merged == [[1, 10], [20, 30]]

    def test_overlap_merged(self):
        merged = extract_sequence.merge_children_intervals(
            [_IntervalStub(1, 15), _IntervalStub(10, 20)]
        )
        assert merged == [[1, 20]]

    def test_unsorted_input_sorted_first(self):
        merged = extract_sequence.merge_children_intervals(
            [_IntervalStub(50, 60), _IntervalStub(1, 10)]
        )
        assert merged == [[1, 10], [50, 60]]

    def test_empty_input(self):
        assert extract_sequence.merge_children_intervals([]) == []


# ---------------------------------------------------------------------------
# get_padding_length
# ---------------------------------------------------------------------------

class TestGetPaddingLength:
    @pytest.mark.parametrize("L,expected", [
        (0, 0), (1, 2), (2, 1), (3, 0), (4, 2), (5, 1), (6, 0),
    ])
    def test_pads_to_codon_multiple(self, L, expected):
        assert extract_sequence.get_padding_length(L) == expected


# ---------------------------------------------------------------------------
# get_dna_sequence
# ---------------------------------------------------------------------------

class TestGetDnaSequence:
    def test_forward_strand_concat(self, fasta_standard):
        fa = Fasta(str(fasta_standard))
        # Two exon-shaped intervals (101..199, 301..399), strand +
        parent = SimpleNamespace(seqid="chr1", strand="+")
        children = [_IntervalStub(101, 199), _IntervalStub(301, 399)]
        seq = extract_sequence.get_dna_sequence(parent, fa, children)
        # First exon starts with ATG, second ends with TAA before any padding
        assert seq.startswith("ATG")
        assert seq[:198].endswith("TAA")  # 99 + 99 nt, no padding (already %3==0)
        # Total length must be multiple of 3 (padded if needed)
        assert len(seq) % 3 == 0

    def test_reverse_strand_complemented(self, fasta_standard):
        fa = Fasta(str(fasta_standard))
        parent = SimpleNamespace(seqid="chr1", strand="-")
        children = [_IntervalStub(101, 199)]
        seq = extract_sequence.get_dna_sequence(parent, fa, children)
        # ATG on + strand becomes CAT on the reverse complement (at the END
        # of the forward exon). Just check we get a different string back.
        fwd = extract_sequence.get_dna_sequence(
            SimpleNamespace(seqid="chr1", strand="+"), fa, children
        )
        assert seq != fwd
        assert len(seq) == len(fwd)

    def test_missing_chromosome_returns_empty(self, fasta_missing_chrom):
        fa = Fasta(str(fasta_missing_chrom))
        parent = SimpleNamespace(seqid="chr1", strand="+")
        children = [_IntervalStub(1, 10)]
        # Empty string returned per legacy contract
        assert extract_sequence.get_dna_sequence(parent, fa, children) == ""


# ---------------------------------------------------------------------------
# get_protein_sequence
# ---------------------------------------------------------------------------

class TestGetProteinSequence:
    def test_translates_well_formed_orf(self, fasta_standard):
        fa = Fasta(str(fasta_standard))
        parent = SimpleNamespace(seqid="chr1", strand="+")
        children = [_IntervalStub(101, 199), _IntervalStub(301, 399)]
        protein = extract_sequence.get_protein_sequence(parent, fa, children)
        assert protein.startswith("M")
        assert protein.endswith("*")
        assert protein.count("*") == 1


# ---------------------------------------------------------------------------
# extract_features (integration with gffutils)
# ---------------------------------------------------------------------------

class TestExtractFeatures:
    def test_standard_gene_yields_dicts(self, gff_standard, fasta_standard):
        ref_db = annotation.Annotation(
            str(gff_standard), False, False, "create_unique", None, True, False,
        )
        fa = Fasta(str(fasta_standard))
        ref_trans, ref_proteins = extract_sequence.extract_features(
            ref_db, ["gene"], fa
        )
        # Recursion descends gene -> mRNA -> exon/CDS, so the populated key
        # is the transcript id ("tx1"), not the gene id ("gene1")
        assert "tx1" in ref_trans
        assert "tx1" in ref_proteins
        assert ref_proteins["tx1"].startswith("M")
        assert ref_proteins["tx1"].endswith("*")

    def test_noncoding_gene_yields_trans_no_protein(
            self, gff_noncoding, fasta_standard):
        ref_db = annotation.Annotation(
            str(gff_noncoding), False, False, "create_unique", None, True, False,
        )
        fa = Fasta(str(fasta_standard))
        ref_trans, ref_proteins = extract_sequence.extract_features(
            ref_db, ["gene"], fa
        )
        assert "nctx" in ref_trans
        assert "nctx" not in ref_proteins


# ---------------------------------------------------------------------------
# Iteration 20 regression: a gene-like type (ncRNA, pseudogene) that ALSO
# occurs as a CHILD of a gene must not be extracted twice. On full Arabidopsis
# RefSeq the ncRNA `rna-NR_143387.1` is a child of `gene-AT1G04425` AND the
# auto-detected gene-like type list includes `ncRNA`, so the old code wrote the
# transcript twice -> pyfaidx `Duplicate key` crash. The fix skips Parent-bearing
# instances at the top-level enumeration.
# ---------------------------------------------------------------------------

@pytest.fixture
def gff_ncrna_child_of_gene(tmp_path: Path) -> Path:
    """A gene whose single transcript is an `ncRNA` child (mirrors
    gene-AT1G04425 -> rna-NR_143387.1). `ncRNA` is therefore both a child
    type AND a candidate gene-like top-level type."""
    fp = tmp_path / "ref_ncrna_child.gff3"
    fp.write_text(
        "##gff-version 3\n"
        "chr1\ttest\tgene\t101\t399\t.\t+\t.\tID=gene1;gene_biotype=ncRNA\n"
        "chr1\ttest\tncRNA\t101\t399\t.\t+\t.\tID=rna1;Parent=gene1\n"
        "chr1\ttest\texon\t101\t199\t.\t+\t.\tID=exon1;Parent=rna1\n"
        "chr1\ttest\texon\t301\t399\t.\t+\t.\tID=exon2;Parent=rna1\n"
    )
    return fp


class TestGeneLikeChildNotDoubleExtracted:
    def _db(self, gff):
        return annotation.Annotation(
            str(gff), False, False, "create_unique", None, True, False,
        )

    def test_streaming_no_duplicate_fasta_key(
            self, gff_ncrna_child_of_gene, fasta_standard, tmp_path):
        """extract_features_to_fasta with features=['gene','ncRNA'] must write
        the child ncRNA transcript exactly once; pyfaidx.Fasta then opens the
        result without raising `Duplicate key` (the crash this fixes)."""
        ref_db = self._db(gff_ncrna_child_of_gene)
        fa = Fasta(str(fasta_standard))
        out_dir = tmp_path / "out"
        trans_path, _ = extract_sequence.extract_features_to_fasta(
            ref_db, ["gene", "ncRNA"], fa, str(out_dir),
        )
        headers = [ln[1:].strip() for ln in open(trans_path)
                   if ln.startswith(">")]
        assert headers.count("rna1") == 1, headers
        # The actual crash was here: pyfaidx indexing a dup-keyed FASTA.
        Fasta(str(trans_path))  # must not raise ValueError("Duplicate key ...")

    def test_streaming_gene_like_superset_equals_gene_only(
            self, gff_ncrna_child_of_gene, fasta_standard, tmp_path):
        """Adding the child type `ncRNA` to the feature list must be a no-op vs
        lifting `gene` alone (the child is already emitted under its gene)."""
        ref_db = self._db(gff_ncrna_child_of_gene)
        fa = Fasta(str(fasta_standard))
        gene_only, _ = extract_sequence.extract_features_to_fasta(
            ref_db, ["gene"], fa, str(tmp_path / "g"))
        superset, _ = extract_sequence.extract_features_to_fasta(
            ref_db, ["gene", "ncRNA"], fa, str(tmp_path / "gs"))
        assert open(gene_only).read() == open(superset).read()

    def test_legacy_dict_path_no_double_process(
            self, gff_ncrna_child_of_gene, fasta_standard):
        """The legacy dict extractor must also be unaffected by listing the
        child type: same {id: seq} whether or not `ncRNA` is in the list."""
        ref_db = self._db(gff_ncrna_child_of_gene)
        fa = Fasta(str(fasta_standard))
        t_gene, _ = extract_sequence.extract_features(ref_db, ["gene"], fa)
        t_super, _ = extract_sequence.extract_features(
            ref_db, ["gene", "ncRNA"], fa)
        assert "rna1" in t_gene
        assert t_gene == t_super
