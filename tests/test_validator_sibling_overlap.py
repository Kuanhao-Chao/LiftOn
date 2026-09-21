"""The validator must catch what it has always claimed to catch.

`gff3_validator`'s module docstring listed "No overlapping CDS within one
transcript" from the day the file was written. No such check existed, at any
severity, and there was none for overlapping exons either -- which is why 27
transcripts of overlapping exons shipped in a published CHM13 annotation that
`gff3-validate` called clean, and why issue #26's second half went unnoticed
through several releases.

`normalize_containment` makes this worse, not better: the overlapping pair used
to share an exon ID, which the mandatory `duplicate_id` check WOULD have caught
as a fatal error, and the collision renumbering gives them distinct IDs first.
"""
import textwrap

import pytest

from lifton import gff3_validator


def _write(tmp_path, body):
    path = tmp_path / "out.gff3"
    path.write_text("##gff-version 3\n" + textwrap.dedent(body).lstrip())
    return str(path)


GENE = "chr1\tLiftOn\tgene\t100\t900\t.\t+\t.\tID=g1\n"
MRNA = "chr1\tLiftOn\tmRNA\t100\t900\t.\t+\t.\tID=tx1;Parent=g1\n"


def _codes(result):
    return sorted({i.check for i in result.issues})


class TestOverlappingExons:
    def test_an_overlapping_exon_pair_is_an_error(self, tmp_path):
        path = _write(tmp_path, GENE + MRNA + (
            "chr1\tLiftOn\texon\t100\t400\t.\t+\t.\tID=e1;Parent=tx1\n"
            "chr1\tLiftOn\texon\t300\t900\t.\t+\t.\tID=e2;Parent=tx1\n"))
        result = gff3_validator.validate_gff3_file(path)
        assert "exon_overlap" in _codes(result)
        assert not result.is_valid

    def test_distinct_ids_do_not_excuse_the_overlap(self, tmp_path):
        """The exact shape normalize_containment leaves behind."""
        path = _write(tmp_path, GENE + MRNA + (
            "chr1\tLiftOn\texon\t100\t400\t.\t+\t.\tID=exon-X-1;Parent=tx1\n"
            "chr1\tLiftOn\texon\t150\t300\t.\t+\t.\tID=exon-X-2;Parent=tx1\n"))
        result = gff3_validator.validate_gff3_file(path)
        assert "exon_overlap" in _codes(result)

    def test_abutting_exons_are_fine(self, tmp_path):
        """400 then 401 is not an overlap -- inclusive coordinates."""
        path = _write(tmp_path, GENE + MRNA + (
            "chr1\tLiftOn\texon\t100\t400\t.\t+\t.\tID=e1;Parent=tx1\n"
            "chr1\tLiftOn\texon\t401\t900\t.\t+\t.\tID=e2;Parent=tx1\n"))
        assert "exon_overlap" not in _codes(gff3_validator.validate_gff3_file(path))

    def test_ordinary_spliced_exons_are_fine(self, tmp_path):
        path = _write(tmp_path, GENE + MRNA + (
            "chr1\tLiftOn\texon\t100\t400\t.\t+\t.\tID=e1;Parent=tx1\n"
            "chr1\tLiftOn\texon\t600\t900\t.\t+\t.\tID=e2;Parent=tx1\n"))
        result = gff3_validator.validate_gff3_file(path)
        assert "exon_overlap" not in _codes(result)
        assert result.is_valid


class TestOverlappingCds:
    def test_the_docstrings_long_standing_claim_now_holds(self, tmp_path):
        path = _write(tmp_path, GENE + MRNA + (
            "chr1\tLiftOn\texon\t100\t900\t.\t+\t.\tID=e1;Parent=tx1\n"
            "chr1\tLiftOn\tCDS\t100\t400\t.\t+\t0\tID=c1;Parent=tx1\n"
            "chr1\tLiftOn\tCDS\t380\t700\t.\t+\t0\tID=c1;Parent=tx1\n"))
        assert "cds_overlap" in _codes(gff3_validator.validate_gff3_file(path))

    def test_a_duplicated_stop_codon_sized_cds_is_caught(self, tmp_path):
        """5 CHM13 transcripts carried an exactly duplicated 3 bp CDS row."""
        path = _write(tmp_path, GENE + MRNA + (
            "chr1\tLiftOn\texon\t100\t900\t.\t+\t.\tID=e1;Parent=tx1\n"
            "chr1\tLiftOn\tCDS\t100\t400\t.\t+\t0\tID=c1;Parent=tx1\n"
            "chr1\tLiftOn\tCDS\t398\t400\t.\t+\t0\tID=c1;Parent=tx1\n"))
        assert "cds_overlap" in _codes(gff3_validator.validate_gff3_file(path))

    def test_a_normal_discontinuous_cds_still_passes(self, tmp_path):
        """Segments of one multi-exon CDS share an ID and must stay valid."""
        path = _write(tmp_path, GENE + MRNA + (
            "chr1\tLiftOn\texon\t100\t400\t.\t+\t.\tID=e1;Parent=tx1\n"
            "chr1\tLiftOn\texon\t600\t900\t.\t+\t.\tID=e2;Parent=tx1\n"
            "chr1\tLiftOn\tCDS\t100\t400\t.\t+\t0\tID=c1;Parent=tx1\n"
            "chr1\tLiftOn\tCDS\t600\t900\t.\t+\t2\tID=c1;Parent=tx1\n"))
        result = gff3_validator.validate_gff3_file(path)
        assert "cds_overlap" not in _codes(result)


class TestTheDocstringIsNoLongerStale:
    def test_the_claim_matches_the_implementation(self):
        assert hasattr(gff3_validator, "_check_sibling_overlap")
        assert "No overlapping exons" in gff3_validator.__doc__
