"""Tests for the OUTPUT GFF3 validator (``lifton.gff3_validator`` — the
``gff3-validate`` console script and the benchmark's validity yardstick).

Covers the v1.0.9 false-positive fixes:
  * ``?`` is a spec-valid strand (stranded but unknown) and must not be flagged;
  * a discontinuous CDS (multiple segments sharing one ID) is exempt from the
    duplicate-ID check (same type + same Parent).
"""
from __future__ import annotations

from lifton.gff3_validator import validate_gff3_file


def _error_checks(result):
    return {i.check for i in result.errors}


def test_question_mark_strand_is_valid(tmp_path):
    # GFF3 permits '?' (stranded but unknown); it appears on real RefSeq
    # organellar features and must NOT be flagged as an error (v1.0.9 fix).
    fp = tmp_path / "qmark.gff3"
    fp.write_text(
        "##gff-version 3\n"
        "chr1\tLiftOn\tgene\t100\t900\t.\t?\t.\tID=gene1\n"
        "chr1\tLiftOn\tmRNA\t100\t900\t.\t?\t.\tID=mrna1;Parent=gene1\n"
        "chr1\tLiftOn\texon\t100\t900\t.\t?\t.\tID=exon1;Parent=mrna1\n"
    )
    assert "strand_valid" not in _error_checks(validate_gff3_file(str(fp)))


def test_bogus_strand_still_flagged(tmp_path):
    fp = tmp_path / "bogus.gff3"
    fp.write_text(
        "##gff-version 3\n"
        "chr1\tLiftOn\tgene\t100\t900\t.\tX\t.\tID=gene1\n"
    )
    assert "strand_valid" in _error_checks(validate_gff3_file(str(fp)))


def test_discontinuous_cds_shared_id_not_duplicate(tmp_path):
    # A multi-segment (discontinuous) CDS legitimately shares one ID across its
    # segments (same type + same Parent) — it must NOT be a duplicate_id error.
    fp = tmp_path / "disc_cds.gff3"
    fp.write_text(
        "##gff-version 3\n"
        "chr1\tLiftOn\tgene\t100\t900\t.\t+\t.\tID=gene1\n"
        "chr1\tLiftOn\tmRNA\t100\t900\t.\t+\t.\tID=mrna1;Parent=gene1\n"
        "chr1\tLiftOn\texon\t100\t400\t.\t+\t.\tID=exon1;Parent=mrna1\n"
        "chr1\tLiftOn\texon\t600\t900\t.\t+\t.\tID=exon2;Parent=mrna1\n"
        "chr1\tLiftOn\tCDS\t100\t400\t.\t+\t0\tID=cds1;Parent=mrna1\n"
        "chr1\tLiftOn\tCDS\t600\t900\t.\t+\t2\tID=cds1;Parent=mrna1\n"
    )
    assert "duplicate_id" not in _error_checks(validate_gff3_file(str(fp)))
