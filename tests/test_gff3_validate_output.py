"""Tests for the OUTPUT GFF3 validator (``lifton.gff3_validator`` — the
``gff3-validate`` console script and the benchmark's validity yardstick).

Covers the v1.0.9 false-positive fixes:
  * ``?`` is a spec-valid strand (stranded but unknown) and must not be flagged;
  * a discontinuous CDS (multiple segments sharing one ID) is exempt from the
    duplicate-ID check (same type + same Parent).
"""
from __future__ import annotations

import pytest

from lifton.gff3_validator import validate_gff3_file, validate_gff3_structure


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


def test_directive_only_file_is_invalid(tmp_path):
    fp = tmp_path / "header_only.gff3"
    fp.write_text("##gff-version 3\n")

    result = validate_gff3_file(str(fp))

    assert not result.is_valid
    assert "features_present" in _error_checks(result)


def test_structural_validator_accepts_parent_first_discontinuous_cds(tmp_path):
    fp = tmp_path / "valid.gff3"
    fp.write_text(
        "##gff-version 3\n"
        "chr1\tLiftOn\tgene\t1\t30\t.\t+\t.\tID=g1\n"
        "chr1\tLiftOn\tmRNA\t1\t30\t.\t+\t.\tID=t1;Parent=g1\n"
        "chr1\tLiftOn\tCDS\t1\t9\t.\t+\t0\tID=c1;Parent=t1\n"
        "chr1\tLiftOn\tCDS\t20\t30\t.\t+\t1\tID=c1;Parent=t1\n"
    )

    assert validate_gff3_structure(str(fp)).is_valid


@pytest.mark.parametrize("root_type", ["ncRNA_gene", "mobile_genetic_element"])
def test_validators_accept_child_bearing_top_level_types(tmp_path, root_type):
    fp = tmp_path / f"{root_type}.gff3"
    fp.write_text(
        "##gff-version 3\n"
        f"chr1\tLiftOn\t{root_type}\t1\t30\t.\t+\t.\tID=root1\n"
        "chr1\tLiftOn\tncRNA\t1\t30\t.\t+\t.\tID=t1;Parent=root1\n"
        "chr1\tLiftOn\texon\t1\t30\t.\t+\t.\tID=e1;Parent=t1\n"
    )

    assert validate_gff3_structure(str(fp)).is_valid
    assert validate_gff3_file(str(fp)).is_valid


def test_validators_accept_child_bearing_parentless_transcript(tmp_path):
    fp = tmp_path / "transcript-root.gff3"
    fp.write_text(
        "##gff-version 3\n"
        "chr1\tLiftOn\tmRNA\t1\t30\t.\t+\t.\tID=t1\n"
        "chr1\tLiftOn\texon\t1\t30\t.\t+\t.\tID=e1;Parent=t1\n"
        "chr1\tLiftOn\tCDS\t1\t30\t.\t+\t0\tID=c1;Parent=t1\n"
    )

    assert validate_gff3_structure(str(fp)).is_valid
    assert validate_gff3_file(str(fp)).is_valid


def test_validators_accept_direct_exon_under_ncrna_gene(tmp_path):
    fp = tmp_path / "direct-ncrna-gene.gff3"
    fp.write_text(
        "##gff-version 3\n"
        "chr1\tLiftOn\tncRNA_gene\t1\t30\t.\t+\t.\tID=root1\n"
        "chr1\tLiftOn\texon\t1\t30\t.\t+\t.\tID=e1;Parent=root1\n"
    )

    assert validate_gff3_structure(str(fp)).is_valid
    assert validate_gff3_file(str(fp)).is_valid


def test_childless_parentless_transcript_remains_invalid(tmp_path):
    fp = tmp_path / "childless-transcript.gff3"
    fp.write_text(
        "##gff-version 3\n"
        "chr1\tLiftOn\tmRNA\t1\t30\t.\t+\t.\tID=t1\n"
    )

    assert "missing_parent" in _error_checks(validate_gff3_structure(str(fp)))
    assert "transcript_no_parent" in _error_checks(validate_gff3_file(str(fp)))


@pytest.mark.parametrize("changed", [
    "chr2\tLiftOn\tCDS\t20\t30\t.\t+\t1\tID=c1;Parent=t1",
    "chr1\tother\tCDS\t20\t30\t.\t+\t1\tID=c1;Parent=t1",
    "chr1\tLiftOn\tCDS\t20\t30\t.\t-\t1\tID=c1;Parent=t1",
    "chr1\tLiftOn\tCDS\t20\t30\t.\t+\t1\tID=c1;Parent=other",
    "chr1\tLiftOn\tCDS\t20\t30\t.\t+\t1\tID=c1",
])
def test_repeated_cds_id_requires_one_complete_signature(tmp_path, changed):
    fp = tmp_path / "collision.gff3"
    fp.write_text(
        "##gff-version 3\n"
        "chr1\tLiftOn\tgene\t1\t30\t.\t+\t.\tID=g1\n"
        "chr1\tLiftOn\tmRNA\t1\t30\t.\t+\t.\tID=t1;Parent=g1\n"
        "chr1\tLiftOn\tCDS\t1\t9\t.\t+\t0\tID=c1;Parent=t1\n"
        f"{changed}\n"
    )

    assert "duplicate_id" in _error_checks(validate_gff3_structure(str(fp)))
    assert "duplicate_id" in _error_checks(validate_gff3_file(str(fp)))


def test_structural_validator_requires_header_and_parent_first_order(tmp_path):
    fp = tmp_path / "unordered.gff3"
    fp.write_text(
        "# comment before the version\n"
        "chr1\tLiftOn\tmRNA\t1\t30\t.\t+\t.\tID=t1;Parent=g1\n"
        "chr1\tLiftOn\tgene\t1\t30\t.\t+\t.\tID=g1\n"
    )

    checks = _error_checks(validate_gff3_structure(str(fp)))
    assert {"gff3_header", "parent_not_declared"} <= checks


@pytest.mark.parametrize("feature_type,phase", [
    ("mRNA", "."),
    ("exon", "."),
    ("CDS", "0"),
])
def test_structural_validator_requires_parent_on_hierarchy_children(
        tmp_path, feature_type, phase):
    fp = tmp_path / "parentless.gff3"
    fp.write_text(
        "##gff-version 3\n"
        "chr1\tLiftOn\tgene\t1\t30\t.\t+\t.\tID=g1\n"
        f"chr1\tLiftOn\t{feature_type}\t1\t30\t.\t+\t{phase}\tID=x1\n"
    )

    assert "missing_parent" in _error_checks(
        validate_gff3_structure(str(fp))
    )


def test_structural_validator_rejects_wrong_parent_types_and_nested_gene(
        tmp_path):
    fp = tmp_path / "wrong-hierarchy.gff3"
    fp.write_text(
        "##gff-version 3\n"
        "chr1\tLiftOn\tgene\t1\t100\t.\t+\t.\tID=g1\n"
        "chr1\tLiftOn\tmisc_feature\t1\t100\t.\t+\t.\tID=p1;Parent=g1\n"
        "chr1\tLiftOn\texon\t1\t20\t.\t+\t.\tID=e1;Parent=p1\n"
        "chr1\tLiftOn\tmRNA\t1\t100\t.\t+\t.\tID=t1;Parent=p1\n"
        "chr1\tLiftOn\tCDS\t1\t20\t.\t+\t0\tID=c1;Parent=p1\n"
        "chr1\tLiftOn\tncRNA_gene\t200\t300\t.\t+\t.\tID=g2;Parent=g1\n"
    )

    checks = _error_checks(validate_gff3_structure(str(fp)))
    assert {
        "exon_parent_type",
        "transcript_parent_type",
        "cds_parent_type",
        "gene_has_parent",
    } <= checks
