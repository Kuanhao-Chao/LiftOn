from __future__ import annotations

import gzip
import os

import pytest
from hypothesis import HealthCheck, given, settings, strategies as st

from lifton.gff3_validator import (
    fasta_sequence_lengths,
    validate_gff3_target_bounds,
)


def _checks(result):
    return {issue.check for issue in result.errors}


def _warning_checks(result):
    return {issue.check for issue in result.warnings}


def test_target_bound_validator_accepts_declared_regions(tmp_path):
    fasta = tmp_path / "target.fa"
    fasta.write_text(">chr1 description\nACGTACGTAA\n>chr2\nNNNN\n")
    gff = tmp_path / "output.gff3"
    gff.write_text(
        "##gff-version 3\n"
        "##sequence-region chr1 1 10\n"
        "chr1\tLiftOn\tgene\t1\t10\t.\t+\t.\tID=g1\n"
    )

    result = validate_gff3_target_bounds(str(gff), str(fasta))

    assert result.is_valid
    assert result.stats == {"gene": 1}


def test_fasta_lengths_stream_plain_and_gzip(tmp_path):
    plain = tmp_path / "target.fa"
    compressed = tmp_path / "target.fa.gz"
    content = ">chr1 full header\nAC GT\nA\n>chr2\nNN\n"
    plain.write_text(content)
    with gzip.open(compressed, "wt") as handle:
        handle.write(content)

    assert fasta_sequence_lengths(str(plain)) == {"chr1": 5, "chr2": 2}
    assert fasta_sequence_lengths(str(compressed)) == {"chr1": 5, "chr2": 2}


def test_fasta_lengths_ignore_newer_index_from_different_fasta(tmp_path):
    fasta = tmp_path / "target.fa"
    fasta.write_text(">chr1\nAAAAAAAAAA\n")
    index = tmp_path / "target.fa.fai"
    # A timestamp-only check accepted this plausible five-column sidecar and
    # treated a 10-base target as 100 bases. Keep it newer to prove that content
    # verification, rather than mtime ordering, controls reuse.
    index.write_text("chr1\t100\t6\t100\t101\n")
    newer = fasta.stat().st_mtime_ns + 1_000_000
    os.utime(index, ns=(newer, newer))

    assert fasta_sequence_lengths(str(fasta)) == {"chr1": 10}

    gff = tmp_path / "output.gff3"
    gff.write_text(
        "##gff-version 3\n"
        "chr1\tLiftOn\tgene\t1\t50\t.\t+\t.\tID=g1\n"
    )
    result = validate_gff3_target_bounds(str(gff), str(fasta))
    assert not result.is_valid
    assert _checks(result) == {"target_coordinate_out_of_bounds"}


@pytest.mark.parametrize(
    ("directive", "feature", "expected"),
    [
        (
            "",
            "chrZ\tLiftOn\tgene\t1\t2\t.\t+\t.\tID=g1",
            {"target_seqid_unknown"},
        ),
        (
            "",
            "chr1\tLiftOn\tgene\t1\t11\t.\t+\t.\tID=g1",
            {"target_coordinate_out_of_bounds"},
        ),
        (
            "##sequence-region chr1 2 8\n",
            "chr1\tLiftOn\tgene\t1\t8\t.\t+\t.\tID=g1",
            set(),
        ),
        (
            "##sequence-region chr1 1 11\n",
            "chr1\tLiftOn\tgene\t1\t8\t.\t+\t.\tID=g1",
            {"sequence_region_out_of_bounds"},
        ),
    ],
)
def test_target_bound_validator_rejects_invalid_target_coordinates(
        tmp_path, directive, feature, expected):
    fasta = tmp_path / "target.fa"
    fasta.write_text(">chr1\nAAAAAAAAAA\n")
    gff = tmp_path / "output.gff3"
    gff.write_text(
        "##gff-version 3\n" + directive + feature + "\n",
    )

    result = validate_gff3_target_bounds(str(gff), str(fasta))
    assert expected <= _checks(result)
    if directive == "##sequence-region chr1 2 8\n":
        assert result.is_valid
        assert _warning_checks(result) == {"sequence_region_containment"}
        assert result.issue_totals == {"sequence_region_containment": 1}
        assert result.severity_totals == {"WARNING": 1}


def test_target_bound_validator_can_strictly_enforce_sequence_regions(tmp_path):
    fasta = tmp_path / "target.fa"
    fasta.write_text(">chr1\nAAAAAAAAAA\n")
    gff = tmp_path / "output.gff3"
    gff.write_text(
        "##gff-version 3\n"
        "##sequence-region chr1 2 8\n"
        "chr1\tLiftOn\tgene\t1\t8\t.\t+\t.\tID=g1\n"
    )

    result = validate_gff3_target_bounds(
        str(gff),
        str(fasta),
        strict_sequence_regions=True,
    )

    assert _checks(result) == {"sequence_region_containment"}


def test_target_bound_validator_checks_sequence_region_order_and_conflicts(
        tmp_path):
    fasta = tmp_path / "target.fa"
    fasta.write_text(">chr1\nAAAAAAAAAA\n")
    gff = tmp_path / "output.gff3"
    gff.write_text(
        "##gff-version 3\n"
        "##sequence-region chr1 1 10\n"
        "chr1\tLiftOn\tgene\t1\t8\t.\t+\t.\tID=g1\n"
        "##sequence-region chr1 2 9\n"
    )

    assert {
        "sequence_region_order",
        "sequence_region_conflict",
    } <= _checks(validate_gff3_target_bounds(str(gff), str(fasta)))


def test_target_bound_validator_fails_closed_on_ambiguous_fasta(tmp_path):
    fasta = tmp_path / "target.fa"
    fasta.write_text(">chr1\nA\n>chr1 duplicate\nA\n")
    gff = tmp_path / "output.gff3"
    gff.write_text(
        "##gff-version 3\n"
        "chr1\tLiftOn\tgene\t1\t1\t.\t+\t.\tID=g1\n"
    )

    result = validate_gff3_target_bounds(str(gff), str(fasta))

    assert _checks(result) == {"target_fasta_readable"}


def test_zero_issue_detail_cap_still_reports_invalid_result(tmp_path):
    fasta = tmp_path / "target.fa"
    fasta.write_text(">chr1\nAAAA\n")
    gff = tmp_path / "output.gff3"
    gff.write_text(
        "##gff-version 3\n"
        "chrZ\tLiftOn\tgene\t1\t2\t.\t+\t.\tID=g1\n"
    )

    result = validate_gff3_target_bounds(
        str(gff), str(fasta), max_issues_per_check=0,
    )

    assert result.issues == []
    assert result.errors == []
    assert result.issue_totals == {"target_seqid_unknown": 1}
    assert result.severity_totals == {"ERROR": 1}
    assert not result.is_valid


@given(
    length=st.integers(min_value=1, max_value=200),
    start=st.integers(min_value=1, max_value=200),
    width=st.integers(min_value=0, max_value=50),
)
@settings(suppress_health_check=[HealthCheck.function_scoped_fixture])
def test_target_bound_property_matches_coordinate_predicate(
        tmp_path, length, start, width):
    end = start + width
    fasta = tmp_path / "target.fa"
    fasta.write_text(">chr1\n" + "A" * length + "\n")
    gff = tmp_path / "output.gff3"
    gff.write_text(
        "##gff-version 3\n"
        f"chr1\tLiftOn\tgene\t{start}\t{end}\t.\t+\t.\tID=g1\n"
    )

    result = validate_gff3_target_bounds(str(gff), str(fasta))

    assert result.is_valid is (end <= length)
