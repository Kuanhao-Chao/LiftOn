from __future__ import annotations

import json

import pytest

from benchmarks.compare import recovery_difference


COLUMNS = (
    "ref_mrna_id", "tool_feature_id", "recovered", "is_coding",
    "protein_identity", "orf_valid",
)


def _write(path, rows):
    lines = ["\t".join(COLUMNS)]
    for row in rows:
        lines.append("\t".join(str(row[column]) for column in COLUMNS))
    path.write_text("\n".join(lines) + "\n")


def _row(identifier, recovered, identity, *, coding=True, orf=True):
    return {
        "ref_mrna_id": identifier,
        "tool_feature_id": f"lifted-{identifier}",
        "recovered": str(recovered).lower(),
        "is_coding": str(coding).lower(),
        "protein_identity": identity,
        "orf_valid": str(orf).lower(),
    }


def _pair(tmp_path):
    candidate = tmp_path / "candidate.transcripts.tsv"
    reference = tmp_path / "reference.transcripts.tsv"
    _write(candidate, [
        _row("t1", True, 0.99),          # both recover
        _row("t2", False, 0.0),          # reference only
        _row("t3", False, 0.0),          # reference only, weak
        _row("t4", True, 0.95),          # candidate only
        _row("nc", True, 0.0, coding=False),
    ])
    _write(reference, [
        _row("t1", True, 0.98),
        _row("t2", True, 0.80, orf=True),
        _row("t3", True, 0.20, orf=False),
        _row("t4", False, 0.0),
        _row("nc", True, 0.0, coding=False),
    ])
    return candidate, reference


def test_partitions_coding_transcripts_by_recovering_arm(tmp_path):
    candidate, reference = _pair(tmp_path)

    result = recovery_difference.compare_arms(candidate, reference)

    assert result["denominator_coding"] == 4, "non-coding rows must be excluded"
    assert result["reference_only"]["n"] == 2
    assert result["candidate_only"]["n"] == 1
    assert result["net_candidate_minus_reference"] == -1


def test_reports_the_quality_of_the_models_only_one_arm_recovered(tmp_path):
    candidate, reference = _pair(tmp_path)

    reference_only = recovery_difference.compare_arms(
        candidate, reference,
    )["reference_only"]

    # The identities are read from the arm that actually recovered them.
    assert reference_only["mean_protein_identity"] == pytest.approx(0.5)
    assert reference_only["median_protein_identity"] == pytest.approx(0.5)
    assert reference_only["n_below_threshold"] == 1
    assert reference_only["fraction_below_threshold"] == pytest.approx(0.5)
    assert reference_only["n_orf_valid"] == 1


def test_threshold_is_configurable(tmp_path):
    candidate, reference = _pair(tmp_path)

    strict = recovery_difference.compare_arms(
        candidate, reference, weak_threshold=0.9,
    )["reference_only"]

    assert strict["n_below_threshold"] == 2


def test_empty_partition_reports_no_statistics(tmp_path):
    candidate = tmp_path / "candidate.transcripts.tsv"
    reference = tmp_path / "reference.transcripts.tsv"
    _write(candidate, [_row("t1", True, 0.9)])
    _write(reference, [_row("t1", True, 0.9)])

    result = recovery_difference.compare_arms(candidate, reference)

    assert result["reference_only"]["n"] == 0
    assert result["reference_only"]["mean_protein_identity"] is None
    assert result["reference_only"]["fraction_below_threshold"] is None


def test_document_is_deterministic_and_campaign_scoped(tmp_path):
    evaluation = (
        tmp_path / "runs" / "r1" / "cells"
        / "paired_e2e__demo__repetition_01" / "pair" / "evaluation"
    )
    evaluation.mkdir(parents=True)
    candidate, reference = _pair(evaluation)
    assert candidate.name == "candidate.transcripts.tsv"
    assert reference.name == "reference.transcripts.tsv"

    first = recovery_difference.build_document(tmp_path)
    second = recovery_difference.build_document(tmp_path)

    assert first == second
    assert json.dumps(first, sort_keys=True) == json.dumps(second, sort_keys=True)
    assert set(first["transfers"]) == {"demo"}
    assert first["method"] == "paired-recovery-difference-v1"


def test_missing_evidence_fails_closed(tmp_path):
    with pytest.raises(ValueError, match="no scored repetitions"):
        recovery_difference.build_document(tmp_path)
