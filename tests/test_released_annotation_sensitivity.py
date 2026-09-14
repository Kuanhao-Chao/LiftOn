from __future__ import annotations

import hashlib
import json

import pytest

from benchmarks.compare import (
    release_evaluation,
    released_annotation_sensitivity as sensitivity,
    target_truth,
    whole_genome_study,
)


def _write_pairs(root):
    for pair_id in whole_genome_study.EXPECTED_PAIR_IDS:
        path = root / pair_id / "pair_result.json"
        path.parent.mkdir(parents=True)
        path.write_text(json.dumps({
            "schema_version": release_evaluation.SCHEMA_VERSION,
            "panel": "e2e",
            "benchmark": pair_id,
            "repetition": 1,
            "modes": {"candidate": "fast", "reference": "fast"},
            "provenance": {
                "candidate": {"sha": sensitivity.CANDIDATE_SHA},
                "reference": {"sha": sensitivity.REFERENCE_SHA},
            },
        }))


def test_repetition_one_pairs_requires_exact_complete_cohort(tmp_path):
    _write_pairs(tmp_path)
    result = sensitivity._repetition_one_pairs(tmp_path)
    assert set(result) == set(whole_genome_study.EXPECTED_PAIR_IDS)

    next(tmp_path.rglob("pair_result.json")).unlink()
    with pytest.raises(sensitivity.SensitivityError, match="lacks one"):
        sensitivity._repetition_one_pairs(tmp_path)


def test_output_path_is_bound_to_execution_fingerprint(tmp_path):
    output = tmp_path / "lifted.gff3"
    output.write_text("##gff-version 3\n")
    version = {
        "output_gff": str(output),
        "fingerprints": {
            "byte_sha256": hashlib.sha256(output.read_bytes()).hexdigest(),
        },
    }
    assert sensitivity._output_path(version, "candidate") == output.resolve()

    output.write_text("changed\n")
    with pytest.raises(sensitivity.SensitivityError, match="changed"):
        sensitivity._output_path(version, "candidate")


def test_primary_projection_rejects_changed_artifact(tmp_path):
    document = {
        "method": target_truth.METHOD,
        "parameters": {"minimum_reciprocal_overlap": 0.5},
        "scope": {},
        "gene": {},
        "transcript": {},
        "structure": {},
    }
    artifact = tmp_path / "target.json"
    artifact.write_text(json.dumps(document))
    version = {
        "summary": {"target_truth": document},
        "evaluation_artifacts": {
            "target_truth": whole_genome_study.file_record(artifact),
        },
    }
    assert sensitivity._primary_from_pair(version, "candidate") == (
        sensitivity._project(document)
    )

    artifact.write_text("{}")
    with pytest.raises(sensitivity.SensitivityError, match="changed"):
        sensitivity._primary_from_pair(version, "candidate")
