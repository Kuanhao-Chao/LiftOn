from __future__ import annotations

import copy
import json

import pytest

from benchmarks.compare import biology_canary, release_evaluation


def _pair(benchmark: str) -> dict:
    fingerprint = {
        "semantic_sha256": "a" * 64,
        "byte_sha256": "b" * 64,
        "feature_counts": {"gene": 1},
    }
    summary = {
        field: None for field in biology_canary.SUMMARY_FIELDS
    }
    summary["target_truth"] = {
        field: None for field in biology_canary.TARGET_FIELDS
    }
    arm = {
        "fingerprints": fingerprint,
        "summary": summary,
        "validation": {
            "is_valid": True,
            "passes": {"full_semantic": {"issue_totals": {}}},
        },
    }
    return {
        "schema_version": release_evaluation.SCHEMA_VERSION,
        "panel": "e2e",
        "benchmark": benchmark,
        "repetition": 1,
        "modes": {"candidate": "fast", "reference": "safe"},
        "provenance": {
            label: {"sha": biology_canary.EXPECTED_SHA}
            for label in ("candidate", "reference")
        },
        "versions": {
            "candidate": copy.deepcopy(arm),
            "reference": copy.deepcopy(arm),
        },
    }


def _root(tmp_path):
    for benchmark in biology_canary.EXPECTED_IDS:
        path = tmp_path / benchmark / "pair_result.json"
        path.parent.mkdir(parents=True)
        path.write_text(json.dumps(_pair(benchmark)), encoding="utf-8")
    return tmp_path


def test_canary_accepts_two_exact_semantically_equivalent_pairs(tmp_path):
    result = biology_canary.audit_canary(_root(tmp_path))

    assert result["fast_supported"] is True
    assert [row["benchmark"] for row in result["rows"]] == list(
        biology_canary.EXPECTED_IDS
    )


def test_canary_rejects_semantic_or_target_metric_drift(tmp_path):
    root = _root(tmp_path)
    path = root / "bee" / "pair_result.json"
    pair = json.loads(path.read_text())
    pair["versions"]["candidate"]["fingerprints"]["semantic_sha256"] = "c" * 64
    pair["versions"]["candidate"]["summary"]["target_truth"]["gene"] = {
        "locus": {"f1": 0.9},
    }
    path.write_text(json.dumps(pair), encoding="utf-8")

    result = biology_canary.audit_canary(root)

    assert result["fast_supported"] is False
    assert result["rows"][0]["semantic_equal"] is False
    assert result["rows"][0]["target_truth_equal"] is False


def test_canary_fails_closed_on_missing_pair_or_wrong_sha(tmp_path):
    root = _root(tmp_path)
    (root / "bee" / "pair_result.json").unlink()
    with pytest.raises(biology_canary.CanaryError, match="expected 2"):
        biology_canary.audit_canary(root)

    root = _root(tmp_path / "wrong")
    path = root / "bee" / "pair_result.json"
    pair = json.loads(path.read_text())
    pair["provenance"]["candidate"]["sha"] = "d" * 40
    path.write_text(json.dumps(pair), encoding="utf-8")
    with pytest.raises(biology_canary.CanaryError, match="exact v1.0.11"):
        biology_canary.audit_canary(root)
