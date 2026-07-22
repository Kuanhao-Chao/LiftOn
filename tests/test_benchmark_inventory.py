"""Regression tests for the deterministic curated benchmark inventory."""
from __future__ import annotations

import copy
import json
from pathlib import Path

import pytest

from benchmarks import inventory


REPO_ROOT = Path(__file__).resolve().parents[1]
RULES_PATH = REPO_ROOT / "benchmarks" / "inventory_rules.json"
JSON_OUTPUT = REPO_ROOT / "benchmarks" / "BENCHMARK_INVENTORY.json"
MARKDOWN_OUTPUT = REPO_ROOT / "benchmarks" / "BENCHMARK_INVENTORY.md"


def _build():
    return inventory.build_inventory(
        REPO_ROOT, inventory.load_rules(RULES_PATH)
    )


def test_committed_inventory_is_current_and_deterministic():
    first = _build()
    second = _build()
    assert first == second
    assert inventory.render_json(first) == JSON_OUTPUT.read_text()
    assert inventory.render_markdown(first) == MARKDOWN_OUTPUT.read_text()
    assert first["inventory_sha256"] == inventory.canonical_sha256({
        key: value for key, value in first.items()
        if key != "inventory_sha256"
    })


def test_inventory_covers_all_curated_classes_and_registry_ids():
    document = _build()
    counts = document["source_summary"]["classification_counts"]
    assert all(counts[name] > 0 for name in inventory.VALID_CLASSIFICATIONS)
    assert document["registry_summary"]["entry_count"] == 35
    assert document["registry_summary"]["duplicate_ids"] == []
    assert document["registry_summary"]["schema_incomplete_ids"] == []
    assert document["canonical_v2_dataset_summary"]["scenario_count"] == 12
    assert document["canonical_v2_dataset_summary"]["source_count"] == 18
    assert document["canonical_v2_dataset_summary"]["remote_source_count"] == 15
    assert len(document["files"]) == document["source_summary"]["file_count"]
    assert len({record["path"] for record in document["files"]}) == len(
        document["files"]
    )


def test_local_run_artifacts_are_not_release_evidence():
    document = _build()
    paths = {record["path"] for record in document["files"]}
    assert not any("/_runs/" in path for path in paths)
    assert "benchmarks/compare/cross_locus_rescue_ab.py" not in paths
    assert "benchmarks/results_rerun/FINAL_RERUN_REPORT.md" not in paths
    assert "benchmarks/compare/_runs/**" in document["ineligible_local_patterns"]


def test_incomplete_soybean_result_is_explicitly_ineligible():
    document = _build()
    record = next(
        item for item in document["files"]
        if item["path"].endswith("_full_runs/t1_soybean_w82_to_lee.json")
    )
    assert record["classification"] == "ineligible"
    assert record["completion_state"] == "incomplete"
    assert record["claim_eligible"] is False


def test_frozen_fourway_hash_is_enforced():
    rules = copy.deepcopy(inventory.load_rules(RULES_PATH))
    rules["frozen_artifacts"][
        "benchmarks/compare/fourway_results.json"
    ] = "0" * 64
    with pytest.raises(inventory.InventoryError, match="frozen artifact changed"):
        inventory.build_inventory(REPO_ROOT, rules)


def test_every_json_artifact_has_a_schema_fingerprint():
    document = json.loads(JSON_OUTPUT.read_text())
    json_records = [
        record for record in document["files"] if record["path"].endswith(".json")
    ]
    assert json_records
    assert all(record["json_schema"] is not None for record in json_records)
    assert all(
        len(record["json_schema"]["shape_sha256"]) == 64
        for record in json_records
    )
