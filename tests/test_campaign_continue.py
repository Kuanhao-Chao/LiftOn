from dataclasses import asdict
from pathlib import Path

from benchmarks.compare import build_controller, continue_campaign


def _plan(tmp_path):
    old_root = tmp_path / "original"
    cells = []
    provenance = {"files": {}}
    provenance["fingerprint"] = build_controller.canonical_hash(provenance)
    policy = build_controller.Policy()
    for benchmark in ("alpha", "beta"):
        cell_dir = old_root / "cells" / f"subset__{benchmark}"
        cell = {
            "id": f"subset__{benchmark}",
            "kind": "subset",
            "stage": "subset",
            "benchmark": benchmark,
            "threads": 8,
            "full_job": False,
            "cell_dir": str(cell_dir),
            "command": ["python", "--cell-dir", str(cell_dir), "--threads", "8"],
            "artifacts": {
                "result_json": str(cell_dir / "result.json"),
                "gff": str(cell_dir / "out.gff3"),
                "manifest": str(cell_dir / "manifest.json"),
            },
        }
        cell["fingerprint"] = build_controller._cell_fingerprint(
            cell, provenance["fingerprint"]
        )
        cells.append(cell)
    plan = {
        "schema_version": build_controller.CONTROLLER_SCHEMA_VERSION,
        "run_id": "original",
        "created_at": build_controller.utc_now(),
        "repo_root": str(tmp_path),
        "run_dir": str(old_root),
        "stage": "subset",
        "ids": ["alpha", "beta"],
        "policy": asdict(policy),
        "inputs": {
            "registry": str(tmp_path / "benchmarks.json"),
            "dataset_registry": str(tmp_path / "datasets.json"),
            "baseline": str(tmp_path / "baseline.json"),
        },
        "provenance": provenance,
        "cells": cells,
    }
    plan["fingerprint"] = build_controller._plan_fingerprint(plan)
    build_controller.validate_plan_integrity(plan)
    return plan


def test_continuation_selects_cells_and_rewrites_immutable_paths(tmp_path):
    original = _plan(tmp_path)
    policy = build_controller.Policy(
        **{
            **original["policy"],
            "max_active": 16,
            "max_full": 4,
            "max_worker_threads": 128,
        }
    )
    continuation = continue_campaign.build_continuation_plan(
        original,
        ["subset__beta"],
        run_dir=tmp_path / "continuation",
        policy=policy,
    )

    assert continuation["ids"] == ["beta"]
    assert continuation["continuation_of"]["plan_fingerprint"] == original[
        "fingerprint"
    ]
    assert [cell["id"] for cell in continuation["cells"]] == ["subset__beta"]
    cell = continuation["cells"][0]
    assert str(tmp_path / "continuation") in cell["cell_dir"]
    assert str(tmp_path / "original") not in str(cell)
    assert cell["threads"] == 8
    build_controller.validate_plan_integrity(continuation)
