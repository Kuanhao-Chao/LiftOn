from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import pytest

from benchmarks.compare import build_controller as controller


def _write_json(path: Path, value) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value), encoding="utf-8")


def _minimal_plan(tmp_path: Path, *, cell_id: str = "gate__unit"):
    cell_dir = tmp_path / "run" / "cells" / cell_id
    cell = {
        "id": cell_id,
        "kind": "gate",
        "benchmark": "unit",
        "mode": "gate",
        "threads": 1,
        "full_job": False,
        "command": [sys.executable, "-c", "pass"],
        "environment": {},
        "artifacts": {},
        "cell_dir": str(cell_dir),
        "fingerprint": "cell-fingerprint",
    }
    plan = {
        "schema_version": controller.CONTROLLER_SCHEMA_VERSION,
        "run_id": "unit-run",
        "created_at": controller.utc_now(),
        "repo_root": str(tmp_path),
        "run_dir": str(tmp_path / "run"),
        "stage": "gates",
        "ids": ["unit"],
        "policy": controller.asdict(controller.Policy()),
        "inputs": {
            "registry": str(tmp_path / "registry.json"),
            "dataset_registry": str(tmp_path / "datasets.json"),
            "baseline": str(tmp_path / "baseline.json"),
        },
        "provenance": {"fingerprint": "provenance"},
        "fingerprint": "plan-fingerprint",
        "cells": [cell],
    }
    return tmp_path / "run", plan, cell


def test_default_policy_is_the_approved_shared_host_envelope():
    policy = controller.Policy()

    assert policy.threads_per_cell == 8
    assert policy.max_active == 4
    assert policy.max_full == 2
    assert policy.max_worker_threads == 32
    assert policy.load1_limit == 32.0
    assert policy.min_available_gib == 256.0
    assert policy.stagger_seconds == 15.0
    assert policy.poll_seconds == 30.0


@pytest.mark.parametrize(
    ("stdout", "stderr", "expected"),
    (
        ("tool 1.2.3\n", "warning: random cache /tmp/cache-123\n", "tool 1.2.3"),
        ("", "tool 4.5.6\n", "tool 4.5.6"),
    ),
)
def test_tool_probe_prefers_stdout_but_supports_stderr_only_versions(
        monkeypatch, stdout, stderr, expected):
    monkeypatch.setattr(
        controller, "_resolve_executable", lambda _candidate: sys.executable,
    )
    monkeypatch.setattr(
        controller,
        "_run_capture",
        lambda *_args, **_kwargs: subprocess.CompletedProcess(
            args=[sys.executable, "--version"], returncode=0,
            stdout=stdout, stderr=stderr,
        ),
    )

    record = controller.probe_tool("example", "example")

    assert record["version"] == expected
    assert record["version_exit_code"] == 0


def test_canonical_matrix_has_34_subsets_and_17_refreshes():
    subsets = controller.select_ids(
        "subset", baseline=controller.DEFAULT_BASELINE,
        dataset_registry=controller.DEFAULT_DATASET_REGISTRY,
    )
    full = controller.select_ids(
        "full", baseline=controller.DEFAULT_BASELINE,
        dataset_registry=controller.DEFAULT_DATASET_REGISTRY,
    )

    assert len(subsets) == 34
    assert len(full) == 17
    assert "t1_soybean_w82_to_lee" not in subsets
    assert controller.select_ids(
        "full-canary", baseline=controller.DEFAULT_BASELINE,
        dataset_registry=controller.DEFAULT_DATASET_REGISTRY,
    ) == ["arabidopsis", "human_to_zebrafish"]


def test_full_cell_redirects_refresh_result_into_run_directory(tmp_path):
    cell_dir = tmp_path / "cells" / "full__bee"
    cell = controller._full_cell("bee", cell_dir, 8)

    assert cell["command"][1].endswith("build_controller.py")
    assert "_run-refresh" in cell["command"]
    assert cell["artifacts"]["result_json"] == str(cell_dir / "bee.json")
    for name in ("result_json", "gff", "manifest"):
        assert Path(cell["artifacts"][name]).is_relative_to(cell_dir)
    assert str(controller.DEFAULT_BASELINE) not in cell["command"]


def test_subset_cell_reads_via_runner_but_writes_all_artifacts_in_cell(tmp_path):
    cell_dir = tmp_path / "cells" / "subset__human_mane"
    cell = controller._subset_cell("human_mane", cell_dir, 8)

    assert "_run-subset" in cell["command"]
    assert cell["environment"] == {}
    for name in ("result_json", "gff", "manifest"):
        assert Path(cell["artifacts"][name]).is_relative_to(cell_dir)
    assert not Path(cell["artifacts"]["gff"]).is_relative_to(controller.HERE / "work")


def test_direct_script_worker_bootstraps_repository_imports(tmp_path):
    result = subprocess.run(
        [
            sys.executable,
            str(Path(controller.__file__).resolve()),
            "_run-subset",
            "--benchmark",
            "definitely-not-a-benchmark",
            "--cell-dir",
            str(tmp_path / "cell"),
            "--threads",
            "1",
        ],
        cwd=tmp_path,
        text=True,
        capture_output=True,
        check=False,
    )

    assert result.returncode == 2
    assert "No module named 'benchmarks'" not in result.stderr
    assert "unknown benchmark id" in result.stderr


def test_plan_layout_rejects_any_output_outside_its_cell(tmp_path):
    run_dir, plan, cell = _minimal_plan(tmp_path)
    cell.update({
        "kind": "subset",
        "artifacts": {
            "result_json": str(Path(cell["cell_dir"]) / "result.json"),
            "gff": str(tmp_path / "canonical" / "stale.gff3"),
            "manifest": str(Path(cell["cell_dir"]) / "run_manifest.json"),
        },
    })

    with pytest.raises(ValueError, match="escapes isolation"):
        controller.initialize_run(run_dir, plan)


def test_refresh_runner_separates_canonical_reads_from_cell_writes(tmp_path, monkeypatch):
    from benchmarks.compare import devel_refresh

    cached = tuple(tmp_path / f"cached-{index}" for index in range(4))
    seen = {}
    monkeypatch.setattr(devel_refresh, "WORK", devel_refresh.WORK)
    monkeypatch.setattr(devel_refresh, "OUT_DIR", devel_refresh.OUT_DIR)
    monkeypatch.setattr(devel_refresh, "_cached_inputs", lambda _benchmark: cached)

    def fake_refresh(benchmark, **kwargs):
        seen["benchmark"] = benchmark
        seen["cached"] = devel_refresh._cached_inputs(benchmark)
        seen["work"] = devel_refresh.WORK
        seen["out"] = devel_refresh.OUT_DIR

    monkeypatch.setattr(devel_refresh, "run_devel_refresh", fake_refresh)
    cell_dir = tmp_path / "cell"

    assert controller._run_refresh("bee", cell_dir, 8) == 0
    assert seen["cached"] == cached
    assert seen["work"] == cell_dir / "work"
    assert seen["out"] == cell_dir


def test_e2e_runner_sanitizes_note_metadata_and_overrides_threads(tmp_path, monkeypatch):
    from benchmarks import run_benchmarks

    registry = tmp_path / "datasets.json"
    _write_json(registry, {
        "datasets": [{
            "_note": "metadata must not reach Dataset(**entry)",
            "id": "bee", "species": "bee", "reference_fa": "ref.fa",
            "target_fa": "target.fa", "reference_gff": "ref.gff",
            "target_gff": "truth.gff",
        }],
        "lifton_flags": ["--stream", "-t", "8"],
        "evaluation_flags": ["-E"],
    })
    fake_here = tmp_path / "benchmarks" / "compare"
    data_dir = fake_here.parent / "data" / "bee"
    data_dir.mkdir(parents=True)
    for name in ("ref.fa", "target.fa", "ref.gff", "truth.gff"):
        (data_dir / name).write_bytes(b"x" * 2048)
    captured = {}

    def fake_run_dataset(_dataset, **kwargs):
        captured.update(kwargs)
        return {
            "dataset": "bee",
            "lift_profile": {"exit_code": 0},
            "eval_profile": {"exit_code": 0},
        }

    monkeypatch.setattr(controller, "HERE", fake_here)
    monkeypatch.setattr(run_benchmarks, "run_dataset", fake_run_dataset)
    monkeypatch.setattr(run_benchmarks, "_ensure_runtime", lambda: {"lifton": "fake"})
    cell_dir = tmp_path / "cell"

    assert controller._run_e2e("bee", cell_dir, registry, 6) == 0
    sanitized = json.loads((cell_dir / "datasets.sanitized.json").read_text())
    assert "_note" not in sanitized["datasets"][0]
    assert sanitized["lifton_flags"] == ["--stream", "-t", "6"]
    assert captured["do_download"] is False
    assert captured["force"] is True


def test_launch_policy_enforces_load_memory_threads_and_full_caps():
    policy = controller.Policy()
    subset = {"threads": 8, "full_job": False}
    full = {"threads": 8, "full_job": True}
    healthy = {"load1": 12.0, "available_gib": 400.0}

    assert controller.launch_allowed(subset, [], healthy, policy) == (True, "ok")
    assert controller.launch_allowed(
        subset, [], {"load1": 32.0, "available_gib": 400.0}, policy,
    )[0] is False
    assert controller.launch_allowed(
        subset, [], {"load1": 10.0, "available_gib": 255.9}, policy,
    )[0] is False
    assert controller.launch_allowed(subset, [subset] * 4, healthy, policy)[0] is False
    assert controller.launch_allowed(full, [full, full], healthy, policy)[0] is False
    assert controller.launch_allowed(
        subset, [{"threads": 32, "full_job": False}], healthy, policy,
    )[0] is False


def test_isolated_performance_retry_waits_for_an_empty_host_slot():
    policy = controller.Policy()
    retry = {"threads": 8, "full_job": False, "isolated_retry": True}
    active = [{"threads": 8, "full_job": False, "isolated_retry": False}]
    healthy = {"load1": 1.0, "available_gib": 900.0}

    assert controller.launch_allowed(retry, active, healthy, policy)[0] is False
    assert controller.launch_allowed(retry, [], healthy, policy)[0] is True


def test_structural_gff_check_rejects_malformed_records(tmp_path):
    valid = tmp_path / "valid.gff3"
    valid.write_text(
        "##gff-version 3\n"
        "chr1\tLiftOn\tgene\t1\t12\t.\t+\t.\tID=gene-1\n",
        encoding="utf-8",
    )
    malformed = tmp_path / "bad.gff3"
    malformed.write_text(
        "chr1\tLiftOn\tCDS\t12\t1\t.\twrong\t.\tID=cds-1\n",
        encoding="utf-8",
    )

    assert controller.validate_gff3_structure(valid) == []
    errors = controller.validate_gff3_structure(malformed)
    assert any("gff-version" in error for error in errors)
    assert any("coordinates" in error for error in errors)
    assert any("strand" in error for error in errors)
    assert any("CDS phase" in error for error in errors)


def test_artifact_gate_requires_result_manifest_and_full_validator(tmp_path, monkeypatch):
    cell_dir = tmp_path / "cell"
    result = cell_dir / "result.json"
    manifest = cell_dir / "run_manifest.json"
    gff = cell_dir / "out.gff3"
    _write_json(result, {
        "subset:demo": {
            "benchmark": "demo", "mode": "subset", "tools": ["lifton_devel"],
            "validity": {"lifton_devel": {"valid": True}},
            "wall_s": {"lifton_devel": 1.0},
            "peak_rss_mb": {"lifton_devel": 2.0},
            "completeness_coding": {"lifton_devel": 0.9},
            "mean_pi": {"lifton_devel": 0.95},
        }
    })
    _write_json(manifest, {
        "schema_version": 1,
        "run": {"status": "success", "finished_at": controller.utc_now()},
    })
    gff.write_text(
        "##gff-version 3\nchr1\tLiftOn\tgene\t1\t12\t.\t+\t.\tID=gene-1\n",
        encoding="utf-8",
    )
    cell = {
        "id": "subset__demo", "kind": "subset", "benchmark": "demo",
        "cell_dir": str(cell_dir),
        "artifacts": {
            "result_json": str(result), "result_key": "subset:demo",
            "manifest": str(manifest), "gff": str(gff),
            "gff_validator": ["validator", "{gff}"],
        },
    }
    monkeypatch.setattr(
        controller, "_run_gff_validator",
        lambda _cell, _gff: ([], {"exit_code": 0, "result": {"is_valid": True}}),
    )

    errors, validation = controller.validate_artifacts(cell, 0)
    assert errors == []
    assert set(validation["artifacts"]) == {"result_json", "gff", "manifest"}

    manifest.unlink()
    errors, _ = controller.validate_artifacts(cell, 0)
    assert any("manifest" in error for error in errors)


def test_e2e_result_schema_rejects_a_failed_evaluation(tmp_path):
    result = tmp_path / "result.json"
    _write_json(result, {
        "rows": [{
            "dataset": "bee",
            "lift_profile": {
                "exit_code": 0, "wall_clock_seconds": 10.0, "peak_rss_mb": 20.0,
            },
            "eval_profile": {"exit_code": 1},
        }]
    })
    cell = {"kind": "end_to_end", "benchmark": "bee"}

    errors = controller.validate_result_schema(cell, result)
    assert errors == ["end-to-end evaluation profile does not have exit_code 0"]


def test_worker_never_writes_success_when_validation_fails(tmp_path, monkeypatch):
    run_dir, plan, cell = _minimal_plan(tmp_path)
    controller.initialize_run(run_dir, plan)
    monkeypatch.setattr(controller, "assert_matching_provenance", lambda _plan: None)
    monkeypatch.setattr(
        controller, "validate_artifacts", lambda _cell, _started: (["invalid output"], {}),
    )

    assert controller.execute_cell(run_dir, cell["id"]) == 1
    assert not (Path(cell["cell_dir"]) / ".success").exists()
    failed = json.loads((Path(cell["cell_dir"]) / ".failed.json").read_text())
    assert failed["errors"] == ["invalid output"]
    assert json.loads((Path(cell["cell_dir"]) / "exit.json").read_text())["returncode"] == 0


def test_worker_writes_atomic_success_after_all_checks(tmp_path, monkeypatch):
    run_dir, plan, cell = _minimal_plan(tmp_path)
    controller.initialize_run(run_dir, plan)
    monkeypatch.setattr(controller, "assert_matching_provenance", lambda _plan: None)
    monkeypatch.setattr(
        controller, "validate_artifacts", lambda _cell, _started: ([], {"checked": True}),
    )

    assert controller.execute_cell(run_dir, cell["id"]) == 0
    success = json.loads((Path(cell["cell_dir"]) / ".success").read_text())
    assert success["fingerprint"] == cell["fingerprint"]
    assert success["validation"] == {"checked": True}
    assert not list(Path(cell["cell_dir"]).glob(".*.tmp"))


def test_performance_regression_gets_one_isolated_retry_then_hard_fails(
        tmp_path, monkeypatch):
    run_dir, plan, cell = _minimal_plan(tmp_path, cell_id="subset__unit")
    cell_dir = Path(cell["cell_dir"])
    cell.update({
        "kind": "subset", "mode": "subset", "benchmark": "unit",
        "artifacts": {
            "result_key": "subset:unit",
            "result_json": str(cell_dir / "result.json"),
            "gff": str(cell_dir / "out.gff3"),
            "manifest": str(cell_dir / "run_manifest.json"),
        },
    })
    controller.initialize_run(run_dir, plan)
    regression = [{
        "metric": "wall_s", "baseline": 100.0, "current": 126.0,
        "ratio": 1.26, "threshold": 1.25,
    }]
    monkeypatch.setattr(controller, "assert_matching_provenance", lambda _plan: None)
    monkeypatch.setattr(
        controller, "validate_artifacts", lambda _cell, _started: ([], {"checked": True}),
    )
    monkeypatch.setattr(
        controller, "performance_regressions", lambda _cell, _baseline: regression,
    )

    assert controller.execute_cell(run_dir, cell["id"]) == 75
    status = json.loads((cell_dir / "status.json").read_text())
    assert status["state"] == "retry_pending"
    assert status["isolated_retry"] is True
    assert not (cell_dir / ".success").exists()
    assert not (cell_dir / ".failed.json").exists()

    assert controller.execute_cell(run_dir, cell["id"]) == 1
    failed = json.loads((cell_dir / ".failed.json").read_text())
    assert "isolated rerun" in failed["errors"][0]


def test_failed_status_replaces_retry_evidence_with_latest_attempt(tmp_path):
    run_dir, plan, cell = _minimal_plan(tmp_path)
    controller.initialize_run(run_dir, plan)
    controller._write_status(
        cell, "retry_pending", attempts=1,
        validation={"attempt": 1},
        performance=[{"current": 126.0}],
    )

    controller._mark_failed(
        cell, ["isolated retry failed"], returncode=0, attempt=2,
        validation={"attempt": 2},
        performance=[{"current": 130.0}],
    )

    status = json.loads((Path(cell["cell_dir"]) / "status.json").read_text())
    assert status["validation"] == {"attempt": 2}
    assert status["performance"] == [{"current": 130.0}]


def test_performance_gate_compares_wall_and_rss_to_canonical_cell(tmp_path):
    baseline = tmp_path / "baseline.json"
    result = tmp_path / "result.json"
    _write_json(baseline, {
        "subset:demo": {
            "wall_s": {"lifton_devel": 100.0},
            "peak_rss_mb": {"lifton_devel": 200.0},
        }
    })
    _write_json(result, {
        "subset:demo": {
            "wall_s": {"lifton_devel": 126.0},
            "peak_rss_mb": {"lifton_devel": 251.0},
        }
    })
    cell = {
        "kind": "subset",
        "artifacts": {"result_key": "subset:demo", "result_json": str(result)},
    }

    regressions = controller.performance_regressions(cell, baseline)
    assert {item["metric"] for item in regressions} == {"wall_s", "peak_rss_mb"}
    assert all(item["ratio"] > 1.25 for item in regressions)


def test_subset_wall_gate_normalizes_to_paired_stable_control(tmp_path):
    baseline = tmp_path / "baseline.json"
    result = tmp_path / "result.json"
    _write_json(baseline, {
        "subset:demo": {
            "wall_s": {"lifton_stable": 200.0, "lifton_devel": 100.0},
            "peak_rss_mb": {"lifton_devel": 200.0},
        }
    })
    current = {
        "subset:demo": {
            "wall_s": {"lifton_stable": 240.0, "lifton_devel": 126.0},
            "peak_rss_mb": {"lifton_devel": 200.0},
        }
    }
    _write_json(result, current)
    cell = {
        "kind": "subset",
        "artifacts": {"result_key": "subset:demo", "result_json": str(result)},
    }

    assert controller.performance_regressions(cell, baseline) == []

    current["subset:demo"]["wall_s"]["lifton_devel"] = 160.0
    _write_json(result, current)
    regressions = controller.performance_regressions(cell, baseline)

    assert len(regressions) == 1
    regression = regressions[0]
    assert regression["metric"] == "wall_s"
    assert regression["comparison"] == "paired_stable_normalized"
    assert regression["raw_ratio"] == pytest.approx(1.6)
    assert regression["control"]["ratio"] == pytest.approx(1.2)
    assert regression["ratio"] == pytest.approx(1.6 / 1.2)


def test_resume_rejects_changed_provenance(monkeypatch):
    plan = {"provenance": {"fingerprint": "expected"}}
    monkeypatch.setattr(
        controller, "collect_current_provenance", lambda _plan: {"fingerprint": "changed"},
    )

    with pytest.raises(RuntimeError, match="start a new run"):
        controller.assert_matching_provenance(plan)
