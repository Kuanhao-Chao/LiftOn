from __future__ import annotations

import json
import os
import subprocess
import sys
import time
from pathlib import Path

import pytest

from benchmarks.compare import build_controller as controller


def _write_json(path: Path, value) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value), encoding="utf-8")


def _e2e_biology() -> dict:
    return {
        "biological_summary": {
            "schema_version": 1,
            "reference_features": 10,
            "mapped_features": 8,
            "lost_features": 2,
            "extra_copy_features": 1,
            "feature_completeness": 0.8,
            "emitted_transcript_records": 9,
            "mapped_transcripts_reported": 8,
            "evaluated_transcript_records": 9,
            "evaluated_coding_records": 7,
            "mean_protein_identity": 0.95,
        },
        "score_summary": {
            "format": "lifton_score_v1",
            "records": 9,
            "coding_records": 7,
            "noncoding_records": 2,
            "malformed_records": 0,
        },
        "evaluation_summary": {
            "format": "lifton_eval_v1",
            "records": 9,
            "coding_records": 7,
            "noncoding_records": 2,
            "malformed_records": 0,
        },
    }


def _paired_configuration(
    tmp_path: Path,
    *,
    panel: str = "subset",
    repetitions: int = 2,
) -> dict:
    candidate_root = tmp_path / "candidate-source"
    reference_root = tmp_path / "reference-source"
    candidate_root.mkdir(parents=True, exist_ok=True)
    reference_root.mkdir(parents=True, exist_ok=True)
    return controller.paired_configuration(
        stage=f"paired-{panel}",
        candidate_root=candidate_root,
        candidate_sha="a" * 40,
        reference_root=reference_root,
        reference_sha="b" * 40,
        repetitions=repetitions,
        lifton_executable=Path(sys.executable),
        candidate_e2e_mode="stream" if panel == "e2e" else "fast",
        reference_e2e_mode="safe" if panel == "e2e" else "fast",
    )


def _paired_result_fixture(tmp_path: Path, *, panel: str = "subset"):
    from benchmarks.compare import release_evaluation

    configuration = _paired_configuration(tmp_path, panel=panel, repetitions=2)
    input_path = tmp_path / "inputs" / "reference.gff3"
    input_path.parent.mkdir(parents=True)
    input_path.write_text(
        "##gff-version 3\n"
        "chr1\tRef\tgene\t1\t12\t.\t+\t.\tID=reference\n",
        encoding="utf-8",
    )
    input_stat = input_path.stat()
    input_record = {
        "path": str(input_path.resolve()),
        "size": input_stat.st_size,
        "sha256": controller.sha256_file(input_path),
        "mtime_ns": input_stat.st_mtime_ns,
        "ctime_ns": input_stat.st_ctime_ns,
        "st_dev": input_stat.st_dev,
        "st_ino": input_stat.st_ino,
    }
    cell = controller._paired_cell(
        "demo",
        tmp_path / "run" / "cells" / "paired",
        8,
        configuration=configuration,
        repetition=1,
        input_fingerprints={"ref_gff": input_record},
    )
    profiles = {
        "candidate": {
            "exit_code": 0,
            "wall_clock_seconds": 2.0,
            "peak_rss_mb": 10.0,
        },
        "reference": {
            "exit_code": 0,
            "wall_clock_seconds": 4.0,
            "peak_rss_mb": 20.0,
        },
    }
    versions = {}
    provenance = {}
    for label in ("candidate", "reference"):
        gff = Path(cell["artifacts"][f"{label}_gff"])
        gff.parent.mkdir(parents=True, exist_ok=True)
        gff.write_text(
            "##gff-version 3\n"
            f"chr1\tLiftOn\tgene\t1\t12\t.\t+\t.\tID={label}\n",
            encoding="utf-8",
        )
        source = {
            "label": label,
            "root": configuration[label]["root"],
            "sha": configuration[label]["sha"],
            "lifton_executable": configuration["lifton_executable"],
        }
        provenance[label] = {
            **source,
            "imported_package": str(
                Path(configuration[label]["root"]) / "lifton" / "__init__.py"
            ),
        }
        profile = profiles[label]
        summary = {
            "tool": label,
            "benchmark": "demo",
            "n_reference_total": 10,
            "completeness_all": 0.8,
            "completeness_coding": 0.8,
            "protein_identity": {"n": 8, "mean": 0.95},
            "profile": {
                "wall_clock_seconds": profile["wall_clock_seconds"],
                "peak_rss_mb": profile["peak_rss_mb"],
            },
        }
        transcript = (
            Path(cell["artifacts"]["result_json"]).parent
            / "evaluation"
            / f"{label}.transcripts.tsv"
        )
        transcript.parent.mkdir(parents=True, exist_ok=True)
        transcript.write_text(
            "ref_mrna_id\trecovered\tis_coding\tprotein_identity\n"
            f"{label}-tx\t1\t1\t0.95\n",
            encoding="utf-8",
        )
        summary["transcripts_tsv"] = str(transcript.resolve())
        versions[label] = {
            "source": source,
            "profile": profile,
            "argv": [sys.executable, "-m", "lifton"],
            "output_gff": str(gff),
            "fingerprints": release_evaluation.gff3_fingerprints(gff),
            "validation": {"is_valid": True, "n_errors": 0},
            "summary": summary,
            "evaluation_artifacts": {
                "transcripts_tsv": {
                    "path": str(transcript.resolve()),
                    "size": transcript.stat().st_size,
                    "sha256": controller.sha256_file(transcript),
                },
            },
        }
        if panel == "e2e":
            versions[label].update({
                "e2e_mode": configuration[label]["e2e_mode"],
                "lifton_flags": ["-t", "8"],
                "evaluation_flags": ["-E"],
                "evaluation_profile": {
                    "exit_code": 0,
                    "wall_clock_seconds": 1.0,
                    "peak_rss_mb": 5.0,
                },
                **_e2e_biology(),
            })
        manifest_path = Path(cell["artifacts"][f"{label}_manifest"])
        native_manifests = {
            "lift": {
                "path": str(
                    manifest_path.parent / "lifton_output" / "run_manifest.json"
                ),
                "present": False,
            },
        }
        versions[label].update({
            "release_manifest": str(manifest_path),
            "native_manifests": native_manifests,
        })
        fingerprints = versions[label]["fingerprints"]
        protocol = (
            {
                "kind": "e2e",
                "mode": configuration[label]["e2e_mode"],
                "lifton_flags": ["-t", "8"],
                "evaluation_flags": ["-E"],
            }
            if panel == "e2e"
            else {"kind": panel, "argv": [sys.executable, "-m", "lifton"]}
        )
        _write_json(manifest_path, {
            "schema_version": release_evaluation.SCHEMA_VERSION,
            "kind": "paired_release_arm",
            "run": {
                "status": "success",
                "finished_at": controller.utc_now(),
            },
            "source": source,
            "protocol": protocol,
            "profile": profile,
            "artifacts": {
                "output_gff": {
                    "path": str(gff.resolve()),
                    "size": gff.stat().st_size,
                    "byte_sha256": fingerprints["byte_sha256"],
                    "semantic_sha256": fingerprints["semantic_sha256"],
                },
                "native_manifests": native_manifests,
            },
            "validation": versions[label]["validation"],
        })
    raw = {
        "schema_version": release_evaluation.SCHEMA_VERSION,
        "panel": panel,
        "benchmark": "demo",
        "repetition": 1,
        "order": ["reference", "candidate"],
        "threads": 8,
        "modes": {
            label: (
                configuration[label]["e2e_mode"] if panel == "e2e" else panel
            )
            for label in ("candidate", "reference")
        },
        "registries": dict(configuration["registries"]),
        "provenance": provenance,
        "inputs": {
            "ref_gff": {
                key: value for key, value in input_record.items()
                if key != "mtime_ns"
            },
        },
        "versions": versions,
        "ratios": {"wall": 0.5, "peak_rss": 0.5},
    }
    result = Path(cell["artifacts"]["result_json"])
    _write_json(result, raw)
    return cell, raw


def _seal_plan(plan) -> None:
    provenance = plan["provenance"]
    provenance["fingerprint"] = controller.canonical_hash({
        key: value for key, value in provenance.items()
        if key != "fingerprint"
    })
    for cell in plan["cells"]:
        cell["fingerprint"] = controller._cell_fingerprint(
            cell, provenance["fingerprint"],
        )
    plan["fingerprint"] = controller._plan_fingerprint(plan)


def _initialize_run(run_dir: Path, plan) -> None:
    _seal_plan(plan)
    controller.initialize_run(run_dir, plan)


def _minimal_plan(tmp_path: Path, *, cell_id: str = "gate__unit",
                  run_name: str = "run"):
    cell_dir = tmp_path / run_name / "cells" / cell_id
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
    }
    plan = {
        "schema_version": controller.CONTROLLER_SCHEMA_VERSION,
        "run_id": f"unit-{run_name}",
        "created_at": controller.utc_now(),
        "repo_root": str(tmp_path),
        "run_dir": str(tmp_path / run_name),
        "stage": "gates",
        "ids": ["unit"],
        "policy": controller.asdict(controller.Policy()),
        "inputs": {
            "registry": str(tmp_path / "registry.json"),
            "dataset_registry": str(tmp_path / "datasets.json"),
            "baseline": str(tmp_path / "baseline.json"),
        },
        "provenance": {"source": "unit-test"},
        "cells": [cell],
    }
    _seal_plan(plan)
    return tmp_path / run_name, plan, cell


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
    assert policy.subset_timeout_seconds == 3 * 60 * 60
    assert policy.full_timeout_seconds == 8 * 60 * 60
    assert policy.e2e_timeout_seconds == 24 * 60 * 60
    assert policy.stall_timeout_seconds == 2 * 60 * 60
    assert policy.watchdog_poll_seconds == 30.0
    assert policy.terminate_grace_seconds == 30.0


def test_legacy_policy_documents_load_with_watchdog_defaults():
    legacy = controller.asdict(controller.Policy())
    for key in (
        "subset_timeout_seconds", "full_timeout_seconds", "e2e_timeout_seconds",
        "stall_timeout_seconds", "watchdog_poll_seconds", "terminate_grace_seconds",
    ):
        legacy.pop(key)

    policy = controller.Policy(**legacy)

    assert policy.subset_timeout_seconds == 3 * 60 * 60
    assert policy.full_timeout_seconds == 8 * 60 * 60
    assert policy.e2e_timeout_seconds == 24 * 60 * 60


def test_watchdog_uses_stage_specific_hard_limits():
    policy = controller.Policy(
        subset_timeout_seconds=11,
        full_timeout_seconds=22,
        e2e_timeout_seconds=33,
    )

    assert controller._hard_timeout_seconds({"kind": "gate"}, policy) == 11
    assert controller._hard_timeout_seconds({"kind": "subset"}, policy) == 11
    assert controller._hard_timeout_seconds({"kind": "full_refresh"}, policy) == 22
    assert controller._hard_timeout_seconds({"kind": "end_to_end"}, policy) == 33
    assert controller._hard_timeout_seconds(
        {"kind": "paired_release", "panel": "subset"}, policy,
    ) == 22
    assert controller._hard_timeout_seconds(
        {"kind": "paired_release", "panel": "full"}, policy,
    ) == 44
    assert controller._hard_timeout_seconds(
        {"kind": "paired_release", "panel": "e2e"}, policy,
    ) == 66


def test_policy_rejects_nonpositive_watchdog_limits():
    with pytest.raises(ValueError, match="stall_timeout_seconds"):
        controller.Policy(stall_timeout_seconds=0).validate()


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


def test_runtime_dependency_manifest_is_deterministic_and_complete(monkeypatch):
    class FakeDistribution:
        version = "1.2.3"
        metadata = {"Name": "Canonical-Name"}

        @staticmethod
        def read_text(filename):
            return "metadata\n" if filename == "METADATA" else None

    requested = []
    monkeypatch.setattr(
        controller.importlib_metadata,
        "distribution",
        lambda name: requested.append(name) or FakeDistribution(),
    )

    first = controller.collect_runtime_dependencies()
    second = controller.collect_runtime_dependencies()

    expected = sorted(
        set(controller.REQUIRED_RUNTIME_DISTRIBUTIONS)
        | set(controller.OPTIONAL_RUNTIME_DISTRIBUTIONS)
    )
    assert requested == expected + expected
    assert first == second
    assert first["python"]["version"] == ".".join(
        str(part) for part in sys.version_info[:3]
    )
    assert sorted(first["distributions"]) == expected
    assert first["distributions"]["duckdb"]["version"] == "1.2.3"
    assert len(
        first["distributions"]["duckdb"]["metadata"]["METADATA"]["sha256"]
    ) == 64


def test_runtime_dependency_manifest_fails_closed_when_required_is_missing(
        monkeypatch):
    def missing_distribution(name):
        raise controller.importlib_metadata.PackageNotFoundError(name)

    monkeypatch.setattr(
        controller.importlib_metadata, "distribution", missing_distribution,
    )

    with pytest.raises(RuntimeError, match="required.*duckdb"):
        controller.collect_runtime_dependencies()


def test_runtime_dependency_drift_invalidates_resume_provenance(
        tmp_path, monkeypatch):
    registry = tmp_path / "benchmarks.json"
    datasets = tmp_path / "datasets.json"
    baseline = tmp_path / "baseline.json"
    for path in (registry, datasets, baseline):
        path.write_text("{}\n", encoding="utf-8")
    runtime = {"python": {"version": "3.11.15"}, "distributions": {"numpy": "1"}}
    monkeypatch.setattr(
        controller, "collect_git_state", lambda _root: {"head": "fixed"},
    )
    monkeypatch.setattr(controller, "_registry_tools", lambda _registry: {})
    monkeypatch.setattr(
        controller, "collect_runtime_dependencies", lambda: runtime,
    )
    expected = controller.collect_provenance(
        repo_root=tmp_path,
        registry=registry,
        dataset_registry=datasets,
        baseline=baseline,
    )
    plan = {
        "repo_root": str(tmp_path),
        "inputs": {
            "registry": str(registry),
            "dataset_registry": str(datasets),
            "baseline": str(baseline),
        },
        "provenance": expected,
    }

    runtime["distributions"]["numpy"] = "2"

    with pytest.raises(RuntimeError, match="start a new run"):
        controller.assert_matching_provenance(plan)


def test_tooling_source_hashes_are_frozen_and_drift_invalidates_resume(
        tmp_path, monkeypatch):
    registry = tmp_path / "benchmarks.json"
    datasets = tmp_path / "datasets.json"
    baseline = tmp_path / "baseline.json"
    tooling = tmp_path / "release_report.py"
    for path in (registry, datasets, baseline):
        path.write_text("{}\n", encoding="utf-8")
    tooling.write_text("alpha\n", encoding="utf-8")
    monkeypatch.setattr(
        controller,
        "PROVENANCE_TOOLING_FILES",
        {"tooling_release_report": tooling},
    )
    monkeypatch.setattr(
        controller, "collect_git_state", lambda _root: {"head": "fixed"},
    )
    monkeypatch.setattr(controller, "_registry_tools", lambda _registry: {})
    monkeypatch.setattr(
        controller,
        "collect_runtime_dependencies",
        lambda: {"python": {"version": "fixed"}, "distributions": {}},
    )
    expected = controller.collect_provenance(
        repo_root=tmp_path,
        registry=registry,
        dataset_registry=datasets,
        baseline=baseline,
    )
    record = expected["files"]["tooling_release_report"]
    assert record == {
        "path": str(tooling.resolve()),
        "size": tooling.stat().st_size,
        "sha256": controller.sha256_file(tooling),
    }
    plan = {
        "repo_root": str(tmp_path),
        "inputs": {
            "registry": str(registry),
            "dataset_registry": str(datasets),
            "baseline": str(baseline),
        },
        "provenance": expected,
    }

    tooling.write_text("bravo\n", encoding="utf-8")

    current = controller.collect_current_provenance(plan)
    assert current["fingerprint"] != expected["fingerprint"]
    with pytest.raises(RuntimeError, match="start a new run"):
        controller.assert_matching_provenance(plan)


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


def test_paired_configuration_requires_exact_sources_and_supports_diagnostic_modes(
        tmp_path):
    configuration = controller.paired_configuration(
        stage="paired-e2e-canary",
        candidate_root=tmp_path / "candidate",
        candidate_sha="A" * 40,
        reference_root=tmp_path / "reference",
        reference_sha="b" * 40,
        repetitions=5,
        lifton_executable=Path(sys.executable),
        candidate_e2e_mode="stream-native",
        reference_e2e_mode="inmemory",
    )

    assert configuration["panel"] == "e2e"
    assert configuration["candidate"]["sha"] == "a" * 40
    assert configuration["candidate"]["e2e_mode"] == "stream-native"
    assert configuration["reference"]["e2e_mode"] == "inmemory"
    with pytest.raises(ValueError, match="exact 40-character"):
        controller.paired_configuration(
            stage="paired-e2e",
            candidate_root=tmp_path,
            candidate_sha="short",
            reference_root=tmp_path,
            reference_sha="b" * 40,
            repetitions=4,
            lifton_executable=Path(sys.executable),
            candidate_e2e_mode="fast",
            reference_e2e_mode="safe",
        )
    with pytest.raises(ValueError, match="even repetition count"):
        controller.paired_configuration(
            stage="paired-e2e",
            candidate_root=tmp_path,
            candidate_sha="a" * 40,
            reference_root=tmp_path,
            reference_sha="b" * 40,
            repetitions=3,
            lifton_executable=Path(sys.executable),
            candidate_e2e_mode="fast",
            reference_e2e_mode="safe",
        )
    diagnostic = controller.paired_configuration(
        stage="paired-e2e-canary",
        candidate_root=tmp_path,
        candidate_sha="a" * 40,
        reference_root=tmp_path,
        reference_sha="b" * 40,
        repetitions=3,
        lifton_executable=Path(sys.executable),
        candidate_e2e_mode="fast",
        reference_e2e_mode="safe",
    )
    assert diagnostic["repetitions"] == 3
    with pytest.raises(ValueError, match="between 1 and 5"):
        controller.paired_configuration(
            stage="paired-e2e",
            candidate_root=tmp_path,
            candidate_sha="a" * 40,
            reference_root=tmp_path,
            reference_sha="b" * 40,
            repetitions=6,
            lifton_executable=Path(sys.executable),
            candidate_e2e_mode="fast",
            reference_e2e_mode="safe",
        )


def test_controller_cli_exposes_exact_paired_sources_repetitions_and_modes(tmp_path):
    benchmark_registry = tmp_path / "benchmarks.json"
    dataset_registry = tmp_path / "datasets.json"
    args = controller.build_parser().parse_args([
        "start",
        "--stage", "paired-e2e-canary",
        "--candidate-root", str(tmp_path / "candidate"),
        "--candidate-sha", "a" * 40,
        "--reference-root", str(tmp_path / "reference"),
        "--reference-sha", "b" * 40,
        "--paired-repetitions", "4",
        "--paired-lifton-executable", sys.executable,
        "--registry", str(benchmark_registry),
        "--dataset-registry", str(dataset_registry),
        "--candidate-e2e-mode", "stream-inmemory",
        "--reference-e2e-mode", "native",
        "--dry-run",
    ])
    configuration = controller._paired_configuration_from_args(args, args.stage)

    assert args.stage == "paired-e2e-canary"
    assert configuration["repetitions"] == 4
    assert configuration["candidate"]["e2e_mode"] == "stream-inmemory"
    assert configuration["reference"]["e2e_mode"] == "native"
    assert configuration["candidate"]["sha"] == "a" * 40
    assert configuration["reference"]["sha"] == "b" * 40
    assert configuration["registries"] == {
        "benchmark": str(benchmark_registry.resolve()),
        "dataset": str(dataset_registry.resolve()),
    }


def test_noncanary_paired_cli_defaults_to_four_balanced_repetitions(tmp_path):
    args = controller.build_parser().parse_args([
        "start",
        "--stage", "paired-subset",
        "--candidate-root", str(tmp_path / "candidate"),
        "--candidate-sha", "a" * 40,
        "--reference-root", str(tmp_path / "reference"),
        "--reference-sha", "b" * 40,
        "--paired-lifton-executable", sys.executable,
        "--dry-run",
    ])

    configuration = controller._paired_configuration_from_args(args, args.stage)

    assert configuration["repetitions"] == 4


def test_paired_plan_is_immutable_resumable_and_one_cell_per_repetition(
        tmp_path, monkeypatch):
    configuration = _paired_configuration(
        tmp_path, panel="subset", repetitions=4,
    )
    benchmark_ids = ["human_mane", "human_to_zebrafish"]
    prepared = {
        "configuration": configuration,
        "sources": {
            "candidate": {"sha": "a" * 40},
            "reference": {"sha": "b" * 40},
        },
        "inputs": {
            benchmark: {
                "ref_gff": {
                    "path": str(tmp_path / f"{benchmark}.gff3"),
                    "size": 100,
                    "sha256": "c" * 64,
                    "mtime_ns": 123,
                    "ctime_ns": 124,
                    "st_dev": 1,
                    "st_ino": 2,
                },
            }
            for benchmark in benchmark_ids
        },
    }
    monkeypatch.setattr(
        controller,
        "collect_provenance",
        lambda **_kwargs: {"base": "fixed", "fingerprint": "d" * 64},
    )
    monkeypatch.setattr(
        controller,
        "_prepare_paired_provenance",
        lambda paired, ids: (
            prepared
            if paired == configuration and list(ids) == benchmark_ids
            else pytest.fail("unexpected paired plan inputs")
        ),
    )

    run_dir, plan = controller.create_plan(
        run_id="paired-unit",
        stage="paired-subset",
        requested_ids=benchmark_ids,
        runs_root=tmp_path / "runs",
        repo_root=tmp_path,
        registry=controller.DEFAULT_REGISTRY,
        dataset_registry=controller.DEFAULT_DATASET_REGISTRY,
        baseline=controller.DEFAULT_BASELINE,
        policy=controller.Policy(),
        paired=configuration,
    )

    assert plan["paired"] == configuration
    assert plan["provenance"]["paired"] == prepared
    assert plan["ids"] == benchmark_ids
    assert len(plan["cells"]) == 8
    assert {cell["repetition"] for cell in plan["cells"]} == {1, 2, 3, 4}
    assert all(cell["exclusive"] is True for cell in plan["cells"])
    assert all(cell["kind"] == "paired_release" for cell in plan["cells"])
    assert all(
        cell["command"][cell["command"].index("--benchmark-registry") + 1]
        == str(controller.DEFAULT_REGISTRY.resolve())
        for cell in plan["cells"]
    )
    assert all(
        cell["command"][cell["command"].index("--dataset-registry") + 1]
        == str(controller.DEFAULT_DATASET_REGISTRY.resolve())
        for cell in plan["cells"]
    )
    assert [cell["expected_order"] for cell in plan["cells"][:4]] == [
        ["reference", "candidate"],
        ["candidate", "reference"],
        ["reference", "candidate"],
        ["candidate", "reference"],
    ]
    assert len({cell["fingerprint"] for cell in plan["cells"]}) == 8
    assert all(
        Path(cell["artifacts"]["result_json"]).is_relative_to(
            Path(cell["cell_dir"])
        )
        for cell in plan["cells"]
    )
    assert run_dir == (tmp_path / "runs" / "paired-unit").resolve()
    controller.validate_plan_layout(plan)


def test_paired_cell_fingerprint_covers_every_release_execution_contract(
        tmp_path):
    cell, _raw = _paired_result_fixture(tmp_path)
    provenance = "a" * 64
    baseline = controller._cell_fingerprint(cell, provenance)
    mutations = [
        lambda item: item["command"].append("--different"),
        lambda item: item["paired"]["candidate"].update({"sha": "c" * 40}),
        lambda item: item["input_fingerprints"]["ref_gff"].update({"size": 999}),
        lambda item: item.update({"expected_order": ["candidate", "reference"]}),
        lambda item: item.update({"panel": "full"}),
        lambda item: item.update({"repetition": 2}),
        lambda item: item["environment"].update({"EXTRA": "1"}),
        lambda item: item["artifacts"].update({"result_key": "different"}),
    ]

    for mutate in mutations:
        changed = json.loads(json.dumps(cell))
        mutate(changed)
        assert controller._cell_fingerprint(changed, provenance) != baseline


def test_plan_integrity_cross_checks_paired_config_inputs_and_order(tmp_path):
    cell, _raw = _paired_result_fixture(tmp_path)
    plan = {
        "schema_version": controller.CONTROLLER_SCHEMA_VERSION,
        "run_id": "paired-integrity",
        "created_at": controller.utc_now(),
        "repo_root": str(tmp_path),
        "run_dir": str(tmp_path / "run"),
        "stage": "paired-subset",
        "ids": ["demo"],
        "policy": controller.asdict(controller.Policy()),
        "inputs": {},
        "paired": cell["paired"],
        "provenance": {
            "source": "unit-test",
            "paired": {
                "configuration": cell["paired"],
                "sources": {},
                "inputs": {"demo": cell["input_fingerprints"]},
            },
        },
        "cells": [cell],
    }
    _seal_plan(plan)
    controller.validate_plan_integrity(plan)

    changed = json.loads(json.dumps(plan))
    changed["paired"]["candidate"]["sha"] = "c" * 40
    _seal_plan(changed)
    with pytest.raises(ValueError, match="frozen provenance"):
        controller.validate_plan_integrity(changed)

    changed = json.loads(json.dumps(plan))
    changed["cells"][0]["input_fingerprints"]["ref_gff"]["size"] += 1
    _seal_plan(changed)
    with pytest.raises(ValueError, match="input fingerprints"):
        controller.validate_plan_integrity(changed)

    changed = json.loads(json.dumps(plan))
    changed["cells"][0]["expected_order"] = ["candidate", "reference"]
    _seal_plan(changed)
    with pytest.raises(ValueError, match="expected order"):
        controller.validate_plan_integrity(changed)

    changed = json.loads(json.dumps(plan))
    changed["policy"]["threads_per_cell"] += 1
    with pytest.raises(ValueError, match="plan fingerprint"):
        controller.validate_plan_integrity(changed)

    changed = json.loads(json.dumps(plan))
    changed["provenance"]["source"] = "tampered"
    with pytest.raises(ValueError, match="provenance fingerprint"):
        controller.validate_plan_integrity(changed)


def test_all_plan_read_paths_reject_tampered_immutable_content(tmp_path):
    run_dir, plan, cell = _minimal_plan(tmp_path)
    _initialize_run(run_dir, plan)
    plan_path = run_dir / "plan.json"
    tampered = json.loads(plan_path.read_text())
    tampered["cells"][0]["command"].append("--tampered")
    _write_json(plan_path, tampered)

    with pytest.raises(ValueError, match="cell fingerprint"):
        controller.load_plan(run_dir)
    with pytest.raises(ValueError, match="cell fingerprint"):
        controller.execute_cell(run_dir, cell["id"])
    with pytest.raises(ValueError, match="cell fingerprint"):
        controller.reconcile_run(run_dir)
    with pytest.raises(SystemExit) as error:
        controller.main(["status", "--run-dir", str(run_dir), "--json"])
    assert error.value.code == controller.STATUS_EXIT_INVALID


def test_schema_one_plan_is_explicitly_rejected(tmp_path):
    run_dir, plan, _cell = _minimal_plan(tmp_path)
    _initialize_run(run_dir, plan)
    plan_path = run_dir / "plan.json"
    legacy = json.loads(plan_path.read_text())
    legacy["schema_version"] = 1
    _write_json(plan_path, legacy)

    with pytest.raises(ValueError, match="unsupported.*plan schema"):
        controller.load_plan(run_dir)


def test_paired_plan_fails_closed_when_plan_time_inputs_are_missing(
        tmp_path, monkeypatch):
    configuration = _paired_configuration(tmp_path)
    monkeypatch.setattr(
        controller,
        "collect_provenance",
        lambda **_kwargs: {"fingerprint": "d" * 64},
    )
    monkeypatch.setattr(
        controller,
        "_prepare_paired_provenance",
        lambda *_args: (_ for _ in ()).throw(
            RuntimeError("missing or empty paired input: ref_gff")
        ),
    )

    with pytest.raises(RuntimeError, match="missing or empty paired input"):
        controller.create_plan(
            run_id="paired-missing",
            stage="paired-subset",
            requested_ids=["human_mane"],
            runs_root=tmp_path / "runs",
            repo_root=tmp_path,
            registry=controller.DEFAULT_REGISTRY,
            dataset_registry=controller.DEFAULT_DATASET_REGISTRY,
            baseline=controller.DEFAULT_BASELINE,
            policy=controller.Policy(),
            paired=configuration,
        )


def test_paired_resume_rechecks_source_and_input_provenance(tmp_path, monkeypatch):
    configuration = _paired_configuration(tmp_path)
    base = {"git": {"head": "fixed"}, "fingerprint": "ignored"}
    paired = {
        "configuration": configuration,
        "sources": {
            "candidate": {"sha": "a" * 40},
            "reference": {"sha": "b" * 40},
        },
        "inputs": {
            "human_mane": {
                "ref_gff": {
                    "path": str(tmp_path / "ref.gff3"),
                    "size": 10,
                    "mtime_ns": 20,
                    "sha256": "c" * 64,
                },
            },
        },
    }
    expected = controller._add_paired_provenance(base, paired)
    plan = {
        "repo_root": str(tmp_path),
        "inputs": {
            "registry": str(tmp_path / "registry.json"),
            "dataset_registry": str(tmp_path / "datasets.json"),
            "baseline": str(tmp_path / "baseline.json"),
        },
        "paired": configuration,
        "provenance": expected,
    }
    monkeypatch.setattr(controller, "collect_provenance", lambda **_kwargs: base)
    monkeypatch.setattr(controller, "_current_paired_provenance", lambda _plan: paired)

    assert controller.collect_current_provenance(plan)["fingerprint"] == (
        expected["fingerprint"]
    )
    changed = json.loads(json.dumps(paired))
    changed["inputs"]["human_mane"]["ref_gff"]["sha256"] = "d" * 64
    monkeypatch.setattr(
        controller, "_current_paired_provenance", lambda _plan: changed,
    )
    with pytest.raises(RuntimeError, match="start a new run"):
        controller.assert_matching_provenance(plan)


def test_paired_retry_archives_prior_workspace_without_discarding_evidence(tmp_path):
    from benchmarks.compare import release_evaluation

    cell, original = _paired_result_fixture(tmp_path)
    pair_dir = Path(cell["artifacts"]["result_json"]).parent

    controller._prepare_attempt_workspace(cell, 2)

    archive = Path(cell["cell_dir"]) / "attempt-01.pair"
    assert not pair_dir.exists()
    assert (archive / "pair_result.json").is_file()
    assert (archive / "candidate" / "release_run_manifest.json").is_file()
    relocated = json.loads((archive / "pair_result.json").read_text())
    for label in ("candidate", "reference"):
        archived_gff = archive / label / f"{label}.gff3"
        version = relocated["versions"][label]
        assert version["output_gff"] == str(archived_gff)
        assert version["release_manifest"] == str(
            archive / label / "release_run_manifest.json"
        )
        assert version["native_manifests"]["lift"]["path"].startswith(
            str(archive / label)
        )
        assert version["summary"]["transcripts_tsv"] == str(
            archive / "evaluation" / f"{label}.transcripts.tsv"
        )
        assert version["evaluation_artifacts"]["transcripts_tsv"]["path"] == str(
            archive / "evaluation" / f"{label}.transcripts.tsv"
        )
        assert version["fingerprints"] == (
            release_evaluation.gff3_fingerprints(archived_gff)
        )
        arm = json.loads(
            (archive / label / "release_run_manifest.json").read_text()
        )
        assert arm["artifacts"]["output_gff"]["path"] == str(archived_gff)
        assert arm["artifacts"]["output_gff"]["byte_sha256"] == (
            version["fingerprints"]["byte_sha256"]
        )
    assert relocated["inputs"] == original["inputs"]
    assert relocated["provenance"] == original["provenance"]
    relocation = json.loads((archive / "relocation.json").read_text())
    assert relocation["original_root"] == str(pair_dir)
    assert relocation["archive_root"] == str(archive)
    assert relocation["gff3_bytes_modified"] is False
    assert len(relocation["documents"]) == 3
    controller._prepare_attempt_workspace(cell, 2)

    pair_dir.mkdir()
    with pytest.raises(RuntimeError, match="archive already exists"):
        controller._prepare_attempt_workspace(cell, 2)


def test_paired_input_sha_is_reused_only_when_full_stat_identity_matches(
        tmp_path, monkeypatch):
    from benchmarks.compare import release_evaluation

    configuration = _paired_configuration(tmp_path)
    input_path = tmp_path / "input.gff3"
    input_path.write_text("immutable input\n", encoding="utf-8")
    stat = input_path.stat()
    record = {
        "path": str(input_path.resolve()),
        "size": stat.st_size,
        "mtime_ns": stat.st_mtime_ns,
        "ctime_ns": stat.st_ctime_ns,
        "st_dev": stat.st_dev,
        "st_ino": stat.st_ino,
        "sha256": "a" * 64,
    }
    plan = {
        "paired": configuration,
        "provenance": {
            "paired": {"inputs": {"demo": {"ref_gff": record}}},
        },
    }
    monkeypatch.setattr(
        release_evaluation,
        "verify_source",
        lambda source: {"label": source.label, "sha": source.sha},
    )
    calls = []
    monkeypatch.setattr(
        release_evaluation,
        "sha256_file",
        lambda path: calls.append(Path(path)) or "b" * 64,
    )

    current = controller._current_paired_provenance(plan)
    assert calls == []
    assert current["inputs"]["demo"]["ref_gff"]["sha256"] == "a" * 64

    plan["provenance"]["paired"]["inputs"]["demo"]["ref_gff"]["ctime_ns"] = -1
    current = controller._current_paired_provenance(plan)
    assert calls == [input_path]
    assert current["inputs"]["demo"]["ref_gff"]["sha256"] == "b" * 64


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
        _initialize_run(run_dir, plan)


def test_paired_plan_layout_checks_both_versions_and_manifests(tmp_path):
    cell, _raw = _paired_result_fixture(tmp_path)
    plan = {
        "run_dir": str(tmp_path / "run"),
        "cells": [cell],
    }
    controller.validate_plan_layout(plan)
    cell["artifacts"]["reference_manifest"] = str(
        tmp_path / "escaped" / "run_manifest.json"
    )

    with pytest.raises(ValueError, match="reference_manifest escapes isolation"):
        controller.validate_plan_layout(plan)


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
    # target_gff is deliberately absent: candidate-side ``-E`` evaluation
    # must not require a separately published comparison annotation.
    for name in ("ref.fa", "target.fa", "ref.gff"):
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


def test_paired_cells_run_exclusively():
    policy = controller.Policy()
    paired = {
        "threads": 8, "full_job": False, "exclusive": True,
    }
    ordinary = {
        "threads": 8, "full_job": False, "exclusive": False,
    }
    healthy = {"load1": 1.0, "available_gib": 900.0}

    assert controller.launch_allowed(paired, [], healthy, policy) == (True, "ok")
    assert controller.launch_allowed(paired, [ordinary], healthy, policy)[0] is False
    assert controller.launch_allowed(ordinary, [paired], healthy, policy)[0] is False


def test_controller_lease_rejects_a_second_stage_in_the_same_runs_root(tmp_path):
    _run_one, plan_one, _cell_one = _minimal_plan(tmp_path, run_name="run-one")
    _run_two, plan_two, _cell_two = _minimal_plan(tmp_path, run_name="run-two")

    with controller._controller_lease(plan_one):
        with pytest.raises(controller.ControllerIsolationError, match="holds"):
            with controller._controller_lease(plan_two):
                pass


def test_controller_lease_rejects_foreign_live_cell_sessions(tmp_path, monkeypatch):
    run_one, plan_one, cell_one = _minimal_plan(tmp_path, run_name="run-one")
    run_two, plan_two, _cell_two = _minimal_plan(tmp_path, run_name="run-two")
    _initialize_run(run_one, plan_one)
    _initialize_run(run_two, plan_two)
    controller._write_status(
        cell_one, "running", attempts=1, started_ns=time.time_ns(),
    )
    live_session = controller.cell_session_name(plan_one, cell_one)
    monkeypatch.setattr(
        controller, "tmux_has_session", lambda name: name == live_session,
    )

    with pytest.raises(controller.ControllerIsolationError, match="foreign benchmark cells"):
        with controller._controller_lease(plan_two):
            pass


@pytest.mark.skipif(not Path("/proc").is_dir(), reason="orphan identity is Linux /proc based")
def test_controller_lease_counts_a_detached_foreign_process(tmp_path, monkeypatch):
    run_one, plan_one, cell_one = _minimal_plan(tmp_path, run_name="run-one")
    run_two, plan_two, _cell_two = _minimal_plan(tmp_path, run_name="run-two")
    _initialize_run(run_one, plan_one)
    _initialize_run(run_two, plan_two)
    process = subprocess.Popen(
        [sys.executable, "-c", "import time;time.sleep(60)"],
        start_new_session=True,
    )
    record = controller._read_proc_stat(process.pid)
    assert record is not None
    controller._write_status(
        cell_one, "running", attempts=1, started_ns=0,
        command_pid=process.pid, command_pgid=process.pid,
        command_start_ticks=record["start_ticks"],
    )
    monkeypatch.setattr(controller, "tmux_has_session", lambda _name: False)
    try:
        with pytest.raises(controller.ControllerIsolationError, match="foreign benchmark cells"):
            with controller._controller_lease(plan_two):
                pass
    finally:
        os.killpg(process.pid, 9)
        process.wait(timeout=2)


@pytest.mark.skipif(not Path("/proc").is_dir(), reason="orphan identity is Linux /proc based")
def test_orphan_cleanup_terminates_the_recorded_process_group(tmp_path, monkeypatch):
    run_dir, plan, cell = _minimal_plan(tmp_path)
    _short_watchdog(plan)
    _initialize_run(run_dir, plan)
    process = subprocess.Popen(
        [sys.executable, "-c", "import time;time.sleep(60)"],
        start_new_session=True,
    )
    record = controller._read_proc_stat(process.pid)
    assert record is not None
    controller._write_status(
        cell, "running", attempts=1, started_ns=0,
        command_pid=process.pid, command_pgid=process.pid,
        command_start_ticks=record["start_ticks"],
    )
    monkeypatch.setattr(controller, "tmux_has_session", lambda _name: False)

    controller._mark_orphans(plan, grace_seconds=0)
    process.wait(timeout=2)

    failed = json.loads((Path(cell["cell_dir"]) / ".failed.json").read_text())
    assert failed["watchdog"]["reason"] == "orphaned_worker"
    assert failed["watchdog"]["cleanup"]["identity_verified"] is True
    assert failed["watchdog"]["cleanup"]["term_sent"] is True


@pytest.mark.skipif(not Path("/proc").is_dir(), reason="orphan identity is Linux /proc based")
def test_orphan_cleanup_refuses_a_recycled_pid_identity():
    process = subprocess.Popen(
        [sys.executable, "-c", "import time;time.sleep(60)"],
        start_new_session=True,
    )
    record = controller._read_proc_stat(process.pid)
    assert record is not None
    try:
        cleanup = controller._terminate_process_group(
            pid=process.pid,
            pgid=process.pid,
            start_ticks=record["start_ticks"] + 1,
            grace_seconds=0.01,
        )

        assert cleanup["identity_verified"] is False
        assert "signal refused" in cleanup["error"]
        assert process.poll() is None
    finally:
        os.killpg(process.pid, 9)
        process.wait(timeout=2)


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
    assert set(validation["artifacts"]) == {
        "result_json", "gff", "manifest", "gff_validation",
    }
    assert all(
        len(record["sha256"]) == 64
        for record in validation["artifacts"].values()
    )

    manifest.unlink()
    errors, _ = controller.validate_artifacts(cell, 0)
    assert any("manifest" in error for error in errors)


def test_paired_result_schema_checks_sources_inputs_profiles_and_outputs(tmp_path):
    cell, raw = _paired_result_fixture(tmp_path)
    result = Path(cell["artifacts"]["result_json"])

    assert controller.validate_result_schema(cell, result) == []

    cases = []
    changed = json.loads(json.dumps(raw))
    changed["provenance"]["candidate"]["sha"] = "c" * 40
    cases.append((changed, "source SHA"))
    changed = json.loads(json.dumps(raw))
    changed["inputs"]["ref_gff"]["sha256"] = "d" * 64
    cases.append((changed, "input 'ref_gff' hash"))
    changed = json.loads(json.dumps(raw))
    changed["registries"]["benchmark"] = str(tmp_path / "other.json")
    cases.append((changed, "registries"))
    changed = json.loads(json.dumps(raw))
    changed["versions"]["candidate"]["profile"]["wall_clock_seconds"] = 0
    cases.append((changed, "wall_clock_seconds"))
    changed = json.loads(json.dumps(raw))
    changed["versions"]["reference"]["validation"]["is_valid"] = False
    cases.append((changed, "reference GFF3 validation"))
    changed = json.loads(json.dumps(raw))
    changed["versions"]["candidate"]["fingerprints"]["byte_sha256"] = "bad"
    cases.append((changed, "byte_sha256"))
    changed = json.loads(json.dumps(raw))
    changed["versions"]["candidate"]["evaluation_artifacts"][
        "transcripts_tsv"
    ]["sha256"] = "bad"
    cases.append((changed, "transcripts_tsv sha256"))
    changed = json.loads(json.dumps(raw))
    changed["versions"]["candidate"]["summary"]["transcripts_tsv"] = str(
        tmp_path / "different.transcripts.tsv"
    )
    cases.append((changed, "summary transcripts_tsv path"))
    changed = json.loads(json.dumps(raw))
    changed["ratios"]["wall"] = float("inf")
    cases.append((changed, "ratio wall"))
    changed = json.loads(json.dumps(raw))
    changed["order"] = ["candidate", "reference"]
    cases.append((changed, "AB/BA order"))

    for changed, message in cases:
        _write_json(result, changed)
        errors = controller.validate_result_schema(cell, result)
        assert errors, message
        assert any(message in error for error in errors), errors


def test_schema_one_paired_result_and_arm_manifest_are_rejected(tmp_path):
    from benchmarks.compare import release_evaluation

    cell, raw = _paired_result_fixture(tmp_path)
    result = Path(cell["artifacts"]["result_json"])
    raw["schema_version"] = 1
    _write_json(result, raw)

    errors = controller.validate_result_schema(cell, result)

    assert any(
        f"schema_version is not {release_evaluation.SCHEMA_VERSION}" in error
        for error in errors
    )

    manifest = Path(cell["artifacts"]["candidate_manifest"])
    document = json.loads(manifest.read_text())
    document["schema_version"] = 1
    _write_json(manifest, document)
    manifest_errors = controller.validate_release_arm_manifest(
        cell,
        "candidate",
        manifest,
        version=raw["versions"]["candidate"],
        observed_fingerprints=raw["versions"]["candidate"]["fingerprints"],
    )
    assert any(
        f"schema_version is not {release_evaluation.SCHEMA_VERSION}" in error
        for error in manifest_errors
    )


def test_paired_e2e_result_requires_meaningful_biology_for_both_versions(tmp_path):
    cell, raw = _paired_result_fixture(tmp_path, panel="e2e")
    result = Path(cell["artifacts"]["result_json"])

    assert controller.validate_result_schema(cell, result) == []

    raw["versions"]["candidate"]["biological_summary"]["lost_features"] = 3
    _write_json(result, raw)
    errors = controller.validate_result_schema(cell, result)
    assert any("candidate E2E biology" in error for error in errors)
    assert any("does not equal" in error for error in errors)


def test_paired_artifact_gate_validates_both_fresh_independent_outputs(
        tmp_path, monkeypatch):
    cell, _raw = _paired_result_fixture(tmp_path)
    validated = []

    def fake_validator(_cell, gff):
        validated.append(Path(gff))
        return [], {"exit_code": 0, "result": {"is_valid": True}}

    monkeypatch.setattr(controller, "_run_gff_validator", fake_validator)

    errors, validation = controller.validate_artifacts(cell, 0)

    assert errors == []
    assert set(validated) == {
        Path(cell["artifacts"]["candidate_gff"]),
        Path(cell["artifacts"]["reference_gff"]),
    }
    assert set(validation["artifacts"]) == {
        "result_json",
        "candidate_gff", "reference_gff",
        "candidate_manifest", "reference_manifest",
        "candidate_gff_validation", "reference_gff_validation",
    }
    assert all(
        len(record["sha256"]) == 64
        for record in validation["artifacts"].values()
    )
    assert set(validation["gff_fingerprints"]) == {"candidate", "reference"}
    assert set(validation["evaluation_artifacts"]) == {
        "candidate", "reference",
    }
    for label in ("candidate", "reference"):
        evidence = validation["evaluation_artifacts"][label]["transcripts_tsv"]
        assert evidence["sha256"] == _raw["versions"][label][
            "evaluation_artifacts"
        ]["transcripts_tsv"]["sha256"]
        assert evidence["size"] > 0
    assert (
        Path(cell["cell_dir"]) / "candidate_gff_validation.json"
    ).is_file()
    assert (
        Path(cell["cell_dir"]) / "reference_gff_validation.json"
    ).is_file()

    candidate_manifest = Path(cell["artifacts"]["candidate_manifest"])
    manifest_document = json.loads(candidate_manifest.read_text())
    manifest_document["profile"]["wall_clock_seconds"] = 99.0
    _write_json(candidate_manifest, manifest_document)
    errors, _ = controller.validate_artifacts(cell, 0)
    assert any(
        "candidate release manifest profile disagrees" in error
        for error in errors
    )
    manifest_document["profile"]["wall_clock_seconds"] = 2.0
    _write_json(candidate_manifest, manifest_document)

    manifest_document["protocol"]["argv"] = ["different-lifton-command"]
    _write_json(candidate_manifest, manifest_document)
    errors, _ = controller.validate_artifacts(cell, 0)
    assert any(
        "candidate release manifest argv disagrees" in error
        for error in errors
    )
    manifest_document["protocol"]["argv"] = [
        sys.executable, "-m", "lifton",
    ]
    _write_json(candidate_manifest, manifest_document)

    Path(cell["artifacts"]["candidate_gff"]).write_text(
        "##gff-version 3\n"
        "chr1\tLiftOn\tgene\t1\t12\t.\t+\t.\tID=candidate-changed\n",
        encoding="utf-8",
    )
    errors, _ = controller.validate_artifacts(cell, 0)
    assert any("candidate GFF3 fingerprints" in error for error in errors)

    Path(cell["artifacts"]["reference_manifest"]).unlink()
    errors, _ = controller.validate_artifacts(cell, 0)
    assert any("reference_manifest is missing" in error for error in errors)


def test_e2e_release_manifest_protocol_must_match_pair_result(tmp_path):
    cell, raw = _paired_result_fixture(tmp_path, panel="e2e")
    manifest = Path(cell["artifacts"]["candidate_manifest"])
    document = json.loads(manifest.read_text())
    document["protocol"]["evaluation_flags"] = ["--different"]
    _write_json(manifest, document)

    errors = controller.validate_release_arm_manifest(
        cell,
        "candidate",
        manifest,
        version=raw["versions"]["candidate"],
        observed_fingerprints=raw["versions"]["candidate"]["fingerprints"],
    )

    assert any(
        "evaluation_flags disagrees with pair result" in error
        for error in errors
    )


def test_paired_artifact_gate_rejects_evaluator_tsv_tampering(tmp_path, monkeypatch):
    cell, raw = _paired_result_fixture(tmp_path)
    monkeypatch.setattr(
        controller,
        "_run_gff_validator",
        lambda *_args: ([], {"exit_code": 0, "result": {"is_valid": True}}),
    )
    transcript = (
        Path(cell["artifacts"]["result_json"]).parent
        / "evaluation"
        / "candidate.transcripts.tsv"
    )
    transcript.write_text("tampered\n", encoding="utf-8")

    errors, validation = controller.validate_artifacts(cell, 0)

    assert any("candidate evaluator TSV size disagrees" in error for error in errors)
    assert any("candidate evaluator TSV hash disagrees" in error for error in errors)
    assert validation["evaluation_artifacts"]["candidate"][
        "transcripts_tsv"
    ]["sha256"] != raw["versions"]["candidate"]["evaluation_artifacts"][
        "transcripts_tsv"
    ]["sha256"]


def test_paired_success_evidence_hashes_present_native_manifests(
        tmp_path, monkeypatch):
    cell, raw = _paired_result_fixture(tmp_path)
    native = (
        Path(cell["artifacts"]["result_json"]).parent
        / "candidate"
        / "lifton_output"
        / "run_manifest.json"
    )
    native.parent.mkdir(parents=True, exist_ok=True)
    native.write_text('{"native": true}\n', encoding="utf-8")
    native_record = {
        "path": str(native.resolve()),
        "present": True,
        "size": native.stat().st_size,
        "sha256": controller.sha256_file(native),
    }
    raw["versions"]["candidate"]["native_manifests"]["lift"] = native_record
    _write_json(Path(cell["artifacts"]["result_json"]), raw)
    manifest = Path(cell["artifacts"]["candidate_manifest"])
    document = json.loads(manifest.read_text())
    document["artifacts"]["native_manifests"]["lift"] = native_record
    _write_json(manifest, document)
    monkeypatch.setattr(
        controller,
        "_run_gff_validator",
        lambda *_args: ([], {"exit_code": 0, "result": {"is_valid": True}}),
    )

    errors, validation = controller.validate_artifacts(cell, 0)

    assert errors == []
    record = validation["artifacts"]["candidate_native_manifest_lift"]
    assert record["path"] == str(native.resolve())
    assert record["sha256"] == controller.sha256_file(native)


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
    assert "end-to-end evaluation profile does not have exit_code 0" in errors
    assert any("biological_summary" in error for error in errors)


def test_e2e_result_schema_accepts_only_consistent_nonempty_biology(tmp_path):
    row = {
        "dataset": "bee",
        "lift_profile": {
            "exit_code": 0,
            "wall_clock_seconds": 10.0,
            "peak_rss_mb": 20.0,
        },
        "eval_profile": {
            "exit_code": 0,
            "wall_clock_seconds": 2.0,
            "peak_rss_mb": 5.0,
        },
        **_e2e_biology(),
    }
    result = tmp_path / "result.json"
    cell = {"kind": "end_to_end", "benchmark": "bee"}
    _write_json(result, {"rows": [row]})

    assert controller.validate_result_schema(cell, result) == []

    mutations = [
        (("biological_summary", "lost_features"), 3, "does not equal"),
        (
            ("biological_summary", "mean_protein_identity"),
            float("nan"),
            "finite",
        ),
        (("evaluation_summary", "records"), 0, "positive integer"),
        (("score_summary", "malformed_records"), 1, "not zero"),
    ]
    for path, value, message in mutations:
        changed = json.loads(json.dumps(row))
        changed[path[0]][path[1]] = value
        _write_json(result, {"rows": [changed]})
        errors = controller.validate_result_schema(cell, result)
        assert any(message in error for error in errors), errors


def _short_watchdog(plan, *, hard=2.0, stall=1.0, poll=0.02, grace=0.05):
    plan["policy"].update({
        "subset_timeout_seconds": hard,
        "full_timeout_seconds": hard,
        "e2e_timeout_seconds": hard,
        "stall_timeout_seconds": stall,
        "watchdog_poll_seconds": poll,
        "terminate_grace_seconds": grace,
    })


def _pid_is_running(pid: int) -> bool:
    try:
        text = (Path("/proc") / str(pid) / "stat").read_text()
    except OSError:
        return False
    state = text[text.rfind(")") + 2:].split()[0]
    return state != "Z"


@pytest.mark.skipif(not Path("/proc").is_dir(), reason="watchdog identity is Linux /proc based")
def test_worker_hard_timeout_kills_descendants_and_preserves_evidence(tmp_path, monkeypatch):
    run_dir, plan, cell = _minimal_plan(tmp_path)
    cell_dir = Path(cell["cell_dir"])
    child_pid_path = cell_dir / "child.pid"
    partial_path = cell_dir / "partial.gff3"
    script = (
        "import pathlib,signal,subprocess,sys,time;"
        "signal.signal(signal.SIGTERM,signal.SIG_IGN);"
        "pathlib.Path(sys.argv[1]).write_text('partial');"
        "p=subprocess.Popen([sys.executable,'-c',"
        "'import signal,time;signal.signal(signal.SIGTERM,signal.SIG_IGN);time.sleep(60)']);"
        "pathlib.Path(sys.argv[2]).write_text(str(p.pid));"
        "time.sleep(60)"
    )
    cell["command"] = [
        sys.executable, "-c", script, str(partial_path), str(child_pid_path),
    ]
    _short_watchdog(plan, hard=0.25, stall=5.0)
    _initialize_run(run_dir, plan)
    monkeypatch.setattr(controller, "assert_matching_provenance", lambda _plan: None)

    try:
        assert controller.execute_cell(run_dir, cell["id"]) == 1
        child_pid = int(child_pid_path.read_text())
        for _ in range(20):
            if not _pid_is_running(child_pid):
                break
            time.sleep(0.05)
        assert not _pid_is_running(child_pid)
    finally:
        if child_pid_path.exists():
            child_pid = int(child_pid_path.read_text())
            if _pid_is_running(child_pid):
                os.kill(child_pid, 9)

    failed = json.loads((cell_dir / ".failed.json").read_text())
    exit_record = json.loads((cell_dir / "attempt-01.exit.json").read_text())
    assert failed["watchdog"]["reason"] == "hard_timeout"
    assert exit_record["watchdog"]["cleanup"]["term_sent"] is True
    assert exit_record["watchdog"]["cleanup"]["kill_sent"] is True
    assert partial_path.read_text() == "partial"
    assert (cell_dir / "attempt-01.stdout.log").exists()
    assert (cell_dir / "attempt-01.stderr.log").exists()
    assert (cell_dir / "attempt-01.watchdog.json").exists()


@pytest.mark.skipif(not Path("/proc").is_dir(), reason="stall watchdog is Linux /proc based")
def test_worker_stall_timeout_requires_no_progress_and_low_cpu(tmp_path, monkeypatch):
    run_dir, plan, cell = _minimal_plan(tmp_path)
    cell["command"] = [sys.executable, "-c", "import time;time.sleep(60)"]
    _short_watchdog(plan, hard=2.0, stall=0.15)
    _initialize_run(run_dir, plan)
    monkeypatch.setattr(controller, "assert_matching_provenance", lambda _plan: None)

    assert controller.execute_cell(run_dir, cell["id"]) == 1

    failed = json.loads((Path(cell["cell_dir"]) / ".failed.json").read_text())
    assert failed["watchdog"]["reason"] == "stall_timeout"


@pytest.mark.skipif(not Path("/proc").is_dir(), reason="stall watchdog is Linux /proc based")
@pytest.mark.parametrize("activity", ["file", "cpu"])
def test_worker_activity_prevents_false_stall_timeout(tmp_path, monkeypatch, activity):
    run_dir, plan, cell = _minimal_plan(tmp_path)
    progress = Path(cell["cell_dir"]) / "progress.log"
    if activity == "file":
        script = (
            "import pathlib,sys,time;"
            "p=pathlib.Path(sys.argv[1]);"
            "[(p.write_text(str(i)),time.sleep(.04)) for i in range(8)]"
        )
        cell["command"] = [sys.executable, "-c", script, str(progress)]
    else:
        script = "import time;t=time.monotonic()+.35\nwhile time.monotonic()<t: pass"
        cell["command"] = [sys.executable, "-c", script]
    _short_watchdog(plan, hard=2.0, stall=0.12)
    _initialize_run(run_dir, plan)
    monkeypatch.setattr(controller, "assert_matching_provenance", lambda _plan: None)

    assert controller.execute_cell(run_dir, cell["id"]) == 0

    watchdog = json.loads(
        (Path(cell["cell_dir"]) / "attempt-01.watchdog.json").read_text()
    )
    assert watchdog["reason"] is None


def test_worker_never_writes_success_when_validation_fails(tmp_path, monkeypatch):
    run_dir, plan, cell = _minimal_plan(tmp_path)
    _initialize_run(run_dir, plan)
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
    _initialize_run(run_dir, plan)
    monkeypatch.setattr(controller, "assert_matching_provenance", lambda _plan: None)
    validation = {
        "evaluation_artifacts": {
            "candidate": {
                "transcripts_tsv": {
                    "path": "/tmp/candidate.transcripts.tsv",
                    "size": 10,
                    "sha256": "a" * 64,
                },
            },
        },
    }
    monkeypatch.setattr(
        controller,
        "validate_artifacts",
        lambda _cell, _started: ([], validation),
    )

    assert controller.execute_cell(run_dir, cell["id"]) == 0
    success = json.loads((Path(cell["cell_dir"]) / ".success").read_text())
    assert success["fingerprint"] == cell["fingerprint"]
    assert success["validation"] == validation
    assert success["validation"]["evaluation_artifacts"]["candidate"][
        "transcripts_tsv"
    ]["sha256"] == "a" * 64
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
    _initialize_run(run_dir, plan)
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
    _initialize_run(run_dir, plan)
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


@pytest.mark.parametrize(
    ("cell_state", "controller_state", "expected"),
    (
        ("pending", "planned", controller.STATUS_EXIT_INCOMPLETE),
        ("running", "running", controller.STATUS_EXIT_INCOMPLETE),
        ("failed", "failed", controller.STATUS_EXIT_FAILED),
        ("pending", "blocked", controller.STATUS_EXIT_FAILED),
        ("success", "success", controller.STATUS_EXIT_SUCCESS),
    ),
)
def test_status_exit_code_is_zero_only_for_fully_successful_runs(
        tmp_path, cell_state, controller_state, expected):
    run_dir, plan, cell = _minimal_plan(
        tmp_path,
        run_name=f"status-{cell_state}-{controller_state}",
    )
    _initialize_run(run_dir, plan)
    controller._write_status(cell, cell_state, attempts=1)
    controller.atomic_write_json(run_dir / "controller.json", {
        "schema_version": controller.CONTROLLER_SCHEMA_VERSION,
        "run_id": plan["run_id"],
        "fingerprint": plan["fingerprint"],
        "state": controller_state,
    })

    result = controller.main([
        "status", "--run-dir", str(run_dir), "--json",
    ])

    assert result == expected


def test_reconcile_detects_same_size_same_mtime_artifact_hash_tampering(
        tmp_path, monkeypatch):
    run_dir, plan, cell = _minimal_plan(tmp_path)
    _initialize_run(run_dir, plan)
    evidence = Path(cell["cell_dir"]) / "evidence.txt"
    evidence.write_text("before\n", encoding="utf-8")
    stat = evidence.stat()
    validation = {
        "artifacts": {
            "evidence": controller._success_artifact_record(evidence),
        },
    }
    controller.atomic_write_json(Path(cell["cell_dir"]) / ".success", {
        "schema_version": controller.CONTROLLER_SCHEMA_VERSION,
        "cell_id": cell["id"],
        "fingerprint": cell["fingerprint"],
        "validation": validation,
    })
    controller._write_status(cell, "success", attempts=1)
    evidence.write_text("after!\n", encoding="utf-8")
    os.utime(evidence, ns=(stat.st_atime_ns, stat.st_mtime_ns))
    monkeypatch.setattr(
        controller, "assert_matching_provenance", lambda _plan: None,
    )

    summary = controller.reconcile_run(run_dir)

    assert summary["counts"]["failed"] == 1
    failed = json.loads(
        (Path(cell["cell_dir"]) / ".failed.json").read_text()
    )
    assert any("hash changed after success" in error for error in failed["errors"])


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
