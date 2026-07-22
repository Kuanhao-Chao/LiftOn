from __future__ import annotations

import json
from pathlib import Path

import pytest

from benchmarks import manifest_tools
from benchmarks.compare import (
    build_controller,
    campaign_orchestrator,
    campaign_profiles,
    protocol_analysis,
)


def _ready_acquisition_evidence(tmp_path: Path, monkeypatch):
    cache = tmp_path / "cache"
    cache.mkdir()
    manifest = tmp_path / "manifest.json"
    lock = cache / "acquisition.lock.json"
    registry = cache / "benchmarks.json"
    dataset_registry = cache / "datasets.json"
    preflight = cache / "preflight.json"
    manifest_document = {
        "profile_id": "canonical-v2",
        "sources": [{
            "id": "source-one",
            "files": [
                {"role": "reference_genome"},
                {"role": "reference_annotation"},
                {"role": "target_truth"},
            ],
        }],
        "scenarios": [
            {
                "id": "synthetic-one",
                "kind": "synthetic",
                "inputs": {
                    "reference_genome": "source-one:reference_genome",
                    "reference_annotation": "source-one:reference_annotation",
                },
                "design": {"seed": "unit-test"},
            },
            {
                "id": "biological-one",
                "kind": "biological",
                "inputs": {"target_truth": "source-one:target_truth"},
                "design": {},
            },
        ],
    }
    lock_document = {"schema_version": 1, "files": []}
    manifest.write_text(json.dumps(manifest_document))
    lock.write_text(json.dumps(lock_document))
    registry.write_text(json.dumps({"benchmarks": []}))
    dataset_registry.write_text(json.dumps({"datasets": []}))

    artifacts = cache / "runtime" / "artifacts"
    artifacts.mkdir(parents=True)

    def artifact(name: str) -> dict:
        path = artifacts / name
        path.write_text(name + "\n", encoding="utf-8")
        return campaign_orchestrator._file_record(path)

    source_provenance = {
        "source-one": {
            role: {
                "source": artifact(f"source-{role}"),
                "runtime": artifact(f"runtime-{role}"),
                "transform": "verified-copy",
            }
            for role in (
                "reference_genome", "reference_annotation", "target_truth",
            )
        },
    }
    synthetic_provenance = {
        "synthetic-one": {
            "design_sha256": manifest_tools.canonical_sha256(
                manifest_document["scenarios"][0]["design"]
            ),
            "source": {
                "fasta": source_provenance["source-one"][
                    "reference_genome"
                ]["runtime"],
                "gff": source_provenance["source-one"][
                    "reference_annotation"
                ]["runtime"],
            },
            "transform": {"kind": "unit-test"},
            "transform_manifest": artifact("synthetic-transform.json"),
            "outputs": {
                "target_fasta": artifact("synthetic-target.fa"),
                "truth_gff": artifact("synthetic-truth.gff3"),
                "ortholog_map": artifact("synthetic-map.json"),
            },
        },
    }
    full_mapping = artifact("full-ortholog-map.json")
    ortholog_mappings = {
        "biological-one": {
            "path": full_mapping["path"],
            "sha256": full_mapping["sha256"],
            "id_policy": "ortholog-map",
            "entries": 1,
        },
    }
    panel_truth_output = artifact("subset-truth.gff3")
    panel_map_output = artifact("subset-ortholog-map.json")
    panel_truth_provenance = {
        "biological-one": {
            "schema_version": 2,
            "method": "canonical-v2-panel-truth-v2",
            "scenario_id": "biological-one",
            "inputs": {
                "subset_manifest": artifact("subset-manifest.json"),
                "source_subset_annotation": artifact("source-subset.gff3"),
                "target_subset_fasta": artifact("target-subset.fa"),
                "full_truth_gff": source_provenance["source-one"][
                    "target_truth"
                ]["runtime"],
                "full_ortholog_map": full_mapping,
            },
            "filter": {
                "inputs": {
                    "truth_gff": source_provenance["source-one"][
                        "target_truth"
                    ]["runtime"],
                },
                "output": panel_truth_output,
            },
            "mapping": {
                "sha256": manifest_tools.canonical_sha256({"mappings": []}),
                "counts": {},
            },
            "outputs": {
                "truth_gff": panel_truth_output,
                "ortholog_map": panel_map_output,
            },
            "manifest": artifact("panel-truth-manifest.json"),
        },
    }
    preflight.write_text(json.dumps({
        "schema_version": 2,
        "manifest_sha256": manifest_tools.canonical_sha256(manifest_document),
        "acquisition_lock_sha256": manifest_tools.canonical_sha256(lock_document),
        "campaign_ready": True,
        "registries_exported": True,
        "blockers": [],
        "remaining_actions": [],
        "registries": {
            "benchmark": campaign_orchestrator._file_record(registry),
            "dataset": campaign_orchestrator._file_record(dataset_registry),
        },
        "source_provenance": source_provenance,
        "synthetic_provenance": synthetic_provenance,
        "ortholog_mappings": ortholog_mappings,
        "panel_truth_provenance": panel_truth_provenance,
    }))
    monkeypatch.setattr(
        manifest_tools, "validate_manifest", lambda document: document,
    )
    monkeypatch.setattr(
        manifest_tools,
        "verify_acquisition_lock",
        lambda *_args, **_kwargs: {"verified": True},
    )
    arguments = {
        "manifest_path": manifest,
        "acquisition_lock": lock,
        "cache_root": cache,
        "acquisition_preflight": preflight,
        "registry": registry,
        "dataset_registry": dataset_registry,
    }
    return arguments, preflight


def test_acquisition_evidence_requires_final_ready_preflight(
        tmp_path, monkeypatch):
    arguments, preflight = _ready_acquisition_evidence(tmp_path, monkeypatch)

    evidence = campaign_orchestrator._acquisition_evidence(
        {"id": "canonical-v2"}, **arguments,
    )
    assert evidence["readiness"]["campaign_ready"] is True
    assert evidence["preflight"] == campaign_orchestrator._file_record(preflight)
    derived = evidence["readiness"]["derived_artifacts"]
    assert derived["source_provenance"]["ids"] == ["source-one"]
    assert derived["synthetic_provenance"]["ids"] == ["synthetic-one"]
    assert derived["ortholog_mappings"]["ids"] == ["biological-one"]
    assert derived["panel_truth_provenance"]["ids"] == ["biological-one"]
    assert any(
        location.endswith("/outputs/truth_gff")
        for location in derived["panel_truth_provenance"]["files"]
    )
    assert any(
        location.endswith("/manifest")
        for location in derived["panel_truth_provenance"]["files"]
    )

    document = json.loads(preflight.read_text())
    document["campaign_ready"] = False
    preflight.write_text(json.dumps(document))
    with pytest.raises(RuntimeError, match="not campaign-ready"):
        campaign_orchestrator._acquisition_evidence(
            {"id": "canonical-v2"}, **arguments,
        )


@pytest.mark.parametrize(
    "field",
    (
        "source_provenance",
        "synthetic_provenance",
        "ortholog_mappings",
        "panel_truth_provenance",
    ),
)
def test_acquisition_evidence_binds_derived_ids_exactly(
        tmp_path, monkeypatch, field):
    arguments, preflight = _ready_acquisition_evidence(tmp_path, monkeypatch)
    document = json.loads(preflight.read_text())
    document[field]["unexpected-id"] = {}
    preflight.write_text(json.dumps(document))

    with pytest.raises(RuntimeError, match=rf"{field} IDs do not match"):
        campaign_orchestrator._acquisition_evidence(
            {"id": "canonical-v2"}, **arguments,
        )


def test_campaign_resume_rejects_mutated_derived_truth_artifact(
        tmp_path, monkeypatch):
    arguments, preflight = _ready_acquisition_evidence(tmp_path, monkeypatch)
    profile = {
        "id": "canonical-v2",
        "digest": "profile-digest",
        "registry": {"path": str(tmp_path / "profiles.json")},
    }
    acquisition = campaign_orchestrator._acquisition_evidence(
        profile, **arguments,
    )
    monkeypatch.setattr(
        campaign_profiles,
        "load_profile",
        lambda *_args, **_kwargs: profile,
    )
    plan = {
        "profile": {
            "id": profile["id"],
            "digest": profile["digest"],
            "registry": profile["registry"],
        },
        "inputs": {
            "registry": str(arguments["registry"]),
            "dataset_registry": str(arguments["dataset_registry"]),
            "acquisition": acquisition,
        },
        "children": [],
    }
    campaign_orchestrator.assert_campaign_provenance(plan)

    document = json.loads(preflight.read_text())
    truth_path = Path(
        document["panel_truth_provenance"]["biological-one"]["outputs"]
        ["truth_gff"]["path"]
    )
    truth_path.write_text("mutated derived truth\n", encoding="utf-8")

    with pytest.raises(RuntimeError, match="changed after preparation"):
        campaign_orchestrator.assert_campaign_provenance(plan)


def _gate_registry(path: Path, cases: int = 1) -> Path:
    document = {
        "schema_version": 1,
        "profiles": [{
            "id": "unit",
            "description": "Unit-test gate campaign.",
            "campaigns": [
                {
                    "id": f"gates-{index}",
                    "stage": "gates",
                    "ids": [
                        "python-suite",
                        "fast-gate",
                        "benchmark-gate",
                    ],
                    "repetitions": 1,
                    "threads": 1,
                    "candidate_mode": None,
                    "reference_mode": None,
                    "baseline_policy": "not_applicable",
                    "truth_policy": "not_applicable",
                }
                for index in range(1, cases + 1)
            ],
        }],
    }
    path.write_text(json.dumps(document))
    return path


def _create_gate_campaign(tmp_path: Path, *, cases: int = 1):
    profile_registry = _gate_registry(
        tmp_path / "profiles.json",
        cases=cases,
    )
    return campaign_orchestrator.create_campaign_plan(
        run_id="unit-campaign",
        profile_id="unit",
        campaigns_root=tmp_path / "campaigns",
        profile_registry=profile_registry,
        candidate_root=tmp_path / "candidate",
        candidate_sha="a" * 40,
        reference_root=tmp_path / "reference",
        reference_sha="b" * 40,
        lifton_executable=Path("/bin/false"),
    )


def test_campaign_plan_initializes_and_reloads_profiled_children(tmp_path):
    campaign_dir, plan, children = _create_gate_campaign(tmp_path)
    campaign_orchestrator.initialize_campaign(campaign_dir, plan, children)

    loaded = campaign_orchestrator.load_campaign_plan(campaign_dir)
    summary = campaign_orchestrator.summarize_campaign(loaded)

    assert loaded["fingerprint"] == plan["fingerprint"]
    assert loaded["full_profile"] is True
    assert loaded["profile"]["selected_campaign_ids"] == ["gates-1"]
    assert summary["state"] == "planned"
    assert summary["counts"] == {"pending": 1}
    assert campaign_orchestrator.status_exit_code(summary) == (
        build_controller.STATUS_EXIT_INCOMPLETE
    )
    child_plan = build_controller.load_plan(
        Path(loaded["children"][0]["root"])
    )
    assert child_plan["campaign_case"]["case"]["id"] == "gates-1"
    assert all(
        cell["campaign_case"] == child_plan["campaign_case"]
        for cell in child_plan["cells"]
    )


def test_campaign_stops_after_first_failed_child(tmp_path, monkeypatch):
    campaign_dir, plan, children = _create_gate_campaign(tmp_path, cases=2)
    campaign_orchestrator.initialize_campaign(campaign_dir, plan, children)
    calls = []
    monkeypatch.setattr(
        campaign_orchestrator,
        "assert_campaign_provenance",
        lambda _plan: None,
    )
    monkeypatch.setattr(
        build_controller,
        "scheduler_loop",
        lambda root: calls.append(Path(root)) or 1,
    )

    assert campaign_orchestrator.run_campaign(campaign_dir) == 1
    assert len(calls) == 1
    state = build_controller.read_json(campaign_dir / "campaign_state.json")
    assert state["state"] == "failed"
    assert state["current_campaign_id"] == "gates-1"
    assert "no automatic retry" in state["error"]


def test_protocol_child_embeds_the_shared_frozen_schedule(tmp_path, monkeypatch):
    profile = campaign_profiles.load_profile("canonical-v2")
    case = campaign_profiles.get_case(
        profile, "v2_protocol_thread_scaling_bee",
    )
    monkeypatch.setattr(
        campaign_orchestrator,
        "_protocol_provenance",
        lambda *_args, **_kwargs: {"source": {}, "inputs": {}},
    )

    child = campaign_orchestrator._protocol_child(
        case,
        campaign_profiles.case_identity(profile, case),
        campaign_dir=tmp_path / "campaign",
        candidate_root=tmp_path / "candidate",
        candidate_sha="a" * 40,
        lifton_executable=Path("/bin/false"),
        registry=build_controller.DEFAULT_REGISTRY,
        dataset_registry=build_controller.DEFAULT_DATASET_REGISTRY,
    )

    expected = protocol_analysis.protocol_cases("thread_scaling")
    observed = [
        {
            key: cell[key] for key in (
                "case_id", "kind", "dataset", "repetition",
                "position", "threads", "mode",
            )
        }
        for cell in child["cells"]
    ]
    assert observed == expected
    campaign_orchestrator.validate_campaign_plan({
        "schema_version": campaign_orchestrator.SCHEMA_VERSION,
        "run_id": "unit",
        "campaign_dir": str((tmp_path / "campaign").resolve()),
        "profile": {
            "selected_campaign_ids": ["v2_protocol_thread_scaling_bee"],
        },
        "children": [{
            "id": "v2_protocol_thread_scaling_bee",
            "kind": "protocol",
            "root": child["root"],
            "fingerprint": child["fingerprint"],
            "plan": child,
        }],
        "fingerprint": "",
    } | {
        "fingerprint": build_controller.canonical_hash({
            "schema_version": campaign_orchestrator.SCHEMA_VERSION,
            "run_id": "unit",
            "campaign_dir": str((tmp_path / "campaign").resolve()),
            "profile": {
                "selected_campaign_ids": ["v2_protocol_thread_scaling_bee"],
            },
            "children": [{
                "id": "v2_protocol_thread_scaling_bee",
                "kind": "protocol",
                "root": child["root"],
                "fingerprint": child["fingerprint"],
                "plan": child,
            }],
        }),
    })


def test_partial_profile_cannot_be_finalized(tmp_path):
    campaign_dir, plan, children = _create_gate_campaign(tmp_path, cases=2)
    plan["full_profile"] = False
    plan["fingerprint"] = campaign_orchestrator._hash_without_fingerprint(plan)
    campaign_orchestrator.initialize_campaign(campaign_dir, plan, children)

    with pytest.raises(ValueError, match="partial profile"):
        campaign_orchestrator.finalize_campaign(campaign_dir)


def test_canonical_v2_requires_verified_acquisition_before_materialization(
        tmp_path):
    with pytest.raises(RuntimeError, match="acquisition-lock"):
        campaign_orchestrator.create_campaign_plan(
            run_id="v2-preflight",
            profile_id="canonical-v2",
            campaigns_root=tmp_path / "campaigns",
            profile_registry=campaign_profiles.DEFAULT_PROFILE_REGISTRY,
            candidate_root=tmp_path / "candidate",
            candidate_sha="a" * 40,
            reference_root=tmp_path / "reference",
            reference_sha="b" * 40,
            lifton_executable=Path("/bin/false"),
        )


def test_materialization_preflight_names_missing_registry_ids(tmp_path):
    benchmark_registry = tmp_path / "benchmarks.json"
    dataset_registry = tmp_path / "datasets.json"
    benchmark_registry.write_text(json.dumps({
        "benchmarks": [{"id": "present"}],
    }))
    dataset_registry.write_text(json.dumps({
        "datasets": [{"id": "present_e2e"}],
    }))
    cases = [
        {"stage": "paired-subset", "ids": ["present", "missing_benchmark"]},
        {"stage": "paired-e2e", "ids": ["present_e2e", "missing_e2e"]},
    ]

    with pytest.raises(RuntimeError, match="missing_benchmark") as error:
        campaign_orchestrator._validate_materialized_profile_inputs(
            cases,
            registry=benchmark_registry,
            dataset_registry=dataset_registry,
        )
    assert "missing_e2e" in str(error.value)
    assert "never invents downloadable registry rows" in str(error.value)


def test_materialization_preflight_rejects_required_e2e_truth_without_map(
        tmp_path):
    benchmark_registry = tmp_path / "benchmarks.json"
    dataset_registry = tmp_path / "datasets.json"
    benchmark_registry.write_text(json.dumps({"benchmarks": []}))
    dataset_registry.write_text(json.dumps({
        "datasets": [{
            "id": "cross",
            "truth_gff": "/data/independent.gff3",
            "truth_id_policy": "ortholog-map",
        }],
    }))
    cases = [{
        "id": "e2e-truth",
        "stage": "paired-e2e",
        "ids": ["cross"],
        "truth_policy": "target_truth_required",
    }]

    with pytest.raises(
        RuntimeError,
        match="target_truth_required needs.*ortholog_map",
    ):
        campaign_orchestrator._validate_materialized_profile_inputs(
            cases,
            registry=benchmark_registry,
            dataset_registry=dataset_registry,
        )


def test_materialization_preflight_selects_panel_specific_truth(tmp_path):
    benchmark_registry = tmp_path / "benchmarks.json"
    benchmark_registry.write_text(json.dumps({
        "benchmarks": [{
            "id": "cross",
            "target_truth_by_panel": {
                "subset": {
                    "gff": "/truth/subset.gff3",
                    "ortholog_map": "/truth/subset-map.json",
                    "id_policy": "ortholog-map",
                },
                "full": {
                    "gff": "/truth/full.gff3",
                    "ortholog_map": "/truth/full-map.json",
                    "id_policy": "ortholog-map",
                },
            },
        }],
    }))
    dataset_registry = tmp_path / "datasets.json"
    dataset_registry.write_text(json.dumps({"datasets": []}))

    campaign_orchestrator._validate_materialized_profile_inputs(
        [{
            "id": "subset-truth",
            "stage": "paired-subset",
            "ids": ["cross"],
            "truth_policy": "target_truth_required",
        }],
        registry=benchmark_registry,
        dataset_registry=dataset_registry,
    )


def test_synthetic_truth_preflight_requires_generated_coordinate_map(
        tmp_path):
    truth = tmp_path / "truth.gff3"
    mapping = tmp_path / "ortholog-map.json"
    truth.write_text("##gff-version 3\n")
    mapping.write_text(json.dumps({
        "schema_version": 1,
        "method": "deterministic-synthetic-coordinate-map-v1",
        "mappings": [{
            "source_id": "g1",
            "truth_ids": [],
            "feature_type": "gene",
            "status": "breakpoint_crossing",
        }],
    }))
    benchmark_registry = tmp_path / "benchmarks.json"
    benchmark_registry.write_text(json.dumps({
        "benchmarks": [{
            "id": "synthetic",
            "target_truth": {
                "gff": str(truth),
                "ortholog_map": str(mapping),
                "id_policy": "ortholog-map",
            },
        }],
    }))
    dataset_registry = tmp_path / "datasets.json"
    dataset_registry.write_text(json.dumps({"datasets": []}))
    cases = [{
        "id": "synthetic-truth",
        "stage": "paired-subset",
        "ids": ["synthetic"],
        "truth_policy": "synthetic_exact_required",
    }]

    campaign_orchestrator._validate_materialized_profile_inputs(
        cases,
        registry=benchmark_registry,
        dataset_registry=dataset_registry,
    )

    document = json.loads(mapping.read_text())
    document["method"] = "unverified-map"
    mapping.write_text(json.dumps(document))
    with pytest.raises(RuntimeError, match="provenance is unsupported"):
        campaign_orchestrator._validate_materialized_profile_inputs(
            cases,
            registry=benchmark_registry,
            dataset_registry=dataset_registry,
        )
