#!/usr/bin/env python3
"""Resumable, fail-closed meta-controller for benchmark campaign profiles.

One immutable meta-plan sequences child ``build_controller`` runs.  Ordinary
stages retain their existing detached-cell execution and global lease;
single-candidate protocol schedules run in their frozen order.  The
orchestrator never retries automatically and never writes the canonical
four-way baseline.
"""
from __future__ import annotations

import argparse
import json
import os
import traceback
from collections import Counter
from dataclasses import asdict, replace
from pathlib import Path
from typing import Any, Mapping, Sequence

from . import build_controller, campaign_profiles, protocol_analysis


SCHEMA_VERSION = 1
DEFAULT_CAMPAIGNS_ROOT = build_controller.DEFAULT_RUNS_ROOT / "_campaigns"
ACTIVE_STATES = {"running", "launching"}
SUCCESS_STATE = "success"
FAILED_STATE = "failed"


def _hash_without_fingerprint(document: Mapping[str, Any]) -> str:
    return build_controller.canonical_hash({
        key: value for key, value in document.items()
        if key != "fingerprint"
    })


def _file_record(path: Path) -> dict[str, Any]:
    path = Path(path).resolve()
    stat = path.stat()
    return {
        "path": str(path),
        "size": stat.st_size,
        "sha256": build_controller.sha256_file(path),
    }


def _exact_provenance_block(
    preflight: Mapping[str, Any],
    field: str,
    expected_ids: set[str],
) -> Mapping[str, Any]:
    block = preflight.get(field)
    if not isinstance(block, Mapping):
        raise RuntimeError(
            f"canonical-v2 acquisition preflight {field} is missing"
        )
    observed_ids = set(block)
    if observed_ids != expected_ids or any(
        not isinstance(identifier, str) for identifier in block
    ):
        raise RuntimeError(
            f"canonical-v2 acquisition preflight {field} IDs do not match "
            f"the manifest; missing={sorted(expected_ids - observed_ids)}, "
            f"unexpected={sorted(observed_ids - expected_ids, key=str)}"
        )
    if any(not isinstance(value, Mapping) for value in block.values()):
        raise RuntimeError(
            f"canonical-v2 acquisition preflight {field} contains a malformed "
            "record"
        )
    return block


def _require_declared_file_record(value: Any, location: str) -> Mapping[str, Any]:
    if (
        not isinstance(value, Mapping)
        or not isinstance(value.get("path"), str)
        or not value["path"]
        or not isinstance(value.get("sha256"), str)
    ):
        raise RuntimeError(
            f"canonical-v2 acquisition preflight file record is malformed: "
            f"{location}"
        )
    return value


def _verify_declared_file_records(
    document: Mapping[str, Any],
    field: str,
) -> dict[str, dict[str, Any]]:
    """Recursively verify and normalize every path/SHA record in one block."""

    verified: dict[str, dict[str, Any]] = {}

    def pointer(parts: Sequence[str]) -> str:
        return "/" + "/".join(
            part.replace("~", "~0").replace("/", "~1") for part in parts
        )

    def visit(value: Any, parts: tuple[str, ...]) -> None:
        if isinstance(value, Mapping):
            if "path" in value:
                location = pointer(parts)
                record = _require_declared_file_record(value, location)
                digest = record["sha256"]
                if (
                    len(digest) != 64
                    or any(character not in "0123456789abcdef" for character in digest)
                ):
                    raise RuntimeError(
                        "canonical-v2 acquisition preflight SHA-256 is malformed: "
                        f"{location}"
                    )
                declared_path = Path(record["path"])
                if not declared_path.is_absolute():
                    raise RuntimeError(
                        "canonical-v2 acquisition preflight artifact path is not "
                        f"absolute: {location}"
                    )
                try:
                    stat = declared_path.stat()
                    if not declared_path.is_file() or stat.st_size <= 0:
                        raise RuntimeError(
                            "canonical-v2 acquisition preflight artifact is empty "
                            f"or not a file: {location}"
                        )
                    live = _file_record(declared_path)
                except OSError as exc:
                    raise RuntimeError(
                        "canonical-v2 acquisition preflight artifact cannot be "
                        f"verified: {location}"
                    ) from exc
                size = record.get("size")
                if (
                    live["path"] != record["path"]
                    or live["sha256"] != digest
                    or (
                        "size" in record
                        and (
                            not isinstance(size, int)
                            or isinstance(size, bool)
                            or size <= 0
                            or live["size"] != size
                        )
                    )
                ):
                    raise RuntimeError(
                        "canonical-v2 acquisition preflight artifact changed after "
                        f"preparation: {location}"
                    )
                verified[location] = live
            for key in sorted(value, key=str):
                visit(value[key], (*parts, str(key)))
        elif isinstance(value, list):
            for index, item in enumerate(value):
                visit(item, (*parts, str(index)))

    visit(document, (field,))
    return verified


def _canonical_v2_derived_evidence(
    manifest: Mapping[str, Any],
    preflight: Mapping[str, Any],
) -> dict[str, Any]:
    """Bind ready derived provenance to the manifest and verify live bytes."""

    sources = manifest.get("sources")
    scenarios = manifest.get("scenarios")
    if not isinstance(sources, list) or not isinstance(scenarios, list):
        raise RuntimeError(
            "canonical-v2 manifest source/scenario identities are malformed"
        )
    source_rows = {
        row.get("id"): row
        for row in sources
        if isinstance(row, Mapping) and isinstance(row.get("id"), str)
    }
    scenario_rows = {
        row.get("id"): row
        for row in scenarios
        if isinstance(row, Mapping) and isinstance(row.get("id"), str)
    }
    if len(source_rows) != len(sources) or len(scenario_rows) != len(scenarios):
        raise RuntimeError(
            "canonical-v2 manifest source/scenario identities are duplicated or "
            "malformed"
        )
    synthetic_rows = {
        identifier: row for identifier, row in scenario_rows.items()
        if row.get("kind") == "synthetic"
    }
    biological_rows = {
        identifier: row for identifier, row in scenario_rows.items()
        if row.get("kind") == "biological"
    }
    source_provenance = _exact_provenance_block(
        preflight, "source_provenance", set(source_rows),
    )
    synthetic_provenance = _exact_provenance_block(
        preflight, "synthetic_provenance", set(synthetic_rows),
    )
    ortholog_mappings = _exact_provenance_block(
        preflight, "ortholog_mappings", set(biological_rows),
    )
    panel_truth_provenance = _exact_provenance_block(
        preflight, "panel_truth_provenance", set(biological_rows),
    )

    for source_id, source in source_rows.items():
        files = source.get("files")
        if not isinstance(files, list):
            raise RuntimeError(f"canonical-v2 manifest source files are malformed: {source_id}")
        roles = {
            item.get("role") for item in files
            if isinstance(item, Mapping) and isinstance(item.get("role"), str)
        }
        record = source_provenance[source_id]
        if len(roles) != len(files) or set(record) != roles:
            raise RuntimeError(
                "canonical-v2 acquisition preflight source roles do not match "
                f"the manifest: {source_id}"
            )
        for role in roles:
            role_record = record[role]
            if not isinstance(role_record, Mapping):
                raise RuntimeError(
                    "canonical-v2 acquisition preflight source role is malformed: "
                    f"{source_id}/{role}"
                )
            for artifact in ("source", "runtime"):
                _require_declared_file_record(
                    role_record.get(artifact),
                    f"source_provenance/{source_id}/{role}/{artifact}",
                )

    for scenario_id, scenario in synthetic_rows.items():
        record = synthetic_provenance[scenario_id]
        if record.get("design_sha256") != build_controller.canonical_hash(
            scenario.get("design")
        ):
            raise RuntimeError(
                "canonical-v2 synthetic provenance design does not match the "
                f"manifest: {scenario_id}"
            )
        source_records = record.get("source")
        outputs = record.get("outputs")
        if (
            not isinstance(source_records, Mapping)
            or set(source_records) != {"fasta", "gff"}
            or not isinstance(outputs, Mapping)
            or set(outputs) != {"target_fasta", "truth_gff", "ortholog_map"}
        ):
            raise RuntimeError(
                "canonical-v2 synthetic provenance inputs/outputs are malformed: "
                f"{scenario_id}"
            )
        bindings = scenario.get("inputs")
        if not isinstance(bindings, Mapping):
            raise RuntimeError(
                f"canonical-v2 synthetic manifest inputs are malformed: {scenario_id}"
            )
        for label, manifest_role in (
            ("fasta", "reference_genome"),
            ("gff", "reference_annotation"),
        ):
            try:
                source_id, source_role = str(bindings[manifest_role]).split(":", 1)
                expected = source_provenance[source_id][source_role]["runtime"]
            except (KeyError, TypeError, ValueError) as exc:
                raise RuntimeError(
                    "canonical-v2 synthetic provenance cannot resolve manifest "
                    f"binding: {scenario_id}/{manifest_role}"
                ) from exc
            if source_records[label] != expected:
                raise RuntimeError(
                    "canonical-v2 synthetic provenance source does not match the "
                    f"manifest binding: {scenario_id}/{manifest_role}"
                )
        _require_declared_file_record(
            record.get("transform_manifest"),
            f"synthetic_provenance/{scenario_id}/transform_manifest",
        )
        for role in ("target_fasta", "truth_gff", "ortholog_map"):
            _require_declared_file_record(
                outputs.get(role),
                f"synthetic_provenance/{scenario_id}/outputs/{role}",
            )

    for scenario_id, scenario in biological_rows.items():
        mapping = ortholog_mappings[scenario_id]
        _require_declared_file_record(
            mapping, f"ortholog_mappings/{scenario_id}",
        )
        if mapping.get("id_policy") != "ortholog-map":
            raise RuntimeError(
                f"canonical-v2 ortholog mapping policy is unsafe: {scenario_id}"
            )
        truth = panel_truth_provenance[scenario_id]
        outputs = truth.get("outputs")
        if (
            truth.get("schema_version") != 2
            or truth.get("method") != "canonical-v2-panel-truth-v2"
            or truth.get("scenario_id") != scenario_id
            or not isinstance(outputs, Mapping)
            or set(outputs) != {"truth_gff", "ortholog_map"}
        ):
            raise RuntimeError(
                "canonical-v2 panel truth provenance is not bound to the "
                f"manifest scenario: {scenario_id}"
            )
        for role in ("truth_gff", "ortholog_map"):
            _require_declared_file_record(
                outputs.get(role),
                f"panel_truth_provenance/{scenario_id}/outputs/{role}",
            )
        _require_declared_file_record(
            truth.get("manifest"),
            f"panel_truth_provenance/{scenario_id}/manifest",
        )
        inputs = truth.get("inputs")
        bindings = scenario.get("inputs")
        if not isinstance(inputs, Mapping) or not isinstance(bindings, Mapping):
            raise RuntimeError(
                f"canonical-v2 panel truth inputs are malformed: {scenario_id}"
            )
        try:
            source_id, source_role = str(bindings["target_truth"]).split(":", 1)
            expected_truth = source_provenance[source_id][source_role]["runtime"]
        except (KeyError, TypeError, ValueError) as exc:
            raise RuntimeError(
                "canonical-v2 panel truth cannot resolve target-truth manifest "
                f"binding: {scenario_id}"
            ) from exc
        if inputs.get("full_truth_gff") != expected_truth:
            raise RuntimeError(
                "canonical-v2 panel truth full GFF does not match the manifest "
                f"binding: {scenario_id}"
            )
        full_mapping = inputs.get("full_ortholog_map")
        if (
            not isinstance(full_mapping, Mapping)
            or full_mapping.get("path") != mapping.get("path")
            or full_mapping.get("sha256") != mapping.get("sha256")
        ):
            raise RuntimeError(
                "canonical-v2 panel truth ortholog map does not match the full "
                f"mapping: {scenario_id}"
            )

    blocks = {
        "source_provenance": source_provenance,
        "synthetic_provenance": synthetic_provenance,
        "ortholog_mappings": ortholog_mappings,
        "panel_truth_provenance": panel_truth_provenance,
    }
    return {
        name: {
            "ids": sorted(block),
            "document_sha256": build_controller.canonical_hash(block),
            "files": _verify_declared_file_records(block, name),
        }
        for name, block in blocks.items()
    }


def _acquisition_evidence(
    profile: Mapping[str, Any],
    *,
    manifest_path: Path | None,
    acquisition_lock: Path | None,
    cache_root: Path | None,
    acquisition_preflight: Path | None,
    registry: Path,
    dataset_registry: Path,
) -> dict[str, Any] | None:
    """Validate and freeze the canonical-v2 content-addressed data lock."""

    if profile["id"] != "canonical-v2":
        if any(
            value is not None
            for value in (
                manifest_path, acquisition_lock, cache_root,
                acquisition_preflight,
            )
        ):
            raise ValueError(
                "dataset acquisition options are only valid for canonical-v2"
            )
        return None
    from benchmarks import manifest_tools

    manifest_path = Path(
        manifest_path or manifest_tools.DEFAULT_MANIFEST
    ).resolve()
    if (
        acquisition_lock is None
        or cache_root is None
        or acquisition_preflight is None
    ):
        raise RuntimeError(
            "canonical-v2 planning requires --acquisition-lock and "
            "--dataset-cache-root and --acquisition-preflight. First review "
            "the immutable request plan "
            "with `python -m benchmarks.manifest_tools acquire --dry-run "
            "--cache-root <cache>`, acquire/materialize the declared bytes, "
            "and create a verified content-addressed lock."
        )
    acquisition_lock = Path(acquisition_lock).resolve()
    cache_root = Path(cache_root).resolve()
    acquisition_preflight = Path(acquisition_preflight).resolve()
    manifest = manifest_tools.validate_manifest(
        manifest_tools.load_json(manifest_path)
    )
    if manifest.get("profile_id") != profile["id"]:
        raise RuntimeError(
            f"dataset manifest profile_id {manifest.get('profile_id')!r} "
            f"does not match {profile['id']!r}"
        )
    verification = manifest_tools.verify_acquisition_lock(
        manifest,
        manifest_tools.load_json(acquisition_lock),
        cache_root,
    )
    preflight = manifest_tools.load_json(acquisition_preflight)
    if preflight.get("schema_version") != 2:
        raise RuntimeError("canonical-v2 acquisition preflight schema must be 2")
    if preflight.get("campaign_ready") is not True:
        raise RuntimeError("canonical-v2 acquisition preflight is not campaign-ready")
    if preflight.get("registries_exported") is not True:
        raise RuntimeError("canonical-v2 acquisition preflight did not export registries")
    if preflight.get("blockers") != [] or preflight.get("remaining_actions") != []:
        raise RuntimeError(
            "canonical-v2 acquisition preflight retains blockers or actions"
        )
    manifest_digest = manifest_tools.canonical_sha256(manifest)
    lock_document = manifest_tools.load_json(acquisition_lock)
    lock_digest = manifest_tools.canonical_sha256(lock_document)
    if preflight.get("manifest_sha256") != manifest_digest:
        raise RuntimeError("acquisition preflight manifest digest is stale")
    if preflight.get("acquisition_lock_sha256") != lock_digest:
        raise RuntimeError("acquisition preflight lock digest is stale")
    try:
        acquisition_preflight.relative_to(cache_root)
    except ValueError as exc:
        raise RuntimeError(
            "acquisition preflight must remain below the dataset cache root"
        ) from exc
    frozen_registries = preflight.get("registries")
    if not isinstance(frozen_registries, Mapping):
        raise RuntimeError("acquisition preflight registries are missing")
    expected_registries = {
        "benchmark": _file_record(Path(registry)),
        "dataset": _file_record(Path(dataset_registry)),
    }
    for name, expected in expected_registries.items():
        if frozen_registries.get(name) != expected:
            raise RuntimeError(
                f"acquisition preflight {name} registry is stale or mismatched"
            )
    verified_registry_records = {}
    for name, record in sorted(frozen_registries.items()):
        if not isinstance(name, str) or not isinstance(record, Mapping):
            raise RuntimeError("acquisition preflight registry record is malformed")
        try:
            live = _file_record(Path(str(record["path"])))
        except (KeyError, OSError) as exc:
            raise RuntimeError(
                f"acquisition preflight registry cannot be verified: {name}"
            ) from exc
        if dict(record) != live:
            raise RuntimeError(
                f"acquisition preflight registry changed after preparation: {name}"
            )
        verified_registry_records[name] = live
    derived_artifacts = _canonical_v2_derived_evidence(manifest, preflight)
    return {
        "manifest": _file_record(manifest_path),
        "lock": _file_record(acquisition_lock),
        "preflight": _file_record(acquisition_preflight),
        "cache_root": str(cache_root),
        "verification": verification,
        "readiness": {
            "schema_version": 2,
            "manifest_sha256": manifest_digest,
            "acquisition_lock_sha256": lock_digest,
            "campaign_ready": True,
            "registries_exported": True,
            "registries": verified_registry_records,
            "derived_artifacts": derived_artifacts,
        },
    }


def _validate_materialized_profile_inputs(
    cases: Sequence[Mapping[str, Any]],
    *,
    registry: Path,
    dataset_registry: Path,
) -> None:
    benchmark_ids = set(build_controller._benchmark_ids(registry))
    dataset_ids = set(build_controller._dataset_ids(dataset_registry))
    benchmark_document = build_controller.read_json(Path(registry))
    dataset_document = build_controller.read_json(Path(dataset_registry))
    benchmark_entries = {
        row.get("id"): row
        for row in benchmark_document.get("benchmarks", [])
        if isinstance(row, Mapping) and isinstance(row.get("id"), str)
    }
    dataset_entries = {
        row.get("id"): row
        for row in dataset_document.get("datasets", [])
        if isinstance(row, Mapping) and isinstance(row.get("id"), str)
    }
    missing_benchmarks = set()
    missing_datasets = set()
    truth_errors = []
    for case in cases:
        base = build_controller._base_stage(case["stage"])
        if base.startswith(("subset", "full")):
            missing_benchmarks.update(set(case["ids"]) - benchmark_ids)
            entries = benchmark_entries
            truth_panel = "subset" if base.startswith("subset") else "full"
        elif base.startswith("e2e") or case["stage"] in {
            "protocol-thread-scaling", "protocol-io-modes",
        }:
            missing_datasets.update(set(case["ids"]) - dataset_ids)
            entries = dataset_entries
            truth_panel = "full"
        else:
            entries = {}
            truth_panel = "full"
        truth_policy = case.get("truth_policy")
        if truth_policy not in {
            "target_truth_required", "synthetic_exact_required",
        }:
            continue
        for item_id in case["ids"]:
            entry = entries.get(item_id)
            if entry is None:
                continue
            target_truth_by_panel = entry.get("target_truth_by_panel")
            target_truth = (
                target_truth_by_panel.get(truth_panel)
                if isinstance(target_truth_by_panel, Mapping)
                else None
            )
            if not isinstance(target_truth, Mapping):
                target_truth = entry.get("target_truth")
            truth = target_truth if isinstance(target_truth, Mapping) else {}
            truth_gff = (
                truth.get("gff")
                or truth.get("truth_gff")
                or entry.get("truth_gff")
            )
            ortholog_map = (
                truth.get("ortholog_map")
                or entry.get("ortholog_map")
            )
            id_policy = str(
                truth.get(
                    "id_policy",
                    entry.get("truth_id_policy", "ortholog-map"),
                )
            ).strip().lower()
            prefix = f"{case['id']}/{item_id}"
            if not truth_gff:
                truth_errors.append(
                    f"{prefix}: missing independent truth_gff"
                )
                continue
            if id_policy != "ortholog-map" or not ortholog_map:
                truth_errors.append(
                    f"{prefix}: {truth_policy} needs "
                    "truth_id_policy=ortholog-map and a non-empty "
                    "ortholog_map"
                )
                continue
            if truth_policy == "synthetic_exact_required":
                map_path = Path(str(ortholog_map)).expanduser()
                if not map_path.is_absolute():
                    map_path = Path(registry).resolve().parent / map_path
                try:
                    map_document = build_controller.read_json(map_path)
                except (OSError, TypeError, ValueError) as exc:
                    truth_errors.append(
                        f"{prefix}: generated synthetic ortholog map is "
                        f"unavailable: {exc}"
                    )
                    continue
                if (
                    map_document.get("schema_version") != 1
                    or map_document.get("method")
                    != "deterministic-synthetic-coordinate-map-v1"
                    or not isinstance(map_document.get("mappings"), list)
                ):
                    truth_errors.append(
                        f"{prefix}: synthetic ortholog map provenance is "
                        "unsupported"
                    )
    if missing_benchmarks or missing_datasets:
        raise RuntimeError(
            "campaign profile inputs are not materialized in the runtime "
            "registries; "
            f"missing benchmark IDs={sorted(missing_benchmarks)}, "
            f"missing E2E dataset IDs={sorted(missing_datasets)}. Verify and "
            "acquire the immutable dataset lock, generate synthetic inputs, "
            "then add explicit path records to benchmarks.json/datasets.json. "
            "The orchestrator never invents downloadable registry rows."
        )
    if truth_errors:
        raise RuntimeError(
            "campaign truth inputs are not safely materialized: "
            + "; ".join(truth_errors)
        )


def _campaign_dir(
    campaigns_root: Path,
    run_id: str,
) -> Path:
    return Path(campaigns_root).resolve() / build_controller.safe_name(
        run_id, limit=100,
    )


def _protocol_kind(stage: str) -> str:
    try:
        return campaign_profiles.PROTOCOL_STAGES[stage]
    except KeyError as exc:
        raise ValueError(f"not a protocol stage: {stage}") from exc


def _paired_configuration(
    case: Mapping[str, Any],
    *,
    candidate_root: Path,
    candidate_sha: str,
    reference_root: Path,
    reference_sha: str,
    lifton_executable: Path,
    registry: Path,
    dataset_registry: Path,
) -> dict[str, Any]:
    return build_controller.paired_configuration(
        stage=case["stage"],
        candidate_root=candidate_root,
        candidate_sha=candidate_sha,
        reference_root=reference_root,
        reference_sha=reference_sha,
        repetitions=case["repetitions"],
        lifton_executable=lifton_executable,
        candidate_e2e_mode=case["candidate_mode"],
        reference_e2e_mode=case["reference_mode"],
        benchmark_registry=registry,
        dataset_registry=dataset_registry,
    )


def _protocol_provenance(
    case: Mapping[str, Any],
    *,
    candidate_root: Path,
    candidate_sha: str,
    lifton_executable: Path,
    registry: Path,
    dataset_registry: Path,
) -> dict[str, Any]:
    from . import release_evaluation

    source = release_evaluation.SourceSpec(
        label="candidate",
        root=Path(candidate_root).resolve(),
        sha=candidate_sha,
        lifton_executable=Path(lifton_executable).resolve(),
    )
    verified = release_evaluation.verify_source(source)
    dataset = case["ids"][0]
    inputs = release_evaluation.resolve_panel_inputs(
        "e2e",
        dataset,
        benchmark_registry=registry,
        dataset_registry=dataset_registry,
    )
    return {
        "source": verified,
        "inputs": release_evaluation.input_fingerprints(inputs),
    }


def _protocol_child(
    case: Mapping[str, Any],
    identity: Mapping[str, Any],
    *,
    campaign_dir: Path,
    candidate_root: Path,
    candidate_sha: str,
    lifton_executable: Path,
    registry: Path,
    dataset_registry: Path,
) -> dict[str, Any]:
    kind = _protocol_kind(case["stage"])
    root = campaign_dir / "protocol" / case["id"]
    schedule = protocol_analysis.protocol_cases(kind)
    cells = []
    for item in schedule:
        cell_dir = root / "cells" / item["case_id"]
        cells.append({
            **item,
            "cell_dir": str(cell_dir),
            "result_json": str(cell_dir / "result.json"),
        })
    child = {
        "schema_version": SCHEMA_VERSION,
        "kind": "protocol",
        "campaign_case": identity,
        "root": str(root),
        "protocol_kind": kind,
        "candidate": {
            "root": str(Path(candidate_root).resolve()),
            "sha": candidate_sha,
            "lifton_executable": str(Path(lifton_executable).resolve()),
        },
        "dataset_registry": str(Path(dataset_registry).resolve()),
        "benchmark_registry": str(Path(registry).resolve()),
        "provenance": _protocol_provenance(
            case,
            candidate_root=candidate_root,
            candidate_sha=candidate_sha,
            lifton_executable=lifton_executable,
            registry=registry,
            dataset_registry=dataset_registry,
        ),
        "cells": cells,
    }
    child["fingerprint"] = _hash_without_fingerprint(child)
    return child


def create_campaign_plan(
    *,
    run_id: str,
    profile_id: str,
    campaigns_root: Path,
    profile_registry: Path,
    candidate_root: Path,
    candidate_sha: str,
    reference_root: Path,
    reference_sha: str,
    lifton_executable: Path,
    registry: Path = build_controller.DEFAULT_REGISTRY,
    dataset_registry: Path = build_controller.DEFAULT_DATASET_REGISTRY,
    baseline: Path = build_controller.DEFAULT_BASELINE,
    base_policy: build_controller.Policy = build_controller.Policy(),
    campaign_ids: Sequence[str] | None = None,
    dataset_manifest: Path | None = None,
    acquisition_lock: Path | None = None,
    dataset_cache_root: Path | None = None,
    acquisition_preflight: Path | None = None,
) -> tuple[Path, dict[str, Any], list[tuple[Path, dict[str, Any]]]]:
    """Create every child plan in memory before writing any run state."""

    profile = campaign_profiles.load_profile(
        profile_id,
        registry=profile_registry,
    )
    selected = set(campaign_ids or [])
    available = set(campaign_profiles.profile_case_ids(profile))
    unknown = selected - available
    if unknown:
        raise ValueError(f"unknown profile campaign ids: {sorted(unknown)}")
    cases = [
        case for case in profile["campaigns"]
        if not selected or case["id"] in selected
    ]
    full_profile = not selected
    acquisition = _acquisition_evidence(
        profile,
        manifest_path=dataset_manifest,
        acquisition_lock=acquisition_lock,
        cache_root=dataset_cache_root,
        acquisition_preflight=acquisition_preflight,
        registry=Path(registry),
        dataset_registry=Path(dataset_registry),
    )
    _validate_materialized_profile_inputs(
        cases,
        registry=Path(registry),
        dataset_registry=Path(dataset_registry),
    )
    campaign_dir = _campaign_dir(campaigns_root, run_id)
    if campaign_dir.exists() and any(campaign_dir.iterdir()):
        raise FileExistsError(
            f"campaign directory already exists and is not empty: {campaign_dir}"
        )
    children = []
    controller_plans: list[tuple[Path, dict[str, Any]]] = []
    for index, case in enumerate(cases, start=1):
        identity = campaign_profiles.case_identity(profile, case)
        if case["stage"] in campaign_profiles.PROTOCOL_STAGES:
            child = _protocol_child(
                case,
                identity,
                campaign_dir=campaign_dir,
                candidate_root=candidate_root,
                candidate_sha=candidate_sha,
                lifton_executable=lifton_executable,
                registry=registry,
                dataset_registry=dataset_registry,
            )
            children.append({
                "id": case["id"],
                "kind": "protocol",
                "root": child["root"],
                "fingerprint": child["fingerprint"],
                "plan": child,
            })
            continue
        paired = None
        if case["stage"] in campaign_profiles.PAIRED_STAGES:
            paired = _paired_configuration(
                case,
                candidate_root=candidate_root,
                candidate_sha=candidate_sha,
                reference_root=reference_root,
                reference_sha=reference_sha,
                lifton_executable=lifton_executable,
                registry=registry,
                dataset_registry=dataset_registry,
            )
        policy = replace(
            base_policy,
            threads_per_cell=case["threads"],
            max_worker_threads=max(
                base_policy.max_worker_threads,
                case["threads"],
            ),
        )
        child_run_id = build_controller.safe_name(
            f"{run_id}--{index:02d}-{case['id']}",
            limit=100,
        )
        child_dir, child_plan = build_controller.create_plan(
            run_id=child_run_id,
            stage=case["stage"],
            requested_ids=case["ids"],
            runs_root=campaign_dir / "runs",
            repo_root=build_controller.REPO_ROOT,
            registry=registry,
            dataset_registry=dataset_registry,
            baseline=baseline,
            policy=policy,
            paired=paired,
            campaign_case=identity,
            profile_registry=profile_registry,
        )
        controller_plans.append((child_dir, child_plan))
        children.append({
            "id": case["id"],
            "kind": "controller",
            "root": str(child_dir),
            "fingerprint": child_plan["fingerprint"],
        })
    plan = {
        "schema_version": SCHEMA_VERSION,
        "run_id": build_controller.safe_name(run_id, limit=100),
        "created_at": build_controller.utc_now(),
        "campaign_dir": str(campaign_dir),
        "full_profile": full_profile,
        "profile": {
            "id": profile["id"],
            "digest": profile["digest"],
            "registry": profile["registry"],
            "selected_campaign_ids": [case["id"] for case in cases],
            "counts": campaign_profiles.campaign_counts({
                **profile,
                "campaigns": cases,
            }),
        },
        "sources": {
            "candidate": {
                "root": str(Path(candidate_root).resolve()),
                "sha": candidate_sha,
            },
            "reference": {
                "root": str(Path(reference_root).resolve()),
                "sha": reference_sha,
            },
            "lifton_executable": str(Path(lifton_executable).resolve()),
        },
        "inputs": {
            "registry": str(Path(registry).resolve()),
            "dataset_registry": str(Path(dataset_registry).resolve()),
            "baseline": str(Path(baseline).resolve()),
            "acquisition": acquisition,
        },
        "policy": asdict(base_policy),
        "children": children,
    }
    plan["fingerprint"] = _hash_without_fingerprint(plan)
    validate_campaign_plan(plan)
    return campaign_dir, plan, controller_plans


def validate_campaign_plan(plan: Mapping[str, Any]) -> None:
    if plan.get("schema_version") != SCHEMA_VERSION:
        raise ValueError("unsupported campaign plan schema")
    campaign_dir = Path(str(plan.get("campaign_dir", ""))).resolve()
    children = plan.get("children")
    if not isinstance(children, list) or not children:
        raise ValueError("campaign plan children are missing")
    ids = []
    for child in children:
        if not isinstance(child, Mapping):
            raise ValueError("campaign plan contains a malformed child")
        ids.append(child.get("id"))
        root = Path(str(child.get("root", ""))).resolve()
        if not root.is_relative_to(campaign_dir):
            raise ValueError(f"campaign child escapes campaign directory: {root}")
        if child.get("kind") == "protocol":
            protocol_plan = child.get("plan")
            if not isinstance(protocol_plan, Mapping):
                raise ValueError("protocol child plan is missing")
            if protocol_plan.get("fingerprint") != _hash_without_fingerprint(
                protocol_plan
            ):
                raise ValueError("protocol child fingerprint is invalid")
            if child.get("fingerprint") != protocol_plan["fingerprint"]:
                raise ValueError("protocol child fingerprint disagrees with meta-plan")
            expected = protocol_analysis.protocol_cases(
                protocol_plan["protocol_kind"]
            )
            observed = [
                {
                    key: cell[key] for key in (
                        "case_id", "kind", "dataset", "repetition",
                        "position", "threads", "mode",
                    )
                }
                for cell in protocol_plan.get("cells", [])
            ]
            if observed != expected:
                raise ValueError("protocol child schedule is not canonical")
        elif child.get("kind") != "controller":
            raise ValueError(f"unknown campaign child kind: {child.get('kind')}")
    if len(ids) != len(set(ids)):
        raise ValueError("campaign plan child ids are duplicated")
    selected = plan.get("profile", {}).get("selected_campaign_ids")
    if ids != selected:
        raise ValueError("campaign child order disagrees with the profile")
    if plan.get("fingerprint") != _hash_without_fingerprint(plan):
        raise ValueError("campaign plan fingerprint is invalid")


def initialize_campaign(
    campaign_dir: Path,
    plan: Mapping[str, Any],
    controller_plans: Sequence[tuple[Path, Mapping[str, Any]]],
) -> None:
    validate_campaign_plan(plan)
    campaign_dir = Path(campaign_dir)
    if campaign_dir.exists() and any(campaign_dir.iterdir()):
        raise FileExistsError(
            f"campaign directory already exists and is not empty: {campaign_dir}"
        )
    campaign_dir.mkdir(parents=True, exist_ok=True)
    build_controller.atomic_write_json(campaign_dir / "campaign_plan.json", plan)
    by_root = {str(root): child for root, child in controller_plans}
    for child in plan["children"]:
        root = Path(child["root"])
        if child["kind"] == "controller":
            build_controller.initialize_run(root, by_root[str(root)])
            continue
        protocol_plan = child["plan"]
        root.mkdir(parents=True, exist_ok=True)
        build_controller.atomic_write_json(root / "plan.json", protocol_plan)
        for cell in protocol_plan["cells"]:
            cell_dir = Path(cell["cell_dir"])
            cell_dir.mkdir(parents=True, exist_ok=True)
            build_controller.atomic_write_json(cell_dir / "status.json", {
                "schema_version": SCHEMA_VERSION,
                "case_id": cell["case_id"],
                "state": "pending",
                "attempts": 0,
                "updated_at": build_controller.utc_now(),
            })
    build_controller.atomic_write_json(campaign_dir / "campaign_state.json", {
        "schema_version": SCHEMA_VERSION,
        "run_id": plan["run_id"],
        "fingerprint": plan["fingerprint"],
        "state": "planned",
        "current_campaign_id": None,
        "updated_at": build_controller.utc_now(),
    })


def load_campaign_plan(campaign_dir: Path) -> dict[str, Any]:
    plan = build_controller.read_json(
        Path(campaign_dir).resolve() / "campaign_plan.json"
    )
    validate_campaign_plan(plan)
    for child in plan["children"]:
        if child["kind"] == "protocol":
            disk_plan = build_controller.read_json(
                Path(child["root"]) / "plan.json"
            )
            if disk_plan != child["plan"]:
                raise ValueError(
                    f"protocol child plan changed: {child['id']}"
                )
            continue
        child_plan = build_controller.load_plan(Path(child["root"]))
        if child_plan["fingerprint"] != child["fingerprint"]:
            raise ValueError(
                f"child plan fingerprint changed: {child['id']}"
            )
    return plan


def _assert_protocol_provenance(protocol_plan: Mapping[str, Any]) -> None:
    case = protocol_plan["campaign_case"]["case"]
    current = _protocol_provenance(
        case,
        candidate_root=Path(protocol_plan["candidate"]["root"]),
        candidate_sha=protocol_plan["candidate"]["sha"],
        lifton_executable=Path(protocol_plan["candidate"]["lifton_executable"]),
        registry=Path(protocol_plan["benchmark_registry"]),
        dataset_registry=Path(protocol_plan["dataset_registry"]),
    )
    if current != protocol_plan["provenance"]:
        raise RuntimeError(
            f"protocol provenance changed: {case['id']}; start a new campaign"
        )


def assert_campaign_provenance(plan: Mapping[str, Any]) -> None:
    profile = campaign_profiles.load_profile(
        plan["profile"]["id"],
        registry=Path(plan["profile"]["registry"]["path"]),
    )
    if (
        profile["digest"] != plan["profile"]["digest"]
        or profile["registry"] != plan["profile"]["registry"]
    ):
        raise RuntimeError(
            "campaign profile registry changed; start a new campaign"
        )
    acquisition = plan.get("inputs", {}).get("acquisition")
    if acquisition is not None:
        current = _acquisition_evidence(
            profile,
            manifest_path=Path(acquisition["manifest"]["path"]),
            acquisition_lock=Path(acquisition["lock"]["path"]),
            cache_root=Path(acquisition["cache_root"]),
            acquisition_preflight=Path(acquisition["preflight"]["path"]),
            registry=Path(plan["inputs"]["registry"]),
            dataset_registry=Path(plan["inputs"]["dataset_registry"]),
        )
        if current != acquisition:
            raise RuntimeError(
                "campaign acquisition lock or content changed; start a new campaign"
            )
    for child in plan["children"]:
        if child["kind"] == "controller":
            build_controller.assert_matching_provenance(
                build_controller.load_plan(Path(child["root"]))
            )
        else:
            _assert_protocol_provenance(child["plan"])


def _protocol_status(cell: Mapping[str, Any]) -> dict[str, Any]:
    return build_controller.read_json(
        Path(cell["cell_dir"]) / "status.json"
    )


def _write_protocol_status(
    cell: Mapping[str, Any],
    state: str,
    **fields: Any,
) -> None:
    current = _protocol_status(cell)
    build_controller.atomic_write_json(
        Path(cell["cell_dir"]) / "status.json",
        {
            **current,
            **fields,
            "schema_version": SCHEMA_VERSION,
            "case_id": cell["case_id"],
            "state": state,
            "updated_at": build_controller.utc_now(),
        },
    )


def _run_protocol_cell(
    protocol_plan: Mapping[str, Any],
    cell: Mapping[str, Any],
) -> int:
    from . import release_evaluation

    status = _protocol_status(cell)
    attempts = int(status.get("attempts", 0)) + 1
    _write_protocol_status(
        cell,
        "running",
        attempts=attempts,
        errors=[],
        started_at=build_controller.utc_now(),
    )
    work = Path(cell["cell_dir"]) / "work"
    try:
        source = release_evaluation.SourceSpec(
            label="candidate",
            root=Path(protocol_plan["candidate"]["root"]),
            sha=protocol_plan["candidate"]["sha"],
            lifton_executable=Path(
                protocol_plan["candidate"]["lifton_executable"]
            ),
        )
        inputs = release_evaluation.resolve_panel_inputs(
            "e2e",
            cell["dataset"],
            benchmark_registry=Path(protocol_plan["benchmark_registry"]),
            dataset_registry=Path(protocol_plan["dataset_registry"]),
        )
        _output, result = release_evaluation._run_e2e_one(
            source,
            inputs,
            work,
            threads=cell["threads"],
            mode=cell["mode"],
            dataset_registry=Path(protocol_plan["dataset_registry"]),
        )
        record = {
            key: cell[key] for key in (
                "case_id", "kind", "dataset", "repetition",
                "position", "threads", "mode",
            )
        }
        record.update({
            "candidate_sha": protocol_plan["candidate"]["sha"],
            "profile": result["profile"],
            "fingerprints": result["fingerprints"],
            "validation": result["validation"],
            "output_gff": result["output_gff"],
        })
        build_controller.atomic_write_json(Path(cell["result_json"]), record)
        _write_protocol_status(
            cell,
            "success",
            attempts=attempts,
            finished_at=build_controller.utc_now(),
            result=_file_record(Path(cell["result_json"])),
        )
        return 0
    except BaseException as exc:
        _write_protocol_status(
            cell,
            "failed",
            attempts=attempts,
            finished_at=build_controller.utc_now(),
            errors=[f"{type(exc).__name__}: {exc}"],
            traceback=traceback.format_exc(),
        )
        return 1


def _run_protocol_child(child: Mapping[str, Any]) -> int:
    protocol_plan = child["plan"]
    _assert_protocol_provenance(protocol_plan)
    for cell in protocol_plan["cells"]:
        state = _protocol_status(cell).get("state")
        if state == "success":
            continue
        if state == "failed":
            return 1
        if _run_protocol_cell(protocol_plan, cell) != 0:
            return 1
    return 0


def _campaign_state(
    plan: Mapping[str, Any],
    state: str,
    **fields: Any,
) -> None:
    build_controller.atomic_write_json(
        Path(plan["campaign_dir"]) / "campaign_state.json",
        {
            "schema_version": SCHEMA_VERSION,
            "run_id": plan["run_id"],
            "fingerprint": plan["fingerprint"],
            "state": state,
            "updated_at": build_controller.utc_now(),
            **fields,
        },
    )


def run_campaign(campaign_dir: Path) -> int:
    plan = load_campaign_plan(campaign_dir)
    assert_campaign_provenance(plan)
    _campaign_state(plan, "running", current_campaign_id=None)
    for child in plan["children"]:
        _campaign_state(
            plan,
            "running",
            current_campaign_id=child["id"],
        )
        if child["kind"] == "controller":
            child_plan = build_controller.load_plan(Path(child["root"]))
            summary = build_controller.summarize_run(child_plan)
            code = build_controller.status_exit_code(summary)
            if code == build_controller.STATUS_EXIT_FAILED:
                _campaign_state(
                    plan,
                    "failed",
                    current_campaign_id=child["id"],
                    error="child controller has failed cells; explicit retry required",
                )
                return 1
            if code != build_controller.STATUS_EXIT_SUCCESS:
                code = build_controller.scheduler_loop(Path(child["root"]))
        else:
            code = _run_protocol_child(child)
        if code != 0:
            _campaign_state(
                plan,
                "failed",
                current_campaign_id=child["id"],
                error="campaign child failed; no automatic retry was attempted",
            )
            return 1
    _campaign_state(plan, "success", current_campaign_id=None)
    return 0


def _protocol_summary(child: Mapping[str, Any]) -> dict[str, Any]:
    rows = []
    for cell in child["plan"]["cells"]:
        status = _protocol_status(cell)
        rows.append({
            "case_id": cell["case_id"],
            "state": status.get("state", "unknown"),
            "attempts": status.get("attempts", 0),
            "errors": status.get("errors", []),
        })
    counts = dict(Counter(row["state"] for row in rows))
    total = len(rows)
    if counts.get("failed", 0):
        state = "failed"
    elif counts.get("success", 0) == total:
        state = "success"
    elif counts.get("running", 0):
        state = "running"
    else:
        state = "pending"
    return {
        "id": child["id"],
        "kind": "protocol",
        "state": state,
        "counts": counts,
        "cells": rows,
    }


def _controller_truth_errors(child: Mapping[str, Any]) -> list[str]:
    """Apply profile truth policy to sealed paired controller results."""

    from . import release_report

    plan = build_controller.load_plan(Path(child["root"]))
    case = plan.get("campaign_case", {}).get("case", {})
    policy = case.get("truth_policy")
    if policy not in {
        "target_truth_required", "synthetic_exact_required",
    }:
        return []
    errors = []
    for cell in plan["cells"]:
        if build_controller._read_status(cell).get("state") != "success":
            continue
        result_path = Path(cell["artifacts"]["result_json"])
        pair = build_controller.read_json(result_path)
        summaries = {}
        for label in ("candidate", "reference"):
            try:
                summaries[label] = release_report._target_truth_summary(
                    result_path,
                    pair,
                    label,
                )
            except ValueError as exc:
                errors.append(f"{cell['id']}:{label}: {exc}")
        if set(summaries) != {"candidate", "reference"} or any(
            value is None for value in summaries.values()
        ):
            errors.append(f"{cell['id']}: mandatory target truth is missing")
            continue
        expected_id_policy = "ortholog-map"
        unsafe_mapping = [
            label
            for label, summary in summaries.items()
            if (
                (summary.get("parameters") or {}).get("id_policy")
                != expected_id_policy
                or (summary.get("parameters") or {}).get(
                    "mapping_requirement_satisfied"
                ) is not True
            )
        ]
        if unsafe_mapping:
            errors.append(
                f"{cell['id']}: {policy} has unsafe target-truth mapping "
                f"policy for {unsafe_mapping}"
            )
            continue
        if policy == "target_truth_required":
            for level in ("gene", "transcript"):
                candidate = release_report._truth_f1(
                    summaries["candidate"], level,
                )
                reference = release_report._truth_f1(
                    summaries["reference"], level,
                )
                if (
                    candidate is None
                    or candidate < release_report.TARGET_TRUTH_LOCUS_F1_FLOOR
                ):
                    errors.append(
                        f"{cell['id']}: candidate {level} locus F1 "
                        f"{candidate!r} is below "
                        f"{release_report.TARGET_TRUTH_LOCUS_F1_FLOOR}"
                    )
                if (
                    candidate is None
                    or reference is None
                    or candidate - reference
                    < release_report.TARGET_TRUTH_DELTA_FLOOR
                ):
                    errors.append(
                        f"{cell['id']}: {level} locus F1 delta fails "
                        f"{release_report.TARGET_TRUTH_DELTA_FLOOR}"
                    )
            continue
        candidate = summaries["candidate"]
        for group, names in (
            ("gene", ("locus", "strand", "copy")),
            ("transcript", ("locus", "strand", "copy")),
            ("structure", ("intron_chain", "intron", "exon", "CDS")),
        ):
            for name in names:
                f1 = candidate[group][name].get("f1")
                if f1 != 1.0:
                    errors.append(
                        f"{cell['id']}: synthetic {group}.{name} F1 is "
                        f"{f1!r}, expected 1.0"
                    )
        for group in ("gene", "transcript"):
            rate = candidate[group]["copy_count_exact"].get("rate")
            if rate != 1.0:
                errors.append(
                    f"{cell['id']}: synthetic {group} copy-count exact rate "
                    f"is {rate!r}, expected 1.0"
                )
    return errors


def summarize_campaign(plan: Mapping[str, Any]) -> dict[str, Any]:
    children = []
    for child in plan["children"]:
        if child["kind"] == "protocol":
            children.append(_protocol_summary(child))
            continue
        try:
            summary = build_controller.summarize_run(
                build_controller.load_plan(Path(child["root"]))
            )
            state = (
                "success"
                if build_controller.status_exit_code(summary)
                == build_controller.STATUS_EXIT_SUCCESS
                else "failed"
                if build_controller.status_exit_code(summary)
                == build_controller.STATUS_EXIT_FAILED
                else "running"
                if summary["controller_state"] in ACTIVE_STATES
                else "pending"
            )
            children.append({
                "id": child["id"],
                "kind": "controller",
                "state": state,
                "counts": summary["counts"],
                "root": child["root"],
            })
        except (OSError, RuntimeError, ValueError) as exc:
            children.append({
                "id": child["id"],
                "kind": "controller",
                "state": "invalid",
                "error": str(exc),
            })
    counts = dict(Counter(child["state"] for child in children))
    try:
        state = build_controller.read_json(
            Path(plan["campaign_dir"]) / "campaign_state.json"
        ).get("state", "unknown")
    except (OSError, TypeError, ValueError):
        state = "unknown"
    return {
        "schema_version": SCHEMA_VERSION,
        "run_id": plan["run_id"],
        "profile_id": plan["profile"]["id"],
        "profile_digest": plan["profile"]["digest"],
        "fingerprint": plan["fingerprint"],
        "state": state,
        "counts": counts,
        "children": children,
        "updated_at": build_controller.utc_now(),
    }


def status_exit_code(summary: Mapping[str, Any]) -> int:
    counts = summary.get("counts")
    if not isinstance(counts, Mapping) or counts.get("invalid", 0):
        return build_controller.STATUS_EXIT_INVALID
    if counts.get("failed", 0) or summary.get("state") == "failed":
        return build_controller.STATUS_EXIT_FAILED
    if (
        counts.get("success", 0) == len(summary.get("children", []))
        and summary.get("state") in {"success", "finalized"}
    ):
        return build_controller.STATUS_EXIT_SUCCESS
    return build_controller.STATUS_EXIT_INCOMPLETE


def retry_campaign(
    campaign_dir: Path,
    *,
    campaign_ids: Sequence[str] | None = None,
    cells: Sequence[str] | None = None,
) -> list[str]:
    plan = load_campaign_plan(campaign_dir)
    assert_campaign_provenance(plan)
    wanted = set(campaign_ids or [])
    unknown = wanted - {child["id"] for child in plan["children"]}
    if unknown:
        raise ValueError(f"unknown campaign ids: {sorted(unknown)}")
    requested_cells = set(cells or [])
    available_cells = set()
    for child in plan["children"]:
        if wanted and child["id"] not in wanted:
            continue
        if child["kind"] == "controller":
            child_plan = build_controller.load_plan(Path(child["root"]))
            available_cells.update(cell["id"] for cell in child_plan["cells"])
        else:
            available_cells.update(
                cell["case_id"] for cell in child["plan"]["cells"]
            )
    if requested_cells - available_cells:
        raise ValueError(
            f"unknown campaign cells: "
            f"{sorted(requested_cells - available_cells)}"
        )
    reset = []
    for child in plan["children"]:
        if wanted and child["id"] not in wanted:
            continue
        if child["kind"] == "controller":
            selected_cells = sorted(
                requested_cells & {
                    cell["id"]
                    for cell in build_controller.load_plan(
                        Path(child["root"])
                    )["cells"]
                }
            )
            if requested_cells and not selected_cells:
                continue
            for cell_id in build_controller.retry_failed(
                Path(child["root"]),
                selected_cells if requested_cells else None,
            ):
                reset.append(f"{child['id']}:{cell_id}")
            continue
        available = {
            cell["case_id"] for cell in child["plan"]["cells"]
        }
        if requested_cells and not requested_cells & available:
            continue
        for cell in child["plan"]["cells"]:
            status = _protocol_status(cell)
            if (
                status.get("state") != "failed"
                or (requested_cells and cell["case_id"] not in requested_cells)
            ):
                continue
            cell_dir = Path(cell["cell_dir"])
            archive = cell_dir / f"attempt-{int(status.get('attempts', 0)):02d}"
            archive.mkdir()
            for name in ("work", "result.json"):
                source = cell_dir / name
                if source.exists() or source.is_symlink():
                    os.replace(source, archive / name)
            _write_protocol_status(
                cell,
                "pending",
                attempts=int(status.get("attempts", 0)),
                errors=[],
                manual_retries=int(status.get("manual_retries", 0)) + 1,
            )
            reset.append(f"{child['id']}:{cell['case_id']}")
    if reset:
        _campaign_state(plan, "planned", current_campaign_id=None)
    return reset


def reconcile_campaign(
    campaign_dir: Path,
    *,
    deep: bool = True,
) -> dict[str, Any]:
    plan = load_campaign_plan(campaign_dir)
    assert_campaign_provenance(plan)
    errors = []
    protocol_summaries = {}
    for child in plan["children"]:
        if child["kind"] == "controller":
            summary = build_controller.reconcile_run(
                Path(child["root"]),
                deep=deep,
            )
            if build_controller.status_exit_code(summary) != 0:
                errors.append(
                    f"{child['id']}: controller is not completely successful"
                )
            errors.extend(
                f"{child['id']}: {error}"
                for error in _controller_truth_errors(child)
            )
            continue
        records = []
        for cell in child["plan"]["cells"]:
            status = _protocol_status(cell)
            if status.get("state") != "success":
                errors.append(
                    f"{child['id']}:{cell['case_id']} is "
                    f"{status.get('state', 'unknown')}"
                )
                continue
            record_path = Path(cell["result_json"])
            evidence = status.get("result")
            try:
                observed_record = _file_record(record_path)
            except OSError as exc:
                errors.append(
                    f"{child['id']}:{cell['case_id']} result is "
                    f"unavailable: {exc}"
                )
                continue
            if not isinstance(evidence, Mapping) or observed_record != {
                "path": evidence.get("path"),
                "size": evidence.get("size"),
                "sha256": evidence.get("sha256"),
            }:
                errors.append(
                    f"{child['id']}:{cell['case_id']} result hash changed"
                )
                continue
            records.append(build_controller.read_json(record_path))
        if len(records) == len(child["plan"]["cells"]):
            try:
                protocol_summaries[child["id"]] = (
                    protocol_analysis.summarize_protocol(
                        records,
                        kind=child["plan"]["protocol_kind"],
                    )
                )
            except ValueError as exc:
                errors.append(f"{child['id']}: {exc}")
    result = {
        "schema_version": SCHEMA_VERSION,
        "run_id": plan["run_id"],
        "deep": deep,
        "passed": not errors,
        "errors": errors,
        "protocol_summaries": protocol_summaries,
        "summary": summarize_campaign(plan),
        "updated_at": build_controller.utc_now(),
    }
    build_controller.atomic_write_json(
        Path(campaign_dir) / "campaign_reconciliation.json",
        result,
    )
    if errors:
        _campaign_state(
            plan,
            "failed",
            current_campaign_id=None,
            error="campaign reconciliation failed",
        )
    return result


def finalize_campaign(
    campaign_dir: Path,
    *,
    output_dir: Path | None = None,
) -> dict[str, Any]:
    plan = load_campaign_plan(campaign_dir)
    if not plan.get("full_profile"):
        raise ValueError("a partial profile campaign cannot be finalized for release")
    reconciliation = reconcile_campaign(campaign_dir, deep=True)
    if not reconciliation["passed"]:
        raise ValueError("campaign reconciliation failed; release bundle not written")
    profile = campaign_profiles.load_profile(
        plan["profile"]["id"],
        registry=Path(plan["profile"]["registry"]["path"]),
    )
    release_spec = campaign_profiles.campaign_spec(profile)
    destination = Path(
        output_dir or Path(campaign_dir) / "publication"
    ).resolve()
    destination.mkdir(parents=True, exist_ok=True)
    build_controller.atomic_write_json(
        destination / "campaign_spec.json",
        release_spec,
    )
    evidence = {
        "campaign_plan": _file_record(
            Path(campaign_dir) / "campaign_plan.json"
        ),
        "campaign_reconciliation": _file_record(
            Path(campaign_dir) / "campaign_reconciliation.json"
        ),
        "campaign_spec": _file_record(destination / "campaign_spec.json"),
        "children": [],
    }
    release_roots = []
    release_case_ids = {
        case["id"] for case in campaign_profiles.release_cases(profile)
    }
    for child in plan["children"]:
        root = Path(child["root"])
        row = {
            "id": child["id"],
            "kind": child["kind"],
            "root": str(root),
            "plan": _file_record(root / "plan.json"),
        }
        if child["kind"] == "controller":
            row["summary"] = _file_record(root / "summary.json")
            row["reconciled_results"] = _file_record(
                root / "reconciled_results.json"
            )
            if child["id"] in release_case_ids:
                release_roots.append(str(root))
        else:
            row["results"] = [
                _file_record(Path(cell["result_json"]))
                for cell in child["plan"]["cells"]
            ]
        evidence["children"].append(row)
    bundle = {
        "schema_version": SCHEMA_VERSION,
        "run_id": plan["run_id"],
        "candidate_sha": plan["sources"]["candidate"]["sha"],
        "reference_sha": plan["sources"]["reference"]["sha"],
        "profile_id": plan["profile"]["id"],
        "profile_digest": plan["profile"]["digest"],
        "campaign_fingerprint": plan["fingerprint"],
        "release_report_roots": release_roots,
        "evidence": evidence,
    }
    bundle["fingerprint"] = _hash_without_fingerprint(bundle)
    build_controller.atomic_write_json(
        destination / "campaign_bundle.json",
        bundle,
    )
    _campaign_state(
        plan,
        "finalized",
        current_campaign_id=None,
        publication_bundle=_file_record(
            destination / "campaign_bundle.json"
        ),
    )
    return bundle


def session_name(plan: Mapping[str, Any]) -> str:
    return build_controller.safe_name(
        f"lifton-campaign-{plan['run_id']}-{plan['fingerprint'][:10]}-ctl",
        limit=100,
    )


def start_tmux(plan: Mapping[str, Any]) -> str:
    name = session_name(plan)
    if build_controller.tmux_has_session(name):
        return name
    build_controller._tmux_new_session(name, build_controller._python_worker_command(
        "-m",
        "benchmarks.compare.campaign_orchestrator",
        "_worker",
        "--campaign-dir",
        plan["campaign_dir"],
    ))
    return name


def _print_summary(summary: Mapping[str, Any]) -> None:
    print(
        f"campaign: {summary['run_id']}  profile: {summary['profile_id']}  "
        f"state: {summary['state']}"
    )
    print("states: " + " ".join(
        f"{key}={value}" for key, value in sorted(summary["counts"].items())
    ))
    for child in summary["children"]:
        print(f"  {child['state']:<10} {child['id']} ({child['kind']})")


def _campaign_argument(args: argparse.Namespace) -> Path:
    if getattr(args, "campaign_dir", None):
        return Path(args.campaign_dir).resolve()
    return _campaign_dir(Path(args.campaigns_root), args.run_id)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="action", required=True)

    plan = subparsers.add_parser("plan", help="freeze and initialize every child plan")
    plan.add_argument("--run-id", required=True)
    plan.add_argument("--campaign-profile", default="canonical-v2")
    plan.add_argument("--campaign-id", action="append")
    plan.add_argument("--profile-registry", default=str(
        campaign_profiles.DEFAULT_PROFILE_REGISTRY
    ))
    plan.add_argument("--campaigns-root", default=str(DEFAULT_CAMPAIGNS_ROOT))
    plan.add_argument("--candidate-root", required=True)
    plan.add_argument("--candidate-sha", required=True)
    plan.add_argument("--reference-root", required=True)
    plan.add_argument("--reference-sha", required=True)
    plan.add_argument("--lifton-executable", required=True)
    plan.add_argument("--registry", default=str(build_controller.DEFAULT_REGISTRY))
    plan.add_argument(
        "--dataset-registry",
        default=str(build_controller.DEFAULT_DATASET_REGISTRY),
    )
    plan.add_argument("--baseline", default=str(build_controller.DEFAULT_BASELINE))
    plan.add_argument("--dataset-manifest")
    plan.add_argument("--acquisition-lock")
    plan.add_argument("--dataset-cache-root")
    plan.add_argument("--acquisition-preflight")
    build_controller._add_policy_arguments(plan)

    for action in ("start", "status", "reconcile", "finalize"):
        command = subparsers.add_parser(action)
        command.add_argument("run_id", nargs="?")
        command.add_argument("--campaign-dir")
        command.add_argument(
            "--campaigns-root", default=str(DEFAULT_CAMPAIGNS_ROOT),
        )
        if action == "start":
            command.add_argument("--foreground", action="store_true")
        elif action == "status":
            command.add_argument("--json", action="store_true")
        elif action == "reconcile":
            command.add_argument("--shallow", action="store_true")
            command.add_argument("--json", action="store_true")
        elif action == "finalize":
            command.add_argument("--output-dir")

    retry = subparsers.add_parser("retry")
    retry.add_argument("run_id", nargs="?")
    retry.add_argument("--campaign-dir")
    retry.add_argument("--campaigns-root", default=str(DEFAULT_CAMPAIGNS_ROOT))
    retry.add_argument("--campaign-id", action="append")
    retry.add_argument("--cells", nargs="+")

    worker = subparsers.add_parser("_worker")
    worker.add_argument("--campaign-dir", required=True)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    if args.action == "plan":
        try:
            campaign_dir, plan, children = create_campaign_plan(
                run_id=args.run_id,
                profile_id=args.campaign_profile,
                campaigns_root=Path(args.campaigns_root),
                profile_registry=Path(args.profile_registry),
                candidate_root=Path(args.candidate_root),
                candidate_sha=args.candidate_sha,
                reference_root=Path(args.reference_root),
                reference_sha=args.reference_sha,
                lifton_executable=Path(args.lifton_executable),
                registry=Path(args.registry),
                dataset_registry=Path(args.dataset_registry),
                baseline=Path(args.baseline),
                base_policy=build_controller._policy_from_args(args),
                campaign_ids=args.campaign_id,
                dataset_manifest=(
                    Path(args.dataset_manifest)
                    if args.dataset_manifest else None
                ),
                acquisition_lock=(
                    Path(args.acquisition_lock)
                    if args.acquisition_lock else None
                ),
                dataset_cache_root=(
                    Path(args.dataset_cache_root)
                    if args.dataset_cache_root else None
                ),
                acquisition_preflight=(
                    Path(args.acquisition_preflight)
                    if args.acquisition_preflight else None
                ),
            )
            initialize_campaign(campaign_dir, plan, children)
        except (OSError, RuntimeError, ValueError) as exc:
            parser.error(str(exc))
        print(f"campaign directory: {campaign_dir}")
        print(f"campaign fingerprint: {plan['fingerprint']}")
        print(json.dumps(plan["profile"]["counts"], indent=2, sort_keys=True))
        return 0

    if args.action == "_worker":
        return run_campaign(Path(args.campaign_dir))
    if not getattr(args, "campaign_dir", None) and not getattr(
        args, "run_id", None
    ):
        parser.error("provide run_id or --campaign-dir")
    campaign_dir = _campaign_argument(args)
    try:
        plan = load_campaign_plan(campaign_dir)
        if args.action == "start":
            assert_campaign_provenance(plan)
            if args.foreground:
                return run_campaign(campaign_dir)
            name = start_tmux(plan)
            print(f"campaign tmux session: {name}")
            return 0
        if args.action == "status":
            summary = summarize_campaign(plan)
            if args.json:
                print(json.dumps(summary, indent=2, sort_keys=True))
            else:
                _print_summary(summary)
            return status_exit_code(summary)
        if args.action == "retry":
            reset = retry_campaign(
                campaign_dir,
                campaign_ids=args.campaign_id,
                cells=args.cells,
            )
            print("reset: " + (", ".join(reset) if reset else "no failed cells"))
            return 0
        if args.action == "reconcile":
            result = reconcile_campaign(
                campaign_dir,
                deep=not args.shallow,
            )
            if args.json:
                print(json.dumps(result, indent=2, sort_keys=True))
            else:
                print(
                    f"reconciliation: {'PASS' if result['passed'] else 'FAIL'}"
                )
                for error in result["errors"]:
                    print(f"  {error}")
            return 0 if result["passed"] else 1
        if args.action == "finalize":
            bundle = finalize_campaign(
                campaign_dir,
                output_dir=(
                    Path(args.output_dir) if args.output_dir else None
                ),
            )
            print(json.dumps(bundle, indent=2, sort_keys=True))
            return 0
    except (OSError, RuntimeError, ValueError) as exc:
        parser.error(str(exc))
    parser.error(f"unsupported action: {args.action}")
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
