#!/usr/bin/env python3
"""Validate and plan acquisition for canonical-v2 benchmark datasets.

This module never downloads data.  ``acquire --dry-run`` emits an exact,
content-addressed plan: each immutable source request has a SHA-256 staging key,
and acquired bytes must ultimately be stored below
``sha256/<prefix>/<content-sha256>/``.  ``verify-lock`` independently hashes
those stored files once a separate acquisition worker has populated them.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
import sys
from pathlib import Path
from typing import Any, Mapping, Sequence


HERE = Path(__file__).resolve().parent
DEFAULT_MANIFEST = HERE / "manifests" / "canonical_v2_datasets.json"
SCHEMA_VERSION = 1
LOCK_SCHEMA_VERSION = 1
SHA256_RE = re.compile(r"^[0-9a-f]{64}$")
ID_RE = re.compile(r"^[a-z][a-z0-9_]*$")
ACCESSION_RE = re.compile(r"^(?:GC[AF]_\d+\.\d+|PRJ[A-Z]{2}\d+)$")
URL_RE = re.compile(r"^https://[^?#]+$")
VALID_TRANSPORTS = {"ncbi_datasets", "https_files", "existing_registry"}
VALID_KINDS = {"biological", "synthetic", "protocol"}
VALID_PANELS = {
    "subset", "full", "truth", "cross_e2e", "synthetic", "protocol",
}
EXPECTED_SCENARIO_IDS = (
    "v2_truth_human_grch38_chm13",
    "v2_dialect_ensembl116_gtf",
    "v2_dialect_flybase_dmel_dere",
    "v2_dialect_wormbase_celegans_cbriggsae",
    "v2_truth_soybean_w82_lee",
    "v2_deep_zebrafish_xenopus",
    "v2_deep_tomato_rice",
    "v2_truth_rat_mouse_e2e",
    "v2_synth_chr22_fragmented",
    "v2_synth_chr22_sv",
    "v2_protocol_thread_scaling_bee",
    "v2_protocol_io_modes_arabidopsis",
)
TOP_LEVEL_KEYS = {
    "schema_version", "manifest_id", "profile_id", "description",
    "content_addressing", "sources", "scenarios",
}
SOURCE_KEYS = {
    "id", "provider", "release", "transport", "identity", "acquisition",
    "files",
}
FILE_KEYS = {"role", "locator", "expected_sha256", "pin_state"}
SCENARIO_KEYS = {"id", "kind", "panels", "description", "inputs", "design"}
FRAGMENTED_SCENARIO_ID = "v2_synth_chr22_fragmented"
SV_SCENARIO_ID = "v2_synth_chr22_sv"
SYNTHETIC_SEED_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]{0,127}$")
FRAGMENT_BREAKPOINT_POLICY = (
    "Cuts are intergenic; every source model remains exactly scorable after "
    "deterministic coordinate remapping."
)
SV_BREAKPOINT_POLICY = (
    "Models crossing a structural breakpoint are reported separately and "
    "excluded from exact-recovery denominators; retained, deleted, inverted, "
    "and duplicated models have explicit expected copy and strand states."
)


class ManifestError(ValueError):
    """Raised when a dataset manifest or acquisition lock is invalid."""


def canonical_json(value: Any) -> str:
    try:
        return json.dumps(
            value, sort_keys=True, separators=(",", ":"), ensure_ascii=True,
            allow_nan=False,
        )
    except (TypeError, ValueError) as exc:
        raise ManifestError(f"value is not finite canonical JSON: {exc}") from exc


def canonical_sha256(value: Any) -> str:
    return hashlib.sha256(canonical_json(value).encode("utf-8")).hexdigest()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_json(path: Path) -> Any:
    try:
        return json.loads(Path(path).read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise ManifestError(f"cannot load JSON {path}: {exc}") from exc


def _require_exact_keys(
    value: Any,
    keys: set[str],
    *,
    label: str,
) -> Mapping[str, Any]:
    if not isinstance(value, dict) or set(value) != keys:
        observed = sorted(value) if isinstance(value, dict) else type(value).__name__
        raise ManifestError(
            f"{label} keys must be exactly {sorted(keys)}; observed {observed}"
        )
    return value


def _require_id(value: Any, *, label: str) -> str:
    if not isinstance(value, str) or ID_RE.fullmatch(value) is None:
        raise ManifestError(f"{label} must match {ID_RE.pattern}")
    return value


def _validate_sha_pin(file_record: Mapping[str, Any], label: str) -> None:
    digest = file_record["expected_sha256"]
    pin_state = file_record["pin_state"]
    if pin_state == "identity_pinned_bytes_pending":
        if digest is not None:
            raise ManifestError(f"{label}: pending bytes must have null SHA-256")
    elif pin_state == "sha256_pinned":
        if not isinstance(digest, str) or SHA256_RE.fullmatch(digest) is None:
            raise ManifestError(f"{label}: pinned bytes require lowercase SHA-256")
    elif pin_state == "runtime_fingerprint_required":
        if digest is not None:
            raise ManifestError(f"{label}: runtime fingerprint must have null SHA")
    else:
        raise ManifestError(f"{label}: invalid pin_state {pin_state!r}")


def _validate_ncbi_source(source: Mapping[str, Any], label: str) -> None:
    identity = _require_exact_keys(
        source["identity"], {"accession", "assembly"}, label=f"{label}.identity"
    )
    accession = identity["accession"]
    if not isinstance(accession, str) or ACCESSION_RE.fullmatch(accession) is None:
        raise ManifestError(f"{label}: invalid NCBI accession {accession!r}")
    if not accession.startswith(("GCF_", "GCA_")):
        raise ManifestError(f"{label}: genome package requires GCF_/GCA_ accession")
    acquisition = _require_exact_keys(
        source["acquisition"], {"include", "filename"},
        label=f"{label}.acquisition",
    )
    include = acquisition["include"]
    if (
        not isinstance(include, list) or not include
        or include != sorted(set(include))
        or any(item not in {"genome", "gff3", "protein"} for item in include)
    ):
        raise ManifestError(f"{label}: include must be a sorted unique NCBI set")
    if (
        not isinstance(acquisition["filename"], str)
        or not acquisition["filename"].endswith(".zip")
        or "/" in acquisition["filename"]
    ):
        raise ManifestError(f"{label}: invalid NCBI archive filename")
    prefix = f"ncbi_dataset/data/{accession}/"
    for file_record in source["files"]:
        locator = _require_exact_keys(
            file_record["locator"], {"archive_member"},
            label=f"{label}.{file_record['role']}.locator",
        )
        member = locator["archive_member"]
        if not isinstance(member, str) or not member.startswith(prefix):
            raise ManifestError(
                f"{label}.{file_record['role']}: archive member must start {prefix}"
            )


def _validate_https_source(source: Mapping[str, Any], label: str) -> None:
    if not isinstance(source["identity"], dict) or not source["identity"]:
        raise ManifestError(f"{label}.identity must be a non-empty object")
    _require_exact_keys(
        source["acquisition"], set(), label=f"{label}.acquisition"
    )
    for file_record in source["files"]:
        locator = _require_exact_keys(
            file_record["locator"], {"url", "filename"},
            label=f"{label}.{file_record['role']}.locator",
        )
        if not isinstance(locator["url"], str) or URL_RE.fullmatch(locator["url"]) is None:
            raise ManifestError(
                f"{label}.{file_record['role']}: URL must be fixed HTTPS without query"
            )
        filename = locator["filename"]
        if (
            not isinstance(filename, str) or not filename
            or "/" in filename or filename in {".", ".."}
        ):
            raise ManifestError(
                f"{label}.{file_record['role']}: invalid destination filename"
            )


def _validate_registry_source(source: Mapping[str, Any], label: str) -> None:
    identity = _require_exact_keys(
        source["identity"], {"registry_path", "benchmark_id"},
        label=f"{label}.identity",
    )
    if identity["registry_path"] != "benchmarks/compare/benchmarks.json":
        raise ManifestError(f"{label}: existing source must use canonical registry")
    _require_id(identity["benchmark_id"], label=f"{label}.benchmark_id")
    _require_exact_keys(
        source["acquisition"], set(), label=f"{label}.acquisition"
    )
    for file_record in source["files"]:
        locator = _require_exact_keys(
            file_record["locator"], {"registry_key"},
            label=f"{label}.{file_record['role']}.locator",
        )
        if locator["registry_key"] not in {
            "ref_genome", "ref_gff", "ref_proteins", "tgt_genome",
        }:
            raise ManifestError(
                f"{label}.{file_record['role']}: invalid registry key"
            )
        if file_record["pin_state"] != "runtime_fingerprint_required":
            raise ManifestError(
                f"{label}.{file_record['role']}: registry bytes require runtime pin"
            )


def _validate_source(source: Any, index: int) -> tuple[str, set[str]]:
    source = _require_exact_keys(source, SOURCE_KEYS, label=f"sources[{index}]")
    source_id = _require_id(source["id"], label=f"sources[{index}].id")
    label = f"source {source_id}"
    if (
        not isinstance(source["provider"], str) or not source["provider"]
        or not isinstance(source["release"], str) or not source["release"]
    ):
        raise ManifestError(f"{label}: provider and release must be non-empty")
    transport = source["transport"]
    if transport not in VALID_TRANSPORTS:
        raise ManifestError(f"{label}: invalid transport {transport!r}")
    if not isinstance(source["files"], list) or not source["files"]:
        raise ManifestError(f"{label}: files must be a non-empty list")
    roles: list[str] = []
    for file_index, file_record in enumerate(source["files"]):
        file_record = _require_exact_keys(
            file_record, FILE_KEYS, label=f"{label}.files[{file_index}]"
        )
        role = _require_id(
            file_record["role"], label=f"{label}.files[{file_index}].role"
        )
        roles.append(role)
        _validate_sha_pin(file_record, f"{label}.{role}")
    if roles != sorted(set(roles)):
        raise ManifestError(f"{label}: file roles must be sorted and unique")
    if transport == "ncbi_datasets":
        _validate_ncbi_source(source, label)
    elif transport == "https_files":
        _validate_https_source(source, label)
    else:
        _validate_registry_source(source, label)
    canonical_json(source)
    return source_id, set(roles)


def _validate_biological_design(
    scenario: Mapping[str, Any],
    label: str,
) -> None:
    design = _require_exact_keys(
        scenario["design"],
        {
            "reference_label", "target_label", "input_format",
            "reference_scope", "target_scope", "truth_policy",
        },
        label=f"{label}.design",
    )
    if design["input_format"] not in {"GFF3", "GTF"}:
        raise ManifestError(f"{label}: biological input format must be GFF3/GTF")
    if design["truth_policy"] != "ortholog_aware_independent_target_annotation":
        raise ManifestError(f"{label}: target truth policy must be ortholog-aware")
    required = {
        "reference_genome", "reference_annotation",
        "target_genome", "target_truth",
    }
    missing = sorted(required - set(scenario["inputs"]))
    if missing:
        raise ManifestError(f"{label}: missing biological inputs {missing}")


def _validate_synthetic_design(
    scenario: Mapping[str, Any],
    label: str,
) -> None:
    scenario_id = scenario["id"]
    if scenario_id == FRAGMENTED_SCENARIO_ID:
        design_keys = {
            "base_benchmark", "chromosome", "coordinate_system",
            "operations", "breakpoint_policy",
        }
        expected_policy = FRAGMENT_BREAKPOINT_POLICY
    elif scenario_id == SV_SCENARIO_ID:
        design_keys = {
            "base_benchmark", "chromosome", "seed", "coordinate_system",
            "operations", "breakpoint_policy",
        }
        expected_policy = SV_BREAKPOINT_POLICY
    else:
        raise ManifestError(f"{label}: unsupported synthetic scenario")

    design = _require_exact_keys(
        scenario["design"],
        design_keys,
        label=f"{label}.design",
    )
    if design["base_benchmark"] != "human_mane" or design["chromosome"] != "chr22":
        raise ManifestError(f"{label}: synthetic truth must use human_mane chr22")
    if (
        scenario_id == SV_SCENARIO_ID
        and (
            not isinstance(design["seed"], str)
            or SYNTHETIC_SEED_RE.fullmatch(design["seed"]) is None
        )
    ):
        raise ManifestError(f"{label}: synthetic seed must be a safe non-empty string")
    if design["coordinate_system"] != "1-based-inclusive":
        raise ManifestError(f"{label}: coordinate system must be 1-based-inclusive")
    if design["breakpoint_policy"] != expected_policy:
        raise ManifestError(f"{label}: breakpoint policy is not canonical")

    operations = design["operations"]
    expected_count = 1 if scenario_id == FRAGMENTED_SCENARIO_ID else 4
    if not isinstance(operations, list) or len(operations) != expected_count:
        raise ManifestError(
            f"{label}: operations must match the exact {scenario_id} schema"
        )
    if scenario_id == FRAGMENTED_SCENARIO_ID:
        operation = _require_exact_keys(
            operations[0],
            {"type", "cut_selection", "cut_after"},
            label=f"{label}.design.operations[0]",
        )
        if operation["type"] != "fragment":
            raise ManifestError(f"{label}: fragment operation type is unsupported")
        if operation["cut_selection"] != "fixed_intergenic":
            raise ManifestError(f"{label}: fragment cut selection is unsupported")
        cuts = operation["cut_after"]
        if (
            not isinstance(cuts, list)
            or len(cuts) != 8
            or any(
                not isinstance(cut, int) or isinstance(cut, bool) or cut <= 0
                for cut in cuts
            )
            or any(left >= right for left, right in zip(cuts, cuts[1:]))
        ):
            raise ManifestError(
                f"{label}: fragment cuts must be eight strictly increasing "
                "positive integers"
            )
    else:
        operation_schemas = (
            ("delete", {"type", "start", "end"}),
            ("invert", {"type", "start", "end"}),
            ("tandem_duplicate", {"type", "start", "end", "copies"}),
            ("insert", {"type", "after", "length", "generator"}),
        )
        for index, (operation, (operation_type, keys)) in enumerate(
            zip(operations, operation_schemas)
        ):
            operation = _require_exact_keys(
                operation,
                keys,
                label=f"{label}.design.operations[{index}]",
            )
            if operation["type"] != operation_type:
                raise ManifestError(
                    f"{label}: operation {index} must be {operation_type}"
                )
        intervals = [
            (operation["start"], operation["end"])
            for operation in operations[:3]
        ]
        if any(
            not isinstance(value, int) or isinstance(value, bool) or value <= 0
            for interval in intervals
            for value in interval
        ) or any(start > end for start, end in intervals):
            raise ManifestError(f"{label}: SV intervals must be positive and ordered")
        insert_after = operations[3]["after"]
        insertion_length = operations[3]["length"]
        if (
            not isinstance(insert_after, int)
            or isinstance(insert_after, bool)
            or insert_after <= 0
            or not isinstance(insertion_length, int)
            or isinstance(insertion_length, bool)
            or insertion_length <= 0
        ):
            raise ManifestError(
                f"{label}: insertion position and length must be positive integers"
            )
        if not (
            intervals[0][1] < intervals[1][0]
            and intervals[1][1] < intervals[2][0]
            and intervals[2][1] < insert_after
        ):
            raise ManifestError(
                f"{label}: SV operations must be ordered and non-overlapping"
            )
        copies = operations[2]["copies"]
        if not isinstance(copies, int) or isinstance(copies, bool) or copies != 2:
            raise ManifestError(f"{label}: tandem duplication must have two copies")
        if operations[3]["generator"] != "sha256_counter_dna_v1":
            raise ManifestError(f"{label}: insertion generator is unsupported")

    expected_inputs = {"reference_genome", "reference_annotation"}
    if set(scenario["inputs"]) != expected_inputs:
        raise ManifestError(
            f"{label}: synthetic inputs must be exactly {sorted(expected_inputs)}"
        )


def synthetic_builder_kwargs(scenario: Mapping[str, Any]) -> dict[str, Any]:
    """Translate a validated synthetic manifest design into builder kwargs."""

    if not isinstance(scenario, dict) or not isinstance(scenario.get("id"), str):
        raise ManifestError("synthetic scenario must be an object with an ID")
    label = f"scenario {scenario['id']}"
    _validate_synthetic_design(scenario, label)
    design = scenario["design"]
    operations = design["operations"]
    if scenario["id"] == FRAGMENTED_SCENARIO_ID:
        return {
            "source_seqid": design["chromosome"],
            "cuts": tuple(operations[0]["cut_after"]),
        }
    return {
        "source_seqid": design["chromosome"],
        "deletion": (operations[0]["start"], operations[0]["end"]),
        "inversion": (operations[1]["start"], operations[1]["end"]),
        "duplication": (operations[2]["start"], operations[2]["end"]),
        "insert_after": operations[3]["after"],
        "insertion_length": operations[3]["length"],
        "seed": design["seed"],
    }


def _deterministic_dna_sha256(length: int, seed: str) -> str:
    """Hash the sequence specified by the sha256_counter_dna_v1 generator."""

    alphabet = "ACGT"
    sequence_digest = hashlib.sha256()
    emitted = 0
    counter = 0
    while emitted < length:
        digest = hashlib.sha256(
            seed.encode("utf-8") + b"\0" + counter.to_bytes(8, "big")
        ).digest()
        bases = "".join(
            alphabet[index]
            for byte in digest
            for index in (
                (byte >> 6) & 3,
                (byte >> 4) & 3,
                (byte >> 2) & 3,
                byte & 3,
            )
        )
        chunk = bases[:length - emitted]
        sequence_digest.update(chunk.encode("ascii"))
        emitted += len(chunk)
        counter += 1
    return sequence_digest.hexdigest()


def expected_synthetic_transform(scenario: Mapping[str, Any]) -> dict[str, Any]:
    """Return the exact transform metadata a synthetic builder must publish."""

    kwargs = synthetic_builder_kwargs(scenario)
    scenario_id = scenario["id"]
    if scenario_id == FRAGMENTED_SCENARIO_ID:
        return {
            "id": scenario_id,
            "kind": "fragmentation",
            "cuts_after_source_coordinate": list(kwargs["cuts"]),
        }
    deletion = kwargs["deletion"]
    return {
        "id": scenario_id,
        "kind": "structural_variation",
        "coordinate_basis": "1-based-inclusive-source-simultaneous",
        "deletion": list(deletion),
        "deleted_length": deletion[1] - deletion[0] + 1,
        "inversion": list(kwargs["inversion"]),
        "tandem_duplication": list(kwargs["duplication"]),
        "insert_after": kwargs["insert_after"],
        "insertion_length": kwargs["insertion_length"],
        "insertion_sha256": _deterministic_dna_sha256(
            kwargs["insertion_length"], kwargs["seed"]
        ),
        "seed": kwargs["seed"],
    }


def _validate_protocol_design(
    scenario: Mapping[str, Any],
    label: str,
) -> None:
    design = _require_exact_keys(
        scenario["design"],
        {
            "benchmark_id", "variable", "values", "repetitions",
            "order_design", "correctness_requirement",
        },
        label=f"{label}.design",
    )
    _require_id(design["benchmark_id"], label=f"{label}.benchmark_id")
    if (
        not isinstance(design["values"], list) or not design["values"]
        or len(design["values"]) != len({
            canonical_json(item) for item in design["values"]
        })
    ):
        raise ManifestError(f"{label}: protocol values must be non-empty/unique")
    if not isinstance(design["repetitions"], int) or design["repetitions"] < 1:
        raise ManifestError(f"{label}: repetitions must be positive")
    if design["order_design"] != "fixed_balanced":
        raise ManifestError(f"{label}: protocol order must be fixed_balanced")
    if design["correctness_requirement"] != "byte_and_semantic_identity":
        raise ManifestError(f"{label}: protocol correctness gate is not strict")
    required = {"reference_genome", "reference_annotation", "target_genome"}
    missing = sorted(required - set(scenario["inputs"]))
    if missing:
        raise ManifestError(f"{label}: missing protocol inputs {missing}")


def _validate_scenario(
    scenario: Any,
    index: int,
    source_roles: Mapping[str, set[str]],
) -> str:
    scenario = _require_exact_keys(
        scenario, SCENARIO_KEYS, label=f"scenarios[{index}]"
    )
    scenario_id = _require_id(scenario["id"], label=f"scenarios[{index}].id")
    label = f"scenario {scenario_id}"
    kind = scenario["kind"]
    if kind not in VALID_KINDS:
        raise ManifestError(f"{label}: invalid kind {kind!r}")
    panels = scenario["panels"]
    if (
        not isinstance(panels, list) or not panels
        or len(panels) != len(set(panels))
        or any(panel not in VALID_PANELS for panel in panels)
    ):
        raise ManifestError(f"{label}: panels must be non-empty and unique")
    if not isinstance(scenario["description"], str) or not scenario["description"]:
        raise ManifestError(f"{label}: description must be non-empty")
    if not isinstance(scenario["inputs"], dict) or not scenario["inputs"]:
        raise ManifestError(f"{label}: inputs must be a non-empty object")
    for role, binding in scenario["inputs"].items():
        _require_id(role, label=f"{label}.input role")
        if not isinstance(binding, str) or binding.count(":") != 1:
            raise ManifestError(f"{label}.{role}: binding must be source_id:file_role")
        source_id, file_role = binding.split(":", 1)
        if source_id not in source_roles or file_role not in source_roles[source_id]:
            raise ManifestError(f"{label}.{role}: unknown binding {binding!r}")
    if kind == "biological":
        _validate_biological_design(scenario, label)
    elif kind == "synthetic":
        _validate_synthetic_design(scenario, label)
    else:
        _validate_protocol_design(scenario, label)
    canonical_json(scenario)
    return scenario_id


def validate_manifest(document: Any) -> dict[str, Any]:
    document = dict(_require_exact_keys(
        document, TOP_LEVEL_KEYS, label="dataset manifest"
    ))
    if document["schema_version"] != SCHEMA_VERSION:
        raise ManifestError(
            f"unsupported manifest schema {document['schema_version']!r}"
        )
    if document["manifest_id"] != "lifton_canonical_v2_datasets":
        raise ManifestError("unexpected manifest_id")
    if document["profile_id"] != "canonical-v2":
        raise ManifestError("unexpected profile_id")
    if not isinstance(document["description"], str) or not document["description"]:
        raise ManifestError("description must be non-empty")
    addressing = _require_exact_keys(
        document["content_addressing"],
        {
            "algorithm", "request_cache_layout", "content_cache_layout",
            "lock_schema_version", "byte_pinning_policy",
        },
        label="content_addressing",
    )
    if addressing != {
        "algorithm": "sha256",
        "request_cache_layout": "requests/{request_sha256}",
        "content_cache_layout": "sha256/{content_sha256_prefix}/{content_sha256}",
        "lock_schema_version": 1,
        "byte_pinning_policy": (
            "Identity is fixed here; acquisition locks the downloaded bytes "
            "before any benchmark may run."
        ),
    }:
        raise ManifestError("content_addressing policy is not canonical")

    if not isinstance(document["sources"], list) or not document["sources"]:
        raise ManifestError("sources must be a non-empty list")
    source_rows = [
        _validate_source(source, index)
        for index, source in enumerate(document["sources"])
    ]
    source_ids = [row[0] for row in source_rows]
    if source_ids != sorted(set(source_ids)):
        raise ManifestError("source IDs must be sorted and unique")
    source_roles = {source_id: roles for source_id, roles in source_rows}

    if not isinstance(document["scenarios"], list):
        raise ManifestError("scenarios must be a list")
    scenario_ids = tuple(
        _validate_scenario(scenario, index, source_roles)
        for index, scenario in enumerate(document["scenarios"])
    )
    if scenario_ids != EXPECTED_SCENARIO_IDS:
        raise ManifestError(
            "scenario IDs/order differ from the approved canonical-v2 design"
        )
    canonical_json(document)
    return document


def source_request(source: Mapping[str, Any]) -> dict[str, Any]:
    """Return immutable request identity, excluding eventual content hashes."""

    return {
        "source_id": source["id"],
        "provider": source["provider"],
        "release": source["release"],
        "transport": source["transport"],
        "identity": source["identity"],
        "acquisition": source["acquisition"],
        "files": [
            {"role": item["role"], "locator": item["locator"]}
            for item in source["files"]
        ],
    }


def source_request_sha256(source: Mapping[str, Any]) -> str:
    return canonical_sha256(source_request(source))


def build_acquisition_plan(
    document: Mapping[str, Any],
    cache_root: Path,
) -> dict[str, Any]:
    manifest = validate_manifest(document)
    cache_root = Path(cache_root)
    steps = []
    for source in manifest["sources"]:
        if source["transport"] == "existing_registry":
            continue
        request_sha = source_request_sha256(source)
        staging = cache_root / "requests" / request_sha
        commands: list[list[str]] = []
        if source["transport"] == "ncbi_datasets":
            acquisition = source["acquisition"]
            commands.append([
                "datasets", "download", "genome", "accession",
                source["identity"]["accession"],
                "--include", ",".join(acquisition["include"]),
                "--filename", str(staging / acquisition["filename"]),
            ])
        else:
            for item in source["files"]:
                locator = item["locator"]
                commands.append([
                    "curl", "--fail", "--location", "--retry", "3",
                    "--output", str(staging / locator["filename"]),
                    locator["url"],
                ])
        steps.append({
            "source_id": source["id"],
            "transport": source["transport"],
            "request_sha256": request_sha,
            "staging_directory": str(staging),
            "commands": commands,
            "content_destination_template": str(
                cache_root / "sha256" / "{content_sha256_prefix}"
                / "{content_sha256}" / "{filename}"
            ),
            "expected_file_roles": [item["role"] for item in source["files"]],
        })
    plan = {
        "schema_version": 1,
        "dry_run": True,
        "manifest_id": manifest["manifest_id"],
        "manifest_sha256": canonical_sha256(manifest),
        "cache_root": str(cache_root),
        "remote_source_count": len(steps),
        "steps": steps,
    }
    plan["plan_sha256"] = canonical_sha256(plan)
    return plan


def verify_acquisition_lock(
    document: Mapping[str, Any],
    lock: Any,
    cache_root: Path,
) -> dict[str, Any]:
    manifest = validate_manifest(document)
    lock = _require_exact_keys(
        lock, {"schema_version", "manifest_sha256", "sources"},
        label="acquisition lock",
    )
    if lock["schema_version"] != LOCK_SCHEMA_VERSION:
        raise ManifestError("unsupported acquisition lock schema")
    manifest_sha = canonical_sha256(manifest)
    if lock["manifest_sha256"] != manifest_sha:
        raise ManifestError("acquisition lock belongs to a different manifest")
    if not isinstance(lock["sources"], dict):
        raise ManifestError("acquisition lock sources must be an object")
    remote = {
        source["id"]: source
        for source in manifest["sources"]
        if source["transport"] != "existing_registry"
    }
    if set(lock["sources"]) != set(remote):
        raise ManifestError("acquisition lock source set is incomplete or unexpected")

    cache_root = Path(cache_root).resolve()
    checked_files = 0
    checked_bytes = 0
    for source_id, source in sorted(remote.items()):
        record = _require_exact_keys(
            lock["sources"][source_id], {"request_sha256", "files"},
            label=f"lock source {source_id}",
        )
        if record["request_sha256"] != source_request_sha256(source):
            raise ManifestError(f"{source_id}: request SHA does not match manifest")
        if not isinstance(record["files"], dict):
            raise ManifestError(f"{source_id}: lock files must be an object")
        expected_files = {item["role"]: item for item in source["files"]}
        if set(record["files"]) != set(expected_files):
            raise ManifestError(f"{source_id}: lock file roles differ from manifest")
        for role, locked_file in sorted(record["files"].items()):
            locked_file = _require_exact_keys(
                locked_file, {"path", "bytes", "sha256"},
                label=f"{source_id}.{role}",
            )
            digest = locked_file["sha256"]
            if not isinstance(digest, str) or SHA256_RE.fullmatch(digest) is None:
                raise ManifestError(f"{source_id}.{role}: invalid locked SHA-256")
            relative = Path(locked_file["path"])
            if relative.is_absolute() or ".." in relative.parts:
                raise ManifestError(f"{source_id}.{role}: lock path must be safe/relative")
            path = (cache_root / relative).resolve()
            try:
                path.relative_to(cache_root)
            except ValueError as exc:
                raise ManifestError(
                    f"{source_id}.{role}: lock path escapes cache"
                ) from exc
            expected_parent = Path("sha256") / digest[:2] / digest
            if relative.parent != expected_parent:
                raise ManifestError(
                    f"{source_id}.{role}: path is not content-addressed by its SHA"
                )
            if not path.is_file():
                raise ManifestError(f"{source_id}.{role}: locked file is missing")
            stat = path.stat()
            if not isinstance(locked_file["bytes"], int) or locked_file["bytes"] < 0:
                raise ManifestError(f"{source_id}.{role}: invalid byte count")
            if stat.st_size != locked_file["bytes"]:
                raise ManifestError(f"{source_id}.{role}: locked byte count changed")
            observed = sha256_file(path)
            if observed != digest:
                raise ManifestError(f"{source_id}.{role}: locked content changed")
            expected_digest = expected_files[role]["expected_sha256"]
            if expected_digest is not None and digest != expected_digest:
                raise ManifestError(f"{source_id}.{role}: expected SHA mismatch")
            checked_files += 1
            checked_bytes += stat.st_size
    return {
        "manifest_sha256": manifest_sha,
        "source_count": len(remote),
        "file_count": checked_files,
        "bytes": checked_bytes,
        "verified": True,
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    verify = subparsers.add_parser("verify", help="validate manifest schema")
    verify.add_argument("manifest", type=Path, nargs="?", default=DEFAULT_MANIFEST)

    acquire = subparsers.add_parser(
        "acquire", help="emit the acquisition plan; never downloads"
    )
    acquire.add_argument("manifest", type=Path, nargs="?", default=DEFAULT_MANIFEST)
    acquire.add_argument("--cache-root", type=Path, required=True)
    acquire.add_argument(
        "--dry-run", action="store_true", required=True,
        help="required safety acknowledgement; no filesystem/network mutation occurs",
    )

    lock = subparsers.add_parser(
        "verify-lock", help="rehash an acquisition lock and its cached bytes"
    )
    lock.add_argument("lock", type=Path)
    lock.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    lock.add_argument("--cache-root", type=Path, required=True)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        manifest = validate_manifest(load_json(args.manifest))
        if args.command == "verify":
            print(json.dumps({
                "valid": True,
                "manifest_id": manifest["manifest_id"],
                "manifest_sha256": canonical_sha256(manifest),
                "sources": len(manifest["sources"]),
                "scenarios": len(manifest["scenarios"]),
            }, sort_keys=True))
        elif args.command == "acquire":
            plan = build_acquisition_plan(manifest, args.cache_root)
            print(json.dumps(plan, indent=2, sort_keys=True))
        else:
            result = verify_acquisition_lock(
                manifest, load_json(args.lock), args.cache_root
            )
            print(json.dumps(result, sort_keys=True))
    except ManifestError as exc:
        print(f"manifest error: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
