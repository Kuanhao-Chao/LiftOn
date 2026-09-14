"""Strict, versioned benchmark campaign-profile loading.

The profile registry is deliberately declarative.  It selects existing
``build_controller`` stages and records release/truth policy without changing
the frozen four-way baseline.  A canonical digest covers the normalized
profile so a controller or meta-campaign can reject registry drift on resume.
"""
from __future__ import annotations

import hashlib
import json
import re
from pathlib import Path
from typing import Any, Mapping, Sequence


HERE = Path(__file__).resolve().parent
DEFAULT_PROFILE_REGISTRY = HERE / "campaign_profiles.json"
PROFILE_REGISTRY_SCHEMA_VERSION = 1
CAMPAIGN_SPEC_SCHEMA_VERSION = 2
LEGACY_PROFILE_ID = "canonical-v1"
DEFAULT_PROFILE_ID = LEGACY_PROFILE_ID

PAIRED_STAGES = {
    "paired-subset-canary", "paired-subset",
    "paired-full-canary", "paired-full",
    "paired-e2e-canary", "paired-e2e",
}
PROTOCOL_STAGES = {
    "protocol-thread-scaling": "thread_scaling",
    "protocol-io-modes": "io_modes",
}
SUPPORTED_STAGES = {"gates", *PAIRED_STAGES, *PROTOCOL_STAGES}
E2E_MODES = {
    "safe", "stream", "inmemory", "native",
    "stream-inmemory", "stream-native", "inmemory-native", "fast",
}
BASELINE_POLICIES = {
    "not_applicable",
    "paired_required",
    "diagnostic",
}
TRUTH_POLICIES = {
    "not_applicable",
    "none",
    "where_available",
    "target_truth_required",
    "synthetic_exact_required",
    "equivalence_required",
}
CASE_FIELDS = {
    "id",
    "stage",
    "ids",
    "repetitions",
    "threads",
    "candidate_mode",
    "reference_mode",
    "baseline_policy",
    "truth_policy",
}
PROFILE_FIELDS = {"id", "description", "campaigns"}
_ID_RE = re.compile(r"[a-z0-9][a-z0-9_.-]*")


def canonical_hash(value: Any) -> str:
    payload = json.dumps(
        value,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=True,
        allow_nan=False,
    ).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _positive_integer(value: Any) -> bool:
    return (
        isinstance(value, int)
        and not isinstance(value, bool)
        and value > 0
    )


def _normalize_case(raw: Any, *, profile_id: str) -> dict[str, Any]:
    if not isinstance(raw, Mapping):
        raise ValueError(f"profile {profile_id!r} contains a non-object campaign")
    unknown = set(raw) - CASE_FIELDS
    missing = CASE_FIELDS - set(raw)
    if unknown or missing:
        raise ValueError(
            f"profile {profile_id!r} campaign fields are invalid; "
            f"missing={sorted(missing)}, unknown={sorted(unknown)}"
        )
    case_id = raw["id"]
    if not isinstance(case_id, str) or _ID_RE.fullmatch(case_id) is None:
        raise ValueError(f"profile {profile_id!r} has invalid campaign id {case_id!r}")
    stage = raw["stage"]
    if stage not in SUPPORTED_STAGES:
        raise ValueError(
            f"profile {profile_id!r}/{case_id} has unsupported stage {stage!r}"
        )
    ids = raw["ids"]
    if (
        not isinstance(ids, list)
        or not ids
        or not all(
            isinstance(value, str) and _ID_RE.fullmatch(value) is not None
            for value in ids
        )
        or len(ids) != len(set(ids))
    ):
        raise ValueError(
            f"profile {profile_id!r}/{case_id} ids must be a non-empty "
            "ordered list of unique identifiers"
        )
    repetitions = raw["repetitions"]
    if not _positive_integer(repetitions) or repetitions > 10:
        raise ValueError(
            f"profile {profile_id!r}/{case_id} repetitions must be between 1 and 10"
        )
    threads = raw["threads"]
    if not _positive_integer(threads) or threads > 32:
        raise ValueError(
            f"profile {profile_id!r}/{case_id} threads must be between 1 and 32"
        )
    baseline_policy = raw["baseline_policy"]
    truth_policy = raw["truth_policy"]
    if baseline_policy not in BASELINE_POLICIES:
        raise ValueError(
            f"profile {profile_id!r}/{case_id} has invalid baseline_policy"
        )
    if truth_policy not in TRUTH_POLICIES:
        raise ValueError(
            f"profile {profile_id!r}/{case_id} has invalid truth_policy"
        )
    candidate_mode = raw["candidate_mode"]
    reference_mode = raw["reference_mode"]
    paired = stage in PAIRED_STAGES
    if paired:
        if candidate_mode not in E2E_MODES or reference_mode not in E2E_MODES:
            raise ValueError(
                f"profile {profile_id!r}/{case_id} paired modes are invalid"
            )
        if not stage.endswith("-canary") and repetitions % 2:
            raise ValueError(
                f"profile {profile_id!r}/{case_id} requires even AB/BA repetitions"
            )
        if baseline_policy == "not_applicable":
            raise ValueError(
                f"profile {profile_id!r}/{case_id} needs a paired baseline policy"
            )
    elif stage == "gates":
        if candidate_mode is not None or reference_mode is not None:
            raise ValueError(
                f"profile {profile_id!r}/{case_id} non-paired modes must be null"
            )
        if repetitions != 1:
            raise ValueError(
                f"profile {profile_id!r}/{case_id} non-paired repetitions must be 1"
            )
        if baseline_policy != "not_applicable" or truth_policy != "not_applicable":
            raise ValueError(
                f"profile {profile_id!r}/{case_id} non-paired policies must "
                "be not_applicable"
            )
    else:
        expected = {
            "protocol-thread-scaling": {
                "ids": ["bee"],
                "repetitions": 6,
                "threads": 32,
                "candidate_mode": "fast",
            },
            "protocol-io-modes": {
                "ids": ["arabidopsis"],
                "repetitions": 4,
                "threads": 8,
                "candidate_mode": "fast",
            },
        }[stage]
        for field, value in expected.items():
            observed = (
                ids if field == "ids" else
                repetitions if field == "repetitions" else
                threads if field == "threads" else
                candidate_mode
            )
            if observed != value:
                raise ValueError(
                    f"profile {profile_id!r}/{case_id} {field} must match "
                    "the frozen protocol schedule"
                )
        if (
            reference_mode is not None
            or baseline_policy != "diagnostic"
            or truth_policy != "equivalence_required"
        ):
            raise ValueError(
                f"profile {profile_id!r}/{case_id} protocol policy is invalid"
            )
    return {
        "id": case_id,
        "stage": stage,
        "ids": list(ids),
        "repetitions": repetitions,
        "threads": threads,
        "candidate_mode": candidate_mode,
        "reference_mode": reference_mode,
        "baseline_policy": baseline_policy,
        "truth_policy": truth_policy,
    }


def load_registry(
    path: Path = DEFAULT_PROFILE_REGISTRY,
) -> dict[str, Any]:
    """Load and strictly normalize the complete profile registry."""

    path = Path(path).resolve()
    try:
        raw = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, TypeError, ValueError) as exc:
        raise ValueError(f"campaign profile registry is unreadable: {path}: {exc}") from exc
    if not isinstance(raw, Mapping) or set(raw) != {"schema_version", "profiles"}:
        raise ValueError(
            "campaign profile registry must contain exactly schema_version and profiles"
        )
    if raw["schema_version"] != PROFILE_REGISTRY_SCHEMA_VERSION:
        raise ValueError(
            "campaign profile registry schema_version must be "
            f"{PROFILE_REGISTRY_SCHEMA_VERSION}"
        )
    profiles = raw["profiles"]
    if not isinstance(profiles, list) or not profiles:
        raise ValueError("campaign profile registry profiles must be non-empty")
    normalized_profiles = []
    seen_profiles: set[str] = set()
    for raw_profile in profiles:
        if not isinstance(raw_profile, Mapping) or set(raw_profile) != PROFILE_FIELDS:
            raise ValueError(
                "each campaign profile must contain exactly id, description, campaigns"
            )
        profile_id = raw_profile["id"]
        if (
            not isinstance(profile_id, str)
            or _ID_RE.fullmatch(profile_id) is None
            or profile_id in seen_profiles
        ):
            raise ValueError(f"campaign profile id is invalid or duplicated: {profile_id!r}")
        description = raw_profile["description"]
        if not isinstance(description, str) or not description.strip():
            raise ValueError(f"profile {profile_id!r} description must be non-empty")
        campaigns = raw_profile["campaigns"]
        if not isinstance(campaigns, list) or not campaigns:
            raise ValueError(f"profile {profile_id!r} campaigns must be non-empty")
        normalized_campaigns = [
            _normalize_case(case, profile_id=profile_id)
            for case in campaigns
        ]
        case_ids = [case["id"] for case in normalized_campaigns]
        if len(case_ids) != len(set(case_ids)):
            raise ValueError(f"profile {profile_id!r} campaign ids are duplicated")
        normalized_profiles.append({
            "id": profile_id,
            "description": description.strip(),
            "campaigns": normalized_campaigns,
        })
        seen_profiles.add(profile_id)
    normalized = {
        "schema_version": PROFILE_REGISTRY_SCHEMA_VERSION,
        "profiles": normalized_profiles,
        "registry": {
            "path": str(path),
            "size": path.stat().st_size,
            "sha256": sha256_file(path),
        },
    }
    return normalized


def load_profile(
    profile_id: str = DEFAULT_PROFILE_ID,
    *,
    registry: Path = DEFAULT_PROFILE_REGISTRY,
) -> dict[str, Any]:
    document = load_registry(registry)
    for profile in document["profiles"]:
        if profile["id"] != profile_id:
            continue
        normalized = {
            "id": profile["id"],
            "description": profile["description"],
            "campaigns": profile["campaigns"],
        }
        return {
            **normalized,
            "digest": canonical_hash(normalized),
            "registry": document["registry"],
        }
    available = [profile["id"] for profile in document["profiles"]]
    raise ValueError(
        f"unknown campaign profile {profile_id!r}; available={available}"
    )


def get_case(profile: Mapping[str, Any], case_id: str) -> dict[str, Any]:
    for case in profile.get("campaigns", []):
        if isinstance(case, Mapping) and case.get("id") == case_id:
            return dict(case)
    available = [
        case.get("id") for case in profile.get("campaigns", [])
        if isinstance(case, Mapping)
    ]
    raise ValueError(
        f"unknown campaign id {case_id!r} in profile "
        f"{profile.get('id')!r}; available={available}"
    )


def release_cases(profile: Mapping[str, Any]) -> list[dict[str, Any]]:
    """Return core publishable paired cases, excluding diagnostics/truth-only runs."""

    return [
        dict(case)
        for case in profile["campaigns"]
        if (
            case["stage"] in PAIRED_STAGES
            and not case["stage"].endswith("-canary")
            and case["baseline_policy"] == "paired_required"
            and case["truth_policy"] != "synthetic_exact_required"
        )
    ]


def campaign_counts(profile: Mapping[str, Any]) -> dict[str, Any]:
    cases = []
    total_cells = 0
    for case in profile["campaigns"]:
        if case["stage"] in PROTOCOL_STAGES:
            from .protocol_analysis import protocol_cases

            cells = len(protocol_cases(PROTOCOL_STAGES[case["stage"]]))
            arms = cells
        else:
            cells = len(case["ids"]) * case["repetitions"]
            arms = cells * (2 if case["stage"] in PAIRED_STAGES else 1)
        cases.append({
            "id": case["id"],
            "stage": case["stage"],
            "ids": len(case["ids"]),
            "repetitions": case["repetitions"],
            "cells": cells,
            "arms": arms,
        })
        total_cells += cells
    return {
        "campaigns": len(cases),
        "cells": total_cells,
        "arms": sum(case["arms"] for case in cases),
        "cases": cases,
    }


def campaign_spec(
    profile: Mapping[str, Any],
    *,
    repetitions: int | None = None,
    legacy_v1: bool = False,
) -> dict[str, Any]:
    """Build the release-report matrix from a normalized profile."""

    if repetitions is not None and (
        not _positive_integer(repetitions) or repetitions % 2
    ):
        raise ValueError("release repetitions must be a positive even integer")
    matrix = []
    panels: dict[str, dict[str, Any]] = {}
    for case in release_cases(profile):
        panel = case["stage"].removeprefix("paired-")
        case_repetitions = repetitions or case["repetitions"]
        row = panels.setdefault(panel, {
            "ids": [],
            "repetitions": case_repetitions,
        })
        if row["repetitions"] != case_repetitions:
            raise ValueError(
                f"profile {profile['id']!r} has inconsistent {panel} repetitions"
            )
        duplicated = set(row["ids"]) & set(case["ids"])
        if duplicated:
            raise ValueError(
                f"profile {profile['id']!r} duplicates {panel} ids: "
                f"{sorted(duplicated)}"
            )
        row["ids"].extend(case["ids"])
        matrix.append({
            **case,
            "panel": panel,
            "repetitions": case_repetitions,
        })
    expected_panels = (
        {"e2e"}
        if profile["id"] == "biology-v1"
        else {"subset", "full", "e2e"}
    )
    if set(panels) != expected_panels:
        raise ValueError(
            f"profile {profile['id']!r} release cases must cover "
            f"{sorted(expected_panels)}"
        )
    if legacy_v1:
        return {
            "schema_version": 1,
            "panels": panels,
        }
    return {
        "schema_version": CAMPAIGN_SPEC_SCHEMA_VERSION,
        "profile_id": profile["id"],
        "profile_digest": profile["digest"],
        "profile_registry_sha256": profile["registry"]["sha256"],
        "matrix": matrix,
        "panels": panels,
    }


def select_profile_case(
    profile_id: str,
    case_id: str,
    *,
    registry: Path = DEFAULT_PROFILE_REGISTRY,
) -> tuple[dict[str, Any], dict[str, Any]]:
    profile = load_profile(profile_id, registry=registry)
    return profile, get_case(profile, case_id)


def case_identity(
    profile: Mapping[str, Any],
    case: Mapping[str, Any],
) -> dict[str, Any]:
    return {
        "profile_id": profile["id"],
        "profile_digest": profile["digest"],
        "registry": dict(profile["registry"]),
        "case": dict(case),
    }


def profile_case_ids(profile: Mapping[str, Any]) -> Sequence[str]:
    return tuple(case["id"] for case in profile["campaigns"])
