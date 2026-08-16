"""Aggregate immutable paired runs into a concise LiftOn release report.

Raw controller artifacts remain below ``benchmarks/compare/_runs``.  This
module emits only the small, reviewable publication set: ``REPORT.md``,
``metrics.json``, and ``manifest.json``.
"""
from __future__ import annotations

import argparse
import csv
import datetime as dt
import hashlib
import json
import math
import os
import platform
import re
import shutil
import statistics
import tempfile
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

from .release_evaluation import (
    DEFAULT_BOOTSTRAP_REPLICATES,
    DEFAULT_BOOTSTRAP_SEED,
    SCHEMA_VERSION,
    STABLE_ID_FEATURE_TYPES,
    STABLE_ID_METHOD,
    _percentile,
    bootstrap_geomean_ratio,
    gff3_fingerprints,
    sha256_file,
    validate_e2e_biology,
)


QUALITY_METRICS = (
    "completeness_coding",
    "completeness_feature_total",
    "covpi",
    "recall_at_0.5",
    "recall_at_0.75",
    "recall_at_0.9",
    "recall_at_0.95",
    "intron_chain_exact_recall",
    "orf_valid_recall",
)
QUALITY_AGGREGATE_FLOOR = -0.001
QUALITY_CELL_FLOOR = -0.010
COMMON_PI_CELL_FLOOR = -0.005
ABSOLUTE_BASELINE_TOLERANCES = {
    "completeness_coding": 0.010,
    "completeness_feature_total": 0.010,
    "covpi": 0.010,
    "mean_pi": 0.010,
}
E2E_FEATURE_COMPLETENESS_FLOOR = 0.50
E2E_MEAN_PI_FLOOR = 0.50
TARGET_TRUTH_LOCUS_F1_FLOOR = 0.50
TARGET_TRUTH_DELTA_FLOOR = -0.01
TARGET_TRUTH_MIN_GROUPS = {
    "subset": 10,
    "full": 100,
    "e2e": 100,
}
PANEL_RATIO_LIMIT = 1.10
CELL_RATIO_LIMIT = 1.25
MEMORY_ENVELOPE_GIB = 192.0
APPROVED_RESOURCE_POLICY = {
    "threads_per_cell": 8,
    "max_active": 4,
    "max_full": 2,
    "max_worker_threads": 32,
    "load1_limit": 32.0,
    "min_available_gib": 256.0,
    "stagger_seconds": 15.0,
    "poll_seconds": 30.0,
}
TOOLING_FILE_KEYS = (
    "tooling_build_controller",
    "tooling_evaluator",
    "tooling_release_evaluation",
    "tooling_release_report",
    "tooling_run_benchmarks",
    "tooling_gff3_validator",
)
PANEL_CONCURRENCY = {
    # Mirrors the approved shared-host controller policy: subset cells may use
    # all four active slots; genome-scale full/E2E cells are capped at two.
    "subset": 4,
    "full": 2,
    "e2e": 2,
}
RELEASE_PANELS = tuple(PANEL_CONCURRENCY)
CANONICAL_PANEL_COUNTS = {
    "subset": 34,
    "full": 17,
    "e2e": 5,
}
CAMPAIGN_SCHEMA_VERSION = 1
DEFAULT_RELEASE_REPETITIONS = 4
_SHA1_RE = re.compile(r"[0-9a-f]{40}")
_SHA256_RE = re.compile(r"[0-9a-f]{64}")


def _controller_panel_concurrency(
    controller_evidence: Mapping[str, Any] | None,
) -> dict[str, int]:
    limits = dict(PANEL_CONCURRENCY)
    for record in (controller_evidence or {}).get("roots", []):
        plan = record.get("plan") or {}
        policy = plan.get("policy") or {}
        stage = plan.get("stage")
        panel = {
            "paired-subset": "subset",
            "paired-full": "full",
            "paired-e2e": "e2e",
        }.get(stage)
        if panel is None:
            continue
        field = "max_active" if panel == "subset" else "max_full"
        value = policy.get(field)
        if isinstance(value, int) and value > 0:
            limits[panel] = max(limits[panel], value)
    return limits


def _live_artifact_record(path: Path) -> dict[str, Any]:
    """Return a stable cryptographic record, rejecting a concurrent rewrite."""

    resolved = Path(path).resolve()
    before = resolved.stat()
    digest = sha256_file(resolved)
    after = resolved.stat()
    if (
        before.st_size != after.st_size
        or before.st_mtime_ns != after.st_mtime_ns
    ):
        raise ValueError(f"artifact changed while it was hashed: {resolved}")
    return {
        "path": str(resolved),
        "size": after.st_size,
        "mtime_ns": after.st_mtime_ns,
        "sha256": digest,
    }


def _finalization_host_snapshot() -> dict[str, Any]:
    cpu_model = None
    physical_ids = set()
    physical_cores = set()
    try:
        processor = {}
        for line in Path("/proc/cpuinfo").read_text().splitlines() + [""]:
            if line and ":" in line:
                key, value = line.split(":", 1)
                processor[key.strip()] = value.strip()
                continue
            if not processor:
                continue
            cpu_model = cpu_model or processor.get("model name")
            physical_id = processor.get("physical id")
            core_id = processor.get("core id")
            if physical_id is not None:
                physical_ids.add(physical_id)
            if physical_id is not None and core_id is not None:
                physical_cores.add((physical_id, core_id))
            processor = {}
    except OSError:
        pass
    memory_gib = None
    try:
        for line in Path("/proc/meminfo").read_text().splitlines():
            if line.startswith("MemTotal:"):
                memory_gib = int(line.split()[1]) / (1024 ** 2)
                break
    except (OSError, TypeError, ValueError):
        pass
    return {
        "classification": "reconstructed_at_report_finalization",
        "observed_at": (
            dt.datetime.now(dt.timezone.utc)
            .isoformat()
            .replace("+00:00", "Z")
        ),
        "hostname": platform.node(),
        "kernel": platform.release(),
        "machine": platform.machine(),
        "cpu_model": cpu_model,
        "logical_cpus": os.cpu_count(),
        "sockets": len(physical_ids) or None,
        "physical_cores": len(physical_cores) or None,
        "memory_gib": memory_gib,
    }


def _verify_artifact_record(
    record: Any,
    *,
    source: Path,
    label: str,
    expected_path: Path | None = None,
) -> dict[str, Any]:
    """Verify one controller record against the live artifact byte-for-byte."""

    if not isinstance(record, Mapping):
        raise ValueError(f"{source}: {label} evidence is malformed")
    raw_path = record.get("path")
    size = record.get("size")
    mtime_ns = record.get("mtime_ns")
    digest = record.get("sha256")
    if (
        not isinstance(raw_path, str)
        or not Path(raw_path).is_absolute()
        or not _integer(size)
        or size <= 0
        or not _integer(mtime_ns)
        or mtime_ns <= 0
        or not isinstance(digest, str)
        or _SHA256_RE.fullmatch(digest) is None
    ):
        raise ValueError(f"{source}: {label} evidence is malformed")
    path = Path(raw_path).resolve()
    if expected_path is not None and path != Path(expected_path).resolve():
        raise ValueError(
            f"{source}: {label} evidence path does not match the controller plan"
        )
    try:
        observed = _live_artifact_record(path)
    except OSError as exc:
        raise ValueError(
            f"{source}: {label} evidence is unavailable: {exc}"
        ) from exc
    if dict(record) != observed:
        raise ValueError(
            f"{source}: {label} changed after success "
            "(controller cryptographic evidence)"
        )
    return observed


def _truthy(value: Any) -> bool:
    return str(value).strip().lower() in {"1", "true", "yes", "t"}


def _binary_flag(value: Any, *, path: Path, line: int, field: str) -> bool:
    normalized = str(value).strip().lower()
    if normalized in {"1", "true", "yes", "t"}:
        return True
    if normalized in {"0", "false", "no", "f"}:
        return False
    raise ValueError(
        f"{path}:{line}: {field} must be a binary flag, got {value!r}"
    )


def _optional_unit_value(
    value: Any,
    *,
    path: Path,
    line: int,
    field: str,
) -> float | None:
    if value in (None, ""):
        return None
    try:
        parsed = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(
            f"{path}:{line}: {field} must be numeric when non-empty"
        ) from exc
    if not math.isfinite(parsed) or not 0.0 <= parsed <= 1.0:
        raise ValueError(
            f"{path}:{line}: {field} must be finite and within [0, 1]"
        )
    return parsed


def _number(value: Any) -> bool:
    return (
        isinstance(value, (int, float))
        and not isinstance(value, bool)
        and math.isfinite(float(value))
    )


def _integer(value: Any) -> bool:
    return isinstance(value, int) and not isinstance(value, bool)


def _normalized_git_sha(value: Any, *, label: str) -> str:
    if (
        not isinstance(value, str)
        or _SHA1_RE.fullmatch(value.lower()) is None
    ):
        raise ValueError(
            f"{label} SHA must be an exact 40-character Git SHA"
        )
    return value.lower()


def _positive_profile_value(
    pair: Mapping[str, Any],
    label: str,
    key: str,
    *,
    source: Path,
) -> float:
    value = (pair.get("versions", {}).get(label, {}).get("profile", {})).get(key)
    if not _number(value) or float(value) <= 0:
        raise ValueError(
            f"{source}: {label} profile {key!r} must be positive and finite"
        )
    return float(value)


def _e2e_biology_errors(version: Any) -> list[str]:
    """Reuse the producer's complete E2E biology contract."""

    if not isinstance(version, Mapping):
        return ["end-to-end version evidence is missing or not an object"]
    try:
        validate_e2e_biology(version)
    except (RuntimeError, TypeError, ValueError) as exc:
        return [str(exc)]
    return []


def _e2e_determinism_document(value: Any) -> Any:
    if isinstance(value, Mapping):
        return {
            key: (
                "<path>"
                if key in {"score_file", "transcripts_tsv"}
                else _e2e_determinism_document(item)
            )
            for key, item in value.items()
        }
    if isinstance(value, list):
        return [_e2e_determinism_document(item) for item in value]
    return value


def _load_transcript_metrics(path: Path, n_reference_coding: Any) -> dict[str, Any]:
    if not _integer(n_reference_coding) or n_reference_coding <= 0:
        raise ValueError(
            f"{path}: n_reference_coding must be a positive integer"
        )
    denominator = n_reference_coding
    recovered_pi: list[float] = []
    recovered_coding = 0
    coding_ids: set[str] = set()
    intron_exact = 0
    orf_valid = 0
    intron_scored = 0
    orf_scored = 0
    exon_sn: list[float] = []
    exon_sp: list[float] = []
    with path.open(encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {
            "ref_mrna_id",
            "recovered",
            "is_coding",
            "protein_identity",
            "dna_identity",
            "intron_chain_exact",
            "orf_valid",
            "exon_sn",
            "exon_sp",
        }
        missing_columns = required - set(reader.fieldnames or ())
        if missing_columns:
            raise ValueError(
                f"{path}: transcript TSV lacks columns "
                f"{sorted(missing_columns)}"
            )
        for line_number, row in enumerate(reader, 2):
            protein_identity = _optional_unit_value(
                row.get("protein_identity"),
                path=path,
                line=line_number,
                field="protein_identity",
            )
            _optional_unit_value(
                row.get("dna_identity"),
                path=path,
                line=line_number,
                field="dna_identity",
            )
            is_coding = _binary_flag(
                row.get("is_coding"),
                path=path,
                line=line_number,
                field="is_coding",
            )
            recovered = _binary_flag(
                row.get("recovered"),
                path=path,
                line=line_number,
                field="recovered",
            )
            if not is_coding:
                continue
            reference_id = str(row.get("ref_mrna_id") or "").strip()
            if not reference_id:
                raise ValueError(
                    f"{path}:{line_number}: coding row lacks ref_mrna_id"
                )
            if reference_id in coding_ids:
                raise ValueError(
                    f"{path}:{line_number}: duplicate coding reference row "
                    f"{reference_id!r}"
                )
            coding_ids.add(reference_id)
            if not recovered:
                continue
            recovered_coding += 1
            if protein_identity is not None:
                recovered_pi.append(protein_identity)
            if row.get("intron_chain_exact") not in (None, ""):
                intron_scored += 1
                intron_exact += int(_binary_flag(
                    row["intron_chain_exact"],
                    path=path,
                    line=line_number,
                    field="intron_chain_exact",
                ))
            if row.get("orf_valid") not in (None, ""):
                orf_scored += 1
                orf_valid += int(_binary_flag(
                    row["orf_valid"],
                    path=path,
                    line=line_number,
                    field="orf_valid",
                ))
            for key, values in (("exon_sn", exon_sn), ("exon_sp", exon_sp)):
                parsed = _optional_unit_value(
                    row.get(key),
                    path=path,
                    line=line_number,
                    field=key,
                )
                if parsed is not None:
                    values.append(parsed)
    if len(coding_ids) != denominator:
        raise ValueError(
            f"{path}: found {len(coding_ids)} unique coding reference rows, "
            f"expected {denominator}"
        )
    structural_counts = {
        "intron_chain_exact": intron_scored,
        "orf_valid": orf_scored,
        "exon_sn": len(exon_sn),
        "exon_sp": len(exon_sp),
    }
    incomplete = {
        field: count
        for field, count in structural_counts.items()
        if count != recovered_coding
    }
    if incomplete:
        raise ValueError(
            f"{path}: recovered coding rows have incomplete structural "
            f"quality values: expected {recovered_coding}, observed={incomplete}"
        )
    result: dict[str, Any] = {
        "n_reference_coding": denominator,
        "n_recovered_coding": recovered_coding,
        "n_recovered_coding_with_pi": len(recovered_pi),
        "sum_protein_identity": sum(recovered_pi),
        "mean_protein_identity_scored": (
            statistics.fmean(recovered_pi) if recovered_pi else None
        ),
        "median_protein_identity_scored": (
            statistics.median(recovered_pi) if recovered_pi else None
        ),
        "covpi": sum(recovered_pi) / denominator if denominator else None,
        "intron_chain_exact_successes": intron_exact,
        "intron_chain_exact_n_scored": intron_scored,
        "intron_chain_exact_recall": intron_exact / denominator if denominator else None,
        "orf_valid_successes": orf_valid,
        "orf_valid_n_scored": orf_scored,
        "orf_valid_recall": orf_valid / denominator if denominator else None,
        "exon_sn_mean": statistics.fmean(exon_sn) if exon_sn else None,
        "exon_sp_mean": statistics.fmean(exon_sp) if exon_sp else None,
    }
    for threshold in (0.5, 0.75, 0.9, 0.95):
        successes = sum(value >= threshold for value in recovered_pi)
        result[f"recall_at_{threshold:g}_successes"] = successes
        result[f"recall_at_{threshold:g}"] = (
            successes / denominator if denominator else None
        )
    return result


def _required_unit_metric(
    value: Any,
    *,
    path: Path,
    label: str,
    field: str,
) -> float:
    if not _number(value) or not 0.0 <= float(value) <= 1.0:
        raise ValueError(
            f"{path}: {label} summary {field!r} must be finite and within [0, 1]"
        )
    return float(value)


def _validate_feature_completeness(
    feature_types: Any,
    *,
    path: Path,
    label: str,
) -> tuple[dict[str, Any], float]:
    if not isinstance(feature_types, Mapping) or not feature_types:
        raise ValueError(
            f"{path}: {label} completeness_by_type is missing or empty"
        )
    normalized = {}
    for feature_type, raw in feature_types.items():
        if (
            not isinstance(feature_type, str)
            or not feature_type
            or not isinstance(raw, Mapping)
        ):
            raise ValueError(
                f"{path}: {label} completeness_by_type is malformed"
            )
        n_reference = raw.get("n_reference")
        n_recovered = raw.get("n_recovered")
        n_extra = raw.get("n_extra_copies")
        n_tool = raw.get("n_tool_features")
        if (
            not _integer(n_reference)
            or n_reference <= 0
            or not _integer(n_recovered)
            or not 0 <= n_recovered <= n_reference
            or not _integer(n_extra)
            or n_extra < 0
            or not _integer(n_tool)
            or n_tool < 0
        ):
            raise ValueError(
                f"{path}: {label} completeness_by_type[{feature_type!r}] "
                "has invalid counts"
            )
        fraction = _required_unit_metric(
            raw.get("pct_recovered"),
            path=path,
            label=label,
            field=f"completeness_by_type.{feature_type}.pct_recovered",
        )
        if not math.isclose(
            fraction,
            n_recovered / n_reference,
            rel_tol=0.0,
            abs_tol=5.1e-6,
        ):
            raise ValueError(
                f"{path}: {label} completeness_by_type[{feature_type!r}] "
                "fraction disagrees with its counts"
            )
        normalized[feature_type] = dict(raw)
    overall = normalized.get("_overall_")
    if not isinstance(overall, Mapping):
        raise ValueError(
            f"{path}: {label} completeness_by_type lacks _overall_"
        )
    return normalized, float(overall["pct_recovered"])


def _validate_stable_id_preservation(
    value: Any,
    *,
    path: Path,
    label: str,
) -> dict[str, dict[str, Any]]:
    """Validate identifier-continuity evidence independently of completeness."""

    if (
        not isinstance(value, Mapping)
        or value.get("method") != STABLE_ID_METHOD
        or not isinstance(value.get("by_type"), Mapping)
    ):
        raise ValueError(
            f"{path}: {label} stable_id_preservation is missing or malformed"
        )
    by_type = value["by_type"]
    if set(by_type) != set(STABLE_ID_FEATURE_TYPES):
        raise ValueError(
            f"{path}: {label} stable_id_preservation must contain exactly "
            f"{list(STABLE_ID_FEATURE_TYPES)}"
        )
    normalized = {}
    for feature_type in STABLE_ID_FEATURE_TYPES:
        row = by_type[feature_type]
        if not isinstance(row, Mapping):
            raise ValueError(
                f"{path}: {label} stable ID evidence for {feature_type} "
                "is malformed"
            )
        n_records = row.get("n_reference_records")
        n_records_with_id = row.get("n_reference_records_with_id")
        n_reference = row.get("n_reference_ids")
        n_preserved = row.get("n_preserved_ids")
        n_output_records = row.get("n_output_records")
        n_output_records_with_id = row.get("n_output_records_with_id")
        n_output = row.get("n_output_ids")
        counts_valid = (
            _integer(n_records)
            and n_records >= 0
            and _integer(n_records_with_id)
            and 0 <= n_records_with_id <= n_records
            and _integer(n_reference)
            and 0 <= n_reference <= n_records_with_id
            and (n_records_with_id == 0) == (n_reference == 0)
            and _integer(n_preserved)
            and 0 <= n_preserved <= n_reference
            and _integer(n_output_records)
            and n_output_records >= 0
            and _integer(n_output_records_with_id)
            and 0 <= n_output_records_with_id <= n_output_records
            and _integer(n_output)
            and 0 <= n_output <= n_output_records_with_id
            and (n_output_records_with_id == 0) == (n_output == 0)
            and n_preserved <= n_output
        )
        if not counts_valid:
            raise ValueError(
                f"{path}: {label} stable ID counts for {feature_type} "
                "are invalid"
            )
        applicable = row.get("applicable")
        expected_applicable = n_reference > 0
        if not isinstance(applicable, bool) or applicable != expected_applicable:
            raise ValueError(
                f"{path}: {label} stable ID applicability for {feature_type} "
                "disagrees with its denominator"
            )
        rate = row.get("preservation_rate")
        reason = row.get("reason")
        if applicable:
            parsed_rate = _required_unit_metric(
                rate,
                path=path,
                label=label,
                field=(
                    f"stable_id_preservation.{feature_type}."
                    "preservation_rate"
                ),
            )
            if reason is not None or not math.isclose(
                parsed_rate,
                n_preserved / n_reference,
                rel_tol=0.0,
                abs_tol=1e-12,
            ):
                raise ValueError(
                    f"{path}: {label} stable ID rate for {feature_type} "
                    "disagrees with its counts"
                )
        else:
            expected_reason = (
                "reference_feature_type_absent"
                if n_records == 0 else "no_declared_reference_ids"
            )
            if rate is not None or reason != expected_reason:
                raise ValueError(
                    f"{path}: {label} non-applicable stable ID evidence for "
                    f"{feature_type} is inconsistent"
                )
        normalized[feature_type] = dict(row)
    return normalized


def _version_quality(
    pair_path: Path,
    pair: Mapping[str, Any],
    label: str,
) -> dict[str, Any]:
    version = pair["versions"][label]
    summary = version.get("summary")
    if not isinstance(summary, Mapping):
        raise ValueError(f"{pair_path}: {label} evaluator summary is malformed")
    transcript_path = _transcript_metrics_path(pair_path, pair, label)
    transcript = _load_transcript_metrics(
        transcript_path,
        summary.get("n_reference_coding"),
    )
    completeness = _required_unit_metric(
        summary.get("completeness_coding"),
        path=pair_path,
        label=label,
        field="completeness_coding",
    )
    expected_completeness = (
        transcript["n_recovered_coding"]
        / transcript["n_reference_coding"]
    )
    if not math.isclose(
        completeness,
        expected_completeness,
        rel_tol=0.0,
        abs_tol=5.1e-6,
    ):
        raise ValueError(
            f"{pair_path}: {label} completeness_coding disagrees with "
            "the evaluator TSV"
        )
    feature_total = _required_unit_metric(
        summary.get("completeness_feature_total"),
        path=pair_path,
        label=label,
        field="completeness_feature_total",
    )
    feature_types, expected_feature_total = _validate_feature_completeness(
        summary.get("completeness_by_type"),
        path=pair_path,
        label=label,
    )
    stable_ids = _validate_stable_id_preservation(
        summary.get("stable_id_preservation"),
        path=pair_path,
        label=label,
    )
    if not math.isclose(
        feature_total,
        expected_feature_total,
        rel_tol=0.0,
        abs_tol=5.1e-6,
    ):
        raise ValueError(
            f"{pair_path}: {label} completeness_feature_total disagrees "
            "with completeness_by_type._overall_"
        )
    if summary.get("n_recovered_coding") != transcript["n_recovered_coding"]:
        raise ValueError(
            f"{pair_path}: {label} n_recovered_coding disagrees with "
            "the evaluator TSV"
        )
    identity = summary.get("protein_identity")
    if not isinstance(identity, Mapping):
        raise ValueError(
            f"{pair_path}: {label} protein_identity summary is malformed"
        )
    if identity.get("n") != transcript["n_recovered_coding_with_pi"]:
        raise ValueError(
            f"{pair_path}: {label} protein_identity.n disagrees with "
            "the evaluator TSV"
        )
    for field, observed_key in (
        ("mean", "mean_protein_identity_scored"),
        ("median", "median_protein_identity_scored"),
    ):
        observed = transcript[observed_key]
        reported = identity.get(field)
        if observed is None:
            if reported is not None:
                raise ValueError(
                    f"{pair_path}: {label} protein_identity.{field} must be null "
                    "when no coding identity was scored"
                )
        else:
            parsed = _required_unit_metric(
                reported,
                path=pair_path,
                label=label,
                field=f"protein_identity.{field}",
            )
            if not math.isclose(
                parsed,
                observed,
                rel_tol=0.0,
                abs_tol=5.1e-6,
            ):
                raise ValueError(
                    f"{pair_path}: {label} protein_identity.{field} disagrees "
                    "with the evaluator TSV"
                )
    return {
        "completeness_coding": completeness,
        "completeness_feature_total": feature_total,
        "mean_pi": identity.get("mean"),
        "median_pi": identity.get("median"),
        "feature_types": feature_types,
        "stable_id_preservation": stable_ids,
        **transcript,
    }


def _warning_inventory(pair: Mapping[str, Any], label: str) -> Counter[str]:
    return Counter(
        str(issue.get("check", "unknown"))
        for issue in pair["versions"][label]["validation"].get("issues", [])
        if str(issue.get("severity", "")).upper().endswith("WARNING")
    )


def _transcript_metrics_path(
    pair_path: Path,
    pair: Mapping[str, Any],
    label: str,
) -> Path:
    version = (pair.get("versions") or {}).get(label)
    overlay = (
        version.get("evaluation_overlay")
        if isinstance(version, Mapping) else None
    )
    if isinstance(overlay, Mapping):
        artifacts = version.get("evaluation_artifacts")
        transcript = (
            artifacts.get("transcripts_tsv")
            if isinstance(artifacts, Mapping) else None
        )
        raw_path = transcript.get("path") if isinstance(transcript, Mapping) else None
        if not isinstance(raw_path, str) or not Path(raw_path).is_absolute():
            raise ValueError(
                f"{pair_path}: {label} evaluation overlay TSV path is malformed"
            )
        return Path(raw_path).resolve()
    return (
        pair_path.parent / "evaluation" / f"{label}.transcripts.tsv"
    ).resolve()


def _verify_transcript_artifact(
    pair_path: Path,
    pair: Mapping[str, Any],
    label: str,
    controller_record: Any,
) -> Path:
    """Verify the live neutral-evaluator TSV against both evidence layers."""

    expected = (
        pair_path.parent / "evaluation" / f"{label}.transcripts.tsv"
    ).resolve()
    version = (pair.get("versions") or {}).get(label)
    artifact_group = (
        version.get("evaluation_artifacts")
        if isinstance(version, Mapping) else None
    )
    pair_record = (
        artifact_group.get("transcripts_tsv")
        if isinstance(artifact_group, Mapping) else None
    )
    if not isinstance(pair_record, Mapping):
        raise ValueError(
            f"{pair_path}: {label} pair result lacks evaluator TSV evidence"
        )
    summary = version.get("summary")
    summary_path = (
        summary.get("transcripts_tsv")
        if isinstance(summary, Mapping) else None
    )
    for source, value in (
        ("pair result", pair_record.get("path")),
        ("evaluator summary", summary_path),
    ):
        if (
            not isinstance(value, str)
            or not Path(value).is_absolute()
            or Path(value).resolve() != expected
        ):
            raise ValueError(
                f"{pair_path}: {label} {source} evaluator TSV path does not "
                "match the canonical cell path"
            )
    pair_size = pair_record.get("size")
    pair_sha = pair_record.get("sha256")
    if (
        not _integer(pair_size)
        or pair_size <= 0
        or not isinstance(pair_sha, str)
        or _SHA256_RE.fullmatch(pair_sha) is None
    ):
        raise ValueError(
            f"{pair_path}: {label} pair-result evaluator TSV evidence is malformed"
        )
    if not isinstance(controller_record, Mapping):
        raise ValueError(
            f"{pair_path}: {label} controller success lacks evaluator TSV evidence"
        )
    controller_path = controller_record.get("path")
    if (
        not isinstance(controller_path, str)
        or not Path(controller_path).is_absolute()
        or Path(controller_path).resolve() != expected
    ):
        raise ValueError(
            f"{pair_path}: {label} controller evaluator TSV path does not "
            "match the canonical cell path"
        )
    if (
        controller_record.get("size") != pair_size
        or controller_record.get("sha256") != pair_sha
        or not _integer(controller_record.get("mtime_ns"))
    ):
        raise ValueError(
            f"{pair_path}: {label} controller and pair-result evaluator TSV "
            "evidence disagree"
        )
    try:
        stat = expected.stat()
        observed_sha = sha256_file(expected)
    except OSError as exc:
        raise ValueError(
            f"{pair_path}: {label} evaluator TSV is unavailable: {exc}"
        ) from exc
    if (
        not expected.is_file()
        or stat.st_size != pair_size
        or stat.st_mtime_ns != controller_record["mtime_ns"]
        or observed_sha != pair_sha
    ):
        raise ValueError(
            f"{pair_path}: live {label} evaluator TSV no longer matches "
            "cryptographic controller evidence"
        )
    return expected


def _target_truth_summary(
    pair_path: Path,
    pair: Mapping[str, Any],
    label: str,
    controller_record: Any = None,
) -> dict[str, Any] | None:
    """Verify and normalize optional independent target-truth evidence."""

    version = (pair.get("versions") or {}).get(label)
    summary = (
        version.get("summary")
        if isinstance(version, Mapping) else None
    )
    metrics = (
        summary.get("target_truth")
        if isinstance(summary, Mapping) else None
    )
    artifacts = (
        version.get("evaluation_artifacts")
        if isinstance(version, Mapping) else None
    )
    record = (
        artifacts.get("target_truth")
        if isinstance(artifacts, Mapping) else None
    )
    if metrics is None and record is None and controller_record is None:
        return None
    if not isinstance(metrics, Mapping) or not isinstance(record, Mapping):
        raise ValueError(
            f"{pair_path}: {label} target-truth summary/artifact is incomplete"
        )
    expected = (
        pair_path.parent / "evaluation" / f"{label}.target_truth.json"
    ).resolve()
    raw_path = record.get("path")
    size = record.get("size")
    digest = record.get("sha256")
    if (
        not isinstance(raw_path, str)
        or not Path(raw_path).is_absolute()
        or Path(raw_path).resolve() != expected
        or not _integer(size)
        or size <= 0
        or not isinstance(digest, str)
        or _SHA256_RE.fullmatch(digest) is None
    ):
        raise ValueError(
            f"{pair_path}: {label} target-truth artifact evidence is malformed"
        )
    if controller_record is not None:
        if (
            not isinstance(controller_record, Mapping)
            or controller_record.get("path") != str(expected)
            or controller_record.get("size") != size
            or controller_record.get("sha256") != digest
            or not _integer(controller_record.get("mtime_ns"))
        ):
            raise ValueError(
                f"{pair_path}: {label} controller and pair target-truth "
                "evidence disagree"
            )
    try:
        stat = expected.stat()
        observed_digest = sha256_file(expected)
        observed = json.loads(expected.read_text())
    except (OSError, TypeError, ValueError) as exc:
        raise ValueError(
            f"{pair_path}: {label} target-truth artifact is unavailable: {exc}"
        ) from exc
    if (
        not expected.is_file()
        or stat.st_size != size
        or observed_digest != digest
        or (
            controller_record is not None
            and stat.st_mtime_ns != controller_record["mtime_ns"]
        )
        or observed != metrics
    ):
        raise ValueError(
            f"{pair_path}: live {label} target-truth artifact no longer "
            "matches sealed evidence"
        )
    if metrics.get("schema_version") != 1:
        raise ValueError(
            f"{pair_path}: {label} target-truth schema is unsupported"
        )
    if metrics.get("method") != "ortholog-scoped-target-coordinate-v1":
        raise ValueError(
            f"{pair_path}: {label} target-truth method is unsupported"
        )
    parameters = metrics.get("parameters")
    if not isinstance(parameters, Mapping):
        raise ValueError(
            f"{pair_path}: {label} target-truth parameters are malformed"
        )
    id_policy = parameters.get("id_policy")
    mapping_required = parameters.get("mapping_required")
    mapping_satisfied = parameters.get("mapping_requirement_satisfied")
    if (
        id_policy not in {"ortholog-map", "exact-id"}
        or mapping_required is not (id_policy == "ortholog-map")
        or mapping_satisfied is not True
    ):
        raise ValueError(
            f"{pair_path}: {label} target-truth mapping policy is unsafe"
        )
    scope = metrics.get("scope")
    if not isinstance(scope, Mapping):
        raise ValueError(f"{pair_path}: {label} target-truth scope is malformed")
    for field in ("gene_groups", "transcript_groups", "mapping_entries"):
        value = scope.get(field)
        if not _integer(value) or value < 0:
            raise ValueError(
                f"{pair_path}: {label} target-truth scope {field} is invalid"
            )
    if mapping_required and scope["mapping_entries"] <= 0:
        raise ValueError(
            f"{pair_path}: {label} target-truth mapping scope is empty"
        )
    for level, group_field in (
        ("gene", "gene_groups"),
        ("transcript", "transcript_groups"),
    ):
        level_scope = scope.get(level)
        required_scope_fields = {
            "groups",
            "predicted_scored",
            "expected_scored",
            "prediction_models_total",
            "truth_models_total",
            "prediction_models_ignored",
            "truth_models_ignored",
        }
        if (
            not isinstance(level_scope, Mapping)
            or set(level_scope) != required_scope_fields
            or any(
                not _integer(level_scope[field]) or level_scope[field] < 0
                for field in required_scope_fields
            )
            or level_scope["groups"] != scope[group_field]
            or (
                level_scope["predicted_scored"]
                + level_scope["prediction_models_ignored"]
                != level_scope["prediction_models_total"]
            )
            or (
                level_scope["expected_scored"]
                + level_scope["truth_models_ignored"]
                != level_scope["truth_models_total"]
            )
        ):
            raise ValueError(
                f"{pair_path}: {label} target-truth {level} scope is malformed"
            )
    for group, names in (
        ("gene", ("locus", "strand", "copy")),
        ("transcript", ("locus", "strand", "copy")),
        ("structure", ("intron_chain", "intron", "exon", "CDS")),
    ):
        rows = metrics.get(group)
        if not isinstance(rows, Mapping):
            raise ValueError(
                f"{pair_path}: {label} target-truth {group} is malformed"
            )
        for name in names:
            row = rows.get(name)
            if not isinstance(row, Mapping):
                raise ValueError(
                    f"{pair_path}: {label} target-truth {group}.{name} "
                    "is malformed"
                )
            for metric_name in ("precision", "recall", "f1"):
                value = row.get(metric_name)
                if value is not None and (
                    not _number(value) or not 0.0 <= float(value) <= 1.0
                ):
                    raise ValueError(
                        f"{pair_path}: {label} target-truth "
                        f"{group}.{name}.{metric_name} is invalid"
                    )
    return dict(metrics)


def _paired_common_pi(
    pair_path: Path,
    pair: Mapping[str, Any],
) -> dict[str, Any] | None:
    paths = {
        label: _transcript_metrics_path(pair_path, pair, label)
        for label in ("candidate", "reference")
    }
    recovered = {}
    for label, path in paths.items():
        values = {}
        with path.open(encoding="utf-8") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                if not _truthy(row.get("is_coding")) or not _truthy(row.get("recovered")):
                    continue
                try:
                    value = float(row["protein_identity"])
                except (KeyError, TypeError, ValueError):
                    continue
                key = str(row["ref_mrna_id"])
                values[key] = max(values.get(key, -math.inf), value)
        recovered[label] = values
    keys = sorted(set(recovered["candidate"]) & set(recovered["reference"]))
    if not keys:
        return None
    deltas = [
        recovered["candidate"][key] - recovered["reference"][key]
        for key in keys
    ]
    return {
        "n_common": len(keys),
        "mean_delta": statistics.fmean(deltas),
        "n_improved": sum(value > 1e-9 for value in deltas),
        "n_regressed": sum(value < -1e-9 for value in deltas),
        "n_tied": sum(abs(value) <= 1e-9 for value in deltas),
    }


def _controller_pair_paths(
    root: Path,
    *,
    verify_content: bool = False,
    allow_incomplete: bool = False,
) -> list[Path] | None:
    """Select current successful cell artifacts from a controller run."""

    plan_path = root / "plan.json"
    if not plan_path.is_file():
        return None
    try:
        plan = json.loads(plan_path.read_text())
    except (OSError, TypeError, ValueError) as exc:
        raise ValueError(f"{plan_path}: controller plan is unreadable: {exc}") from exc
    from . import build_controller

    if plan.get("schema_version") != build_controller.CONTROLLER_SCHEMA_VERSION:
        raise ValueError(f"{plan_path}: unsupported controller plan schema")
    if verify_content:
        try:
            build_controller.validate_plan_integrity(plan)
        except (KeyError, TypeError, ValueError) as exc:
            raise ValueError(
                f"{plan_path}: controller plan integrity validation failed: {exc}"
            ) from exc
    cells = [
        cell for cell in plan.get("cells", [])
        if isinstance(cell, Mapping) and cell.get("kind") == "paired_release"
    ]
    if not cells:
        raise ValueError(f"{plan_path}: controller plan has no paired release cells")
    selected = []
    for cell in cells:
        cell_dir = Path(str(cell.get("cell_dir", ""))).resolve()
        expected_root = (root / "cells").resolve()
        if not cell_dir.is_relative_to(expected_root):
            raise ValueError(
                f"{plan_path}: paired cell escapes the controller run: {cell_dir}"
            )
        status_path = cell_dir / "status.json"
        success_path = cell_dir / ".success"
        try:
            status = json.loads(status_path.read_text())
        except (OSError, TypeError, ValueError) as exc:
            raise ValueError(
                f"{cell_dir}: controller status is unreadable: {exc}"
            ) from exc
        if status.get("state") != "success":
            if allow_incomplete:
                continue
            raise ValueError(
                f"{cell_dir}: controller cell state is {status.get('state')!r}, "
                "not 'success'"
            )
        try:
            success = json.loads(success_path.read_text())
        except (OSError, TypeError, ValueError) as exc:
            raise ValueError(
                f"{cell_dir}: successful controller evidence is unreadable: {exc}"
            ) from exc
        if (
            status.get("schema_version")
            != build_controller.CONTROLLER_SCHEMA_VERSION
            or success.get("schema_version")
            != build_controller.CONTROLLER_SCHEMA_VERSION
        ):
            raise ValueError(
                f"{cell_dir}: unsupported controller success/status schema"
            )
        fingerprint = cell.get("fingerprint")
        if (
            not isinstance(fingerprint, str)
            or status.get("fingerprint") != fingerprint
            or success.get("fingerprint") != fingerprint
        ):
            raise ValueError(
                f"{cell_dir}: controller success fingerprint does not match the plan"
            )
        recorded_artifacts = (
            success.get("validation", {}).get("artifacts", {})
            if isinstance(success.get("validation"), Mapping)
            else {}
        )
        if not isinstance(recorded_artifacts, Mapping) or not recorded_artifacts:
            raise ValueError(
                f"{cell_dir}: controller success lacks validated artifacts"
            )
        for name, metadata in recorded_artifacts.items():
            verified = _verify_artifact_record(
                metadata,
                source=success_path,
                label=f"validated artifact {name!r}",
            )
            if not Path(verified["path"]).is_relative_to(cell_dir):
                raise ValueError(
                    f"{cell_dir}: validated artifact {name!r} escapes its cell"
                )
        result_path = Path(
            str((cell.get("artifacts") or {}).get("result_json", ""))
        ).resolve()
        if not result_path.is_relative_to(cell_dir):
            raise ValueError(
                f"{cell_dir}: paired result escapes its controller cell"
            )
        if not result_path.is_file():
            raise ValueError(f"{cell_dir}: paired result is missing: {result_path}")
        if verify_content:
            required_names = (
                "result_json",
                "candidate_gff",
                "reference_gff",
                "candidate_manifest",
                "reference_manifest",
                "candidate_gff_validation",
                "reference_gff_validation",
            )
            configured = cell.get("artifacts")
            if not isinstance(configured, Mapping):
                raise ValueError(f"{cell_dir}: controller artifacts are malformed")
            expected_artifacts = {
                name: Path(str(configured.get(name, ""))).resolve()
                for name in required_names[:5]
            }
            expected_artifacts.update({
                "candidate_gff_validation": (
                    cell_dir / "candidate_gff_validation.json"
                ).resolve(),
                "reference_gff_validation": (
                    cell_dir / "reference_gff_validation.json"
                ).resolve(),
            })
            for name in required_names:
                metadata = recorded_artifacts.get(name)
                if not isinstance(metadata, Mapping):
                    raise ValueError(
                        f"{cell_dir}: controller success lacks {name!r} evidence"
                    )
                _verify_artifact_record(
                    metadata,
                    source=success_path,
                    label=f"controller success {name!r}",
                    expected_path=expected_artifacts[name],
                )

            try:
                pair = json.loads(result_path.read_text())
            except (OSError, TypeError, ValueError) as exc:
                raise ValueError(
                    f"{result_path}: paired result is unreadable: {exc}"
                ) from exc
            if pair.get("schema_version") != SCHEMA_VERSION:
                raise ValueError(f"{result_path}: unsupported paired schema")
            expected_identity = {
                "panel": cell.get("panel"),
                "benchmark": cell.get("benchmark"),
                "repetition": cell.get("repetition"),
                "order": cell.get("expected_order"),
                "threads": cell.get("threads"),
            }
            if any(
                pair.get(name) != expected
                for name, expected in expected_identity.items()
            ):
                raise ValueError(
                    f"{result_path}: pair identity/protocol does not match "
                    "the controller cell"
                )
            pair_inputs = pair.get("inputs")
            controller_inputs = cell.get("input_fingerprints")
            if (
                not isinstance(pair_inputs, Mapping)
                or not isinstance(controller_inputs, Mapping)
                or set(pair_inputs) != set(controller_inputs)
                or any(
                    not isinstance(pair_inputs.get(name), Mapping)
                    or not isinstance(controller_inputs.get(name), Mapping)
                    or Path(str(pair_inputs[name].get("path", ""))).resolve()
                    != Path(str(controller_inputs[name].get("path", ""))).resolve()
                    or pair_inputs[name].get("size")
                    != controller_inputs[name].get("size")
                    or pair_inputs[name].get("sha256")
                    != controller_inputs[name].get("sha256")
                    for name in pair_inputs
                )
            ):
                raise ValueError(
                    f"{result_path}: pair inputs do not match controller "
                    "fingerprint evidence"
                )
            paired_configuration = plan.get("paired")
            expected_modes = {
                label: (
                    (paired_configuration.get(label) or {}).get("e2e_mode")
                    if cell.get("panel") == "e2e"
                    else cell.get("panel")
                )
                for label in ("candidate", "reference")
            }
            if pair.get("modes") != expected_modes:
                raise ValueError(
                    f"{result_path}: pair modes do not match the controller plan"
                )
            paired_provenance = (plan.get("provenance") or {}).get("paired")
            frozen_sources = (
                paired_provenance.get("sources")
                if isinstance(paired_provenance, Mapping) else None
            )
            if pair.get("provenance") != frozen_sources:
                raise ValueError(
                    f"{result_path}: pair source provenance does not match "
                    "the frozen controller plan"
                )
            if pair.get("registries") != paired_configuration.get("registries"):
                raise ValueError(
                    f"{result_path}: pair registries do not match the frozen "
                    "controller plan"
                )
            success_validation = success.get("validation")
            success_fingerprints = (
                success_validation.get("gff_fingerprints", {})
                if isinstance(success_validation, Mapping)
                else {}
            )
            success_evaluation = (
                success_validation.get("evaluation_artifacts", {})
                if isinstance(success_validation, Mapping)
                else {}
            )
            for label in ("candidate", "reference"):
                gff_path = Path(str(configured[f"{label}_gff"])).resolve()
                observed = gff3_fingerprints(gff_path)
                if success_fingerprints.get(label) != observed:
                    raise ValueError(
                        f"{cell_dir}: live {label} GFF3 no longer matches "
                        "controller success evidence"
                    )
                version = (pair.get("versions") or {}).get(label)
                if not isinstance(version, Mapping):
                    raise ValueError(
                        f"{result_path}: {label} version evidence is malformed"
                    )
                if Path(str(version.get("output_gff", ""))).resolve() != gff_path:
                    raise ValueError(
                        f"{result_path}: {label} output path does not match "
                        "the controller plan"
                    )
                manifest_path = Path(
                    str(configured[f"{label}_manifest"])
                ).resolve()
                if (
                    Path(str(version.get("release_manifest", ""))).resolve()
                    != manifest_path
                ):
                    raise ValueError(
                        f"{result_path}: {label} release-manifest path does "
                        "not match the controller plan"
                    )
                if version.get("fingerprints") != observed:
                    raise ValueError(
                        f"{result_path}: {label} GFF3 fingerprints disagree "
                        "with the live artifact"
                    )
                label_evaluation = (
                    success_evaluation.get(label)
                    if isinstance(success_evaluation, Mapping) else None
                )
                transcript_record = (
                    label_evaluation.get("transcripts_tsv")
                    if isinstance(label_evaluation, Mapping) else None
                )
                transcript_path = _verify_transcript_artifact(
                    result_path,
                    pair,
                    label,
                    transcript_record,
                )
                truth_record = (
                    label_evaluation.get("target_truth")
                    if isinstance(label_evaluation, Mapping) else None
                )
                truth_metrics = (
                    (version.get("summary") or {}).get("target_truth")
                    if isinstance(version.get("summary"), Mapping) else None
                )
                if truth_record is not None or truth_metrics is not None:
                    _target_truth_summary(
                        result_path,
                        pair,
                        label,
                        truth_record,
                    )
                if (
                    transcript_path.stat().st_mtime_ns
                    > result_path.stat().st_mtime_ns
                ):
                    raise ValueError(
                        f"{cell_dir}: evaluator TSV changed after pair_result.json "
                        f"was published: {transcript_path}"
                    )
        selected.append(result_path)
    return selected


def load_pairs(
    roots: Iterable[Path], *, allow_incomplete: bool = False,
) -> list[tuple[Path, dict[str, Any]]]:
    pairs = []
    seen = set()
    for root in roots:
        root = Path(root).resolve()
        controller_paths = _controller_pair_paths(
            root, allow_incomplete=allow_incomplete,
        )
        paths = (
            controller_paths
            if controller_paths is not None
            else sorted(root.rglob("pair_result.json"))
        )
        for path in paths:
            resolved = path.resolve()
            if resolved in seen:
                continue
            raw = json.loads(path.read_text())
            if raw.get("schema_version") != SCHEMA_VERSION:
                raise ValueError(f"unsupported paired schema in {path}")
            pairs.append((path, raw))
            seen.add(resolved)
    return pairs


def _canonical_release_ids(
    profile_id: str = "canonical-v1",
    *,
    profile_registry: Path | None = None,
) -> dict[str, tuple[str, ...]]:
    """Load explicit release panels without honoring partial ``--ids``."""

    from . import build_controller, campaign_profiles

    profile = campaign_profiles.load_profile(
        profile_id,
        registry=Path(
            profile_registry or campaign_profiles.DEFAULT_PROFILE_REGISTRY
        ),
    )
    specification = campaign_profiles.campaign_spec(
        profile,
        legacy_v1=profile_id == campaign_profiles.LEGACY_PROFILE_ID,
    )
    selected = {
        panel: tuple(specification["panels"][panel]["ids"])
        for panel in RELEASE_PANELS
    }
    if profile_id == campaign_profiles.LEGACY_PROFILE_ID:
        malformed = {
            panel: len(ids)
            for panel, ids in selected.items()
            if (
                len(ids) != CANONICAL_PANEL_COUNTS[panel]
                or len(ids) != len(set(ids))
            )
        }
        baseline_selected = {
            panel: tuple(build_controller.select_ids(
                f"paired-{panel}",
                baseline=build_controller.DEFAULT_BASELINE,
                dataset_registry=build_controller.DEFAULT_DATASET_REGISTRY,
            ))
            for panel in RELEASE_PANELS
        }
        if malformed or selected != baseline_selected:
            raise RuntimeError(
                "frozen canonical-v1 profile disagrees with repository "
                f"baseline panels: malformed={malformed}"
            )
    return selected


def _canonical_quality_baselines() -> tuple[dict[str, Any], dict[str, Any]]:
    """Load frozen per-cell stable quality, with explicit legacy fallback."""

    from . import build_controller

    path = Path(build_controller.DEFAULT_BASELINE).resolve()
    document = json.loads(path.read_text())
    if not isinstance(document, Mapping):
        raise ValueError(f"{path}: canonical quality baseline is not an object")
    baselines = {}
    for panel in ("subset", "full"):
        for benchmark in _canonical_release_ids()[panel]:
            key = f"{panel}:{benchmark}"
            row = document.get(key)
            if not isinstance(row, Mapping):
                raise ValueError(f"{path}: canonical baseline lacks {key}")
            metrics = {}
            for metric, source_path in (
                ("completeness_coding", ("completeness_coding",)),
                (
                    "completeness_feature_total",
                    ("completeness_feature_total",),
                ),
                ("mean_pi", ("mean_pi",)),
                ("covpi", ("joint", "covpi")),
            ):
                current: Any = row
                for field in source_path:
                    current = (
                        current.get(field)
                        if isinstance(current, Mapping) else None
                    )
                stable = (
                    current.get("lifton_stable")
                    if isinstance(current, Mapping) else None
                )
                tool = "lifton_stable"
                if not _number(stable):
                    stable = (
                        current.get("lifton_devel")
                        if isinstance(current, Mapping) else None
                    )
                    tool = "lifton_devel_frozen_fallback"
                if not _number(stable) or not 0.0 <= float(stable) <= 1.0:
                    raise ValueError(
                        f"{path}: {key} has no finite canonical {metric}"
                    )
                tolerance = ABSOLUTE_BASELINE_TOLERANCES[metric]
                metrics[metric] = {
                    "baseline": float(stable),
                    "tool": tool,
                    "tolerance": tolerance,
                    "floor": max(0.0, float(stable) - tolerance),
                }
            baselines[key] = metrics
    return baselines, {
        "path": str(path),
        "size": path.stat().st_size,
        "sha256": sha256_file(path),
    }


def canonical_campaign_spec(
    repetitions: int = DEFAULT_RELEASE_REPETITIONS,
    *,
    profile_id: str = "canonical-v1",
    profile_registry: Path | None = None,
) -> dict[str, Any]:
    """Return a reviewable template for the complete canonical campaign."""

    from . import campaign_profiles

    profile = campaign_profiles.load_profile(
        profile_id,
        registry=Path(
            profile_registry or campaign_profiles.DEFAULT_PROFILE_REGISTRY
        ),
    )
    return campaign_profiles.campaign_spec(
        profile,
        repetitions=repetitions,
        legacy_v1=profile_id == campaign_profiles.LEGACY_PROFILE_ID,
    )


def _normalize_profile_campaign_spec(document: Mapping[str, Any]) -> dict[str, Any]:
    from . import campaign_profiles

    required = {
        "schema_version",
        "profile_id",
        "profile_digest",
        "profile_registry_sha256",
        "matrix",
        "panels",
    }
    if set(document) != required:
        raise ValueError(
            "profile campaign specification fields are invalid; "
            f"missing={sorted(required - set(document))}, "
            f"unknown={sorted(set(document) - required)}"
        )
    for name in ("profile_digest", "profile_registry_sha256"):
        value = document[name]
        if not isinstance(value, str) or _SHA256_RE.fullmatch(value) is None:
            raise ValueError(f"profile campaign {name} must be a SHA-256 digest")
    if not isinstance(document["profile_id"], str) or not document["profile_id"]:
        raise ValueError("profile campaign profile_id must be non-empty")
    matrix = document["matrix"]
    if not isinstance(matrix, list) or not matrix:
        raise ValueError("profile campaign matrix must be non-empty")
    normalized_matrix = []
    matrix_panels: dict[str, dict[str, Any]] = {}
    for raw in matrix:
        if not isinstance(raw, Mapping):
            raise ValueError("profile campaign matrix rows must be objects")
        expected_fields = campaign_profiles.CASE_FIELDS | {"panel"}
        if set(raw) != expected_fields:
            raise ValueError("profile campaign matrix row fields are invalid")
        case = campaign_profiles._normalize_case(
            {key: raw[key] for key in campaign_profiles.CASE_FIELDS},
            profile_id=document["profile_id"],
        )
        panel = raw["panel"]
        if (
            panel not in RELEASE_PANELS
            or case["stage"] != f"paired-{panel}"
            or case["baseline_policy"] != "paired_required"
        ):
            raise ValueError(
                f"profile campaign matrix case {case['id']!r} is not publishable"
            )
        row = matrix_panels.setdefault(panel, {
            "ids": [],
            "repetitions": case["repetitions"],
        })
        if row["repetitions"] != case["repetitions"]:
            raise ValueError(
                f"profile campaign panel {panel!r} mixes repetition counts"
            )
        duplicated = set(row["ids"]) & set(case["ids"])
        if duplicated:
            raise ValueError(
                f"profile campaign panel {panel!r} duplicates ids "
                f"{sorted(duplicated)}"
            )
        row["ids"].extend(case["ids"])
        normalized_matrix.append({**case, "panel": panel})
    if set(matrix_panels) != set(RELEASE_PANELS):
        raise ValueError(
            "profile campaign matrix must cover subset, full, and e2e"
        )
    panels = document["panels"]
    if not isinstance(panels, Mapping) or dict(panels) != matrix_panels:
        raise ValueError(
            "profile campaign panels do not match the ordered campaign matrix"
        )
    return {
        "schema_version": campaign_profiles.CAMPAIGN_SPEC_SCHEMA_VERSION,
        "profile_id": document["profile_id"],
        "profile_digest": document["profile_digest"],
        "profile_registry_sha256": document["profile_registry_sha256"],
        "matrix": normalized_matrix,
        "panels": matrix_panels,
    }


def _normalize_campaign_spec(document: Any) -> dict[str, Any]:
    if not isinstance(document, Mapping):
        raise ValueError("campaign specification must be a JSON object")
    if document.get("schema_version") == 2:
        return _normalize_profile_campaign_spec(document)
    if document.get("schema_version") != CAMPAIGN_SCHEMA_VERSION:
        raise ValueError(
            f"campaign specification schema_version must be "
            f"{CAMPAIGN_SCHEMA_VERSION}"
        )
    panels = document.get("panels")
    if not isinstance(panels, Mapping) or set(panels) != set(RELEASE_PANELS):
        raise ValueError(
            "campaign specification must define exactly subset, full, and e2e"
        )
    normalized: dict[str, Any] = {
        "schema_version": CAMPAIGN_SCHEMA_VERSION,
        "panels": {},
    }
    canonical_ids = _canonical_release_ids()
    for panel in RELEASE_PANELS:
        row = panels.get(panel)
        if not isinstance(row, Mapping):
            raise ValueError(f"campaign panel {panel!r} must be an object")
        ids = row.get("ids")
        if (
            not isinstance(ids, list)
            or not ids
            or not all(isinstance(value, str) and value for value in ids)
            or len(ids) != len(set(ids))
        ):
            raise ValueError(
                f"campaign panel {panel!r} ids must be a non-empty unique list"
            )
        if ids != list(canonical_ids[panel]):
            missing = sorted(set(canonical_ids[panel]) - set(ids))
            unexpected = sorted(set(ids) - set(canonical_ids[panel]))
            raise ValueError(
                f"campaign panel {panel!r} must enumerate the canonical "
                f"{len(canonical_ids[panel])} IDs in canonical order; "
                f"missing={missing}, unexpected={unexpected}"
            )
        repetitions = row.get("repetitions")
        if (
            not _integer(repetitions)
            or repetitions <= 0
            or repetitions % 2
        ):
            raise ValueError(
                f"campaign panel {panel!r} repetitions must be a positive "
                "even integer so AB/BA order is balanced"
            )
        normalized["panels"][panel] = {
            "ids": list(ids),
            "repetitions": repetitions,
        }
    return normalized


def _expected_campaign_keys(
    campaign: Mapping[str, Any],
) -> set[tuple[str, str, int]]:
    return {
        (panel, benchmark, repetition)
        for panel, row in campaign["panels"].items()
        for benchmark in row["ids"]
        for repetition in range(1, row["repetitions"] + 1)
    }


def _verify_source_file_record(
    record: Any,
    *,
    source: Path,
    label: str,
    expected_path: Path | None = None,
) -> dict[str, Any]:
    """Verify a provenance file record whose contract omits mutable mtimes."""

    if not isinstance(record, Mapping):
        raise ValueError(f"{source}: {label} provenance is malformed")
    raw_path = record.get("path")
    size = record.get("size")
    digest = record.get("sha256")
    if (
        not isinstance(raw_path, str)
        or not Path(raw_path).is_absolute()
        or not _integer(size)
        or size <= 0
        or not isinstance(digest, str)
        or _SHA256_RE.fullmatch(digest) is None
    ):
        raise ValueError(f"{source}: {label} provenance is malformed")
    path = Path(raw_path).resolve()
    if expected_path is not None and path != Path(expected_path).resolve():
        raise ValueError(
            f"{source}: {label} provenance path is not the expected source file"
        )
    try:
        observed = _live_artifact_record(path)
    except OSError as exc:
        raise ValueError(
            f"{source}: {label} provenance file is unavailable: {exc}"
        ) from exc
    normalized = {
        "path": observed["path"],
        "size": observed["size"],
        "sha256": observed["sha256"],
    }
    if dict(record) != normalized:
        raise ValueError(
            f"{source}: {label} provenance file changed after planning"
        )
    return normalized


def _current_reporter_git_state() -> dict[str, Any]:
    from . import build_controller

    return build_controller.collect_git_state(build_controller.REPO_ROOT)


def _validate_resource_policy(policy: Mapping[str, Any], plan_path: Path) -> None:
    for name, expected in APPROVED_RESOURCE_POLICY.items():
        value = policy.get(name)
        if (
            not _number(value)
            or (
                isinstance(expected, int)
                and (not _integer(value) or value != expected)
            )
            or (
                not isinstance(expected, int)
                and float(value) != expected
            )
        ):
            raise ValueError(
                f"{plan_path}: controller resource policy {name!r} is "
                f"{value!r}, expected {expected!r}"
            )


def _expected_tooling_paths() -> dict[str, Path]:
    from . import build_controller

    root = Path(build_controller.REPO_ROOT).resolve()
    return {
        "tooling_build_controller": Path(build_controller.__file__).resolve(),
        "tooling_evaluator": (
            root / "benchmarks" / "compare" / "evaluator.py"
        ).resolve(),
        "tooling_release_evaluation": (
            root / "benchmarks" / "compare" / "release_evaluation.py"
        ).resolve(),
        "tooling_release_report": Path(__file__).resolve(),
        "tooling_run_benchmarks": (
            root / "benchmarks" / "run_benchmarks.py"
        ).resolve(),
        "tooling_gff3_validator": (
            root / "lifton" / "gff3_validator.py"
        ).resolve(),
    }


def _controller_publication_evidence(
    roots: Sequence[Path],
    campaign: Mapping[str, Any],
    *,
    candidate_sha: str,
    reference_sha: str,
) -> tuple[list[Path], dict[str, Any]]:
    """Require complete controller roots for every declared release case."""

    from . import build_controller

    profile_matrix = campaign.get("matrix")
    expected_root_count = (
        len(profile_matrix)
        if isinstance(profile_matrix, list)
        else len(RELEASE_PANELS)
    )
    if len(roots) != expected_root_count:
        raise ValueError(
            f"release publication requires exactly {expected_root_count} "
            "controller roots for the declared campaign"
        )
    by_panel: dict[str, dict[str, Any]] = {}
    by_case: dict[str, dict[str, Any]] = {}
    pair_paths: list[Path] = []
    compatibility_signature: str | None = None
    compatibility_sha256: str | None = None
    reporter_git = _current_reporter_git_state()
    clean_diff = hashlib.sha256(b"").hexdigest()
    if (
        reporter_git.get("dirty_tracked_paths") != []
        or reporter_git.get("tracked_diff_sha256") != clean_diff
    ):
        raise ValueError(
            "release reporting requires the current reporter checkout to be "
            "clean and committed"
        )
    tooling_paths = _expected_tooling_paths()
    for raw_root in roots:
        root = Path(raw_root).resolve()
        plan_path = root / "plan.json"
        try:
            plan = json.loads(plan_path.read_text())
        except (OSError, TypeError, ValueError) as exc:
            raise ValueError(
                f"{plan_path}: release publication requires a readable "
                f"controller plan: {exc}"
            ) from exc
        if plan.get("schema_version") != build_controller.CONTROLLER_SCHEMA_VERSION:
            raise ValueError(f"{plan_path}: unsupported controller plan schema")
        stage = plan.get("stage")
        if isinstance(stage, str) and stage.endswith("-canary"):
            raise ValueError(
                f"{plan_path}: canary stage {stage!r} is diagnostic and cannot "
                "be published as a full release campaign"
            )
        try:
            build_controller.validate_plan_integrity(plan)
        except (KeyError, TypeError, ValueError) as exc:
            raise ValueError(
                f"{plan_path}: controller plan integrity validation failed: {exc}"
            ) from exc
        if not isinstance(stage, str) or stage not in {
            f"paired-{panel}" for panel in RELEASE_PANELS
        }:
            raise ValueError(
                f"{plan_path}: stage must be paired-subset, paired-full, "
                "or paired-e2e"
            )
        provenance = plan.get("provenance")
        policy = plan.get("policy")
        if not isinstance(provenance, Mapping) or not isinstance(policy, Mapping):
            raise ValueError(
                f"{plan_path}: controller provenance or policy is malformed"
            )
        plan_cells = plan.get("cells")
        if (
            not isinstance(plan_cells, list)
            or not plan_cells
            or not all(isinstance(cell, Mapping) for cell in plan_cells)
        ):
            raise ValueError(
                f"{plan_path}: controller cells are missing or malformed"
            )
        _validate_resource_policy(policy, plan_path)
        panel = stage.removeprefix("paired-")
        expected_case = None
        if isinstance(profile_matrix, list):
            identity = plan.get("campaign_case")
            case = identity.get("case") if isinstance(identity, Mapping) else None
            if (
                not isinstance(identity, Mapping)
                or identity.get("profile_id") != campaign.get("profile_id")
                or identity.get("profile_digest") != campaign.get(
                    "profile_digest"
                )
                or not isinstance(identity.get("registry"), Mapping)
                or identity["registry"].get("sha256")
                != campaign.get("profile_registry_sha256")
                or not isinstance(case, Mapping)
            ):
                raise ValueError(
                    f"{plan_path}: controller campaign profile identity "
                    "does not match the report"
                )
            matches = [
                row for row in profile_matrix
                if row.get("id") == case.get("id")
            ]
            if len(matches) != 1 or case != {
                key: matches[0][key]
                for key in matches[0] if key != "panel"
            }:
                raise ValueError(
                    f"{plan_path}: controller campaign case does not match "
                    "the release matrix"
                )
            expected_case = matches[0]
            if expected_case["id"] in by_case:
                raise ValueError(
                    f"duplicate controller root for campaign case "
                    f"{expected_case['id']!r}"
                )
            expected = {
                "ids": expected_case["ids"],
                "repetitions": expected_case["repetitions"],
            }
        else:
            if panel in by_panel:
                raise ValueError(
                    f"duplicate controller root for panel {panel!r}"
                )
            expected = campaign["panels"][panel]
        if plan.get("ids") != expected["ids"]:
            raise ValueError(
                f"{plan_path}: controller ids do not match the expected "
                f"{panel!r} campaign"
            )
        paired = plan.get("paired")
        if not isinstance(paired, Mapping):
            raise ValueError(f"{plan_path}: paired configuration is missing")
        if (
            paired.get("panel") != panel
            or paired.get("repetitions") != expected["repetitions"]
        ):
            raise ValueError(
                f"{plan_path}: paired panel/repetition configuration does not "
                "match the campaign specification"
            )
        inputs = plan.get("inputs")
        registries = paired.get("registries")
        if (
            not isinstance(inputs, Mapping)
            or not isinstance(registries, Mapping)
            or Path(str(inputs.get("registry", ""))).resolve()
            != Path(str(registries.get("benchmark", ""))).resolve()
            or Path(str(inputs.get("dataset_registry", ""))).resolve()
            != Path(str(registries.get("dataset", ""))).resolve()
        ):
            raise ValueError(
                f"{plan_path}: controller input registries disagree with "
                "the frozen paired configuration"
            )
        for label, expected_sha in (
            ("candidate", candidate_sha),
            ("reference", reference_sha),
        ):
            source = paired.get(label)
            if not isinstance(source, Mapping) or source.get("sha") != expected_sha:
                raise ValueError(
                    f"{plan_path}: {label} SHA does not match the report"
                )

        cells = [
            cell for cell in plan_cells
            if cell.get("kind") == "paired_release"
        ]
        if len(cells) != len(plan_cells):
            raise ValueError(
                f"{plan_path}: release panel contains non-paired cells"
            )
        expected_threads = (
            expected_case["threads"]
            if expected_case is not None
            else APPROVED_RESOURCE_POLICY["threads_per_cell"]
        )
        if any(cell.get("threads") != expected_threads for cell in cells):
            raise ValueError(
                f"{plan_path}: a paired cell violates the approved thread policy"
            )
        if expected_case is not None and (
            paired["candidate"].get("e2e_mode")
            != expected_case["candidate_mode"]
            or paired["reference"].get("e2e_mode")
            != expected_case["reference_mode"]
        ):
            raise ValueError(
                f"{plan_path}: paired modes do not match the campaign case"
            )
        observed_keys = {
            (panel, cell.get("benchmark"), cell.get("repetition"))
            for cell in cells
        }
        expected_keys = {
            (panel, benchmark, repetition)
            for benchmark in expected["ids"]
            for repetition in range(1, expected["repetitions"] + 1)
        }
        if observed_keys != expected_keys or len(cells) != len(expected_keys):
            raise ValueError(
                f"{plan_path}: controller cells do not exactly match the "
                "expected campaign matrix"
            )
        selected = _controller_pair_paths(root, verify_content=True)
        if selected is None or len(selected) != len(expected_keys):
            raise ValueError(
                f"{plan_path}: controller did not provide every expected result"
            )
        pair_paths.extend(selected)

        required_provenance = {
            name: provenance.get(name)
            for name in ("git", "files", "tools", "runtime")
        }
        if any(
            not isinstance(value, Mapping) or not value
            for value in required_provenance.values()
        ):
            raise ValueError(
                f"{plan_path}: release controller provenance lacks git, files, "
                "tools, or Python runtime evidence"
            )
        git_state = required_provenance["git"]
        if (
            git_state.get("dirty_tracked_paths") != []
            or git_state.get("tracked_diff_sha256") != clean_diff
        ):
            raise ValueError(
                f"{plan_path}: release controller was planned from a dirty "
                "tracked working tree"
            )
        if git_state != reporter_git:
            raise ValueError(
                f"{plan_path}: current reporter Git HEAD/state does not match "
                "the controller plan"
            )
        files = required_provenance["files"]
        for name in TOOLING_FILE_KEYS:
            _verify_source_file_record(
                files.get(name),
                source=plan_path,
                label=name,
                expected_path=tooling_paths[name],
            )
        for name, expected_path in (
            ("benchmark_registry", Path(registries["benchmark"])),
            ("dataset_registry", Path(registries["dataset"])),
            ("baseline", Path(inputs.get("baseline", ""))),
        ):
            _verify_source_file_record(
                files.get(name),
                source=plan_path,
                label=name,
                expected_path=expected_path,
            )
        paired_provenance = provenance.get("paired")
        paired_sources = (
            paired_provenance.get("sources")
            if isinstance(paired_provenance, Mapping) else None
        )
        if not isinstance(paired_sources, Mapping) or set(paired_sources) != {
            "candidate", "reference",
        }:
            raise ValueError(
                f"{plan_path}: paired source provenance is missing or malformed"
            )
        common = {
            **required_provenance,
            "paired_sources": paired_sources,
            "policy": dict(policy),
            "lifton_executable": paired.get("lifton_executable"),
            "registries": paired.get("registries"),
            "repetitions": paired.get("repetitions"),
            "sources": {
                label: {
                    "root": (paired.get(label) or {}).get("root"),
                    "sha": (paired.get(label) or {}).get("sha"),
                }
                for label in ("candidate", "reference")
            },
        }
        signature = _canonical_json(
            common,
            source=plan_path,
            field="cross-root controller provenance",
        )
        if compatibility_signature is None:
            compatibility_signature = signature
            compatibility_sha256 = hashlib.sha256(signature.encode()).hexdigest()
        elif signature != compatibility_signature:
            raise ValueError(
                f"{plan_path}: controller toolchain, registry, source, or "
                "thread provenance differs across release roots"
            )
        root_record = {
            "root": str(root),
            "plan": _live_artifact_record(plan_path),
            "stage": stage,
            "ids": list(plan["ids"]),
            "repetitions": paired["repetitions"],
            "resource_policy": {
                name: policy[name] for name in APPROVED_RESOURCE_POLICY
            },
        }
        if expected_case is not None:
            root_record["campaign_id"] = expected_case["id"]
            by_case[expected_case["id"]] = root_record
        else:
            by_panel[panel] = root_record
    if isinstance(profile_matrix, list):
        expected_cases = [row["id"] for row in profile_matrix]
        if set(by_case) != set(expected_cases):
            raise ValueError(
                "release controller roots do not cover every campaign case"
            )
        return pair_paths, {
            "validated": True,
            "compatibility_sha256": compatibility_sha256,
            "profile_id": campaign["profile_id"],
            "profile_digest": campaign["profile_digest"],
            "campaigns": {
                case_id: by_case[case_id]
                for case_id in expected_cases
            },
        }
    if set(by_panel) != set(RELEASE_PANELS):
        raise ValueError("release controller roots do not cover every panel")
    return pair_paths, {
        "validated": True,
        "compatibility_sha256": compatibility_sha256,
        "panels": {
            panel: by_panel[panel]
            for panel in RELEASE_PANELS
        },
    }


def _controller_execution_environment(
    provenance: Mapping[str, Any],
) -> dict[str, Any]:
    runtime = provenance.get("runtime")
    tools = provenance.get("tools")
    runtime = runtime if isinstance(runtime, Mapping) else {}
    tools = tools if isinstance(tools, Mapping) else {}
    python = runtime.get("python")
    distributions = runtime.get("distributions")
    python = python if isinstance(python, Mapping) else {}
    distributions = (
        distributions if isinstance(distributions, Mapping) else {}
    )
    package_names = (
        "duckdb", "gffutils", "mappy", "numpy", "parasail",
        "pyfaidx", "pysam",
    )
    tool_names = {
        "minimap2": "minimap2_bin",
        "miniprot": "miniprot_bin",
    }
    return {
        "python": python.get("version"),
        "tools": {
            label: (tools.get(key) or {}).get("version")
            for label, key in tool_names.items()
            if isinstance(tools.get(key), Mapping)
        },
        "packages": {
            name: distributions[name].get("version")
            for name in package_names
            if isinstance(distributions.get(name), Mapping)
        },
    }


def _diagnostic_controller_evidence(
    roots: Sequence[Path],
) -> dict[str, Any]:
    """Record partial controller plans and execution policies for diagnostics."""

    from . import build_controller

    records = []
    for raw_root in roots:
        root = Path(raw_root).resolve()
        plan_path = root / "plan.json"
        record: dict[str, Any] = {
            "root": str(root),
            "plan": None,
            "errors": [],
        }
        if not plan_path.is_file():
            records.append(record)
            continue
        try:
            plan = json.loads(plan_path.read_text())
            build_controller.validate_plan_integrity(plan)
            cells = plan.get("cells", [])
            counts = Counter()
            for cell in cells:
                status_path = Path(cell["cell_dir"]) / "status.json"
                status = json.loads(status_path.read_text())
                counts[status.get("state", "unknown")] += 1
            selected = _controller_pair_paths(root, allow_incomplete=True)
            controller_state = "unknown"
            controller_path = root / "controller.json"
            if controller_path.is_file():
                controller_state = json.loads(controller_path.read_text()).get(
                    "state", "unknown"
                )
            record["plan"] = {
                "run_id": plan.get("run_id"),
                "stage": plan.get("stage"),
                "fingerprint": plan.get("fingerprint"),
                "sha256": sha256_file(plan_path),
                "ids": list(plan.get("ids", [])),
                "repetitions": (plan.get("paired") or {}).get("repetitions"),
                "policy": dict(plan.get("policy") or {}),
                "counts": dict(sorted(counts.items())),
                "controller_state": controller_state,
                "successful_pair_results": len(selected or []),
                "execution_environment": _controller_execution_environment(
                    plan.get("provenance") or {},
                ),
            }
        except (OSError, TypeError, ValueError, KeyError) as exc:
            record["errors"].append(str(exc))
        records.append(record)
    return {
        "validated": False,
        "kind": "diagnostic_controller_roots",
        "roots": records,
    }


def _controller_artifact_evidence(
    roots: Sequence[Path],
) -> dict[str, Any]:
    """Rehash the complete controller evidence set for publication."""

    cells = []
    tooling = None
    provenance_files = None
    for raw_root in roots:
        root = Path(raw_root).resolve()
        plan_path = root / "plan.json"
        plan = json.loads(plan_path.read_text())
        files = plan["provenance"]["files"]
        current_tooling = {
            name: _verify_source_file_record(
                files[name],
                source=plan_path,
                label=name,
                expected_path=_expected_tooling_paths()[name],
            )
            for name in TOOLING_FILE_KEYS
        }
        current_provenance_files = {
            name: _verify_source_file_record(
                record,
                source=plan_path,
                label=name,
            )
            for name, record in sorted(files.items())
        }
        if tooling is None:
            tooling = current_tooling
            provenance_files = current_provenance_files
        elif tooling != current_tooling:
            raise ValueError(
                f"{plan_path}: tooling evidence differs across controller roots"
            )
        elif provenance_files != current_provenance_files:
            raise ValueError(
                f"{plan_path}: source-file evidence differs across controller roots"
            )
        for cell in plan["cells"]:
            if cell.get("kind") != "paired_release":
                continue
            cell_dir = Path(cell["cell_dir"]).resolve()
            status_path = cell_dir / "status.json"
            success_path = cell_dir / ".success"
            success = json.loads(success_path.read_text())
            artifacts = success["validation"]["artifacts"]
            evaluations = success["validation"]["evaluation_artifacts"]
            validated_artifacts = {
                name: _verify_artifact_record(
                    record,
                    source=success_path,
                    label=f"validated artifact {name!r}",
                )
                for name, record in sorted(artifacts.items())
            }
            evaluator_artifacts = {
                label: {
                    name: _verify_artifact_record(
                        record,
                        source=success_path,
                        label=f"{label} evaluator artifact {name!r}",
                    )
                    for name, record in sorted(evidence.items())
                }
                for label, evidence in sorted(evaluations.items())
                if isinstance(evidence, Mapping)
            }
            cells.append({
                "panel": cell["panel"],
                "benchmark": cell["benchmark"],
                "repetition": cell["repetition"],
                "status": _live_artifact_record(status_path),
                "success": _live_artifact_record(success_path),
                "validated_artifacts": validated_artifacts,
                "evaluator_artifacts": evaluator_artifacts,
            })
    return {
        "tooling": tooling or {},
        "provenance_files": provenance_files or {},
        "cells": cells,
    }


def _canonical_json(value: Any, *, source: Path, field: str) -> str:
    try:
        return json.dumps(
            value,
            sort_keys=True,
            separators=(",", ":"),
            allow_nan=False,
        )
    except (TypeError, ValueError) as exc:
        raise ValueError(
            f"{source}: {field} is not finite canonical JSON: {exc}"
        ) from exc


def validate_campaign(
    pairs: Sequence[tuple[Path, Mapping[str, Any]]],
    *,
    candidate_sha: str | None = None,
    reference_sha: str | None = None,
    expected_campaign: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """Fail closed on provenance, repetition identity, and raw measurements."""

    if not pairs:
        raise ValueError("paired campaign is empty")

    def pair_sha(path: Path, pair: Mapping[str, Any], label: str) -> str:
        value = (pair.get("provenance", {}).get(label, {})).get("sha")
        if (
            not isinstance(value, str)
            or _SHA1_RE.fullmatch(value.lower()) is None
        ):
            raise ValueError(
                f"{path}: {label} provenance SHA must be an exact Git SHA"
            )
        return value.lower()

    first_path, first_pair = pairs[0]
    expected_candidate = (
        candidate_sha
        if candidate_sha is not None
        else pair_sha(first_path, first_pair, "candidate")
    )
    expected_reference = (
        reference_sha
        if reference_sha is not None
        else pair_sha(first_path, first_pair, "reference")
    )
    expected_candidate = _normalized_git_sha(
        expected_candidate, label="candidate",
    )
    expected_reference = _normalized_git_sha(
        expected_reference, label="reference",
    )

    seen_repetitions: dict[tuple[str, str, int], Path] = {}
    cell_signatures: dict[tuple[str, str], tuple[str, int, str]] = {}
    panel_protocols: dict[str, tuple[int, str]] = {}
    for path, pair in pairs:
        panel = pair.get("panel")
        benchmark = pair.get("benchmark")
        repetition = pair.get("repetition")
        if not isinstance(panel, str) or not panel:
            raise ValueError(f"{path}: panel must be a non-empty string")
        if not isinstance(benchmark, str) or not benchmark:
            raise ValueError(f"{path}: benchmark must be a non-empty string")
        if not _integer(repetition) or repetition <= 0:
            raise ValueError(f"{path}: repetition must be a positive integer")

        key = (panel, benchmark, repetition)
        if key in seen_repetitions:
            raise ValueError(
                f"duplicate paired repetition {key}: "
                f"{seen_repetitions[key]} and {path}"
            )
        seen_repetitions[key] = path

        order = pair.get("order")
        if (
            not isinstance(order, (list, tuple))
            or len(order) != 2
            or sorted(order) != ["candidate", "reference"]
        ):
            raise ValueError(
                f"{path}: order must contain candidate and reference exactly once"
            )
        expected_order = (
            ["reference", "candidate"]
            if repetition % 2 else ["candidate", "reference"]
        )
        if list(order) != expected_order:
            raise ValueError(
                f"{path}: repetition {repetition} violates alternating AB/BA "
                f"order; expected {expected_order}"
            )

        actual_candidate = pair_sha(path, pair, "candidate")
        actual_reference = pair_sha(path, pair, "reference")
        if actual_candidate != expected_candidate:
            raise ValueError(
                f"{path}: candidate provenance SHA {actual_candidate} does not "
                f"match report SHA {expected_candidate}"
            )
        if actual_reference != expected_reference:
            raise ValueError(
                f"{path}: reference provenance SHA {actual_reference} does not "
                f"match report SHA {expected_reference}"
            )

        threads = pair.get("threads")
        if not _integer(threads) or threads <= 0:
            raise ValueError(f"{path}: threads must be a positive integer")
        inputs = pair.get("inputs")
        if not isinstance(inputs, Mapping) or not inputs:
            raise ValueError(f"{path}: inputs must be a non-empty fingerprint object")
        for name, record in inputs.items():
            if not isinstance(name, str) or not name or not isinstance(record, Mapping):
                raise ValueError(f"{path}: input fingerprint entries are malformed")
            fingerprint = record.get("sha256")
            if (
                not isinstance(fingerprint, str)
                or _SHA256_RE.fullmatch(fingerprint.lower()) is None
                or not _integer(record.get("size"))
                or record["size"] <= 0
                or not isinstance(record.get("path"), str)
                or not record["path"]
            ):
                raise ValueError(
                    f"{path}: input fingerprint {name!r} is incomplete or malformed"
                )
        modes = pair.get("modes")
        if not isinstance(modes, Mapping):
            raise ValueError(f"{path}: modes must be an object")
        if set(modes) != {"candidate", "reference"} or not all(
            isinstance(modes[label], str) and modes[label]
            for label in ("candidate", "reference")
        ):
            raise ValueError(
                f"{path}: modes must contain non-empty candidate/reference values"
            )
        signature = (
            _canonical_json(inputs, source=path, field="inputs"),
            threads,
            _canonical_json(modes, source=path, field="modes"),
        )
        cell_key = (panel, benchmark)
        previous = cell_signatures.setdefault(cell_key, signature)
        if previous != signature:
            raise ValueError(
                f"{path}: inputs, threads, or modes changed across repetitions "
                f"of {cell_key}"
            )
        panel_protocol = (threads, signature[2])
        previous_protocol = panel_protocols.setdefault(panel, panel_protocol)
        if previous_protocol != panel_protocol:
            raise ValueError(
                f"{path}: threads or modes differ across {panel!r} panel cells"
            )

        raw_profiles = {}
        versions = pair.get("versions")
        if not isinstance(versions, Mapping):
            raise ValueError(f"{path}: versions must be an object")
        for label in ("candidate", "reference"):
            version = versions.get(label)
            if not isinstance(version, Mapping):
                raise ValueError(f"{path}: {label} version evidence is malformed")
            wall = _positive_profile_value(
                pair, label, "wall_clock_seconds", source=path,
            )
            rss = _positive_profile_value(
                pair, label, "peak_rss_mb", source=path,
            )
            raw_profiles[label] = {"wall": wall, "peak_rss": rss}
            profile = version.get("profile")
            if not isinstance(profile, Mapping) or profile.get("exit_code") != 0:
                raise ValueError(
                    f"{path}: {label} profile must be a successful profile object"
                )
            for required in ("user_cpu_seconds", "sys_cpu_seconds"):
                value = profile.get(required)
                if not _number(value) or float(value) < 0:
                    raise ValueError(
                        f"{path}: {label} profile {required!r} must be "
                        "non-negative and finite"
                    )
            validation = version.get("validation")
            if (
                not isinstance(validation, Mapping)
                or not isinstance(validation.get("is_valid"), bool)
                or not _integer(validation.get("n_errors"))
                or validation["n_errors"] < 0
                or validation["is_valid"] != (validation["n_errors"] == 0)
            ):
                raise ValueError(
                    f"{path}: {label} GFF3 validation evidence is malformed"
                )
            fingerprints = version.get("fingerprints")
            if not isinstance(fingerprints, Mapping) or any(
                not isinstance(fingerprints.get(name), str)
                or _SHA256_RE.fullmatch(fingerprints[name].lower()) is None
                for name in ("byte_sha256", "semantic_sha256")
            ):
                raise ValueError(
                    f"{path}: {label} GFF3 fingerprints are malformed"
                )
            if panel == "e2e":
                evaluation_profile = version.get("evaluation_profile")
                if (
                    not isinstance(evaluation_profile, Mapping)
                    or evaluation_profile.get("exit_code") != 0
                ):
                    raise ValueError(
                        f"{path}: {label} E2E evaluation profile is unsuccessful"
                    )
                evaluation_method = version.get("evaluation_method")
                neutral_evaluation = (
                    evaluation_method == "neutral_evaluator_v1"
                    or evaluation_profile.get("method") == "neutral_evaluator_v1"
                )
                if neutral_evaluation:
                    if (
                        evaluation_method != "neutral_evaluator_v1"
                        or evaluation_profile.get("method")
                        != "neutral_evaluator_v1"
                    ):
                        raise ValueError(
                            f"{path}: {label} neutral E2E evaluation method "
                            "is inconsistent"
                        )
                    evaluation_summary = version.get("evaluation_summary")
                    if (
                        not isinstance(evaluation_summary, Mapping)
                        or evaluation_summary.get("format")
                        != "neutral_evaluator_v1"
                    ):
                        raise ValueError(
                            f"{path}: {label} neutral E2E evaluation summary "
                            "is missing or has the wrong format"
                        )
                    required_profile_fields = (
                        "wall_clock_seconds",
                        "peak_rss_mb",
                    )
                else:
                    required_profile_fields = (
                        "wall_clock_seconds",
                        "peak_rss_mb",
                        "user_cpu_seconds",
                        "sys_cpu_seconds",
                    )
                for name in (
                    *required_profile_fields,
                ):
                    value = evaluation_profile.get(name)
                    minimum = 0.0 if name.endswith("cpu_seconds") else 0.0
                    if (
                        not _number(value)
                        or float(value) < minimum
                        or (
                            name in ("wall_clock_seconds", "peak_rss_mb")
                            and float(value) <= 0
                        )
                    ):
                        raise ValueError(
                            f"{path}: {label} E2E evaluation profile {name!r} "
                            "is invalid"
                        )

        ratios = pair.get("ratios")
        if not isinstance(ratios, Mapping):
            raise ValueError(f"{path}: ratios must be an object")
        expected_ratios = {
            "wall": (
                raw_profiles["candidate"]["wall"]
                / raw_profiles["reference"]["wall"]
            ),
            "peak_rss": (
                raw_profiles["candidate"]["peak_rss"]
                / raw_profiles["reference"]["peak_rss"]
            ),
        }
        for name, expected in expected_ratios.items():
            value = ratios.get(name)
            if not _number(value) or float(value) <= 0:
                raise ValueError(
                    f"{path}: raw ratio {name!r} must be positive and finite"
                )
            if not math.isclose(
                float(value), expected, rel_tol=1e-9, abs_tol=1e-12,
            ):
                raise ValueError(
                    f"{path}: raw ratio {name!r} does not match profile data"
                )

    for panel, benchmark in cell_signatures:
        repetitions = sorted(
            repetition
            for seen_panel, seen_benchmark, repetition in seen_repetitions
            if (seen_panel, seen_benchmark) == (panel, benchmark)
        )
        expected = list(range(1, max(repetitions) + 1))
        if repetitions != expected:
            raise ValueError(
                f"{panel}/{benchmark}: repetitions must be contiguous from 1; "
                f"found {repetitions}"
            )
    matrix_complete = False
    normalized_campaign = None
    if expected_campaign is not None:
        normalized_campaign = _normalize_campaign_spec(expected_campaign)
        observed_keys = set(seen_repetitions)
        expected_keys = _expected_campaign_keys(normalized_campaign)
        if observed_keys != expected_keys:
            missing = sorted(expected_keys - observed_keys)
            unexpected = sorted(observed_keys - expected_keys)
            raise ValueError(
                "paired results do not match the expected campaign matrix; "
                f"missing={missing}, unexpected={unexpected}"
            )
        matrix_complete = True

    return {
        "candidate_sha": expected_candidate,
        "reference_sha": expected_reference,
        "n_pairs": len(pairs),
        "n_cells": len(cell_signatures),
        "matrix_complete": matrix_complete,
        "expected_campaign": normalized_campaign,
    }


def _bootstrap_delta(
    values: Sequence[float],
    *,
    seed: int,
    replicates: int,
) -> dict[str, Any] | None:
    clean = [float(value) for value in values if _number(value)]
    if not clean:
        return None
    if len(clean) == 1:
        return {
            "estimate": clean[0],
            "low": clean[0],
            "high": clean[0],
            "replicates": 1,
            "seed": seed,
        }
    import numpy as np

    array = np.asarray(clean, dtype=float)
    rng = np.random.default_rng(seed)
    means = []
    remaining = max(1, int(replicates))
    while remaining:
        batch = min(256, remaining)
        indices = rng.integers(0, len(array), size=(batch, len(array)))
        means.extend(array[indices].mean(axis=1).tolist())
        remaining -= batch
    return {
        "estimate": statistics.fmean(clean),
        "low": _percentile(means, 0.025),
        "high": _percentile(means, 0.975),
        "replicates": len(means),
        "seed": seed,
    }


def _median(values: Iterable[float]) -> float | None:
    clean = [float(value) for value in values if _number(value)]
    return statistics.median(clean) if clean else None


def _log_ratio_cv(values: Sequence[float]) -> float | None:
    """Return the geometric CV of positive repetition-level ratios."""
    clean = [float(value) for value in values if _number(value) and value > 0]
    if len(clean) < 2:
        return None
    logs = [math.log(value) for value in clean]
    spread = statistics.pstdev(logs)
    return math.sqrt(math.expm1(spread * spread))


def _performance_outliers(cells: Sequence[Mapping[str, Any]]) -> list[dict[str, Any]]:
    outliers = []
    for cell in cells:
        wall = [float(value) for value in cell.get("wall_ratios", [])]
        rss = [float(value) for value in cell.get("rss_ratios", [])]
        wall_geomean = math.exp(statistics.fmean(math.log(value) for value in wall))
        rss_geomean = math.exp(statistics.fmean(math.log(value) for value in rss))
        if (
            wall_geomean > 1.10
            or rss_geomean > 1.10
            or max(wall, default=0.0) > 1.25
            or max(rss, default=0.0) > 1.25
            or (cell.get("wall_ratio_cv") or 0.0) > 0.15
            or (cell.get("rss_ratio_cv") or 0.0) > 0.15
        ):
            outliers.append({
                "panel": cell.get("panel"),
                "benchmark": cell.get("benchmark"),
                "repetitions": cell.get("repetitions"),
                "wall_geomean": wall_geomean,
                "rss_geomean": rss_geomean,
                "wall_ratio_cv": cell.get("wall_ratio_cv"),
                "rss_ratio_cv": cell.get("rss_ratio_cv"),
            })
    return outliers


def _group_cells(
    pairs: Sequence[tuple[Path, Mapping[str, Any]]],
) -> dict[tuple[str, str], list[tuple[Path, Mapping[str, Any]]]]:
    grouped: dict[tuple[str, str], list[tuple[Path, Mapping[str, Any]]]] = defaultdict(list)
    for path, pair in pairs:
        grouped[(str(pair["panel"]), str(pair["benchmark"]))].append((path, pair))
    for rows in grouped.values():
        rows.sort(key=lambda item: int(item[1]["repetition"]))
    return grouped


def _quality_is_deterministic(
    rows: Sequence[Mapping[str, Any]],
    *,
    source: Path,
) -> bool:
    return len({
        _canonical_json(row, source=source, field="quality evidence")
        for row in rows
    }) == 1


def _aggregate_quality_repetitions(
    rows: Sequence[Mapping[str, Any]],
) -> dict[str, Any]:
    result = dict(rows[0])
    numeric_fields = {
        *QUALITY_METRICS,
        "mean_pi",
        "median_pi",
        "exon_sn_mean",
        "exon_sp_mean",
    }
    for field in numeric_fields:
        values = [row.get(field) for row in rows]
        if all(_number(value) for value in values):
            result[field] = statistics.median(float(value) for value in values)
    result["repetition_metrics"] = [dict(row) for row in rows]
    return result


def _aggregate_common_pi(
    values: Sequence[Mapping[str, Any] | None],
) -> dict[str, Any] | None:
    present = [value for value in values if value is not None]
    if not present:
        return None
    result = dict(present[0])
    result["mean_delta"] = statistics.median(
        float(value["mean_delta"]) for value in present
    )
    result["repetitions_scored"] = len(present)
    result["repetition_metrics"] = [
        dict(value) if isinstance(value, Mapping) else None
        for value in values
    ]
    return result


def _paired_stable_id_preservation(
    candidate: Mapping[str, Any],
    reference: Mapping[str, Any],
    *,
    source: Path,
) -> dict[str, dict[str, Any]]:
    result = {}
    candidate_types = candidate["stable_id_preservation"]
    reference_types = reference["stable_id_preservation"]
    denominator_fields = (
        "n_reference_records",
        "n_reference_records_with_id",
        "n_reference_ids",
        "applicable",
        "reason",
    )
    for feature_type in STABLE_ID_FEATURE_TYPES:
        candidate_row = candidate_types[feature_type]
        reference_row = reference_types[feature_type]
        mismatched = [
            field for field in denominator_fields
            if candidate_row.get(field) != reference_row.get(field)
        ]
        if mismatched:
            raise ValueError(
                f"{source}: candidate/reference stable ID denominators for "
                f"{feature_type} disagree: {mismatched}"
            )
        applicable = candidate_row["applicable"]
        result[feature_type] = {
            "applicable": applicable,
            "reason": candidate_row["reason"],
            "n_reference_records": candidate_row["n_reference_records"],
            "n_reference_records_with_id": (
                candidate_row["n_reference_records_with_id"]
            ),
            "n_reference_ids": candidate_row["n_reference_ids"],
            "candidate_n_preserved_ids": candidate_row["n_preserved_ids"],
            "candidate_n_output_records": (
                candidate_row["n_output_records"]
            ),
            "candidate_n_output_records_with_id": (
                candidate_row["n_output_records_with_id"]
            ),
            "candidate_n_output_ids": candidate_row["n_output_ids"],
            "candidate_preservation_rate": (
                candidate_row["preservation_rate"]
            ),
            "reference_n_preserved_ids": reference_row["n_preserved_ids"],
            "reference_n_output_records": (
                reference_row["n_output_records"]
            ),
            "reference_n_output_records_with_id": (
                reference_row["n_output_records_with_id"]
            ),
            "reference_n_output_ids": reference_row["n_output_ids"],
            "reference_preservation_rate": (
                reference_row["preservation_rate"]
            ),
            "rate_delta": (
                candidate_row["preservation_rate"]
                - reference_row["preservation_rate"]
                if applicable else None
            ),
        }
    return result


def _campaign_policy_map(
    campaign: Mapping[str, Any],
) -> dict[tuple[str, str], dict[str, str]]:
    specification = campaign.get("expected_campaign")
    if not isinstance(specification, Mapping) or not isinstance(
        specification.get("matrix"), list,
    ):
        return {}
    result = {}
    for case in specification["matrix"]:
        for benchmark in case["ids"]:
            key = (case["panel"], benchmark)
            if key in result:
                raise ValueError(
                    f"campaign policy matrix duplicates {key[0]}/{key[1]}"
                )
            result[key] = {
                "campaign_id": case["id"],
                "baseline_policy": case["baseline_policy"],
                "truth_policy": case["truth_policy"],
            }
    return result


def _truth_f1(
    metrics: Mapping[str, Any] | None,
    level: str,
) -> float | None:
    rows = metrics.get(level) if isinstance(metrics, Mapping) else None
    locus = rows.get("locus") if isinstance(rows, Mapping) else None
    value = locus.get("f1") if isinstance(locus, Mapping) else None
    return float(value) if _number(value) else None


def _target_truth_determinism_projection(
    metrics: Mapping[str, Any],
) -> dict[str, Any]:
    """Retain target-truth semantics while omitting cell-local input paths."""

    projected = dict(metrics)
    inputs = metrics.get("inputs")
    if not isinstance(inputs, Mapping):
        return projected
    projected_inputs = dict(inputs)
    for field in ("prediction_gff", "truth_gff", "source_gff"):
        projected_inputs.pop(field, None)
    mapping = inputs.get("mapping")
    if isinstance(mapping, Mapping):
        projected_mapping = dict(mapping)
        projected_mapping.pop("path", None)
        projected_inputs["mapping"] = projected_mapping
    projected["inputs"] = projected_inputs
    return projected


def aggregate_pairs(
    pairs: Sequence[tuple[Path, Mapping[str, Any]]],
    *,
    seed: int = DEFAULT_BOOTSTRAP_SEED,
    replicates: int = DEFAULT_BOOTSTRAP_REPLICATES,
    candidate_sha: str | None = None,
    reference_sha: str | None = None,
    expected_campaign: Mapping[str, Any] | None = None,
    publication_mode: str = "diagnostic",
    controller_evidence: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    if not _integer(replicates) or replicates <= 0:
        raise ValueError("bootstrap replicates must be a positive integer")
    if not _integer(seed):
        raise ValueError("bootstrap seed must be an integer")
    if publication_mode not in {"diagnostic", "release"}:
        raise ValueError("publication_mode must be diagnostic or release")
    if publication_mode == "release" and (
        expected_campaign is None
        or not isinstance(controller_evidence, Mapping)
        or controller_evidence.get("validated") is not True
    ):
        raise ValueError(
            "release publication requires an expected campaign and validated "
            "controller evidence"
        )
    campaign = validate_campaign(
        pairs,
        candidate_sha=candidate_sha,
        reference_sha=reference_sha,
        expected_campaign=expected_campaign,
    )
    quality_baselines, quality_baseline_artifact = (
        _canonical_quality_baselines()
    )
    policy_map = _campaign_policy_map(campaign)
    panel_concurrency = _controller_panel_concurrency(controller_evidence)
    grouped = _group_cells(pairs)
    cells = []
    warning_totals = {
        "candidate": Counter(),
        "reference": Counter(),
    }
    for (panel, benchmark), rows in sorted(grouped.items()):
        if panel not in PANEL_CONCURRENCY:
            raise ValueError(
                f"no release-report concurrency policy for panel {panel!r}"
            )
        candidate_wall_seconds = [
            _positive_profile_value(
                pair, "candidate", "wall_clock_seconds", source=path,
            )
            for path, pair in rows
        ]
        reference_wall_seconds = [
            _positive_profile_value(
                pair, "reference", "wall_clock_seconds", source=path,
            )
            for path, pair in rows
        ]
        candidate_rss_mb = [
            _positive_profile_value(
                pair, "candidate", "peak_rss_mb", source=path,
            )
            for path, pair in rows
        ]
        reference_rss_mb = [
            _positive_profile_value(
                pair, "reference", "peak_rss_mb", source=path,
            )
            for path, pair in rows
        ]
        wall = [
            candidate / reference
            for candidate, reference in zip(
                candidate_wall_seconds, reference_wall_seconds,
            )
        ]
        rss = [
            candidate / reference
            for candidate, reference in zip(candidate_rss_mb, reference_rss_mb)
        ]
        # Validate every evaluator TSV before selecting the first repetition's
        # deterministic quality values for the cell-level comparison.
        repetition_quality = {
            label: [
                _version_quality(path, pair, label)
                for path, pair in rows
            ]
            for label in ("candidate", "reference")
        }
        first_path, first = rows[0]
        candidate = _aggregate_quality_repetitions(
            repetition_quality["candidate"]
        )
        reference = _aggregate_quality_repetitions(
            repetition_quality["reference"]
        )
        candidate_quality_deterministic = _quality_is_deterministic(
            repetition_quality["candidate"],
            source=first_path,
        )
        reference_quality_deterministic = _quality_is_deterministic(
            repetition_quality["reference"],
            source=first_path,
        )
        stable_id_preservation = _paired_stable_id_preservation(
            candidate,
            reference,
            source=first_path,
        )
        deltas = {
            metric: statistics.median(
                candidate_row[metric] - reference_row[metric]
                for candidate_row, reference_row in zip(
                    repetition_quality["candidate"],
                    repetition_quality["reference"],
                )
            )
            for metric in QUALITY_METRICS
        }
        candidate_hashes = {
            pair["versions"]["candidate"]["fingerprints"]["byte_sha256"]
            for _, pair in rows
        }
        candidate_semantic = {
            pair["versions"]["candidate"]["fingerprints"]["semantic_sha256"]
            for _, pair in rows
        }
        reference_hashes = {
            pair["versions"]["reference"]["fingerprints"]["byte_sha256"]
            for _, pair in rows
        }
        reference_semantic = {
            pair["versions"]["reference"]["fingerprints"]["semantic_sha256"]
            for _, pair in rows
        }
        for label in warning_totals:
            for _, pair in rows[:1]:
                warning_totals[label].update(_warning_inventory(pair, label))
        common_pi_repetitions = [
            _paired_common_pi(path, pair)
            for path, pair in rows
        ]
        common_pi = _aggregate_common_pi(common_pi_repetitions)
        common_pi_deterministic = len({
            _canonical_json(
                value,
                source=first_path,
                field="common-set protein identity",
            )
            for value in common_pi_repetitions
        }) == 1
        target_truth_results = {}
        for label in ("candidate", "reference"):
            repetition_truth = [
                _target_truth_summary(path, pair, label)
                for path, pair in rows
            ]
            present = [value for value in repetition_truth if value is not None]
            if present and len(present) != len(repetition_truth):
                raise ValueError(
                    f"{panel}/{benchmark}: {label} target-truth evidence is "
                    "missing from some repetitions"
                )
            signatures = {
                _canonical_json(
                    _target_truth_determinism_projection(value),
                    source=first_path,
                    field=f"{label} target-truth metrics",
                )
                for value in present
            }
            target_truth_results[label] = {
                "available": bool(present),
                "deterministic": bool(present) and len(signatures) == 1,
                "summary": dict(present[0]) if present else None,
            }
        candidate_truth = target_truth_results["candidate"]["summary"]
        reference_truth = target_truth_results["reference"]["summary"]
        target_truth_results["locus_f1_delta"] = {
            level: (
                _truth_f1(candidate_truth, level)
                - _truth_f1(reference_truth, level)
                if (
                    _truth_f1(candidate_truth, level) is not None
                    and _truth_f1(reference_truth, level) is not None
                )
                else None
            )
            for level in ("gene", "transcript")
        }
        e2e_biology = None
        if panel == "e2e":
            e2e_biology = {}
            for label in ("candidate", "reference"):
                summaries = [
                    pair["versions"][label]
                    for _, pair in rows
                ]
                errors = []
                for (_, pair), summary in zip(rows, summaries):
                    errors.extend(
                        f"repetition {pair['repetition']}: {error}"
                        for error in _e2e_biology_errors(summary)
                    )
                biological_summaries = [
                    summary.get("biological_summary")
                    if isinstance(summary, Mapping) else None
                    for summary in summaries
                ]
                biology_signatures = set()
                deterministic = True
                for summary in summaries:
                    if not isinstance(summary, Mapping):
                        deterministic = False
                        continue
                    try:
                        biology_signatures.add(_canonical_json(
                            _e2e_determinism_document({
                                "biological_summary": summary.get(
                                    "biological_summary"
                                ),
                                "score_summary": summary.get("score_summary"),
                                "evaluation_summary": summary.get(
                                    "evaluation_summary"
                                ),
                            }),
                            source=first_path,
                            field=f"{label} E2E biology",
                        ))
                    except ValueError:
                        deterministic = False
                deterministic = deterministic and len(biology_signatures) == 1
                e2e_biology[label] = {
                    "valid": not errors,
                    "errors": errors,
                    "deterministic": deterministic,
                    "summary": (
                        dict(biological_summaries[0])
                        if isinstance(biological_summaries[0], Mapping)
                        else None
                    ),
                    "repetition_summaries": [
                        {
                            "biological_summary": summary.get(
                                "biological_summary"
                            ),
                            "score_summary": summary.get("score_summary"),
                            "evaluation_summary": summary.get(
                                "evaluation_summary"
                            ),
                        }
                        if isinstance(summary, Mapping) else None
                        for summary in summaries
                    ],
                }
        absolute_quality = {
            "baseline": quality_baselines.get(f"{panel}:{benchmark}"),
        }
        for label, quality in (
            ("candidate", candidate),
            ("reference", reference),
        ):
            absolute_quality[label] = {
                "n_reference_coding": quality["n_reference_coding"],
                "n_recovered_coding": quality["n_recovered_coding"],
                "n_recovered_coding_with_pi": (
                    quality["n_recovered_coding_with_pi"]
                ),
                "completeness_coding": quality["completeness_coding"],
                "completeness_feature_total": (
                    quality["completeness_feature_total"]
                ),
                "covpi": quality["covpi"],
                "mean_pi": quality["mean_pi"],
            }
        cell = {
            "panel": panel,
            "benchmark": benchmark,
            "repetitions": len(rows),
            "orders": [pair["order"] for _, pair in rows],
            "modes": dict(first["modes"]),
            "candidate_wall_seconds_median": statistics.median(
                candidate_wall_seconds
            ),
            "reference_wall_seconds_median": statistics.median(
                reference_wall_seconds
            ),
            "wall_ratios": wall,
            "rss_ratios": rss,
            "wall_ratio_median": statistics.median(wall),
            "wall_ratio_min": min(wall),
            "wall_ratio_max": max(wall),
            "wall_ratio_cv": _log_ratio_cv(wall),
            "rss_ratio_median": statistics.median(rss),
            "rss_ratio_min": min(rss),
            "rss_ratio_max": max(rss),
            "rss_ratio_cv": _log_ratio_cv(rss),
            "candidate_peak_rss_gib": max(candidate_rss_mb) / 1024.0,
            "candidate": candidate,
            "reference": reference,
            "quality_deltas": deltas,
            "common_pi": common_pi,
            "stable_id_preservation": stable_id_preservation,
            "common_pi_deterministic": common_pi_deterministic,
            "e2e_biology": e2e_biology,
            "candidate_quality_deterministic": candidate_quality_deterministic,
            "reference_quality_deterministic": reference_quality_deterministic,
            "candidate_byte_deterministic": len(candidate_hashes) == 1,
            "candidate_semantic_deterministic": len(candidate_semantic) == 1,
            "reference_byte_deterministic": len(reference_hashes) == 1,
            "reference_semantic_deterministic": len(reference_semantic) == 1,
            "candidate_valid": all(
                pair["versions"]["candidate"]["validation"]["is_valid"]
                for _, pair in rows
            ),
            "reference_valid": all(
                pair["versions"]["reference"]["validation"]["is_valid"]
                for _, pair in rows
            ),
            "absolute_quality": absolute_quality,
            "campaign_policy": policy_map.get(
                (panel, benchmark),
                {
                    "campaign_id": None,
                    "baseline_policy": "paired_required",
                    "truth_policy": "none",
                },
            ),
            "target_truth": target_truth_results,
        }
        cells.append(cell)

    panels: dict[str, Any] = {}
    for panel in sorted({cell["panel"] for cell in cells}):
        selected = [cell for cell in cells if cell["panel"] == panel]
        wall_values = [cell["wall_ratio_median"] for cell in selected]
        rss_values = [cell["rss_ratio_median"] for cell in selected]
        candidate_wall_total = sum(
            cell["candidate_wall_seconds_median"] for cell in selected
        )
        reference_wall_total = sum(
            cell["reference_wall_seconds_median"] for cell in selected
        )
        concurrency = panel_concurrency[panel]
        memory_cells = sorted(
            selected,
            key=lambda cell: (
                -cell["candidate_peak_rss_gib"],
                cell["benchmark"],
            ),
        )[:concurrency]
        memory_contributors = [
            {
                "benchmark": cell["benchmark"],
                "candidate_peak_rss_gib": cell["candidate_peak_rss_gib"],
            }
            for cell in memory_cells
        ]
        quality = {}
        for metric in QUALITY_METRICS:
            deltas = [
                cell["quality_deltas"][metric]
                for cell in selected
                if _number(cell["quality_deltas"][metric])
            ]
            quality[metric] = _bootstrap_delta(
                deltas,
                seed=seed + sum(ord(char) for char in f"{panel}:{metric}"),
                replicates=replicates,
            )
        stable_ids = {}
        for feature_type in STABLE_ID_FEATURE_TYPES:
            applicable = [
                cell["stable_id_preservation"][feature_type]
                for cell in selected
                if cell["stable_id_preservation"][feature_type]["applicable"]
            ]
            denominator = sum(
                row["n_reference_ids"] for row in applicable
            )
            candidate_preserved = sum(
                row["candidate_n_preserved_ids"] for row in applicable
            )
            reference_preserved = sum(
                row["reference_n_preserved_ids"] for row in applicable
            )
            stable_ids[feature_type] = {
                "n_applicable_cells": len(applicable),
                "n_reference_ids": denominator,
                "candidate_n_preserved_ids": candidate_preserved,
                "reference_n_preserved_ids": reference_preserved,
                "candidate_preservation_rate": (
                    candidate_preserved / denominator
                    if denominator else None
                ),
                "reference_preservation_rate": (
                    reference_preserved / denominator
                    if denominator else None
                ),
                "rate_delta": (
                    (candidate_preserved - reference_preserved) / denominator
                    if denominator else None
                ),
            }
        panels[panel] = {
            "n_cells": len(selected),
            "wall_ratio": bootstrap_geomean_ratio(
                wall_values, replicates=replicates, seed=seed,
            ),
            "rss_ratio": bootstrap_geomean_ratio(
                rss_values, replicates=replicates, seed=seed + 1,
            ),
            "candidate_wall_seconds_total_of_cell_medians": candidate_wall_total,
            "reference_wall_seconds_total_of_cell_medians": reference_wall_total,
            "wall_total_ratio": candidate_wall_total / reference_wall_total,
            "wall_ratio_median": _median(wall_values),
            "wall_ratio_p95": _percentile(wall_values, 0.95),
            "wall_ratio_worst": max(wall_values) if wall_values else None,
            "rss_ratio_median": _median(rss_values),
            "rss_ratio_p95": _percentile(rss_values, 0.95),
            "rss_ratio_worst": max(rss_values) if rss_values else None,
            "candidate_peak_rss_gib": max(
                (cell["candidate_peak_rss_gib"] for cell in selected),
                default=None,
            ),
            "candidate_concurrency_limit": concurrency,
            "candidate_concurrent_memory_envelope_gib": sum(
                item["candidate_peak_rss_gib"]
                for item in memory_contributors
            ),
            "candidate_concurrent_memory_contributors": memory_contributors,
            "modes": dict(selected[0]["modes"]),
            "quality": quality,
            "stable_id_preservation": stable_ids,
        }
    release_evidence_valid = (
        publication_mode == "release"
        and campaign["matrix_complete"]
        and isinstance(controller_evidence, Mapping)
        and controller_evidence.get("validated") is True
    )
    payload = {
        "schema_version": SCHEMA_VERSION,
        "campaign": campaign,
        "publication": {
            "mode": publication_mode,
            "release_evidence_valid": release_evidence_valid,
            "controller_evidence": (
                dict(controller_evidence)
                if isinstance(controller_evidence, Mapping)
                else None
            ),
        },
        "bootstrap": {"seed": seed, "replicates": replicates},
        "measurement_semantics": {
            "wall_clock_seconds": (
                "Elapsed wall-clock time reported by GNU /usr/bin/time -v for "
                "the profiled command."
            ),
            "peak_rss_mb": (
                "Maximum resident set size reported by GNU /usr/bin/time -v "
                "for the profiled command; this is not a sampled simultaneous "
                "RSS peak for the complete process tree or host."
            ),
            "candidate_concurrent_memory_envelope_gib": (
                "Sum of the largest observed candidate command-level maximum-"
                "RSS values up to the panel admission-slot limit. This is an "
                "engineering admission proxy, not a measured host-memory peak "
                "or a guaranteed process-tree upper bound."
            ),
        },
        "policy": {
            "panel_concurrency": panel_concurrency,
            "candidate_concurrent_memory_limit_gib": MEMORY_ENVELOPE_GIB,
        },
        "quality_baseline_artifact": quality_baseline_artifact,
        "cells": cells,
        "panels": panels,
        "performance_outliers": _performance_outliers(cells),
        "warnings": {
            label: dict(sorted(counter.items()))
            for label, counter in warning_totals.items()
        },
    }
    payload["verdict"] = evaluate_verdict(payload)
    return payload


def evaluate_verdict(metrics: Mapping[str, Any]) -> dict[str, Any]:
    checks = []
    reference_diagnostics = []

    def add(name: str, passed: bool, actual: Any, limit: Any) -> None:
        checks.append({
            "name": name,
            "passed": bool(passed),
            "actual": actual,
            "limit": limit,
        })

    def add_reference_diagnostic(
        name: str, passed: bool, actual: Any, limit: Any,
    ) -> None:
        reference_diagnostics.append({
            "name": name,
            "passed": bool(passed),
            "actual": actual,
            "limit": limit,
        })

    cells = list(metrics.get("cells", []))
    publication = metrics.get("publication")
    release_evidence_valid = (
        isinstance(publication, Mapping)
        and publication.get("mode") == "release"
        and publication.get("release_evidence_valid") is True
    )
    add(
        "complete_immutable_release_campaign",
        release_evidence_valid,
        (
            publication.get("mode")
            if isinstance(publication, Mapping) else "diagnostic"
        ),
        "release mode with complete controller-backed evidence",
    )
    add(
        "all_candidate_outputs_valid",
        bool(cells) and all(cell.get("candidate_valid") is True for cell in cells),
        sum(cell.get("candidate_valid") is True for cell in cells),
        len(cells),
    )
    add(
        "candidate_quality_byte_deterministic",
        bool(cells) and all(
            cell.get("candidate_quality_deterministic") is True
            and cell.get("common_pi_deterministic") is True
            for cell in cells
        ),
        sum(
            cell.get("candidate_quality_deterministic") is True
            and cell.get("common_pi_deterministic") is True
            for cell in cells
        ),
        len(cells),
    )
    add(
        "candidate_outputs_byte_deterministic",
        bool(cells) and all(
            cell.get("candidate_byte_deterministic") is True
            for cell in cells
        ),
        sum(
            cell.get("candidate_byte_deterministic") is True
            for cell in cells
        ),
        len(cells),
    )
    add(
        "reference_quality_byte_deterministic",
        bool(cells) and all(
            cell.get("reference_quality_deterministic") is True
            and cell.get("common_pi_deterministic") is True
            for cell in cells
        ),
        sum(
            cell.get("reference_quality_deterministic") is True
            and cell.get("common_pi_deterministic") is True
            for cell in cells
        ),
        len(cells),
    )
    add_reference_diagnostic(
        "all_reference_outputs_valid",
        bool(cells) and all(cell.get("reference_valid") is True for cell in cells),
        sum(cell.get("reference_valid") is True for cell in cells),
        len(cells),
    )
    add(
        "reference_outputs_byte_deterministic",
        bool(cells) and all(
            cell.get("reference_byte_deterministic") is True
            for cell in cells
        ),
        sum(
            cell.get("reference_byte_deterministic") is True
            for cell in cells
        ),
        len(cells),
    )
    for panel, summary in metrics.get("panels", {}).items():
        wall = summary.get("wall_ratio") or {}
        rss = summary.get("rss_ratio") or {}
        add(
            f"{panel}.wall_gmr_upper",
            _number(wall.get("high")) and wall["high"] <= PANEL_RATIO_LIMIT,
            wall.get("high"),
            PANEL_RATIO_LIMIT,
        )
        add(
            f"{panel}.rss_gmr_upper",
            _number(rss.get("high")) and rss["high"] <= PANEL_RATIO_LIMIT,
            rss.get("high"),
            PANEL_RATIO_LIMIT,
        )
        add(
            f"{panel}.wall_total_ratio",
            _number(summary.get("wall_total_ratio"))
            and summary["wall_total_ratio"] <= PANEL_RATIO_LIMIT,
            summary.get("wall_total_ratio"),
            PANEL_RATIO_LIMIT,
        )
        add(
            f"{panel}.candidate_concurrent_memory_envelope_gib",
            _number(summary.get("candidate_concurrent_memory_envelope_gib"))
            and summary["candidate_concurrent_memory_envelope_gib"]
            <= MEMORY_ENVELOPE_GIB,
            summary.get("candidate_concurrent_memory_envelope_gib"),
            MEMORY_ENVELOPE_GIB,
        )
        add(
            f"{panel}.worst_wall_cell",
            _number(summary.get("wall_ratio_worst"))
            and summary["wall_ratio_worst"] <= CELL_RATIO_LIMIT,
            summary.get("wall_ratio_worst"),
            CELL_RATIO_LIMIT,
        )
        add(
            f"{panel}.worst_rss_cell",
            _number(summary.get("rss_ratio_worst"))
            and summary["rss_ratio_worst"] <= CELL_RATIO_LIMIT,
            summary.get("rss_ratio_worst"),
            CELL_RATIO_LIMIT,
        )
        for metric, interval in summary.get("quality", {}).items():
            if interval is None:
                continue
            add(
                f"{panel}.{metric}.lower",
                interval["low"] >= QUALITY_AGGREGATE_FLOOR,
                interval["low"],
                QUALITY_AGGREGATE_FLOOR,
            )
    for cell in cells:
        absolute = cell.get("absolute_quality") or {}
        baseline = absolute.get("baseline")
        truth_policy = (
            cell.get("campaign_policy", {}).get("truth_policy", "none")
        )
        truth = cell.get("target_truth") or {}
        truth_required = truth_policy == "target_truth_required"
        truth_present = any(
            (truth.get(label) or {}).get("available") is True
            for label in ("candidate", "reference")
        )
        if truth_required or truth_present:
            for label in ("candidate", "reference"):
                result = truth.get(label) or {}
                add(
                    f"{cell['panel']}.{cell['benchmark']}.{label}."
                    "target_truth_available",
                    result.get("available") is True,
                    result.get("available"),
                    True,
                )
                add(
                    f"{cell['panel']}.{cell['benchmark']}.{label}."
                    "target_truth_deterministic",
                    result.get("deterministic") is True,
                    result.get("deterministic"),
                    True,
                )
            candidate_truth = (
                truth.get("candidate", {}).get("summary")
                if isinstance(truth.get("candidate"), Mapping) else None
            )
            for level in ("gene", "transcript"):
                truth_scope = (
                    candidate_truth.get("scope")
                    if isinstance(candidate_truth, Mapping) else None
                )
                scope_groups = (
                    truth_scope.get(f"{level}_groups")
                    if isinstance(truth_scope, Mapping) else None
                )
                minimum_groups = TARGET_TRUTH_MIN_GROUPS.get(
                    cell["panel"], 1,
                )
                add(
                    f"{cell['panel']}.{cell['benchmark']}.candidate."
                    f"target_truth_{level}_scope_groups",
                    _integer(scope_groups)
                    and scope_groups >= minimum_groups,
                    scope_groups,
                    minimum_groups,
                )
                candidate_f1 = _truth_f1(candidate_truth, level)
                delta = (
                    truth.get("locus_f1_delta", {}).get(level)
                    if isinstance(truth.get("locus_f1_delta"), Mapping)
                    else None
                )
                add(
                    f"{cell['panel']}.{cell['benchmark']}.candidate."
                    f"target_truth_{level}_locus_f1",
                    _number(candidate_f1)
                    and candidate_f1 >= TARGET_TRUTH_LOCUS_F1_FLOOR,
                    candidate_f1,
                    TARGET_TRUTH_LOCUS_F1_FLOOR,
                )
                add(
                    f"{cell['panel']}.{cell['benchmark']}."
                    f"target_truth_{level}_locus_f1_delta",
                    _number(delta) and delta >= TARGET_TRUTH_DELTA_FLOOR,
                    delta,
                    TARGET_TRUTH_DELTA_FLOOR,
                )
        for label in ("candidate", "reference"):
            observed = absolute.get(label) or {}
            add(
                f"{cell['panel']}.{cell['benchmark']}.{label}."
                "positive_coding_recovery",
                _integer(observed.get("n_reference_coding"))
                and observed["n_reference_coding"] > 0
                and _integer(observed.get("n_recovered_coding"))
                and observed["n_recovered_coding"] > 0,
                {
                    "reference": observed.get("n_reference_coding"),
                    "recovered": observed.get("n_recovered_coding"),
                },
                "positive recovered coding count",
            )
            add(
                f"{cell['panel']}.{cell['benchmark']}.{label}.valid_pi",
                _integer(observed.get("n_recovered_coding_with_pi"))
                and observed["n_recovered_coding_with_pi"] > 0
                and _number(observed.get("mean_pi"))
                and observed["mean_pi"] > 0
                and _number(observed.get("covpi"))
                and observed["covpi"] > 0,
                {
                    "n": observed.get("n_recovered_coding_with_pi"),
                    "mean_pi": observed.get("mean_pi"),
                    "covpi": observed.get("covpi"),
                },
                "positive scored coding count and finite PI",
            )
            if isinstance(baseline, Mapping):
                for metric, threshold in baseline.items():
                    floor = (
                        threshold.get("floor")
                        if isinstance(threshold, Mapping) else None
                    )
                    value = observed.get(metric)
                    add(
                        f"{cell['panel']}.{cell['benchmark']}.{label}."
                        f"absolute_{metric}",
                        _number(value)
                        and _number(floor)
                        and float(value) >= float(floor),
                        value,
                        threshold,
                    )
        if cell.get("panel") == "e2e":
            biology = cell.get("e2e_biology") or {}
            modes = cell.get("modes") or {}
            add(
                f"e2e.{cell['benchmark']}.candidate_reference_modes_match",
                modes.get("candidate") == modes.get("reference")
                and bool(modes.get("candidate")),
                modes,
                "identical candidate/reference E2E modes",
            )
            for label in ("candidate", "reference"):
                result = biology.get(label) or {}
                add(
                    f"e2e.{cell['benchmark']}.{label}.biological_summary",
                    result.get("valid") is True,
                    result.get("errors", ["missing validation result"]),
                    "schema-valid, non-empty biological summary",
                )
                add(
                    f"e2e.{cell['benchmark']}.{label}.biology_deterministic",
                    result.get("deterministic") is True,
                    result.get("deterministic"),
                    True,
                )
                summary = result.get("summary") or {}
                add(
                    f"e2e.{cell['benchmark']}.{label}."
                    "absolute_feature_completeness",
                    _number(summary.get("feature_completeness"))
                    and summary["feature_completeness"]
                    >= E2E_FEATURE_COMPLETENESS_FLOOR,
                    summary.get("feature_completeness"),
                    E2E_FEATURE_COMPLETENESS_FLOOR,
                )
                add(
                    f"e2e.{cell['benchmark']}.{label}."
                    "absolute_mean_protein_identity",
                    _number(summary.get("mean_protein_identity"))
                    and summary["mean_protein_identity"]
                    >= E2E_MEAN_PI_FLOOR,
                    summary.get("mean_protein_identity"),
                    E2E_MEAN_PI_FLOOR,
                )
        for feature_type, result in (
            cell.get("stable_id_preservation") or {}
        ).items():
            if result.get("applicable") is not True:
                continue
            add(
                f"{cell['panel']}.{cell['benchmark']}."
                f"stable_id_preservation.{feature_type}.no_regression",
                _integer(result.get("candidate_n_preserved_ids"))
                and _integer(result.get("reference_n_preserved_ids"))
                and result["candidate_n_preserved_ids"]
                >= result["reference_n_preserved_ids"],
                {
                    "candidate": result.get("candidate_n_preserved_ids"),
                    "reference": result.get("reference_n_preserved_ids"),
                    "n_reference_ids": result.get("n_reference_ids"),
                },
                "candidate preserved-ID count >= reference preserved-ID count",
            )
        for metric, delta in cell.get("quality_deltas", {}).items():
            if _number(delta):
                add(
                    f"{cell['panel']}.{cell['benchmark']}.{metric}",
                    delta >= QUALITY_CELL_FLOOR,
                    delta,
                    QUALITY_CELL_FLOOR,
                )
        common = cell.get("common_pi")
        add(
            f"{cell['panel']}.{cell['benchmark']}.common_pi_present",
            isinstance(common, Mapping)
            and _integer(common.get("n_common"))
            and common["n_common"] > 0,
            common.get("n_common") if isinstance(common, Mapping) else None,
            "positive common coding set",
        )
        if common is not None:
            add(
                f"{cell['panel']}.{cell['benchmark']}.common_pi",
                common["mean_delta"] >= COMMON_PI_CELL_FLOOR,
                common["mean_delta"],
                COMMON_PI_CELL_FLOOR,
            )
    diagnostic_checks = [
        check for check in checks
        if check["name"] != "complete_immutable_release_campaign"
    ]
    return {
        "passed": bool(checks) and all(check["passed"] for check in checks),
        "diagnostic_passed": (
            bool(diagnostic_checks)
            and all(check["passed"] for check in diagnostic_checks)
        ),
        "reference_diagnostics_passed": (
            bool(reference_diagnostics)
            and all(check["passed"] for check in reference_diagnostics)
        ),
        "reference_diagnostics": reference_diagnostics,
        "checks": checks,
        "n_passed": sum(check["passed"] for check in checks),
        "n_failed": sum(not check["passed"] for check in checks),
    }


def _format_ratio(value: Any) -> str:
    return f"{float(value):.3f}×" if _number(value) else "n/a"


def _format_gib(value: Any) -> str:
    return f"{float(value):.2f} GiB" if _number(value) else "n/a"


def _format_ratio_interval(interval: Any, n_cells: Any) -> str:
    if not isinstance(interval, Mapping):
        return "n/a"
    estimate = _format_ratio(interval.get("estimate"))
    if not _integer(n_cells) or n_cells < 2:
        return f"{estimate} (descriptive)"
    return (
        f"{estimate} [{_format_ratio(interval.get('low'))}, "
        f"{_format_ratio(interval.get('high'))}]"
    )


def _format_preserved(count: Any, total: Any, rate: Any) -> str:
    if not _integer(count) or not _integer(total) or not _number(rate):
        return "n/a"
    return f"{count}/{total} ({float(rate):.5f})"


def render_markdown(metrics: Mapping[str, Any], manifest: Mapping[str, Any]) -> str:
    verdict = metrics["verdict"]
    publication = metrics.get("publication") or {}
    diagnostic = publication.get("mode") != "release"
    panel_concurrency = (
        metrics.get("policy", {}).get("panel_concurrency")
        or PANEL_CONCURRENCY
    )
    concurrency_text = ", ".join(
        f"{panel}={panel_concurrency.get(panel, 'n/a')}"
        for panel in RELEASE_PANELS
    )
    diagnostic_text = (
        "CHECKS PASS" if verdict.get("diagnostic_passed") else "CHECKS FAIL"
    )
    if not verdict.get("reference_diagnostics_passed", True):
        diagnostic_text += "; REFERENCE DIAGNOSTICS FAIL"
    verdict_text = (
        "DIAGNOSTIC ONLY "
        f"({diagnostic_text})"
        if diagnostic
        else ("PASS" if verdict["passed"] else "FAIL")
    )
    lines = [
        "# LiftOn v1.0.11 Release Evaluation",
        "",
        f"**Verdict:** {verdict_text}  ",
        f"**Candidate:** `{manifest['candidate_sha']}`  ",
        f"**Reference:** `{manifest['reference_sha']}`",
        (
            "**Evaluation:** immutable corrected overlay"
            if manifest.get("evaluation_overlays")
            else "**Evaluation:** sealed in-run evaluator artifacts"
        ),
        "",
        "## Performance",
        "",
        "| Panel | Candidate mode | Reference mode | Cells | Wall GMR (95% CI) | "
        "Total wall | RSS GMR (95% CI) | Worst wall | Worst RSS | "
        "Candidate command-RSS slot-sum proxy |",
        "|---|---|---|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for panel, summary in sorted(metrics["panels"].items()):
        wall = summary.get("wall_ratio") or {}
        rss = summary.get("rss_ratio") or {}
        wall_ci = _format_ratio_interval(wall, summary.get("n_cells"))
        rss_ci = _format_ratio_interval(rss, summary.get("n_cells"))
        modes = summary.get("modes") or {}
        lines.append(
            f"| {panel} | {modes.get('candidate', 'n/a')} | "
            f"{modes.get('reference', 'n/a')} | {summary['n_cells']} | {wall_ci} | "
            f"{_format_ratio(summary.get('wall_total_ratio'))} | {rss_ci} | "
            f"{_format_ratio(summary.get('wall_ratio_worst'))} | "
            f"{_format_ratio(summary.get('rss_ratio_worst'))} | "
            f"{_format_gib(summary.get('candidate_concurrent_memory_envelope_gib'))} "
            f"({summary.get('candidate_concurrency_limit', 0)} slots) |"
        )
    lines.extend([
        "",
        "Ratios below 1.0 favor the candidate. Timings are paired and alternate "
        "reference/candidate execution order. Total wall is the sum of candidate "
        "per-cell median seconds divided by the corresponding reference sum. "
        f"The command-RSS slot-sum proxy uses the largest observed GNU-time "
        f"maximum-RSS values up to the configured panel slots "
        f"({concurrency_text}) and must remain at or below "
        f"{MEMORY_ENVELOPE_GIB:.0f} GiB. It is an admission metric, not a "
        f"sampled simultaneous process-tree or host-memory peak.",
        (
            "This is an ad-hoc diagnostic aggregation, not a publishable release "
            "verdict. A release PASS requires an explicit complete campaign "
            "matrix and immutable successful controller evidence."
            if diagnostic else
            "The release verdict is backed by the explicit complete campaign "
            "matrix and immutable successful controller evidence."
        ),
        "",
        "## Annotation correctness",
        "",
        "| Panel | Metric | Delta (95% panel bootstrap CI) |",
        "|---|---|---:|",
    ])
    for panel, summary in sorted(metrics["panels"].items()):
        for metric, interval in summary.get("quality", {}).items():
            if interval is None:
                continue
            lines.append(
                f"| {panel} | {metric} | {interval['estimate']:+.5f} "
                f"[{interval['low']:+.5f}, {interval['high']:+.5f}] |"
            )
    lines.extend([
        "",
        "## Stable feature-ID preservation",
        "",
        "This measures continuity of explicit GFF3 `ID` values for CDS and "
        "exon features. It is copy-suffix-aware and is not a measure of "
        "biological completeness, coordinate accuracy, or sequence identity.",
        "",
        "| Panel | Feature type | Applicable cells | Candidate preserved | "
        "Reference preserved | Rate delta |",
        "|---|---|---:|---:|---:|---:|",
    ])
    for panel, summary in sorted(metrics["panels"].items()):
        stable_ids = summary.get("stable_id_preservation") or {}
        for feature_type in STABLE_ID_FEATURE_TYPES:
            row = stable_ids.get(feature_type) or {}
            delta = row.get("rate_delta")
            delta_text = f"{float(delta):+.5f}" if _number(delta) else "n/a"
            lines.append(
                f"| {panel} | {feature_type} | "
                f"{row.get('n_applicable_cells', 0)} | "
                f"{_format_preserved(row.get('candidate_n_preserved_ids'), row.get('n_reference_ids'), row.get('candidate_preservation_rate'))} | "
                f"{_format_preserved(row.get('reference_n_preserved_ids'), row.get('n_reference_ids'), row.get('reference_preservation_rate'))} | "
                f"{delta_text} |"
            )
    e2e_cells = [
        cell for cell in metrics.get("cells", [])
        if cell.get("panel") == "e2e"
    ]
    if e2e_cells:
        lines.extend([
            "",
            "## End-to-end biological validation",
            "",
            "| Dataset | Version | Status | Reference features | Feature "
            "completeness | Evaluated coding records | Mean protein identity |",
            "|---|---|---|---:|---:|---:|---:|",
        ])
        for cell in e2e_cells:
            for label in ("candidate", "reference"):
                biology = (cell.get("e2e_biology") or {}).get(label) or {}
                summary = biology.get("summary") or {}
                status = "PASS" if biology.get("valid") else "FAIL"
                completeness = summary.get("feature_completeness")
                completeness_text = (
                    f"{float(completeness):.5f}"
                    if _number(completeness) else "n/a"
                )
                identity = summary.get("mean_protein_identity")
                identity_text = (
                    f"{float(identity):.5f}" if _number(identity) else "n/a"
                )
                lines.append(
                    f"| {cell['benchmark']} | {label} | {status} | "
                    f"{summary.get('reference_features', 'n/a')} | "
                    f"{completeness_text} | "
                    f"{summary.get('evaluated_coding_records', 'n/a')} | "
                    f"{identity_text} |"
                )
    failures = [check for check in verdict["checks"] if not check["passed"]]
    lines.extend([
        "",
        "## Gate results",
        "",
        f"- Passed: {verdict['n_passed']}",
        f"- Failed: {verdict['n_failed']}",
    ])
    if failures:
        lines.extend(["", "Failed checks:"])
        for check in failures:
            lines.append(
                f"- `{check['name']}`: actual `{check['actual']}`, "
                f"limit `{check['limit']}`"
            )
    reference_failures = [
        check for check in verdict.get("reference_diagnostics", [])
        if not check["passed"]
    ]
    lines.extend([
        "",
        "## Reference control diagnostics",
        "",
        "Reference validity is reported for baseline transparency and does not "
        "block a candidate-gated release. Reference byte and quality "
        "determinism remain blocking controls.",
        "",
        f"- Status: {'PASS' if verdict.get('reference_diagnostics_passed') else 'FAIL'}",
    ])
    if reference_failures:
        lines.append("- Non-blocking diagnostic failures:")
        for check in reference_failures:
            lines.append(
                f"  - `{check['name']}`: actual `{check['actual']}`, "
                f"limit `{check['limit']}`"
            )
    lines.extend([
        "",
        "The publication manifest hashes every controller success marker, "
        "validated GFF3, arm manifest, GFF3 validation record, evaluator TSV, "
        "tooling source, and paired result used for this report. Other raw logs "
        "and profiles remain in the controller run directories.",
        "",
    ])
    return "\n".join(lines)


def _write_durable_text(path: Path, text: str) -> None:
    with path.open("w", encoding="utf-8") as handle:
        handle.write(text)
        handle.flush()
        os.fsync(handle.fileno())


def _remove_publication_path(path: Path) -> None:
    if path.is_dir() and not path.is_symlink():
        shutil.rmtree(path)
    else:
        path.unlink()


def _validate_publication_destination(
    output_dir: Path,
    roots: Sequence[Path],
) -> None:
    """Reject unsafe replacement targets before moving or deleting anything."""

    from . import build_controller

    repo_root = Path(build_controller.REPO_ROOT).resolve()
    run_roots = [Path(root).resolve() for root in roots]
    if repo_root == output_dir or repo_root.is_relative_to(output_dir):
        raise ValueError(
            f"output directory may not be an ancestor of the repository: "
            f"{output_dir}"
        )
    if any(
        run_root == output_dir
        or run_root.is_relative_to(output_dir)
        or output_dir.is_relative_to(run_root)
        for run_root in run_roots
    ):
        raise ValueError(
            f"output directory may not overlap a controller run root: "
            f"{output_dir}"
        )
    if not (output_dir.exists() or output_dir.is_symlink()):
        return
    if output_dir.is_symlink() or not output_dir.is_dir():
        raise ValueError(
            f"refusing to replace a non-directory or symlink output: {output_dir}"
        )
    expected_names = {"REPORT.md", "manifest.json", "metrics.json"}
    observed_names = {path.name for path in output_dir.iterdir()}
    if observed_names != expected_names or any(
        not (output_dir / name).is_file() for name in expected_names
    ):
        raise ValueError(
            f"refusing to replace unrecognized publication directory: {output_dir}"
        )
    manifest_path = output_dir / "manifest.json"
    try:
        manifest = json.loads(manifest_path.read_text())
    except (OSError, TypeError, ValueError) as exc:
        raise ValueError(
            f"refusing to replace publication with invalid manifest: {exc}"
        ) from exc
    required_manifest_fields = {
        "schema_version",
        "candidate_sha",
        "reference_sha",
        "publication_mode",
        "expected_campaign",
        "controller_evidence",
        "publication_evidence",
        "evaluation_overlays",
        "finalization_host",
        "quality_baseline_artifact",
        "bootstrap",
        "run_roots",
        "pair_results",
        "controller_plans",
    }
    bootstrap = manifest.get("bootstrap") if isinstance(manifest, Mapping) else None
    if (
        not isinstance(manifest, Mapping)
        or set(manifest) != required_manifest_fields
        or manifest.get("schema_version") != SCHEMA_VERSION
        or manifest.get("publication_mode") not in {"diagnostic", "release"}
        or any(
            not isinstance(manifest.get(name), str)
            or _SHA1_RE.fullmatch(manifest[name]) is None
            for name in ("candidate_sha", "reference_sha")
        )
        or not isinstance(bootstrap, Mapping)
        or not _integer(bootstrap.get("seed"))
        or not _integer(bootstrap.get("replicates"))
        or bootstrap["replicates"] <= 0
        or not isinstance(manifest.get("run_roots"), list)
        or not manifest["run_roots"]
        or not all(
            isinstance(path, str) and Path(path).is_absolute()
            for path in manifest["run_roots"]
        )
        or not isinstance(manifest.get("pair_results"), list)
        or not manifest["pair_results"]
        or not isinstance(manifest.get("controller_plans"), list)
        or not (
            manifest.get("evaluation_overlays") is None
            or isinstance(manifest.get("evaluation_overlays"), Mapping)
        )
        or not isinstance(manifest.get("finalization_host"), Mapping)
    ):
        raise ValueError(
            f"refusing to replace publication with invalid manifest: "
            f"{manifest_path}"
        )
    try:
        metrics = json.loads((output_dir / "metrics.json").read_text())
    except (OSError, TypeError, ValueError) as exc:
        raise ValueError(
            f"refusing to replace publication with invalid metrics: {exc}"
        ) from exc
    if (
        not isinstance(metrics, Mapping)
        or metrics.get("schema_version") != SCHEMA_VERSION
        or not isinstance(metrics.get("verdict"), Mapping)
        or not isinstance(metrics["verdict"].get("passed"), bool)
        or not (output_dir / "REPORT.md").read_text().startswith(
            "# LiftOn v1.0.11 Release Evaluation\n"
        )
    ):
        raise ValueError(
            f"refusing to replace publication with invalid report/metrics: "
            f"{output_dir}"
        )


def _write_report_once(
    roots: Sequence[Path],
    output_dir: Path,
    *,
    candidate_sha: str,
    reference_sha: str,
    seed: int = DEFAULT_BOOTSTRAP_SEED,
    replicates: int = DEFAULT_BOOTSTRAP_REPLICATES,
    expected_campaign: Mapping[str, Any] | None = None,
    diagnostic: bool = False,
    evaluation_roots: Sequence[Path] = (),
) -> dict[str, Any]:
    candidate_sha = _normalized_git_sha(candidate_sha, label="candidate")
    reference_sha = _normalized_git_sha(reference_sha, label="reference")
    controller_evidence = None
    normalized_campaign = None
    if diagnostic:
        if expected_campaign is not None:
            normalized_campaign = _normalize_campaign_spec(expected_campaign)
        controller_evidence = _diagnostic_controller_evidence(roots)
        pairs = load_pairs(roots, allow_incomplete=True)
    else:
        if expected_campaign is None:
            raise ValueError(
                "release publication requires an explicit campaign specification; "
                "use diagnostic=True for a non-release report"
            )
        normalized_campaign = _normalize_campaign_spec(expected_campaign)
        pair_paths, controller_evidence = _controller_publication_evidence(
            roots,
            normalized_campaign,
            candidate_sha=candidate_sha,
            reference_sha=reference_sha,
        )
        pairs = []
        for path in pair_paths:
            raw = json.loads(path.read_text())
            if raw.get("schema_version") != SCHEMA_VERSION:
                raise ValueError(f"unsupported paired schema in {path}")
            pairs.append((path, raw))
    if not pairs:
        raise RuntimeError("no pair_result.json files found")
    evaluation_evidence = None
    if evaluation_roots:
        from . import release_rescore

        pairs, evaluation_evidence = release_rescore.apply_overlays(
            pairs, evaluation_roots,
        )
    metrics = aggregate_pairs(
        pairs,
        seed=seed,
        replicates=replicates,
        candidate_sha=candidate_sha,
        reference_sha=reference_sha,
        expected_campaign=(None if diagnostic else normalized_campaign),
        publication_mode="diagnostic" if diagnostic else "release",
        controller_evidence=controller_evidence,
    )
    metrics["publication"]["evaluation_overlays"] = evaluation_evidence
    finalization_host = _finalization_host_snapshot()
    metrics["environment"] = {
        "finalization_host": finalization_host,
    }
    if diagnostic and normalized_campaign is not None:
        metrics["campaign"]["planned_campaign"] = normalized_campaign
        observed_keys = {
            (pair["panel"], pair["benchmark"], int(pair["repetition"]))
            for _, pair in pairs
        }
        expected_keys = _expected_campaign_keys(normalized_campaign)
        metrics["campaign"]["missing_keys"] = [
            list(key) for key in sorted(expected_keys - observed_keys)
        ]
        metrics["campaign"]["unexpected_keys"] = [
            list(key) for key in sorted(observed_keys - expected_keys)
        ]
    publication_evidence = (
        None if diagnostic else _controller_artifact_evidence(roots)
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    manifest = {
        "schema_version": SCHEMA_VERSION,
        "candidate_sha": candidate_sha,
        "reference_sha": reference_sha,
        "publication_mode": "diagnostic" if diagnostic else "release",
        "expected_campaign": None if diagnostic else normalized_campaign,
        "controller_evidence": controller_evidence,
        "publication_evidence": publication_evidence,
        "evaluation_overlays": evaluation_evidence,
        "finalization_host": finalization_host,
        "quality_baseline_artifact": metrics["quality_baseline_artifact"],
        "bootstrap": {"seed": seed, "replicates": replicates},
        "run_roots": [str(Path(root).resolve()) for root in roots],
        "pair_results": [
            {
                "path": str(path.resolve()),
                "sha256": sha256_file(path),
                "panel": pair["panel"],
                "benchmark": pair["benchmark"],
                "repetition": pair["repetition"],
                "transcript_metrics": {
                    label: {
                        "path": str(_transcript_metrics_path(path, pair, label)),
                        "sha256": sha256_file(
                            _transcript_metrics_path(path, pair, label)
                        ),
                    }
                    for label in ("candidate", "reference")
                },
                "outputs": {
                    label: {
                        "path": str(
                            Path(pair["versions"][label].get("output_gff", ""))
                            .resolve()
                        ),
                        "present": Path(
                            pair["versions"][label].get("output_gff", "")
                        ).is_file(),
                        "sha256": (
                            sha256_file(Path(
                                pair["versions"][label]["output_gff"]
                            ))
                            if Path(
                                pair["versions"][label].get("output_gff", "")
                            ).is_file()
                            else None
                        ),
                    }
                    for label in ("candidate", "reference")
                },
            }
            for path, pair in pairs
        ],
        "controller_plans": [
            {
                "path": str((Path(root) / "plan.json").resolve()),
                "sha256": sha256_file(Path(root) / "plan.json"),
            }
            for root in roots
            if (Path(root) / "plan.json").is_file()
        ],
    }
    _write_durable_text(
        output_dir / "metrics.json",
        json.dumps(metrics, indent=2) + "\n",
    )
    _write_durable_text(
        output_dir / "manifest.json",
        json.dumps(manifest, indent=2) + "\n",
    )
    _write_durable_text(
        output_dir / "REPORT.md",
        render_markdown(metrics, manifest),
    )
    return {"metrics": metrics, "manifest": manifest}


def write_report(
    roots: Sequence[Path],
    output_dir: Path,
    *,
    candidate_sha: str,
    reference_sha: str,
    seed: int = DEFAULT_BOOTSTRAP_SEED,
    replicates: int = DEFAULT_BOOTSTRAP_REPLICATES,
    expected_campaign: Mapping[str, Any] | None = None,
    diagnostic: bool = False,
    evaluation_roots: Sequence[Path] = (),
) -> dict[str, Any]:
    """Build and atomically replace a stale-PASS-safe publication directory."""

    output_dir = Path(output_dir).resolve()
    _validate_publication_destination(
        output_dir, [*roots, *evaluation_roots],
    )
    parent = output_dir.parent
    parent.mkdir(parents=True, exist_ok=True)
    quarantined = None
    if output_dir.exists() or output_dir.is_symlink():
        quarantined = Path(tempfile.mkdtemp(
            prefix=f".{output_dir.name}.stale-",
            dir=parent,
        ))
        quarantined.rmdir()
        os.replace(output_dir, quarantined)
    staging = Path(tempfile.mkdtemp(
        prefix=f".{output_dir.name}.tmp-",
        dir=parent,
    ))
    published = False
    try:
        result = _write_report_once(
            roots,
            staging,
            candidate_sha=candidate_sha,
            reference_sha=reference_sha,
            seed=seed,
            replicates=replicates,
            expected_campaign=expected_campaign,
            diagnostic=diagnostic,
            evaluation_roots=evaluation_roots,
        )
        directory_fd = os.open(staging, os.O_RDONLY)
        try:
            os.fsync(directory_fd)
        finally:
            os.close(directory_fd)
        os.replace(staging, output_dir)
        published = True
        parent_fd = os.open(parent, os.O_RDONLY)
        try:
            os.fsync(parent_fd)
        finally:
            os.close(parent_fd)
    except BaseException:
        if published and (output_dir.exists() or output_dir.is_symlink()):
            try:
                os.replace(output_dir, staging)
            except OSError:
                pass
        if staging.exists() or staging.is_symlink():
            _remove_publication_path(staging)
        raise
    if quarantined is not None:
        try:
            _remove_publication_path(quarantined)
        except OSError:
            # The canonical path is already the new durable publication. A
            # best-effort stale-backup cleanup must not invalidate it.
            pass
    return result


def _positive_int_argument(value: str) -> int:
    parsed = int(value)
    if parsed <= 0:
        raise argparse.ArgumentTypeError("must be a positive integer")
    return parsed


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--runs-root", action="append", required=True)
    parser.add_argument(
        "--evaluation-root",
        action="append",
        default=[],
        help=(
            "immutable corrected-evaluation overlay; coverage must exactly "
            "match the selected successful pair results"
        ),
    )
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--candidate-sha", required=True)
    parser.add_argument("--reference-sha", required=True)
    parser.add_argument(
        "--campaign-spec",
        help=(
            "JSON declaring exact subset/full/e2e ids and repetitions; required "
            "for a release verdict"
        ),
    )
    parser.add_argument(
        "--campaign-profile",
        help=(
            "load an exact release matrix from campaign_profiles.json; "
            "mutually exclusive with --campaign-spec"
        ),
    )
    parser.add_argument(
        "--profile-registry",
        help="alternate strict campaign profile registry",
    )
    parser.add_argument(
        "--diagnostic",
        action="store_true",
        help="allow partial/ad-hoc roots but never emit a release PASS",
    )
    parser.add_argument("--seed", type=int, default=DEFAULT_BOOTSTRAP_SEED)
    parser.add_argument(
        "--replicates",
        type=_positive_int_argument,
        default=DEFAULT_BOOTSTRAP_REPLICATES,
    )
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    campaign_inputs = int(bool(args.campaign_spec)) + int(
        bool(args.campaign_profile)
    )
    if not args.diagnostic and campaign_inputs != 1:
        parser.error(
            "provide exactly one of --campaign-spec or --campaign-profile "
            "unless --diagnostic is used"
        )
    if args.diagnostic and campaign_inputs > 1:
        parser.error(
            "provide at most one campaign specification with --diagnostic"
        )
    if args.profile_registry and not args.campaign_profile:
        parser.error("--profile-registry requires --campaign-profile")
    if args.campaign_spec:
        campaign = json.loads(Path(args.campaign_spec).read_text())
    elif args.campaign_profile:
        from . import campaign_profiles

        profile = campaign_profiles.load_profile(
            args.campaign_profile,
            registry=(
                Path(args.profile_registry)
                if args.profile_registry
                else campaign_profiles.DEFAULT_PROFILE_REGISTRY
            ),
        )
        campaign = campaign_profiles.campaign_spec(
            profile,
            legacy_v1=False,
        )
    else:
        campaign = None
    result = write_report(
        [Path(root) for root in args.runs_root],
        Path(args.output_dir),
        candidate_sha=args.candidate_sha,
        reference_sha=args.reference_sha,
        seed=args.seed,
        replicates=args.replicates,
        expected_campaign=campaign,
        diagnostic=args.diagnostic,
        evaluation_roots=[Path(root) for root in args.evaluation_root],
    )
    if not args.diagnostic and not result["metrics"]["verdict"]["passed"]:
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
