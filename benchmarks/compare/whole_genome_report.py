#!/usr/bin/env python3
"""Reduce the biology-first exact-release campaign into frozen evidence."""
from __future__ import annotations

import argparse
import json
import math
import os
import statistics
from collections import defaultdict
from pathlib import Path
from typing import Any, Mapping, Sequence

from . import release_evaluation, release_report, whole_genome_study


SCHEMA_VERSION = 1
METHOD = "lifton-v1.0.11-biology-first-reducer-v1"
SOURCE_METRICS = (
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
TARGET_METRICS = {
    "gene_locus_precision": ("gene", "locus", "precision"),
    "gene_locus_recall": ("gene", "locus", "recall"),
    "gene_locus_f1": ("gene", "locus", "f1"),
    "transcript_locus_precision": ("transcript", "locus", "precision"),
    "transcript_locus_recall": ("transcript", "locus", "recall"),
    "transcript_locus_f1": ("transcript", "locus", "f1"),
    "intron_chain_precision": ("structure", "intron_chain", "precision"),
    "intron_chain_recall": ("structure", "intron_chain", "recall"),
    "intron_chain_f1": ("structure", "intron_chain", "f1"),
    "exon_precision": ("structure", "exon", "precision"),
    "exon_recall": ("structure", "exon", "recall"),
    "cds_precision": ("structure", "CDS", "precision"),
    "cds_recall": ("structure", "CDS", "recall"),
}
SEQUENCE_METRICS = {
    "exact_protein_rate": ("exact_protein", "rate"),
    "mean_protein_identity": ("protein_identity", "mean"),
    "coverage_weighted_protein_identity": (
        "protein_identity", "coverage_weighted",
    ),
    "mean_reciprocal_length_coverage": (
        "protein_identity", "mean_reciprocal_length_coverage",
    ),
    "prediction_orf_valid_rate": ("orf", "prediction_valid_rate"),
    "truth_orf_valid_rate": ("orf", "truth_valid_rate"),
    "both_orf_valid_rate": ("orf", "both_valid_rate"),
}


class ReportError(RuntimeError):
    """Raised when campaign evidence cannot support the report."""


def _fingerprint(document: Mapping[str, Any], label: str) -> None:
    material = dict(document)
    observed = material.pop("fingerprint", None)
    if observed != whole_genome_study.canonical_sha256(material):
        raise ReportError(f"{label} fingerprint is invalid")


def _artifact(
    record: Any,
    label: str,
    *,
    expected_document: Mapping[str, Any] | None = None,
) -> Path:
    if not isinstance(record, Mapping):
        raise ReportError(f"{label}: artifact record is missing")
    path = Path(str(record.get("path", ""))).resolve()
    try:
        observed = whole_genome_study.file_record(path)
    except OSError as exc:
        raise ReportError(f"{label}: artifact cannot be read: {exc}") from exc
    if (
        observed["size"] != record.get("size")
        or observed["sha256"] != record.get("sha256")
    ):
        raise ReportError(f"{label}: artifact changed")
    if expected_document is not None:
        try:
            document = json.loads(path.read_text(encoding="utf-8"))
        except json.JSONDecodeError as exc:
            raise ReportError(f"{label}: artifact JSON is invalid") from exc
        if document != expected_document:
            raise ReportError(f"{label}: embedded and artifact evidence disagree")
    return path


def _process_group_trace(path: Path, label: str) -> dict[str, float | int]:
    samples = 0
    peak_rss = 0.0
    peak_processes = 0
    try:
        with path.open("r", encoding="utf-8", errors="strict") as handle:
            for line_number, raw in enumerate(handle, start=1):
                try:
                    row = json.loads(raw)
                except json.JSONDecodeError as exc:
                    raise ReportError(
                        f"{label}: invalid trace row {line_number}"
                    ) from exc
                rss = row.get("rss_mb") if isinstance(row, Mapping) else None
                processes = (
                    row.get("processes") if isinstance(row, Mapping) else None
                )
                if (
                    isinstance(rss, bool)
                    or not isinstance(rss, (int, float))
                    or not math.isfinite(float(rss))
                    or float(rss) < 0
                    or not isinstance(processes, int)
                    or isinstance(processes, bool)
                    or processes < 0
                ):
                    raise ReportError(
                        f"{label}: malformed trace row {line_number}"
                    )
                samples += 1
                peak_rss = max(peak_rss, float(rss))
                peak_processes = max(peak_processes, processes)
    except OSError as exc:
        raise ReportError(f"{label}: process-group trace cannot be read") from exc
    return {
        "samples": samples,
        "peak_rss_mb": peak_rss,
        "peak_processes": peak_processes,
    }


def _value(document: Mapping[str, Any], path: Sequence[str]) -> Any:
    current: Any = document
    for key in path:
        if not isinstance(current, Mapping):
            return None
        current = current.get(key)
    return current


def _unit(value: Any, label: str) -> float:
    if (
        isinstance(value, bool)
        or not isinstance(value, (int, float))
        or not math.isfinite(float(value))
        or not 0 <= float(value) <= 1
    ):
        raise ReportError(f"{label} must be finite and within [0, 1]")
    return float(value)


def _positive(value: Any, label: str) -> float:
    if (
        isinstance(value, bool)
        or not isinstance(value, (int, float))
        or not math.isfinite(float(value))
        or float(value) <= 0
    ):
        raise ReportError(f"{label} must be finite and positive")
    return float(value)


def _median(values: Sequence[float]) -> float:
    if not values:
        raise ReportError("cannot summarize an empty metric")
    return float(statistics.median(values))


def _geomean(values: Sequence[float]) -> float:
    if not values or any(value <= 0 for value in values):
        raise ReportError("geometric mean needs positive observations")
    return math.exp(statistics.fmean(math.log(value) for value in values))


def _metric_projection(
    document: Mapping[str, Any],
    fields: Mapping[str, Sequence[str]],
    *,
    label: str,
) -> dict[str, float | None]:
    result = {}
    for name, path in fields.items():
        value = _value(document, path)
        result[name] = None if value is None else _unit(value, f"{label}.{name}")
    return result


def _transcript_metrics(version: Mapping[str, Any], label: str) -> dict[str, Any]:
    summary = version.get("summary")
    artifacts = version.get("evaluation_artifacts")
    if not isinstance(summary, Mapping) or not isinstance(artifacts, Mapping):
        raise ReportError(f"{label}: neutral evaluation evidence is missing")
    path = _artifact(
        artifacts.get("transcripts_tsv"), f"{label}.transcripts_tsv",
    )
    metrics = release_report._load_transcript_metrics(
        path, summary.get("n_reference_coding"),
    )
    metrics["completeness_coding"] = _unit(
        summary.get("completeness_coding"),
        f"{label}.completeness_coding",
    )
    metrics["completeness_feature_total"] = _unit(
        summary.get("completeness_feature_total"),
        f"{label}.completeness_feature_total",
    )
    if not math.isclose(
        metrics["completeness_coding"],
        metrics["n_recovered_coding"] / metrics["n_reference_coding"],
        rel_tol=0,
        abs_tol=1e-9,
    ):
        raise ReportError(f"{label}: coding completeness disagrees with TSV")
    return metrics


def _target_metrics(version: Mapping[str, Any], label: str) -> dict[str, Any]:
    target = _value(version, ("summary", "target_truth"))
    if not isinstance(target, Mapping):
        raise ReportError(f"{label}: target-annotation concordance is missing")
    if (
        target.get("method") != "ortholog-scoped-target-coordinate-v1"
        or target.get("parameters", {}).get("id_policy") != "ortholog-map"
        or not target.get("parameters", {}).get("mapping_source_scope_validated")
    ):
        raise ReportError(
            f"{label}: target-annotation concordance is not scope-validated"
        )
    artifacts = version.get("evaluation_artifacts")
    if not isinstance(artifacts, Mapping):
        raise ReportError(f"{label}: target-annotation artifact is missing")
    _artifact(
        artifacts.get("target_truth"),
        f"{label}.target_annotation",
        expected_document=target,
    )
    return {
        "scope": {
            "gene_groups": int(target.get("scope", {}).get("gene_groups", 0)),
            "transcript_groups": int(
                target.get("scope", {}).get("transcript_groups", 0)
            ),
        },
        "coordinate": _metric_projection(
            target, TARGET_METRICS, label=f"{label}.target",
        ),
        "sequence": (
            _metric_projection(
                target["target_sequence"],
                SEQUENCE_METRICS,
                label=f"{label}.target_sequence",
            )
            if isinstance(target.get("target_sequence"), Mapping)
            else None
        ),
    }


def _stable_projection(version: Mapping[str, Any], label: str) -> dict[str, Any]:
    stable = _value(version, ("summary", "stable_id_preservation", "by_type"))
    if not isinstance(stable, Mapping):
        raise ReportError(f"{label}: stable-ID evidence is missing")
    result = {}
    for feature_type in ("CDS", "exon"):
        row = stable.get(feature_type)
        if not isinstance(row, Mapping):
            raise ReportError(f"{label}: {feature_type} stable-ID row is missing")
        result[feature_type] = {
            "applicable": bool(row.get("applicable")),
            "n_reference_ids": int(row.get("n_reference_ids", 0)),
            "n_preserved_ids": int(row.get("n_preserved_ids", 0)),
            "preservation_rate": (
                _unit(
                    row["preservation_rate"],
                    f"{label}.{feature_type}.preservation_rate",
                )
                if row.get("applicable")
                else None
            ),
        }
    return result


def _performance(version: Mapping[str, Any], label: str) -> dict[str, float]:
    profile = version.get("profile")
    if not isinstance(profile, Mapping) or profile.get("exit_code") != 0:
        raise ReportError(f"{label}: successful profile evidence is missing")
    peak_group_rss = _positive(
        profile.get("peak_process_group_rss_mb"),
        f"{label}.process_group_rss",
    )
    trace = Path(str(profile.get("process_group_trace_path", ""))).resolve()
    try:
        trace_record = whole_genome_study.file_record(trace)
    except OSError as exc:
        raise ReportError(f"{label}: process-group trace is missing") from exc
    trace_summary = _process_group_trace(trace, label)
    if (
        trace_record["sha256"] != profile.get("process_group_trace_sha256")
        or not isinstance(profile.get("process_group_sample_count"), int)
        or isinstance(profile.get("process_group_sample_count"), bool)
        or profile["process_group_sample_count"] < 1
        or profile.get("process_group_sample_interval_seconds") != 1.0
        or not isinstance(profile.get("peak_process_group_processes"), int)
        or profile["peak_process_group_processes"] < 1
        or trace_summary["samples"] != profile["process_group_sample_count"]
        or not math.isclose(
            float(trace_summary["peak_rss_mb"]),
            peak_group_rss,
            rel_tol=0,
            abs_tol=1e-9,
        )
        or trace_summary["peak_processes"]
        != profile["peak_process_group_processes"]
    ):
        raise ReportError(f"{label}: process-group trace metadata is invalid")
    return {
        "wall_seconds": _positive(profile.get("wall_clock_seconds"), f"{label}.wall"),
        "command_peak_rss_mb": _positive(
            profile.get("peak_rss_mb"), f"{label}.command_rss",
        ),
        "process_group_peak_rss_mb": _positive(
            peak_group_rss, f"{label}.process_group_rss",
        ),
    }


def _version_row(version: Mapping[str, Any], label: str) -> dict[str, Any]:
    validation = version.get("validation")
    fingerprints = version.get("fingerprints")
    if not isinstance(validation, Mapping) or not isinstance(fingerprints, Mapping):
        raise ReportError(f"{label}: validation or fingerprint evidence is missing")
    for name in ("semantic_sha256", "byte_sha256"):
        value = fingerprints.get(name)
        if (
            not isinstance(value, str)
            or len(value) != 64
            or any(character not in "0123456789abcdef" for character in value)
        ):
            raise ReportError(f"{label}: {name} is malformed")
    if not isinstance(fingerprints.get("feature_counts"), Mapping):
        raise ReportError(f"{label}: feature-count fingerprint is missing")
    return {
        "source": _transcript_metrics(version, label),
        "target": _target_metrics(version, label),
        "stable_ids": _stable_projection(version, label),
        "performance": _performance(version, label),
        "valid": bool(validation.get("is_valid")),
        "validation_errors": int(validation.get("n_errors", 0)),
        "semantic_sha256": fingerprints.get("semantic_sha256"),
        "byte_sha256": fingerprints.get("byte_sha256"),
        "feature_counts": fingerprints.get("feature_counts"),
    }


def _pair_results(campaign_root: str | Path) -> dict[str, list[tuple[Path, dict]]]:
    grouped: dict[str, list[tuple[Path, dict]]] = defaultdict(list)
    for path in sorted(Path(campaign_root).resolve().rglob("pair_result.json")):
        pair = json.loads(path.read_text(encoding="utf-8"))
        if pair.get("schema_version") != release_evaluation.SCHEMA_VERSION:
            raise ReportError(f"{path}: unsupported pair-result schema")
        if pair.get("panel") != "e2e":
            raise ReportError(f"{path}: result is not a raw whole-genome transfer")
        grouped[str(pair.get("benchmark"))].append((path, pair))
    return grouped


def _determinism(rows: Sequence[dict[str, Any]]) -> dict[str, Any]:
    semantic = {row["semantic_sha256"] for row in rows}
    byte = {row["byte_sha256"] for row in rows}
    source = {
        json.dumps(row["source"], sort_keys=True, allow_nan=False)
        for row in rows
    }
    target_coordinate = {
        json.dumps(
            {
                "scope": row["target"]["scope"],
                "coordinate": row["target"]["coordinate"],
            },
            sort_keys=True,
            allow_nan=False,
        )
        for row in rows
    }
    stable_ids = {
        json.dumps(row["stable_ids"], sort_keys=True, allow_nan=False)
        for row in rows
    }
    feature_counts = {
        json.dumps(row["feature_counts"], sort_keys=True, allow_nan=False)
        for row in rows
    }
    return {
        "semantic": len(semantic) == 1,
        "byte": len(byte) == 1,
        "source_metrics": len(source) == 1,
        "target_coordinate_metrics": len(target_coordinate) == 1,
        "stable_id_metrics": len(stable_ids) == 1,
        "feature_counts": len(feature_counts) == 1,
        "passed": all((
            len(semantic) == 1,
            len(byte) == 1,
            len(source) == 1,
            len(target_coordinate) == 1,
            len(stable_ids) == 1,
            len(feature_counts) == 1,
        )),
    }


def _arm_summary(rows: Sequence[dict[str, Any]]) -> dict[str, Any]:
    source = {
        metric: _median([float(row["source"][metric]) for row in rows])
        for metric in SOURCE_METRICS
    }
    coordinate = {
        metric: _median([
            float(row["target"]["coordinate"][metric]) for row in rows
            if row["target"]["coordinate"][metric] is not None
        ])
        for metric in TARGET_METRICS
    }
    sequence_rows = [
        row["target"]["sequence"] for row in rows
        if row["target"]["sequence"] is not None
    ]
    if len(sequence_rows) != 1:
        raise ReportError("exactly one repetition must carry target-sequence metrics")
    performance = {
        key: {
            "median": _median([row["performance"][key] for row in rows]),
            "minimum": min(row["performance"][key] for row in rows),
            "maximum": max(row["performance"][key] for row in rows),
        }
        for key in (
            "wall_seconds", "command_peak_rss_mb", "process_group_peak_rss_mb",
        )
    }
    return {
        "source": source,
        "target_coordinate": coordinate,
        "target_sequence": sequence_rows[0],
        "target_scope": rows[0]["target"]["scope"],
        "stable_ids": rows[0]["stable_ids"],
        "performance": performance,
        "valid_repetitions": sum(row["valid"] for row in rows),
        "validation_errors": sum(row["validation_errors"] for row in rows),
        "determinism": _determinism(rows),
    }


def _cluster_bootstrap(
    values: Mapping[str, Sequence[float]],
    *,
    kind: str,
    seed: int,
    replicates: int,
) -> dict[str, Any]:
    if kind not in {"difference", "ratio"}:
        raise ValueError("bootstrap kind must be difference or ratio")
    pair_values = {
        pair_id: (
            statistics.fmean(observations)
            if kind == "difference"
            else _geomean(observations)
        )
        for pair_id, observations in values.items()
        if observations
    }
    if not pair_values:
        raise ReportError("bootstrap has no pair observations")
    identifiers = sorted(pair_values)
    estimates = [pair_values[pair_id] for pair_id in identifiers]
    estimate = (
        statistics.fmean(estimates)
        if kind == "difference" else _geomean(estimates)
    )
    import numpy as np

    rng = np.random.default_rng(seed)
    array = np.asarray(estimates, dtype=float)
    sampled = []
    remaining = replicates
    while remaining:
        batch = min(256, remaining)
        indices = rng.integers(0, len(array), size=(batch, len(array)))
        if kind == "difference":
            sampled.extend(array[indices].mean(axis=1).tolist())
        else:
            sampled.extend(
                np.exp(np.log(array[indices]).mean(axis=1)).tolist()
            )
        remaining -= batch
    return {
        "estimate": estimate,
        "low": release_evaluation._percentile(sampled, 0.025),
        "high": release_evaluation._percentile(sampled, 0.975),
        "pairs": len(estimates),
        "replicates": replicates,
        "seed": seed,
        "method": "equal-pair cluster bootstrap",
    }


def _summarize_pair(
    pair: Mapping[str, Any],
    rows: Sequence[tuple[Path, dict]],
    repetitions: int,
) -> tuple[dict[str, Any], dict[str, dict[str, list[float]]]]:
    rows = sorted(rows, key=lambda item: item[1].get("repetition", 0))
    if len(rows) != repetitions:
        raise ReportError(
            f"{pair['id']}: expected {repetitions} repetitions, observed {len(rows)}"
        )
    repetitions_seen = sorted(document.get("repetition") for _, document in rows)
    if repetitions_seen != list(range(1, repetitions + 1)):
        raise ReportError(f"{pair['id']}: repetition keys are incomplete")
    arm_rows = {"candidate": [], "reference": []}
    raw_values = {
        "source_deltas": defaultdict(list),
        "target_deltas": defaultdict(list),
        "ratios": defaultdict(list),
    }
    orders = []
    for path, document in rows:
        if (
            document.get("benchmark") != pair["id"]
            or document.get("threads") != 8
            or document.get("modes") != {
                "candidate": "fast", "reference": "fast",
            }
        ):
            raise ReportError(f"{path}: execution contract disagrees with study")
        orders.append(document.get("order"))
        provenance = document.get("provenance", {})
        if (
            provenance.get("candidate", {}).get("sha")
            != "c623f0bddc5b5051a3670a3b0064c52cc61bb719"
            or provenance.get("reference", {}).get("sha")
            != "e503643d8346c600fedabcd3a4dff5c0873a4a37"
        ):
            raise ReportError(f"{path}: exact release SHAs do not match")
        projected = {
            label: _version_row(
                document["versions"][label], f"{pair['id']}.{label}",
            )
            for label in ("candidate", "reference")
        }
        for label in projected:
            arm_rows[label].append(projected[label])
        for metric in SOURCE_METRICS:
            raw_values["source_deltas"][metric].append(
                projected["candidate"]["source"][metric]
                - projected["reference"]["source"][metric]
            )
        for metric in ("gene_locus_f1", "transcript_locus_f1"):
            raw_values["target_deltas"][metric].append(
                projected["candidate"]["target"]["coordinate"][metric]
                - projected["reference"]["target"]["coordinate"][metric]
            )
        raw_values["ratios"]["wall"].append(
            projected["candidate"]["performance"]["wall_seconds"]
            / projected["reference"]["performance"]["wall_seconds"]
        )
        raw_values["ratios"]["process_group_rss"].append(
            projected["candidate"]["performance"]["process_group_peak_rss_mb"]
            / projected["reference"]["performance"]["process_group_peak_rss_mb"]
        )
    expected_orders = [
        ["reference", "candidate"] if repetition % 2 else ["candidate", "reference"]
        for repetition in range(1, repetitions + 1)
    ]
    if orders != expected_orders:
        raise ReportError(f"{pair['id']}: AB/BA order is not complete")
    candidate = _arm_summary(arm_rows["candidate"])
    reference = _arm_summary(arm_rows["reference"])
    result = {
        "id": pair["id"],
        "public_label": pair["public_label"],
        "biological_class": pair["biological_class"],
        "target_annotation": {
            key: pair["target_annotation"][key]
            for key in (
                "provider", "assembly_accession", "assembly_name", "release",
                "release_date", "reported_gene_count",
            )
        },
        "repetitions": repetitions,
        "orders": orders,
        "pair_results": [
            whole_genome_study.file_record(path) for path, _document in rows
        ],
        "candidate": candidate,
        "reference": reference,
        "source_deltas": {
            metric: _median(values)
            for metric, values in raw_values["source_deltas"].items()
        },
        "target_deltas": {
            metric: _median(values)
            for metric, values in raw_values["target_deltas"].items()
        },
        "performance_ratios": {
            metric: {
                "median": _median(values),
                "minimum": min(values),
                "maximum": max(values),
            }
            for metric, values in raw_values["ratios"].items()
        },
    }
    return result, raw_values


def reduce_campaign(
    study_path: str | Path,
    preflight_path: str | Path,
    campaign_root: str | Path,
    canary_path: str | Path,
    *,
    provider_lock: str | Path | None = None,
) -> dict[str, Any]:
    study = whole_genome_study.load_study(study_path)
    preflight = json.loads(Path(preflight_path).read_text(encoding="utf-8"))
    if (
        preflight.get("schema_version") != whole_genome_study.SCHEMA_VERSION
        or preflight.get("kind")
        != "lifton-v1.0.11-biology-study-preflight"
        or not preflight.get("campaign_ready")
        or preflight.get("study") != whole_genome_study.file_record(study_path)
    ):
        raise ReportError("study preflight is not valid and campaign-ready")
    _fingerprint(preflight, "study preflight")
    canary = json.loads(Path(canary_path).read_text(encoding="utf-8"))
    canary_material = dict(canary)
    canary_fingerprint = canary_material.pop("fingerprint", None)
    if (
        canary.get("method") != "lifton-v1.0.11-fast-safe-equivalence-v1"
        or canary.get("exact_sha")
        != "c623f0bddc5b5051a3670a3b0064c52cc61bb719"
        or not canary.get("fast_supported")
        or canary_fingerprint
        != whole_genome_study.canonical_sha256(canary_material)
    ):
        raise ReportError("fast-mode canary did not establish equivalence")
    canary_rows = canary.get("rows")
    if (
        not isinstance(canary_rows, list)
        or not all(isinstance(row, Mapping) for row in canary_rows)
        or {row.get("benchmark") for row in canary_rows}
        != set(("bee", "t2_human_to_gorilla"))
        or not all(row.get("equivalent") is True for row in canary_rows)
    ):
        raise ReportError("fast-mode canary coverage is incomplete")
    grouped = _pair_results(campaign_root)
    expected_ids = [pair["id"] for pair in study["pairs"]]
    unexpected = set(grouped) - set(expected_ids)
    if unexpected:
        raise ReportError(f"campaign contains unexpected IDs: {sorted(unexpected)}")
    pair_results = []
    all_values = {
        "source_deltas": defaultdict(dict),
        "target_deltas": defaultdict(dict),
        "ratios": defaultdict(dict),
    }
    incomplete = []
    repetitions = study["execution"]["repetitions"]
    for pair in study["pairs"]:
        rows = grouped.get(pair["id"], [])
        if len(rows) != repetitions:
            incomplete.append({"id": pair["id"], "observed": len(rows)})
            continue
        result, values = _summarize_pair(pair, rows, repetitions)
        pair_results.append(result)
        for family, metrics in values.items():
            for metric, observations in metrics.items():
                all_values[family][metric][pair["id"]] = observations
    seed = study["statistics"]["bootstrap_seed"]
    replicates = study["statistics"]["bootstrap_replicates"]
    aggregate = {
        "source_deltas": {
            metric: _cluster_bootstrap(
                pair_values,
                kind="difference",
                seed=seed + index,
                replicates=replicates,
            )
            for index, (metric, pair_values) in enumerate(
                all_values["source_deltas"].items()
            )
        },
        "target_deltas": {
            metric: _cluster_bootstrap(
                pair_values,
                kind="difference",
                seed=seed + 100 + index,
                replicates=replicates,
            )
            for index, (metric, pair_values) in enumerate(
                all_values["target_deltas"].items()
            )
        },
        "performance_ratios": {
            metric: _cluster_bootstrap(
                pair_values,
                kind="ratio",
                seed=seed + 200 + index,
                replicates=replicates,
            )
            for index, (metric, pair_values) in enumerate(
                all_values["ratios"].items()
            )
        },
    }
    if pair_results:
        candidate_wall = sum(
            row["candidate"]["performance"]["wall_seconds"]["median"]
            for row in pair_results
        )
        reference_wall = sum(
            row["reference"]["performance"]["wall_seconds"]["median"]
            for row in pair_results
        )
        candidate_rss = sorted((
            row["candidate"]["performance"]["process_group_peak_rss_mb"][
                "maximum"
            ]
            for row in pair_results
        ), reverse=True)
        aggregate["total_work"] = {
            "candidate_wall_seconds": candidate_wall,
            "reference_wall_seconds": reference_wall,
            "wall_ratio": candidate_wall / reference_wall,
        }
        aggregate["candidate_concurrent_memory_proxy_gib"] = (
            sum(candidate_rss[:2]) / 1024.0
        )
    else:
        aggregate["total_work"] = None
        aggregate["candidate_concurrent_memory_proxy_gib"] = None
    limits = study["qualification"]
    gates = []

    def gate(name: str, passed: bool, observed: Any, threshold: Any) -> None:
        gates.append({
            "name": name,
            "passed": bool(passed),
            "observed": observed,
            "threshold": threshold,
        })

    gate(
        "complete paired cohort",
        not incomplete and len(pair_results) == len(expected_ids),
        len(pair_results),
        len(expected_ids),
    )
    gate(
        "candidate validity",
        all(row["candidate"]["valid_repetitions"] == repetitions for row in pair_results),
        sum(row["candidate"]["valid_repetitions"] for row in pair_results),
        len(expected_ids) * repetitions,
    )
    gate(
        "candidate determinism",
        all(row["candidate"]["determinism"]["passed"] for row in pair_results),
        sum(row["candidate"]["determinism"]["passed"] for row in pair_results),
        len(expected_ids),
    )
    for metric, summary in aggregate["source_deltas"].items():
        gate(
            f"aggregate source {metric}",
            summary["low"] >= limits["minimum_aggregate_source_delta_lower_bound"],
            summary["low"],
            limits["minimum_aggregate_source_delta_lower_bound"],
        )
    for row in pair_results:
        gate(
            f"{row['id']} source quality",
            min(row["source_deltas"].values())
            >= limits["minimum_per_pair_source_delta"],
            min(row["source_deltas"].values()),
            limits["minimum_per_pair_source_delta"],
        )
        gate(
            f"{row['id']} candidate target locus F1",
            min(
                row["candidate"]["target_coordinate"]["gene_locus_f1"],
                row["candidate"]["target_coordinate"]["transcript_locus_f1"],
            ) >= limits["minimum_absolute_locus_f1"],
            min(
                row["candidate"]["target_coordinate"]["gene_locus_f1"],
                row["candidate"]["target_coordinate"]["transcript_locus_f1"],
            ),
            limits["minimum_absolute_locus_f1"],
        )
        gate(
            f"{row['id']} target locus delta",
            min(row["target_deltas"].values())
            >= limits["minimum_per_pair_target_delta"],
            min(row["target_deltas"].values()),
            limits["minimum_per_pair_target_delta"],
        )
        gate(
            f"{row['id']} wall ratio",
            row["performance_ratios"]["wall"]["median"]
            <= limits["maximum_per_pair_performance_ratio"],
            row["performance_ratios"]["wall"]["median"],
            limits["maximum_per_pair_performance_ratio"],
        )
        gate(
            f"{row['id']} process-group RSS ratio",
            row["performance_ratios"]["process_group_rss"]["median"]
            <= limits["maximum_per_pair_performance_ratio"],
            row["performance_ratios"]["process_group_rss"]["median"],
            limits["maximum_per_pair_performance_ratio"],
        )
    if pair_results:
        gate(
            "aggregate wall ratio upper confidence bound",
            aggregate["performance_ratios"]["wall"]["high"]
            <= limits["maximum_aggregate_performance_ratio"],
            aggregate["performance_ratios"]["wall"]["high"],
            limits["maximum_aggregate_performance_ratio"],
        )
        gate(
            "aggregate process-group RSS ratio upper confidence bound",
            aggregate["performance_ratios"]["process_group_rss"]["high"]
            <= limits["maximum_aggregate_performance_ratio"],
            aggregate["performance_ratios"]["process_group_rss"]["high"],
            limits["maximum_aggregate_performance_ratio"],
        )
        gate(
            "concurrent candidate memory proxy",
            aggregate["candidate_concurrent_memory_proxy_gib"]
            <= limits["maximum_concurrent_candidate_rss_gib"],
            aggregate["candidate_concurrent_memory_proxy_gib"],
            limits["maximum_concurrent_candidate_rss_gib"],
        )
    provider = None
    if provider_lock is not None:
        provider_document = json.loads(
            Path(provider_lock).read_text(encoding="utf-8")
        )
        if (
            provider_document.get("schema_version") != 1
            or provider_document.get("kind")
            != "lifton-v1.0.11-provider-ortholog-lock"
            or provider_document.get("study")
            != whole_genome_study.file_record(study_path)
            or provider_document.get("preflight")
            != whole_genome_study.file_record(preflight_path)
        ):
            raise ReportError("provider ortholog lock is not bound to this study")
        _fingerprint(provider_document, "provider ortholog lock")
        provider_maps = provider_document.get("maps")
        if not isinstance(provider_maps, Mapping) or set(provider_maps) != set(
            expected_ids
        ):
            raise ReportError("provider ortholog lock has incomplete maps")
        for pair_id, row in provider_maps.items():
            if pair_id not in expected_ids or not isinstance(row, Mapping):
                raise ReportError("provider ortholog lock has unexpected maps")
            if row.get("available"):
                whole_genome_study._resolved_object(
                    Path(preflight_path).resolve().parent,
                    row.get("mapping", {}),
                )
        provider = whole_genome_study.file_record(provider_lock)
    document = {
        "schema_version": SCHEMA_VERSION,
        "method": METHOD,
        "created_at": whole_genome_study.utc_now(),
        "study": whole_genome_study.file_record(study_path),
        "preflight": whole_genome_study.file_record(preflight_path),
        "canary": whole_genome_study.file_record(canary_path),
        "provider_ortholog_lock": provider,
        "coverage": {
            "planned_pairs": len(expected_ids) * repetitions,
            "observed_pairs": len(pair_results) * repetitions,
            "planned_genome_transfers": len(expected_ids),
            "successful_genome_transfers": len(pair_results),
            "incomplete": incomplete,
        },
        "pairs": pair_results,
        "aggregate": aggregate,
        "qualification": {
            "scope": "this seven-transfer cohort only",
            "status": "PASS" if gates and all(row["passed"] for row in gates) else "FAIL",
            "whole_release_claim": "DIAGNOSTIC ONLY",
            "gates": gates,
        },
    }
    material = dict(document)
    material.pop("created_at")
    document["fingerprint"] = whole_genome_study.canonical_sha256(material)
    return document


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--study", type=Path, default=whole_genome_study.DEFAULT_STUDY)
    parser.add_argument("--preflight", type=Path, required=True)
    parser.add_argument("--campaign-root", type=Path, required=True)
    parser.add_argument("--canary", type=Path, required=True)
    parser.add_argument("--provider-lock", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    arguments = build_parser().parse_args(argv)
    try:
        result = reduce_campaign(
            arguments.study,
            arguments.preflight,
            arguments.campaign_root,
            arguments.canary,
            provider_lock=arguments.provider_lock,
        )
        whole_genome_study._atomic_json(arguments.output, result)
    except (OSError, ValueError, ReportError) as exc:
        print(f"whole-genome-report: {exc}", file=os.sys.stderr)
        return 2
    print(json.dumps(result["coverage"], indent=2, sort_keys=True))
    print(json.dumps(result["qualification"], indent=2, sort_keys=True))
    return 0 if result["qualification"]["status"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
