#!/usr/bin/env python
"""Render and verify result-derived tables for the LiftOn technical report.

Legacy comparative tables continue to read ``fourway_results.json`` and
``version_compare.results.json``. Exact-release qualification tables read the
publication bundle's ``metrics.json``, ``manifest.json``, and selected failure
audit. Generated report regions use
markers such as::

    <!-- BEGIN GENERATED: release-performance -->
    ...
    <!-- END GENERATED: release-performance -->

``--check`` fails when a marked region differs from its source evidence, while
``--update`` replaces only marked regions.
"""
from __future__ import annotations

import argparse
from collections import Counter
import difflib
import json
import math
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable, Mapping

HERE = Path(__file__).resolve().parent
FOURWAY = HERE / "fourway_results.json"
VERSION_CMP = HERE / "version_compare.results.json"
CANDIDATE_SHA = "c623f0bddc5b5051a3670a3b0064c52cc61bb719"
REFERENCE_SHA = "e503643d8346c600fedabcd3a4dff5c0873a4a37"

TOOLS = ("liftoff", "miniprot", "lifton_stable", "lifton_devel")
PANEL_ORDER = {"subset": 0, "full": 1, "e2e": 2}
START_MARKER = "<!-- BEGIN GENERATED: {name} -->"
END_MARKER = "<!-- END GENERATED: {name} -->"
MDX_START_MARKER = "{{/* BEGIN GENERATED: {name} */}}"
MDX_END_MARKER = "{{/* END GENERATED: {name} */}}"
MARKER_RE = re.compile(
    r"^(?:<!-- BEGIN GENERATED: (?P<html_name>[a-z0-9-]+) -->|"
    r"\{\/\* BEGIN GENERATED: (?P<mdx_name>[a-z0-9-]+) \*\/\})$",
    re.MULTILINE,
)
PLACEHOLDER_RE = re.compile(r"@@[A-Z0-9_]+@@")

TIERS = [
    ("same_species", "same"),
    ("cross_species", "cross"),
    ("close_cross_species", "close"),
    ("distant_cross_species", "distant"),
    ("very_distant_cross_species", "very distant"),
]

PRETTY = {
    "arabidopsis": "arabidopsis",
    "bee": "bee",
    "rice": "rice",
    "t1_maize_b73_to_mo17": "maize",
    "t1_tomato_microtom_to_heinz": "tomato",
    "drosophila": "drosophila",
    "t2_human_to_gorilla": "human→gorilla",
    "t2_mouse_to_caroli": "mouse→caroli",
    "t2_tomato_to_potato": "tomato→potato",
    "t3_dog_to_cat": "dog→cat",
    "t3_human_to_macaque": "human→macaque",
    "t3_human_to_marmoset": "human→marmoset",
    "arabidopsis_to_rice": "arabidopsis→rice",
    "human_to_zebrafish": "human→zebrafish",
    "t4_human_to_chicken": "human→chicken",
    "t4_human_to_xenopus": "human→xenopus",
    "t4_drosophila_to_bee": "drosophila→bee",
}

LEGACY_FULL_FAILURES = (
    "arabidopsis",
    "rice",
    "t1_maize_b73_to_mo17",
    "t1_tomato_microtom_to_heinz",
    "t2_tomato_to_potato",
)


@dataclass(frozen=True)
class Labels:
    candidate: str = "LiftOn v1.0.11"
    reference: str = "LiftOn v1.0.8"
    legacy_candidate: str = "LiftOn candidate (7fc8419)"


def load(path: Path) -> Any:
    if not path.is_file():
        raise ValueError(f"missing results document: {path}")
    try:
        return json.loads(path.read_text())
    except json.JSONDecodeError as exc:
        raise ValueError(f"invalid JSON in {path}: {exc}") from exc


def full_cells(db: Mapping[str, Any]) -> list[tuple[str, Mapping[str, Any]]]:
    """Return full-genome cells in report divergence order."""

    cells = {
        key.split(":", 1)[1]: value
        for key, value in db.items()
        if key.startswith("full:")
    }
    order = []
    for divergence_class, _label in TIERS:
        order.extend(sorted(
            (
                benchmark
                for benchmark, cell in cells.items()
                if cell.get("divergence_class") == divergence_class
            ),
            key=lambda benchmark: PRETTY.get(benchmark, benchmark),
        ))
    order.extend(benchmark for benchmark in cells if benchmark not in order)
    return [(benchmark, cells[benchmark]) for benchmark in order]


def _fmt(value: Any, digits: int = 5, dash: str = "—") -> str:
    if value is None:
        return dash
    if isinstance(value, bool):
        return "yes" if value else "no"
    if isinstance(value, float):
        return f"{value:.{digits}f}"
    return f"{value:,}" if isinstance(value, int) else str(value)


def _signed(value: Any, digits: int = 5) -> str:
    return "—" if not _number(value) else f"{float(value):+.{digits}f}"


def _ratio(value: Any) -> str:
    return "—" if not _number(value) else f"{float(value):.3f}×"


def _inverse_ratio(value: Any) -> str:
    if not _number(value) or float(value) <= 0:
        return "—"
    return f"{1.0 / float(value):.2f}×"


def _preserved(count: Any, denominator: Any, rate: Any) -> str:
    if not isinstance(count, int) or not isinstance(denominator, int) or not _number(rate):
        return "—"
    return f"{count:,}/{denominator:,} ({float(rate):.5f})"


def _ci(interval: Any, n_units: Any = None) -> str:
    if not isinstance(interval, Mapping) or not all(
        _number(interval.get(key)) for key in ("estimate", "low", "high")
    ):
        return "—"
    estimate = f"{float(interval['estimate']):.3f}×"
    if isinstance(n_units, int) and not isinstance(n_units, bool) and n_units < 2:
        return f"{estimate} *(descriptive; n={n_units} cell)*"
    return (
        f"{estimate} "
        f"[{float(interval['low']):.3f}, {float(interval['high']):.3f}]"
    )


def _delta_ci(interval: Any, n_units: Any = None) -> str:
    if not isinstance(interval, Mapping) or not all(
        _number(interval.get(key)) for key in ("estimate", "low", "high")
    ):
        return "—"
    estimate = f"{float(interval['estimate']):+.5f}"
    if isinstance(n_units, int) and not isinstance(n_units, bool) and n_units < 2:
        return f"{estimate} *(descriptive; n={n_units} cell)*"
    return (
        f"{estimate} "
        f"[{float(interval['low']):+.5f}, {float(interval['high']):+.5f}]"
    )


def _number(value: Any) -> bool:
    return (
        isinstance(value, (int, float))
        and not isinstance(value, bool)
        and math.isfinite(float(value))
    )


def _tier(cell: Mapping[str, Any]) -> str:
    divergence_class = cell.get("divergence_class")
    return dict(TIERS).get(divergence_class, divergence_class or "—")


def _legacy_header(labels: Labels) -> tuple[str, str]:
    return labels.reference, labels.legacy_candidate


def table2(db: Mapping[str, Any], _vc: Mapping[str, Any], _labels: Labels) -> str:
    rows = [
        "| Whole genome | Divergence | n coding | n features |",
        "|---|---|---:|---:|",
    ]
    for benchmark, cell in full_cells(db):
        rows.append(
            f"| {PRETTY.get(benchmark, benchmark)} | {_tier(cell)} | "
            f"{_fmt(cell.get('n_reference_coding'))} | "
            f"{_fmt(cell.get('n_reference_total'))} |"
        )
    return "\n".join(rows)


def table3(db: Mapping[str, Any], _vc: Mapping[str, Any], labels: Labels) -> str:
    reference, candidate = _legacy_header(labels)
    rows = [
        f"| Full genome | {reference} recovered | {candidate} recovered | "
        f"gain | {reference} completeness | {candidate} completeness |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    full_by_benchmark = dict(full_cells(db))
    for benchmark in LEGACY_FULL_FAILURES:
        cell = full_by_benchmark.get(benchmark)
        if cell is None:
            continue
        recovered = cell.get("n_recovered_coding") or {}
        completeness = cell.get("completeness_coding") or {}
        stable = recovered.get("lifton_stable")
        development = recovered.get("lifton_devel")
        gain = (
            f"{development - stable:+,}"
            if isinstance(stable, int) and isinstance(development, int)
            else "—"
        )
        stable_cell = _fmt(stable) if stable is not None else "crash (no output)"
        stable_completeness = completeness.get("lifton_stable")
        stable_completeness_cell = (
            "crash (no output)"
            if stable_completeness is None
            else f"{stable_completeness * 100:.2f}% (partial output)"
        )
        development_completeness_value = completeness.get("lifton_devel")
        development_completeness = (
            f"{development_completeness_value * 100:.2f}%"
            if _number(development_completeness_value)
            else "—"
        )
        rows.append(
            f"| {PRETTY.get(benchmark, benchmark)} | {stable_cell} | "
            f"{_fmt(development)} | {gain} | {stable_completeness_cell} | "
            f"{development_completeness} |"
        )
    return "\n".join(rows)


def table4(_db: Mapping[str, Any], vc: Mapping[str, Any], labels: Labels) -> str:
    reference, candidate = _legacy_header(labels)
    rows = [
        f"| Benchmark | wall {reference} (s) | wall {candidate} (s) | speedup | "
        f"RSS {reference} (MB) | RSS {candidate} (MB) | RSS reduction |",
        "|---|---:|---:|---:|---:|---:|---:|",
    ]
    for key, record in vc.items():
        if not key.startswith("controlled:"):
            continue
        benchmark = key.split(":", 1)[1]
        wall = record.get("wall_s") or {}
        rss = record.get("peak_rss_mb") or {}
        stable_wall, development_wall = wall.get("stable"), wall.get("devel")
        stable_rss, development_rss = rss.get("stable"), rss.get("devel")
        if all(value is not None for value in (
            stable_wall, development_wall, stable_rss, development_rss,
        )):
            rows.append(
                f"| {benchmark} | {stable_wall:.1f} | {development_wall:.1f} | "
                f"{stable_wall / development_wall:.2f}× | {stable_rss:,.0f} | "
                f"{development_rss:,.0f} | {stable_rss / development_rss:.2f}× |"
            )
        else:
            rows.append(f"| {benchmark} | — | — | — | — | — | — |")
    return "\n".join(rows)


def table5(db: Mapping[str, Any], _vc: Mapping[str, Any], labels: Labels) -> str:
    reference, candidate = _legacy_header(labels)
    rows = [
        "| Whole genome | Divergence | Liftoff | *miniprot (evidence)* | "
        f"{reference} | **{candidate}** | candidate − Liftoff |",
        "|---|---|---:|---:|---:|---:|---:|",
    ]
    for benchmark, cell in full_cells(db):
        identity = cell.get("mean_pi") or {}
        liftoff = identity.get("liftoff")
        development = identity.get("lifton_devel")
        lead = _signed(development - liftoff) if liftoff is not None and development is not None else "—"
        rows.append(
            f"| {PRETTY.get(benchmark, benchmark)} | {_tier(cell)} | {_fmt(liftoff)} | "
            f"{_fmt(identity.get('miniprot'))} | {_fmt(identity.get('lifton_stable'))} | "
            f"**{_fmt(development)}** | {lead} |"
        )
    return "\n".join(rows)


def table7(db: Mapping[str, Any], _vc: Mapping[str, Any], labels: Labels) -> str:
    reference, candidate = _legacy_header(labels)
    rows = [
        "| Whole genome | Divergence | Liftoff | *miniprot (evidence)* | "
        f"{reference} | **{candidate}** |",
        "|---|---|---:|---:|---:|---:|",
    ]
    for benchmark, cell in full_cells(db):
        completeness = cell.get("completeness_coding") or {}
        values = []
        for tool in TOOLS:
            value = completeness.get(tool)
            if value is None:
                values.append("crash")
            elif tool == "lifton_devel":
                values.append(f"**{value * 100:.2f}%**")
            else:
                values.append(f"{value * 100:.2f}%")
        rows.append(
            f"| {PRETTY.get(benchmark, benchmark)} | {_tier(cell)} | "
            + " | ".join(values)
            + " |"
        )
    return "\n".join(rows)


def table9(db: Mapping[str, Any], _vc: Mapping[str, Any], labels: Labels) -> str:
    reference, candidate = _legacy_header(labels)
    rows = [
        "| Whole genome | Divergence | Liftoff | *miniprot (evidence)* | "
        f"{reference} | **{candidate}** |",
        "|---|---|---:|---:|---:|---:|",
    ]
    for benchmark, cell in full_cells(db):
        validity = cell.get("validity") or {}
        values = []
        for tool in TOOLS:
            errors = (validity.get(tool) or {}).get("n_errors")
            if errors is None:
                values.append("—")
            elif tool == "lifton_devel":
                values.append(f"**{errors}**")
            else:
                values.append(str(errors))
        rows.append(
            f"| {PRETTY.get(benchmark, benchmark)} | {_tier(cell)} | "
            + " | ".join(values)
            + " |"
        )
    return "\n".join(rows)


TABLES: dict[str, tuple[str, Callable[[Mapping[str, Any], Mapping[str, Any], Labels], str]]] = {
    "2": ("The seventeen whole-genome lift-overs", table2),
    "3": ("Full-genome recovery where v1.0.8 fails", table3),
    "4": ("Legacy controlled comparison", table4),
    "5": ("Legacy whole-genome mean protein identity", table5),
    "7": ("Legacy whole-genome coding completeness", table7),
    "9": ("Legacy validator error counts", table9),
}


def _validate_release_metrics(metrics: Any) -> Mapping[str, Any]:
    if not isinstance(metrics, Mapping) or metrics.get("schema_version") != 4:
        raise ValueError("release metrics must be a schema-version 4 object")
    for key in ("campaign", "panels", "cells", "verdict"):
        if key not in metrics:
            raise ValueError(f"release metrics are missing {key!r}")
    if not isinstance(metrics["panels"], Mapping) or not metrics["panels"]:
        raise ValueError("release metrics contain no panel summaries")
    if not isinstance(metrics["cells"], list) or not metrics["cells"]:
        raise ValueError("release metrics contain no cells")
    campaign = metrics["campaign"]
    if not isinstance(campaign, Mapping):
        raise ValueError("release metrics campaign is malformed")
    if campaign.get("n_cells") != len(metrics["cells"]):
        raise ValueError("release metrics n_cells disagrees with cells")
    seen = set()
    for cell in metrics["cells"]:
        if not isinstance(cell, Mapping):
            raise ValueError("release metrics contain a malformed cell")
        key = (cell.get("panel"), cell.get("benchmark"))
        repetitions = cell.get("repetitions")
        if (
            key in seen
            or key[0] not in metrics["panels"]
            or not isinstance(key[1], str)
            or not isinstance(repetitions, int)
            or isinstance(repetitions, bool)
            or repetitions < 1
        ):
            raise ValueError(f"release metrics contain malformed or duplicate cell {key!r}")
        seen.add(key)
    return metrics


def _planned_pair_keys(metrics: Mapping[str, Any]) -> set[tuple[str, str, int]]:
    planned = (metrics["campaign"].get("planned_campaign") or {}).get("panels")
    if not isinstance(planned, Mapping) or not planned:
        raise ValueError("release metrics contain no planned panel matrix")
    keys = set()
    for panel, row in planned.items():
        if not isinstance(row, Mapping):
            raise ValueError(f"planned panel {panel!r} is malformed")
        ids = row.get("ids")
        repetitions = row.get("repetitions")
        if (
            not isinstance(ids, list)
            or not ids
            or len(set(ids)) != len(ids)
            or not all(isinstance(benchmark, str) and benchmark for benchmark in ids)
            or not isinstance(repetitions, int)
            or isinstance(repetitions, bool)
            or repetitions < 1
        ):
            raise ValueError(f"planned panel {panel!r} has invalid ids/repetitions")
        keys.update(
            (panel, benchmark, repetition)
            for benchmark in ids
            for repetition in range(1, repetitions + 1)
        )
    return keys


def _observed_pair_keys(metrics: Mapping[str, Any]) -> set[tuple[str, str, int]]:
    return {
        (cell["panel"], cell["benchmark"], repetition)
        for cell in metrics["cells"]
        for repetition in range(1, int(cell["repetitions"]) + 1)
    }


def validate_failure_audit(document: Any, metrics: Mapping[str, Any]) -> Mapping[str, Any]:
    """Validate the sealed deterministic-control-failure classification."""

    metrics = _validate_release_metrics(metrics)
    campaign = metrics["campaign"]
    if campaign.get("candidate_sha") != CANDIDATE_SHA:
        raise ValueError("failure-audit tables require exact v1.0.11 candidate SHA")
    if campaign.get("reference_sha") != REFERENCE_SHA:
        raise ValueError("failure-audit tables require exact v1.0.8 reference SHA")
    if not isinstance(document, Mapping) or document.get("schema_version") != 1:
        raise ValueError("failure audit must be a schema-version 1 object")
    if document.get("classification") != "deterministic historical-reference failure":
        raise ValueError("failure audit has an unexpected classification")
    roots = document.get("selected_roots")
    cells = document.get("cells")
    invariants = document.get("invariants")
    if (
        not isinstance(roots, list)
        or not roots
        or len(set(roots)) != len(roots)
        or not all(isinstance(root, str) and root for root in roots)
    ):
        raise ValueError("failure audit selected roots are malformed")
    if not isinstance(cells, list) or not cells or not isinstance(invariants, Mapping):
        raise ValueError("failure audit cells or invariants are malformed")
    audit_keys = set()
    counts = Counter()
    candidate_first = 0
    for cell in cells:
        if not isinstance(cell, Mapping):
            raise ValueError("failure audit contains a malformed cell")
        key = (cell.get("panel"), cell.get("benchmark"), cell.get("repetition"))
        if (
            key in audit_keys
            or key[0] not in {"full", "e2e"}
            or not isinstance(key[1], str)
            or not isinstance(key[2], int)
            or isinstance(key[2], bool)
        ):
            raise ValueError(f"failure audit contains malformed or duplicate key {key!r}")
        audit_keys.add(key)
        counts[(key[0], key[1])] += 1
        if cell.get("root") not in roots:
            raise ValueError(f"failure audit key {key!r} names an unselected root")
        if (
            cell.get("returncode") != 1
            or cell.get("watchdog_reason") is not None
            or cell.get("feature_not_found") is not True
            or cell.get("exception_string_typeerror") is not True
            or not _number(cell.get("elapsed_seconds"))
        ):
            raise ValueError(f"failure audit key {key!r} is not a deterministic v1.0.8 failure")
        candidate_completed = cell.get("candidate_completed_before_control_failure")
        if not isinstance(candidate_completed, bool):
            raise ValueError(f"failure audit key {key!r} has malformed candidate-first state")
        candidate_first += int(candidate_completed)
    if document.get("failed_pair_count") != len(cells):
        raise ValueError("failure audit failed_pair_count disagrees with cells")
    missing = campaign.get("missing_keys")
    if not isinstance(missing, list) or any(
        not isinstance(key, list) or len(key) != 3 for key in missing
    ):
        raise ValueError("release metrics missing_keys are malformed")
    if audit_keys != {tuple(key) for key in missing}:
        raise ValueError("failure audit keys disagree with release metrics missing_keys")
    planned_keys = _planned_pair_keys(metrics)
    observed_keys = _observed_pair_keys(metrics)
    if observed_keys & audit_keys or observed_keys | audit_keys != planned_keys:
        raise ValueError("observed and failed pairs do not partition the planned matrix")
    if campaign.get("n_pairs") != len(observed_keys):
        raise ValueError("release metrics n_pairs disagrees with observed repetitions")
    expected_counts = [
        {"panel": panel, "benchmark": benchmark, "failed_repetitions": count}
        for (panel, benchmark), count in sorted(counts.items())
    ]
    if document.get("failed_case_counts") != expected_counts:
        raise ValueError("failure audit case counts disagree with cells")
    for name in (
        "all_returncode_1",
        "all_watchdog_reason_null",
        "all_feature_not_found",
        "all_exception_string_typeerror",
    ):
        if invariants.get(name) is not True:
            raise ValueError(f"failure audit invariant {name!r} does not hold")
    if invariants.get("candidate_first_outputs_preserved") != candidate_first:
        raise ValueError("failure audit candidate-first count disagrees with cells")
    if len(cells) != 28 or len(counts) != 7 or candidate_first != 14:
        raise ValueError("failure audit disagrees with the sealed exact-release boundary")
    return document


def validate_release_manifest(document: Any, metrics: Mapping[str, Any]) -> Mapping[str, Any]:
    """Cross-check the report manifest without rewriting or rehashing sealed files."""

    metrics = _validate_release_metrics(metrics)
    if not isinstance(document, Mapping) or document.get("schema_version") != 4:
        raise ValueError("release manifest must be a schema-version 4 object")
    campaign = metrics["campaign"]
    if document.get("candidate_sha") != campaign.get("candidate_sha"):
        raise ValueError("release manifest candidate SHA disagrees with metrics")
    if document.get("reference_sha") != campaign.get("reference_sha"):
        raise ValueError("release manifest reference SHA disagrees with metrics")
    if document.get("publication_mode") != (metrics.get("publication") or {}).get("mode"):
        raise ValueError("release manifest publication mode disagrees with metrics")
    pair_results = document.get("pair_results")
    if not isinstance(pair_results, list):
        raise ValueError("release manifest pair results are malformed")
    manifest_keys = {
        (record.get("panel"), record.get("benchmark"), record.get("repetition"))
        for record in pair_results
        if isinstance(record, Mapping)
    }
    if len(manifest_keys) != len(pair_results) or manifest_keys != _observed_pair_keys(metrics):
        raise ValueError("release manifest pair results disagree with metrics cells")
    if len(pair_results) != campaign.get("n_pairs"):
        raise ValueError("release manifest pair count disagrees with metrics")
    if document.get("bootstrap") != metrics.get("bootstrap"):
        raise ValueError("release manifest bootstrap settings disagree with metrics")
    if document.get("quality_baseline_artifact") != metrics.get("quality_baseline_artifact"):
        raise ValueError("release manifest quality baseline disagrees with metrics")
    manifest_overlays = document.get("evaluation_overlays") or {}
    metric_overlays = (metrics.get("publication") or {}).get("evaluation_overlays") or {}
    if (
        manifest_overlays.get("validated") != metric_overlays.get("validated")
        or manifest_overlays.get("pair_count") != metric_overlays.get("pair_count")
    ):
        raise ValueError("release manifest overlays disagree with metrics")
    manifest_roots = set(document.get("run_roots") or [])
    metric_roots = {
        row.get("root")
        for row in ((metrics.get("publication") or {}).get("controller_evidence") or {}).get("roots", [])
        if isinstance(row, Mapping)
    }
    if manifest_roots != metric_roots:
        raise ValueError("release manifest controller roots disagree with metrics")
    if document.get("finalization_host") != (metrics.get("environment") or {}).get("finalization_host"):
        raise ValueError("release manifest finalization host disagrees with metrics")
    return document


def release_provenance(metrics: Mapping[str, Any], _db: Mapping[str, Any], _vc: Mapping[str, Any],
                       labels: Labels) -> str:
    campaign = metrics["campaign"]
    expected = (
        campaign.get("expected_campaign")
        or campaign.get("planned_campaign")
        or {}
    )
    bootstrap = metrics.get("bootstrap") or {}
    verdict = metrics.get("verdict") or {}
    publication = metrics.get("publication") or {}
    publication_verdict = (
        "PASS"
        if publication.get("mode") == "release" and verdict.get("passed")
        else "FAIL"
        if publication.get("mode") == "release"
        else "DIAGNOSTIC ONLY"
    )
    rows = [
        "| Evidence | Value |",
        "|---|---|",
        f"| Candidate | {labels.candidate} `{campaign.get('candidate_sha', '—')}` |",
        f"| Reference | {labels.reference} `{campaign.get('reference_sha', '—')}` |",
        f"| Profile | `{expected.get('profile_id', '—')}` |",
        f"| Profile digest | `{expected.get('profile_digest', '—')}` |",
        f"| Paired repetitions | {campaign.get('n_pairs', '—')} AB/BA pairs across "
        f"{campaign.get('n_cells', '—')} benchmark cells |",
        f"| Matrix complete | {'yes' if campaign.get('matrix_complete') else 'no'} |",
        f"| Publication mode | `{publication.get('mode', '—')}` |",
        f"| Bootstrap | seed {bootstrap.get('seed', '—')}; "
        f"{bootstrap.get('replicates', '—'):,} replicates |"
        if isinstance(bootstrap.get("replicates"), int)
        else f"| Bootstrap | seed {bootstrap.get('seed', '—')}; — replicates |",
        f"| Release verdict | **{publication_verdict}** |",
    ]
    evaluation = publication.get("evaluation_overlays")
    if isinstance(evaluation, Mapping) and evaluation.get("validated") is True:
        roots = evaluation.get("roots") or []
        evaluator_sha = (
            roots[0].get("evaluator_sha256")
            if roots and isinstance(roots[0], Mapping) else None
        )
        cutoff = max(
            (
                str(root.get("created_at"))
                for root in roots
                if isinstance(root, Mapping) and root.get("created_at")
            ),
            default="—",
        )
        rows.extend([
            f"| Corrected evaluation | immutable overlay; "
            f"{evaluation.get('pair_count', '—')} paired repetitions; "
            f"evaluator `{str(evaluator_sha)[:16] if evaluator_sha else '—'}…` |",
            f"| Evaluation cutoff | `{cutoff}` |",
        ])
    else:
        rows.append("| Corrected evaluation | unavailable; in-run scores only |")
    host = (metrics.get("environment") or {}).get("finalization_host")
    if isinstance(host, Mapping):
        memory = host.get("memory_gib")
        memory_text = f"{float(memory):.1f} GiB" if _number(memory) else "—"
        rows.append(
            f"| Finalization host *(not a historical run record)* | "
            f"`{host.get('hostname', '—')}`; {host.get('cpu_model', '—')}; "
            f"{host.get('sockets', '—')} sockets; "
            f"{host.get('logical_cpus', '—')} logical CPUs; {memory_text}; "
            f"kernel `{host.get('kernel', '—')}` |"
        )
    controller_evidence = publication.get("controller_evidence") or {}
    diagnostic_roots = controller_evidence.get("roots", [])
    if diagnostic_roots:
        environments = {
            json.dumps(
                (root.get("plan") or {}).get("execution_environment"),
                sort_keys=True,
            )
            for root in diagnostic_roots
            if isinstance(
                (root.get("plan") or {}).get("execution_environment"),
                Mapping,
            )
        }
        if len(environments) == 1:
            environment = json.loads(next(iter(environments)))
            tools = environment.get("tools") or {}
            packages = environment.get("packages") or {}
            rows.append(
                "| Frozen execution software | "
                f"Python {environment.get('python', '—')}; "
                f"minimap2 {tools.get('minimap2', '—')}; "
                f"miniprot {tools.get('miniprot', '—')}; "
                f"gffutils {packages.get('gffutils', '—')}; "
                f"parasail {packages.get('parasail', '—')} |"
            )
        elif environments:
            rows.append(
                f"| Frozen execution software | mixed across "
                f"{len(environments)} controller environments |"
            )
        rows.append(
            f"| Execution roots | {len(diagnostic_roots)} diagnostic controller roots |"
        )
        for root_record in diagnostic_roots:
            plan = root_record.get("plan") or {}
            policy = plan.get("policy") or {}
            if not plan:
                continue
            rows.append(
                f"| Policy `{plan.get('run_id', 'unknown')}` | active "
                f"{policy.get('max_active', '—')}; full "
                f"{policy.get('max_full', '—')}; worker threads "
                f"{policy.get('max_worker_threads', '—')}; states "
                f"`{plan.get('counts', {})}` |"
            )
    return "\n".join(rows)


def release_performance(metrics: Mapping[str, Any], _db: Mapping[str, Any], _vc: Mapping[str, Any],
                        _labels: Labels) -> str:
    rows = [
        "| Panel | Modes (candidate/reference) | Cells | Wall ratio GMR (95% CI when n≥2) | "
        "Total wall ratio | RSS ratio GMR (95% CI when n≥2) | Worst wall ratio | Worst RSS ratio | "
        "Candidate command-RSS slot-sum proxy |",
        "|---|---|---:|---:|---:|---:|---:|---:|---:|",
    ]
    panels = sorted(PANEL_ORDER, key=PANEL_ORDER.get)
    for panel in panels:
        summary = metrics["panels"].get(panel)
        if summary is None:
            rows.append(f"| {panel} | —/— | 0 | — | — | — | — | — | — |")
            continue
        modes = summary.get("modes") or {}
        rows.append(
            f"| {panel} | {modes.get('candidate', '—')}/{modes.get('reference', '—')} | "
            f"{summary.get('n_cells', 0)} | "
            f"{_ci(summary.get('wall_ratio'), summary.get('n_cells'))} | "
            f"{_ratio(summary.get('wall_total_ratio'))} | "
            f"{_ci(summary.get('rss_ratio'), summary.get('n_cells'))} | "
            f"{_ratio(summary.get('wall_ratio_worst'))} | "
            f"{_ratio(summary.get('rss_ratio_worst'))} | "
            f"{_fmt(summary.get('candidate_concurrent_memory_envelope_gib'), 2)} GiB "
            f"({summary.get('candidate_concurrency_limit', 0)} slots) |"
        )
    return "\n".join(rows)


def release_correctness(metrics: Mapping[str, Any], _db: Mapping[str, Any], _vc: Mapping[str, Any],
                        _labels: Labels) -> str:
    cells = metrics["cells"]
    total = len(cells)
    e2e_cells = [cell for cell in cells if cell.get("panel") == "e2e"]
    rows = [
        "| Check | Candidate | Reference/control |",
        "|---|---:|---:|",
        f"| GFF3 valid | {sum(cell.get('candidate_valid') is True for cell in cells)}/{total} | "
        f"{sum(cell.get('reference_valid') is True for cell in cells)}/{total} *(diagnostic)* |",
        f"| Byte deterministic | "
        f"{sum(cell.get('candidate_byte_deterministic') is True for cell in cells)}/{total} | "
        f"{sum(cell.get('reference_byte_deterministic') is True for cell in cells)}/{total} |",
        f"| Semantic deterministic | "
        f"{sum(cell.get('candidate_semantic_deterministic') is True for cell in cells)}/{total} | "
        f"{sum(cell.get('reference_semantic_deterministic') is True for cell in cells)}/{total} |",
        f"| Quality deterministic | "
        f"{sum(cell.get('candidate_quality_deterministic') is True for cell in cells)}/{total} | "
        f"{sum(cell.get('reference_quality_deterministic') is True for cell in cells)}/{total} |",
        f"| Common-set PI deterministic | "
        f"{sum(cell.get('common_pi_deterministic') is True for cell in cells)}/{total} | same paired evidence |",
    ]
    if e2e_cells:
        e2e_total = len(e2e_cells)
        for field, label in (
            ("valid", "E2E biological summary valid"),
            ("deterministic", "E2E biological summary deterministic"),
        ):
            rows.append(
                f"| {label} | "
                f"{sum(((cell.get('e2e_biology') or {}).get('candidate') or {}).get(field) is True for cell in e2e_cells)}/{e2e_total} | "
                f"{sum(((cell.get('e2e_biology') or {}).get('reference') or {}).get(field) is True for cell in e2e_cells)}/{e2e_total} |"
            )
    rows.extend([
        "",
        "| Panel | Quality metric | Candidate − reference (95% bootstrap CI when n≥2) |",
        "|---|---|---:|",
    ])
    panels = sorted(metrics["panels"].items(), key=lambda item: PANEL_ORDER.get(item[0], 99))
    for panel, summary in panels:
        for metric, interval in (summary.get("quality") or {}).items():
            if interval is not None:
                rows.append(
                    f"| {panel} | `{metric}` | "
                    f"{_delta_ci(interval, summary.get('n_cells'))} |"
                )
    rows.extend([
        "",
        "| Panel | Stable ID type | Applicable cells | Candidate preserved | "
        "Reference preserved | Rate delta |",
        "|---|---|---:|---:|---:|---:|",
    ])
    for panel, summary in panels:
        for feature_type, record in (summary.get("stable_id_preservation") or {}).items():
            rows.append(
                f"| {panel} | `{feature_type}` | {record.get('n_applicable_cells', 0)} | "
                f"{_preserved(record.get('candidate_n_preserved_ids'), record.get('n_reference_ids'), record.get('candidate_preservation_rate'))} | "
                f"{_preserved(record.get('reference_n_preserved_ids'), record.get('n_reference_ids'), record.get('reference_preservation_rate'))} | "
                f"{_signed(record.get('rate_delta'))} |"
            )
    return "\n".join(rows)


def release_coverage(metrics: Mapping[str, Any], audit: Mapping[str, Any], _labels: Labels) -> str:
    planned = metrics["campaign"]["planned_campaign"]
    panels = planned["panels"]
    matrix = {
        row.get("panel"): row
        for row in planned.get("matrix", [])
        if isinstance(row, Mapping)
    }
    failed_by_panel = Counter(cell["panel"] for cell in audit["cells"])
    rows = [
        "| Panel | Endpoint modes | Plan thread budget | Planned cases | Observed cases | "
        "Planned AB/BA pairs | Observed pairs | Missing v1.0.8 controls |",
        "|---|---|---:|---:|---:|---:|---:|---:|",
    ]
    totals = Counter()
    for panel in sorted(PANEL_ORDER, key=PANEL_ORDER.get):
        row = panels[panel]
        protocol = matrix.get(panel) or {}
        planned_cases = len(row["ids"])
        observed_cells = [cell for cell in metrics["cells"] if cell["panel"] == panel]
        observed_cases = len(observed_cells)
        planned_pairs = planned_cases * int(row["repetitions"])
        observed_pairs = sum(int(cell["repetitions"]) for cell in observed_cells)
        missing = failed_by_panel[panel]
        totals.update({
            "planned_cases": planned_cases,
            "observed_cases": observed_cases,
            "planned_pairs": planned_pairs,
            "observed_pairs": observed_pairs,
            "missing": missing,
        })
        rows.append(
            f"| {panel} | {protocol.get('candidate_mode', '—')}/"
            f"{protocol.get('reference_mode', '—')} | {protocol.get('threads', '—')} | "
            f"{planned_cases} | {observed_cases} | {planned_pairs} | "
            f"{observed_pairs} | {missing} |"
        )
    rows.extend([
        f"| **Total** | — | — | **{totals['planned_cases']}** | "
        f"**{totals['observed_cases']}** | **{totals['planned_pairs']}** | "
        f"**{totals['observed_pairs']}** | **{totals['missing']}** |",
        "",
        "| Evidence-integrity status | Result |",
        "|---|---:|",
        f"| Candidate GFF3 validation | "
        f"{sum(cell.get('candidate_valid') is True for cell in metrics['cells'])}/"
        f"{len(metrics['cells'])} cells |",
        f"| Candidate output + score determinism | "
        f"{sum(_candidate_deterministic(cell) for cell in metrics['cells'])}/"
        f"{len(metrics['cells'])} cells |",
        f"| Deterministic historical-reference failures | "
        f"{audit['failed_pair_count']} pairs across {len(audit['failed_case_counts'])} cases |",
        f"| Candidate-first outputs retained but excluded from pairing | "
        f"{audit['invariants']['candidate_first_outputs_preserved']} |",
    ])
    return "\n".join(rows)


def _candidate_deterministic(cell: Mapping[str, Any]) -> bool:
    return bool(
        cell.get("candidate_byte_deterministic") is True
        and cell.get("candidate_semantic_deterministic") is True
        and cell.get("candidate_quality_deterministic") is True
        and cell.get("common_pi_deterministic") is True
        and (
            cell.get("panel") != "e2e"
            or ((cell.get("e2e_biology") or {}).get("candidate") or {}).get(
                "deterministic"
            ) is True
        )
    )


def release_comprehensiveness(metrics: Mapping[str, Any], _audit: Mapping[str, Any],
                              _labels: Labels) -> str:
    rows = [
        "| Panel | Cells | Δ coding completeness (95% CI) | "
        "Δ all-feature completeness (95% CI) | Candidate/reference CDS IDs | "
        "Candidate/reference exon IDs |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for panel, summary in sorted(
        metrics["panels"].items(), key=lambda item: PANEL_ORDER.get(item[0], 99)
    ):
        quality = summary.get("quality") or {}
        stable = summary.get("stable_id_preservation") or {}
        cds = stable.get("CDS") or {}
        exon = stable.get("exon") or {}
        rows.append(
            f"| {panel} | {summary.get('n_cells', 0)} | "
            f"{_delta_ci(quality.get('completeness_coding'), summary.get('n_cells'))} | "
            f"{_delta_ci(quality.get('completeness_feature_total'), summary.get('n_cells'))} | "
            f"{_preserved(cds.get('candidate_n_preserved_ids'), cds.get('n_reference_ids'), cds.get('candidate_preservation_rate'))} / "
            f"{_preserved(cds.get('reference_n_preserved_ids'), cds.get('n_reference_ids'), cds.get('reference_preservation_rate'))} | "
            f"{_preserved(exon.get('candidate_n_preserved_ids'), exon.get('n_reference_ids'), exon.get('candidate_preservation_rate'))} / "
            f"{_preserved(exon.get('reference_n_preserved_ids'), exon.get('n_reference_ids'), exon.get('reference_preservation_rate'))} |"
        )
    return "\n".join(rows)


def _arm_quality(cell: Mapping[str, Any], arm: str) -> Mapping[str, Any]:
    absolute = cell.get("absolute_quality") or {}
    value = absolute.get(arm) if isinstance(absolute, Mapping) else None
    detailed = cell.get(arm)
    merged = dict(detailed) if isinstance(detailed, Mapping) else {}
    if isinstance(value, Mapping):
        merged.update(value)
    return merged


def _paired_metric(candidate: Any, reference: Any, digits: int = 5) -> str:
    return f"{_fmt(candidate, digits)} / {_fmt(reference, digits)}"


def release_e2e_outcomes(metrics: Mapping[str, Any], _audit: Mapping[str, Any],
                         _labels: Labels) -> str:
    rows = [
        "| E2E genome | Coding recovered candidate/reference (denominator) | "
        "Coding completeness | All-feature completeness | Mean PI | CovPI | "
        "Exact intron chain | ORF valid |",
        "|---|---:|---:|---:|---:|---:|---:|---:|",
    ]
    e2e_cells = sorted(
        (cell for cell in metrics["cells"] if cell.get("panel") == "e2e"),
        key=lambda cell: cell.get("benchmark", ""),
    )
    for cell in e2e_cells:
        candidate = _arm_quality(cell, "candidate")
        reference = _arm_quality(cell, "reference")
        denominator = candidate.get("n_reference_coding")
        rows.append(
            f"| {PRETTY.get(cell.get('benchmark'), str(cell.get('benchmark')).title())} | "
            f"{_fmt(candidate.get('n_recovered_coding'))} / "
            f"{_fmt(reference.get('n_recovered_coding'))} ({_fmt(denominator)}) | "
            f"{_paired_metric(candidate.get('completeness_coding'), reference.get('completeness_coding'))} | "
            f"{_paired_metric(candidate.get('completeness_feature_total'), reference.get('completeness_feature_total'))} | "
            f"{_paired_metric(candidate.get('mean_pi'), reference.get('mean_pi'))} | "
            f"{_paired_metric(candidate.get('covpi'), reference.get('covpi'), 6)} | "
            f"{_paired_metric(candidate.get('intron_chain_exact_recall'), reference.get('intron_chain_exact_recall'))} | "
            f"{_paired_metric(candidate.get('orf_valid_recall'), reference.get('orf_valid_recall'))} |"
        )
    if not e2e_cells:
        raise ValueError("release E2E outcomes require observed end-to-end cells")
    return "\n".join(rows)


def _gate_value(value: Any) -> str:
    if _number(value):
        return f"{float(value):.6g}"
    if isinstance(value, str):
        return f"`{value}`"
    return f"`{json.dumps(value, sort_keys=True)}`"


def release_failed_gates(metrics: Mapping[str, Any], _audit: Mapping[str, Any],
                         _labels: Labels) -> str:
    checks = [
        check
        for check in (metrics.get("verdict") or {}).get("checks", [])
        if isinstance(check, Mapping) and check.get("passed") is False
    ]
    if not checks:
        return "| Failed gate | Observed | Required |\n|---|---:|---:|\n| None | — | — |"
    rows = [
        "| Failed gate | Observed | Required | Interpretation |",
        "|---|---:|---:|---|",
    ]
    for check in checks:
        name = str(check.get("name", "unnamed"))
        if name == "complete_immutable_release_campaign":
            interpretation = "Incomplete controller-backed release matrix"
        elif name == "all_candidate_outputs_valid":
            interpretation = "Exact-tag structural validity defect"
        elif name.endswith(("worst_wall_cell", "worst_rss_cell")):
            interpretation = "Worst-cell performance cap exceeded"
        elif name.endswith(".lower"):
            interpretation = "Bootstrap lower confidence bound below floor"
        else:
            interpretation = "Prespecified qualification check failed"
        rows.append(
            f"| `{name}` | {_gate_value(check.get('actual'))} | "
            f"{_gate_value(check.get('limit'))} | {interpretation} |"
        )
    return "\n".join(rows)


def release_failure_audit(metrics: Mapping[str, Any], audit: Mapping[str, Any],
                          _labels: Labels) -> str:
    missing_by_case = Counter(
        (cell["panel"], cell["benchmark"])
        for cell in audit["cells"]
    )
    candidate_first = Counter(
        (cell["panel"], cell["benchmark"])
        for cell in audit["cells"]
        if cell["candidate_completed_before_control_failure"]
    )
    elapsed = {
        (cell["panel"], cell["benchmark"]): []
        for cell in audit["cells"]
    }
    for cell in audit["cells"]:
        elapsed[(cell["panel"], cell["benchmark"])].append(float(cell["elapsed_seconds"]))
    rows = [
        "| Panel | Missing case | Failed control repetitions | Candidate-first outputs retained | "
        "Control elapsed range (s) | Classification |",
        "|---|---|---:|---:|---:|---|",
    ]
    for panel, benchmark in sorted(missing_by_case, key=lambda key: (PANEL_ORDER.get(key[0], 99), key[1])):
        values = elapsed[(panel, benchmark)]
        rows.append(
            f"| {panel} | {PRETTY.get(benchmark, benchmark)} | "
            f"{missing_by_case[(panel, benchmark)]} | "
            f"{candidate_first[(panel, benchmark)]} | "
            f"{min(values):.1f}–{max(values):.1f} | "
            "v1.0.8 `FeatureNotFoundError` + exception-string `TypeError` |"
        )
    observed = _observed_pair_keys(metrics)
    rows.append(
        f"| **Total** | **7 cases** | **{audit['failed_pair_count']}** | "
        f"**{audit['invariants']['candidate_first_outputs_preserved']}** | — | "
        f"{audit['classification']}; excluded from {len(observed)} observed pairs |"
    )
    return "\n".join(rows)


def release_cells(metrics: Mapping[str, Any], _db: Mapping[str, Any], _vc: Mapping[str, Any],
                  _labels: Labels) -> str:
    rows = [
        "| Panel | Benchmark | Reps | Candidate wall (s) | Reference wall (s) | "
        "Wall factor (reference/candidate) | Candidate peak RSS | "
        "RSS factor (reference/candidate) | Δ coding | Δ all-feature | Δ CovPI | "
        "Δ exact intron | Δ ORF valid | Δ common-set PI | Candidate valid |",
        "|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    cells = sorted(
        metrics["cells"],
        key=lambda cell: (PANEL_ORDER.get(cell.get("panel"), 99), cell.get("benchmark", "")),
    )
    for cell in cells:
        common = cell.get("common_pi") or {}
        quality = cell.get("quality_deltas") or {}
        rows.append(
            f"| {cell.get('panel', '—')} | {PRETTY.get(cell.get('benchmark'), cell.get('benchmark', '—'))} | "
            f"{cell.get('repetitions', 0)} | "
            f"{_fmt(cell.get('candidate_wall_seconds_median'), 2)} | "
            f"{_fmt(cell.get('reference_wall_seconds_median'), 2)} | "
            f"{_inverse_ratio(cell.get('wall_ratio_median'))} | "
            f"{_fmt(cell.get('candidate_peak_rss_gib'), 2)} GiB | "
            f"{_inverse_ratio(cell.get('rss_ratio_median'))} | "
            f"{_signed(quality.get('completeness_coding'))} | "
            f"{_signed(quality.get('completeness_feature_total'))} | "
            f"{_signed(quality.get('covpi'))} | "
            f"{_signed(quality.get('intron_chain_exact_recall'))} | "
            f"{_signed(quality.get('orf_valid_recall'))} | "
            f"{_signed(common.get('mean_delta'))} | "
            f"{'yes' if cell.get('candidate_valid') is True else 'no'} |"
        )
    return "\n".join(rows)


BlockRenderer = Callable[
    [Mapping[str, Any], Mapping[str, Any], Mapping[str, Any], Labels],
    str,
]

BLOCKS: dict[str, BlockRenderer] = {
    "release-provenance": release_provenance,
    "release-performance": release_performance,
    "release-correctness": release_correctness,
    "release-cells": release_cells,
    "legacy-recovery": lambda _metrics, db, vc, labels: table3(db, vc, labels),
    "legacy-performance": lambda _metrics, db, vc, labels: table4(db, vc, labels),
    "legacy-accuracy": lambda _metrics, db, vc, labels: table5(db, vc, labels),
    "legacy-completeness": lambda _metrics, db, vc, labels: table7(db, vc, labels),
}

AuditBlockRenderer = Callable[
    [Mapping[str, Any], Mapping[str, Any], Labels],
    str,
]

AUDIT_BLOCKS: dict[str, AuditBlockRenderer] = {
    "release-coverage": release_coverage,
    "release-comprehensiveness": release_comprehensiveness,
    "release-e2e-outcomes": release_e2e_outcomes,
    "release-failed-gates": release_failed_gates,
    "release-failure-audit": release_failure_audit,
}
BLOCK_NAMES = set(BLOCKS) | set(AUDIT_BLOCKS)


def render_block(name: str, *, metrics: Mapping[str, Any] | None,
                 fourway: Mapping[str, Any] | None, version_compare: Mapping[str, Any] | None,
                 labels: Labels, failure_audit: Mapping[str, Any] | None = None,
                 manifest: Mapping[str, Any] | None = None) -> str:
    if name not in BLOCK_NAMES:
        raise ValueError(f"unknown generated block {name!r}")
    if name.startswith("release-"):
        source_metrics = _validate_release_metrics(metrics)
    else:
        source_metrics = metrics or {}
    if name.startswith("legacy-") and (fourway is None or version_compare is None):
        raise ValueError(f"generated block {name!r} requires legacy result documents")
    if failure_audit is not None:
        validated_audit = validate_failure_audit(failure_audit, source_metrics)
    else:
        validated_audit = None
    if manifest is not None:
        validate_release_manifest(manifest, source_metrics)
    if name in AUDIT_BLOCKS:
        if validated_audit is None:
            raise ValueError(f"generated block {name!r} requires --failure-audit")
        return AUDIT_BLOCKS[name](source_metrics, validated_audit, labels).rstrip()
    return BLOCKS[name](
        source_metrics,
        fourway or {},
        version_compare or {},
        labels,
    ).rstrip()


def wrapped_block(name: str, content: str, *, mdx: bool = False) -> str:
    start_marker = MDX_START_MARKER if mdx else START_MARKER
    end_marker = MDX_END_MARKER if mdx else END_MARKER
    return "\n".join([
        start_marker.format(name=name),
        content.rstrip(),
        end_marker.format(name=name),
    ])


def update_generated_blocks(report: str, *, metrics: Mapping[str, Any] | None,
                            fourway: Mapping[str, Any] | None,
                            version_compare: Mapping[str, Any] | None,
                            labels: Labels,
                            failure_audit: Mapping[str, Any] | None = None,
                            manifest: Mapping[str, Any] | None = None) -> str:
    matches = list(MARKER_RE.finditer(report))
    names = [match.group("html_name") or match.group("mdx_name") for match in matches]
    if not names:
        raise ValueError("report contains no generated table markers")
    duplicates = sorted({name for name in names if names.count(name) > 1})
    if duplicates:
        raise ValueError(f"report contains duplicate generated blocks: {duplicates}")
    updated = report
    for match, name in zip(matches, names):
        mdx = match.group("mdx_name") is not None
        start_marker = MDX_START_MARKER if mdx else START_MARKER
        end_marker = MDX_END_MARKER if mdx else END_MARKER
        start = start_marker.format(name=name)
        end = end_marker.format(name=name)
        start_index = updated.find(start)
        end_index = updated.find(end, start_index + len(start))
        if end_index < 0:
            raise ValueError(f"generated block {name!r} has no end marker")
        rendered = wrapped_block(
            name,
            render_block(
                name,
                metrics=metrics,
                fourway=fourway,
                version_compare=version_compare,
                labels=labels,
                failure_audit=failure_audit,
                manifest=manifest,
            ),
            mdx=mdx,
        )
        updated = updated[:start_index] + rendered + updated[end_index + len(end):]
    return updated


def update_report_placeholders(
    report: str,
    metrics: Mapping[str, Any],
) -> str:
    """Replace narrative tokens with values derived from release metrics."""

    campaign = metrics.get("campaign") or {}
    publication = metrics.get("publication") or {}
    cells = list(metrics.get("cells") or [])
    planned = campaign.get("planned_campaign") or {}
    planned_count = sum(
        len(row.get("ids", [])) * int(row.get("repetitions", 0))
        for row in (planned.get("panels") or {}).values()
    )
    observed_count = int(campaign.get("n_cells", len(cells)))
    observed_pairs = int(campaign.get("n_pairs", 0))
    panel_counts = Counter(cell.get("panel", "unknown") for cell in cells)
    panel_text = ", ".join(
        f"{panel}={panel_counts.get(panel, 0)}"
        for panel in sorted(PANEL_ORDER, key=PANEL_ORDER.get)
    )
    missing_keys = campaign.get("missing_keys") or []
    missing_count = len(missing_keys)
    missing_by_panel = Counter(
        key[0] for key in missing_keys if isinstance(key, (list, tuple)) and key
    )
    missing_panel_text = ", ".join(
        f"{panel}={count}" for panel, count in sorted(missing_by_panel.items())
    ) or "none"
    outlier_count = len(metrics.get("performance_outliers") or [])
    candidate_valid = sum(cell.get("candidate_valid") is True for cell in cells)
    candidate_deterministic = sum(
        cell.get("candidate_byte_deterministic") is True
        and cell.get("candidate_semantic_deterministic") is True
        and cell.get("candidate_quality_deterministic") is True
        and cell.get("common_pi_deterministic") is True
        and (
            cell.get("panel") != "e2e"
            or ((cell.get("e2e_biology") or {}).get("candidate") or {}).get(
                "deterministic"
            ) is True
        )
        for cell in cells
    )
    mode = publication.get("mode", "diagnostic")
    verdict = metrics.get("verdict") or {}
    failed_checks = [
        check.get("name", "unnamed check")
        for check in verdict.get("checks", [])
        if isinstance(check, Mapping)
        and check.get("passed") is False
        and check.get("name") != "complete_immutable_release_campaign"
    ]
    performance_suffixes = (
        ".wall_gmr_upper",
        ".rss_gmr_upper",
        ".wall_total_ratio",
        ".candidate_concurrent_memory_envelope_gib",
        ".worst_wall_cell",
        ".worst_rss_cell",
    )
    failed_performance = [
        name for name in failed_checks
        if name.endswith(performance_suffixes)
    ]
    failed_correctness = [
        name for name in failed_checks
        if name not in failed_performance
    ]
    verdict_label = (
        "PASS" if mode == "release" and verdict.get("passed") is True
        else "FAIL" if mode == "release"
        else "DIAGNOSTIC ONLY"
    )
    blocking_text = (
        f" {len(failed_checks)} additional blocking check(s) fail in the "
        "observed matrix."
        if failed_checks else ""
    )
    panel_labels = {
        "subset": "subset",
        "full": "full genome",
        "e2e": "end to end",
    }
    abstract_ratios = []
    performance_details = []
    quality_details = []
    for panel in sorted(PANEL_ORDER, key=PANEL_ORDER.get):
        summary = (metrics.get("panels") or {}).get(panel)
        if not isinstance(summary, Mapping):
            continue
        n_cells = int(summary.get("n_cells", 0))
        wall = summary.get("wall_ratio") or {}
        rss = summary.get("rss_ratio") or {}
        wall_estimate = wall.get("estimate")
        rss_estimate = rss.get("estimate")
        if _number(wall_estimate) and _number(rss_estimate):
            abstract_ratios.append(
                f"{panel_labels[panel]} {float(wall_estimate):.3f}× wall/"
                f"{float(rss_estimate):.3f}× RSS"
            )
            if n_cells >= 2 and all(
                _number(interval.get(bound))
                for interval in (wall, rss)
                for bound in ("low", "high")
            ):
                interval_text = (
                    f"wall {float(wall_estimate):.3f}× "
                    f"[95% CI {float(wall['low']):.3f}, "
                    f"{float(wall['high']):.3f}], RSS "
                    f"{float(rss_estimate):.3f}× "
                    f"[95% CI {float(rss['low']):.3f}, "
                    f"{float(rss['high']):.3f}]"
                )
            else:
                interval_text = (
                    f"wall {float(wall_estimate):.3f}× and RSS "
                    f"{float(rss_estimate):.3f}× (descriptive)"
                )
            performance_details.append(
                f"{panel_labels[panel]}: {interval_text}, total wall "
                f"{_ratio(summary.get('wall_total_ratio'))}, worst wall "
                f"{_ratio(summary.get('wall_ratio_worst'))}, and worst RSS "
                f"{_ratio(summary.get('rss_ratio_worst'))}"
            )
        quality = summary.get("quality") or {}
        quality_values = []
        for metric, label in (
            ("completeness_coding", "coding completeness"),
            ("completeness_feature_total", "feature completeness"),
            ("covpi", "CovPI"),
        ):
            interval = quality.get(metric) or {}
            if _number(interval.get("estimate")):
                quality_values.append(
                    f"Δ {label} {_signed(interval['estimate'])}"
                )
        if quality_values:
            quality_details.append(
                f"{panel_labels[panel]} ({', '.join(quality_values)})"
            )
    ratio_text = "; ".join(abstract_ratios) or "unavailable"
    performance_text = "; ".join(performance_details) or "No panel estimates are available"
    quality_text = "; ".join(quality_details) or "no finite panel deltas"
    failed_performance_text = (
        " Named performance gates failing: "
        + ", ".join(f"`{name}`" for name in failed_performance)
        + "."
        if failed_performance
        else " All named performance gates pass in the observed matrix."
    )
    failed_correctness_text = (
        f" {len(failed_correctness)} named correctness/control gate(s) fail "
        "in the observed matrix."
        if failed_correctness
        else " All named correctness/control gates pass in the observed matrix."
    )
    replacements = {
        "@@FINAL_ABSTRACT_RESULTS@@": (
            f"The selected immutable roots contain {observed_pairs} of "
            f"{planned_count} planned paired executions across "
            f"{observed_count} observed benchmark cells ({panel_text}); "
            f"{missing_count} planned executions are unavailable. "
            "Candidate/reference geometric-mean ratios are "
            f"{ratio_text}; values below 1 favor v1.0.11. Candidate outputs "
            f"pass current validation in {candidate_valid}/{observed_count} "
            "cells, and output plus score evidence is deterministic in "
            f"{candidate_deterministic}/{observed_count}. The publication "
            f"verdict is {verdict_label}."
        ),
        "@@FINAL_PERFORMANCE_INTERPRETATION@@": (
            f"The paired AB/BA panel estimates are {performance_text}. Ratios "
            "below 1 favor v1.0.11."
            f"{failed_performance_text} {outlier_count} observed cell(s) meet "
            "the variability or ratio screening criteria; a screening flag is "
            "not itself a release gate. Execution-policy epochs remain explicit, "
            "and isolated reruns are diagnostic rather than replacements for "
            "canonical measurements."
        ),
        "@@FINAL_CORRECTNESS_INTERPRETATION@@": (
            f"Corrected candidate-minus-reference panel estimates are {quality_text}. "
            f"Among {observed_count} observed cells, candidate outputs pass "
            f"current validation in {candidate_valid} and output plus score "
            f"evidence is deterministic in {candidate_deterministic} ({panel_text})."
            f"{failed_correctness_text} Missing control pairs are not converted "
            "into scores."
        ),
        "@@FINAL_RELEASE_VERDICT@@": (
            f"The exact-tag result is **{verdict_label}**: "
            f"{observed_count} cells were observed and {missing_count} planned "
            f"paired executions are unavailable in the selected immutable roots "
            f"({missing_panel_text}). A release PASS is not emitted "
            "because the immutable campaign evidence is not a complete passing "
            f"release publication.{blocking_text}"
        ),
        "@@FINAL_CONCLUSION@@": (
            "The current evidence distinguishes exact-release behavior, historical "
            "comparisons, and contention-aware throughput. All reported values "
            "are generated from the immutable result documents, with missing or "
            "failed control cells retained as explicit limitations."
        ),
    }
    updated = report
    for token, replacement in replacements.items():
        updated = updated.replace(token, replacement)
    return updated


def check_generated_blocks(path: Path, *, metrics: Mapping[str, Any] | None,
                           fourway: Mapping[str, Any] | None,
                           version_compare: Mapping[str, Any] | None,
                           labels: Labels,
                           failure_audit: Mapping[str, Any] | None = None,
                           manifest: Mapping[str, Any] | None = None) -> tuple[bool, str]:
    observed = path.read_text()
    expected = update_generated_blocks(
        observed,
        metrics=metrics,
        fourway=fourway,
        version_compare=version_compare,
        labels=labels,
        failure_audit=failure_audit,
        manifest=manifest,
    )
    if metrics is not None:
        expected = update_report_placeholders(expected, metrics)
    if observed == expected:
        return True, ""
    diff = "".join(difflib.unified_diff(
        observed.splitlines(keepends=True),
        expected.splitlines(keepends=True),
        fromfile=str(path),
        tofile=f"{path} (generated)",
    ))
    return False, diff


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    operation = parser.add_mutually_exclusive_group()
    operation.add_argument("--table", choices=sorted(TABLES, key=int))
    operation.add_argument("--all", action="store_true", help="render every legacy table")
    operation.add_argument("--block", choices=sorted(BLOCK_NAMES))
    operation.add_argument("--blocks", action="store_true", help="render every generated block")
    operation.add_argument("--list", action="store_true", help="list tables and blocks")
    operation.add_argument("--check", type=Path, metavar="REPORT")
    operation.add_argument("--update", type=Path, metavar="REPORT")
    parser.add_argument("--release-metrics", type=Path)
    parser.add_argument(
        "--whole-genome-metrics", type=Path,
        help="seven-transfer biology-first metrics",
    )
    parser.add_argument(
        "--annotation-sensitivity", type=Path,
        help="released target-annotation sensitivity evidence",
    )
    parser.add_argument(
        "--recovery-difference", type=Path,
        help="paired coding-transcript recovery-difference evidence",
    )
    parser.add_argument("--release-manifest", type=Path)
    parser.add_argument("--failure-audit", type=Path)
    parser.add_argument("--fourway", type=Path, default=FOURWAY)
    parser.add_argument("--version-compare", type=Path, default=VERSION_CMP)
    parser.add_argument("--candidate-label", default=Labels.candidate)
    parser.add_argument("--reference-label", default=Labels.reference)
    parser.add_argument("--legacy-candidate-label", default=Labels.legacy_candidate)
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    labels = Labels(
        candidate=args.candidate_label,
        reference=args.reference_label,
        legacy_candidate=args.legacy_candidate_label,
    )
    if args.list:
        for number, (title, _renderer) in sorted(TABLES.items(), key=lambda item: int(item[0])):
            print(f"Table {number}: {title}")
        for name in sorted(BLOCK_NAMES):
            print(f"Block {name}")
        return 0

    if args.whole_genome_metrics:
        if not (args.check or args.update):
            parser.error(
                "--whole-genome-metrics requires --check or --update"
            )
        if args.release_metrics or args.release_manifest or args.failure_audit:
            parser.error(
                "whole-genome metrics cannot be combined with the earlier "
                "release bundle"
            )
        from . import whole_genome_assets

        try:
            study_metrics = whole_genome_assets.load_metrics(
                args.whole_genome_metrics
            )
            sensitivity = (
                whole_genome_assets.load_sensitivity(
                    args.annotation_sensitivity, study_metrics,
                )
                if args.annotation_sensitivity else None
            )
            recovery = (
                whole_genome_assets.load_recovery_difference(
                    args.recovery_difference, study_metrics,
                )
                if args.recovery_difference else None
            )
            report_path = args.check or args.update
            passed, diff = whole_genome_assets.check_report(
                report_path, study_metrics, sensitivity, recovery,
            )
            if args.update:
                if not passed:
                    report_path.write_text(
                        whole_genome_assets.update_report(
                            report_path.read_text(encoding="utf-8"),
                            study_metrics,
                            sensitivity,
                            recovery,
                        ),
                        encoding="utf-8",
                    )
                    print(f"updated {report_path}")
                else:
                    print(f"up to date: {report_path}")
                return 0
            if passed:
                print(f"generated tables match: {report_path}")
                return 0
            sys.stderr.write(diff)
            return 1
        except (OSError, TypeError, ValueError, RuntimeError) as exc:
            parser.error(str(exc))

    if args.recovery_difference:
        parser.error("--recovery-difference requires --whole-genome-metrics")

    needs_legacy = bool(
        args.table
        or args.all
        or args.blocks
        or (args.block and args.block.startswith("legacy-"))
        or args.check
        or args.update
    )
    needs_metrics = bool(
        args.blocks
        or (args.block and args.block.startswith("release-"))
        or args.check
        or args.update
    )
    try:
        fourway = load(args.fourway) if needs_legacy else None
        version_compare = load(args.version_compare) if needs_legacy else None
        metrics = load(args.release_metrics) if args.release_metrics else None
        manifest = load(args.release_manifest) if args.release_manifest else None
        failure_audit = load(args.failure_audit) if args.failure_audit else None
        if needs_metrics and metrics is None:
            raise ValueError("release table operation requires --release-metrics")

        if args.table:
            title, renderer = TABLES[args.table]
            print(f"**Table {args.table}.** {title}\n")
            print(renderer(fourway, version_compare, labels))
            return 0
        if args.all:
            for number, (title, renderer) in sorted(TABLES.items(), key=lambda item: int(item[0])):
                print(f"\n**Table {number}.** {title}\n")
                print(renderer(fourway, version_compare, labels))
            return 0
        if args.block:
            print(wrapped_block(
                args.block,
                render_block(
                    args.block,
                    metrics=metrics,
                    fourway=fourway,
                    version_compare=version_compare,
                    labels=labels,
                    failure_audit=failure_audit,
                    manifest=manifest,
                ),
            ))
            return 0
        if args.blocks:
            names = sorted(BLOCK_NAMES if failure_audit is not None else BLOCKS)
            for name in names:
                print(wrapped_block(
                    name,
                    render_block(
                        name,
                        metrics=metrics,
                        fourway=fourway,
                        version_compare=version_compare,
                        labels=labels,
                        failure_audit=failure_audit,
                        manifest=manifest,
                    ),
                ))
                print()
            return 0
        if args.check or args.update:
            report_path = args.check or args.update
            passed, diff = check_generated_blocks(
                report_path,
                metrics=metrics,
                fourway=fourway,
                version_compare=version_compare,
                labels=labels,
                failure_audit=failure_audit,
                manifest=manifest,
            )
            if args.update:
                if not passed:
                    updated = update_generated_blocks(
                        report_path.read_text(),
                        metrics=metrics,
                        fourway=fourway,
                        version_compare=version_compare,
                        labels=labels,
                        failure_audit=failure_audit,
                        manifest=manifest,
                    )
                    report_path.write_text(update_report_placeholders(updated, metrics))
                    print(f"updated {report_path}")
                else:
                    print(f"up to date: {report_path}")
                return 0
            if passed:
                print(f"generated tables match: {report_path}")
                return 0
            sys.stderr.write(diff)
            return 1
    except (OSError, TypeError, ValueError) as exc:
        parser.error(str(exc))

    parser.print_help()
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
