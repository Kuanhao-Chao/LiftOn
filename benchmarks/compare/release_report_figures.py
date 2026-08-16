#!/usr/bin/env python
"""Generate evidence-bound figures for the LiftOn v1.0.11 report.

The historical two-figure interface remains available when only ``--metrics``
is supplied.  Passing ``--failure-audit`` selects the exact-release report
inventory and verifies the complete metrics/audit correspondence before any
figure is written.  ``--manifest`` and ``--legacy-results`` add bundle-level
cross-checks and the explicitly separated historical cross-tool appendix.
"""
from __future__ import annotations

import argparse
from collections import Counter
import json
import math
from pathlib import Path
from typing import Any, Mapping, Sequence

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

CANDIDATE_SHA = "c623f0bddc5b5051a3670a3b0064c52cc61bb719"
REFERENCE_SHA = "e503643d8346c600fedabcd3a4dff5c0873a4a37"
LEGACY_FUNCTIONAL_SHA = "7fc84192d85044c946eee5b90e0683eed722f9ea"

PANEL_ORDER = {"subset": 0, "full": 1, "e2e": 2}
PANEL_LABELS = {
    "subset": "Controlled subset",
    "full": "Controlled whole genome",
    "e2e": "End to end",
}
PANEL_COLORS = {
    "subset": "#4c78a8",
    "full": "#7a5195",
    "e2e": "#e45756",
}
TOOL_COLORS = {
    "reference": "#9d755d",
    "candidate": "#54a24b",
    "liftoff": "#4c78a8",
    "miniprot": "#f58518",
}
QUALITY_LABELS = {
    "completeness_coding": "Coding completeness",
    "completeness_feature_total": "All-feature completeness",
    "covpi": "CovPI",
    "recall_at_0.5": "Recall at PI ≥ 0.50",
    "recall_at_0.75": "Recall at PI ≥ 0.75",
    "recall_at_0.9": "Recall at PI ≥ 0.90",
    "recall_at_0.95": "Recall at PI ≥ 0.95",
    "intron_chain_exact_recall": "Exact intron-chain recall",
    "orf_valid_recall": "ORF-valid recall",
}
PRIMARY_QUALITY = (
    "completeness_coding",
    "completeness_feature_total",
    "covpi",
    "intron_chain_exact_recall",
    "orf_valid_recall",
)
ACCURACY_QUALITY = (
    "covpi",
    "recall_at_0.5",
    "recall_at_0.75",
    "recall_at_0.9",
    "recall_at_0.95",
    "intron_chain_exact_recall",
    "orf_valid_recall",
)

STUDY_FIGURE = "rfig_release_study_design.png"
COMPREHENSIVENESS_FIGURE = "rfig_release_comprehensiveness.png"
ACCURACY_FIGURE = "rfig_release_accuracy.png"
PERFORMANCE_FIGURE = "rfig_release_performance.png"
CORRECTNESS_FIGURE = "rfig_release_correctness.png"
HISTORICAL_FIGURE = "rfig_historical_cross_tool.png"
EXACT_FIGURE_INVENTORY = {
    STUDY_FIGURE,
    COMPREHENSIVENESS_FIGURE,
    ACCURACY_FIGURE,
    PERFORMANCE_FIGURE,
}


def _figure_context(metrics: Mapping[str, Any]) -> tuple[str, str]:
    publication = metrics.get("publication")
    if isinstance(publication, Mapping) and publication.get("mode") == "release":
        return "exact-release qualification", "PASS" if metrics["verdict"]["passed"] else "FAIL"
    return "paired benchmark diagnostic", "DIAGNOSTIC"


def _number(value: Any) -> bool:
    return (
        isinstance(value, (int, float))
        and not isinstance(value, bool)
        and math.isfinite(float(value))
    )


def _load_json(path: Path, label: str) -> Any:
    if not path.is_file():
        raise ValueError(f"missing {label}: {path}")
    try:
        return json.loads(path.read_text())
    except json.JSONDecodeError as exc:
        raise ValueError(f"invalid {label} JSON: {exc}") from exc


def load_metrics(path: Path) -> Mapping[str, Any]:
    return validate_metrics(_load_json(path, "release metrics"))


def load_failure_audit(path: Path) -> Mapping[str, Any]:
    document = _load_json(path, "failure audit")
    if not isinstance(document, Mapping):
        raise ValueError("failure audit must be a JSON object")
    return document


def load_manifest(path: Path) -> Mapping[str, Any]:
    document = _load_json(path, "release manifest")
    if not isinstance(document, Mapping):
        raise ValueError("release manifest must be a JSON object")
    return document


def validate_metrics(document: Any) -> Mapping[str, Any]:
    if not isinstance(document, Mapping) or document.get("schema_version") != 4:
        raise ValueError("release figures require schema-version 4 metrics")
    panels = document.get("panels")
    cells = document.get("cells")
    verdict = document.get("verdict")
    campaign = document.get("campaign")
    if not isinstance(campaign, Mapping):
        raise ValueError("release metrics contain no campaign")
    if not isinstance(panels, Mapping) or not panels:
        raise ValueError("release metrics contain no panels")
    if not isinstance(cells, list) or not cells:
        raise ValueError("release metrics contain no cells")
    if not isinstance(verdict, Mapping) or not isinstance(verdict.get("passed"), bool):
        raise ValueError("release metrics contain no boolean verdict")
    if campaign.get("n_cells") != len(cells):
        raise ValueError("release metrics campaign n_cells disagrees with cells")
    seen_cells = set()
    for cell in cells:
        if not isinstance(cell, Mapping):
            raise ValueError("release metrics contain a malformed cell")
        key = (cell.get("panel"), cell.get("benchmark"))
        if key in seen_cells:
            raise ValueError(f"release metrics contain duplicate cell {key!r}")
        seen_cells.add(key)
        if key[0] not in panels or not isinstance(key[1], str):
            raise ValueError(f"release metrics contain malformed cell key {key!r}")
        repetitions = cell.get("repetitions")
        if not isinstance(repetitions, int) or isinstance(repetitions, bool) or repetitions < 1:
            raise ValueError(f"release metrics cell {key!r} has invalid repetitions")
    for panel, summary in panels.items():
        if not isinstance(summary, Mapping):
            raise ValueError(f"panel {panel!r} summary is malformed")
        observed = sum(cell.get("panel") == panel for cell in cells)
        if summary.get("n_cells") != observed:
            raise ValueError(f"panel {panel!r} n_cells disagrees with cells")
        for name in ("wall_ratio", "rss_ratio"):
            interval = summary.get(name)
            if not isinstance(interval, Mapping) or not all(
                _number(interval.get(key)) for key in ("estimate", "low", "high")
            ):
                raise ValueError(f"panel {panel!r} has malformed {name}")
    return document


def _planned_pair_keys(metrics: Mapping[str, Any]) -> set[tuple[str, str, int]]:
    planned = (metrics["campaign"].get("planned_campaign") or {}).get("panels")
    if not isinstance(planned, Mapping) or not planned:
        raise ValueError("exact-release metrics contain no planned panel matrix")
    keys = set()
    for panel, row in planned.items():
        if panel not in PANEL_ORDER or not isinstance(row, Mapping):
            raise ValueError(f"planned panel {panel!r} is malformed")
        ids = row.get("ids")
        repetitions = row.get("repetitions")
        if (
            not isinstance(ids, list)
            or not ids
            or len(set(ids)) != len(ids)
            or not all(isinstance(item, str) and item for item in ids)
            or not isinstance(repetitions, int)
            or isinstance(repetitions, bool)
            or repetitions < 1
        ):
            raise ValueError(f"planned panel {panel!r} has invalid ids/repetitions")
        keys.update((panel, benchmark, repetition) for benchmark in ids for repetition in range(1, repetitions + 1))
    return keys


def _observed_pair_keys(metrics: Mapping[str, Any]) -> set[tuple[str, str, int]]:
    return {
        (cell["panel"], cell["benchmark"], repetition)
        for cell in metrics["cells"]
        for repetition in range(1, int(cell["repetitions"]) + 1)
    }


def validate_failure_audit(document: Any, metrics: Mapping[str, Any]) -> Mapping[str, Any]:
    metrics = validate_metrics(metrics)
    if metrics["campaign"].get("candidate_sha") != CANDIDATE_SHA:
        raise ValueError("failure-audit figures require exact v1.0.11 candidate SHA")
    if metrics["campaign"].get("reference_sha") != REFERENCE_SHA:
        raise ValueError("failure-audit figures require exact v1.0.8 reference SHA")
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
        raise ValueError("failure audit selected_roots are malformed")
    if not isinstance(cells, list) or not cells:
        raise ValueError("failure audit contains no failed cells")
    if not isinstance(invariants, Mapping):
        raise ValueError("failure audit contains no invariants")
    audited_keys = set()
    candidate_first = 0
    counts = Counter()
    for cell in cells:
        if not isinstance(cell, Mapping):
            raise ValueError("failure audit contains a malformed cell")
        key = (cell.get("panel"), cell.get("benchmark"), cell.get("repetition"))
        if (
            key in audited_keys
            or key[0] not in {"full", "e2e"}
            or not isinstance(key[1], str)
            or not isinstance(key[2], int)
            or isinstance(key[2], bool)
        ):
            raise ValueError(f"failure audit contains malformed or duplicate key {key!r}")
        audited_keys.add(key)
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
        if not isinstance(cell.get("candidate_completed_before_control_failure"), bool):
            raise ValueError(f"failure audit key {key!r} has malformed candidate-first state")
        candidate_first += int(cell["candidate_completed_before_control_failure"])
    if document.get("failed_pair_count") != len(cells):
        raise ValueError("failure audit failed_pair_count disagrees with cells")
    missing_keys = metrics["campaign"].get("missing_keys")
    if not isinstance(missing_keys, list) or any(
        not isinstance(key, list) or len(key) != 3 for key in missing_keys
    ):
        raise ValueError("release metrics missing_keys are malformed")
    metrics_missing = {tuple(key) for key in missing_keys}
    if audited_keys != metrics_missing:
        raise ValueError("failure audit keys disagree with release metrics missing_keys")
    planned_keys = _planned_pair_keys(metrics)
    observed_keys = _observed_pair_keys(metrics)
    if observed_keys & audited_keys or observed_keys | audited_keys != planned_keys:
        raise ValueError("observed and failed pair keys do not partition the planned matrix")
    if metrics["campaign"].get("n_pairs") != len(observed_keys):
        raise ValueError("release metrics n_pairs disagrees with observed repetitions")
    expected_counts = [
        {"panel": panel, "benchmark": benchmark, "failed_repetitions": count}
        for (panel, benchmark), count in sorted(counts.items())
    ]
    if document.get("failed_case_counts") != expected_counts:
        raise ValueError("failure audit failed_case_counts disagree with cells")
    required_true = (
        "all_returncode_1",
        "all_watchdog_reason_null",
        "all_feature_not_found",
        "all_exception_string_typeerror",
    )
    if any(invariants.get(name) is not True for name in required_true):
        raise ValueError("failure audit deterministic-failure invariants do not all hold")
    if invariants.get("candidate_first_outputs_preserved") != candidate_first:
        raise ValueError("failure audit candidate-first invariant disagrees with cells")
    if len(cells) != 28 or len(counts) != 7 or candidate_first != 14:
        raise ValueError("failure audit disagrees with the sealed exact-release boundary")
    return document


def validate_manifest(document: Any, metrics: Mapping[str, Any]) -> Mapping[str, Any]:
    metrics = validate_metrics(metrics)
    if not isinstance(document, Mapping) or document.get("schema_version") != 4:
        raise ValueError("release manifest must be a schema-version 4 object")
    campaign = metrics["campaign"]
    for manifest_key, metrics_key in (
        ("candidate_sha", "candidate_sha"),
        ("reference_sha", "reference_sha"),
    ):
        if document.get(manifest_key) != campaign.get(metrics_key):
            raise ValueError(f"release manifest {manifest_key} disagrees with metrics")
    mode = (metrics.get("publication") or {}).get("mode")
    if document.get("publication_mode") != mode:
        raise ValueError("release manifest publication mode disagrees with metrics")
    pair_results = document.get("pair_results")
    if not isinstance(pair_results, list):
        raise ValueError("release manifest pair_results are malformed")
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
    metrics_overlays = (metrics.get("publication") or {}).get("evaluation_overlays") or {}
    if (
        manifest_overlays.get("validated") != metrics_overlays.get("validated")
        or manifest_overlays.get("pair_count") != metrics_overlays.get("pair_count")
    ):
        raise ValueError("release manifest evaluation overlays disagree with metrics")
    manifest_roots = set(document.get("run_roots") or [])
    metric_roots = {
        root.get("root")
        for root in ((metrics.get("publication") or {}).get("controller_evidence") or {}).get("roots", [])
        if isinstance(root, Mapping)
    }
    if manifest_roots != metric_roots:
        raise ValueError("release manifest controller roots disagree with metrics")
    if document.get("finalization_host") != (metrics.get("environment") or {}).get("finalization_host"):
        raise ValueError("release manifest finalization host disagrees with metrics")
    return document


def validate_legacy_results(document: Any) -> Mapping[str, Any]:
    if not isinstance(document, Mapping):
        raise ValueError("legacy result document must be an object")
    full = [value for key, value in document.items() if key.startswith("full:")]
    if len(full) != 17 or any(not isinstance(cell, Mapping) for cell in full):
        raise ValueError("historical appendix requires the frozen 17-genome result set")
    return document


def _panels(metrics: Mapping[str, Any]) -> list[tuple[str, Mapping[str, Any]]]:
    return sorted(metrics["panels"].items(), key=lambda item: PANEL_ORDER.get(item[0], 99))


def _save(fig: Any, path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(
        path,
        format="png",
        dpi=300,
        bbox_inches="tight",
        facecolor="white",
        metadata={"Software": "LiftOn release_report_figures.py"},
    )
    plt.close(fig)
    return path


def _ratio_interval(summary: Mapping[str, Any], field: str) -> tuple[float, float | None, float | None]:
    interval = summary[field]
    estimate = float(interval["estimate"])
    if int(summary.get("n_cells", 0)) < 2:
        return estimate, None, None
    return estimate, estimate - float(interval["low"]), float(interval["high"]) - estimate


def _delta_error(interval: Mapping[str, Any]) -> tuple[float, float, float]:
    estimate = float(interval["estimate"])
    return estimate, estimate - float(interval["low"]), float(interval["high"]) - estimate


def _quality_rows(metrics: Mapping[str, Any]) -> list[tuple[str, str, Mapping[str, Any], int]]:
    rows = []
    for panel, summary in _panels(metrics):
        quality = summary.get("quality") or {}
        for metric in PRIMARY_QUALITY:
            interval = quality.get(metric)
            if isinstance(interval, Mapping):
                rows.append((panel, metric, interval, int(summary.get("n_cells", 0))))
    return rows


def figure_study_design(metrics: Mapping[str, Any], audit: Mapping[str, Any], output: Path) -> Path:
    planned = metrics["campaign"]["planned_campaign"]["panels"]
    panels = sorted(PANEL_ORDER, key=PANEL_ORDER.get)
    planned_cases = [len(planned[panel]["ids"]) for panel in panels]
    observed_cases = [sum(cell["panel"] == panel for cell in metrics["cells"]) for panel in panels]
    planned_pairs = [len(planned[panel]["ids"]) * int(planned[panel]["repetitions"]) for panel in panels]
    observed_pairs = [
        sum(int(cell["repetitions"]) for cell in metrics["cells"] if cell["panel"] == panel)
        for panel in panels
    ]
    missing_pairs = [planned_value - observed_value for planned_value, observed_value in zip(planned_pairs, observed_pairs)]
    positions = np.arange(len(panels))
    labels = [PANEL_LABELS[panel] for panel in panels]
    colors = [PANEL_COLORS[panel] for panel in panels]

    fig, axes = plt.subplots(1, 3, figsize=(14.6, 4.8), gridspec_kw={"width_ratios": [1.0, 1.08, 1.0]})
    ax = axes[0]
    ax.barh(positions, planned_cases, color="#dedede", edgecolor="#888888", label="planned")
    observed_bars = ax.barh(positions, observed_cases, color=colors, label="observed")
    ax.bar_label(observed_bars, labels=[f"{value}/{planned}" for value, planned in zip(observed_cases, planned_cases)], padding=3)
    ax.set_yticks(positions, labels)
    ax.invert_yaxis()
    ax.set_xlabel("logical benchmark cases")
    ax.set_title("A. Planned and observed cases", loc="left", fontweight="bold")
    ax.legend(frameon=False, fontsize=8)
    ax.grid(axis="x", alpha=0.2)

    ax = axes[1]
    ax.barh(positions, observed_pairs, color=colors, label="successful paired repetitions")
    missing_bars = ax.barh(positions, missing_pairs, left=observed_pairs, color=TOOL_COLORS["reference"], alpha=0.8, label="missing after v1.0.8 failure")
    for position, observed, missing, planned_count in zip(positions, observed_pairs, missing_pairs, planned_pairs):
        ax.text(observed / 2, position, str(observed), ha="center", va="center", color="white", fontsize=9, fontweight="bold")
        if missing:
            ax.text(observed + missing / 2, position, str(missing), ha="center", va="center", color="white", fontsize=9, fontweight="bold")
        ax.text(planned_count + 2, position, f"of {planned_count}", va="center", fontsize=8)
    ax.set_yticks(positions, labels)
    ax.invert_yaxis()
    ax.set_xlim(0, max(planned_pairs) * 1.18)
    ax.set_xlabel("AB/BA paired repetitions")
    ax.set_title("B. Observed repetitions and controls", loc="left", fontweight="bold")
    ax.legend(frameon=False, fontsize=7.5, loc="lower right")
    ax.grid(axis="x", alpha=0.2)

    ax = axes[2]
    total = len(metrics["cells"])
    validity = sum(cell.get("candidate_valid") is True for cell in metrics["cells"])
    deterministic = sum(
        cell.get("candidate_byte_deterministic") is True
        and cell.get("candidate_semantic_deterministic") is True
        and cell.get("candidate_quality_deterministic") is True
        and cell.get("common_pi_deterministic") is True
        for cell in metrics["cells"]
    )
    first = int(audit["invariants"]["candidate_first_outputs_preserved"])
    status_labels = ["Candidate valid", "Candidate deterministic", "Candidate-first\noutputs preserved"]
    values = [validity, deterministic, first]
    denominators = [total, total, int(audit["failed_pair_count"])]
    status_positions = np.arange(3)
    status_bars = ax.barh(status_positions, [100 * value / denominator for value, denominator in zip(values, denominators)], color=[TOOL_COLORS["candidate"], TOOL_COLORS["candidate"], "#72b7b2"])
    ax.bar_label(status_bars, labels=[f"{value}/{denominator}" for value, denominator in zip(values, denominators)], padding=3)
    ax.axvline(100, color="#333333", linewidth=1)
    ax.set_xlim(0, 112)
    ax.set_yticks(status_positions, status_labels)
    ax.invert_yaxis()
    ax.set_xlabel("percent of applicable evidence")
    ax.set_title("C. Integrity and retained evidence", loc="left", fontweight="bold")
    ax.grid(axis="x", alpha=0.2)
    ax.text(0.02, -0.28, "28/28 missing controls: deterministic historical-reference failures", transform=ax.transAxes, fontsize=8, color=TOOL_COLORS["reference"])

    fig.suptitle("Exact v1.0.11 diagnostic design and evidence status", fontsize=14, fontweight="bold", y=1.03)
    fig.tight_layout()
    return _save(fig, output)


def _arm_quality(cell: Mapping[str, Any], arm: str) -> Mapping[str, Any]:
    absolute = cell.get("absolute_quality") or {}
    value = absolute.get(arm) if isinstance(absolute, Mapping) else None
    if isinstance(value, Mapping):
        return value
    value = cell.get(arm)
    return value if isinstance(value, Mapping) else {}


def figure_comprehensiveness(metrics: Mapping[str, Any], output: Path,
                             candidate_label: str, reference_label: str) -> Path:
    panels = _panels(metrics)
    fig = plt.figure(figsize=(15.2, 8.0))
    grid = fig.add_gridspec(2, 2, height_ratios=[0.9, 1.1], width_ratios=[1.0, 1.3], hspace=0.38, wspace=0.27)
    delta_axis = fig.add_subplot(grid[0, 0])
    stable_axis = fig.add_subplot(grid[1, 0])
    e2e_axis = fig.add_subplot(grid[:, 1])

    metric_names = ("completeness_coding", "completeness_feature_total")
    offsets = np.linspace(-0.20, 0.20, len(panels))
    positions = np.arange(len(metric_names))
    for offset, (panel, summary) in zip(offsets, panels):
        estimates, lower, upper = [], [], []
        for metric in metric_names:
            estimate, low, high = _delta_error(summary["quality"][metric])
            estimates.append(estimate)
            lower.append(low)
            upper.append(high)
        delta_axis.errorbar(estimates, positions + offset, xerr=[lower, upper], fmt="o", color=PANEL_COLORS[panel], capsize=3, label=PANEL_LABELS[panel])
    delta_axis.axvline(0, color="#333333", linewidth=1)
    delta_axis.axvline(-0.001, color="#d62728", linestyle="--", linewidth=1, label="−0.001 panel floor")
    delta_axis.set_yticks(positions, [QUALITY_LABELS[name] for name in metric_names])
    delta_axis.invert_yaxis()
    delta_axis.set_xlabel(f"{candidate_label} − {reference_label}")
    delta_axis.set_title("A. Panel-level comprehensiveness deltas", loc="left", fontweight="bold")
    delta_axis.legend(frameon=False, fontsize=7.2, loc="lower right")
    delta_axis.grid(axis="x", alpha=0.2)

    stable_rows = []
    for panel, summary in panels:
        for feature_type in ("CDS", "exon"):
            record = (summary.get("stable_id_preservation") or {}).get(feature_type)
            if isinstance(record, Mapping):
                stable_rows.append((panel, feature_type, record))
    stable_positions = np.arange(len(stable_rows))
    height = 0.34
    candidate_values = [100 * float(record["candidate_preservation_rate"]) for _panel, _feature, record in stable_rows]
    reference_values = [100 * float(record["reference_preservation_rate"]) for _panel, _feature, record in stable_rows]
    stable_axis.barh(stable_positions - height / 2, candidate_values, height, color=TOOL_COLORS["candidate"], label=candidate_label)
    stable_axis.barh(stable_positions + height / 2, reference_values, height, color=TOOL_COLORS["reference"], label=reference_label)
    stable_axis.set_yticks(stable_positions, [f"{PANEL_LABELS[panel]} · {feature}" for panel, feature, _record in stable_rows], fontsize=8)
    stable_axis.invert_yaxis()
    stable_axis.set_xlabel("stable-ID preservation (%)")
    stable_axis.set_title("B. CDS and exon stable identifiers", loc="left", fontweight="bold")
    stable_axis.set_xlim(0, 105)
    stable_axis.legend(frameon=False, fontsize=8, loc="lower right")
    stable_axis.grid(axis="x", alpha=0.2)

    e2e_cells = sorted((cell for cell in metrics["cells"] if cell["panel"] == "e2e"), key=lambda cell: cell["benchmark"])
    e2e_metrics = ("completeness_coding", "completeness_feature_total")
    group_positions = np.arange(len(e2e_cells))
    bar_width = 0.18
    for metric_index, metric in enumerate(e2e_metrics):
        for arm_index, (arm, color, hatch) in enumerate((
            ("candidate", TOOL_COLORS["candidate"], None),
            ("reference", TOOL_COLORS["reference"], "//"),
        )):
            shift = (metric_index * 2 + arm_index - 1.5) * bar_width
            values = [100 * float(_arm_quality(cell, arm)[metric]) for cell in e2e_cells]
            label = f"{candidate_label if arm == 'candidate' else reference_label} · {QUALITY_LABELS[metric]}"
            e2e_axis.bar(group_positions + shift, values, bar_width, color=color, alpha=1.0 if metric_index == 0 else 0.58, hatch=hatch, label=label)
    e2e_axis.set_xticks(group_positions, [cell["benchmark"].title() for cell in e2e_cells])
    e2e_axis.set_ylabel("absolute recovery (%)")
    e2e_axis.set_ylim(0, 106)
    e2e_axis.set_title("C. End-to-end absolute recovery", loc="left", fontweight="bold")
    e2e_axis.legend(
        frameon=False,
        fontsize=7.2,
        ncol=2,
        loc="upper center",
        bbox_to_anchor=(0.5, -0.08),
    )
    e2e_axis.grid(axis="y", alpha=0.2)
    for index, cell in enumerate(e2e_cells):
        candidate = _arm_quality(cell, "candidate")
        reference = _arm_quality(cell, "reference")
        recovered = f"coding {int(candidate['n_recovered_coding']):,}/{int(reference['n_recovered_coding']):,}\nof {int(candidate['n_reference_coding']):,}"
        e2e_axis.text(index, 3, recovered, ha="center", va="bottom", fontsize=7)

    fig.suptitle("Comprehensiveness and identifier preservation", fontsize=14, fontweight="bold", y=0.99)
    return _save(fig, output)


def figure_accuracy(metrics: Mapping[str, Any], output: Path,
                    candidate_label: str, reference_label: str) -> Path:
    panels = _panels(metrics)
    fig = plt.figure(figsize=(15.2, 8.2))
    grid = fig.add_gridspec(2, 2, height_ratios=[1.2, 0.8], width_ratios=[1.35, 1.0], hspace=0.36, wspace=0.30)
    delta_axis = fig.add_subplot(grid[:, 0])
    tradeoff_axis = fig.add_subplot(grid[0, 1])
    common_axis = fig.add_subplot(grid[1, 1])

    positions = np.arange(len(ACCURACY_QUALITY))
    offsets = np.linspace(-0.22, 0.22, len(panels))
    for offset, (panel, summary) in zip(offsets, panels):
        estimates, lower, upper = [], [], []
        for metric in ACCURACY_QUALITY:
            estimate, low, high = _delta_error(summary["quality"][metric])
            estimates.append(estimate)
            lower.append(low)
            upper.append(high)
        delta_axis.errorbar(estimates, positions + offset, xerr=[lower, upper], fmt="o", markersize=5, color=PANEL_COLORS[panel], capsize=3, label=PANEL_LABELS[panel])
    delta_axis.axvline(0, color="#333333", linewidth=1)
    delta_axis.axvline(-0.001, color="#d62728", linestyle="--", linewidth=1, label="−0.001 panel floor")
    delta_axis.set_yticks(positions, [QUALITY_LABELS[name] for name in ACCURACY_QUALITY], fontsize=8.5)
    delta_axis.invert_yaxis()
    delta_axis.set_xlabel(f"{candidate_label} − {reference_label}")
    delta_axis.set_title("A. Accuracy and structural-quality deltas", loc="left", fontweight="bold")
    delta_axis.legend(frameon=False, fontsize=7.5, loc="lower right")
    delta_axis.grid(axis="x", alpha=0.2)

    e2e_cells = sorted((cell for cell in metrics["cells"] if cell["panel"] == "e2e"), key=lambda cell: cell["benchmark"])
    for cell in e2e_cells:
        candidate = _arm_quality(cell, "candidate")
        reference = _arm_quality(cell, "reference")
        x0, y0 = 100 * float(reference["completeness_coding"]), 100 * float(reference["mean_pi"])
        x1, y1 = 100 * float(candidate["completeness_coding"]), 100 * float(candidate["mean_pi"])
        tradeoff_axis.annotate("", xy=(x1, y1), xytext=(x0, y0), arrowprops={"arrowstyle": "->", "color": "#777777", "lw": 1.2})
        tradeoff_axis.scatter(x0, y0, color=TOOL_COLORS["reference"], marker="s", s=45, edgecolor="white", linewidth=0.5)
        tradeoff_axis.scatter(x1, y1, color=TOOL_COLORS["candidate"], marker="o", s=55, edgecolor="white", linewidth=0.5)
        tradeoff_axis.annotate(cell["benchmark"].title(), (x1, y1), xytext=(5, 4), textcoords="offset points", fontsize=8)
    tradeoff_axis.scatter([], [], color=TOOL_COLORS["reference"], marker="s", label=reference_label)
    tradeoff_axis.scatter([], [], color=TOOL_COLORS["candidate"], marker="o", label=candidate_label)
    tradeoff_axis.set_xlabel("coding recovery (%)")
    tradeoff_axis.set_ylabel("mean PI among scored recovery (%)")
    tradeoff_axis.set_title("B. E2E accuracy–recovery tradeoff", loc="left", fontweight="bold")
    tradeoff_axis.legend(frameon=False, fontsize=8)
    tradeoff_axis.grid(alpha=0.2)

    common_values = {
        panel: [float((cell.get("common_pi") or {}).get("mean_delta")) for cell in metrics["cells"] if cell["panel"] == panel and _number((cell.get("common_pi") or {}).get("mean_delta"))]
        for panel in PANEL_ORDER
    }
    for position, panel in enumerate(sorted(PANEL_ORDER, key=PANEL_ORDER.get)):
        values = common_values[panel]
        jitter = np.linspace(-0.09, 0.09, len(values)) if values else []
        common_axis.scatter(np.full(len(values), position) + jitter, values, color=PANEL_COLORS[panel], alpha=0.75, s=25, edgecolor="white", linewidth=0.35)
        if values:
            common_axis.plot([position - 0.22, position + 0.22], [np.median(values)] * 2, color="#222222", linewidth=2)
    common_axis.axhline(0, color="#333333", linewidth=1)
    common_axis.axhline(-0.005, color="#d62728", linestyle="--", linewidth=1, label="−0.005 cell floor")
    common_axis.set_xticks(np.arange(3), [PANEL_LABELS[panel] for panel in sorted(PANEL_ORDER, key=PANEL_ORDER.get)], rotation=12, ha="right", fontsize=8)
    common_axis.set_ylabel("common-recovered-set mean PI delta")
    common_axis.set_title("C. Dataset-level common-set PI", loc="left", fontweight="bold")
    common_axis.legend(frameon=False, fontsize=8)
    common_axis.grid(axis="y", alpha=0.2)

    fig.suptitle("Reference-content retention and sequence/structure plausibility", fontsize=14, fontweight="bold", y=0.99)
    return _save(fig, output)


def figure_performance(metrics: Mapping[str, Any], output: Path,
                       candidate_label: str, reference_label: str) -> Path:
    panels = _panels(metrics)
    panel_names = [panel for panel, _summary in panels]
    fig = plt.figure(figsize=(15.0, 8.2))
    grid = fig.add_gridspec(2, 2, hspace=0.37, wspace=0.26)
    wall_axis = fig.add_subplot(grid[0, 0])
    rss_axis = fig.add_subplot(grid[0, 1])
    aggregate_axis = fig.add_subplot(grid[1, 0])
    memory_axis = fig.add_subplot(grid[1, 1])
    outliers = {(row.get("panel"), row.get("benchmark")) for row in metrics.get("performance_outliers") or []}

    def distribution(axis: Any, field: str, summary_field: str, title: str) -> None:
        for position, (panel, summary) in enumerate(panels):
            values = [float(cell[field]) for cell in metrics["cells"] if cell["panel"] == panel]
            if len(values) > 1:
                parts = axis.violinplot(values, positions=[position], widths=0.68, showextrema=False)
                for body in parts["bodies"]:
                    body.set_facecolor(PANEL_COLORS[panel])
                    body.set_edgecolor(PANEL_COLORS[panel])
                    body.set_alpha(0.17)
            jitter = np.linspace(-0.16, 0.16, len(values)) if values else []
            cells = [cell for cell in metrics["cells"] if cell["panel"] == panel]
            normal = [(x, value) for x, value, cell in zip(jitter, values, cells) if (panel, cell["benchmark"]) not in outliers]
            flagged = [(x, value) for x, value, cell in zip(jitter, values, cells) if (panel, cell["benchmark"]) in outliers]
            if normal:
                axis.scatter([position + x for x, _value in normal], [value for _x, value in normal], s=24, color=PANEL_COLORS[panel], edgecolor="white", linewidth=0.35, alpha=0.85)
            if flagged:
                axis.scatter([position + x for x, _value in flagged], [value for _x, value in flagged], s=38, facecolor=PANEL_COLORS[panel], edgecolor="#d62728", linewidth=1.2, alpha=0.95, label="screened shared-host outlier" if position == 0 else None)
            estimate, low_error, high_error = _ratio_interval(summary, summary_field)
            if low_error is not None:
                axis.errorbar(position, estimate, yerr=[[low_error], [high_error]], fmt="D", color="#111111", markersize=5, capsize=4, linewidth=1.5, zorder=5)
        axis.axhline(1.0, color="#333333", linewidth=1)
        axis.axhline(1.10, color="#d62728", linestyle="--", linewidth=1, label="1.10 aggregate upper gate")
        axis.axhline(1.25, color="#d62728", linestyle=":", linewidth=1, label="1.25 worst-cell gate")
        axis.set_xticks(np.arange(len(panels)), [PANEL_LABELS[panel] for panel in panel_names], rotation=10, ha="right")
        axis.set_ylabel(f"{candidate_label} / {reference_label}")
        axis.set_title(title, loc="left", fontweight="bold")
        axis.grid(axis="y", alpha=0.2)
        axis.legend(frameon=False, fontsize=7.2, loc="upper right")

    distribution(wall_axis, "wall_ratio_median", "wall_ratio", "A. Per-cell wall ratios and panel GMR")
    distribution(rss_axis, "rss_ratio_median", "rss_ratio", "B. Per-cell command-RSS ratios and panel GMR")

    positions = np.arange(len(panels))
    width = 0.34
    gmr = [float(summary["wall_ratio"]["estimate"]) for _panel, summary in panels]
    total = [float(summary["wall_total_ratio"]) for _panel, summary in panels]
    aggregate_axis.bar(positions - width / 2, gmr, width, color=[PANEL_COLORS[panel] for panel in panel_names], label="cell-level wall GMR")
    aggregate_axis.bar(positions + width / 2, total, width, color=[PANEL_COLORS[panel] for panel in panel_names], alpha=0.45, hatch="//", label="total work ratio")
    aggregate_axis.axhline(1.0, color="#333333", linewidth=1)
    aggregate_axis.axhline(1.10, color="#d62728", linestyle="--", linewidth=1)
    aggregate_axis.set_xticks(positions, [PANEL_LABELS[panel] for panel in panel_names], rotation=10, ha="right")
    aggregate_axis.set_ylabel(f"{candidate_label} / {reference_label}")
    aggregate_axis.set_title("C. Aggregate and total-work wall ratios", loc="left", fontweight="bold")
    aggregate_axis.legend(frameon=False, fontsize=8)
    aggregate_axis.grid(axis="y", alpha=0.2)
    for position, (_panel, summary) in enumerate(panels):
        candidate_seconds = summary.get("candidate_wall_seconds_total_of_cell_medians")
        reference_seconds = summary.get("reference_wall_seconds_total_of_cell_medians")
        if _number(candidate_seconds) and _number(reference_seconds):
            candidate_hours = float(candidate_seconds) / 3600
            reference_hours = float(reference_seconds) / 3600
            aggregate_axis.text(position, 0.03, f"{candidate_hours:.1f}/{reference_hours:.1f} h", ha="center", va="bottom", fontsize=7)

    envelopes = [float(summary.get("candidate_concurrent_memory_envelope_gib", 0.0)) for _panel, summary in panels]
    peaks = [float(summary.get("candidate_peak_rss_gib", 0.0)) for _panel, summary in panels]
    memory_axis.barh(positions - 0.18, envelopes, 0.34, color=[PANEL_COLORS[panel] for panel in panel_names], label="slot-sum proxy")
    memory_axis.barh(positions + 0.18, peaks, 0.34, color=[PANEL_COLORS[panel] for panel in panel_names], alpha=0.42, hatch="//", label="largest command RSS")
    memory_axis.axvline(192.0, color="#d62728", linewidth=1.2, linestyle="--", label="192 GiB admission gate")
    memory_axis.set_yticks(positions, [PANEL_LABELS[panel] for panel in panel_names])
    memory_axis.invert_yaxis()
    memory_axis.set_xlabel("GiB (GNU-time command-RSS proxy)")
    memory_axis.set_title("D. Concurrency memory proxies", loc="left", fontweight="bold")
    memory_axis.set_xlim(0, 205)
    memory_axis.legend(frameon=False, fontsize=8, loc="lower right")
    memory_axis.grid(axis="x", alpha=0.2)
    for position, value in enumerate(envelopes):
        memory_axis.text(value + 2, position - 0.18, f"{value:.1f}", va="center", fontsize=8)

    context, verdict = _figure_context(metrics)
    fig.suptitle(f"LiftOn v1.0.11 {context}: performance — {verdict}", fontsize=14, fontweight="bold", y=0.99)
    return _save(fig, output)


def figure_correctness(metrics: Mapping[str, Any], output: Path,
                       candidate_label: str, reference_label: str) -> Path:
    """Retain the original compact correctness figure for old callers."""
    cells = metrics["cells"]
    total = len(cells)
    fig = plt.figure(figsize=(14.2, 5.8))
    grid = fig.add_gridspec(1, 3, width_ratios=[1.45, 0.9, 0.95], wspace=0.38)
    quality_axis = fig.add_subplot(grid[0, 0])
    evidence_axis = fig.add_subplot(grid[0, 1])
    stable_axis = fig.add_subplot(grid[0, 2])

    quality_rows = _quality_rows(metrics)
    quality_positions = np.arange(len(quality_rows))
    for position, (panel, _metric, interval, n_cells) in enumerate(quality_rows):
        estimate = float(interval["estimate"])
        if n_cells < 2:
            quality_axis.plot(estimate, position, marker="D", linestyle="none", color=PANEL_COLORS.get(panel, "#777777"), markersize=5)
        else:
            quality_axis.errorbar(estimate, position, xerr=[[estimate - float(interval["low"])], [float(interval["high"]) - estimate]], fmt="o", color=PANEL_COLORS.get(panel, "#777777"), ecolor=PANEL_COLORS.get(panel, "#777777"), capsize=3, markersize=5)
    quality_axis.axvline(0.0, color="#333333", linewidth=1.0)
    quality_axis.axvline(-0.001, color="#d62728", linewidth=1.0, linestyle="--")
    quality_axis.set_yticks(quality_positions)
    quality_axis.set_yticklabels([f"{panel}: {QUALITY_LABELS.get(metric, metric)}" for panel, metric, _interval, _n_cells in quality_rows], fontsize=8)
    quality_axis.invert_yaxis()
    quality_axis.set_xlabel(f"{candidate_label} − {reference_label}")
    quality_axis.set_title("A. Quality deltas (95% CI where n≥2)", loc="left", fontweight="bold")
    quality_axis.grid(axis="x", alpha=0.25)

    evidence_names = ["GFF3 valid", "byte deterministic", "semantic deterministic", "quality deterministic"]
    candidate_counts = [
        sum(cell.get("candidate_valid") is True for cell in cells),
        sum(cell.get("candidate_byte_deterministic") is True for cell in cells),
        sum(cell.get("candidate_semantic_deterministic") is True for cell in cells),
        sum(cell.get("candidate_quality_deterministic") is True for cell in cells),
    ]
    reference_counts = [
        sum(cell.get("reference_valid") is True for cell in cells),
        sum(cell.get("reference_byte_deterministic") is True for cell in cells),
        sum(cell.get("reference_semantic_deterministic") is True for cell in cells),
        sum(cell.get("reference_quality_deterministic") is True for cell in cells),
    ]
    evidence_positions = np.arange(len(evidence_names))
    height = 0.36
    candidate_bars = evidence_axis.barh(evidence_positions - height / 2, candidate_counts, height, color=TOOL_COLORS["candidate"], label=candidate_label)
    reference_bars = evidence_axis.barh(evidence_positions + height / 2, reference_counts, height, color=TOOL_COLORS["reference"], label=reference_label)
    evidence_axis.bar_label(candidate_bars, labels=[f"{value}/{total}" for value in candidate_counts], fontsize=8)
    evidence_axis.bar_label(reference_bars, labels=[f"{value}/{total}" for value in reference_counts], fontsize=8)
    evidence_axis.set_yticks(evidence_positions, evidence_names, fontsize=8)
    evidence_axis.invert_yaxis()
    evidence_axis.set_xlim(0, max(total * 1.22, 1))
    evidence_axis.set_xlabel("cells")
    evidence_axis.set_title("B. Evidence integrity", loc="left", fontweight="bold")
    evidence_axis.legend(fontsize=7.5, loc="lower right", frameon=False)
    evidence_axis.grid(axis="x", alpha=0.25)

    stable_rows = []
    for panel, summary in _panels(metrics):
        for feature_type, row in (summary.get("stable_id_preservation") or {}).items():
            value = row.get("rate_delta")
            if _number(value):
                stable_rows.append((panel, feature_type, float(value)))
    stable_positions = np.arange(len(stable_rows))
    stable_axis.barh(stable_positions, [row[2] for row in stable_rows], color=[PANEL_COLORS.get(row[0], "#777777") for row in stable_rows])
    stable_axis.axvline(0.0, color="#333333", linewidth=1.0)
    stable_axis.set_yticks(stable_positions, [f"{panel}: {feature_type}" for panel, feature_type, _value in stable_rows], fontsize=8)
    stable_axis.invert_yaxis()
    stable_axis.set_xlabel("candidate − reference preservation rate")
    stable_axis.set_title("C. Stable feature IDs", loc="left", fontweight="bold")
    stable_axis.grid(axis="x", alpha=0.25)

    context, verdict = _figure_context(metrics)
    fig.suptitle(f"LiftOn v1.0.11 {context} — {verdict}", fontsize=13, fontweight="bold", y=1.03)
    return _save(fig, output)


def _historical_cells(document: Mapping[str, Any]) -> list[tuple[str, Mapping[str, Any]]]:
    tiers = {"same_species": 0, "cross_species": 1, "close_cross_species": 2, "distant_cross_species": 3, "very_distant_cross_species": 4}
    cells = [(key.split(":", 1)[1], value) for key, value in document.items() if key.startswith("full:")]
    return sorted(cells, key=lambda item: (tiers.get(item[1].get("divergence_class"), 99), item[0]))


def figure_historical_cross_tool(document: Mapping[str, Any], output: Path) -> Path:
    cells = _historical_cells(validate_legacy_results(document))
    labels = [name.replace("_to_", "→").replace("t1_", "").replace("t2_", "").replace("t3_", "").replace("t4_", "").replace("_", " ") for name, _cell in cells]
    positions = np.arange(len(cells))
    tools = (
        ("liftoff", "Liftoff"),
        ("miniprot", "miniprot"),
        ("lifton_stable", "LiftOn v1.0.8"),
        ("lifton_devel", "LiftOn functional SHA 7fc8419"),
    )
    fig, axes = plt.subplots(1, 2, figsize=(14.8, 9.2), sharey=True, gridspec_kw={"wspace": 0.08})
    for axis, field, title, scale in (
        (axes[0], "mean_pi", "A. Mean protein identity", 1.0),
        (axes[1], "completeness_coding", "B. Coding completeness", 100.0),
    ):
        for tool, label in tools:
            values = []
            for _name, cell in cells:
                value = (cell.get(field) or {}).get(tool)
                values.append(np.nan if not _number(value) else float(value) * scale)
            axis.plot(values, positions, marker="o", linewidth=1.2, markersize=4.2, color=TOOL_COLORS.get(tool.replace("lifton_devel", "candidate").replace("lifton_stable", "reference"), "#777777"), label=label)
        axis.set_title(title, loc="left", fontweight="bold")
        axis.grid(axis="x", alpha=0.22)
        axis.invert_yaxis()
    axes[0].set_yticks(positions, labels, fontsize=8)
    axes[0].set_xlabel("protein identity")
    axes[0].set_xlim(0, 1.03)
    axes[1].set_xlabel("reference coding transcripts recovered (%)")
    axes[1].set_xlim(0, 103)
    axes[1].legend(frameon=False, fontsize=8, loc="lower right")
    fig.suptitle("Historical 17-genome cross-tool study — measured LiftOn arm: functional SHA 7fc8419", fontsize=13.5, fontweight="bold", y=0.98)
    fig.text(0.5, 0.01, "Historical context only; these values are not exact-tag v1.0.11 qualification measurements.", ha="center", fontsize=9, color="#555555")
    return _save(fig, output)


def generate_figures(metrics: Mapping[str, Any], output_dir: Path,
                     candidate_label: str = "LiftOn v1.0.11",
                     reference_label: str = "LiftOn v1.0.8",
                     failure_audit: Mapping[str, Any] | None = None,
                     manifest: Mapping[str, Any] | None = None,
                     legacy_results: Mapping[str, Any] | None = None) -> list[Path]:
    validated = validate_metrics(metrics)
    output_dir = Path(output_dir)
    if failure_audit is None:
        if manifest is not None or legacy_results is not None:
            raise ValueError("--manifest and --legacy-results require --failure-audit")
        return [
            figure_performance(validated, output_dir / PERFORMANCE_FIGURE, candidate_label, reference_label),
            figure_correctness(validated, output_dir / CORRECTNESS_FIGURE, candidate_label, reference_label),
        ]
    audit = validate_failure_audit(failure_audit, validated)
    if manifest is not None:
        validate_manifest(manifest, validated)
    paths = [
        figure_study_design(validated, audit, output_dir / STUDY_FIGURE),
        figure_comprehensiveness(validated, output_dir / COMPREHENSIVENESS_FIGURE, candidate_label, reference_label),
        figure_accuracy(validated, output_dir / ACCURACY_FIGURE, candidate_label, reference_label),
        figure_performance(validated, output_dir / PERFORMANCE_FIGURE, candidate_label, reference_label),
    ]
    if legacy_results is not None:
        paths.append(figure_historical_cross_tool(legacy_results, output_dir / HISTORICAL_FIGURE))
    return paths


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    evidence = parser.add_mutually_exclusive_group(required=True)
    evidence.add_argument("--metrics", type=Path)
    evidence.add_argument(
        "--whole-genome-metrics", type=Path,
        help="seven-transfer biology-first metrics",
    )
    parser.add_argument("--failure-audit", type=Path)
    parser.add_argument("--manifest", type=Path)
    parser.add_argument("--legacy-results", type=Path)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--candidate-label", default="LiftOn v1.0.11")
    parser.add_argument("--reference-label", default="LiftOn v1.0.8")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        if args.whole_genome_metrics:
            if args.failure_audit or args.manifest or args.legacy_results:
                raise ValueError(
                    "whole-genome metrics cannot be combined with legacy "
                    "release-bundle options"
                )
            from . import whole_genome_assets

            paths = whole_genome_assets.generate_figures(
                whole_genome_assets.load_metrics(args.whole_genome_metrics),
                args.outdir,
            )
        else:
            paths = generate_figures(
                load_metrics(args.metrics),
                args.outdir,
                candidate_label=args.candidate_label,
                reference_label=args.reference_label,
                failure_audit=load_failure_audit(args.failure_audit) if args.failure_audit else None,
                manifest=load_manifest(args.manifest) if args.manifest else None,
                legacy_results=_load_json(args.legacy_results, "legacy results") if args.legacy_results else None,
            )
    except (OSError, TypeError, ValueError, RuntimeError) as exc:
        parser.error(str(exc))
    for path in paths:
        print(path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
