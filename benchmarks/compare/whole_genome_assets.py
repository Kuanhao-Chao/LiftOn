#!/usr/bin/env python3
"""Generate figures and evidence-derived tables for the seven-transfer study."""
from __future__ import annotations

import argparse
import difflib
import json
import math
import os
import re
from pathlib import Path
from typing import Any, Callable, Mapping, Sequence

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from . import whole_genome_report, whole_genome_study


CANDIDATE = "LiftOn v1.0.11"
REFERENCE = "LiftOn v1.0.8"
COLORS = {"candidate": "#54a24b", "reference": "#9d755d"}
SHORT_LABELS = {
    "bee": "Honey bee",
    "drosophila": "Fruit fly",
    "t2_human_to_gorilla": "Human → gorilla",
    "t2_mouse_to_caroli": "Mouse → M. caroli",
    "t3_dog_to_cat": "Dog → cat",
    "t3_human_to_macaque": "Human → macaque",
    "t3_human_to_marmoset": "Human → marmoset",
}
FIGURES = {
    "study": "rfig_release_study_design.png",
    "comprehensiveness": "rfig_release_comprehensiveness.png",
    "accuracy": "rfig_release_accuracy.png",
    "performance": "rfig_release_performance.png",
}
START = "{{/* BEGIN GENERATED: {name} */}}"
END = "{{/* END GENERATED: {name} */}}"
MARKER = re.compile(
    r"^\{\/\* BEGIN GENERATED: (?P<name>biology-[a-z0-9-]+) \*\/\}$",
    re.MULTILINE,
)


class AssetError(RuntimeError):
    """Raised when report evidence or a generated region is inconsistent."""


def _file_record(record: Any, label: str) -> Path:
    if not isinstance(record, Mapping):
        raise AssetError(f"{label} file record is missing")
    path = Path(str(record.get("path", ""))).resolve()
    try:
        observed = whole_genome_study.file_record(path)
    except OSError as exc:
        raise AssetError(f"{label} cannot be read: {exc}") from exc
    if observed != record:
        raise AssetError(f"{label} changed after reduction")
    return path


def _load_json(path: Path, label: str) -> Mapping[str, Any]:
    try:
        document = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise AssetError(f"cannot read {label}: {exc}") from exc
    if not isinstance(document, Mapping):
        raise AssetError(f"{label} must be a JSON object")
    return document


def _validate_fingerprint(document: Mapping[str, Any], label: str) -> None:
    material = dict(document)
    observed = material.pop("fingerprint", None)
    material.pop("created_at", None)
    if observed != whole_genome_study.canonical_sha256(material):
        raise AssetError(f"{label} fingerprint is invalid")


def validate_metrics(document: Any) -> Mapping[str, Any]:
    if (
        not isinstance(document, Mapping)
        or document.get("schema_version") != whole_genome_report.SCHEMA_VERSION
        or document.get("method") != whole_genome_report.METHOD
    ):
        raise AssetError("study metrics use an unsupported schema or method")
    _validate_fingerprint(document, "study metrics")
    study_path = _file_record(document.get("study"), "study manifest")
    study = whole_genome_study.load_study(study_path)
    if (
        study["candidate"]["sha"]
        != "c623f0bddc5b5051a3670a3b0064c52cc61bb719"
        or study["reference"]["sha"]
        != "e503643d8346c600fedabcd3a4dff5c0873a4a37"
    ):
        raise AssetError("study manifest does not name the exact releases")
    _file_record(document.get("preflight"), "study preflight")
    _file_record(document.get("canary"), "fast/safe canary")
    if document.get("provider_ortholog_lock") is not None:
        _file_record(document["provider_ortholog_lock"], "provider lock")
    pairs = document.get("pairs")
    if (
        not isinstance(pairs, list)
        or [row.get("id") for row in pairs]
        != list(whole_genome_study.EXPECTED_PAIR_IDS)
    ):
        raise AssetError("study metrics do not contain the exact seven transfers")
    public_labels = {pair["id"]: pair["public_label"] for pair in study["pairs"]}
    for row in pairs:
        if row.get("public_label") != public_labels[row["id"]]:
            raise AssetError(f"{row['id']}: public label disagrees with study")
        records = row.get("pair_results")
        if not isinstance(records, list) or len(records) != 4:
            raise AssetError(f"{row['id']}: paired result records are incomplete")
        for index, record in enumerate(records, start=1):
            _file_record(record, f"{row['id']} repetition {index}")
    coverage = document.get("coverage", {})
    if (
        coverage.get("planned_pairs") != 28
        or coverage.get("planned_genome_transfers") != 7
        or coverage.get("observed_pairs") != len(pairs) * 4
        or coverage.get("successful_genome_transfers") != len(pairs)
    ):
        raise AssetError("study coverage disagrees with paired results")
    qualification = document.get("qualification", {})
    if (
        qualification.get("scope") != "this seven-transfer cohort only"
        or qualification.get("whole_release_claim") != "DIAGNOSTIC ONLY"
        or qualification.get("status") not in {"PASS", "FAIL"}
        or not isinstance(qualification.get("gates"), list)
        or not qualification["gates"]
    ):
        raise AssetError("study qualification boundary is missing or invalid")
    return document


def validate_sensitivity(
    document: Any,
    metrics: Mapping[str, Any],
) -> Mapping[str, Any]:
    if (
        not isinstance(document, Mapping)
        or document.get("schema_version") != 1
        or document.get("method")
        != "released-target-annotation-sensitivity-v1"
    ):
        raise AssetError("annotation sensitivity uses an unsupported schema")
    _validate_fingerprint(document, "annotation sensitivity")
    for key in ("study", "preflight", "provider_ortholog_lock"):
        if document.get(key) != metrics.get(
            "provider_ortholog_lock" if key == "provider_ortholog_lock" else key
        ):
            raise AssetError(f"annotation sensitivity {key} binding disagrees")
        _file_record(document[key], f"annotation sensitivity {key}")
    if [row.get("id") for row in document.get("pairs", [])] != list(
        whole_genome_study.EXPECTED_PAIR_IDS
    ):
        raise AssetError("annotation sensitivity transfer inventory disagrees")
    return document


def load_metrics(path: Path) -> Mapping[str, Any]:
    return validate_metrics(_load_json(path, "study metrics"))


def load_sensitivity(
    path: Path,
    metrics: Mapping[str, Any],
) -> Mapping[str, Any]:
    return validate_sensitivity(
        _load_json(path, "annotation sensitivity"), metrics,
    )


def _setup() -> None:
    plt.rcParams.update({
        "font.size": 9.2,
        "axes.titlesize": 11,
        "axes.labelsize": 9.5,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "figure.facecolor": "white",
        "axes.facecolor": "white",
        "savefig.facecolor": "white",
    })


def _save(fig: Any, path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(
        path,
        dpi=210,
        bbox_inches="tight",
        metadata={"Software": "LiftOn whole_genome_assets.py"},
    )
    plt.close(fig)
    if not path.is_file() or path.stat().st_size < 10_000:
        raise AssetError(f"generated figure is missing or unexpectedly small: {path}")
    return path


def _paired_bars(
    ax: Any,
    pairs: Sequence[Mapping[str, Any]],
    getter: Callable[[Mapping[str, Any], str], float],
    *,
    title: str,
    xlabel: str,
    limits: tuple[float, float] | None = (0, 1.01),
) -> None:
    y = np.arange(len(pairs))
    height = 0.34
    for offset, arm in ((-height / 2, "reference"), (height / 2, "candidate")):
        ax.barh(
            y + offset,
            [getter(row, arm) for row in pairs],
            height=height,
            color=COLORS[arm],
            label=REFERENCE if arm == "reference" else CANDIDATE,
        )
    ax.set_yticks(y, [SHORT_LABELS[row["id"]] for row in pairs])
    ax.invert_yaxis()
    ax.set_title(title, loc="left", fontweight="bold")
    ax.set_xlabel(xlabel)
    if limits is not None:
        ax.set_xlim(*limits)
    ax.grid(axis="x", color="#d7d7d7", linewidth=0.6, alpha=0.8)


def figure_study(metrics: Mapping[str, Any], output: Path) -> Path:
    metrics = validate_metrics(metrics)
    _setup()
    pairs = metrics["pairs"]
    fig = plt.figure(figsize=(15.2, 7.7))
    grid = fig.add_gridspec(1, 3, width_ratios=(1.25, 1.2, 1.05), wspace=0.34)
    ax = fig.add_subplot(grid[0, 0])
    ax.axis("off")
    ax.set_title("A  Biological range", loc="left", fontweight="bold")
    rows = []
    for pair in pairs:
        rows.append([
            SHORT_LABELS[pair["id"]],
            pair["biological_class"].replace(" transfer", ""),
        ])
    table = ax.table(
        cellText=rows,
        colLabels=["Transfer", "Biological rationale"],
        cellLoc="left",
        colLoc="left",
        loc="upper left",
        colWidths=(0.40, 0.60),
    )
    table.auto_set_font_size(False)
    table.set_fontsize(8.3)
    table.scale(1, 1.75)
    for (row, _column), cell in table.get_celld().items():
        cell.set_edgecolor("#d5d5d5")
        if row == 0:
            cell.set_facecolor("#eef4ee")
            cell.set_text_props(weight="bold")
    ax.text(
        0, 0.04,
        "One same-species assembly comparison, three close comparative\n"
        "transfers, and three moderate-divergence mammalian transfers.",
        transform=ax.transAxes,
        color="#4a4a4a",
        va="bottom",
    )

    ax = fig.add_subplot(grid[0, 1])
    repetitions = np.arange(1, 5)
    matrix = np.ones((len(pairs), 4))
    incomplete = {row["id"]: row["observed"] for row in metrics["coverage"]["incomplete"]}
    for row_index, pair in enumerate(pairs):
        observed = incomplete.get(pair["id"], 4)
        matrix[row_index, observed:] = 0
    ax.imshow(matrix, cmap=matplotlib.colors.ListedColormap(["#efefef", "#54a24b"]), aspect="auto")
    ax.set_xticks(repetitions - 1, [f"Rep {item}" for item in repetitions])
    ax.set_yticks(np.arange(len(pairs)), [SHORT_LABELS[row["id"]] for row in pairs])
    ax.set_title("B  Paired execution matrix", loc="left", fontweight="bold")
    for row_index in range(len(pairs)):
        for column in range(4):
            order = "BA" if (column + 1) % 2 else "AB"
            ax.text(column, row_index, order, ha="center", va="center", color="white", fontweight="bold")
    ax.set_xlabel("AB: v1.0.11 first · BA: v1.0.8 first")
    ax.tick_params(length=0)

    ax = fig.add_subplot(grid[0, 2])
    ax.axis("off")
    valid = sum(row["candidate"]["valid_repetitions"] for row in pairs)
    deterministic = sum(row["candidate"]["determinism"]["passed"] for row in pairs)
    status = metrics["qualification"]["status"]
    cards = [
        ("Transfers completed", f"{len(pairs)} / 7"),
        ("Paired repetitions", f"{metrics['coverage']['observed_pairs']} / 28"),
        ("Valid v1.0.11 outputs", f"{valid} / 28"),
        ("Deterministic transfers", f"{deterministic} / 7"),
        ("Cohort qualification", status),
    ]
    ax.set_title("C  Evidence status", loc="left", fontweight="bold")
    for index, (label, value) in enumerate(cards):
        top = 0.95 - index * 0.18
        ax.add_patch(plt.Rectangle(
            (0.02, top - 0.13), 0.96, 0.14,
            transform=ax.transAxes,
            facecolor="#f5f7f5" if index < 4 else "#eef4ee",
            edgecolor="#ccd6cc",
        ))
        ax.text(0.06, top - 0.035, label, transform=ax.transAxes, color="#555")
        ax.text(
            0.94, top - 0.06, value,
            transform=ax.transAxes, ha="right", va="center",
            fontsize=15, fontweight="bold",
            color=COLORS["candidate"] if value != "FAIL" else "#b23a33",
        )
    ax.text(
        0.02, 0.015,
        "The verdict applies only to this prespecified cohort; it is not a\n"
        "blanket qualification of every supported genome pair.",
        transform=ax.transAxes,
        color="#4a4a4a",
    )
    fig.suptitle(
        "Exact-release genome-to-genome annotation study",
        x=0.02, ha="left", fontsize=15, fontweight="bold",
    )
    return _save(fig, output)


def figure_comprehensiveness(metrics: Mapping[str, Any], output: Path) -> Path:
    metrics = validate_metrics(metrics)
    _setup()
    pairs = metrics["pairs"]
    fig, axes = plt.subplots(2, 2, figsize=(15.2, 10.2))
    _paired_bars(
        axes[0, 0], pairs,
        lambda row, arm: row[arm]["source"]["completeness_coding"],
        title="A  Reference coding transcripts recovered",
        xlabel="Fraction of source coding transcripts",
    )
    _paired_bars(
        axes[0, 1], pairs,
        lambda row, arm: row[arm]["source"]["completeness_feature_total"],
        title="B  All annotated source features recovered",
        xlabel="Fraction of counted source features",
    )
    _paired_bars(
        axes[1, 0], pairs,
        lambda row, arm: row[arm]["target_coordinate"]["gene_locus_recall"],
        title="C  Released target-annotation gene recall",
        xlabel="Ortholog-scoped locus recall",
    )
    y = np.arange(len(pairs))
    width = 0.20
    entries = [
        ("reference", "CDS", -1.5 * width, "v1.0.8 CDS", "#c9a995"),
        ("candidate", "CDS", -0.5 * width, "v1.0.11 CDS", "#78b872"),
        ("reference", "exon", 0.5 * width, "v1.0.8 exon", COLORS["reference"]),
        ("candidate", "exon", 1.5 * width, "v1.0.11 exon", COLORS["candidate"]),
    ]
    for arm, feature, offset, label, color in entries:
        values = [
            row[arm]["stable_ids"][feature]["preservation_rate"]
            if row[arm]["stable_ids"][feature]["applicable"] else 0
            for row in pairs
        ]
        axes[1, 1].barh(y + offset, values, width, label=label, color=color)
    axes[1, 1].set_yticks(y, [SHORT_LABELS[row["id"]] for row in pairs])
    axes[1, 1].invert_yaxis()
    axes[1, 1].set_xlim(0, 1.01)
    axes[1, 1].set_title("D  Stable child-feature identifiers", loc="left", fontweight="bold")
    axes[1, 1].set_xlabel("Preservation fraction when source IDs are present")
    axes[1, 1].grid(axis="x", color="#d7d7d7", linewidth=0.6)
    axes[1, 1].legend(frameon=False, ncol=2, fontsize=8)
    axes[0, 0].legend(frameon=False, ncol=2, loc="lower right")
    fig.suptitle(
        "Comprehensiveness across completed annotation transfers",
        x=0.02, ha="left", fontsize=15, fontweight="bold",
    )
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    return _save(fig, output)


def figure_accuracy(metrics: Mapping[str, Any], output: Path) -> Path:
    metrics = validate_metrics(metrics)
    _setup()
    pairs = metrics["pairs"]
    fig, axes = plt.subplots(2, 2, figsize=(15.2, 10.4))
    _paired_bars(
        axes[0, 0], pairs,
        lambda row, arm: row[arm]["source"]["covpi"],
        title="A  Coverage-weighted reference-protein identity",
        xlabel="CovPI (unrecovered coding transcripts contribute zero)",
    )
    y = np.arange(len(pairs))
    source_metrics = (
        ("intron_chain_exact_recall", "Exact intron chain"),
        ("orf_valid_recall", "ORF valid"),
    )
    for index, (metric, label) in enumerate(source_metrics):
        values = [row["source_deltas"][metric] for row in pairs]
        axes[0, 1].scatter(
            values, y + (index - 0.5) * 0.18,
            label=label,
            color=("#4c78a8", "#e45756")[index],
            s=44,
        )
    axes[0, 1].axvline(0, color="#555", linewidth=0.9)
    axes[0, 1].set_yticks(y, [SHORT_LABELS[row["id"]] for row in pairs])
    axes[0, 1].invert_yaxis()
    axes[0, 1].set_title("B  Source-structure change", loc="left", fontweight="bold")
    axes[0, 1].set_xlabel("v1.0.11 − v1.0.8 recall")
    axes[0, 1].grid(axis="x", color="#d7d7d7", linewidth=0.6)
    axes[0, 1].legend(frameon=False)
    _paired_bars(
        axes[1, 0], pairs,
        lambda row, arm: row[arm]["target_coordinate"]["transcript_locus_f1"],
        title="C  Released target-annotation transcript-locus agreement",
        xlabel="Ortholog-scoped F1 at 50% reciprocal overlap",
    )
    metrics_to_show = (
        ("intron_chain_f1", "Intron chain"),
        ("exon_recall", "Exon recall"),
        ("cds_recall", "CDS recall"),
    )
    x = np.arange(len(pairs))
    width = 0.25
    for index, (metric, label) in enumerate(metrics_to_show):
        axes[1, 1].bar(
            x + (index - 1) * width,
            [row["candidate"]["target_coordinate"][metric] for row in pairs],
            width,
            label=label,
            color=("#4c78a8", "#72b7b2", "#54a24b")[index],
        )
    axes[1, 1].set_xticks(
        x, [SHORT_LABELS[row["id"]] for row in pairs], rotation=28, ha="right",
    )
    axes[1, 1].set_ylim(0, 1.01)
    axes[1, 1].set_title("D  v1.0.11 target-structure concordance", loc="left", fontweight="bold")
    axes[1, 1].set_ylabel("Recall or F1")
    axes[1, 1].grid(axis="y", color="#d7d7d7", linewidth=0.6)
    axes[1, 1].legend(frameon=False, ncol=3, fontsize=8)
    axes[0, 0].legend(frameon=False, ncol=2, loc="lower right")
    fig.suptitle(
        "Sequence and structural quality of completed transfers",
        x=0.02, ha="left", fontsize=15, fontweight="bold",
    )
    fig.text(
        0.02, 0.006,
        "Target metrics quantify concordance with released annotations within fixed ortholog scopes; they do not establish biological locus truth.",
        color="#4a4a4a",
    )
    fig.tight_layout(rect=(0, 0.025, 1, 0.96))
    return _save(fig, output)


def figure_performance(metrics: Mapping[str, Any], output: Path) -> Path:
    metrics = validate_metrics(metrics)
    _setup()
    pairs = metrics["pairs"]
    fig, axes = plt.subplots(2, 2, figsize=(15.2, 10.2))
    y = np.arange(len(pairs))
    for ax, metric, title, gate in (
        (axes[0, 0], "wall", "A  Wall-time ratio", 1.25),
        (
            axes[0, 1], "process_group_rss",
            "B  Process-group RSS ratio", 1.25,
        ),
    ):
        medians = np.asarray([row["performance_ratios"][metric]["median"] for row in pairs])
        low = np.asarray([row["performance_ratios"][metric]["minimum"] for row in pairs])
        high = np.asarray([row["performance_ratios"][metric]["maximum"] for row in pairs])
        ax.errorbar(
            medians, y,
            xerr=np.vstack((medians - low, high - medians)),
            fmt="o", color=COLORS["candidate"], ecolor="#8fba8b",
            capsize=3,
        )
        ax.axvline(1, color="#555", linewidth=0.9, linestyle="--")
        if gate is not None:
            ax.axvline(gate, color="#b23a33", linewidth=0.9, linestyle=":")
        ax.set_yticks(y, [SHORT_LABELS[row["id"]] for row in pairs])
        ax.invert_yaxis()
        ax.set_xscale("log")
        ax.set_title(title, loc="left", fontweight="bold")
        ax.set_xlabel("v1.0.11 / v1.0.8 (range across four repetitions)")
        ax.grid(axis="x", color="#d7d7d7", linewidth=0.6)
    _paired_bars(
        axes[1, 0], pairs,
        lambda row, arm: row[arm]["performance"]["wall_seconds"]["median"] / 3600,
        title="C  Median elapsed time",
        xlabel="Hours per transfer",
        limits=None,
    )
    _paired_bars(
        axes[1, 1], pairs,
        lambda row, arm: row[arm]["performance"]["process_group_peak_rss_mb"]["median"] / 1024,
        title="D  Median process-group memory proxy",
        xlabel="GiB (summed RSS; shared pages may be counted more than once)",
        limits=None,
    )
    axes[1, 0].legend(frameon=False, ncol=2)
    aggregate = metrics["aggregate"]
    wall = aggregate["performance_ratios"]["wall"]
    rss = aggregate["performance_ratios"]["process_group_rss"]
    total_work = aggregate["total_work"]
    fig.suptitle(
        "Runtime and memory across completed annotation transfers",
        x=0.02, ha="left", fontsize=15, fontweight="bold",
    )
    fig.text(
        0.02, 0.006,
        f"Equal-transfer geometric means: wall {wall['estimate']:.3f}× "
        f"(95% CI {wall['low']:.3f}–{wall['high']:.3f}); process-group RSS "
        f"{rss['estimate']:.3f}× (95% CI {rss['low']:.3f}–{rss['high']:.3f}). "
        f"Two-concurrent-transfer v1.0.11 proxy: "
        f"{aggregate['candidate_concurrent_memory_proxy_gib']:.1f} GiB.",
        color="#4a4a4a",
    )
    fig.tight_layout(rect=(0, 0.03, 1, 0.96))
    return _save(fig, output)


def generate_figures(metrics: Mapping[str, Any], output_dir: Path) -> list[Path]:
    metrics = validate_metrics(metrics)
    functions = {
        "study": figure_study,
        "comprehensiveness": figure_comprehensiveness,
        "accuracy": figure_accuracy,
        "performance": figure_performance,
    }
    return [
        functions[name](metrics, output_dir / filename)
        for name, filename in FIGURES.items()
    ]


def _fmt(value: Any, digits: int = 4) -> str:
    if value is None:
        return "—"
    if isinstance(value, bool):
        return "yes" if value else "no"
    if isinstance(value, int):
        return f"{value:,}"
    if isinstance(value, float) and math.isfinite(value):
        return f"{value:.{digits}f}"
    return str(value)


def _pair_value(row: Mapping[str, Any], path: Sequence[str]) -> Any:
    current: Any = row
    for key in path:
        if not isinstance(current, Mapping):
            return None
        current = current.get(key)
    return current


def _two(row: Mapping[str, Any], path: Sequence[str], digits: int = 4) -> str:
    return " / ".join(
        _fmt(_pair_value(row[arm], path), digits)
        for arm in ("candidate", "reference")
    )


def table_cohort(metrics: Mapping[str, Any], _sensitivity: Mapping[str, Any] | None) -> str:
    rows = [
        "| Biological transfer | Rationale | Released target annotation | Reported genes | Protein-RBH gene / transcript groups |",
        "|---|---|---|---:|---:|",
    ]
    for row in metrics["pairs"]:
        target = row["target_annotation"]
        groups = row["candidate"]["target_scope"]
        rows.append(
            f"| {row['public_label']} | {row['biological_class']} | "
            f"{target['provider']}; `{target['assembly_accession']}`; "
            f"{target['release']} (evidence date {target['release_date']}) | "
            f"{target['reported_gene_count']:,} | "
            f"{groups['gene_groups']:,} / {groups['transcript_groups']:,} |"
        )
    return "\n".join(rows)


def headline_summary(
    metrics: Mapping[str, Any],
    _sensitivity: Mapping[str, Any] | None,
) -> str:
    pairs = metrics["pairs"]
    aggregate = metrics["aggregate"]
    wall = aggregate["performance_ratios"]["wall"]
    rss = aggregate["performance_ratios"]["process_group_rss"]
    total_work = aggregate["total_work"]
    valid = sum(row["candidate"]["valid_repetitions"] for row in pairs)
    deterministic = sum(
        row["candidate"]["determinism"]["passed"] for row in pairs
    )
    failed = sum(
        gate["passed"] is False
        for gate in metrics["qualification"]["gates"]
    )
    return (
        f"All {len(pairs)} prespecified transfers completed all "
        f"{metrics['coverage']['observed_pairs']} paired repetitions. "
        f"LiftOn v1.0.11 produced {valid}/28 valid outputs and was "
        f"deterministic for {deterministic}/7 transfers. Its equal-transfer "
        f"geometric-mean wall ratio was {wall['estimate']:.3f}× "
        f"(95% CI {wall['low']:.3f}–{wall['high']:.3f}) and its "
        f"process-group RSS ratio was {rss['estimate']:.3f}× "
        f"(95% CI {rss['low']:.3f}–{rss['high']:.3f}); values below one "
        f"favor v1.0.11. The summed median-duration ratio was "
        f"{total_work['wall_ratio']:.3f}×, and the two-concurrent-transfer "
        f"candidate memory proxy was "
        f"{aggregate['candidate_concurrent_memory_proxy_gib']:.1f} GiB. "
        f"The prespecified cohort verdict is "
        f"**{metrics['qualification']['status']}**, with {failed} failed "
        f"gate{'s' if failed != 1 else ''}. This verdict is restricted to "
        "the seven transfers and is diagnostic for the release as a whole."
    )


def table_comprehensiveness(metrics: Mapping[str, Any], _sensitivity: Mapping[str, Any] | None) -> str:
    rows = [
        "| Transfer | Coding completeness v1.0.11 / v1.0.8 | All-feature completeness | Target gene recall | Target transcript recall | CDS ID preservation | Exon ID preservation |",
        "|---|---:|---:|---:|---:|---:|---:|",
    ]
    for row in metrics["pairs"]:
        rows.append(
            f"| {SHORT_LABELS[row['id']]} | "
            f"{_two(row, ('source', 'completeness_coding'))} | "
            f"{_two(row, ('source', 'completeness_feature_total'))} | "
            f"{_two(row, ('target_coordinate', 'gene_locus_recall'))} | "
            f"{_two(row, ('target_coordinate', 'transcript_locus_recall'))} | "
            f"{_two(row, ('stable_ids', 'CDS', 'preservation_rate'))} | "
            f"{_two(row, ('stable_ids', 'exon', 'preservation_rate'))} |"
        )
    return "\n".join(rows)


def table_accuracy(metrics: Mapping[str, Any], _sensitivity: Mapping[str, Any] | None) -> str:
    rows = [
        "| Transfer | Source CovPI v1.0.11 / v1.0.8 | Exact intron recall | ORF-valid recall | Target gene F1 | Target transcript F1 | Target intron-chain F1 | Target protein identity |",
        "|---|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in metrics["pairs"]:
        rows.append(
            f"| {SHORT_LABELS[row['id']]} | {_two(row, ('source', 'covpi'))} | "
            f"{_two(row, ('source', 'intron_chain_exact_recall'))} | "
            f"{_two(row, ('source', 'orf_valid_recall'))} | "
            f"{_two(row, ('target_coordinate', 'gene_locus_f1'))} | "
            f"{_two(row, ('target_coordinate', 'transcript_locus_f1'))} | "
            f"{_two(row, ('target_coordinate', 'intron_chain_f1'))} | "
            f"{_two(row, ('target_sequence', 'mean_protein_identity'))} |"
        )
    labels = {
        "completeness_coding": "Source coding completeness",
        "completeness_feature_total": "Source all-feature completeness",
        "covpi": "Source CovPI",
        "recall_at_0.5": "Source recall at PI ≥ 0.50",
        "recall_at_0.75": "Source recall at PI ≥ 0.75",
        "recall_at_0.9": "Source recall at PI ≥ 0.90",
        "recall_at_0.95": "Source recall at PI ≥ 0.95",
        "intron_chain_exact_recall": "Source exact-intron recall",
        "orf_valid_recall": "Source ORF-valid recall",
        "gene_locus_f1": "Released-target gene-locus F1",
        "transcript_locus_f1": "Released-target transcript-locus F1",
    }
    rows.extend([
        "",
        "| Equal-transfer paired delta (v1.0.11 − v1.0.8) | Estimate (95% CI) |",
        "|---|---:|",
    ])
    for family in ("source_deltas", "target_deltas"):
        for metric, interval in metrics["aggregate"][family].items():
            rows.append(
                f"| {labels[metric]} | {interval['estimate']:+.5f} "
                f"({interval['low']:+.5f}–{interval['high']:+.5f}) |"
            )
    return "\n".join(rows)


def table_performance(metrics: Mapping[str, Any], _sensitivity: Mapping[str, Any] | None) -> str:
    rows = [
        "| Transfer | Median wall h v1.0.11 / v1.0.8 | Wall ratio | Median process-group RSS GiB | RSS ratio | Valid v1.0.11 reps | Deterministic v1.0.11 |",
        "|---|---:|---:|---:|---:|---:|---:|",
    ]
    for row in metrics["pairs"]:
        wall = " / ".join(
            _fmt(row[arm]["performance"]["wall_seconds"]["median"] / 3600, 2)
            for arm in ("candidate", "reference")
        )
        rss = " / ".join(
            _fmt(row[arm]["performance"]["process_group_peak_rss_mb"]["median"] / 1024, 2)
            for arm in ("candidate", "reference")
        )
        rows.append(
            f"| {SHORT_LABELS[row['id']]} | {wall} | "
            f"{_fmt(row['performance_ratios']['wall']['median'], 3)}× | {rss} | "
            f"{_fmt(row['performance_ratios']['process_group_rss']['median'], 3)}× | "
            f"{row['candidate']['valid_repetitions']}/4 | "
            f"{'yes' if row['candidate']['determinism']['passed'] else 'no'} |"
        )
    aggregate = metrics["aggregate"]
    wall = aggregate["performance_ratios"]["wall"]
    rss = aggregate["performance_ratios"]["process_group_rss"]
    rows.extend([
        "",
        "| Cohort performance summary | Observed |",
        "|---|---:|",
        (
            "| Equal-transfer wall ratio (95% CI) | "
            f"{wall['estimate']:.3f}× "
            f"({wall['low']:.3f}–{wall['high']:.3f}) |"
        ),
        (
            "| Summed median-duration ratio | "
            f"{aggregate['total_work']['wall_ratio']:.3f}× |"
        ),
        (
            "| Equal-transfer process-group RSS ratio (95% CI) | "
            f"{rss['estimate']:.3f}× "
            f"({rss['low']:.3f}–{rss['high']:.3f}) |"
        ),
        (
            "| Two-concurrent-transfer v1.0.11 memory proxy | "
            f"{aggregate['candidate_concurrent_memory_proxy_gib']:.1f} GiB |"
        ),
    ])
    return "\n".join(rows)


def table_qualification(metrics: Mapping[str, Any], _sensitivity: Mapping[str, Any] | None) -> str:
    rows = [
        "| Prespecified cohort gate | Result | Observed | Required |",
        "|---|---:|---:|---:|",
    ]
    for gate in metrics["qualification"]["gates"]:
        rows.append(
            f"| {gate['name']} | {'PASS' if gate['passed'] else '**FAIL**'} | "
            f"{_fmt(gate['observed'], 6)} | {_fmt(gate['threshold'], 6)} |"
        )
    return "\n".join(rows)


def table_sensitivity(metrics: Mapping[str, Any], sensitivity: Mapping[str, Any] | None) -> str:
    if sensitivity is None:
        raise AssetError("biology-sensitivity requires annotation sensitivity")
    rows = [
        "| Transfer | Protein-RBH gene F1 at overlap 0.25 / 0.50 / 0.80 (v1.0.11) | v1.0.8 at 0.50 | NCBI Gene groups | Provider-defined gene F1 v1.0.11 / v1.0.8 |",
        "|---|---:|---:|---:|---:|",
    ]
    by_id = {row["id"]: row for row in sensitivity["pairs"]}
    for metric_row in metrics["pairs"]:
        row = by_id[metric_row["id"]]
        primary = row["primary_protein_rbh"]
        candidate_values = " / ".join(
            _fmt(primary[str(threshold)]["candidate"]["gene"]["locus"]["f1"])
            for threshold in (0.25, 0.5, 0.8)
        )
        reference = _fmt(primary["0.5"]["reference"]["gene"]["locus"]["f1"])
        provider = row["provider_gene_relationships"]
        provider_value = "—"
        if provider["available"]:
            provider_value = " / ".join(
                _fmt(provider["arms"][arm]["gene"]["locus"]["f1"])
                for arm in ("candidate", "reference")
            )
        rows.append(
            f"| {SHORT_LABELS[row['id']]} | {candidate_values} | {reference} | "
            f"{provider['groups']:,} | {provider_value} |"
        )
    return "\n".join(rows)


def table_provenance(metrics: Mapping[str, Any], sensitivity: Mapping[str, Any] | None) -> str:
    rows = [
        "| Evidence object | SHA-256 / fingerprint |",
        "|---|---|",
        f"| Exact study manifest | `{metrics['study']['sha256']}` |",
        f"| Annotation preflight | `{metrics['preflight']['sha256']}` |",
        f"| Fast/safe canary | `{metrics['canary']['sha256']}` |",
        f"| Reduced study metrics | `{metrics['fingerprint']}` |",
    ]
    if metrics.get("provider_ortholog_lock"):
        rows.append(
            f"| NCBI Gene relationship lock | "
            f"`{metrics['provider_ortholog_lock']['sha256']}` |"
        )
    if sensitivity is not None:
        rows.append(
            f"| Released-annotation sensitivity | `{sensitivity['fingerprint']}` |"
        )
    return "\n".join(rows)


BLOCKS: dict[str, Callable[[Mapping[str, Any], Mapping[str, Any] | None], str]] = {
    "biology-headline": headline_summary,
    "biology-cohort": table_cohort,
    "biology-comprehensiveness": table_comprehensiveness,
    "biology-accuracy": table_accuracy,
    "biology-performance": table_performance,
    "biology-qualification": table_qualification,
    "biology-sensitivity": table_sensitivity,
    "biology-provenance": table_provenance,
}


def render_block(
    name: str,
    metrics: Mapping[str, Any],
    sensitivity: Mapping[str, Any] | None,
) -> str:
    if name not in BLOCKS:
        raise AssetError(f"unknown biology report block {name!r}")
    return BLOCKS[name](metrics, sensitivity).rstrip()


def update_report(
    text: str,
    metrics: Mapping[str, Any],
    sensitivity: Mapping[str, Any] | None,
) -> str:
    matches = list(MARKER.finditer(text))
    names = [match.group("name") for match in matches]
    if set(names) != set(BLOCKS) or len(names) != len(BLOCKS):
        raise AssetError(
            "report biology markers must occur exactly once: "
            f"expected={sorted(BLOCKS)}, observed={sorted(names)}"
        )
    updated = text
    for name in names:
        start = START.format(name=name)
        end = END.format(name=name)
        start_index = updated.find(start)
        end_index = updated.find(end, start_index + len(start))
        if end_index < 0:
            raise AssetError(f"generated block {name!r} has no end marker")
        rendered = "\n".join((start, render_block(name, metrics, sensitivity), end))
        updated = updated[:start_index] + rendered + updated[end_index + len(end):]
    return updated


def check_report(
    path: Path,
    metrics: Mapping[str, Any],
    sensitivity: Mapping[str, Any] | None,
) -> tuple[bool, str]:
    observed = path.read_text(encoding="utf-8")
    expected = update_report(observed, metrics, sensitivity)
    if observed == expected:
        return True, ""
    return False, "".join(difflib.unified_diff(
        observed.splitlines(keepends=True),
        expected.splitlines(keepends=True),
        fromfile=str(path),
        tofile=f"{path} (generated)",
    ))


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metrics", type=Path, required=True)
    parser.add_argument("--sensitivity", type=Path)
    operation = parser.add_mutually_exclusive_group(required=True)
    operation.add_argument("--figures", type=Path, metavar="DIR")
    operation.add_argument("--update", type=Path, metavar="REPORT")
    operation.add_argument("--check", type=Path, metavar="REPORT")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    arguments = build_parser().parse_args(argv)
    try:
        metrics = load_metrics(arguments.metrics)
        sensitivity = (
            load_sensitivity(arguments.sensitivity, metrics)
            if arguments.sensitivity else None
        )
        if arguments.figures:
            for path in generate_figures(metrics, arguments.figures):
                print(path)
            return 0
        report = arguments.update or arguments.check
        passed, diff = check_report(report, metrics, sensitivity)
        if arguments.update:
            if not passed:
                report.write_text(
                    update_report(report.read_text(encoding="utf-8"), metrics, sensitivity),
                    encoding="utf-8",
                )
                print(f"updated {report}")
            else:
                print(f"up to date: {report}")
            return 0
        if passed:
            print(f"generated tables match: {report}")
            return 0
        os.sys.stderr.write(diff)
        return 1
    except (OSError, TypeError, ValueError, AssetError) as exc:
        print(f"whole-genome-assets: {exc}", file=os.sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
