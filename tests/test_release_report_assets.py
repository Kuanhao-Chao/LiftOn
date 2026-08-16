from copy import deepcopy
import json

import pytest

from benchmarks.compare import release_report_figures, report_figures, report_tables


def _interval(estimate, low=None, high=None):
    return {
        "estimate": estimate,
        "low": estimate if low is None else low,
        "high": estimate if high is None else high,
    }


def _metrics():
    panels = {}
    cells = []
    for index, panel in enumerate(("subset", "full", "e2e"), start=1):
        panels[panel] = {
            "n_cells": 1,
            "modes": {"candidate": "fast", "reference": "fast"},
            "wall_ratio": _interval(0.8 + index * 0.02, 0.78, 0.91),
            "rss_ratio": _interval(0.2 + index * 0.01, 0.18, 0.25),
            "wall_total_ratio": 0.83,
            "wall_ratio_worst": 0.91,
            "rss_ratio_worst": 0.25,
            "candidate_concurrent_memory_envelope_gib": 8.0 * index,
            "candidate_concurrency_limit": 4 if panel == "subset" else 2,
            "quality": {
                "completeness_coding": _interval(0.001, 0.0, 0.002),
                "completeness_feature_total": _interval(0.002, 0.001, 0.003),
                "covpi": _interval(0.003, 0.002, 0.004),
                "intron_chain_exact_recall": _interval(0.001, 0.0, 0.002),
                "orf_valid_recall": _interval(0.001, 0.0, 0.002),
            },
            "stable_id_preservation": {
                "CDS": {"rate_delta": 0.001},
                "exon": {"rate_delta": 0.002},
            },
        }
        cells.append({
            "panel": panel,
            "benchmark": f"demo_{panel}",
            "repetitions": 4,
            "candidate_wall_seconds_median": 8.0,
            "reference_wall_seconds_median": 10.0,
            "wall_ratio_median": 0.8,
            "rss_ratio_median": 0.25,
            "candidate_peak_rss_gib": 2.0,
            "common_pi": {"mean_delta": 0.001},
            "candidate_valid": True,
            "reference_valid": False,
            "candidate_byte_deterministic": True,
            "reference_byte_deterministic": True,
            "candidate_semantic_deterministic": True,
            "reference_semantic_deterministic": True,
            "candidate_quality_deterministic": True,
            "reference_quality_deterministic": True,
            "common_pi_deterministic": True,
            "e2e_biology": (
                {
                    "candidate": {"valid": True, "deterministic": True},
                    "reference": {"valid": True, "deterministic": True},
                }
                if panel == "e2e" else None
            ),
        })
    return {
        "schema_version": 4,
        "campaign": {
            "candidate_sha": "c" * 40,
            "reference_sha": "r" * 40,
            "n_pairs": 12,
            "n_cells": 3,
            "matrix_complete": True,
            "expected_campaign": {
                "profile_id": "canonical-v1",
                "profile_digest": "d" * 64,
            },
        },
        "bootstrap": {"seed": 17, "replicates": 1000},
        "publication": {"mode": "release"},
        "panels": panels,
        "cells": cells,
        "verdict": {"passed": True},
    }


def _exact_metrics():
    """Small-valued but matrix-complete analogue of the sealed diagnostic."""

    base = _metrics()
    panel_ids = {
        "subset": [f"subset_{index:02d}" for index in range(34)],
        "full": [
            "bee", "drosophila", "arabidopsis_to_rice", "human_to_zebrafish",
            "t2_human_to_gorilla", "t2_mouse_to_caroli", "t3_dog_to_cat",
            "t3_human_to_macaque", "t3_human_to_marmoset",
            "t4_drosophila_to_bee", "t4_human_to_chicken", "t4_human_to_xenopus",
        ],
        "e2e": ["bee", "human", "mouse"],
    }
    missing_ids = {
        "full": [
            "arabidopsis", "rice", "t1_maize_b73_to_mo17",
            "t1_tomato_microtom_to_heinz", "t2_tomato_to_potato",
        ],
        "e2e": ["arabidopsis", "rice"],
    }
    e2e_values = {
        "bee": (23471, 23365, 23364, 0.99548, 0.99544, 0.94570, 0.94493,
                0.99051, 0.99024, 0.986034, 0.985728),
        "human": (144329, 129408, 129795, 0.89662, 0.89930, 0.40655, 0.37493,
                  0.99679, 0.99592, 0.893704, 0.895538),
        "mouse": (96192, 95430, 95546, 0.99208, 0.99328, 0.87002, 0.82091,
                  0.99506, 0.99463, 0.986781, 0.987551),
    }
    template = base["cells"][0]
    cells = []
    for panel, benchmarks in panel_ids.items():
        for index, benchmark in enumerate(benchmarks):
            cell = deepcopy(template)
            cell.update({
                "panel": panel,
                "benchmark": benchmark,
                "repetitions": 4,
                "candidate_valid": not (panel == "full" and benchmark == "arabidopsis_to_rice"),
                "quality_deltas": {
                    "completeness_coding": -0.0001,
                    "completeness_feature_total": 0.002,
                    "covpi": 0.001,
                    "recall_at_0.5": 0.001,
                    "recall_at_0.75": 0.001,
                    "recall_at_0.9": 0.001,
                    "recall_at_0.95": 0.001,
                    "intron_chain_exact_recall": 0.001,
                    "orf_valid_recall": 0.001,
                },
                "wall_ratio_median": 0.75 + index / 1000,
                "rss_ratio_median": 0.20 + index / 1000,
                "candidate_peak_rss_gib": 2.0 + index / 100,
                "common_pi": {"mean_delta": 0.001, "n_common": 100},
                "candidate": {
                    "n_reference_coding": 100,
                    "n_recovered_coding": 99,
                    "completeness_coding": 0.99,
                    "completeness_feature_total": 0.95,
                    "mean_pi": 0.98,
                    "covpi": 0.97,
                    "intron_chain_exact_recall": 0.96,
                    "orf_valid_recall": 0.97,
                },
                "reference": {
                    "n_reference_coding": 100,
                    "n_recovered_coding": 98,
                    "completeness_coding": 0.98,
                    "completeness_feature_total": 0.94,
                    "mean_pi": 0.97,
                    "covpi": 0.96,
                    "intron_chain_exact_recall": 0.95,
                    "orf_valid_recall": 0.96,
                },
            })
            if panel == "e2e":
                (denominator, candidate_recovered, reference_recovered,
                 candidate_coding, reference_coding, candidate_feature,
                 reference_feature, candidate_pi, reference_pi,
                 candidate_covpi, reference_covpi) = e2e_values[benchmark]
                cell["candidate"].update({
                    "n_reference_coding": denominator,
                    "n_recovered_coding": candidate_recovered,
                    "completeness_coding": candidate_coding,
                    "completeness_feature_total": candidate_feature,
                    "mean_pi": candidate_pi,
                    "covpi": candidate_covpi,
                })
                cell["reference"].update({
                    "n_reference_coding": denominator,
                    "n_recovered_coding": reference_recovered,
                    "completeness_coding": reference_coding,
                    "completeness_feature_total": reference_feature,
                    "mean_pi": reference_pi,
                    "covpi": reference_covpi,
                })
                cell["absolute_quality"] = {
                    "candidate": dict(cell["candidate"]),
                    "reference": dict(cell["reference"]),
                }
            else:
                cell["e2e_biology"] = None
            cells.append(cell)
    planned_panels = {
        panel: {"ids": ids + missing_ids.get(panel, []), "repetitions": 4}
        for panel, ids in panel_ids.items()
    }
    matrix = [
        {
            "panel": panel,
            "candidate_mode": "fast",
            "reference_mode": "fast",
            "threads": 8,
        }
        for panel in ("subset", "full", "e2e")
    ]
    base["cells"] = cells
    base["campaign"].update({
        "candidate_sha": report_tables.CANDIDATE_SHA,
        "reference_sha": report_tables.REFERENCE_SHA,
        "n_pairs": 196,
        "n_cells": 49,
        "matrix_complete": False,
        "planned_campaign": {"profile_id": "canonical-v1", "panels": planned_panels, "matrix": matrix},
        "missing_keys": [
            [panel, benchmark, repetition]
            for panel, benchmarks in missing_ids.items()
            for benchmark in benchmarks
            for repetition in range(1, 5)
        ],
    })
    base["publication"] = {
        "mode": "diagnostic",
        "evaluation_overlays": {"validated": True, "pair_count": 196},
        "controller_evidence": {"roots": [{"root": "/sealed/root"}]},
    }
    base["bootstrap"] = {"seed": 20260717, "replicates": 10000}
    base["quality_baseline_artifact"] = {"path": "/sealed/fourway.json", "sha256": "a" * 64}
    base["environment"] = {"finalization_host": {"hostname": "sealed-host"}}
    for panel, summary in base["panels"].items():
        summary["n_cells"] = len(panel_ids[panel])
        summary["candidate_wall_seconds_total_of_cell_medians"] = 100.0
        summary["reference_wall_seconds_total_of_cell_medians"] = 150.0
        summary["stable_id_preservation"] = {
            feature: {
                "n_applicable_cells": len(panel_ids[panel]),
                "n_reference_ids": 1000,
                "candidate_n_preserved_ids": 900,
                "reference_n_preserved_ids": 700,
                "candidate_preservation_rate": 0.9,
                "reference_preservation_rate": 0.7,
                "rate_delta": 0.2,
            }
            for feature in ("CDS", "exon")
        }
        summary["quality"].update({
            name: _interval(0.001, 0.0005, 0.0015)
            for name in (
                "recall_at_0.5", "recall_at_0.75", "recall_at_0.9",
                "recall_at_0.95",
            )
        })
    base["performance_outliers"] = [{"panel": "full", "benchmark": "arabidopsis_to_rice"}]
    base["verdict"] = {
        "passed": False,
        "checks": [
            {"name": "complete_immutable_release_campaign", "passed": False,
             "actual": "diagnostic", "limit": "release mode with complete controller-backed evidence"},
            {"name": "all_candidate_outputs_valid", "passed": False, "actual": 48, "limit": 49},
        ],
    }
    return base


def _failure_audit(metrics):
    roots = ["/sealed/root"]
    cells = []
    for panel, benchmark, repetition in metrics["campaign"]["missing_keys"]:
        cells.append({
            "root": roots[0],
            "cell_id": f"{panel}__{benchmark}__{repetition:02d}",
            "panel": panel,
            "benchmark": benchmark,
            "repetition": repetition,
            "returncode": 1,
            "elapsed_seconds": 100.0 + repetition,
            "watchdog_reason": None,
            "feature_not_found": True,
            "exception_string_typeerror": True,
            "candidate_completed_before_control_failure": repetition % 2 == 0,
        })
    counts = {}
    for cell in cells:
        counts[(cell["panel"], cell["benchmark"])] = counts.get(
            (cell["panel"], cell["benchmark"]), 0
        ) + 1
    return {
        "schema_version": 1,
        "classification": "deterministic historical-reference failure",
        "selected_roots": roots,
        "failed_pair_count": 28,
        "failed_case_counts": [
            {"panel": panel, "benchmark": benchmark, "failed_repetitions": count}
            for (panel, benchmark), count in sorted(counts.items())
        ],
        "invariants": {
            "all_returncode_1": True,
            "all_watchdog_reason_null": True,
            "all_feature_not_found": True,
            "all_exception_string_typeerror": True,
            "candidate_first_outputs_preserved": 14,
        },
        "cells": cells,
    }


def _manifest(metrics):
    return {
        "schema_version": 4,
        "candidate_sha": metrics["campaign"]["candidate_sha"],
        "reference_sha": metrics["campaign"]["reference_sha"],
        "publication_mode": metrics["publication"]["mode"],
        "evaluation_overlays": metrics["publication"]["evaluation_overlays"],
        "finalization_host": metrics["environment"]["finalization_host"],
        "quality_baseline_artifact": metrics["quality_baseline_artifact"],
        "bootstrap": metrics["bootstrap"],
        "run_roots": ["/sealed/root"],
        "pair_results": [
            {"panel": panel, "benchmark": benchmark, "repetition": repetition}
            for panel, benchmark, repetition in sorted(
                report_tables._observed_pair_keys(metrics)
            )
        ],
    }


def test_release_blocks_render_exact_provenance_and_performance():
    metrics = _metrics()
    metrics["environment"] = {
        "finalization_host": {
            "hostname": "report-host",
            "cpu_model": "Example CPU",
            "sockets": 2,
            "logical_cpus": 64,
            "memory_gib": 512.0,
            "kernel": "example-kernel",
        },
    }
    metrics["publication"]["controller_evidence"] = {
        "roots": [{
            "plan": {
                "execution_environment": {
                    "python": "3.11.15",
                    "tools": {
                        "minimap2": "2.28-r1209",
                        "miniprot": "0.13-r248",
                    },
                    "packages": {
                        "gffutils": "0.14",
                        "parasail": "1.3.4",
                    },
                },
            },
        }],
    }
    labels = report_tables.Labels()

    provenance = report_tables.render_block(
        "release-provenance",
        metrics=metrics,
        fourway=None,
        version_compare=None,
        labels=labels,
    )
    performance = report_tables.render_block(
        "release-performance",
        metrics=metrics,
        fourway=None,
        version_compare=None,
        labels=labels,
    )

    assert "LiftOn v1.0.11" in provenance
    assert "`canonical-v1`" in provenance
    assert "**PASS**" in provenance
    assert "12 AB/BA pairs across 3 benchmark cells" in provenance
    assert "Frozen execution software" in provenance
    assert "Finalization host *(not a historical run record)*" in provenance
    assert "| subset | fast/fast | 1 |" in performance
    assert "descriptive; n=1 cell" in performance
    assert "24.00 GiB (2 slots)" in performance

    correctness = report_tables.render_block(
        "release-correctness",
        metrics=metrics,
        fourway=None,
        version_compare=None,
        labels=labels,
    )
    assert "| E2E biological summary valid | 1/1 | 1/1 |" in correctness
    assert "| subset | `CDS` | 0 | — | — | +0.00100 |" in correctness


def test_legacy_recovery_table_only_lists_documented_failures():
    table = report_tables.table3(
        report_tables.load(report_tables.FOURWAY),
        {},
        report_tables.Labels(),
    )

    rows = [line for line in table.splitlines()[2:] if line.startswith("|")]
    assert len(rows) == 5
    assert [line.split("|")[1].strip() for line in rows] == [
        "arabidopsis",
        "rice",
        "maize",
        "tomato",
        "tomato→potato",
    ]
    assert "+34,589" in table
    assert "+9,594" in table
    assert table.count("crash (no output)") == 6
    assert "+-" not in table


def test_generated_block_check_detects_and_repairs_drift(tmp_path):
    metrics = _metrics()
    labels = report_tables.Labels()
    report = tmp_path / "report.mdx"
    report.write_text(
        "before\n"
        "<!-- BEGIN GENERATED: release-performance -->\n"
        "stale\n"
        "<!-- END GENERATED: release-performance -->\n"
        "after\n"
    )

    passed, diff = report_tables.check_generated_blocks(
        report,
        metrics=metrics,
        fourway={},
        version_compare={},
        labels=labels,
    )
    assert passed is False
    assert "-stale" in diff

    report.write_text(report_tables.update_generated_blocks(
        report.read_text(),
        metrics=metrics,
        fourway={},
        version_compare={},
        labels=labels,
    ))
    passed, diff = report_tables.check_generated_blocks(
        report,
        metrics=metrics,
        fourway={},
        version_compare={},
        labels=labels,
    )
    assert passed is True
    assert diff == ""


def test_report_placeholders_are_derived_from_metrics():
    metrics = _metrics()
    metrics["campaign"]["planned_campaign"] = {
        "panels": {
            "subset": {"ids": ["demo_subset"], "repetitions": 4},
            "full": {"ids": ["demo_full"], "repetitions": 4},
            "e2e": {"ids": ["demo_e2e"], "repetitions": 4},
        },
    }
    metrics["campaign"]["missing_keys"] = [["full", "demo_full", 4]]
    metrics["publication"] = {"mode": "diagnostic"}

    report = report_tables.update_report_placeholders(
        "@@FINAL_ABSTRACT_RESULTS@@\n@@FINAL_RELEASE_VERDICT@@",
        metrics,
    )

    assert "@@" not in report
    assert "12 of 12 planned paired executions" in report
    assert "subset 0.820× wall/0.210× RSS" in report
    assert "pass current validation in 3/3 cells" in report
    assert "DIAGNOSTIC ONLY" in report


def test_report_tables_support_mdx_generated_markers(tmp_path):
    metrics = _metrics()
    report = tmp_path / "report.mdx"
    report.write_text(
        "{/* BEGIN GENERATED: release-performance */}\n"
        "stale\n"
        "{/* END GENERATED: release-performance */}\n"
    )

    common = {
        "metrics": metrics,
        "fourway": {},
        "version_compare": {},
        "labels": report_tables.Labels(),
    }
    updated = report_tables.update_generated_blocks(report.read_text(), **common)
    report.write_text(updated)

    passed, diff = report_tables.check_generated_blocks(report, **common)
    assert passed is True
    assert diff == ""
    assert "{/* BEGIN GENERATED: release-performance */}" in updated
    assert "<!-- BEGIN GENERATED" not in updated


def test_report_tables_cli_updates_then_checks(tmp_path):
    metrics_path = tmp_path / "metrics.json"
    metrics_path.write_text(json.dumps(_metrics()))
    report = tmp_path / "report.mdx"
    report.write_text(
        "<!-- BEGIN GENERATED: release-provenance -->\n"
        "stale\n"
        "<!-- END GENERATED: release-provenance -->\n"
    )

    common = ["--release-metrics", str(metrics_path)]
    assert report_tables.main(["--update", str(report), *common]) == 0
    assert report_tables.main(["--check", str(report), *common]) == 0
    assert "stale" not in report.read_text()


def test_generated_blocks_reject_duplicate_or_unclosed_markers():
    metrics = _metrics()
    labels = report_tables.Labels()
    duplicate = (
        "<!-- BEGIN GENERATED: release-performance -->\n"
        "x\n<!-- END GENERATED: release-performance -->\n"
        "<!-- BEGIN GENERATED: release-performance -->\n"
        "y\n<!-- END GENERATED: release-performance -->\n"
    )
    with pytest.raises(ValueError, match="duplicate generated blocks"):
        report_tables.update_generated_blocks(
            duplicate,
            metrics=metrics,
            fourway={},
            version_compare={},
            labels=labels,
        )
    with pytest.raises(ValueError, match="no end marker"):
        report_tables.update_generated_blocks(
            "<!-- BEGIN GENERATED: release-performance -->\n",
            metrics=metrics,
            fourway={},
            version_compare={},
            labels=labels,
        )


def test_release_report_figures_write_nonempty_pngs(tmp_path):
    paths = release_report_figures.generate_figures(_metrics(), tmp_path)

    assert {path.name for path in paths} == {
        release_report_figures.PERFORMANCE_FIGURE,
        release_report_figures.CORRECTNESS_FIGURE,
    }
    assert all(path.stat().st_size > 10_000 for path in paths)


def test_release_report_figures_label_diagnostic_metrics():
    metrics = _metrics()
    metrics["publication"] = {"mode": "diagnostic"}

    assert release_report_figures._figure_context(metrics) == (
        "paired benchmark diagnostic", "DIAGNOSTIC",
    )


def test_release_report_figures_reject_wrong_schema():
    with pytest.raises(ValueError, match="schema-version 4"):
        release_report_figures.validate_metrics({"schema_version": 3})


def test_legacy_figure_labels_are_configurable():
    arguments = report_figures.build_parser().parse_args([
        "/tmp/report-figures",
        "--candidate-label", "candidate SHA",
        "--reference-label", "reference SHA",
    ])

    assert arguments.candidate_label == "candidate SHA"
    assert arguments.reference_label == "reference SHA"


def test_failure_audit_and_manifest_cross_validate_exact_matrix():
    metrics = _exact_metrics()
    audit = _failure_audit(metrics)
    manifest = _manifest(metrics)

    assert report_tables.validate_failure_audit(audit, metrics) is audit
    assert report_tables.validate_release_manifest(manifest, metrics) is manifest
    assert release_report_figures.validate_failure_audit(audit, metrics) is audit
    assert release_report_figures.validate_manifest(manifest, metrics) is manifest


@pytest.mark.parametrize(
    ("mutation", "message"),
    [
        (lambda metrics, audit: audit.update(schema_version=2), "schema-version 1"),
        (
            lambda metrics, audit: audit["cells"][0].update(returncode=137),
            "not a deterministic v1.0.8 failure",
        ),
        (
            lambda metrics, audit: metrics["campaign"]["missing_keys"].pop(),
            "disagree with release metrics missing_keys",
        ),
        (
            lambda metrics, audit: audit["invariants"].update(
                candidate_first_outputs_preserved=13
            ),
            "candidate-first",
        ),
    ],
)
def test_failure_audit_rejects_malformed_or_mismatched_documents(mutation, message):
    metrics = _exact_metrics()
    audit = _failure_audit(metrics)
    mutation(metrics, audit)

    with pytest.raises(ValueError, match=message):
        report_tables.validate_failure_audit(audit, metrics)
    with pytest.raises(ValueError, match=message):
        release_report_figures.validate_failure_audit(audit, metrics)


def test_failure_audit_rejects_candidate_or_reference_sha_drift():
    metrics = _exact_metrics()
    audit = _failure_audit(metrics)
    metrics["campaign"]["candidate_sha"] = "f" * 40

    with pytest.raises(ValueError, match="exact v1.0.11 candidate SHA"):
        report_tables.validate_failure_audit(audit, metrics)
    with pytest.raises(ValueError, match="exact v1.0.11 candidate SHA"):
        release_report_figures.validate_failure_audit(audit, metrics)


def test_release_manifest_rejects_sha_and_pair_inventory_drift():
    metrics = _exact_metrics()
    manifest = _manifest(metrics)
    manifest["reference_sha"] = "f" * 40
    with pytest.raises(ValueError, match="reference SHA"):
        report_tables.validate_release_manifest(manifest, metrics)

    manifest = _manifest(metrics)
    manifest["pair_results"].pop()
    with pytest.raises(ValueError, match="pair results"):
        release_report_figures.validate_manifest(manifest, metrics)


def test_failure_audit_blocks_are_generated_and_drift_checked(tmp_path):
    metrics = _exact_metrics()
    audit = _failure_audit(metrics)
    manifest = _manifest(metrics)
    report = tmp_path / "report.mdx"
    report.write_text(
        "{/* BEGIN GENERATED: release-coverage */}\n"
        "stale\n"
        "{/* END GENERATED: release-coverage */}\n\n"
        "{/* BEGIN GENERATED: release-e2e-outcomes */}\n"
        "stale\n"
        "{/* END GENERATED: release-e2e-outcomes */}\n"
    )
    common = {
        "metrics": metrics,
        "fourway": {},
        "version_compare": {},
        "labels": report_tables.Labels(),
        "failure_audit": audit,
        "manifest": manifest,
    }

    updated = report_tables.update_generated_blocks(report.read_text(), **common)
    assert "196" in updated
    assert "23,365 / 23,364 (23,471)" in updated
    report.write_text(updated)
    passed, diff = report_tables.check_generated_blocks(report, **common)
    assert passed is True
    assert diff == ""

    report.write_text(updated.replace("23,365 / 23,364", "23,366 / 23,364"))
    passed, diff = report_tables.check_generated_blocks(report, **common)
    assert passed is False
    assert "23,366" in diff


def test_exact_release_figure_inventory_and_historical_separation(tmp_path):
    metrics = _exact_metrics()
    paths = release_report_figures.generate_figures(
        metrics,
        tmp_path,
        failure_audit=_failure_audit(metrics),
        manifest=_manifest(metrics),
        legacy_results=report_tables.load(report_tables.FOURWAY),
    )

    assert {path.name for path in paths} == (
        release_report_figures.EXACT_FIGURE_INVENTORY
        | {release_report_figures.HISTORICAL_FIGURE}
    )
    assert all(path.stat().st_size > 10_000 for path in paths)


def test_failure_derived_block_requires_explicit_audit():
    with pytest.raises(ValueError, match="requires --failure-audit"):
        report_tables.render_block(
            "release-coverage",
            metrics=_exact_metrics(),
            fourway={},
            version_compare={},
            labels=report_tables.Labels(),
        )
