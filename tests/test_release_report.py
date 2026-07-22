from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path

import pytest

from benchmarks.compare import build_controller, release_report


_DEFAULT = object()
CANDIDATE_SHA = "c" * 40
REFERENCE_SHA = "b" * 40
GIT_STATE = {
    "head": "d" * 40,
    "tracked_diff_sha256": hashlib.sha256(b"").hexdigest(),
    "dirty_tracked_paths": [],
}


TSV_HEADER = (
    "ref_mrna_id\ttool_feature_id\trecovered\tis_coding\tcopy_index\t"
    "protein_identity\tdna_identity\tref_prot_len\tlifted_prot_len\t"
    "ref_dna_len\tlifted_dna_len\tdna_basis\tstatus\tseqid\t"
    "n_cds_lifted\tn_cds_ref\tintron_chain_exact\texon_sn\texon_sp\t"
    "orf_start_ok\torf_stop_ok\torf_valid\n"
)


def _sha256(value: str) -> str:
    return hashlib.sha256(value.encode()).hexdigest()


def _source_file_record(path: Path) -> dict:
    path = Path(path).resolve()
    return {
        "path": str(path),
        "size": path.stat().st_size,
        "sha256": release_report.sha256_file(path),
    }


@pytest.fixture(autouse=True)
def _stable_reporter_git_state(monkeypatch):
    monkeypatch.setattr(
        release_report,
        "_current_reporter_git_state",
        lambda: dict(GIT_STATE),
    )


def _tsv_row(
    reference_id: str,
    *,
    recovered: str = "1",
    is_coding: str = "1",
    protein_identity: str = "1.0",
    dna_identity: str = "1.0",
    status: str = "recovered",
) -> str:
    return (
        f"{reference_id}\ttool-{reference_id}\t{recovered}\t{is_coding}\t0\t"
        f"{protein_identity}\t{dna_identity}\t10\t10\t30\t30\ttranscript\t"
        f"{status}\tchr1\t1\t1\t1\t1.0\t1.0\t1\t1\t1\n"
    )


def _write_tsv(path: Path, identities: list[float]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    rows = [TSV_HEADER]
    for index, identity in enumerate(identities, 1):
        rows.append(_tsv_row(f"tx{index}", protein_identity=str(identity)))
    path.write_text("".join(rows))


def _write_pair(
    root: Path,
    *,
    repetition: int,
    candidate_wall: float,
    reference_wall: float,
    panel: str = "subset",
    benchmark: str = "demo",
    candidate_rss_mb: float = 100.0,
    reference_rss_mb: float = 200.0,
    candidate_biology=_DEFAULT,
    reference_biology=_DEFAULT,
    candidate_sha: str = CANDIDATE_SHA,
    reference_sha: str = REFERENCE_SHA,
    threads: int = 8,
    input_tag: str = "shared-input",
    candidate_mode: str | None = None,
    reference_mode: str | None = None,
    order=None,
    candidate_valid: bool = True,
    reference_valid: bool = True,
    candidate_byte_sha: str = "candidate-bytes",
    reference_byte_sha: str = "reference-bytes",
) -> Path:
    cell = root / f"{panel}-{benchmark}" / f"rep-{repetition}"
    evaluation = cell / "evaluation"
    _write_tsv(evaluation / "candidate.transcripts.tsv", [1.0, 1.0])
    _write_tsv(evaluation / "reference.transcripts.tsv", [1.0, 1.0])
    summary = {
        "n_reference_coding": 2,
        "n_recovered_coding": 2,
        "completeness_coding": 1.0,
        "completeness_feature_total": 1.0,
        "protein_identity": {"n": 2, "mean": 1.0, "median": 1.0},
        "stable_id_preservation": {
            "method": release_report.STABLE_ID_METHOD,
            "by_type": {
                feature_type: {
                    "applicable": True,
                    "reason": None,
                    "n_reference_records": 2,
                    "n_reference_records_with_id": 2,
                    "n_reference_ids": 2,
                    "n_preserved_ids": 2,
                    "n_output_records": 2,
                    "n_output_records_with_id": 2,
                    "n_output_ids": 2,
                    "preservation_rate": 1.0,
                }
                for feature_type in release_report.STABLE_ID_FEATURE_TYPES
            },
        },
        "completeness_by_type": {
            "gene": {
                "n_reference": 2,
                "n_recovered": 2,
                "pct_recovered": 1.0,
                "n_extra_copies": 0,
                "n_tool_features": 2,
            },
            "_overall_": {
                "n_reference": 2,
                "n_recovered": 2,
                "pct_recovered": 1.0,
                "n_extra_copies": 0,
                "n_tool_features": 2,
            },
        },
    }
    payload = {
        "schema_version": release_report.SCHEMA_VERSION,
        "panel": panel,
        "benchmark": benchmark,
        "repetition": repetition,
        "order": order if order is not None else (
            ["reference", "candidate"]
            if repetition % 2 else ["candidate", "reference"]
        ),
        "threads": threads,
        "modes": {
            "candidate": candidate_mode or panel,
            "reference": reference_mode or panel,
        },
        "provenance": {
            "candidate": {"sha": candidate_sha},
            "reference": {"sha": reference_sha},
        },
        "inputs": {
            "ref_fa": {
                "path": "/immutable-inputs/ref.fa",
                "sha256": _sha256(input_tag),
                "size": 1,
            }
        },
        "evaluation_inputs": {},
        "versions": {
            "candidate": {
                "profile": {
                    "wall_clock_seconds": candidate_wall,
                    "peak_rss_mb": candidate_rss_mb,
                    "user_cpu_seconds": 1.0,
                    "sys_cpu_seconds": 0.5,
                    "exit_code": 0,
                },
                "summary": json.loads(json.dumps(summary)),
                "validation": {
                    "is_valid": candidate_valid,
                    "n_errors": 0 if candidate_valid else 1,
                    "issues": [],
                },
                "fingerprints": {
                    "byte_sha256": _sha256(candidate_byte_sha),
                    "semantic_sha256": _sha256("candidate-semantic"),
                },
            },
            "reference": {
                "profile": {
                    "wall_clock_seconds": reference_wall,
                    "peak_rss_mb": reference_rss_mb,
                    "user_cpu_seconds": 1.0,
                    "sys_cpu_seconds": 0.5,
                    "exit_code": 0,
                },
                "summary": json.loads(json.dumps(summary)),
                "validation": {
                    "is_valid": reference_valid,
                    "n_errors": 0 if reference_valid else 1,
                    "issues": [],
                },
                "fingerprints": {
                    "byte_sha256": _sha256(reference_byte_sha),
                    "semantic_sha256": _sha256("reference-semantic"),
                },
            },
        },
        "ratios": {
            "wall": candidate_wall / reference_wall,
            "peak_rss": candidate_rss_mb / reference_rss_mb,
        },
    }
    if panel == "e2e":
        default_biology = {
            "schema_version": 1,
            "reference_features": 10,
            "mapped_features": 9,
            "lost_features": 1,
            "extra_copy_features": 1,
            "feature_completeness": 0.9,
            "emitted_transcript_records": 12,
            "mapped_transcripts_reported": 11,
            "evaluated_transcript_records": 12,
            "evaluated_coding_records": 10,
            "mean_protein_identity": 0.95,
        }
        payload["versions"]["candidate"]["biological_summary"] = (
            default_biology
            if candidate_biology is _DEFAULT else candidate_biology
        )
        payload["versions"]["reference"]["biological_summary"] = (
            default_biology
            if reference_biology is _DEFAULT else reference_biology
        )
        for version in payload["versions"].values():
            version["evaluation_profile"] = {
                "wall_clock_seconds": 1.0,
                "peak_rss_mb": 10.0,
                "user_cpu_seconds": 0.5,
                "sys_cpu_seconds": 0.25,
                "exit_code": 0,
            }
            version["score_summary"] = {
                "format": "lifton_score_v1",
                "records": 12,
                "coding_records": 10,
                "noncoding_records": 2,
                "malformed_records": 0,
            }
            version["evaluation_summary"] = {
                "format": "lifton_eval_v1",
                "records": 12,
                "coding_records": 10,
                "noncoding_records": 2,
                "malformed_records": 0,
            }
    for label in ("candidate", "reference"):
        transcript = (
            evaluation / f"{label}.transcripts.tsv"
        ).resolve()
        payload["versions"][label]["summary"]["transcripts_tsv"] = str(
            transcript
        )
        payload["versions"][label]["evaluation_artifacts"] = {
            "transcripts_tsv": {
                "path": str(transcript),
                "size": transcript.stat().st_size,
                "sha256": release_report.sha256_file(transcript),
            },
        }
    path = cell / "pair_result.json"
    path.write_text(json.dumps(payload))
    return path


def _target_truth_document(
    pair_path: Path,
    label: str,
    *,
    mapping_sha: str,
) -> dict:
    level_scope = {
        "groups": 1,
        "predicted_scored": 1,
        "expected_scored": 1,
        "prediction_models_total": 1,
        "truth_models_total": 1,
        "prediction_models_ignored": 0,
        "truth_models_ignored": 0,
    }
    score = {"precision": 1.0, "recall": 1.0, "f1": 1.0}
    cell = pair_path.parent.resolve()
    return {
        "schema_version": 1,
        "method": "ortholog-scoped-target-coordinate-v1",
        "inputs": {
            "prediction_gff": str(cell / label / f"{label}.gff3"),
            "truth_gff": str(cell / "score-inputs" / "truth.gff3"),
            "source_gff": str(cell / "score-inputs" / "source.gff3"),
            "mapping": {
                "kind": "ortholog-map",
                "path": str(cell / "score-inputs" / "ortholog-map.json"),
                "entries": 1,
                "sha256": mapping_sha,
            },
        },
        "parameters": {
            "minimum_reciprocal_overlap": 0.5,
            "id_policy": "ortholog-map",
            "mapping_required": True,
            "mapping_requirement_satisfied": True,
            "mapping_source_scope_validated": True,
        },
        "scope": {
            "gene_groups": 1,
            "transcript_groups": 1,
            "mapping_entries": 1,
            "mapping_status_counts": {"retained": 1},
            "gene": dict(level_scope),
            "transcript": dict(level_scope),
        },
        "gene": {
            name: dict(score) for name in ("locus", "strand", "copy")
        },
        "transcript": {
            name: dict(score) for name in ("locus", "strand", "copy")
        },
        "structure": {
            name: dict(score)
            for name in ("intron_chain", "intron", "exon", "CDS")
        },
    }


def _attach_target_truth_evidence(
    pair_path: Path,
    *,
    candidate_mapping_sha: str,
    reference_mapping_sha: str,
) -> None:
    payload = json.loads(pair_path.read_text())
    mapping_shas = {
        "candidate": candidate_mapping_sha,
        "reference": reference_mapping_sha,
    }
    for label in ("candidate", "reference"):
        document = _target_truth_document(
            pair_path, label, mapping_sha=mapping_shas[label],
        )
        artifact = (
            pair_path.parent / "evaluation" / f"{label}.target_truth.json"
        ).resolve()
        artifact.write_text(json.dumps(document, sort_keys=True))
        payload["versions"][label]["summary"]["target_truth"] = document
        payload["versions"][label]["evaluation_artifacts"]["target_truth"] = {
            "path": str(artifact),
            "size": artifact.stat().st_size,
            "sha256": release_report.sha256_file(artifact),
        }
    pair_path.write_text(json.dumps(payload))


def _campaign_spec(repetitions: int = 2) -> dict:
    return release_report.canonical_campaign_spec(repetitions)


def _write_controller_run(
    root: Path,
    panel: str,
    *,
    repetitions: int = 2,
    stage: str | None = None,
    compatibility_tag: str = "shared-toolchain",
    candidate_mode: str | None = None,
    reference_mode: str | None = None,
) -> Path:
    benchmarks = _campaign_spec(repetitions)["panels"][panel]["ids"]
    root = Path(root).resolve()
    benchmark_registry = Path(build_controller.DEFAULT_REGISTRY).resolve()
    dataset_registry = Path(build_controller.DEFAULT_DATASET_REGISTRY).resolve()
    baseline = Path(build_controller.DEFAULT_BASELINE).resolve()
    lifton_executable = Path("/env/bin/lifton").resolve()
    paired = {
        "panel": panel,
        "repetitions": repetitions,
        "lifton_executable": str(lifton_executable),
        "registries": {
            "benchmark": str(benchmark_registry),
            "dataset": str(dataset_registry),
        },
        "candidate": {
            "root": str(Path("/sources/candidate").resolve()),
            "sha": CANDIDATE_SHA,
            "e2e_mode": candidate_mode or "fast",
        },
        "reference": {
            "root": str(Path("/sources/reference").resolve()),
            "sha": REFERENCE_SHA,
            "e2e_mode": reference_mode or "fast",
        },
    }
    source_records = {
        label: {
            "label": label,
            "root": paired[label]["root"],
            "sha": paired[label]["sha"],
            "lifton_executable": str(lifton_executable),
            "imported_package": str(
                Path(paired[label]["root"]) / "lifton" / "__init__.py"
            ),
        }
        for label in ("candidate", "reference")
    }
    cells = []
    pending_success = []
    frozen_inputs = {}
    frozen_evaluation_inputs = {}
    for benchmark in benchmarks:
        for repetition in range(1, repetitions + 1):
            result_path = _write_pair(
                root / "cells" / f"{benchmark}-{repetition}",
                panel=panel,
                benchmark=benchmark,
                repetition=repetition,
                candidate_wall=8.0,
                reference_wall=10.0,
                candidate_mode=candidate_mode,
                reference_mode=reference_mode,
            )
            cell_dir = result_path.parent
            raw = json.loads(result_path.read_text())
            raw["registries"] = dict(paired["registries"])
            raw["provenance"] = json.loads(json.dumps(source_records))
            artifacts = {"result_json": str(result_path)}
            observed = {}
            for label in ("candidate", "reference"):
                gff = cell_dir / label / f"{label}.gff3"
                gff.parent.mkdir()
                gff.write_text(
                    "##gff-version 3\n"
                    f"chr1\tLiftOn\tgene\t1\t2\t.\t+\t.\t"
                    f"ID={label}-{benchmark}\n"
                )
                manifest = cell_dir / label / "release_run_manifest.json"
                manifest.write_text(json.dumps({
                    "schema_version": release_report.SCHEMA_VERSION,
                    "status": "success",
                }) + "\n")
                fingerprints = release_report.gff3_fingerprints(gff)
                observed[label] = fingerprints
                raw["versions"][label]["output_gff"] = str(gff)
                raw["versions"][label]["fingerprints"] = fingerprints
                raw["versions"][label]["release_manifest"] = str(manifest)
                artifacts[f"{label}_gff"] = str(gff)
                artifacts[f"{label}_manifest"] = str(manifest)
            result_path.write_text(json.dumps(raw))
            for label in ("candidate", "reference"):
                validation = cell_dir / f"{label}_gff_validation.json"
                validation.write_text(json.dumps({
                    "schema_version": 1,
                    "is_valid": True,
                    "errors": [],
                }) + "\n")
                artifacts[f"{label}_gff_validation"] = str(validation)
            evaluation_artifacts = {}
            for label in ("candidate", "reference"):
                transcript = (
                    cell_dir / "evaluation" / f"{label}.transcripts.tsv"
                ).resolve()
                pair_record = raw["versions"][label][
                    "evaluation_artifacts"
                ]["transcripts_tsv"]
                evaluation_artifacts[label] = {
                    "transcripts_tsv": {
                        **pair_record,
                        "mtime_ns": transcript.stat().st_mtime_ns,
                    },
                }
            cell = {
                "id": f"{panel}-{benchmark}-{repetition}",
                "kind": "paired_release",
                "mode": f"paired_{panel}",
                "panel": panel,
                "benchmark": benchmark,
                "repetition": repetition,
                "expected_order": raw["order"],
                "threads": raw["threads"],
                "full_job": panel in {"full", "e2e"},
                "exclusive": True,
                "command": ["python", "-m", "release_evaluation", benchmark],
                "environment": {},
                "paired": json.loads(json.dumps(paired)),
                "cell_dir": str(cell_dir),
                "input_fingerprints": raw["inputs"],
                "evaluation_input_fingerprints": raw["evaluation_inputs"],
                "artifacts": {
                    name: value
                    for name, value in artifacts.items()
                    if not name.endswith("_gff_validation")
                },
            }
            cells.append(cell)
            frozen_inputs[benchmark] = raw["inputs"]
            frozen_evaluation_inputs[benchmark] = raw["evaluation_inputs"]
            pending_success.append(
                (cell, artifacts, observed, evaluation_artifacts)
            )
    tooling_paths = {
        "tooling_build_controller": Path(build_controller.__file__),
        "tooling_release_evaluation": (
            Path(build_controller.REPO_ROOT)
            / "benchmarks" / "compare" / "release_evaluation.py"
        ),
        "tooling_release_report": Path(release_report.__file__),
        "tooling_run_benchmarks": (
            Path(build_controller.REPO_ROOT) / "benchmarks" / "run_benchmarks.py"
        ),
        "tooling_gff3_validator": (
            Path(build_controller.REPO_ROOT) / "lifton" / "gff3_validator.py"
        ),
    }
    provenance = {
        "git": dict(GIT_STATE),
        "files": {
            "benchmark_registry": _source_file_record(benchmark_registry),
            "dataset_registry": _source_file_record(dataset_registry),
            "baseline": _source_file_record(baseline),
            **{
                name: _source_file_record(path)
                for name, path in tooling_paths.items()
            },
        },
        "tools": {
            "python": {
                "path": "/env/bin/python",
                "version": compatibility_tag,
            },
        },
        "runtime": {
            "python": {
                "implementation": "cpython",
                "version": "3.11.15",
            },
            "distributions": {"numpy": {"version": compatibility_tag}},
        },
        "paired": {
            "configuration": json.loads(json.dumps(paired)),
            "sources": json.loads(json.dumps(source_records)),
            "inputs": json.loads(json.dumps(frozen_inputs)),
            "evaluation_inputs": json.loads(json.dumps(
                frozen_evaluation_inputs
            )),
        },
    }
    provenance["fingerprint"] = build_controller.canonical_hash(provenance)
    for cell, artifacts, observed, evaluation_artifacts in pending_success:
        fingerprint = build_controller._cell_fingerprint(
            cell, provenance["fingerprint"],
        )
        cell["fingerprint"] = fingerprint
        cell_dir = Path(cell["cell_dir"])
        recorded = {
            name: release_report._live_artifact_record(Path(path))
            for name, path in artifacts.items()
        }
        (cell_dir / "status.json").write_text(json.dumps({
            "schema_version": build_controller.CONTROLLER_SCHEMA_VERSION,
            "state": "success",
            "fingerprint": fingerprint,
        }))
        (cell_dir / ".success").write_text(json.dumps({
            "schema_version": build_controller.CONTROLLER_SCHEMA_VERSION,
            "fingerprint": fingerprint,
            "validation": {
                "artifacts": recorded,
                "gff_fingerprints": observed,
                "evaluation_artifacts": evaluation_artifacts,
            },
        }))
    plan = {
        "schema_version": build_controller.CONTROLLER_SCHEMA_VERSION,
        "run_id": f"test-{panel}",
        "created_at": "2026-07-17T00:00:00+00:00",
        "repo_root": str(Path(build_controller.REPO_ROOT).resolve()),
        "run_dir": str(root),
        "stage": stage or f"paired-{panel}",
        "ids": benchmarks,
        "policy": {
            name: getattr(build_controller.Policy(), name)
            for name in build_controller.Policy.__dataclass_fields__
        },
        "inputs": {
            "registry": str(benchmark_registry),
            "dataset_registry": str(dataset_registry),
            "baseline": str(baseline),
        },
        "paired": paired,
        "provenance": provenance,
        "cells": cells,
    }
    plan["fingerprint"] = build_controller._plan_fingerprint(plan)
    build_controller.validate_plan_integrity(plan)
    root.mkdir(parents=True, exist_ok=True)
    (root / "plan.json").write_text(json.dumps(plan))
    return root


def _write_release_roots(
    tmp_path: Path,
    *,
    compatibility_tags: dict[str, str] | None = None,
    e2e_modes: tuple[str, str] = ("fast", "fast"),
) -> list[Path]:
    compatibility_tags = compatibility_tags or {}
    return [
        _write_controller_run(
            tmp_path / panel,
            panel,
            compatibility_tag=compatibility_tags.get(
                panel, "shared-toolchain"
            ),
            candidate_mode=e2e_modes[0] if panel == "e2e" else None,
            reference_mode=e2e_modes[1] if panel == "e2e" else None,
        )
        for panel in ("subset", "full", "e2e")
    ]


def test_release_report_aggregates_pairs_and_writes_publication_set(tmp_path):
    runs = tmp_path / "runs"
    _write_pair(runs, repetition=1, candidate_wall=8.0, reference_wall=10.0)
    _write_pair(runs, repetition=2, candidate_wall=9.0, reference_wall=10.0)
    _write_pair(runs, repetition=3, candidate_wall=8.5, reference_wall=10.0)
    output = tmp_path / "publication"

    result = release_report.write_report(
        [runs],
        output,
        candidate_sha="c" * 40,
        reference_sha=REFERENCE_SHA,
        seed=17,
        replicates=200,
        diagnostic=True,
    )

    assert set(path.name for path in output.iterdir()) == {
        "REPORT.md",
        "manifest.json",
        "metrics.json",
    }
    panel = result["metrics"]["panels"]["subset"]
    assert panel["n_cells"] == 1
    assert panel["wall_ratio"]["estimate"] == 0.85
    assert panel["rss_ratio"]["estimate"] == 0.5
    assert panel["candidate_wall_seconds_total_of_cell_medians"] == 8.5
    assert panel["reference_wall_seconds_total_of_cell_medians"] == 10.0
    assert panel["wall_total_ratio"] == 0.85
    assert panel["candidate_concurrent_memory_envelope_gib"] == pytest.approx(
        100.0 / 1024.0
    )
    cell = result["metrics"]["cells"][0]
    assert cell["candidate_wall_seconds_median"] == 8.5
    assert cell["reference_wall_seconds_median"] == 10.0
    assert result["metrics"]["cells"][0]["candidate_byte_deterministic"] is True
    assert result["metrics"]["verdict"]["passed"] is False
    assert result["metrics"]["verdict"]["diagnostic_passed"] is True
    transcript_metrics = result["manifest"]["pair_results"][0][
        "transcript_metrics"
    ]
    assert set(transcript_metrics) == {"candidate", "reference"}
    assert all(
        len(metadata["sha256"]) == 64
        for metadata in transcript_metrics.values()
    )
    report = (output / "REPORT.md").read_text()
    assert "**Verdict:** DIAGNOSTIC ONLY (CHECKS PASS)" in report
    assert "Total wall" in report
    assert "Concurrent candidate RSS" in report
    assert "`" + "c" * 40 + "`" in report


def test_complete_controller_campaign_can_publish_release_pass(tmp_path):
    roots = _write_release_roots(tmp_path / "runs")
    output = tmp_path / "publication"

    result = release_report.write_report(
        roots,
        output,
        candidate_sha=CANDIDATE_SHA,
        reference_sha=REFERENCE_SHA,
        expected_campaign=_campaign_spec(),
        seed=17,
        replicates=100,
    )

    assert result["metrics"]["verdict"]["passed"] is True
    assert result["metrics"]["campaign"]["matrix_complete"] is True
    assert result["metrics"]["publication"]["release_evidence_valid"] is True
    assert result["manifest"]["publication_mode"] == "release"
    report = (output / "REPORT.md").read_text()
    assert "**Verdict:** PASS" in report
    assert "| e2e | fast | fast |" in report


def test_release_publication_rejects_canary_or_incomplete_matrix(tmp_path):
    roots = _write_release_roots(tmp_path / "runs")
    subset_plan = roots[0] / "plan.json"
    plan = json.loads(subset_plan.read_text())
    plan["stage"] = "paired-subset-canary"
    subset_plan.write_text(json.dumps(plan))

    with pytest.raises(ValueError, match="canary stage"):
        release_report.write_report(
            roots,
            tmp_path / "canary-publication",
            candidate_sha=CANDIDATE_SHA,
            reference_sha=REFERENCE_SHA,
            expected_campaign=_campaign_spec(),
            replicates=10,
        )

    roots = _write_release_roots(tmp_path / "other-runs")
    incomplete = _campaign_spec()
    incomplete["panels"]["full"]["ids"] = incomplete["panels"]["full"]["ids"][:-1]
    with pytest.raises(ValueError, match="must enumerate the canonical 17 IDs"):
        release_report.write_report(
            roots,
            tmp_path / "incomplete-publication",
            candidate_sha=CANDIDATE_SHA,
            reference_sha=REFERENCE_SHA,
            expected_campaign=incomplete,
            replicates=10,
        )


def test_release_campaign_requires_even_balanced_repetitions():
    campaign = _campaign_spec()
    for panel in campaign["panels"].values():
        panel["repetitions"] = 3

    with pytest.raises(ValueError, match="positive even integer"):
        release_report._normalize_campaign_spec(campaign)


def test_release_campaign_rejects_explicit_but_partial_id_matrix():
    campaign = _campaign_spec()
    assert {
        panel: len(row["ids"])
        for panel, row in campaign["panels"].items()
    } == {"subset": 34, "full": 17, "e2e": 5}
    campaign["panels"]["subset"]["ids"] = campaign["panels"]["subset"]["ids"][:2]

    with pytest.raises(ValueError, match="canonical 34 IDs"):
        release_report._normalize_campaign_spec(campaign)


def test_verdict_rejects_invalid_or_nondeterministic_candidate():
    metrics = {
        "cells": [{
            "panel": "subset",
            "benchmark": "demo",
            "candidate_valid": False,
            "candidate_byte_deterministic": False,
            "quality_deltas": {},
            "common_pi": None,
        }],
        "panels": {},
    }

    verdict = release_report.evaluate_verdict(metrics)

    assert verdict["passed"] is False
    failed = {
        check["name"] for check in verdict["checks"]
        if not check["passed"]
    }
    assert {
        "complete_immutable_release_campaign",
        "all_candidate_outputs_valid",
        "candidate_outputs_byte_deterministic",
        "subset.demo.candidate.positive_coding_recovery",
        "subset.demo.candidate.valid_pi",
        "subset.demo.common_pi_present",
    } <= failed


def test_wall_total_ratio_uses_sums_of_per_cell_arm_medians(tmp_path):
    runs = tmp_path / "runs"
    for repetition, candidate, reference in (
        (1, 8.0, 10.0),
        (2, 12.0, 10.0),
        (3, 10.0, 10.0),
    ):
        _write_pair(
            runs, panel="subset", benchmark="short",
            repetition=repetition, candidate_wall=candidate,
            reference_wall=reference,
        )
    for repetition, candidate, reference in (
        (1, 9.0, 90.0),
        (2, 11.0, 110.0),
        (3, 10.0, 100.0),
    ):
        _write_pair(
            runs, panel="subset", benchmark="long",
            repetition=repetition, candidate_wall=candidate,
            reference_wall=reference,
        )

    metrics = release_report.aggregate_pairs(
        release_report.load_pairs([runs]), seed=17, replicates=100,
    )

    cells = {
        cell["benchmark"]: cell
        for cell in metrics["cells"]
    }
    assert cells["short"]["candidate_wall_seconds_median"] == 10.0
    assert cells["short"]["reference_wall_seconds_median"] == 10.0
    assert cells["long"]["candidate_wall_seconds_median"] == 10.0
    assert cells["long"]["reference_wall_seconds_median"] == 100.0
    panel = metrics["panels"]["subset"]
    assert panel["candidate_wall_seconds_total_of_cell_medians"] == 20.0
    assert panel["reference_wall_seconds_total_of_cell_medians"] == 110.0
    assert panel["wall_total_ratio"] == pytest.approx(20.0 / 110.0)
    # This must not regress to either an equal-cell mean (0.55) or a
    # reference-annotation-size-weighted mean of per-cell ratios.
    assert panel["wall_total_ratio"] != pytest.approx(0.55)


def test_concurrent_memory_envelope_uses_panel_scheduler_caps(tmp_path):
    runs = tmp_path / "runs"
    for index, gib in enumerate((10, 20, 30, 40, 50)):
        _write_pair(
            runs, panel="subset", benchmark=f"s{index}",
            repetition=1, candidate_wall=10.0, reference_wall=10.0,
            candidate_rss_mb=gib * 1024.0,
            reference_rss_mb=gib * 1024.0,
        )
    for index, gib in enumerate((100, 100, 1)):
        _write_pair(
            runs, panel="full", benchmark=f"f{index}",
            repetition=1, candidate_wall=10.0, reference_wall=10.0,
            candidate_rss_mb=gib * 1024.0,
            reference_rss_mb=gib * 1024.0,
        )

    metrics = release_report.aggregate_pairs(
        release_report.load_pairs([runs]), seed=17, replicates=100,
    )

    subset = metrics["panels"]["subset"]
    assert subset["candidate_concurrency_limit"] == 4
    assert subset["candidate_concurrent_memory_envelope_gib"] == 140
    assert [
        item["benchmark"]
        for item in subset["candidate_concurrent_memory_contributors"]
    ] == ["s4", "s3", "s2", "s1"]
    full = metrics["panels"]["full"]
    assert full["candidate_concurrency_limit"] == 2
    assert full["candidate_concurrent_memory_envelope_gib"] == 200
    memory_checks = {
        check["name"]: check
        for check in metrics["verdict"]["checks"]
        if "concurrent_memory_envelope" in check["name"]
    }
    assert memory_checks[
        "subset.candidate_concurrent_memory_envelope_gib"
    ]["passed"] is True
    assert memory_checks[
        "full.candidate_concurrent_memory_envelope_gib"
    ]["passed"] is False
    assert metrics["verdict"]["passed"] is False


def test_e2e_biology_is_validated_for_every_repetition_and_reported(tmp_path):
    runs = tmp_path / "runs"
    _write_pair(
        runs, panel="e2e", benchmark="bee", repetition=1,
        candidate_wall=10.0, reference_wall=10.0,
        candidate_rss_mb=1024.0, reference_rss_mb=1024.0,
    )
    malformed_reference = {
        "schema_version": 1,
        "reference_features": 10,
        "mapped_features": 9,
        "lost_features": 0,
        "extra_copy_features": 0,
        "feature_completeness": "0.9",
        "emitted_transcript_records": 12,
        "mapped_transcripts_reported": 11,
        "evaluated_transcript_records": 12,
        "evaluated_coding_records": 13,
        "mean_protein_identity": float("nan"),
    }
    _write_pair(
        runs, panel="e2e", benchmark="bee", repetition=2,
        candidate_wall=10.0, reference_wall=10.0,
        candidate_rss_mb=1024.0, reference_rss_mb=1024.0,
        candidate_biology=None,
        reference_biology=malformed_reference,
    )

    metrics = release_report.aggregate_pairs(
        release_report.load_pairs([runs]), seed=17, replicates=100,
    )

    biology = metrics["cells"][0]["e2e_biology"]
    assert len(biology["candidate"]["repetition_summaries"]) == 2
    assert len(biology["reference"]["repetition_summaries"]) == 2
    assert biology["candidate"]["valid"] is False
    assert biology["reference"]["valid"] is False
    assert any(
        "repetition 2: end-to-end result has no biological_summary"
        in error for error in biology["candidate"]["errors"]
    )
    assert any(
        "does not equal reference_features" in error
        for error in biology["reference"]["errors"]
    )
    assert any(
        "feature_completeness is not finite" in error
        for error in biology["reference"]["errors"]
    )
    assert any(
        "mean_protein_identity is not finite" in error
        for error in biology["reference"]["errors"]
    )
    checks = {
        check["name"]: check
        for check in metrics["verdict"]["checks"]
    }
    assert checks["e2e.bee.candidate.biological_summary"]["passed"] is False
    assert checks["e2e.bee.reference.biological_summary"]["passed"] is False
    assert metrics["verdict"]["passed"] is False
    markdown = release_report.render_markdown(
        metrics,
        {"candidate_sha": "c" * 40, "reference_sha": "r" * 40},
    )
    assert "## End-to-end biological validation" in markdown
    assert "| bee | candidate | FAIL |" in markdown


def test_valid_e2e_biology_and_two_cell_memory_envelope_pass(tmp_path):
    runs = tmp_path / "runs"
    for benchmark, gib in (("bee", 30), ("rice", 20), ("mouse", 10)):
        _write_pair(
            runs, panel="e2e", benchmark=benchmark, repetition=1,
            candidate_wall=10.0, reference_wall=10.0,
            candidate_rss_mb=gib * 1024.0,
            reference_rss_mb=gib * 1024.0,
        )

    metrics = release_report.aggregate_pairs(
        release_report.load_pairs([runs]), seed=17, replicates=100,
    )

    assert all(
        version["valid"]
        for cell in metrics["cells"]
        for version in cell["e2e_biology"].values()
    )
    panel = metrics["panels"]["e2e"]
    assert panel["candidate_concurrency_limit"] == 2
    assert panel["candidate_concurrent_memory_envelope_gib"] == 50
    assert metrics["verdict"]["passed"] is False
    assert metrics["verdict"]["diagnostic_passed"] is True


def test_report_rejects_mismatched_cli_provenance_before_publication(tmp_path):
    runs = tmp_path / "runs"
    _write_pair(
        runs, repetition=1, candidate_wall=8.0, reference_wall=10.0,
        candidate_sha="a" * 40,
    )
    output = tmp_path / "publication"

    with pytest.raises(ValueError, match="candidate provenance SHA"):
        release_report.write_report(
            [runs],
            output,
            candidate_sha=CANDIDATE_SHA,
            reference_sha=REFERENCE_SHA,
            seed=17,
            replicates=100,
            diagnostic=True,
        )

    assert not output.exists()


def test_campaign_rejects_duplicate_panel_benchmark_repetition(tmp_path):
    first = tmp_path / "runs-a"
    second = tmp_path / "runs-b"
    _write_pair(
        first, panel="subset", benchmark="demo", repetition=1,
        candidate_wall=8.0, reference_wall=10.0,
    )
    _write_pair(
        second, panel="subset", benchmark="demo", repetition=1,
        candidate_wall=8.0, reference_wall=10.0,
    )

    with pytest.raises(ValueError, match="duplicate paired repetition"):
        release_report.aggregate_pairs(
            release_report.load_pairs([first, second]),
            seed=17,
            replicates=100,
        )


def test_controller_discovery_uses_only_current_successful_cell_result(tmp_path):
    run = tmp_path / "controller-run"
    cell_dir = run / "cells" / "paired-cell"
    current = _write_pair(
        cell_dir / "pair",
        panel="subset",
        benchmark="demo",
        repetition=1,
        candidate_wall=8.0,
        reference_wall=10.0,
    )
    _write_pair(
        cell_dir / "attempt-archive",
        panel="subset",
        benchmark="demo",
        repetition=1,
        candidate_wall=99.0,
        reference_wall=10.0,
    )
    fingerprint = "f" * 64
    cell = {
        "id": "paired-cell",
        "kind": "paired_release",
        "cell_dir": str(cell_dir),
        "fingerprint": fingerprint,
        "artifacts": {"result_json": str(current)},
    }
    run.mkdir(parents=True, exist_ok=True)
    (run / "plan.json").write_text(json.dumps({
        "schema_version": build_controller.CONTROLLER_SCHEMA_VERSION,
        "cells": [cell],
    }))
    (cell_dir / "status.json").write_text(json.dumps({
        "schema_version": build_controller.CONTROLLER_SCHEMA_VERSION,
        "state": "success",
        "fingerprint": fingerprint,
    }))
    (cell_dir / ".success").write_text(json.dumps({
        "schema_version": build_controller.CONTROLLER_SCHEMA_VERSION,
        "fingerprint": fingerprint,
        "validation": {
            "artifacts": {
                "result_json": release_report._live_artifact_record(current),
            }
        },
    }))

    pairs = release_report.load_pairs([run])

    assert [path for path, _pair in pairs] == [current]
    current.write_text(current.read_text() + "\n")
    with pytest.raises(ValueError, match="changed after success"):
        release_report.load_pairs([run])


def test_controller_discovery_rejects_incomplete_cells(tmp_path):
    run = tmp_path / "controller-run"
    cell_dir = run / "cells" / "paired-cell"
    current = _write_pair(
        cell_dir / "pair",
        repetition=1,
        candidate_wall=8.0,
        reference_wall=10.0,
    )
    fingerprint = "f" * 64
    run.mkdir(parents=True, exist_ok=True)
    (run / "plan.json").write_text(json.dumps({
        "schema_version": build_controller.CONTROLLER_SCHEMA_VERSION,
        "cells": [{
            "id": "paired-cell",
            "kind": "paired_release",
            "cell_dir": str(cell_dir),
            "fingerprint": fingerprint,
            "artifacts": {"result_json": str(current)},
        }],
    }))
    (cell_dir / "status.json").write_text(json.dumps({
        "schema_version": build_controller.CONTROLLER_SCHEMA_VERSION,
        "state": "failed",
        "fingerprint": fingerprint,
    }))
    (cell_dir / ".success").write_text(json.dumps({
        "schema_version": build_controller.CONTROLLER_SCHEMA_VERSION,
        "fingerprint": fingerprint,
        "validation": {
            "artifacts": {
                "result_json": release_report._live_artifact_record(current),
            }
        },
    }))

    with pytest.raises(ValueError, match="state is 'failed'"):
        release_report.load_pairs([run])


@pytest.mark.parametrize(
    "changed",
    (
        {"input_tag": "different-input"},
        {"threads": 4},
        {"candidate_mode": "safe"},
    ),
    ids=("inputs", "threads", "modes"),
)
def test_campaign_rejects_inconsistent_repetition_configuration(
        tmp_path, changed):
    runs = tmp_path / "runs"
    _write_pair(
        runs, panel="subset", benchmark="demo", repetition=1,
        candidate_wall=8.0, reference_wall=10.0,
    )
    _write_pair(
        runs, panel="subset", benchmark="demo", repetition=2,
        candidate_wall=8.0, reference_wall=10.0,
        **changed,
    )

    with pytest.raises(ValueError, match="changed across repetitions"):
        release_report.aggregate_pairs(
            release_report.load_pairs([runs]), seed=17, replicates=100,
        )


def test_campaign_rejects_malformed_order_and_nonfinite_raw_profile(tmp_path):
    malformed_order = tmp_path / "bad-order"
    _write_pair(
        malformed_order, repetition=1,
        candidate_wall=8.0, reference_wall=10.0,
        order=["candidate", "candidate"],
    )
    with pytest.raises(ValueError, match="order must contain"):
        release_report.aggregate_pairs(
            release_report.load_pairs([malformed_order]),
            seed=17,
            replicates=100,
        )

    nonfinite = tmp_path / "nonfinite"
    _write_pair(
        nonfinite, repetition=1,
        candidate_wall=float("nan"), reference_wall=10.0,
    )
    with pytest.raises(ValueError, match="must be positive and finite"):
        release_report.aggregate_pairs(
            release_report.load_pairs([nonfinite]),
            seed=17,
            replicates=100,
        )


def test_campaign_enforces_alternation_and_contiguous_repetitions(tmp_path):
    wrong_order = tmp_path / "wrong-order"
    _write_pair(
        wrong_order, repetition=1,
        candidate_wall=8.0, reference_wall=10.0,
        order=["candidate", "reference"],
    )
    with pytest.raises(ValueError, match="violates alternating AB/BA"):
        release_report.aggregate_pairs(
            release_report.load_pairs([wrong_order]),
            seed=17,
            replicates=100,
        )

    missing_first = tmp_path / "missing-first"
    _write_pair(
        missing_first, repetition=2,
        candidate_wall=8.0, reference_wall=10.0,
    )
    with pytest.raises(ValueError, match="contiguous from 1"):
        release_report.aggregate_pairs(
            release_report.load_pairs([missing_first]),
            seed=17,
            replicates=100,
        )


def test_campaign_rejects_mixed_protocols_within_a_panel(tmp_path):
    runs = tmp_path / "runs"
    _write_pair(
        runs, panel="e2e", benchmark="bee", repetition=1,
        candidate_wall=8.0, reference_wall=10.0,
        candidate_mode="fast", reference_mode="safe",
    )
    _write_pair(
        runs, panel="e2e", benchmark="rice", repetition=1,
        candidate_wall=8.0, reference_wall=10.0,
        candidate_mode="stream", reference_mode="safe",
    )

    with pytest.raises(ValueError, match="modes differ across"):
        release_report.aggregate_pairs(
            release_report.load_pairs([runs]),
            seed=17,
            replicates=100,
        )


def test_reference_validity_and_byte_determinism_are_control_gates(tmp_path):
    runs = tmp_path / "runs"
    _write_pair(
        runs, repetition=1, candidate_wall=8.0, reference_wall=10.0,
    )
    _write_pair(
        runs, repetition=2, candidate_wall=8.0, reference_wall=10.0,
        reference_valid=False,
        reference_byte_sha="changed-reference-bytes",
    )

    metrics = release_report.aggregate_pairs(
        release_report.load_pairs([runs]), seed=17, replicates=100,
    )

    cell = metrics["cells"][0]
    assert cell["candidate_byte_deterministic"] is True
    assert cell["reference_byte_deterministic"] is False
    assert cell["reference_valid"] is False
    checks = {
        check["name"]: check
        for check in metrics["verdict"]["checks"]
    }
    assert checks["all_reference_outputs_valid"]["passed"] is False
    assert checks["reference_outputs_byte_deterministic"]["passed"] is False
    assert metrics["verdict"]["passed"] is False


def test_transcript_metrics_count_recovered_rows_with_blank_pi(tmp_path):
    path = tmp_path / "transcripts.tsv"
    path.write_text(
        TSV_HEADER
        + _tsv_row("tx1", protein_identity="", status="map_failed")
        + _tsv_row("tx2", protein_identity="0.5")
    )

    metrics = release_report._load_transcript_metrics(path, 2)

    assert metrics["n_recovered_coding"] == 2
    assert metrics["n_recovered_coding_with_pi"] == 1
    assert metrics["sum_protein_identity"] == 0.5
    assert metrics["covpi"] == 0.25
    assert metrics["recall_at_0.5"] == 0.5


@pytest.mark.parametrize(
    ("rows", "message"),
    (
        (
            _tsv_row("tx1") + _tsv_row("tx1"),
            "duplicate coding reference row",
        ),
        (
            _tsv_row("tx1"),
            "found 1 unique coding reference rows, expected 2",
        ),
        (
            _tsv_row("tx1", protein_identity="nan") + _tsv_row("tx2"),
            "protein_identity must be finite and within",
        ),
        (
            _tsv_row("tx1", dna_identity="1.1") + _tsv_row("tx2"),
            "dna_identity must be finite and within",
        ),
    ),
    ids=("duplicate", "missing", "nonfinite-pi", "out-of-range-dna"),
)
def test_transcript_metrics_reject_malformed_tsv(tmp_path, rows, message):
    path = tmp_path / "transcripts.tsv"
    path.write_text(TSV_HEADER + rows)

    with pytest.raises(ValueError, match=message):
        release_report._load_transcript_metrics(path, 2)


def test_campaign_validates_transcript_tsv_for_every_repetition(tmp_path):
    runs = tmp_path / "runs"
    _write_pair(
        runs, repetition=1, candidate_wall=8.0, reference_wall=10.0,
    )
    second = _write_pair(
        runs, repetition=2, candidate_wall=8.0, reference_wall=10.0,
    )
    (second.parent / "evaluation" / "reference.transcripts.tsv").write_text(
        TSV_HEADER + _tsv_row("tx1") + _tsv_row("tx1")
    )

    with pytest.raises(ValueError, match="duplicate coding reference row"):
        release_report.aggregate_pairs(
            release_report.load_pairs([runs]), seed=17, replicates=100,
        )


def test_release_rehashes_gff_and_rejects_changed_evaluator_tsv(tmp_path):
    roots = _write_release_roots(tmp_path / "gff-runs")
    subset_plan = json.loads((roots[0] / "plan.json").read_text())
    gff = Path(subset_plan["cells"][0]["artifacts"]["candidate_gff"])
    stat = gff.stat()
    gff.write_text(gff.read_text().replace("\tgene\t", "\txene\t"))
    os.utime(gff, ns=(stat.st_atime_ns, stat.st_mtime_ns))

    with pytest.raises(ValueError, match="changed after success"):
        release_report.write_report(
            roots,
            tmp_path / "gff-publication",
            candidate_sha=CANDIDATE_SHA,
            reference_sha=REFERENCE_SHA,
            expected_campaign=_campaign_spec(),
            replicates=10,
        )

    roots = _write_release_roots(tmp_path / "tsv-runs")
    plan = json.loads((roots[0] / "plan.json").read_text())
    result_path = Path(plan["cells"][0]["artifacts"]["result_json"])
    transcript = (
        result_path.parent / "evaluation" / "candidate.transcripts.tsv"
    )
    transcript_stat = transcript.stat()
    transcript.write_text(
        transcript.read_text().replace("tx1\ttool-tx1\t1", "tx1\ttool-tx1\t0", 1)
    )
    os.utime(
        transcript,
        ns=(transcript_stat.st_atime_ns, transcript_stat.st_mtime_ns),
    )

    with pytest.raises(ValueError, match="cryptographic controller evidence"):
        release_report.write_report(
            roots,
            tmp_path / "tsv-publication",
            candidate_sha=CANDIDATE_SHA,
            reference_sha=REFERENCE_SHA,
            expected_campaign=_campaign_spec(),
            replicates=10,
        )


@pytest.mark.parametrize("missing_from", ("pair", "controller"))
def test_release_requires_both_evaluator_tsv_evidence_layers(
        tmp_path, missing_from):
    roots = _write_release_roots(tmp_path / "runs")
    plan = json.loads((roots[0] / "plan.json").read_text())
    cell = plan["cells"][0]
    result_path = Path(cell["artifacts"]["result_json"])
    if missing_from == "pair":
        raw = json.loads(result_path.read_text())
        raw["versions"]["candidate"].pop("evaluation_artifacts")
        result_path.write_text(json.dumps(raw))
        success_path = Path(cell["cell_dir"]) / ".success"
        success = json.loads(success_path.read_text())
        result_record = success["validation"]["artifacts"]["result_json"]
        result_record["size"] = result_path.stat().st_size
        result_record["mtime_ns"] = result_path.stat().st_mtime_ns
        result_record["sha256"] = release_report.sha256_file(result_path)
        success_path.write_text(json.dumps(success))
        message = "pair result lacks evaluator TSV evidence"
    else:
        success_path = Path(cell["cell_dir"]) / ".success"
        success = json.loads(success_path.read_text())
        success["validation"]["evaluation_artifacts"]["candidate"].pop(
            "transcripts_tsv"
        )
        success_path.write_text(json.dumps(success))
        message = "controller success lacks evaluator TSV evidence"

    with pytest.raises(ValueError, match=message):
        release_report.write_report(
            roots,
            tmp_path / "publication",
            candidate_sha=CANDIDATE_SHA,
            reference_sha=REFERENCE_SHA,
            expected_campaign=_campaign_spec(),
            replicates=10,
        )


def test_release_rejects_cross_root_toolchain_drift(tmp_path):
    roots = _write_release_roots(
        tmp_path / "runs",
        compatibility_tags={"full": "different-python"},
    )

    with pytest.raises(ValueError, match="differs across release roots"):
        release_report.write_report(
            roots,
            tmp_path / "publication",
            candidate_sha=CANDIDATE_SHA,
            reference_sha=REFERENCE_SHA,
            expected_campaign=_campaign_spec(),
            replicates=10,
        )


def test_release_discloses_and_gates_mismatched_e2e_modes(tmp_path):
    roots = _write_release_roots(
        tmp_path / "runs",
        e2e_modes=("fast", "safe"),
    )

    result = release_report.write_report(
        roots,
        tmp_path / "publication",
        candidate_sha=CANDIDATE_SHA,
        reference_sha=REFERENCE_SHA,
        expected_campaign=_campaign_spec(),
        replicates=10,
    )

    assert result["metrics"]["verdict"]["passed"] is False
    checks = {
        check["name"]: check
        for check in result["metrics"]["verdict"]["checks"]
    }
    assert checks[
        "e2e.bee.candidate_reference_modes_match"
    ]["passed"] is False
    report = (tmp_path / "publication" / "REPORT.md").read_text()
    assert "| e2e | fast | safe |" in report


def test_required_quality_metrics_fail_closed(tmp_path):
    runs = tmp_path / "runs"
    path = _write_pair(
        runs, repetition=1, candidate_wall=8.0, reference_wall=10.0,
    )
    raw = json.loads(path.read_text())
    raw["versions"]["candidate"]["summary"].pop("completeness_coding")
    path.write_text(json.dumps(raw))

    with pytest.raises(ValueError, match="completeness_coding"):
        release_report.aggregate_pairs(
            release_report.load_pairs([runs]), replicates=10,
        )

    transcript = tmp_path / "missing-structure.tsv"
    transcript.write_text(
        TSV_HEADER.replace("intron_chain_exact\t", "")
        + _tsv_row("tx1")
        + _tsv_row("tx2")
    )
    with pytest.raises(ValueError, match="lacks columns"):
        release_report._load_transcript_metrics(transcript, 2)

    incomplete = tmp_path / "missing-structural-values.tsv"
    incomplete.write_text(
        TSV_HEADER
        + _tsv_row("tx1").replace("\t1\t1.0\t1.0\t1\t1\t1\n", "\t\t\t\t1\t1\t\n")
        + _tsv_row("tx2")
    )
    with pytest.raises(ValueError, match="incomplete structural quality"):
        release_report._load_transcript_metrics(incomplete, 2)


def test_feature_completeness_must_match_required_summary(tmp_path):
    runs = tmp_path / "runs"
    path = _write_pair(
        runs, repetition=1, candidate_wall=8.0, reference_wall=10.0,
    )
    raw = json.loads(path.read_text())
    overall = raw["versions"]["candidate"]["summary"][
        "completeness_by_type"
    ]["_overall_"]
    overall["n_recovered"] = 1
    path.write_text(json.dumps(raw))

    with pytest.raises(ValueError, match="fraction disagrees with its counts"):
        release_report.aggregate_pairs(
            release_report.load_pairs([runs]), replicates=10,
        )


def test_stable_id_preservation_is_reported_and_gates_regression(tmp_path):
    runs = tmp_path / "runs"
    path = _write_pair(
        runs, repetition=1, candidate_wall=8.0, reference_wall=10.0,
    )
    raw = json.loads(path.read_text())
    candidate_cds = raw["versions"]["candidate"]["summary"][
        "stable_id_preservation"
    ]["by_type"]["CDS"]
    candidate_cds["n_preserved_ids"] = 1
    candidate_cds["preservation_rate"] = 0.5
    path.write_text(json.dumps(raw))

    metrics = release_report.aggregate_pairs(
        release_report.load_pairs([runs]), replicates=10,
    )

    cell = metrics["cells"][0]["stable_id_preservation"]["CDS"]
    assert cell["candidate_n_preserved_ids"] == 1
    assert cell["candidate_n_output_records"] == 2
    assert cell["candidate_n_output_records_with_id"] == 2
    assert cell["reference_n_preserved_ids"] == 2
    checks = {
        check["name"]: check
        for check in metrics["verdict"]["checks"]
    }
    check = checks[
        "subset.demo.stable_id_preservation.CDS.no_regression"
    ]
    assert check["passed"] is False
    panel = metrics["panels"]["subset"]["stable_id_preservation"]["CDS"]
    assert panel["candidate_preservation_rate"] == 0.5
    assert panel["reference_preservation_rate"] == 1.0
    markdown = release_report.render_markdown(
        metrics,
        {"candidate_sha": CANDIDATE_SHA, "reference_sha": REFERENCE_SHA},
    )
    assert "## Stable feature-ID preservation" in markdown
    assert "not a measure of biological completeness" in markdown
    assert "| subset | CDS | 1 | 1/2 (0.50000) | 2/2 (1.00000) | -0.50000 |" in markdown


def test_stable_id_preservation_skips_types_without_declared_reference_ids(
        tmp_path):
    runs = tmp_path / "runs"
    path = _write_pair(
        runs, repetition=1, candidate_wall=8.0, reference_wall=10.0,
    )
    raw = json.loads(path.read_text())
    for version in raw["versions"].values():
        for feature_type, row in version["summary"][
            "stable_id_preservation"
        ]["by_type"].items():
            row.update({
                "applicable": False,
                "reason": (
                    "reference_feature_type_absent"
                    if feature_type == "CDS"
                    else "no_declared_reference_ids"
                ),
                "n_reference_records": 0 if feature_type == "CDS" else 2,
                "n_reference_records_with_id": 0,
                "n_reference_ids": 0,
                "n_preserved_ids": 0,
                "n_output_records": 0 if feature_type == "CDS" else 2,
                "n_output_records_with_id": 0,
                "n_output_ids": 0,
                "preservation_rate": None,
            })
    path.write_text(json.dumps(raw))

    metrics = release_report.aggregate_pairs(
        release_report.load_pairs([runs]), replicates=10,
    )

    assert not any(
        "stable_id_preservation" in check["name"]
        for check in metrics["verdict"]["checks"]
    )
    panel = metrics["panels"]["subset"]["stable_id_preservation"]
    assert panel["CDS"]["n_applicable_cells"] == 0
    assert panel["exon"]["candidate_preservation_rate"] is None


def test_stable_id_preservation_rejects_arm_denominator_mismatch(tmp_path):
    runs = tmp_path / "runs"
    path = _write_pair(
        runs, repetition=1, candidate_wall=8.0, reference_wall=10.0,
    )
    raw = json.loads(path.read_text())
    candidate_cds = raw["versions"]["candidate"]["summary"][
        "stable_id_preservation"
    ]["by_type"]["CDS"]
    candidate_cds.update({
        "n_reference_records": 3,
        "n_reference_records_with_id": 3,
        "n_reference_ids": 3,
        "n_preserved_ids": 2,
        "preservation_rate": 2 / 3,
    })
    path.write_text(json.dumps(raw))

    with pytest.raises(ValueError, match="stable ID denominators"):
        release_report.aggregate_pairs(
            release_report.load_pairs([runs]), replicates=10,
        )


def test_stable_id_preservation_rejects_rate_count_disagreement(tmp_path):
    runs = tmp_path / "runs"
    path = _write_pair(
        runs, repetition=1, candidate_wall=8.0, reference_wall=10.0,
    )
    raw = json.loads(path.read_text())
    raw["versions"]["candidate"]["summary"]["stable_id_preservation"][
        "by_type"
    ]["exon"]["preservation_rate"] = 0.5
    path.write_text(json.dumps(raw))

    with pytest.raises(ValueError, match="stable ID rate"):
        release_report.aggregate_pairs(
            release_report.load_pairs([runs]), replicates=10,
        )


@pytest.mark.parametrize(
    "updates",
    (
        {
            "applicable": False,
            "reason": "no_declared_reference_ids",
            "n_reference_ids": 0,
            "n_preserved_ids": 0,
            "preservation_rate": None,
        },
        {"n_output_records": 1},
        {"n_preserved_ids": 0, "n_output_ids": 0},
        {"n_output_ids": 3},
    ),
)
def test_stable_id_preservation_rejects_impossible_counts(
        tmp_path, updates):
    runs = tmp_path / "runs"
    path = _write_pair(
        runs, repetition=1, candidate_wall=8.0, reference_wall=10.0,
    )
    raw = json.loads(path.read_text())
    raw["versions"]["candidate"]["summary"]["stable_id_preservation"][
        "by_type"
    ]["CDS"].update(updates)
    path.write_text(json.dumps(raw))

    with pytest.raises(ValueError, match="stable ID counts"):
        release_report.aggregate_pairs(
            release_report.load_pairs([runs]), replicates=10,
        )


def test_load_pairs_rejects_legacy_release_schema_2(tmp_path):
    runs = tmp_path / "runs"
    path = _write_pair(
        runs, repetition=1, candidate_wall=8.0, reference_wall=10.0,
    )
    raw = json.loads(path.read_text())
    raw["schema_version"] = 2
    path.write_text(json.dumps(raw))

    assert release_report.SCHEMA_VERSION == 4
    with pytest.raises(ValueError, match="unsupported paired schema"):
        release_report.load_pairs([runs])


def test_quality_aggregates_and_gates_drift_across_repetitions(tmp_path):
    runs = tmp_path / "runs"
    _write_pair(
        runs, repetition=1, candidate_wall=8.0, reference_wall=10.0,
    )
    second = _write_pair(
        runs, repetition=2, candidate_wall=8.0, reference_wall=10.0,
    )
    candidate_tsv = (
        second.parent / "evaluation" / "candidate.transcripts.tsv"
    )
    candidate_tsv.write_text(
        TSV_HEADER
        + _tsv_row("tx1", protein_identity="0.0")
        + _tsv_row("tx2", protein_identity="0.0")
    )
    raw = json.loads(second.read_text())
    raw["versions"]["candidate"]["summary"]["protein_identity"] = {
        "n": 2,
        "mean": 0.0,
        "median": 0.0,
    }
    second.write_text(json.dumps(raw))

    metrics = release_report.aggregate_pairs(
        release_report.load_pairs([runs]), replicates=10,
    )

    cell = metrics["cells"][0]
    assert cell["candidate"]["covpi"] == pytest.approx(0.5)
    assert len(cell["candidate"]["repetition_metrics"]) == 2
    assert len(cell["common_pi"]["repetition_metrics"]) == 2
    assert cell["candidate_quality_deterministic"] is False
    assert cell["common_pi_deterministic"] is False
    assert metrics["verdict"]["diagnostic_passed"] is False


def test_target_truth_determinism_ignores_paths_but_retains_mapping_hashes(
        tmp_path):
    runs = tmp_path / "runs"
    mapping_sha = _sha256("shared-ortholog-map")
    paths = [
        _write_pair(
            runs, repetition=repetition,
            candidate_wall=8.0, reference_wall=10.0,
        )
        for repetition in (1, 2)
    ]
    for path in paths:
        _attach_target_truth_evidence(
            path,
            candidate_mapping_sha=mapping_sha,
            reference_mapping_sha=mapping_sha,
        )

    metrics = release_report.aggregate_pairs(
        release_report.load_pairs([runs]), replicates=10,
    )

    truth = metrics["cells"][0]["target_truth"]
    assert truth["candidate"]["deterministic"] is True
    assert truth["reference"]["deterministic"] is True
    assert (
        truth["candidate"]["summary"]["inputs"]["prediction_gff"]
        != json.loads(paths[1].read_text())["versions"]["candidate"][
            "summary"
        ]["target_truth"]["inputs"]["prediction_gff"]
    )

    _attach_target_truth_evidence(
        paths[1],
        candidate_mapping_sha=_sha256("changed-ortholog-map"),
        reference_mapping_sha=mapping_sha,
    )
    changed = release_report.aggregate_pairs(
        release_report.load_pairs([runs]), replicates=10,
    )["cells"][0]["target_truth"]

    assert changed["candidate"]["deterministic"] is False
    assert changed["reference"]["deterministic"] is True


def test_bootstrap_replicates_must_be_positive(tmp_path):
    runs = tmp_path / "runs"
    _write_pair(
        runs, repetition=1, candidate_wall=8.0, reference_wall=10.0,
    )

    with pytest.raises(ValueError, match="replicates must be"):
        release_report.aggregate_pairs(
            release_report.load_pairs([runs]), replicates=0,
        )
    with pytest.raises(SystemExit):
        release_report.build_parser().parse_args([
            "--runs-root", str(runs),
            "--output-dir", str(tmp_path / "out"),
            "--candidate-sha", CANDIDATE_SHA,
            "--reference-sha", REFERENCE_SHA,
            "--diagnostic",
            "--replicates", "0",
        ])


@pytest.mark.parametrize(
    ("diagnostic", "passed", "expected"),
    (
        (False, True, 0),
        (False, False, 1),
        (True, False, 0),
    ),
)
def test_report_cli_exit_status_reflects_release_verdict(
    monkeypatch, tmp_path, diagnostic, passed, expected,
):
    campaign = tmp_path / "campaign.json"
    campaign.write_text("{}")
    monkeypatch.setattr(
        release_report,
        "write_report",
        lambda *args, **kwargs: {
            "metrics": {"verdict": {"passed": passed}},
            "manifest": {},
        },
    )
    argv = [
        "--runs-root", str(tmp_path / "runs"),
        "--output-dir", str(tmp_path / "out"),
        "--candidate-sha", CANDIDATE_SHA,
        "--reference-sha", REFERENCE_SHA,
    ]
    if diagnostic:
        argv.append("--diagnostic")
    else:
        argv.extend(["--campaign-spec", str(campaign)])

    assert release_report.main(argv) == expected


@pytest.mark.parametrize(
    ("old", "new"),
    (
        ('"wall_clock_seconds": 8.0', '"wall_clock_seconds": 9.0'),
        ('"wall": 0.8', '"wall": 0.9'),
    ),
    ids=("profile", "ratio"),
)
def test_release_rejects_same_stat_pair_result_tampering(
        tmp_path, old, new):
    roots = _write_release_roots(tmp_path / "runs")
    plan = json.loads((roots[0] / "plan.json").read_text())
    result_path = Path(plan["cells"][0]["artifacts"]["result_json"])
    stat = result_path.stat()
    original = result_path.read_text()
    changed = original.replace(old, new, 1)
    assert changed != original
    assert len(changed) == len(original)
    result_path.write_text(changed)
    os.utime(result_path, ns=(stat.st_atime_ns, stat.st_mtime_ns))

    with pytest.raises(ValueError, match="changed after success"):
        release_report.write_report(
            roots,
            tmp_path / "publication",
            candidate_sha=CANDIDATE_SHA,
            reference_sha=REFERENCE_SHA,
            expected_campaign=_campaign_spec(),
            replicates=10,
        )


def test_release_rejects_custom_resource_policy_even_when_refingerprinted(
        tmp_path):
    roots = _write_release_roots(tmp_path / "runs")
    plan_path = roots[0] / "plan.json"
    plan = json.loads(plan_path.read_text())
    plan["policy"]["poll_seconds"] = 31.0
    plan["fingerprint"] = build_controller._plan_fingerprint(plan)
    plan_path.write_text(json.dumps(plan))

    with pytest.raises(ValueError, match="resource policy 'poll_seconds'"):
        release_report.write_report(
            roots,
            tmp_path / "publication",
            candidate_sha=CANDIDATE_SHA,
            reference_sha=REFERENCE_SHA,
            expected_campaign=_campaign_spec(),
            replicates=10,
        )


def test_release_rejects_plan_configuration_not_frozen_in_provenance(tmp_path):
    roots = _write_release_roots(tmp_path / "runs")
    plan_path = roots[0] / "plan.json"
    plan = json.loads(plan_path.read_text())
    plan["paired"]["candidate"]["root"] = "/different/candidate"
    for cell in plan["cells"]:
        cell["paired"] = json.loads(json.dumps(plan["paired"]))
        cell["fingerprint"] = build_controller._cell_fingerprint(
            cell, plan["provenance"]["fingerprint"],
        )
    plan["fingerprint"] = build_controller._plan_fingerprint(plan)
    plan_path.write_text(json.dumps(plan))

    with pytest.raises(ValueError, match="frozen provenance"):
        release_report.write_report(
            roots,
            tmp_path / "publication",
            candidate_sha=CANDIDATE_SHA,
            reference_sha=REFERENCE_SHA,
            expected_campaign=_campaign_spec(),
            replicates=10,
        )


def test_release_rejects_pair_registry_not_frozen_in_controller(tmp_path):
    roots = _write_release_roots(tmp_path / "runs")
    plan = json.loads((roots[0] / "plan.json").read_text())
    cell = plan["cells"][0]
    result_path = Path(cell["artifacts"]["result_json"])
    raw = json.loads(result_path.read_text())
    raw["registries"]["dataset"] = "/different/datasets.json"
    result_path.write_text(json.dumps(raw))
    success_path = Path(cell["cell_dir"]) / ".success"
    success = json.loads(success_path.read_text())
    success["validation"]["artifacts"]["result_json"] = (
        release_report._live_artifact_record(result_path)
    )
    success_path.write_text(json.dumps(success))

    with pytest.raises(ValueError, match="pair registries do not match"):
        release_report.write_report(
            roots,
            tmp_path / "publication",
            candidate_sha=CANDIDATE_SHA,
            reference_sha=REFERENCE_SHA,
            expected_campaign=_campaign_spec(),
            replicates=10,
        )


def test_release_requires_current_reporter_git_head(
        tmp_path, monkeypatch):
    roots = _write_release_roots(tmp_path / "runs")
    monkeypatch.setattr(
        release_report,
        "_current_reporter_git_state",
        lambda: {**GIT_STATE, "head": "e" * 40},
    )

    with pytest.raises(ValueError, match="reporter Git HEAD/state"):
        release_report.write_report(
            roots,
            tmp_path / "publication",
            candidate_sha=CANDIDATE_SHA,
            reference_sha=REFERENCE_SHA,
            expected_campaign=_campaign_spec(),
            replicates=10,
        )


def test_absolute_quality_gates_recovery_baseline_and_e2e_floors(tmp_path):
    subset = tmp_path / "subset"
    _write_pair(
        subset,
        panel="subset",
        benchmark="human_mane",
        repetition=1,
        candidate_wall=8.0,
        reference_wall=10.0,
    )
    metrics = release_report.aggregate_pairs(
        release_report.load_pairs([subset]), replicates=10,
    )
    cell = metrics["cells"][0]
    baseline = cell["absolute_quality"]["baseline"]["mean_pi"]
    cell["absolute_quality"]["candidate"]["n_recovered_coding"] = 0
    cell["absolute_quality"]["candidate"]["mean_pi"] = baseline["floor"] - 0.001
    checks = {
        check["name"]: check
        for check in release_report.evaluate_verdict(metrics)["checks"]
    }
    assert checks[
        "subset.human_mane.candidate.positive_coding_recovery"
    ]["passed"] is False
    assert checks[
        "subset.human_mane.candidate.absolute_mean_pi"
    ]["passed"] is False

    e2e = tmp_path / "e2e"
    _write_pair(
        e2e,
        panel="e2e",
        benchmark="bee",
        repetition=1,
        candidate_wall=8.0,
        reference_wall=10.0,
    )
    metrics = release_report.aggregate_pairs(
        release_report.load_pairs([e2e]), replicates=10,
    )
    summary = metrics["cells"][0]["e2e_biology"]["candidate"]["summary"]
    summary["feature_completeness"] = 0.49
    summary["mean_protein_identity"] = 0.49
    checks = {
        check["name"]: check
        for check in release_report.evaluate_verdict(metrics)["checks"]
    }
    assert checks[
        "e2e.bee.candidate.absolute_feature_completeness"
    ]["passed"] is False
    assert checks[
        "e2e.bee.candidate.absolute_mean_protein_identity"
    ]["passed"] is False


def test_release_manifest_hashes_complete_controller_evidence(tmp_path):
    roots = _write_release_roots(tmp_path / "runs")
    result = release_report.write_report(
        roots,
        tmp_path / "publication",
        candidate_sha=CANDIDATE_SHA,
        reference_sha=REFERENCE_SHA,
        expected_campaign=_campaign_spec(),
        replicates=10,
    )

    evidence = result["manifest"]["publication_evidence"]
    assert set(evidence["tooling"]) == set(release_report.TOOLING_FILE_KEYS)
    assert len(evidence["cells"]) == (34 + 17 + 5) * 2
    first = evidence["cells"][0]
    assert len(first["success"]["sha256"]) == 64
    assert {
        "result_json",
        "candidate_gff",
        "reference_gff",
        "candidate_manifest",
        "reference_manifest",
        "candidate_gff_validation",
        "reference_gff_validation",
    } <= set(first["validated_artifacts"])
    assert set(first["evaluator_artifacts"]) == {"candidate", "reference"}
    assert len(
        result["manifest"]["quality_baseline_artifact"]["sha256"]
    ) == 64


def test_failed_republication_quarantines_prior_pass_and_cleans_staging(
        tmp_path):
    roots = _write_release_roots(tmp_path / "runs")
    output = tmp_path / "publication"
    first = release_report.write_report(
        roots,
        output,
        candidate_sha=CANDIDATE_SHA,
        reference_sha=REFERENCE_SHA,
        expected_campaign=_campaign_spec(),
        replicates=10,
    )
    assert first["metrics"]["verdict"]["passed"] is True
    plan = json.loads((roots[0] / "plan.json").read_text())
    result_path = Path(plan["cells"][0]["artifacts"]["result_json"])
    stat = result_path.stat()
    result_path.write_text(
        result_path.read_text().replace('"wall": 0.8', '"wall": 0.9', 1)
    )
    os.utime(result_path, ns=(stat.st_atime_ns, stat.st_mtime_ns))

    with pytest.raises(ValueError, match="changed after success"):
        release_report.write_report(
            roots,
            output,
            candidate_sha=CANDIDATE_SHA,
            reference_sha=REFERENCE_SHA,
            expected_campaign=_campaign_spec(),
            replicates=10,
        )

    assert not output.exists()
    stale = list(tmp_path.glob(".publication.stale-*"))
    assert len(stale) == 1
    assert "**Verdict:** PASS" in (stale[0] / "REPORT.md").read_text()
    assert not list(tmp_path.glob(".publication.tmp-*"))


def test_successful_republication_is_atomic_and_removes_backups(tmp_path):
    runs = tmp_path / "runs"
    _write_pair(
        runs,
        repetition=1,
        candidate_wall=8.0,
        reference_wall=10.0,
    )
    output = tmp_path / "publication"
    arguments = {
        "candidate_sha": CANDIDATE_SHA,
        "reference_sha": REFERENCE_SHA,
        "replicates": 10,
        "diagnostic": True,
    }
    release_report.write_report([runs], output, **arguments)
    release_report.write_report([runs], output, **arguments)

    assert set(path.name for path in output.iterdir()) == {
        "REPORT.md", "manifest.json", "metrics.json",
    }
    assert not list(tmp_path.glob(".publication.stale-*"))
    assert not list(tmp_path.glob(".publication.tmp-*"))


def test_unrecognized_output_directory_is_never_modified(tmp_path):
    output = tmp_path / "unrelated"
    output.mkdir()
    marker = output / "keep-me.txt"
    marker.write_text("important user data\n")

    with pytest.raises(ValueError, match="unrecognized publication"):
        release_report.write_report(
            [tmp_path / "runs"],
            output,
            candidate_sha=CANDIDATE_SHA,
            reference_sha=REFERENCE_SHA,
            diagnostic=True,
            replicates=10,
        )

    assert marker.read_text() == "important user data\n"
    assert set(output.iterdir()) == {marker}


def test_output_directory_may_not_overlap_a_run_root(tmp_path):
    run_root = tmp_path / "runs" / "controller"
    run_root.mkdir(parents=True)
    output = run_root / "publication"

    with pytest.raises(ValueError, match="overlap a controller run root"):
        release_report.write_report(
            [run_root],
            output,
            candidate_sha=CANDIDATE_SHA,
            reference_sha=REFERENCE_SHA,
            diagnostic=True,
            replicates=10,
        )

    assert not output.exists()
