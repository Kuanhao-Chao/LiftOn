from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pytest

from benchmarks.compare import (
    biology_canary,
    release_evaluation,
    whole_genome_assets,
    whole_genome_report,
    whole_genome_study,
)


HEADER = (
    "ref_mrna_id\ttool_feature_id\trecovered\tis_coding\tcopy_index\t"
    "protein_identity\tdna_identity\tref_prot_len\tlifted_prot_len\t"
    "ref_dna_len\tlifted_dna_len\tdna_basis\tstatus\tseqid\t"
    "n_cds_lifted\tn_cds_ref\tintron_chain_exact\texon_sn\texon_sp\t"
    "orf_start_ok\torf_stop_ok\torf_valid\n"
)
ROW = (
    "tx{index}\ttool-tx{index}\t1\t1\t0\t1.0\t1.0\t10\t10\t30\t30\t"
    "transcript\trecovered\tchr1\t1\t1\t1\t1.0\t1.0\t1\t1\t1\n"
)


def _write_json(path: Path, document: dict) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(document, indent=2, sort_keys=True) + "\n")
    return path


def _record(path: Path) -> dict:
    return whole_genome_study.file_record(path)


def _target_document(*, sequence: bool) -> dict:
    score = {"precision": 1.0, "recall": 1.0, "f1": 1.0}
    result = {
        "schema_version": 1,
        "method": "ortholog-scoped-target-coordinate-v1",
        "parameters": {
            "id_policy": "ortholog-map",
            "mapping_source_scope_validated": True,
        },
        "scope": {"gene_groups": 200, "transcript_groups": 200},
        "gene": {"locus": dict(score)},
        "transcript": {"locus": dict(score)},
        "structure": {
            name: dict(score) for name in ("intron_chain", "exon", "CDS")
        },
    }
    if sequence:
        result["target_sequence"] = {
            "exact_protein": {"rate": 1.0},
            "protein_identity": {
                "mean": 1.0,
                "coverage_weighted": 1.0,
                "mean_reciprocal_length_coverage": 1.0,
            },
            "orf": {
                "prediction_valid_rate": 1.0,
                "truth_valid_rate": 1.0,
                "both_valid_rate": 1.0,
            },
        }
    return result


def _arm(
    cell: Path,
    label: str,
    *,
    repetition: int,
    process_group_rss: float,
) -> dict:
    evaluation = cell / "evaluation"
    transcript_path = evaluation / f"{label}.transcripts.tsv"
    transcript_path.parent.mkdir(parents=True, exist_ok=True)
    transcript_path.write_text(
        HEADER + ROW.format(index=1) + ROW.format(index=2),
        encoding="utf-8",
    )
    target = _target_document(sequence=repetition == 1)
    target_path = _write_json(
        evaluation / f"{label}.target-annotation.json", target,
    )
    trace = cell / label / "process-group.jsonl"
    trace.parent.mkdir(parents=True, exist_ok=True)
    trace.write_text(json.dumps({
        "rss_mb": process_group_rss,
        "processes": 3,
    }) + "\n", encoding="utf-8")
    stable = {
        feature_type: {
            "applicable": True,
            "n_reference_ids": 2,
            "n_preserved_ids": 2,
            "preservation_rate": 1.0,
        }
        for feature_type in ("CDS", "exon")
    }
    return {
        "profile": {
            "exit_code": 0,
            "wall_clock_seconds": 8.0 if label == "candidate" else 10.0,
            "peak_rss_mb": process_group_rss / 2,
            "peak_process_group_rss_mb": process_group_rss,
            "peak_process_group_processes": 3,
            "process_group_sample_count": 1,
            "process_group_sample_interval_seconds": 1.0,
            "process_group_trace_path": str(trace.resolve()),
            "process_group_trace_sha256": hashlib.sha256(
                trace.read_bytes()
            ).hexdigest(),
        },
        "summary": {
            "n_reference_coding": 2,
            "completeness_coding": 1.0,
            "completeness_feature_total": 1.0,
            "stable_id_preservation": {"by_type": stable},
            "target_truth": target,
        },
        "evaluation_artifacts": {
            "transcripts_tsv": _record(transcript_path),
            "target_truth": _record(target_path),
        },
        "validation": {"is_valid": True, "n_errors": 0},
        "fingerprints": {
            "semantic_sha256": hashlib.sha256(
                f"{label}-semantic".encode()
            ).hexdigest(),
            "byte_sha256": hashlib.sha256(
                f"{label}-byte".encode()
            ).hexdigest(),
            "feature_counts": {"gene": 2, "mRNA": 2},
        },
    }


def _pair(root: Path, pair_id: str, repetition: int) -> Path:
    cell = root / pair_id / f"repetition-{repetition:02d}"
    document = {
        "schema_version": release_evaluation.SCHEMA_VERSION,
        "panel": "e2e",
        "benchmark": pair_id,
        "repetition": repetition,
        "order": (
            ["reference", "candidate"]
            if repetition % 2 else ["candidate", "reference"]
        ),
        "threads": 8,
        "modes": {"candidate": "fast", "reference": "fast"},
        "provenance": {
            "candidate": {
                "sha": "c623f0bddc5b5051a3670a3b0064c52cc61bb719"
            },
            "reference": {
                "sha": "e503643d8346c600fedabcd3a4dff5c0873a4a37"
            },
        },
        "versions": {
            "candidate": _arm(
                cell, "candidate", repetition=repetition,
                process_group_rss=100.0,
            ),
            "reference": _arm(
                cell, "reference", repetition=repetition,
                process_group_rss=120.0,
            ),
        },
    }
    return _write_json(cell / "pair_result.json", document)


def _evidence(tmp_path: Path) -> dict[str, Path]:
    study = tmp_path / "study.json"
    source_study = whole_genome_study.load_study(
        whole_genome_study.DEFAULT_STUDY
    )
    _write_json(study, source_study)
    preflight_document = {
        "schema_version": whole_genome_study.SCHEMA_VERSION,
        "kind": "lifton-v1.0.11-biology-study-preflight",
        "campaign_ready": True,
        "study": _record(study),
    }
    preflight_document["fingerprint"] = (
        whole_genome_study.canonical_sha256(preflight_document)
    )
    preflight = _write_json(tmp_path / "preflight.json", preflight_document)
    canary_document = {
        "schema_version": 1,
        "method": biology_canary.METHOD,
        "exact_sha": biology_canary.EXPECTED_SHA,
        "fast_supported": True,
        "rows": [
            {"benchmark": benchmark, "equivalent": True}
            for benchmark in biology_canary.EXPECTED_IDS
        ],
    }
    canary_document["fingerprint"] = whole_genome_study.canonical_sha256(
        canary_document
    )
    canary = _write_json(tmp_path / "canary.json", canary_document)
    campaign = tmp_path / "campaign"
    for pair_id in whole_genome_study.EXPECTED_PAIR_IDS:
        for repetition in range(1, 5):
            _pair(campaign, pair_id, repetition)
    return {
        "study": study,
        "preflight": preflight,
        "canary": canary,
        "campaign": campaign,
    }


def _reduce(paths: dict[str, Path]) -> dict:
    return whole_genome_report.reduce_campaign(
        paths["study"], paths["preflight"], paths["campaign"],
        paths["canary"],
    )


def test_reducer_accepts_complete_exact_release_cohort(tmp_path):
    result = _reduce(_evidence(tmp_path))

    assert result["coverage"] == {
        "planned_pairs": 28,
        "observed_pairs": 28,
        "planned_genome_transfers": 7,
        "successful_genome_transfers": 7,
        "incomplete": [],
    }
    assert result["qualification"]["status"] == "PASS"
    assert result["qualification"]["whole_release_claim"] == "DIAGNOSTIC ONLY"
    assert result["aggregate"]["performance_ratios"]["wall"]["estimate"] == (
        pytest.approx(0.8)
    )
    gate_names = {
        gate["name"] for gate in result["qualification"]["gates"]
    }
    assert "bee process-group RSS ratio" in gate_names
    assert (
        "aggregate process-group RSS ratio upper confidence bound"
        in gate_names
    )


def test_reducer_rejects_changed_target_artifact_or_memory_peak(tmp_path):
    paths = _evidence(tmp_path)
    pair_path = next(paths["campaign"].rglob("pair_result.json"))
    pair = json.loads(pair_path.read_text())
    artifact = Path(
        pair["versions"]["candidate"]["evaluation_artifacts"][
            "target_truth"
        ]["path"]
    )
    artifact.write_text("{}\n", encoding="utf-8")
    with pytest.raises(whole_genome_report.ReportError, match="artifact changed"):
        _reduce(paths)

    paths = _evidence(tmp_path / "memory")
    pair_path = next(paths["campaign"].rglob("pair_result.json"))
    pair = json.loads(pair_path.read_text())
    pair["versions"]["candidate"]["profile"][
        "peak_process_group_rss_mb"
    ] = 101.0
    _write_json(pair_path, pair)
    with pytest.raises(
        whole_genome_report.ReportError, match="trace metadata is invalid"
    ):
        _reduce(paths)


def test_reducer_rejects_incomplete_canary_or_wrong_release_sha(tmp_path):
    paths = _evidence(tmp_path)
    canary = json.loads(paths["canary"].read_text())
    canary["rows"].pop()
    canary.pop("fingerprint")
    canary["fingerprint"] = whole_genome_study.canonical_sha256(canary)
    _write_json(paths["canary"], canary)
    with pytest.raises(whole_genome_report.ReportError, match="coverage"):
        _reduce(paths)

    paths = _evidence(tmp_path / "sha")
    pair_path = next(paths["campaign"].rglob("pair_result.json"))
    pair = json.loads(pair_path.read_text())
    pair["provenance"]["candidate"]["sha"] = "f" * 40
    _write_json(pair_path, pair)
    with pytest.raises(whole_genome_report.ReportError, match="exact release"):
        _reduce(paths)


def test_reducer_marks_missing_repetition_as_failed_coverage(tmp_path):
    paths = _evidence(tmp_path)
    missing = sorted(paths["campaign"].rglob("pair_result.json"))[0]
    missing.unlink()

    result = _reduce(paths)

    assert result["coverage"]["observed_pairs"] == 24
    assert result["coverage"]["successful_genome_transfers"] == 6
    assert result["qualification"]["status"] == "FAIL"


def _with_sensitivity(metrics: dict, tmp_path: Path) -> tuple[dict, dict]:
    provider = _write_json(tmp_path / "provider-lock.json", {"locked": True})
    metrics = json.loads(json.dumps(metrics))
    metrics["provider_ortholog_lock"] = _record(provider)
    material = dict(metrics)
    material.pop("fingerprint")
    material.pop("created_at")
    metrics["fingerprint"] = whole_genome_study.canonical_sha256(material)
    score = {
        "gene": {"locus": {"precision": 1.0, "recall": 1.0, "f1": 1.0}},
        "transcript": {},
        "structure": {},
        "scope": {},
        "parameters": {},
    }
    sensitivity = {
        "schema_version": 1,
        "method": "released-target-annotation-sensitivity-v1",
        "created_at": "2026-08-13T00:00:00Z",
        "study": metrics["study"],
        "preflight": metrics["preflight"],
        "provider_ortholog_lock": metrics["provider_ortholog_lock"],
        "pairs": [
            {
                "id": pair_id,
                "primary_protein_rbh": {
                    str(threshold): {
                        arm: json.loads(json.dumps(score))
                        for arm in ("candidate", "reference")
                    }
                    for threshold in (0.25, 0.5, 0.8)
                },
                "provider_gene_relationships": {
                    "available": False,
                    "reason": "test",
                    "groups": 0,
                },
            }
            for pair_id in whole_genome_study.EXPECTED_PAIR_IDS
        ],
    }
    sensitivity["fingerprint"] = whole_genome_study.canonical_sha256({
        key: value for key, value in sensitivity.items()
        if key != "created_at"
    })
    return metrics, sensitivity


def test_report_assets_generate_expected_figures_and_tables(tmp_path):
    metrics = _reduce(_evidence(tmp_path / "evidence"))
    metrics, sensitivity = _with_sensitivity(metrics, tmp_path)

    paths = whole_genome_assets.generate_figures(metrics, tmp_path / "figures")

    assert {path.name for path in paths} == set(whole_genome_assets.FIGURES.values())
    assert all(path.stat().st_size > 10_000 for path in paths)
    report = "\n\n".join(
        "\n".join((
            whole_genome_assets.START.format(name=name),
            "stale",
            whole_genome_assets.END.format(name=name),
        ))
        for name in whole_genome_assets.BLOCKS
    ) + "\n"
    updated = whole_genome_assets.update_report(report, metrics, sensitivity)
    assert "stale" not in updated
    assert "Honey bee" in updated
    assert "Equal-transfer paired delta" in updated
    assert "Summed median-duration ratio" in updated
    assert "Two-concurrent-transfer v1.0.11 memory proxy" in updated
    for excluded in (
        "arabidopsis", "rice", "zebrafish", "chicken", "xenopus",
    ):
        assert excluded not in updated.lower()
    report_path = tmp_path / "report.mdx"
    report_path.write_text(updated)
    assert whole_genome_assets.check_report(
        report_path, metrics, sensitivity,
    ) == (True, "")


def test_report_assets_reject_metric_or_sensitivity_drift(tmp_path):
    metrics = _reduce(_evidence(tmp_path / "evidence"))
    metrics, sensitivity = _with_sensitivity(metrics, tmp_path)
    changed = json.loads(json.dumps(metrics))
    changed["pairs"][0]["public_label"] = "changed"
    with pytest.raises(whole_genome_assets.AssetError, match="fingerprint"):
        whole_genome_assets.validate_metrics(changed)

    changed = json.loads(json.dumps(metrics))
    changed["qualification"]["whole_release_claim"] = "PASS"
    material = dict(changed)
    material.pop("fingerprint")
    material.pop("created_at")
    changed["fingerprint"] = whole_genome_study.canonical_sha256(material)
    with pytest.raises(
        whole_genome_assets.AssetError, match="qualification boundary"
    ):
        whole_genome_assets.validate_metrics(changed)

    sensitivity["pairs"].pop()
    with pytest.raises(whole_genome_assets.AssetError, match="fingerprint"):
        whole_genome_assets.validate_sensitivity(sensitivity, metrics)


def _transcripts_tsv(path: Path, *, total: int, recovered: int) -> Path:
    """A per-transcript table with a realistic, non-round recovery fraction."""

    path.parent.mkdir(parents=True, exist_ok=True)
    rows = [HEADER]
    for index in range(total):
        row = ROW.format(index=index)
        if index >= recovered:
            columns = row.rstrip("\n").split("\t")
            columns[2] = "0"          # recovered
            columns[5] = "0.0"        # protein_identity
            columns[12] = "map_failed"
            row = "\t".join(columns) + "\n"
        rows.append(row)
    path.write_text("".join(rows))
    return path


def _version_for_completeness(tmp_path: Path, *, total: int, recovered: int) -> dict:
    tsv = _transcripts_tsv(tmp_path / "candidate.transcripts.tsv",
                           total=total, recovered=recovered)
    return {
        "summary": {
            "n_reference_coding": total,
            "n_recovered_coding": recovered,
            # The neutral evaluator rounds this to five decimals, exactly as it
            # does on real data.
            "completeness_coding": round(recovered / total, 5),
            "completeness_feature_total": 0.94570,
        },
        "evaluation_artifacts": {"transcripts_tsv": _record(tsv)},
    }


def test_rounded_completeness_agrees_with_its_own_transcript_table(tmp_path):
    """The stored ratio is rounded; the cross-check must allow for that.

    Comparing a five-decimal value against a full-precision recomputation at
    1e-9 cannot hold on any real cohort -- 23365/23471 stores as 0.99548 but
    recomputes as 0.9954837885. The synthetic fixtures hid this because their
    completeness is exactly 1.0.
    """

    version = _version_for_completeness(tmp_path, total=23471, recovered=23365)

    metrics = whole_genome_report._transcript_metrics(version, "bee.candidate")

    assert metrics["n_recovered_coding"] == 23365
    assert metrics["n_reference_coding"] == 23471
    assert metrics["completeness_coding"] == 0.99548


def test_recovered_count_disagreement_still_fails_closed(tmp_path):
    version = _version_for_completeness(tmp_path, total=23471, recovered=23365)
    version["summary"]["n_recovered_coding"] = 23364

    with pytest.raises(
        whole_genome_report.ReportError, match="recovered coding count disagrees"
    ):
        whole_genome_report._transcript_metrics(version, "bee.candidate")


def test_completeness_inconsistent_with_the_table_still_fails_closed(tmp_path):
    version = _version_for_completeness(tmp_path, total=23471, recovered=23365)
    version["summary"].pop("n_recovered_coding")
    version["summary"]["completeness_coding"] = 0.5

    with pytest.raises(
        whole_genome_report.ReportError, match="coding completeness disagrees"
    ):
        whole_genome_report._transcript_metrics(version, "bee.candidate")
