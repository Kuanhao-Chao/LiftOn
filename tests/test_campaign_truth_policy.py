from __future__ import annotations

import json
from pathlib import Path

from benchmarks.compare import release_evaluation, release_report


def _prf(value: float = 0.8):
    return {
        "true_positive": 8,
        "false_positive": 2,
        "false_negative": 2,
        "predicted": 10,
        "expected": 10,
        "precision": value,
        "recall": value,
        "f1": value,
    }


def _truth(value: float = 0.8):
    level = {
        "locus": _prf(value),
        "strand": _prf(value),
        "copy": _prf(value),
        "copy_count_exact": {
            "groups_exact": 8,
            "groups_total": 10,
            "rate": value,
        },
    }
    return {
        "schema_version": 1,
        "method": "ortholog-scoped-target-coordinate-v1",
        "parameters": {
            "minimum_reciprocal_overlap": 0.5,
            "id_policy": "ortholog-map",
            "mapping_required": True,
            "mapping_requirement_satisfied": True,
        },
        "scope": {
            "gene_groups": 10,
            "transcript_groups": 10,
            "mapping_entries": 20,
            "gene": {
                "groups": 10,
                "predicted_scored": 10,
                "expected_scored": 10,
                "prediction_models_total": 10,
                "truth_models_total": 10,
                "prediction_models_ignored": 0,
                "truth_models_ignored": 0,
            },
            "transcript": {
                "groups": 10,
                "predicted_scored": 10,
                "expected_scored": 10,
                "prediction_models_total": 10,
                "truth_models_total": 10,
                "prediction_models_ignored": 0,
                "truth_models_ignored": 0,
            },
        },
        "gene": level,
        "transcript": level,
        "structure": {
            "intron_chain": _prf(value),
            "intron": _prf(value),
            "exon": _prf(value),
            "CDS": _prf(value),
        },
    }


def test_target_truth_artifact_is_checked_against_pair_and_controller(tmp_path):
    pair_path = tmp_path / "pair_result.json"
    expected = tmp_path / "evaluation" / "candidate.target_truth.json"
    expected.parent.mkdir()
    metrics = _truth()
    expected.write_text(json.dumps(metrics))
    record = {
        "path": str(expected.resolve()),
        "size": expected.stat().st_size,
        "sha256": release_report.sha256_file(expected),
    }
    pair = {
        "versions": {
            "candidate": {
                "summary": {"target_truth": metrics},
                "evaluation_artifacts": {"target_truth": record},
            },
        },
    }
    controller = {
        **record,
        "mtime_ns": expected.stat().st_mtime_ns,
    }

    assert release_report._target_truth_summary(
        pair_path, pair, "candidate", controller,
    ) == metrics

    expected.write_text(json.dumps(_truth(0.7)))
    try:
        release_report._target_truth_summary(
            pair_path, pair, "candidate", controller,
        )
    except ValueError as error:
        assert "sealed evidence" in str(error)
    else:  # pragma: no cover - cryptographic tampering must fail closed
        raise AssertionError("target-truth tampering was accepted")


def test_target_truth_policy_adds_absolute_and_nonregression_gates():
    candidate = _truth(0.8)
    reference = _truth(0.805)
    metrics = {
        "publication": {
            "mode": "release",
            "release_evidence_valid": True,
        },
        "panels": {},
        "cells": [{
            "panel": "subset",
            "benchmark": "truth-demo",
            "campaign_policy": {
                "truth_policy": "target_truth_required",
            },
            "target_truth": {
                "candidate": {
                    "available": True,
                    "deterministic": True,
                    "summary": candidate,
                },
                "reference": {
                    "available": True,
                    "deterministic": True,
                    "summary": reference,
                },
                "locus_f1_delta": {
                    "gene": -0.005,
                    "transcript": -0.005,
                },
            },
            "candidate_valid": True,
            "reference_valid": True,
            "candidate_quality_deterministic": True,
            "reference_quality_deterministic": True,
            "common_pi_deterministic": True,
            "candidate_byte_deterministic": True,
            "reference_byte_deterministic": True,
            "absolute_quality": {
                "candidate": {},
                "reference": {},
            },
            "quality_deltas": {},
            "common_pi": None,
        }],
    }

    verdict = release_report.evaluate_verdict(metrics)
    checks = {
        row["name"]: row
        for row in verdict["checks"]
        if "target_truth" in row["name"]
    }

    assert checks[
        "subset.truth-demo.candidate.target_truth_gene_locus_f1"
    ]["passed"] is True
    assert checks[
        "subset.truth-demo.candidate.target_truth_gene_scope_groups"
    ]["passed"] is True
    assert checks[
        "subset.truth-demo.target_truth_gene_locus_f1_delta"
    ]["passed"] is True
    assert all(row["passed"] for row in checks.values())


def test_target_truth_policy_rejects_trivially_small_scope():
    candidate = _truth(1.0)
    reference = _truth(1.0)
    for metrics in (candidate, reference):
        metrics["scope"]["gene_groups"] = 1
        metrics["scope"]["gene"]["groups"] = 1
    document = {
        "publication": {
            "mode": "release",
            "release_evidence_valid": True,
        },
        "panels": {},
        "cells": [{
            "panel": "subset",
            "benchmark": "tiny-truth",
            "campaign_policy": {
                "truth_policy": "target_truth_required",
            },
            "target_truth": {
                "candidate": {
                    "available": True,
                    "deterministic": True,
                    "summary": candidate,
                },
                "reference": {
                    "available": True,
                    "deterministic": True,
                    "summary": reference,
                },
                "locus_f1_delta": {"gene": 0.0, "transcript": 0.0},
            },
            "candidate_valid": True,
            "reference_valid": True,
            "candidate_quality_deterministic": True,
            "reference_quality_deterministic": True,
            "common_pi_deterministic": True,
            "candidate_byte_deterministic": True,
            "reference_byte_deterministic": True,
            "absolute_quality": {"candidate": {}, "reference": {}},
            "quality_deltas": {},
            "common_pi": None,
        }],
    }

    checks = {
        row["name"]: row
        for row in release_report.evaluate_verdict(document)["checks"]
    }
    assert checks[
        "subset.tiny-truth.candidate.target_truth_gene_scope_groups"
    ]["passed"] is False


def test_full_raw_input_mode_does_not_require_historical_lmtp_cache(
        tmp_path, monkeypatch):
    ref_fa = tmp_path / "ref.fa"
    ref_gff = tmp_path / "ref.gff3"
    target_fa = tmp_path / "target.fa"
    ref_fa.write_text(">chr1\nAAAA\n")
    ref_gff.write_text(
        "##gff-version 3\n"
        "chr1\tt\tgene\t1\t4\t.\t+\t.\tID=g1\n"
    )
    target_fa.write_text(">chr1\nAAAA\n")
    registry = tmp_path / "benchmarks.json"
    registry.write_text(json.dumps({
        "benchmarks": [{
            "id": "v2_raw",
            "species": "Synthetic raw full",
            "cross_species": False,
            "annotation_database": "RefSeq",
            "ref_genome": str(ref_fa),
            "ref_gff": str(ref_gff),
            "tgt_genome": str(target_fa),
            "full_input_mode": "raw",
        }],
    }))
    monkeypatch.setattr(
        release_evaluation.devel_refresh,
        "_cached_inputs",
        lambda _benchmark: (_ for _ in ()).throw(
            AssertionError("raw mode consulted historical cache")
        ),
    )
    monkeypatch.setattr(
        release_evaluation,
        "source_cli_options",
        lambda _source: frozenset(),
    )

    inputs = release_evaluation.resolve_panel_inputs(
        "full",
        "v2_raw",
        benchmark_registry=registry,
    )
    assert inputs.liftoff_gff is None
    assert inputs.miniprot_gff is None
    assert inputs.transcripts_fa is None
    assert inputs.proteins_fa is None
    argv = release_evaluation.build_lifton_argv(
        release_evaluation.SourceSpec(
            "candidate", tmp_path, "a" * 40, Path("/env/bin/lifton"),
        ),
        inputs,
        tmp_path / "out.gff3",
        threads=2,
    )
    assert not {"-L", "-M", "-T", "-P"} & set(argv)
