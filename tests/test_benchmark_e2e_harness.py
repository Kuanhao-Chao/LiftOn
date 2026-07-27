"""Focused regression tests for the Phase-16 biological E2E harness."""
from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest


_HARNESS_DIR = Path(__file__).resolve().parent.parent / "benchmarks"
sys.path.insert(0, str(_HARNESS_DIR))

import run_benchmarks as harness  # noqa: E402


def _profile(exit_code: int = 0) -> harness.ProfileResult:
    return harness.ProfileResult(
        wall_clock_seconds=2.5,
        peak_rss_mb=32.0,
        user_cpu_seconds=1.0,
        sys_cpu_seconds=0.5,
        exit_code=exit_code,
        stdout_path="stdout.log",
        stderr_path="stderr.log",
        time_log_path="time.log",
    )


@pytest.fixture(autouse=True)
def _stable_lift_provenance(monkeypatch):
    """Keep unit tests independent of the host's installed LiftOn stack."""

    monkeypatch.setattr(
        harness,
        "_lift_tool_provenance",
        lambda: {
            "tools": {"lifton": {"sha256": "1" * 64}},
            "lifton_source": {"sha256": "2" * 64},
            "distributions": {"duckdb": {"version": "test"}},
        },
    )


def test_registry_strips_inline_underscore_metadata(tmp_path):
    registry = tmp_path / "datasets.json"
    registry.write_text(json.dumps({
        "datasets": [
            {"_note": "comment-only entry"},
            {
                "_note": "inline metadata",
                "id": "demo",
                "species": "Demo",
                "reference_fa": "ref.fa",
                "target_fa": "target.fa",
                "reference_gff": "ref.gff3",
            },
        ],
        "lifton_flags": ["--stream"],
        "evaluation_flags": ["-E"],
    }))

    loaded = harness.load_registry(registry)

    assert [dataset.id for dataset in loaded.datasets] == ["demo"]
    assert loaded.datasets[0].species == "Demo"


def test_parse_score_txt_uses_the_real_eight_column_schema(tmp_path):
    score = tmp_path / "score.txt"
    score.write_text(
        "tx1\t1.0\t0.9\t0.98\t0.95\tLiftoff\tsynonymous\tchr1:1-9\n"
        "tx2\t0\t0\t0.80\t0\tLiftoff\tnon_coding\tchr1:20-30\n"
        "legacy\tchr1\t1\t9\tmapped\t0.8\t0.7\n"
        "bad_identity\t1\t1\t0.9\tnan\tLiftoff\tsynonymous\tchr1:40-50\n"
    )

    summary = harness.parse_score_txt(score)

    assert summary.format == "lifton_score_v1"
    assert summary.mapped == 2
    assert summary.records == 2
    assert summary.coding_records == 1
    assert summary.noncoding_records == 1
    assert summary.malformed_records == 2
    assert summary.avg_identity == pytest.approx(0.95)
    assert summary.lost is None


def test_parse_eval_txt_uses_the_real_five_column_schema(tmp_path):
    evaluation = tmp_path / "eval.txt"
    evaluation.write_text(
        "tx1\t0.97\t0.91\tsynonymous\tchr1:1-9\n"
        "tx2\t0.88\t0\tnon_coding\tchr1:20-30\n"
        "wrong\tcolumn\tcount\n"
        "bad\t0.8\tinf\tsynonymous\tchr1:40-50\n"
    )

    summary = harness.parse_eval_txt(evaluation)

    assert summary.format == "lifton_eval_v1"
    assert summary.mapped == 2
    assert summary.records == 2
    assert summary.coding_records == 1
    assert summary.noncoding_records == 1
    assert summary.malformed_records == 2
    assert summary.avg_identity == pytest.approx(0.91)
    assert summary.lost is None


def test_no_protein_status_remains_coding_in_score_and_evaluation(tmp_path):
    score = tmp_path / "score.txt"
    score.write_text(
        "tx1\t1\t0\t0.8\t0\tLiftoff\tno_protein\tchr1:1-9\n"
        "tx2\t0\t0\t0.7\t0\tLiftoff\tnon_coding\tchr1:20-30\n"
    )
    evaluation = tmp_path / "eval.txt"
    evaluation.write_text(
        "tx1\t0.8\t0\tno_protein\tchr1:1-9\n"
        "tx2\t0.7\t0\tnon_coding\tchr1:20-30\n"
    )

    score_summary = harness.parse_score_txt(score)
    evaluation_summary = harness.parse_eval_txt(evaluation)

    for summary in (score_summary, evaluation_summary):
        assert summary.records == 2
        assert summary.coding_records == 1
        assert summary.noncoding_records == 1
        assert summary.avg_identity == 0.0


def _valid_biological_result() -> dict:
    return {
        "biological_summary": {
            "schema_version": 1,
            "reference_features": 3,
            "mapped_features": 2,
            "lost_features": 1,
            "extra_copy_features": 1,
            "feature_completeness": 2 / 3,
            "emitted_transcript_records": 2,
            "mapped_transcripts_reported": 1,
            "evaluated_transcript_records": 2,
            "evaluated_coding_records": 1,
            "mean_protein_identity": 0.9,
        },
        "score_summary": {
            "format": "lifton_score_v1",
            "records": 2,
            "coding_records": 1,
            "noncoding_records": 1,
            "malformed_records": 0,
        },
        "evaluation_summary": {
            "format": "lifton_eval_v1",
            "records": 2,
            "coding_records": 1,
            "noncoding_records": 1,
            "malformed_records": 0,
        },
    }


def test_strict_biological_contract_accepts_consistent_results():
    errors = harness.validate_biological_result(
        _valid_biological_result(),
        require_evaluation=True,
        require_identity=True,
    )

    assert errors == []


@pytest.mark.parametrize(
    ("section", "field", "value", "message"),
    [
        ("biological_summary", "lost_features", 2, "does not equal"),
        ("biological_summary", "feature_completeness", 0.5, "inconsistent"),
        ("biological_summary", "mean_protein_identity", float("inf"), "finite"),
        ("score_summary", "malformed_records", 1, "not zero"),
        ("evaluation_summary", "records", 0, "positive integer"),
        ("evaluation_summary", "coding_records", 2, "disagrees"),
    ],
)
def test_strict_biological_contract_rejects_silent_bad_results(
        section, field, value, message):
    result = _valid_biological_result()
    result[section][field] = value

    errors = harness.validate_biological_result(
        result,
        require_evaluation=True,
        require_identity=True,
    )

    assert any(message in error for error in errors)


def _prepare_inputs(data_dir: Path) -> tuple[harness.Dataset, dict[str, bytes]]:
    dataset = harness.Dataset(
        id="demo",
        species="Demo species",
        reference_fa="https://example.test/ref.fa",
        target_fa="https://example.test/target.fa",
        reference_gff="https://example.test/ref.gff3",
        target_gff="https://example.test/published.gff3",
    )
    payloads = {
        "ref.fa": b">ref\nACGT\n",
        "target.fa": b">target\nACGT\n",
        "ref.gff3": b"##gff-version 3\n",
        "published.gff3": b"cached comparison annotation; do not modify\n",
    }
    dataset_dir = data_dir / dataset.id
    dataset_dir.mkdir(parents=True)
    for name, payload in payloads.items():
        (dataset_dir / name).write_bytes(payload)
    return dataset, payloads


def _write_lift_artifacts(out_gff: Path) -> None:
    out_gff.write_text(
        "##gff-version 3\n"
        "chr1\tLiftOn\tgene\t1\t9\t.\t+\t.\tID=gene1\n"
    )
    output_dir = out_gff.parent / "lifton_output"
    stats_dir = output_dir / "stats"
    stats_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "score.txt").write_text(
        "tx1\t1\t1\t0.9\t0.8\tLiftoff\tsynonymous\tchr1:1-9\n"
        "tx2\t0\t0\t0.7\t0\tLiftoff\tnon_coding\tchr1:20-30\n"
    )
    (stats_dir / "mapped_feature.txt").write_text(
        "gene1\t1\tcoding\n"
        "gene2\t1\tnon-coding\n"
    )
    (stats_dir / "unmapped_features.txt").write_text("gene3\tcoding\n")
    (stats_dir / "extra_copy_features.txt").write_text("gene2\t2\tnon-coding\n")
    (stats_dir / "mapped_transcript.txt").write_text("tx1\t1\tcoding\n")


def test_run_dataset_evaluates_isolated_candidate_and_returns_biology(
        tmp_path, monkeypatch):
    data_dir = tmp_path / "data"
    results_dir = tmp_path / "results"
    dataset, payloads = _prepare_inputs(data_dir)
    calls = []

    def fake_profile(argv, *, label, log_dir, log):
        calls.append((label, list(argv)))
        output = Path(argv[argv.index("-o") + 1])
        if label == "lift":
            _write_lift_artifacts(output)
            return _profile()
        assert label == "evaluation"
        assert "-E" in argv
        assert output.name == "candidate.gff3"
        assert output.parent == results_dir / "demo" / "evaluation"
        assert output.read_bytes() == (
            results_dir / "demo" / "lifton.gff3"
        ).read_bytes()
        eval_outdir = output.parent / "lifton_output"
        eval_outdir.mkdir(parents=True, exist_ok=True)
        (eval_outdir / "eval.txt").write_text(
            "tx1\t0.99\t0.90\tsynonymous\tchr1:1-9\n"
            "tx2\t0.80\t0\tnon_coding\tchr1:20-30\n"
        )
        (eval_outdir / "run_manifest.json").write_text("{}\n")
        return _profile()

    monkeypatch.setattr(harness, "run_profiled", fake_profile)

    result = harness.run_dataset(
        dataset,
        data_dir=data_dir,
        results_dir=results_dir,
        lifton_flags=["--stream"],
        evaluation_flags=["-E"],
        do_download=False,
        do_lift=True,
        do_evaluate=True,
        force=True,
        log=lambda *_args, **_kwargs: None,
    )

    assert "error" not in result
    assert [label for label, _argv in calls] == ["lift", "evaluation"]
    assert (
        data_dir / "demo" / "published.gff3"
    ).read_bytes() == payloads["published.gff3"]
    eval_argv = calls[1][1]
    assert eval_argv[eval_argv.index("-o") + 1] != str(
        data_dir / "demo" / "published.gff3"
    )
    assert result["score_summary"]["records"] == 2
    assert result["evaluation_summary"]["records"] == 2
    assert result["evaluation_summary"]["coding_records"] == 1
    biological = result["biological_summary"]
    assert biological["reference_features"] == 3
    assert biological["mapped_features"] == 2
    assert biological["lost_features"] == 1
    assert biological["extra_copy_features"] == 1
    assert biological["feature_completeness"] == pytest.approx(2 / 3)
    assert biological["evaluated_transcript_records"] == 2
    assert biological["mean_protein_identity"] == pytest.approx(0.90)
    assert result["eval_summary"] == {
        "mapped": 2,
        "lost": 1,
        "extra_copies": 1,
        "avg_identity": pytest.approx(0.90),
        "score_file": str(
            results_dir / "demo" / "evaluation" / "lifton_output" / "eval.txt"
        ),
    }


def test_run_dataset_evaluates_candidate_without_published_target_gff(
        tmp_path, monkeypatch):
    data_dir = tmp_path / "data"
    results_dir = tmp_path / "results"
    dataset, _payloads = _prepare_inputs(data_dir)
    dataset.target_gff = None
    calls = []

    def fake_profile(argv, *, label, log_dir, log):
        calls.append(label)
        output = Path(argv[argv.index("-o") + 1])
        if label == "lift":
            _write_lift_artifacts(output)
        else:
            eval_outdir = output.parent / "lifton_output"
            eval_outdir.mkdir(parents=True, exist_ok=True)
            (eval_outdir / "eval.txt").write_text(
                "tx1\t0.99\t0.90\tsynonymous\tchr1:1-9\n"
                "tx2\t0.80\t0\tnon_coding\tchr1:20-30\n"
            )
            (eval_outdir / "run_manifest.json").write_text("{}\n")
        return _profile()

    monkeypatch.setattr(harness, "run_profiled", fake_profile)

    result = harness.run_dataset(
        dataset,
        data_dir=data_dir,
        results_dir=results_dir,
        lifton_flags=[],
        evaluation_flags=["-E"],
        do_download=False,
        do_lift=True,
        do_evaluate=True,
        force=True,
        log=lambda *_args, **_kwargs: None,
    )

    assert calls == ["lift", "evaluation"]
    assert "error" not in result
    assert result["evaluation_summary"]["records"] == 2
    assert result["input_paths"].get("target_gff") is None


def test_run_dataset_propagates_nonzero_evaluation(tmp_path, monkeypatch):
    data_dir = tmp_path / "data"
    results_dir = tmp_path / "results"
    dataset, payloads = _prepare_inputs(data_dir)

    def fake_profile(argv, *, label, log_dir, log):
        output = Path(argv[argv.index("-o") + 1])
        if label == "lift":
            _write_lift_artifacts(output)
            return _profile()
        return _profile(exit_code=9)

    monkeypatch.setattr(harness, "run_profiled", fake_profile)

    result = harness.run_dataset(
        dataset,
        data_dir=data_dir,
        results_dir=results_dir,
        lifton_flags=[],
        evaluation_flags=["-E"],
        do_download=False,
        do_lift=True,
        do_evaluate=True,
        force=True,
        log=lambda *_args, **_kwargs: None,
    )

    assert result["error"] == "evaluation non-zero exit"
    assert result["eval_profile"]["exit_code"] == 9
    assert (
        data_dir / "demo" / "published.gff3"
    ).read_bytes() == payloads["published.gff3"]


def test_lift_cache_reuses_only_exact_key_and_real_profile(
        tmp_path, monkeypatch):
    data_dir = tmp_path / "data"
    results_dir = tmp_path / "results"
    dataset, _payloads = _prepare_inputs(data_dir)
    calls = []
    provenance = {"source": "source-a", "duckdb": "1.5.2"}

    def fake_profile(argv, *, label, log_dir, log):
        calls.append((label, list(argv)))
        _write_lift_artifacts(Path(argv[argv.index("-o") + 1]))
        return _profile()

    monkeypatch.setattr(harness, "run_profiled", fake_profile)
    monkeypatch.setattr(
        harness,
        "_lift_tool_provenance",
        lambda: {
            "source": provenance["source"],
            "distributions": {"duckdb": provenance["duckdb"]},
        },
    )
    arguments = {
        "data_dir": data_dir,
        "results_dir": results_dir,
        "lifton_flags": ["--stream"],
        "evaluation_flags": [],
        "do_download": False,
        "do_lift": True,
        "do_evaluate": False,
        "force": False,
        "log": lambda *_args, **_kwargs: None,
    }

    first = harness.run_dataset(dataset, **arguments)
    second = harness.run_dataset(dataset, **arguments)

    assert [label for label, _argv in calls] == ["lift"]
    assert first["lift_profile"] == second["lift_profile"]
    assert second["lift_profile"]["wall_clock_seconds"] == 2.5
    assert second["lift_profile"]["peak_rss_mb"] == 32.0
    sentinel = json.loads(
        (results_dir / "demo" / ".lifton.done").read_text()
    )
    assert sentinel["kind"] == "lifton_lift_cache"
    assert sentinel["profile"]["wall_clock_seconds"] == 2.5

    target = data_dir / "demo" / "target.fa"
    target.write_bytes(target.read_bytes().replace(b"ACGT", b"TGCA"))
    harness.run_dataset(dataset, **arguments)
    provenance["source"] = "source-b"
    harness.run_dataset(dataset, **arguments)
    provenance["duckdb"] = "1.5.3"
    harness.run_dataset(dataset, **arguments)
    changed_flags = {**arguments, "lifton_flags": ["--native"]}
    harness.run_dataset(dataset, **changed_flags)

    assert [label for label, _argv in calls] == [
        "lift", "lift", "lift", "lift", "lift",
    ]


def test_legacy_or_malformed_lift_cache_is_stale(tmp_path, monkeypatch):
    data_dir = tmp_path / "data"
    results_dir = tmp_path / "results"
    dataset, _payloads = _prepare_inputs(data_dir)
    dataset_results = results_dir / "demo"
    dataset_results.mkdir(parents=True)
    sentinel = dataset_results / ".lifton.done"
    sentinel.write_text("2026-07-17T00:00:00Z")
    calls = []

    def fake_profile(argv, *, label, log_dir, log):
        calls.append(label)
        _write_lift_artifacts(Path(argv[argv.index("-o") + 1]))
        return _profile()

    monkeypatch.setattr(harness, "run_profiled", fake_profile)
    monkeypatch.setattr(
        harness, "_lift_tool_provenance", lambda: {"fingerprint": "fixed"},
    )

    result = harness.run_dataset(
        dataset,
        data_dir=data_dir,
        results_dir=results_dir,
        lifton_flags=[],
        evaluation_flags=[],
        do_download=False,
        do_lift=True,
        do_evaluate=False,
        force=False,
        log=lambda *_args, **_kwargs: None,
    )

    assert calls == ["lift"]
    assert result["lift_profile"]["wall_clock_seconds"] == 2.5
    assert json.loads(sentinel.read_text())["kind"] == "lifton_lift_cache"


def test_lift_cache_refuses_zero_or_synthetic_profile(tmp_path):
    output = tmp_path / "lifton.gff3"
    output.write_text("##gff-version 3\n")
    sentinel = tmp_path / ".lifton.done"
    profile = _profile()
    profile.wall_clock_seconds = 0.0

    with pytest.raises(ValueError, match="zero"):
        harness._write_lift_cache(
            sentinel,
            request={"argv": ["lifton"], "inputs": {}, "provenance": {}},
            profile=profile,
            out_gff=output,
        )

    assert not sentinel.exists()
