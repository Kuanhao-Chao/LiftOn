from __future__ import annotations

import json
from pathlib import Path

import pytest

from benchmarks.compare import release_evaluation as release


def _annotation(gene_id="g1", transcript_id="t1"):
    return (
        "##gff-version 3\n"
        f"chr1\tLiftOn\tgene\t1\t90\t.\t+\t.\tID={gene_id}\n"
        f"chr1\tLiftOn\tmRNA\t1\t90\t.\t+\t."
        f"\tID={transcript_id};Parent={gene_id}\n"
        f"chr1\tLiftOn\texon\t1\t30\t.\t+\t."
        f"\tID={transcript_id}.e1;Parent={transcript_id}\n"
        f"chr1\tLiftOn\texon\t61\t90\t.\t+\t."
        f"\tID={transcript_id}.e2;Parent={transcript_id}\n"
        f"chr1\tLiftOn\tCDS\t1\t30\t.\t+\t0"
        f"\tID={transcript_id}.c;Parent={transcript_id}\n"
        f"chr1\tLiftOn\tCDS\t61\t90\t.\t+\t0"
        f"\tID={transcript_id}.c;Parent={transcript_id}\n"
    )


def _inputs(tmp_path: Path) -> release.PanelInputs:
    tmp_path.mkdir(parents=True, exist_ok=True)
    ref_fa = tmp_path / "ref.fa"
    tgt_fa = tmp_path / "target.fa"
    ref_gff = tmp_path / "ref.gff3"
    truth_gff = tmp_path / "truth.gff3"
    ortholog_map = tmp_path / "ortholog.json"
    ref_fa.write_text(">chr1\n" + "A" * 100 + "\n")
    tgt_fa.write_text(">chr1\n" + "A" * 100 + "\n")
    ref_gff.write_text(_annotation())
    truth_gff.write_text(_annotation("truth-gene", "truth-tx"))
    ortholog_map.write_text(json.dumps({"mappings": [
        {
            "source_id": "g1",
            "truth_ids": ["truth-gene"],
            "feature_type": "gene",
        },
        {
            "source_id": "t1",
            "truth_ids": ["truth-tx"],
            "feature_type": "transcript",
        },
    ]}))
    return release.PanelInputs(
        benchmark="truth-demo",
        panel="subset",
        species="Example",
        cross_species=True,
        annotation_database="RefSeq",
        ref_fa=ref_fa,
        ref_gff=ref_gff,
        tgt_fa=tgt_fa,
        truth_gff=truth_gff,
        ortholog_map=ortholog_map,
    )


def test_release_validation_document_includes_target_fasta_gate(tmp_path):
    fasta = tmp_path / "target.fa"
    output = tmp_path / "output.gff3"
    fasta.write_text(">chr1\n" + "A" * 80 + "\n")
    output.write_text(_annotation())

    result = release._validation_document(output, fasta)

    assert result["is_valid"] is False
    assert result["passes"]["full_semantic"]["is_valid"] is True
    assert result["passes"]["streaming_structure"]["is_valid"] is True
    assert result["passes"]["target_fasta_bounds"]["is_valid"] is False
    assert "target_coordinate_out_of_bounds" in {
        issue["check"] for issue in result["issues"]
        if issue["validator"] == "target_fasta_bounds"
    }


def test_release_validation_rejects_stale_sequence_region(
        tmp_path):
    fasta = tmp_path / "target.fa"
    output = tmp_path / "output.gff3"
    fasta.write_text(">chr1\n" + "A" * 100 + "\n")
    output.write_text(
        _annotation().replace(
            "##gff-version 3\n",
            "##gff-version 3\n##sequence-region chr1 1 80\n",
        )
    )

    result = release._validation_document(output, fasta)
    target_pass = result["passes"]["target_fasta_bounds"]

    assert result["is_valid"] is False
    assert target_pass["n_errors"] == 4
    assert target_pass["n_warnings"] == 0
    assert target_pass["issue_totals"] == {
        "sequence_region_containment": 4,
    }
    assert {
        issue["check"] for issue in result["issues"]
        if issue["validator"] == "target_fasta_bounds"
    } == {"sequence_region_containment"}


def test_truth_inputs_are_not_lifton_execution_inputs(tmp_path, monkeypatch):
    monkeypatch.setattr(
        release, "source_cli_options", lambda _source: frozenset(),
    )
    inputs = _inputs(tmp_path)
    argv = release.build_lifton_argv(
        release.SourceSpec(
            "candidate", tmp_path, "a" * 40, Path("/env/bin/lifton"),
        ),
        inputs,
        tmp_path / "output.gff3",
        threads=2,
    )

    assert set(inputs.evaluation_paths()) == {"truth_gff", "ortholog_map"}
    assert "truth_gff" not in inputs.required_paths()
    assert str(inputs.truth_gff) not in argv
    assert str(inputs.ortholog_map) not in argv


def test_registry_truth_configuration_fails_closed_without_mapping(tmp_path):
    truth = tmp_path / "truth.gff3"
    truth.write_text(_annotation())

    with pytest.raises(ValueError, match="requires a non-empty ortholog_map"):
        release._evaluation_input_paths(
            {"target_truth": {"gff": truth.name}},
            tmp_path / "benchmarks.json",
        )

    config = release._evaluation_input_paths(
        {
            "target_truth": {
                "gff": truth.name,
                "id_policy": "exact-id",
            },
        },
        tmp_path / "benchmarks.json",
    )
    assert config == {
        "truth_gff": truth.resolve(),
        "truth_id_policy": "exact-id",
    }


def test_panel_specific_truth_is_exact_and_never_uses_legacy_fallback(
        tmp_path):
    paths = {
        name: tmp_path / name
        for name in (
            "legacy.gff3",
            "subset.gff3",
            "subset.map.json",
            "full.gff3",
            "full.map.json",
        )
    }
    for path in paths.values():
        path.write_text(path.name + "\n")
    bench = {
        "target_truth": {
            "gff": paths["legacy.gff3"].name,
            "id_policy": "exact-id",
        },
        "target_truth_by_panel": {
            "subset": {
                "gff": paths["subset.gff3"].name,
                "ortholog_map": paths["subset.map.json"].name,
                "id_policy": "ortholog-map",
            },
            "full": {
                "gff": paths["full.gff3"].name,
                "ortholog_map": paths["full.map.json"].name,
                "id_policy": "ortholog-map",
            },
        },
    }

    subset = release._evaluation_input_paths(
        bench, tmp_path / "benchmarks.json", "subset",
    )
    full = release._evaluation_input_paths(
        bench, tmp_path / "benchmarks.json", "full",
    )

    assert subset == {
        "truth_gff": paths["subset.gff3"].resolve(),
        "ortholog_map": paths["subset.map.json"].resolve(),
        "truth_id_policy": "ortholog-map",
    }
    assert full["truth_gff"] == paths["full.gff3"].resolve()
    assert full["ortholog_map"] == paths["full.map.json"].resolve()
    assert paths["legacy.gff3"].resolve() not in subset.values()
    with pytest.raises(ValueError, match="no exact 'e2e' selection"):
        release._evaluation_input_paths(
            bench, tmp_path / "benchmarks.json", "e2e",
        )
    with pytest.raises(ValueError, match="panel is required"):
        release._evaluation_input_paths(
            bench, tmp_path / "benchmarks.json",
        )


def test_neutral_pair_scoring_publishes_target_truth_evidence(
        tmp_path, monkeypatch):
    inputs = _inputs(tmp_path / "inputs")
    outputs = {}
    documents = {}
    for label in ("candidate", "reference"):
        output = tmp_path / f"{label}.gff3"
        output.write_text(_annotation())
        outputs[label] = output
        documents[label] = {
            "profile": {
                "wall_clock_seconds": 1.0,
                "peak_rss_mb": 2.0,
            },
        }

    monkeypatch.setattr(
        release.evaluator,
        "build_reference",
        lambda *_args, **_kwargs: ({}, {}),
    )

    def evaluate(tool, *_args, out_dir=None, **_kwargs):
        # out_dir is positional in the production evaluator interface.
        actual_out_dir = Path(_args[4]) if out_dir is None else Path(out_dir)
        transcript = actual_out_dir / f"{tool}.transcripts.tsv"
        transcript.parent.mkdir(parents=True, exist_ok=True)
        transcript.write_text("ref_mrna_id\trecovered\n")
        return {"transcripts_tsv": str(transcript)}

    monkeypatch.setattr(release.evaluator, "evaluate_tool", evaluate)

    release._score_pair(
        inputs,
        outputs,
        documents,
        tmp_path / "cell",
        threads=2,
    )

    for label in ("candidate", "reference"):
        metrics = documents[label]["summary"]["target_truth"]
        assert metrics["gene"]["locus"]["f1"] == 1.0
        assert metrics["transcript"]["strand"]["f1"] == 1.0
        assert metrics["structure"]["exon"]["f1"] == 1.0
        assert metrics["parameters"]["mapping_required"] is True
        assert metrics["parameters"]["mapping_requirement_satisfied"] is True
        assert metrics["parameters"]["mapping_source_scope_validated"] is True
        assert documents[label]["target_truth_evidence"]["id_policy"] == (
            "ortholog-map"
        )
        artifact = documents[label]["evaluation_artifacts"]["target_truth"]
        assert Path(artifact["path"]).is_file()
        assert artifact["sha256"] == release.sha256_file(
            Path(artifact["path"])
        )
