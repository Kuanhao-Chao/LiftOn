from __future__ import annotations

import json

import pytest
from hypothesis import HealthCheck, given, settings, strategies as st

from benchmarks.compare.target_truth import (
    filter_truth_to_target_fasta,
    main,
    parse_annotation,
    score_target_truth,
    write_target_truth_metrics,
)


def _annotation(
    gene_id,
    transcript_id,
    *,
    seqid="chr1",
    start=100,
    strand="+",
    parent_id=None,
):
    end = start + 100
    parent_id = parent_id or gene_id
    return (
        f"{seqid}\tt\tgene\t{start}\t{end}\t.\t{strand}\t.\tID={gene_id}\n"
        f"{seqid}\tt\tmRNA\t{start}\t{end}\t.\t{strand}\t."
        f"\tID={transcript_id};Parent={parent_id}\n"
        f"{seqid}\tt\texon\t{start}\t{start + 39}\t.\t{strand}\t."
        f"\tID={transcript_id}.e1;Parent={transcript_id}\n"
        f"{seqid}\tt\texon\t{start + 60}\t{end}\t.\t{strand}\t."
        f"\tID={transcript_id}.e2;Parent={transcript_id}\n"
        f"{seqid}\tt\tCDS\t{start + 3}\t{start + 39}\t.\t{strand}\t0"
        f"\tID={transcript_id}.cds;Parent={transcript_id}\n"
        f"{seqid}\tt\tCDS\t{start + 60}\t{end - 4}\t.\t{strand}\t2"
        f"\tID={transcript_id}.cds;Parent={transcript_id}\n"
    )


def _write(path, body):
    path.write_text("##gff-version 3\n" + body)
    return path


def test_ortholog_truth_perfect_gene_and_transcript_mapping(tmp_path):
    prediction = _write(
        tmp_path / "prediction.gff3",
        _annotation("source_gene", "source_tx"),
    )
    truth = _write(
        tmp_path / "truth.gff3",
        _annotation("target_gene", "target_tx"),
    )
    mapping = tmp_path / "mapping.json"
    mapping.write_text(json.dumps({"mappings": [
        {
            "source_id": "source_gene",
            "truth_ids": ["target_gene"],
            "feature_type": "gene",
        },
        {
            "source_id": "source_tx",
            "truth_ids": ["target_tx"],
            "feature_type": "transcript",
        },
    ]}))

    result = score_target_truth(
        prediction, truth, ortholog_map=mapping,
    )

    for level in ("gene", "transcript"):
        for metric in ("locus", "strand", "copy"):
            assert result[level][metric]["f1"] == 1.0
    for metric in ("intron_chain", "intron", "exon", "CDS"):
        assert result["structure"][metric]["f1"] == 1.0


def test_source_scope_validation_rejects_unknown_status_and_target_ids(
        tmp_path):
    source = _write(
        tmp_path / "source.gff3",
        _annotation("source_gene", "source_tx"),
    )
    prediction = _write(
        tmp_path / "prediction.gff3",
        _annotation("source_gene", "source_tx"),
    )
    truth = _write(
        tmp_path / "truth.gff3",
        _annotation("target_gene", "target_tx"),
    )

    unknown_source = tmp_path / "unknown-source.json"
    unknown_source.write_text(json.dumps({"mappings": [{
        "source_id": "fabricated",
        "truth_ids": [],
        "feature_type": "gene",
        "status": "unscored",
    }]}))
    with pytest.raises(ValueError, match="unknown source gene ID"):
        score_target_truth(
            prediction,
            truth,
            ortholog_map=unknown_source,
            source_gff=source,
        )

    unknown_target = tmp_path / "unknown-target.json"
    unknown_target.write_text(json.dumps({"mappings": [{
        "source_id": "source_gene",
        "truth_ids": ["fabricated"],
        "feature_type": "gene",
        "status": "retained",
    }]}))
    with pytest.raises(ValueError, match="missing target truth ID"):
        score_target_truth(
            prediction,
            truth,
            ortholog_map=unknown_target,
            source_gff=source,
        )

    invalid_status = tmp_path / "invalid-status.json"
    invalid_status.write_text(json.dumps({"mappings": [{
        "source_id": "source_gene",
        "truth_ids": ["target_gene"],
        "feature_type": "gene",
        "status": "guessed",
    }]}))
    with pytest.raises(ValueError, match="unsupported mapping status"):
        score_target_truth(
            prediction,
            truth,
            ortholog_map=invalid_status,
            source_gff=source,
        )


def test_target_seqid_filter_is_exact_streaming_and_deterministic(tmp_path):
    truth = _write(
        tmp_path / "whole-truth.gff3",
        _annotation("g1", "t1", seqid="chr1")
        + _annotation("g2", "t2", seqid="chr2", start=300),
    )
    target = tmp_path / "subset.fa"
    target.write_text(">chr2 selected\n" + "A" * 500 + "\n")

    reports = []
    for name in ("first.gff3", "second.gff3"):
        reports.append(filter_truth_to_target_fasta(
            truth,
            target,
            tmp_path / name,
        ))

    first = tmp_path / "first.gff3"
    second = tmp_path / "second.gff3"
    assert first.read_bytes() == second.read_bytes()
    assert "chr1\t" not in first.read_text()
    assert "##sequence-region chr2 1 500\n" in first.read_text()
    filtered = parse_annotation(first)
    assert set(filtered.genes) == {"g2"}
    assert set(filtered.transcripts) == {"t2"}
    assert reports[0]["output"]["feature_seqids"] == ["chr2"]
    assert reports[0]["output"]["feature_records"] == 6
    assert reports[0]["output"]["dropped_feature_records"] == 6


def test_gene_ortholog_mapping_infers_transcript_scope(tmp_path):
    prediction = _write(
        tmp_path / "prediction.gff3",
        _annotation("source_gene", "source_tx"),
    )
    truth = _write(
        tmp_path / "truth.gff3",
        _annotation("target_gene", "different_target_tx"),
    )
    mapping = tmp_path / "mapping.tsv"
    mapping.write_text(
        "source_id\ttruth_id\tfeature_type\tstatus\n"
        "source_gene\ttarget_gene\tgene\tretained\n"
    )

    result = score_target_truth(
        prediction, truth, ortholog_map=mapping,
    )

    assert result["scope"]["transcript_groups"] == 1
    assert result["transcript"]["locus"]["f1"] == 1.0
    assert result["structure"]["exon"]["f1"] == 1.0


def test_truth_metrics_separate_locus_strand_copy_and_structure(tmp_path):
    prediction = _write(
        tmp_path / "prediction.gff3",
        _annotation("g", "t", strand="-").replace(
            "chr1\tt\texon\t160\t200", "chr1\tt\texon\t170\t200"
        ),
    )
    truth = _write(tmp_path / "truth.gff3", _annotation("g", "t"))

    result = score_target_truth(
        prediction, truth, id_policy="exact-id",
    )

    assert result["transcript"]["locus"]["f1"] == 1.0
    assert result["transcript"]["strand"]["f1"] == 0.0
    assert result["transcript"]["copy"]["f1"] == 1.0
    assert result["structure"]["intron_chain"]["f1"] == 0.0
    assert result["structure"]["exon"]["true_positive"] == 0


def test_independent_truth_requires_nonempty_ortholog_map(tmp_path):
    prediction = _write(
        tmp_path / "prediction.gff3",
        _annotation("g", "t"),
    )
    truth = _write(tmp_path / "truth.gff3", _annotation("g", "t"))
    empty_mapping = tmp_path / "empty.json"
    empty_mapping.write_text('{"mappings": []}\n')

    with pytest.raises(ValueError, match="explicit non-empty ortholog map"):
        score_target_truth(prediction, truth)
    with pytest.raises(ValueError, match="contains no mapping entries"):
        score_target_truth(
            prediction, truth, ortholog_map=empty_mapping,
        )

    result = score_target_truth(
        prediction, truth, id_policy="exact-id",
    )
    assert result["parameters"]["mapping_required"] is False
    assert result["parameters"]["mapping_requirement_satisfied"] is True


def test_gff3_ids_are_percent_decoded_after_raw_comma_splitting(tmp_path):
    prediction = _write(
        tmp_path / "prediction.gff3",
        "chr1\tt\tgene\t1\t10\t.\t+\t.\tID=source%3Ag%2Cone\n",
    )
    truth = _write(
        tmp_path / "truth.gff3",
        "chr1\tt\tgene\t1\t10\t.\t+\t.\tID=target%3Ag%2Cone\n",
    )
    mapping = tmp_path / "mapping.json"
    mapping.write_text(json.dumps({"mappings": [{
        "source_id": "source:g,one",
        "truth_ids": ["target:g,one"],
        "feature_type": "gene",
    }]}))

    result = score_target_truth(
        prediction, truth, ortholog_map=mapping,
    )

    assert result["gene"]["locus"]["f1"] == 1.0


def test_deleted_and_breakpoint_models_have_explicit_scope(tmp_path):
    prediction = _write(
        tmp_path / "prediction.gff3",
        _annotation("deleted_gene", "deleted_tx", start=10)
        + _annotation("excluded_gene", "excluded_tx", start=300),
    )
    truth = _write(
        tmp_path / "truth.gff3",
        _annotation("other_gene", "other_tx", start=600),
    )
    mapping = tmp_path / "mapping.json"
    mapping.write_text(json.dumps({"mappings": [
        {
            "source_id": "deleted_gene",
            "truth_ids": [],
            "feature_type": "gene",
            "status": "deleted",
        },
        {
            "source_id": "excluded_gene",
            "truth_ids": [],
            "feature_type": "gene",
            "status": "breakpoint_crossing",
        },
    ]}))

    result = score_target_truth(
        prediction, truth, ortholog_map=mapping,
    )

    assert result["scope"]["gene_groups"] == 1
    assert result["gene"]["copy"]["predicted"] == 1
    assert result["gene"]["copy"]["expected"] == 0
    assert result["gene"]["copy"]["false_positive"] == 1
    assert result["scope"]["gene"]["prediction_models_ignored"] == 1
    assert result["scope"]["gene"]["truth_models_ignored"] == 1


def test_ortholog_map_never_falls_back_to_identical_ids(tmp_path):
    prediction = _write(
        tmp_path / "prediction.gff3",
        _annotation("same_gene", "same_tx"),
    )
    truth = _write(
        tmp_path / "truth.gff3",
        _annotation("same_gene", "same_tx"),
    )
    mapping = tmp_path / "mapping.json"
    mapping.write_text(json.dumps({"mappings": [{
        "source_id": "same_tx",
        "truth_ids": [],
        "feature_type": "transcript",
        "status": "unscored",
    }]}))

    result = score_target_truth(
        prediction, truth, ortholog_map=mapping,
    )

    assert result["scope"]["gene_groups"] == 0
    assert result["scope"]["transcript_groups"] == 0
    assert result["scope"]["gene"]["prediction_models_ignored"] == 1
    assert result["scope"]["transcript"]["prediction_models_ignored"] == 1
    assert result["gene"]["locus"]["f1"] is None
    assert result["transcript"]["locus"]["f1"] is None


def test_copy_suffixes_match_one_source_ortholog_group(tmp_path):
    prediction = _write(
        tmp_path / "prediction.gff3",
        _annotation("g", "t", start=100)
        + _annotation("g_1", "t_1", start=500),
    )
    truth = _write(
        tmp_path / "truth.gff3",
        _annotation("tg1", "tt1", start=100)
        + _annotation("tg2", "tt2", start=500),
    )
    mapping = tmp_path / "mapping.json"
    mapping.write_text(json.dumps({"mappings": [
        {
            "source_id": "g",
            "truth_ids": ["tg1", "tg2"],
            "feature_type": "gene",
        },
        {
            "source_id": "t",
            "truth_ids": ["tt1", "tt2"],
            "feature_type": "transcript",
        },
    ]}))

    result = score_target_truth(
        prediction, truth, ortholog_map=mapping,
    )

    assert result["gene"]["copy"]["f1"] == 1.0
    assert result["gene"]["copy_count_exact"]["rate"] == 1.0
    assert result["transcript"]["locus"]["f1"] == 1.0


def test_copy_matching_maximizes_one_to_one_locus_recovery(tmp_path):
    prediction = _write(
        tmp_path / "prediction.gff3",
        "chr1\tt\tgene\t25\t125\t.\t+\t.\tID=g\n"
        "chr1\tt\tgene\t1\t80\t.\t+\t.\tID=g_1\n",
    )
    truth = _write(
        tmp_path / "truth.gff3",
        "chr1\tt\tgene\t1\t100\t.\t+\t.\tID=t1\n"
        "chr1\tt\tgene\t51\t150\t.\t+\t.\tID=t2\n",
    )
    mapping = tmp_path / "mapping.json"
    mapping.write_text(json.dumps({"mappings": [{
        "source_id": "g",
        "truth_ids": ["t1", "t2"],
        "feature_type": "gene",
    }]}))

    result = score_target_truth(
        prediction, truth, ortholog_map=mapping,
    )

    assert result["gene"]["locus"]["true_positive"] == 2
    assert result["gene"]["locus"]["f1"] == 1.0


def test_metric_writer_is_byte_deterministic(tmp_path):
    document = {"b": 2, "a": {"value": 1.0}}
    first = write_target_truth_metrics(tmp_path / "first.json", document)
    second = write_target_truth_metrics(tmp_path / "second.json", document)

    assert first.read_bytes() == second.read_bytes()


def test_target_truth_cli_writes_machine_readable_result(tmp_path):
    prediction = _write(
        tmp_path / "prediction.gff3",
        _annotation("g", "t"),
    )
    truth = _write(tmp_path / "truth.gff3", _annotation("g", "t"))
    output = tmp_path / "metrics.json"

    assert main([
        str(prediction),
        str(truth),
        "--exact-id",
        "--output",
        str(output),
    ]) == 0
    assert json.loads(output.read_text())["gene"]["locus"]["f1"] == 1.0


@given(
    offset=st.integers(min_value=-200, max_value=200),
    width=st.integers(min_value=10, max_value=100),
)
@settings(suppress_health_check=[HealthCheck.function_scoped_fixture])
def test_locus_property_is_bounded_and_symmetric(tmp_path, offset, width):
    truth = _write(
        tmp_path / "truth.gff3",
        "chr1\tt\tgene\t500\t600\t.\t+\t.\tID=g\n",
    )
    start = max(1, 500 + offset)
    prediction = _write(
        tmp_path / "prediction.gff3",
        f"chr1\tt\tgene\t{start}\t{start + width}\t.\t+\t.\tID=g\n",
    )

    metric = score_target_truth(
        prediction, truth, id_policy="exact-id",
    )["gene"]["locus"]

    for key in ("precision", "recall", "f1"):
        assert metric[key] is None or 0.0 <= metric[key] <= 1.0
