from __future__ import annotations

import json
import subprocess
from pathlib import Path

import pytest

from benchmarks.compare import ortholog_scope, target_truth


def _write(path: Path, text: str) -> Path:
    path.write_text(text, encoding="utf-8")
    return path


def _gff3(prefix: str, transcripts: tuple[str, ...]) -> str:
    rows = ["##gff-version 3"]
    for index, transcript_id in enumerate(transcripts, start=1):
        gene_id = f"{prefix}g{index}"
        start = index * 100
        rows.extend((
            (
                f"chr1\tfixture\tgene\t{start}\t{start + 80}\t.\t+\t.\t"
                f"ID={gene_id}"
            ),
            (
                f"chr1\tfixture\tmRNA\t{start}\t{start + 80}\t.\t+\t.\t"
                f"ID={transcript_id};Parent={gene_id}"
            ),
            (
                f"chr1\tfixture\tCDS\t{start}\t{start + 80}\t.\t+\t0\t"
                f"Parent={transcript_id}"
            ),
        ))
    return "\n".join(rows) + "\n"


def _gtf(prefix: str, transcripts: tuple[str, ...]) -> str:
    rows = []
    for index, transcript_id in enumerate(transcripts, start=1):
        gene_id = f"{prefix}g{index}"
        start = index * 100
        attributes = (
            f'gene_id "{gene_id}"; transcript_id "{transcript_id}";'
        )
        rows.extend((
            (
                f"chr1\tfixture\tgene\t{start}\t{start + 80}\t.\t+\t.\t"
                f'gene_id "{gene_id}";'
            ),
            (
                f"chr1\tfixture\ttranscript\t{start}\t{start + 80}\t.\t+\t."
                f"\t{attributes}"
            ),
            (
                f"chr1\tfixture\tCDS\t{start}\t{start + 80}\t.\t+\t0\t"
                f"{attributes}"
            ),
        ))
    return "\n".join(rows) + "\n"


def _inputs(tmp_path: Path):
    source_gff = _write(
        tmp_path / "source.gff3",
        _gff3("s", ("s1", "s2", "s3")),
    )
    target_gff = _write(
        tmp_path / "target.gff3",
        _gff3("t", ("t1", "t2")),
    )
    source_fa = _write(tmp_path / "source.fa", ">chr1\n" + "A" * 500 + "\n")
    target_fa = _write(tmp_path / "target.fa", ">chr1\n" + "A" * 500 + "\n")
    return source_gff, source_fa, target_gff, target_fa


def test_frozen_record_allows_filesystem_metadata_drift(tmp_path):
    path = _write(tmp_path / "input.fa", ">chr1\nACGT\n")
    expected = ortholog_scope._input_record(path)

    path.chmod(0o600)

    assert ortholog_scope._input_record(path)["sha256"] == expected["sha256"]
    ortholog_scope._assert_frozen_records({"source_genome": expected})


def test_frozen_record_rejects_byte_drift_with_diagnostics(tmp_path):
    path = _write(tmp_path / "input.fa", ">chr1\nACGT\n")
    expected = ortholog_scope._input_record(path)
    _write(path, ">chr1\nTGCA\n")

    with pytest.raises(
        ortholog_scope.ScopeBuildError,
        match=r"source_genome .*sha256: expected",
    ):
        ortholog_scope._assert_frozen_records({"source_genome": expected})


def test_hierarchy_excludes_ambiguous_multi_part_models(tmp_path):
    annotation = _write(
        tmp_path / "multipart.gff3",
        "##gff-version 3\n"
        "chr1\tRefSeq\tgene\t1\t50\t.\t+\t.\tID=g1;part=1\n"
        "chr1\tRefSeq\tgene\t60\t90\t.\t-\t.\tID=g1;part=2\n"
        "chr1\tRefSeq\tmRNA\t1\t50\t.\t+\t.\tID=t1;Parent=g1\n"
        "chr1\tRefSeq\tgene\t100\t190\t.\t+\t.\tID=g2\n"
        "chr1\tRefSeq\tmRNA\t100\t190\t.\t+\t.\tID=t2;Parent=g2\n",
    )

    hierarchy = ortholog_scope.parse_hierarchy(annotation)

    assert hierarchy.genes == ("g2",)
    assert hierarchy.transcripts == ("t2",)
    assert hierarchy.excluded_genes == ("g1",)
    assert hierarchy.excluded_transcripts == ("t1",)


def test_hierarchy_ignores_top_level_non_protein_rnas(tmp_path):
    annotation = _write(
        tmp_path / "organellar.gff3",
        "##gff-version 3\n"
        "chrM\tRefSeq\ttRNA\t1\t70\t.\t+\t.\tID=trna1\n"
        "chr1\tRefSeq\tgene\t100\t190\t.\t+\t.\tID=g1\n"
        "chr1\tRefSeq\tmRNA\t100\t190\t.\t+\t.\tID=t1;Parent=g1\n",
    )

    hierarchy = ortholog_scope.parse_hierarchy(annotation)

    assert hierarchy.genes == ("g1",)
    assert hierarchy.transcripts == ("t1",)


def test_protein_normalization_excludes_annotation_translation_failures(
        tmp_path):
    hierarchy = ortholog_scope.Hierarchy(
        genes=("g1", "g2"),
        transcript_to_gene={"t1": "g1", "t2": "g2"},
        format_counts={"gff3": 2},
    )
    raw = _write(
        tmp_path / "raw.fa",
        ">t1\nMA*AA\n>t2\nMBBBB*\n>pseudogene-direct-cds\nMCCCC*\n",
    )
    output = tmp_path / "normalized.fa"

    statistics = ortholog_scope.normalize_protein_fasta(
        raw, output, hierarchy,
    )

    assert output.read_text() == ">t2\nMBBBB\n"
    assert statistics["excluded_invalid_protein_ids"] == ["t1"]
    assert statistics["ignored_outside_hierarchy_protein_ids"] == [
        "pseudogene-direct-cds"
    ]
    assert statistics["missing_transcript_proteins"] == ["t1"]


def _tools(tmp_path: Path) -> tuple[Path, Path]:
    gffread = _write(tmp_path / "gffread", "fixture executable\n")
    mmseqs = _write(tmp_path / "mmseqs", "fixture executable\n")
    gffread.chmod(0o755)
    mmseqs.chmod(0o755)
    return gffread, mmseqs


class FakeRunner:
    def __init__(
        self,
        *,
        source_proteins: str,
        target_proteins: str,
        hits: str,
    ):
        self.source_proteins = source_proteins
        self.target_proteins = target_proteins
        self.hits = hits
        self.calls = []

    def __call__(self, argv, **kwargs):
        self.calls.append((list(argv), dict(kwargs)))
        executable = Path(argv[0]).name
        if executable == "gffread" and argv[1:] == ["--version"]:
            return subprocess.CompletedProcess(
                argv, 0, stdout="gffread 0.12.7\n", stderr="",
            )
        if executable == "mmseqs" and argv[1:] == ["version"]:
            return subprocess.CompletedProcess(
                argv, 0, stdout="15-6f452\n", stderr="",
            )
        if executable == "gffread":
            output = Path(argv[argv.index("-y") + 1])
            annotation = Path(argv[-1]).name
            output.write_text(
                self.source_proteins
                if annotation.startswith("source")
                else self.target_proteins,
                encoding="utf-8",
            )
            return subprocess.CompletedProcess(argv, 0, stdout="", stderr="")
        if executable == "mmseqs" and argv[1] == "easy-rbh":
            Path(argv[4]).write_text(self.hits, encoding="utf-8")
            return subprocess.CompletedProcess(argv, 0, stdout="", stderr="")
        raise AssertionError(f"unexpected command: {argv}")


def _runner() -> FakeRunner:
    return FakeRunner(
        source_proteins=(
            ">s2 gffread metadata\nMBBBB*\n"
            ">s1\nMAAAA*\n"
        ),
        target_proteins=(
            ">t2\nMCCCC*\n"
            ">t1\nMDDDD*\n"
        ),
        hits=(
            "s1\tt1\t0.90\t90\t0.80\t0.80\t1e-20\t100\n"
            "s2\tt1\t0.95\t90\t0.80\t0.80\t1e-30\t90\n"
            "s2\tt2\t0.20\t90\t0.80\t0.80\t1e-20\t200\n"
        ),
    )


def test_build_scope_is_deterministic_provenance_rich_and_explicit(
        tmp_path):
    source_gff, source_fa, target_gff, target_fa = _inputs(tmp_path)
    gffread, mmseqs = _tools(tmp_path)
    runner = _runner()

    for name in ("first", "second"):
        ortholog_scope.build_ortholog_scope(
            source_gff,
            source_fa,
            target_gff,
            target_fa,
            tmp_path / name,
            gffread=gffread,
            mmseqs=mmseqs,
            minimum_identity=0.30,
            minimum_reciprocal_coverage=0.50,
            maximum_evalue=1e-5,
            threads=2,
            minimum_gene_groups=1,
            minimum_transcript_groups=1,
            runner=runner,
        )

    first = tmp_path / "first"
    second = tmp_path / "second"
    assert (first / "ortholog_map.json").read_bytes() == (
        second / "ortholog_map.json"
    ).read_bytes()
    assert (first / "ortholog_scope.manifest.json").read_bytes() == (
        second / "ortholog_scope.manifest.json"
    ).read_bytes()
    assert (first / "source.proteins.fa").read_text() == (
        ">s1\nMAAAA\n>s2\nMBBBB\n"
    )

    mapping = json.loads((first / "ortholog_map.json").read_text())
    entries = {
        (entry["feature_type"], entry["source_id"]): entry
        for entry in mapping["mappings"]
    }
    assert entries[("gene", "sg1")]["truth_ids"] == ["tg1"]
    assert entries[("transcript", "s1")]["truth_ids"] == ["t1"]
    assert entries[("gene", "sg2")]["status"] == "unscored"
    assert entries[("transcript", "s2")]["reason"] == (
        "transcript_pair_conflict"
    )
    assert entries[("gene", "sg3")]["reason"] == "no_extracted_protein"
    assert entries[("transcript", "s3")]["reason"] == (
        "no_extracted_protein"
    )
    assert len(entries) == 6
    assert mapping["metadata"]["counts"]["gene_groups_retained"] == 1
    assert mapping["metadata"]["counts"]["rbh_hits_raw"] == 3
    assert target_truth.load_ortholog_map(
        first / "ortholog_map.json"
    )[1]["entries"] == 6
    scored = target_truth.score_target_truth(
        source_gff,
        target_gff,
        ortholog_map=first / "ortholog_map.json",
        source_gff=source_gff,
    )
    assert scored["parameters"]["mapping_source_scope_validated"] is True
    assert scored["inputs"]["mapping"]["semantic_validation"][
        "input_fingerprints_validated"
    ] is True
    with pytest.raises(ValueError, match="requires source_gff"):
        target_truth.score_target_truth(
            source_gff,
            target_gff,
            ortholog_map=first / "ortholog_map.json",
        )

    manifest = json.loads(
        (first / "ortholog_scope.manifest.json").read_text()
    )
    assert manifest["tools"]["gffread"]["version"] == "gffread 0.12.7"
    assert manifest["tools"]["mmseqs"]["version"] == "15-6f452"
    assert manifest["parameters"]["threads"] == 2
    assert all(
        len(record["sha256"]) == 64
        for record in manifest["inputs"].values()
    )
    assert all(command["shell"] is False for command in manifest["commands"])
    assert all(isinstance(argv, list) for argv, _kwargs in runner.calls)
    assert all("shell" not in kwargs for _argv, kwargs in runner.calls)


def test_semantic_validation_rejects_unknown_ids_counts_and_stale_inputs(
        tmp_path):
    source_gff, source_fa, target_gff, target_fa = _inputs(tmp_path)
    gffread, mmseqs = _tools(tmp_path)
    output = tmp_path / "scope"
    ortholog_scope.build_ortholog_scope(
        source_gff,
        source_fa,
        target_gff,
        target_fa,
        output,
        gffread=gffread,
        mmseqs=mmseqs,
        minimum_gene_groups=1,
        minimum_transcript_groups=1,
        runner=_runner(),
    )
    original = json.loads((output / "ortholog_map.json").read_text())

    unknown_source = json.loads(json.dumps(original))
    unknown_source["mappings"][0]["source_id"] = "fabricated-gene"
    with pytest.raises(
        ortholog_scope.ScopeBuildError, match="unknown source gene ID",
    ):
        ortholog_scope.validate_mapping_against_annotations(
            unknown_source,
            source_annotation=source_gff,
            target_annotation=target_gff,
        )

    unknown_target = json.loads(json.dumps(original))
    retained = next(
        entry for entry in unknown_target["mappings"]
        if entry["feature_type"] == "transcript"
        and entry["status"] == "retained"
    )
    retained["truth_ids"] = ["fabricated-transcript"]
    with pytest.raises(
        ortholog_scope.ScopeBuildError, match="unknown target transcript ID",
    ):
        ortholog_scope.validate_mapping_against_annotations(
            unknown_target,
            source_annotation=source_gff,
            target_annotation=target_gff,
        )

    stale_count = json.loads(json.dumps(original))
    stale_count["metadata"]["counts"]["gene_groups_retained"] += 1
    with pytest.raises(
        ortholog_scope.ScopeBuildError, match="count gene_groups_retained is stale",
    ):
        ortholog_scope.validate_mapping_against_annotations(
            stale_count,
            source_annotation=source_gff,
            target_annotation=target_gff,
        )

    source_gff.write_text(source_gff.read_text() + "# changed bytes\n")
    with pytest.raises(
        ortholog_scope.ScopeBuildError,
        match="source annotation fingerprint does not match",
    ):
        ortholog_scope.validate_mapping_against_annotations(
            original,
            source_annotation=source_gff,
            target_annotation=target_gff,
        )


def test_registry_finalization_is_deterministic_and_semantically_validated(
        tmp_path):
    inputs = tmp_path / "inputs"
    inputs.mkdir()
    source_gff, source_fa, target_gff, target_fa = _inputs(inputs)
    gffread, mmseqs = _tools(tmp_path)
    mapping_root = tmp_path / "ortholog-scopes"
    ortholog_scope.build_ortholog_scope(
        source_gff,
        source_fa,
        target_gff,
        target_fa,
        mapping_root / "biological-demo",
        gffread=gffread,
        mmseqs=mmseqs,
        minimum_gene_groups=1,
        minimum_transcript_groups=1,
        runner=_runner(),
    )
    manifest = _write(
        tmp_path / "manifest.json",
        json.dumps({
            "schema_version": 1,
            "scenarios": [{
                "id": "biological-demo",
                "kind": "biological",
            }],
        }),
    )
    first_path = tmp_path / "registry.first.json"
    second_path = tmp_path / "registry.second.json"

    first = ortholog_scope.finalize_mapping_registry(
        manifest, mapping_root, first_path,
    )
    assert ortholog_scope.main([
        "finalize-registry",
        str(manifest),
        str(mapping_root),
        str(second_path),
    ]) == 0
    second = json.loads(second_path.read_text())

    assert first == second
    assert first_path.read_bytes() == second_path.read_bytes()
    record = first["mappings"]["biological-demo"]
    assert record["id_policy"] == "ortholog-map"
    assert len(record["sha256"]) == 64
    assert (first_path.parent / record["path"]).resolve() == (
        mapping_root / "biological-demo" / "ortholog_map.json"
    ).resolve()
    bundle_path = mapping_root / "biological-demo"
    validated_bundle = ortholog_scope.validate_scope_bundle(bundle_path)
    assert validated_bundle["mapping_sha256"] == record["sha256"]
    unexpected = bundle_path / "partial.tmp"
    unexpected.write_text("unexpected")
    with pytest.raises(
        ortholog_scope.ScopeBuildError, match="unexpected artifact",
    ):
        ortholog_scope.validate_scope_bundle(bundle_path)
    unexpected.unlink()

    mapping_path = mapping_root / "biological-demo" / "ortholog_map.json"
    mapping_path.write_text(mapping_path.read_text() + "\n")
    with pytest.raises(
        ortholog_scope.ScopeBuildError, match="manifest is stale or tampered",
    ):
        ortholog_scope.finalize_mapping_registry(
            manifest,
            mapping_root,
            tmp_path / "registry.tampered.json",
        )


def test_prepare_publishes_normalized_proteins_without_running_mmseqs(
        tmp_path):
    source_gff, source_fa, target_gff, target_fa = _inputs(tmp_path)
    gffread, _mmseqs = _tools(tmp_path)
    runner = _runner()
    output = tmp_path / "prepared"

    manifest = ortholog_scope.prepare_proteins(
        source_gff,
        source_fa,
        target_gff,
        target_fa,
        output,
        gffread=gffread,
        runner=runner,
    )

    assert manifest["method"] == ortholog_scope.PREPARE_METHOD
    assert (output / "prepare.manifest.json").is_file()
    assert (output / "source.proteins.fa").read_text().startswith(">s1\n")
    assert all(
        not (Path(argv[0]).name == "mmseqs" and "easy-rbh" in argv)
        for argv, _kwargs in runner.calls
    )


def test_gtf_hierarchy_and_subset_filter_remove_out_of_scope_groups(
        tmp_path):
    source_full = _write(
        tmp_path / "source-full.gtf",
        _gtf("s", ("s1", "s2")),
    )
    target_full = _write(
        tmp_path / "target-full.gtf",
        _gtf("t", ("t1", "t2")),
    )
    source_subset = _write(
        tmp_path / "source-subset.gtf",
        _gtf("s", ("s1",)),
    )
    target_subset = _write(
        tmp_path / "target-subset.gtf",
        _gtf("t", ("t1",)),
    )
    source_hierarchy = ortholog_scope.parse_hierarchy(source_subset)
    target_hierarchy = ortholog_scope.parse_hierarchy(target_subset)
    assert source_hierarchy.transcript_to_gene == {"s1": "sg1"}
    assert target_hierarchy.transcript_to_gene == {"t1": "tg1"}

    full = {
        "schema_version": ortholog_scope.SCHEMA_VERSION,
        "method": ortholog_scope.METHOD,
        "metadata": {
            "scope": "full",
            "counts": {
                "source_genes": 2,
                "source_transcripts": 2,
                "target_genes": 2,
                "target_transcripts": 2,
                "gene_groups_retained": 2,
                "transcript_groups_retained": 2,
            },
            "provenance": {
                "inputs": {
                    "source_annotation": ortholog_scope._input_record(
                        source_full
                    ),
                    "target_annotation": ortholog_scope._input_record(
                        target_full
                    ),
                },
            },
        },
        "mappings": [
            {
                "source_id": "sg1",
                "truth_ids": ["tg1"],
                "feature_type": "gene",
                "status": "retained",
            },
            {
                "source_id": "sg2",
                "truth_ids": ["tg2"],
                "feature_type": "gene",
                "status": "retained",
            },
            {
                "source_id": "s1",
                "truth_ids": ["t1"],
                "feature_type": "transcript",
                "status": "retained",
            },
            {
                "source_id": "s2",
                "truth_ids": ["t2"],
                "feature_type": "transcript",
                "status": "retained",
            },
        ],
    }
    filtered = ortholog_scope.filter_scope_document(
        full,
        source_subset_annotation=source_subset,
        target_subset_annotation=target_subset,
    )

    assert [
        (entry["feature_type"], entry["source_id"], entry["truth_ids"])
        for entry in filtered["mappings"]
    ] == [
        ("gene", "sg1", ["tg1"]),
        ("transcript", "s1", ["t1"]),
    ]
    assert filtered["metadata"]["scope"] == "subset"
    assert filtered["metadata"]["counts"]["gene_groups_retained"] == 1

    other_target = _write(
        tmp_path / "other-target.gtf",
        _gtf("x", ("x1",)),
    )
    unscored = ortholog_scope.filter_scope_document(
        full,
        source_subset_annotation=source_subset,
        target_subset_annotation=other_target,
        minimum_gene_groups=0,
        minimum_transcript_groups=0,
    )
    assert all(
        entry["status"] == "unscored"
        and entry["reason"] == "target_outside_subset"
        for entry in unscored["mappings"]
    )


def test_gene_resolution_maximizes_one_to_one_group_count():
    source = ortholog_scope.Hierarchy(
        genes=("sg1", "sg2"),
        transcript_to_gene={
            "s1a": "sg1",
            "s1b": "sg1",
            "s2": "sg2",
        },
        format_counts={"gff3": 1},
    )
    target = ortholog_scope.Hierarchy(
        genes=("tg1", "tg2"),
        transcript_to_gene={
            "t1a": "tg1",
            "t1b": "tg1",
            "t2": "tg2",
        },
        format_counts={"gff3": 1},
    )

    def hit(source_id, target_id, bits):
        return ortholog_scope.Hit(
            source_id,
            target_id,
            0.9,
            100,
            0.9,
            0.9,
            1e-20,
            bits,
            1,
        )

    selected, _evidence = ortholog_scope._resolve_gene_pairs(
        [
            hit("s1a", "t1a", 100),
            hit("s1b", "t2", 90),
            hit("s2", "t1b", 80),
        ],
        source,
        target,
    )

    assert selected == [("sg1", "tg2"), ("sg2", "tg1")]


def test_minimum_group_failure_leaves_no_partial_output(tmp_path):
    source_gff, source_fa, target_gff, target_fa = _inputs(tmp_path)
    gffread, mmseqs = _tools(tmp_path)
    runner = FakeRunner(
        source_proteins=">s1\nMAAAA\n",
        target_proteins=">t1\nMBBBB\n",
        hits="s1\tt1\t0.1\t90\t0.9\t0.9\t1e-20\t100\n",
    )
    output = tmp_path / "failed"

    with pytest.raises(
        ortholog_scope.ScopeBuildError,
        match="retained gene groups 0",
    ):
        ortholog_scope.build_ortholog_scope(
            source_gff,
            source_fa,
            target_gff,
            target_fa,
            output,
            gffread=gffread,
            mmseqs=mmseqs,
            runner=runner,
        )

    assert not output.exists()
    assert not list(tmp_path.glob(".failed.tmp-*"))


def test_build_dry_run_hashes_inputs_but_executes_no_tools(
        tmp_path, capsys):
    source_gff, source_fa, target_gff, target_fa = _inputs(tmp_path)
    gffread, mmseqs = _tools(tmp_path)
    output = tmp_path / "dry-output"

    code = ortholog_scope.main([
        "build",
        str(source_gff),
        str(source_fa),
        str(target_gff),
        str(target_fa),
        str(output),
        "--gffread",
        str(gffread),
        "--mmseqs",
        str(mmseqs),
        "--dry-run",
    ])

    assert code == 0
    assert not output.exists()
    plan = json.loads(capsys.readouterr().out)
    assert plan["dry_run"] is True
    assert plan["action"] == "build"
    assert len(plan["commands"]) == 3
    assert all(command["shell"] is False for command in plan["commands"])
    assert all(
        len(record["sha256"]) == 64 for record in plan["inputs"].values()
    )


def test_malformed_mmseqs_output_fails_closed(tmp_path):
    source_gff, source_fa, target_gff, target_fa = _inputs(tmp_path)
    gffread, mmseqs = _tools(tmp_path)
    runner = FakeRunner(
        source_proteins=">s1\nMAAAA\n",
        target_proteins=">t1\nMBBBB\n",
        hits="s1\tt1\tnot-a-number\n",
    )

    with pytest.raises(ortholog_scope.ScopeBuildError, match="expected 8"):
        ortholog_scope.build_ortholog_scope(
            source_gff,
            source_fa,
            target_gff,
            target_fa,
            tmp_path / "malformed",
            gffread=gffread,
            mmseqs=mmseqs,
            runner=runner,
        )
