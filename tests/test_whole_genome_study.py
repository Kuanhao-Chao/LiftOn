from __future__ import annotations

import json
import zipfile
from pathlib import Path

import pytest

from benchmarks.compare import whole_genome_study


def _write_fasta(path: Path, records: dict[str, str]) -> None:
    path.write_text(
        "".join(f">{identifier}\n{sequence}\n" for identifier, sequence in records.items()),
        encoding="utf-8",
    )


def test_zip_validation_rejects_traversal_and_symlinks(tmp_path):
    traversal = tmp_path / "traversal.zip"
    with zipfile.ZipFile(traversal, "w") as handle:
        handle.writestr("../escape.gff", "unsafe")
    with pytest.raises(whole_genome_study.StudyError, match="unsafe ZIP member"):
        whole_genome_study.validate_zip_archive(traversal)

    symlink = tmp_path / "symlink.zip"
    member = zipfile.ZipInfo("annotation.gff")
    member.external_attr = 0o120777 << 16
    with zipfile.ZipFile(symlink, "w") as handle:
        handle.writestr(member, "target")
    with pytest.raises(whole_genome_study.StudyError, match="unsafe ZIP member"):
        whole_genome_study.validate_zip_archive(symlink)


def test_bee_normalization_rewrites_only_seqid_and_validates_lengths(tmp_path):
    target = tmp_path / "target.fa"
    raw = tmp_path / "raw.gff3"
    output = tmp_path / "normalized.gff3"
    _write_fasta(target, {"CM1.1": "A" * 10, "CM2.1": "C" * 7})
    raw.write_text(
        "##gff-version 3\n"
        "##sequence-region LG1 1 10\n"
        "##sequence-region LG2 1 7\n"
        "LG1\tMAKER\tgene\t1\t9\t.\t+\t.\tID=gene-LG1\n"
        "LG2\tMAKER\tmRNA\t2\t7\t.\t-\t.\tID=tx1;Parent=gene-LG1\n"
        "##FASTA\n>LG1\nAAAAAAAAAA\n>LG2\nCCCCCCC\n",
        encoding="utf-8",
    )

    audit = whole_genome_study.normalize_bee_annotation(
        raw,
        target,
        output,
        {"LG1": "CM1.1", "LG2": "CM2.1"},
        {"LG1": 10, "LG2": 7},
    )

    assert audit["feature_rows"] == 2
    assert audit["feature_counts"] == {"gene": 1, "mRNA": 1}
    rows = [line for line in output.read_text().splitlines() if not line.startswith("#")]
    assert rows == [
        "CM1.1\tMAKER\tgene\t1\t9\t.\t+\t.\tID=gene-LG1",
        "CM2.1\tMAKER\tmRNA\t2\t7\t.\t-\t.\tID=tx1;Parent=gene-LG1",
    ]

    with pytest.raises(whole_genome_study.StudyError, match="length mismatch"):
        whole_genome_study.normalize_bee_annotation(
            raw,
            target,
            tmp_path / "bad.gff3",
            {"LG1": "CM1.1", "LG2": "CM2.1"},
            {"LG1": 11, "LG2": 7},
        )


def test_fasta_lengths_rejects_duplicate_final_identifier(tmp_path):
    fasta = tmp_path / "duplicate.fa"
    fasta.write_text(">one\nA\n>one\nC\n", encoding="utf-8")

    with pytest.raises(whole_genome_study.StudyError, match="duplicate FASTA ID"):
        whole_genome_study.fasta_lengths(fasta)


def test_model_view_keeps_gene_models_and_removes_alignment_evidence(tmp_path):
    source = tmp_path / "released.gff3"
    output = tmp_path / "models.gff3"
    source.write_text(
        "##gff-version 3\n"
        "##sequence-region chr1 1 100\n"
        "chr1\tprovider\tgene\t1\t90\t.\t+\t.\tID=g1\n"
        "chr1\tprovider\tmRNA\t1\t90\t.\t+\t.\tID=t1;Parent=g1\n"
        "chr1\tprovider\texon\t1\t90\t.\t+\t.\tParent=t1\n"
        "chr1\tprovider\tCDS\t1\t90\t.\t+\t0\tParent=t1\n"
        "chr1\tprovider\tgene\t91\t94\t.\t+\t.\tID=split\n"
        "chr1\tprovider\tgene\t95\t98\t.\t-\t.\tID=split\n"
        "chr1\tprovider\tmRNA\t91\t98\t.\t?\t.\tID=bad;Parent=split\n"
        "chr1\tprovider\tCDS\t91\t98\t.\t?\t0\tParent=bad\n"
        "chr1\tprovider\tprotein_match\t1\t80\t.\t+\t.\tID=hit1\n",
        encoding="utf-8",
    )

    audit = whole_genome_study.build_model_view(source, output)

    text = output.read_text(encoding="utf-8")
    assert "protein_match" not in text
    assert "ID=split" not in text
    assert "ID=bad" not in text
    assert audit["feature_counts"] == {
        "CDS": 1, "exon": 1, "gene": 1, "mRNA": 1,
    }
    assert audit["excluded_feature_counts"] == {
        "CDS": 1, "gene": 2, "mRNA": 1,
    }
    assert audit["sequence_regions"] == 1


def test_annotation_lock_is_bound_to_every_content_addressed_object(tmp_path):
    study = whole_genome_study.load_study()
    root = tmp_path / "study-cache"
    artifact = root / "objects" / "sha256" / "aa" / "artifact"
    artifact.parent.mkdir(parents=True)
    artifact.write_text("annotation bytes\n", encoding="utf-8")
    record = whole_genome_study.file_record(artifact)
    record["path"] = str(artifact.relative_to(root))
    targets = {
        pair["id"]: {
            "identity": pair["target_annotation"],
            "archive": dict(record),
            "artifacts": {"truth_gff": dict(record)},
        }
        for pair in study["pairs"]
    }
    lock = {
        "schema_version": whole_genome_study.SCHEMA_VERSION,
        "kind": "lifton-v1.0.11-biology-study-annotation-lock",
        "created_at": "2026-08-13T00:00:00Z",
        "study": whole_genome_study.file_record(
            whole_genome_study.DEFAULT_STUDY
        ),
        "study_sha256": whole_genome_study.canonical_sha256(study),
        "benchmark_registry": whole_genome_study.file_record(
            whole_genome_study.DEFAULT_BENCHMARK_REGISTRY
        ),
        "targets": targets,
    }
    lock["fingerprint"] = whole_genome_study.canonical_sha256(lock)
    lock_path = root / whole_genome_study.LOCK_NAME
    lock_path.write_text(json.dumps(lock), encoding="utf-8")

    result = whole_genome_study.verify_annotation_lock(
        whole_genome_study.DEFAULT_STUDY, lock_path, root,
    )
    assert result["verified"] is True
    assert result["targets"] == 7

    artifact.write_text("corrupted bytes\n", encoding="utf-8")
    with pytest.raises(whole_genome_study.StudyError, match="object changed"):
        whole_genome_study.verify_annotation_lock(
            whole_genome_study.DEFAULT_STUDY, lock_path, root,
        )
