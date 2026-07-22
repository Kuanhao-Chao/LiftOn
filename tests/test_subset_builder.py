from __future__ import annotations

import json
import subprocess
from pathlib import Path

from benchmarks.compare import subset_builder


def _write(path: Path, text: str) -> Path:
    path.write_text(text, encoding="utf-8")
    return path


def _fasta(seqid: str) -> str:
    return f">{seqid}\n" + "A" * 90 + "\n"


def test_subset_gff_preserves_gff3_directives_and_attributes(tmp_path):
    source = _write(
        tmp_path / "source.gff3",
        "##gff-version 3\n"
        "##sequence-region chr1 1 90\n"
        "##sequence-region chr2 1 90\n"
        "chr1\tt\tmRNA\t1\t90\t.\t+\t.\tID=t1;Parent=g1\n"
        "chr1\tt\tCDS\t1\t90\t.\t+\t0\tParent=t1;protein_id=p1\n"
        "chr2\tt\tmRNA\t1\t90\t.\t+\t.\tID=t2;Parent=g2\n",
    )
    output = tmp_path / "subset.gff3"

    info = subset_builder.subset_gff(str(source), "chr1", output)

    text = output.read_text(encoding="utf-8")
    assert text.count("##gff-version 3") == 1
    assert "##sequence-region chr1" in text
    assert "##sequence-region chr2" not in text
    assert "chr2\t" not in text
    assert info["annotation_format"] == "GFF3"
    assert info["protein_acc_to_mrna"] == {"p1": "t1"}


def test_subset_gff_preserves_raw_gtf_dialect(tmp_path):
    source = _write(
        tmp_path / "source.gtf",
        "#!genome-build GRCh38\n"
        'chr1\tt\ttranscript\t1\t90\t.\t+\t.\tgene_id "g1"; '
        'transcript_id "t1";\n'
        'chr1\tt\tCDS\t1\t90\t.\t+\t0\tgene_id "g1"; '
        'transcript_id "t1"; protein_id "p1";\n'
        'chr2\tt\ttranscript\t1\t90\t.\t+\t.\tgene_id "g2"; '
        'transcript_id "t2";\n',
    )
    output = tmp_path / "subset.gtf"

    info = subset_builder.subset_gff(
        str(source), "chr1", output, annotation_format="GTF",
    )

    text = output.read_text(encoding="utf-8")
    assert text.startswith("#!genome-build GRCh38\n")
    assert "##gff-version" not in text
    assert 'transcript_id "t1";' in text
    assert "chr2\t" not in text
    assert info["annotation_format"] == "GTF"
    assert info["protein_acc_to_mrna"] == {"p1": "t1"}


def test_build_subset_reuses_only_content_verified_manifest(
    tmp_path,
    monkeypatch,
):
    reference = _write(tmp_path / "reference.fa", _fasta("chr1"))
    target = _write(tmp_path / "target.fa", _fasta("chrT"))
    annotation = _write(
        tmp_path / "reference.gtf",
        'chr1\tt\tgene\t1\t90\t.\t+\t.\tgene_id "g1";\n'
        'chr1\tt\ttranscript\t1\t90\t.\t+\t.\tgene_id "g1"; '
        'transcript_id "t1";\n'
        'chr1\tt\tCDS\t1\t90\t.\t+\t0\tgene_id "g1"; '
        'transcript_id "t1"; protein_id "p1";\n',
    )
    proteins = _write(tmp_path / "reference.faa", ">p1\nMPEPTIDE\n")
    samtools = _write(tmp_path / "samtools", "#!/bin/sh\nexit 0\n")
    samtools.chmod(0o755)
    calls = []

    def fake_run(argv, **kwargs):
        calls.append([str(value) for value in argv])
        assert argv[1] == "faidx"
        fasta = Path(argv[2])
        header = next(
            line[1:].split()[0]
            for line in fasta.read_text(encoding="utf-8").splitlines()
            if line.startswith(">")
        )
        Path(str(fasta) + ".fai").write_text(
            f"{header}\t90\t6\t90\t91\n", encoding="ascii",
        )
        return subprocess.CompletedProcess(argv, 0, stdout="", stderr="")

    monkeypatch.setattr(subset_builder.subprocess, "run", fake_run)
    bench = {
        "id": "gtf_fixture",
        "species": "fixture",
        "cross_species": False,
        "ref_genome": str(reference),
        "ref_gff": str(annotation),
        "ref_proteins": str(proteins),
        "tgt_genome": str(target),
        "ref_chrom": "chr1",
        "tgt_chrom": "chrT",
        "ready_subset": True,
        "miniprot_target_space": "protein",
        "input_format": "GTF",
    }
    tools = {
        "samtools_bin": str(samtools),
        "minimap2_bin": "/bin/true",
    }
    work = tmp_path / "work"

    first = subset_builder.build_subset(bench, work, tools, threads=2)
    first_call_count = len(calls)
    assert first["schema_version"] == subset_builder.SUBSET_MANIFEST_SCHEMA_VERSION
    assert first["annotation_format"] == "GTF"
    assert Path(first["paths"]["ref_gff"]).suffix == ".gtf"
    assert first["request"]["inputs"]["ref_gff"]["sha256"]
    assert all(set(record) == {"path", "size", "sha256"}
               for record in first["outputs"].values())

    second = subset_builder.build_subset(bench, work, tools, threads=2)
    assert second == first
    assert len(calls) == first_call_count

    subset_manifest = work / "subset" / "subset.manifest.json"
    incomplete = json.loads(subset_manifest.read_text(encoding="utf-8"))
    del incomplete["outputs"]["ref_gff"]
    subset_manifest.write_text(json.dumps(incomplete), encoding="utf-8")
    rebuilt_incomplete = subset_builder.build_subset(
        bench, work, tools, threads=2,
    )
    assert len(calls) > first_call_count

    incomplete_call_count = len(calls)
    samtools.write_text("#!/bin/sh\nexit 0\n# changed\n", encoding="utf-8")
    rebuilt_tool = subset_builder.build_subset(bench, work, tools, threads=2)
    assert len(calls) > incomplete_call_count
    assert rebuilt_tool["request"]["tools"]["samtools"]["sha256"] != (
        rebuilt_incomplete["request"]["tools"]["samtools"]["sha256"]
    )

    tool_call_count = len(calls)
    subset_annotation = Path(rebuilt_tool["paths"]["ref_gff"])
    subset_annotation.write_text("tampered\n", encoding="utf-8")
    rebuilt = subset_builder.build_subset(bench, work, tools, threads=2)
    assert len(calls) > tool_call_count
    assert "tampered" not in subset_annotation.read_text(encoding="utf-8")
    assert rebuilt["outputs"]["ref_gff"]["sha256"] == (
        subset_builder._sha256_file(subset_annotation)
    )
