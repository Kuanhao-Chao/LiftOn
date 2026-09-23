"""A gene-like type that also occurs as a gene's CHILD is lifted once, quietly.

The gene-like lift auto-detects every type with a top-level instance -- on dog
RefSeq that includes `tRNA` and `rRNA` -- and Step 7 then visited EVERY
instance of those types as a locus of its own, including the ones that are
children of a lifted gene and had already been emitted under it. Each such
child failed to resolve a reference gene (`ref_db[None]`, KeyError) and was
recorded as a pipeline failure: 395 of the 396 failures on the dog -> cat
whole-genome lift, all of them rows already in the output. The run reported
`partial_success`, and the one real loss among them was invisible.

This is the Iteration-20 shape -- detect from top-level instances, apply to
all of them -- at the one enumeration site that fix did not reach. The rule is
the same: skip an instance whose parent is itself of a lifted type.
"""

from __future__ import annotations

import json
import textwrap

import pytest


def _wrap(seq: str) -> str:
    return "\n".join(textwrap.wrap(seq, 60)) + "\n"


ORF = "ATG" + "GCT" * 31 + "TAA"


@pytest.fixture
def hermetic_pipeline(monkeypatch):
    from lifton import lifton_utils as _lu, run_liftoff, run_miniprot

    monkeypatch.setattr(_lu, "check_miniprot_installed", lambda: None)
    monkeypatch.setattr(run_miniprot, "check_miniprot_installed", lambda: True)

    def _fail(*args, **kwargs):
        raise AssertionError("External liftoff/miniprot must NOT be invoked in tests")

    monkeypatch.setattr(run_liftoff, "run_liftoff", _fail)
    monkeypatch.setattr(run_miniprot, "run_miniprot", _fail)


ROWS = (
    "chr1\t{src}\tgene\t101\t199\t.\t+\t.\tID=gene1;gene_biotype=protein_coding\n"
    "chr1\t{src}\tmRNA\t101\t199\t.\t+\t.\tID=tx1;Parent=gene1\n"
    "chr1\t{src}\texon\t101\t199\t.\t+\t.\tID=exon1;Parent=tx1\n"
    "chr1\t{src}\tCDS\t101\t199\t.\t+\t0\tID=cds1;Parent=tx1\n"
    # a tRNA gene: the tRNA is the gene's child
    "chr1\t{src}\tgene\t301\t372\t.\t+\t.\tID=gene-trnA;gene_biotype=tRNA\n"
    "chr1\t{src}\ttRNA\t301\t372\t.\t+\t.\tID=rna-trnA;Parent=gene-trnA\n"
    "chr1\t{src}\texon\t301\t372\t.\t+\t.\tID=exon-trnA;Parent=rna-trnA\n"
    # a top-level tRNA, which is what makes `tRNA` a gene-like type
    "chr1\t{src}\ttRNA\t401\t472\t.\t+\t.\tID=rna-trnB\n"
    "chr1\t{src}\texon\t401\t472\t.\t+\t.\tID=exon-trnB;Parent=rna-trnB\n"
)


@pytest.fixture
def workspace(tmp_path):
    work = tmp_path / "work"
    work.mkdir()
    chrom = ["A"] * 600
    for i, ch in enumerate(ORF):
        chrom[100 + i] = ch
    seq = "".join(chrom)
    (work / "ref.fa").write_text(">chr1\n" + _wrap(seq))
    (work / "tgt.fa").write_text(">chr1\n" + _wrap(seq))
    (work / "ref.gff3").write_text("##gff-version 3\n" + ROWS.format(src="test"))
    (work / "liftoff.gff3").write_text("##gff-version 3\n" + ROWS.format(src="Liftoff"))
    (work / "miniprot.gff3").write_text(
        "##gff-version 3\n"
        "chr1\tminiprot\tmRNA\t101\t199\t.\t+\t.\tID=MP1;Target=tx1 1 33\n"
        "chr1\tminiprot\tCDS\t101\t199\t.\t+\t0\tID=MP1.cds1;Parent=MP1\n"
    )
    (work / "out").mkdir()
    return work


def _run(work, name, *extra):
    from lifton import lifton as lifton_main
    out_gff = work / "out" / name
    argv = [str(work / "tgt.fa"), str(work / "ref.fa"),
            "-g", str(work / "ref.gff3"), "-L", str(work / "liftoff.gff3"),
            "-M", str(work / "miniprot.gff3"), "-o", str(out_gff),
            "-ad", "RefSeq", "--force", "-dir", str(work / f"dir_{name}"), *extra]
    lifton_main.run_all_lifton_steps(lifton_main.parse_args(argv))
    manifest = json.loads((work / f"dir_{name}" / "run_manifest.json").read_text())
    return out_gff.read_text(), manifest


@pytest.mark.parametrize("threads", [(), ("-t", "2", "--locus-pipeline")])
def test_a_child_trna_is_lifted_once_and_is_not_a_failure(workspace, hermetic_pipeline,
                                                          threads):
    text, manifest = _run(workspace, "out.gff3", *threads)
    ids = [line.split("ID=")[1].split(";")[0] for line in text.splitlines()
           if line and not line.startswith("#") and "\ttRNA\t" in line]
    assert sorted(ids) == ["rna-trnA", "rna-trnB"], ids
    assert manifest["failures"] == [], manifest["failures"]
    assert manifest["run"]["status"] == "success"
