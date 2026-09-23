"""A transcript whose lifted CDS encodes no protein must not cost its gene.

Liftoff can lift a CDS down to a fragment shorter than one codon after its
phase -- on dog -> cat a 3-bp CDS at phase 2 -- which translates to the empty
string. ``lifton_parasail_align`` returned ``None`` for a ``None`` protein, the
contract every caller handles, but passed ``""`` on to
``parasail_align_protein_base``, whose guard raises ``LiftOnAlignmentError``
"so the caller can attribute the cause". No caller caught it: the error
reached the per-locus handler and the WHOLE gene -- every transcript, the good
ones included -- was left out of the annotation. The run said only
``partial_success``. 396 genes on the dog -> cat whole-genome lift, in every
release since the guard arrived in v1.0.9.
"""

from __future__ import annotations

import json
import textwrap

import pytest

from lifton import align


def _wrap(seq: str) -> str:
    return "\n".join(textwrap.wrap(seq, 60)) + "\n"


ORF = "ATG" + "GCT" * 31 + "TAA"          # 99 nt: M + A*31 + stop, at 101-199


# ---------------------------------------------------------------------------
# The contract, in isolation
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("lifted,reference", [("", "MA"), ("MA", "")])
def test_an_empty_protein_is_not_aligned(monkeypatch, lifted, reference):
    """``None``, as for a protein that cannot be built -- what callers handle
    and what ``Lifton_TRANS.align_coding_seq`` already returns for the same
    sequences -- rather than an exception no caller catches."""
    monkeypatch.setattr(align, "LiftOn_translate",
                        lambda *args: (reference, lifted, [3], []))
    result = align.lifton_parasail_align(
        object(), None, None, {"rna-ref": reference}, "rna-ref")
    assert result is None


def test_the_aligner_guard_itself_still_refuses_empty_input():
    """The guard stays: an empty sequence reaching the aligner is still a bug."""
    from lifton.exceptions import LiftOnAlignmentError
    with pytest.raises(LiftOnAlignmentError):
        align.parasail_align_protein_base("", "MA")


# ---------------------------------------------------------------------------
# End to end, through the real pipeline
# ---------------------------------------------------------------------------

@pytest.fixture
def hermetic_pipeline(monkeypatch):
    from lifton import lifton_utils as _lu, run_liftoff, run_miniprot

    monkeypatch.setattr(_lu, "check_miniprot_installed", lambda: None)
    monkeypatch.setattr(run_miniprot, "check_miniprot_installed", lambda: True)

    def _fail(*args, **kwargs):
        raise AssertionError("External liftoff/miniprot must NOT be invoked in tests")

    monkeypatch.setattr(run_liftoff, "run_liftoff", _fail)
    monkeypatch.setattr(run_miniprot, "run_miniprot", _fail)


@pytest.fixture
def truncated_workspace(tmp_path):
    """gene1 has two transcripts. Liftoff lifts tx1 whole; tx2's CDS comes out
    as 3 bp at phase 2, which leaves no codon to translate."""
    work = tmp_path / "work"
    work.mkdir()
    chrom = ["A"] * 600
    for i, ch in enumerate(ORF):
        chrom[100 + i] = ch
    seq = "".join(chrom)
    (work / "ref.fa").write_text(">chr1\n" + _wrap(seq))
    (work / "tgt.fa").write_text(">chr1\n" + _wrap(seq))
    (work / "ref.gff3").write_text(
        "##gff-version 3\n"
        "chr1\ttest\tgene\t101\t199\t.\t+\t.\tID=gene1;gene_biotype=protein_coding\n"
        "chr1\ttest\tmRNA\t101\t199\t.\t+\t.\tID=tx1;Parent=gene1\n"
        "chr1\ttest\texon\t101\t199\t.\t+\t.\tID=exon1;Parent=tx1\n"
        "chr1\ttest\tCDS\t101\t199\t.\t+\t0\tID=cds1;Parent=tx1\n"
        "chr1\ttest\tmRNA\t101\t199\t.\t+\t.\tID=tx2;Parent=gene1\n"
        "chr1\ttest\texon\t101\t199\t.\t+\t.\tID=exon2;Parent=tx2\n"
        "chr1\ttest\tCDS\t101\t199\t.\t+\t0\tID=cds2;Parent=tx2\n"
    )
    (work / "liftoff.gff3").write_text(
        "##gff-version 3\n"
        "chr1\tLiftoff\tgene\t101\t199\t.\t+\t.\tID=gene1;gene_biotype=protein_coding\n"
        "chr1\tLiftoff\tmRNA\t101\t199\t.\t+\t.\tID=tx1;Parent=gene1\n"
        "chr1\tLiftoff\texon\t101\t199\t.\t+\t.\tID=exon1;Parent=tx1\n"
        "chr1\tLiftoff\tCDS\t101\t199\t.\t+\t0\tID=cds1;Parent=tx1\n"
        "chr1\tLiftoff\tmRNA\t101\t199\t.\t+\t.\tID=tx2;Parent=gene1\n"
        "chr1\tLiftoff\texon\t101\t199\t.\t+\t.\tID=exon2;Parent=tx2\n"
        "chr1\tLiftoff\tCDS\t101\t103\t.\t+\t2\tID=cds2;Parent=tx2\n"
    )
    # One ordinary miniprot hit for tx1, which Step 8 suppresses (it overlaps
    # a lifted gene); an empty miniprot file is rejected at intake.
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
def test_the_gene_and_both_transcripts_are_emitted(truncated_workspace,
                                                   hermetic_pipeline, threads):
    text, manifest = _run(truncated_workspace, "out.gff3", *threads)
    rows = [line.split("\t") for line in text.splitlines()
            if line and not line.startswith("#")]
    genes = [r for r in rows if r[2] == "gene"]
    mrnas = {r[8].split("ID=")[1].split(";")[0] for r in rows if r[2] == "mRNA"}
    assert len(genes) == 1, "the gene was dropped for one unscorable transcript"
    assert mrnas == {"tx1", "tx2"}
    assert manifest["failures"] == []
    assert manifest["run"]["status"] == "success"
