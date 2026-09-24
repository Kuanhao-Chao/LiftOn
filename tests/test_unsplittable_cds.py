"""A CDS that cannot be split at exon boundaries must not cost its gene.

v1.0.14 made ``Lifton_TRANS.add_cds`` refuse a CDS that overlaps several exons
it cannot be split across unambiguously (overlapping exons, spans out of order,
an exon that already carries a CDS) instead of attaching it whole. The refusal
is right; letting it escape was not. In Step 7 ``process_liftoff`` did not
catch it, so the error reached the per-locus handler and the whole gene -- its
good transcripts included -- was left out of the annotation. ``-E`` evaluation
did not catch it either, and its serial loop aborted the run.
"""

from __future__ import annotations

import json
import textwrap

import pytest

from tests.test_run_evaluation import _args, _capture_orf, _db, _evaluate, _ref_features, _row


def _wrap(seq: str) -> str:
    return "\n".join(textwrap.wrap(seq, 60)) + "\n"


ORF = "ATG" + "GCT" * 31 + "TAA"          # 99 nt: M + A*31 + stop, at 101-199


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
def unsplittable_workspace(tmp_path):
    """gene1 has two transcripts. Liftoff lifts tx1 whole; tx2 comes out with
    two overlapping exons and one CDS across both, which no split resolves."""
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
        "chr1\tLiftoff\texon\t101\t160\t.\t+\t.\tID=exon2a;Parent=tx2\n"
        "chr1\tLiftoff\texon\t150\t199\t.\t+\t.\tID=exon2b;Parent=tx2\n"
        "chr1\tLiftoff\tCDS\t101\t199\t.\t+\t0\tID=cds2;Parent=tx2\n"
    )
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
def test_only_the_unsplittable_transcript_is_left_out(unsplittable_workspace,
                                                     hermetic_pipeline, threads):
    text, manifest = _run(unsplittable_workspace, "out.gff3", *threads)
    rows = [line.split("\t") for line in text.splitlines()
            if line and not line.startswith("#")]
    genes = [r for r in rows if r[2] == "gene"]
    mrnas = {r[8].split("ID=")[1].split(";")[0] for r in rows if r[2] == "mRNA"}
    assert len(genes) == 1, "the gene was dropped for one unsplittable transcript"
    assert mrnas == {"tx1"}
    assert manifest["failures"] == []
    assert manifest["run"]["status"] == "success"
    assert manifest["counts"]["dropped_cds_spanning_exons"] == 1


# ---------------------------------------------------------------------------
# -E evaluation
# ---------------------------------------------------------------------------

def test_evaluation_skips_an_unsplittable_transcript_and_scores_the_rest(monkeypatch):
    ref_db = _db([
        _row("gene", 1, 200, "ID=gene1;gene_biotype=protein_coding"),
        _row("mRNA", 1, 90, "ID=bad;Parent=gene1"),
        _row("mRNA", 101, 200, "ID=good;Parent=gene1"),
    ])
    tgt_db = _db([
        _row("gene", 1, 200, "ID=gene1"),
        _row("mRNA", 1, 90, "ID=bad;Parent=gene1"),
        _row("exon", 1, 60, "ID=bad-ex1;Parent=bad"),
        _row("exon", 50, 90, "ID=bad-ex2;Parent=bad"),
        _row("CDS", 1, 90, "ID=bad-cds;Parent=bad", frame="0"),
        _row("mRNA", 101, 200, "ID=good;Parent=gene1"),
        _row("exon", 101, 200, "ID=good-ex;Parent=good"),
        _row("CDS", 101, 199, "ID=good-cds;Parent=good", frame="0"),
    ])
    calls = _capture_orf(monkeypatch)

    _, score = _evaluate(ref_db, tgt_db, _ref_features("gene1"), _args())

    assert [call["ref_trans_id"] for call in calls] == ["good"]
    assert score.startswith("good\t") and "bad" not in score


def test_evaluation_keeps_the_cds_a_nested_stop_codon_row_duplicates(monkeypatch):
    """A ``stop_codon`` row inside the terminal CDS -- miniprot writes one, and
    older LiftOn outputs carry them -- used to replace that CDS on its exon,
    so the evaluated protein was a 3-bp fragment."""
    from lifton import lifton_class

    ref_db = _db([
        _row("gene", 1, 90, "ID=gene1;gene_biotype=protein_coding"),
        _row("mRNA", 1, 90, "ID=tx1;Parent=gene1"),
    ])
    tgt_db = _db([
        _row("gene", 1, 90, "ID=gene1"),
        _row("mRNA", 1, 90, "ID=tx1;Parent=gene1"),
        _row("exon", 1, 90, "ID=ex1;Parent=tx1"),
        _row("CDS", 1, 90, "ID=cds1;Parent=tx1", frame="0"),
        _row("stop_codon", 88, 90, "ID=stop1;Parent=tx1", frame="0"),
    ])
    seen = []

    def fake_orf(self, trans_id, ref_trans_id, _fai, _proteins, _transcripts,
                 status, eval_only=False, eval_liftoff_chm13=False):
        seen.append([(e.cds.entry.start, e.cds.entry.end)
                     for e in self.transcripts[trans_id].exons if e.cds is not None])
        status.lifton_dna = status.lifton_aa = 1.0
        status.status = ["identical"]
        return None, None

    monkeypatch.setattr(lifton_class.Lifton_GENE, "orf_search_protein", fake_orf)
    _evaluate(ref_db, tgt_db, _ref_features("gene1"), _args())
    assert seen == [[(1, 90)]]


# ---------------------------------------------------------------------------
# A CDS inside one of two overlapping exons is not unsplittable
# ---------------------------------------------------------------------------

@pytest.fixture
def shared_base_workspace(unsplittable_workspace):
    """tx2's two Liftoff exons share one base (Liftoff writes such pairs;
    C. elegans -> C. briggsae W07G4.3), and each CDS row equals its exon. Each
    CDS lies wholly inside one exon and merely touches the other."""
    work = unsplittable_workspace
    text = (work / "liftoff.gff3").read_text()
    text = text.replace(
        "chr1\tLiftoff\texon\t101\t160\t.\t+\t.\tID=exon2a;Parent=tx2\n"
        "chr1\tLiftoff\texon\t150\t199\t.\t+\t.\tID=exon2b;Parent=tx2\n"
        "chr1\tLiftoff\tCDS\t101\t199\t.\t+\t0\tID=cds2;Parent=tx2\n",
        "chr1\tLiftoff\texon\t101\t150\t.\t+\t.\tID=exon2a;Parent=tx2\n"
        "chr1\tLiftoff\texon\t150\t199\t.\t+\t.\tID=exon2b;Parent=tx2\n"
        "chr1\tLiftoff\tCDS\t101\t150\t.\t+\t0\tID=cds2;Parent=tx2\n"
        "chr1\tLiftoff\tCDS\t150\t199\t.\t+\t2\tID=cds2;Parent=tx2\n")
    assert "101\t150" in text
    (work / "liftoff.gff3").write_text(text)
    return work


def test_add_cds_attaches_a_contained_cds_to_its_own_exon():
    from types import SimpleNamespace
    from lifton import lifton_class

    def feature(kind, start, end, fid, frame="."):
        return SimpleNamespace(seqid="chr1", strand="+", start=start, end=end, id=fid,
                               featuretype=kind, frame=frame,
                               attributes={"ID": [fid], "Parent": ["tx"]})

    trans = lifton_class.Lifton_TRANS.__new__(lifton_class.Lifton_TRANS)
    trans.entry = feature("mRNA", 101, 199, "tx")
    trans._transl_table = 1
    trans._cds_attr_template = None
    trans.exons = [lifton_class.Lifton_EXON(feature("exon", 101, 150, "e1")),
                   lifton_class.Lifton_EXON(feature("exon", 150, 199, "e2"))]
    trans.add_cds(feature("CDS", 101, 150, "c", "0"))
    trans.add_cds(feature("CDS", 150, 199, "c", "2"))
    assert [(e.cds.entry.start, e.cds.entry.end) for e in trans.exons] == [(101, 150), (150, 199)]


@pytest.mark.parametrize("threads", [(), ("-t", "2", "--locus-pipeline")])
def test_a_cds_inside_one_of_two_overlapping_exons_keeps_its_transcript(
        shared_base_workspace, hermetic_pipeline, threads):
    text, manifest = _run(shared_base_workspace, "out.gff3", *threads)
    rows = [line.split("\t") for line in text.splitlines()
            if line and not line.startswith("#")]
    mrnas = {r[8].split("ID=")[1].split(";")[0] for r in rows if r[2] == "mRNA"}
    assert mrnas == {"tx1", "tx2"}
    assert manifest["counts"]["dropped_cds_spanning_exons"] == 0
    assert manifest["run"]["status"] == "success"
