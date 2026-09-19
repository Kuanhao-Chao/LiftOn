"""Isoform-aware miniprot-only rescue (v1.0.12 A2).

The rescue deduplicates by reference gene, so it emitted one transcript per
rescued gene. With ``--rescue-isoforms`` each rescued gene also receives the
other transcripts of the same reference gene whose miniprot hits sit at its
locus, after every gene has been placed. These tests pin that the isoform is
added, that placement never changes, that a widening into another gene is
refused (leaving the output byte-identical to OFF), and that the identity floor
applies to isoforms too.
"""
from __future__ import annotations

import errno
import textwrap
import types

import pytest
from intervaltree import IntervalTree

from lifton import lifton, miniprot_rescue
from lifton.intervals import _make_interval
from tests.test_integration_pipeline import hermetic_pipeline  # noqa: F401


GENE1_EXON1 = "ATG" + "GCT" * 32
GENE1_EXON2 = "GCT" * 32 + "TAA"
GENE2_EXON1 = "ATG" + "GCT" * 16          # 51 nt
GENE2_EXON2 = "GCT" * 15 + "TAA"          # 48 nt: M + 31 A + stop
GENE3 = "ATG" + "GCT" * 8 + "TAA"         # 30 nt ORF at 1300-1329


def _wrap(sequence):
    return "\n".join(textwrap.wrap(sequence, 60)) + "\n"


def _chromosome(placements, length=1500):
    chromosome = ["A"] * length
    for start, sequence in placements:
        chromosome[start - 1:start - 1 + len(sequence)] = sequence
    return "".join(chromosome)


def _build_workspace(work, *, neighbour=False, isoform_exon2=GENE2_EXON2):
    """gene2 has isoforms tx2a (last exon 801-848) and tx2b (951-998) in the
    reference; its CDS span is 398. In the target the two last exons sit at
    1231-1278 and 1331-1378, so miniprot's hits span 678 and 778 (ratios 1.70
    and 1.96): outside Step 8's band, inside the rescue band. Liftoff lifts
    gene1 (and gene3 at 1300-1329 when ``neighbour``) but not gene2."""
    work.mkdir(parents=True, exist_ok=True)
    extra = [(1300, GENE3)] if neighbour else []
    (work / "ref.fa").write_text(">chr1\n" + _wrap(_chromosome([
        (101, GENE1_EXON1), (301, GENE1_EXON2), (601, GENE2_EXON1),
        (801, GENE2_EXON2), (951, GENE2_EXON2), *extra])))
    (work / "tgt.fa").write_text(">chr1\n" + _wrap(_chromosome([
        (101, GENE1_EXON1), (301, GENE1_EXON2), (601, GENE2_EXON1),
        (1231, GENE2_EXON2), (1331, isoform_exon2), *extra])))
    gene1 = (
        "chr1\t{src}\tgene\t101\t399\t.\t+\t.\tID=gene1;gene_biotype=protein_coding\n"
        "chr1\t{src}\tmRNA\t101\t399\t.\t+\t.\tID=tx1;Parent=gene1\n"
        "chr1\t{src}\texon\t101\t199\t.\t+\t.\tID=exon1;Parent=tx1\n"
        "chr1\t{src}\texon\t301\t399\t.\t+\t.\tID=exon2;Parent=tx1\n"
        "chr1\t{src}\tCDS\t101\t199\t.\t+\t0\tID=cds1;Parent=tx1\n"
        "chr1\t{src}\tCDS\t301\t399\t.\t+\t0\tID=cds1;Parent=tx1\n")
    gene3 = (
        "chr1\t{src}\tgene\t1300\t1329\t.\t+\t.\tID=gene3;gene_biotype=protein_coding\n"
        "chr1\t{src}\tmRNA\t1300\t1329\t.\t+\t.\tID=tx3;Parent=gene3\n"
        "chr1\t{src}\texon\t1300\t1329\t.\t+\t.\tID=exon9;Parent=tx3\n"
        "chr1\t{src}\tCDS\t1300\t1329\t.\t+\t0\tID=cds3;Parent=tx3\n"
        if neighbour else "")
    gene2 = (
        "chr1\ttest\tgene\t601\t998\t.\t+\t.\tID=gene2;gene_biotype=protein_coding\n"
        "chr1\ttest\tmRNA\t601\t848\t.\t+\t.\tID=tx2a;Parent=gene2\n"
        "chr1\ttest\texon\t601\t651\t.\t+\t.\tID=exon3;Parent=tx2a\n"
        "chr1\ttest\texon\t801\t848\t.\t+\t.\tID=exon4;Parent=tx2a\n"
        "chr1\ttest\tCDS\t601\t651\t.\t+\t0\tID=cds2a;Parent=tx2a\n"
        "chr1\ttest\tCDS\t801\t848\t.\t+\t0\tID=cds2a;Parent=tx2a\n"
        "chr1\ttest\tmRNA\t601\t998\t.\t+\t.\tID=tx2b;Parent=gene2\n"
        "chr1\ttest\texon\t601\t651\t.\t+\t.\tID=exon5;Parent=tx2b\n"
        "chr1\ttest\texon\t951\t998\t.\t+\t.\tID=exon6;Parent=tx2b\n"
        "chr1\ttest\tCDS\t601\t651\t.\t+\t0\tID=cds2b;Parent=tx2b\n"
        "chr1\ttest\tCDS\t951\t998\t.\t+\t0\tID=cds2b;Parent=tx2b\n")
    (work / "ref.gff3").write_text("##gff-version 3\n" + gene1.format(src="test")
                                   + gene2 + gene3.format(src="test"))
    (work / "liftoff.gff3").write_text(
        "##gff-version 3\n" + gene1.format(src="Liftoff")
        + gene3.format(src="Liftoff"))
    (work / "miniprot.gff3").write_text(
        "##gff-version 3\n"
        "chr1\tminiprot\tmRNA\t101\t399\t.\t+\t.\tID=MP1;Rank=1;Identity=1.0000;Target=tx1 1 65\n"
        "chr1\tminiprot\tCDS\t101\t199\t.\t+\t0\tID=MP1.c1;Parent=MP1\n"
        "chr1\tminiprot\tCDS\t301\t399\t.\t+\t0\tID=MP1.c2;Parent=MP1\n"
        "chr1\tminiprot\tmRNA\t601\t1278\t.\t+\t.\tID=MP2;Rank=1;Identity=1.0000;Target=tx2a 1 32\n"
        "chr1\tminiprot\tCDS\t601\t651\t.\t+\t0\tID=MP2.c1;Parent=MP2\n"
        "chr1\tminiprot\tCDS\t1231\t1278\t.\t+\t0\tID=MP2.c2;Parent=MP2\n"
        "chr1\tminiprot\tmRNA\t601\t1378\t.\t+\t.\tID=MP3;Rank=1;Identity=1.0000;Target=tx2b 1 32\n"
        "chr1\tminiprot\tCDS\t601\t651\t.\t+\t0\tID=MP3.c1;Parent=MP3\n"
        "chr1\tminiprot\tCDS\t1331\t1378\t.\t+\t0\tID=MP3.c2;Parent=MP3\n")
    (work / "out").mkdir()
    return work


def _run(work, *flags):
    output = work / "out" / "lifton.gff3"
    argv = [str(work / "tgt.fa"), str(work / "ref.fa"), "-g", str(work / "ref.gff3"),
            "-L", str(work / "liftoff.gff3"), "-M", str(work / "miniprot.gff3"),
            "-o", str(output), "-ad", "RefSeq", "--force", *flags]
    lifton.run_all_lifton_steps(lifton.parse_args(argv))
    return output.read_text()


def _rows(body, feature_type):
    return [row.split("\t") for row in body.splitlines()
            if row and not row.startswith("#") and row.split("\t")[2] == feature_type]


def _attribute(row, key):
    for field in row[8].split(";"):
        if field.startswith(key + "="):
            return field.split("=", 1)[1]
    return None


@pytest.fixture(autouse=True)
def _clean_env(monkeypatch):
    for key in ("LIFTON_RESCUE_ISOFORMS", "LIFTON_RESCUE_COVERAGE_GATE",
                "LIFTON_MINIPROT_RESCUE", "LIFTON_MINIPROT_RESCUE_MIN_ID",
                "LIFTON_MINIPROT_RESCUE_LEN", "LIFTON_RESCUE_ADAPTIVE_FLOOR"):
        monkeypatch.delenv(key, raising=False)


class TestIsoformPipeline:
    def test_off_emits_one_transcript(self, tmp_path, hermetic_pipeline):
        body = _run(_build_workspace(tmp_path / "w"), "--no-rescue-isoforms")
        gene2 = [r for r in _rows(body, "mRNA") if _attribute(r, "Parent") == "gene2"]
        assert [_attribute(r, "ID") for r in gene2] == ["tx2a"]

    def test_on_attaches_the_second_isoform(self, tmp_path, hermetic_pipeline):
        off = _run(_build_workspace(tmp_path / "off"), "--no-rescue-isoforms")
        on = _run(_build_workspace(tmp_path / "on"), "--rescue-isoforms")
        mrnas = [r for r in _rows(on, "mRNA") if _attribute(r, "Parent") == "gene2"]
        assert sorted(_attribute(r, "ID") for r in mrnas) == ["tx2a", "tx2b"]
        isoform = next(r for r in mrnas if _attribute(r, "ID") == "tx2b")
        assert _attribute(isoform, "rescue_isoform") == "true"
        assert _attribute(isoform, "lifton_rescue") == "miniprot_only"
        assert _attribute(isoform, "protein_identity") == "1.000"
        assert (isoform[3], isoform[4]) == ("601", "1378")
        # Placement is unchanged: the same genes, and gene2 only grew.
        assert ([_attribute(r, "ID") for r in _rows(off, "gene")]
                == [_attribute(r, "ID") for r in _rows(on, "gene")])
        gene2_off = next(r for r in _rows(off, "gene") if _attribute(r, "ID") == "gene2")
        gene2_on = next(r for r in _rows(on, "gene") if _attribute(r, "ID") == "gene2")
        assert (gene2_off[3], gene2_off[4]) == ("601", "1278")
        assert (gene2_on[3], gene2_on[4]) == ("601", "1378")
        # Every transcript of the OFF output is still there.
        off_ids = {_attribute(r, "ID") for r in _rows(off, "mRNA")}
        on_ids = {_attribute(r, "ID") for r in _rows(on, "mRNA")}
        assert off_ids < on_ids

    def test_environment_enables_it(self, tmp_path, hermetic_pipeline, monkeypatch):
        monkeypatch.setenv("LIFTON_RESCUE_ISOFORMS", "1")
        body = _run(_build_workspace(tmp_path / "w"))
        assert "rescue_isoform=true" in body

    def test_widening_into_a_neighbour_is_refused(self, tmp_path, hermetic_pipeline):
        # tx2b's hit would stretch gene2 over gene3 (1300-1329): not attached,
        # and with nothing attached the ON output is byte-identical to OFF.
        off = _run(_build_workspace(tmp_path / "off", neighbour=True),
                   "--no-rescue-isoforms")
        on = _run(_build_workspace(tmp_path / "on", neighbour=True),
                  "--rescue-isoforms")
        assert "ID=gene3" in on and "ID=tx2a" in on
        assert "ID=tx2b" not in on
        assert on == off

    def test_isoforms_meet_the_identity_floor(self, tmp_path, hermetic_pipeline,
                                              monkeypatch):
        monkeypatch.setenv("LIFTON_MINIPROT_RESCUE_MIN_ID", "0.95")
        monkeypatch.setenv("LIFTON_RESCUE_ADAPTIVE_FLOOR", "0")
        work = _build_workspace(tmp_path / "w",
                                isoform_exon2="CCC" * 15 + "TAA")
        body = _run(work, "--rescue-isoforms")
        assert "ID=tx2a" in body
        assert "ID=tx2b" not in body


class TestExtensionCollides:
    def _tree(self):
        tree = IntervalTree()
        tree.add(_make_interval(100, 200, "self"))
        tree.add(_make_interval(300, 400, "other"))
        return {"chr1": tree}

    def test_growth_inside_free_space(self):
        assert not miniprot_rescue._extension_collides(
            self._tree(), "chr1", "self", 100, 200, 90, 299)

    def test_growth_into_another_gene(self):
        assert miniprot_rescue._extension_collides(
            self._tree(), "chr1", "self", 100, 200, 100, 300)

    def test_no_growth_never_collides(self):
        assert not miniprot_rescue._extension_collides(
            self._tree(), "chr1", "self", 100, 200, 100, 200)

    def test_unknown_sequence(self):
        assert not miniprot_rescue._extension_collides(
            {}, "chr9", "self", 1, 2, 1, 50)


class TestResolution:
    def _parse(self, *extra):
        return lifton.parse_args(["t.fa", "r.fa", "-g", "r.gff3", "-o", "o.gff3",
                                  *extra])

    def test_flags_and_default(self, monkeypatch):
        unset = lifton.resolve_miniprot_rescue_args(self._parse())
        assert unset.rescue_isoforms is None
        for default in (False, True):
            monkeypatch.setattr(miniprot_rescue, "ISOFORMS_DEFAULT", default)
            assert miniprot_rescue._rescue_isoforms_on(unset) is default
        assert miniprot_rescue._rescue_isoforms_on(
            lifton.resolve_miniprot_rescue_args(self._parse("--rescue-isoforms")))
        assert not miniprot_rescue._rescue_isoforms_on(
            lifton.resolve_miniprot_rescue_args(self._parse("--no-rescue-isoforms")))

    def test_environment_wins(self, monkeypatch):
        monkeypatch.setenv("LIFTON_RESCUE_ISOFORMS", "0")
        args = lifton.resolve_miniprot_rescue_args(self._parse("--rescue-isoforms"))
        assert not miniprot_rescue._rescue_isoforms_on(args)
        assert not miniprot_rescue._rescue_isoforms_on(types.SimpleNamespace())


class TestIsoformWorkers:
    def test_process_pool_matches_inline_scoring(self, tmp_path, hermetic_pipeline,
                                                 monkeypatch):
        monkeypatch.setenv("LIFTON_RESCUE_ISOFORM_WORKERS", "0")
        inline = _run(_build_workspace(tmp_path / "inline"), "--rescue-isoforms")
        monkeypatch.setenv("LIFTON_RESCUE_ISOFORM_WORKERS", "2")
        pooled = _run(_build_workspace(tmp_path / "pooled"), "--rescue-isoforms")
        assert "rescue_isoform=true" in inline
        assert pooled == inline

    def test_enomem_starting_the_pool_falls_back_to_inline_scoring(
            self, tmp_path, hermetic_pipeline, monkeypatch):
        # fork() can fail with memory free: under strict overcommit the kernel
        # charges each child the parent's whole address space. Losing the pool
        # must cost speed, not the rescue.
        monkeypatch.setenv("LIFTON_RESCUE_ISOFORM_WORKERS", "0")
        inline = _run(_build_workspace(tmp_path / "inline"), "--rescue-isoforms")
        assert "rescue_isoform=true" in inline

        import multiprocessing

        def _enomem_pool(*_args, **_kwargs):
            raise OSError(errno.ENOMEM, "Cannot allocate memory")

        real_get_context = multiprocessing.get_context

        def _fake_get_context(*args, **kwargs):
            context = real_get_context(*args, **kwargs)
            return types.SimpleNamespace(Pool=_enomem_pool, _real=context)

        monkeypatch.setattr(multiprocessing, "get_context", _fake_get_context)
        monkeypatch.setenv("LIFTON_RESCUE_ISOFORM_WORKERS", "2")
        fallback = _run(_build_workspace(tmp_path / "fallback"), "--rescue-isoforms")
        assert fallback == inline
        assert miniprot_rescue._ISOFORM_SHARED == {}

    def test_worker_count(self, monkeypatch):
        monkeypatch.delenv("LIFTON_RESCUE_ISOFORM_WORKERS", raising=False)
        args = types.SimpleNamespace(threads=8)
        assert miniprot_rescue._isoform_workers(args, 10) == 0      # too little work
        assert miniprot_rescue._isoform_workers(args, 500) == 8
        assert miniprot_rescue._isoform_workers(types.SimpleNamespace(), 500) == 1
        monkeypatch.setenv("LIFTON_RESCUE_ISOFORM_WORKERS", "3")
        assert miniprot_rescue._isoform_workers(args, 10) == 3
