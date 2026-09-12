"""Terminal-stop completion of miniprot-derived models (v1.0.12).

miniprot's CDS ends at the last aligned codon, so a miniprot-derived model is
one codon short of the reference convention every other LiftOn model follows —
and the ORF search cannot add it, because such a model has no UTR to search.
``lifton.orf_completion`` grows the terminal CDS and its exon over a downstream
stop codon. These tests pin when it fires, when it refuses, that it is
strand-correct, and that it reaches a real rescued model end to end.
"""
from __future__ import annotations

import textwrap
import types

import gffutils
import pytest
from Bio.Seq import Seq
from pyfaidx import Fasta

from lifton import lifton, orf_completion
from tests.test_integration_pipeline import hermetic_pipeline  # noqa: F401


def _feature(seqid, start, end, strand, featuretype="CDS", identifier="f"):
    return gffutils.feature.feature_from_line(
        f"{seqid}\tt\t{featuretype}\t{start}\t{end}\t.\t{strand}\t0\tID={identifier}")


def _transcript(seqid, blocks, strand):
    """A CDS-only model: one exon per block, each exon flush with its CDS."""
    exons = []
    for index, (start, end) in enumerate(blocks):
        exons.append(types.SimpleNamespace(
            entry=_feature(seqid, start, end, strand, "exon", f"e{index}"),
            cds=types.SimpleNamespace(
                entry=_feature(seqid, start, end, strand, "CDS", f"c{index}"))))
    exons.sort(key=lambda e: e.entry.start)
    return types.SimpleNamespace(
        entry=_feature(seqid, exons[0].entry.start, exons[-1].entry.end, strand,
                       "mRNA", "tx"),
        exons=exons)


@pytest.fixture
def genome(tmp_path):
    """chr1 positions 1-99 are an ORF-free filler; the model occupies 10-33 and
    the codon at 34-36 is TAA. On the minus strand the model occupies 40-63 and
    the codon at 37-39 is TTA, the reverse complement of a stop."""
    sequence = list("A" * 120)
    sequence[9:33] = "ATG" + "GCT" * 7
    sequence[33:36] = "TAA"
    sequence[36:39] = "TTA"
    sequence[39:63] = list(str(Seq("ATG" + "GCT" * 7).reverse_complement()))
    path = tmp_path / "g.fa"
    path.write_text(">chr1\n" + "\n".join(textwrap.wrap("".join(sequence), 60)) + "\n")
    return Fasta(str(path))


class TestCompleteTerminalStop:
    def test_extends_over_a_downstream_stop(self, genome):
        transcript = _transcript("chr1", [(10, 33)], "+")
        assert orf_completion.complete_terminal_stop(transcript, genome) is True
        assert transcript.exons[-1].cds.entry.end == 36
        assert transcript.exons[-1].entry.end == 36

    def test_minus_strand_extends_downstream_in_transcript_order(self, genome):
        transcript = _transcript("chr1", [(40, 63)], "-")
        assert orf_completion.complete_terminal_stop(transcript, genome) is True
        assert transcript.exons[0].cds.entry.start == 37
        assert transcript.exons[0].entry.start == 37

    def test_multi_exon_extends_only_the_terminal_exon(self, genome):
        transcript = _transcript("chr1", [(10, 21), (22, 33)], "+")
        assert orf_completion.complete_terminal_stop(transcript, genome) is True
        assert [(e.entry.start, e.entry.end) for e in transcript.exons] == [
            (10, 21), (22, 36)]

    def test_refuses_when_the_next_codon_is_not_a_stop(self, genome):
        transcript = _transcript("chr1", [(10, 30)], "+")   # next codon is GCT
        assert orf_completion.complete_terminal_stop(transcript, genome) is False
        assert transcript.exons[-1].cds.entry.end == 30

    def test_refuses_when_the_model_already_ends_in_a_stop(self, genome):
        transcript = _transcript("chr1", [(10, 36)], "+")   # 34-36 is TAA
        assert orf_completion.complete_terminal_stop(transcript, genome) is False

    def test_refuses_a_partial_codon(self, genome):
        transcript = _transcript("chr1", [(10, 32)], "+")   # 23 bases
        assert orf_completion.complete_terminal_stop(transcript, genome) is False

    def test_refuses_when_the_cds_is_not_flush_with_its_exon(self, genome):
        transcript = _transcript("chr1", [(10, 33)], "+")
        transcript.exons[-1].entry.end = 45                 # the model has a UTR
        assert orf_completion.complete_terminal_stop(transcript, genome) is False

    def test_refuses_past_the_end_of_the_sequence(self, genome):
        transcript = _transcript("chr1", [(97, 120)], "+")
        assert orf_completion.complete_terminal_stop(transcript, genome) is False

    def test_refuses_a_model_without_coding_exons(self, genome):
        transcript = _transcript("chr1", [(10, 33)], "+")
        transcript.exons[0].cds = None
        assert orf_completion.complete_terminal_stop(transcript, genome) is False

    def test_refuses_an_unstranded_model(self, genome):
        transcript = _transcript("chr1", [(10, 33)], "+")
        transcript.entry.strand = "."
        assert orf_completion.complete_terminal_stop(transcript, genome) is False

    def test_refuses_an_unknown_sequence(self, genome):
        transcript = _transcript("chrZ", [(10, 33)], "+")
        assert orf_completion.complete_terminal_stop(transcript, genome) is False


class TestSwitch:
    @pytest.fixture(autouse=True)
    def _clean_env(self, monkeypatch):
        monkeypatch.delenv("LIFTON_ORF_STOP_COMPLETION", raising=False)

    def _args(self, *extra):
        return lifton.resolve_miniprot_rescue_args(lifton.parse_args(
            ["t.fa", "r.fa", "-g", "r.gff3", "-o", "o.gff3", *extra]))

    def test_default_is_on(self):
        assert orf_completion.enabled(self._args()) is True

    def test_flags(self):
        assert orf_completion.enabled(self._args("--no-orf-stop-completion")) is False
        assert orf_completion.enabled(self._args("--orf-stop-completion")) is True

    @pytest.mark.parametrize("value, expected", [("1", True), ("yes", True),
                                                 ("0", False), ("false", False),
                                                 ("OFF", False)])
    def test_environment_wins(self, monkeypatch, value, expected):
        monkeypatch.setenv("LIFTON_ORF_STOP_COMPLETION", value)
        flag = "--no-orf-stop-completion" if expected else "--orf-stop-completion"
        assert orf_completion.enabled(self._args(flag)) is expected

    def test_missing_attribute_defaults_on(self):
        assert orf_completion.enabled(types.SimpleNamespace()) is True


# --------------------------------------------------------------------------- #
# End to end (pre-baked Liftoff/miniprot; external tools never run)
# --------------------------------------------------------------------------- #
GENE1_EXON1 = "ATG" + "GCT" * 32
GENE1_EXON2 = "GCT" * 32 + "TAA"
GENE2_CDS = "ATG" + "GCT" * 31            # 96 nt, no stop: miniprot's convention
GENE2_REF = GENE2_CDS + "TAA"             # the reference CDS carries the stop


def _wrap(sequence):
    return "\n".join(textwrap.wrap(sequence, 60)) + "\n"


def _chromosome(placements, length=1200):
    chromosome = ["A"] * length
    for start, sequence in placements:
        chromosome[start - 1:start - 1 + len(sequence)] = sequence
    return "".join(chromosome)


def _build_workspace(work, *, target_stop=True):
    """Liftoff lifts gene1 and misses gene2. miniprot places gene2 at 601-696,
    ending at its last aligned codon; the target carries the stop codon at
    697-699 only when ``target_stop``."""
    work.mkdir(parents=True, exist_ok=True)
    (work / "ref.fa").write_text(">chr1\n" + _wrap(_chromosome([
        (101, GENE1_EXON1), (301, GENE1_EXON2), (601, GENE2_REF)])))
    target_gene2 = GENE2_CDS + ("TAA" if target_stop else "GCT")
    (work / "tgt.fa").write_text(">chr1\n" + _wrap(_chromosome([
        (101, GENE1_EXON1), (301, GENE1_EXON2), (601, target_gene2)])))
    gene1 = (
        "chr1\t{src}\tgene\t101\t399\t.\t+\t.\tID=gene1;gene_biotype=protein_coding\n"
        "chr1\t{src}\tmRNA\t101\t399\t.\t+\t.\tID=tx1;Parent=gene1\n"
        "chr1\t{src}\texon\t101\t199\t.\t+\t.\tID=exon1;Parent=tx1\n"
        "chr1\t{src}\texon\t301\t399\t.\t+\t.\tID=exon2;Parent=tx1\n"
        "chr1\t{src}\tCDS\t101\t199\t.\t+\t0\tID=cds1;Parent=tx1\n"
        "chr1\t{src}\tCDS\t301\t399\t.\t+\t0\tID=cds2;Parent=tx1\n")
    gene2 = (
        "chr1\ttest\tgene\t601\t699\t.\t+\t.\tID=gene2;gene_biotype=protein_coding\n"
        "chr1\ttest\tmRNA\t601\t699\t.\t+\t.\tID=tx2;Parent=gene2\n"
        "chr1\ttest\texon\t601\t699\t.\t+\t.\tID=exon3;Parent=tx2\n"
        "chr1\ttest\tCDS\t601\t699\t.\t+\t0\tID=cds3;Parent=tx2\n")
    (work / "ref.gff3").write_text("##gff-version 3\n" + gene1.format(src="test")
                                   + gene2)
    (work / "liftoff.gff3").write_text("##gff-version 3\n"
                                       + gene1.format(src="Liftoff"))
    (work / "miniprot.gff3").write_text(
        "##gff-version 3\n"
        "chr1\tminiprot\tmRNA\t101\t399\t.\t+\t.\tID=MP1;Rank=1;Identity=1.0000;Target=tx1 1 65\n"
        "chr1\tminiprot\tCDS\t101\t199\t.\t+\t0\tID=MP1.c1;Parent=MP1\n"
        "chr1\tminiprot\tCDS\t301\t399\t.\t+\t0\tID=MP1.c2;Parent=MP1\n"
        "chr1\tminiprot\tmRNA\t601\t696\t.\t+\t.\tID=MP2;Rank=1;Identity=1.0000;Target=tx2 1 32\n"
        "chr1\tminiprot\tCDS\t601\t696\t.\t+\t0\tID=MP2.c1;Parent=MP2\n")
    (work / "out").mkdir()
    return work


def _run(work, *flags):
    output = work / "out" / "lifton.gff3"
    argv = [str(work / "tgt.fa"), str(work / "ref.fa"), "-g", str(work / "ref.gff3"),
            "-L", str(work / "liftoff.gff3"), "-M", str(work / "miniprot.gff3"),
            "-o", str(output), "-ad", "RefSeq", "--force", *flags]
    lifton.run_all_lifton_steps(lifton.parse_args(argv))
    return output.read_text()


def _rows(body, feature_type, parent=None):
    rows = [r.split("\t") for r in body.splitlines()
            if r and not r.startswith("#") and r.split("\t")[2] == feature_type]
    if parent is None:
        return rows
    return [r for r in rows if f"Parent={parent}" in r[8]]


@pytest.fixture(autouse=True)
def _clean_rescue_env(monkeypatch):
    for key in ("LIFTON_ORF_STOP_COMPLETION", "LIFTON_RESCUE_ISOFORMS",
                "LIFTON_RESCUE_COVERAGE_GATE", "LIFTON_MINIPROT_RESCUE"):
        monkeypatch.delenv(key, raising=False)


class TestEndToEnd:
    def test_rescued_model_gains_the_stop_codon(self, tmp_path, hermetic_pipeline):
        on = _run(_build_workspace(tmp_path / "on"))
        cds = _rows(on, "CDS", parent="tx2")
        assert [(r[3], r[4]) for r in cds] == [("601", "699")]
        assert [(r[3], r[4]) for r in _rows(on, "exon", parent="tx2")] == [
            ("601", "699")]
        mrna = _rows(on, "mRNA", parent="gene2")[0]
        assert (mrna[3], mrna[4]) == ("601", "699")
        assert "orf_stop_completed=true" in mrna[8]

    def test_flag_off_keeps_the_miniprot_boundary(self, tmp_path, hermetic_pipeline):
        off = _run(_build_workspace(tmp_path / "off"), "--no-orf-stop-completion")
        assert [(r[3], r[4]) for r in _rows(off, "CDS", parent="tx2")] == [
            ("601", "696")]
        assert "orf_stop_completed" not in "".join(
            r[8] for r in _rows(off, "mRNA", parent="gene2"))

    def test_no_downstream_stop_leaves_the_model_alone(self, tmp_path,
                                                       hermetic_pipeline):
        body = _run(_build_workspace(tmp_path / "nostop", target_stop=False))
        assert [(r[3], r[4]) for r in _rows(body, "CDS", parent="tx2")] == [
            ("601", "696")]

    def test_lifted_genes_are_untouched(self, tmp_path, hermetic_pipeline):
        on = _run(_build_workspace(tmp_path / "a"))
        off = _run(_build_workspace(tmp_path / "b"), "--no-orf-stop-completion")

        def gene1_block(body):
            return [r for r in body.splitlines()
                    if "gene1" in r or "tx1" in r or "exon1" in r or "exon2" in r]

        assert gene1_block(on) == gene1_block(off)


# --------------------------------------------------------------------------- #
# The A/B gate that proves no other coordinate moved
# --------------------------------------------------------------------------- #
from benchmarks.compare.rescue_extension_ab import (  # noqa: E402
    _is_terminal_stop_extension)


class TestStopExtensionShape:
    PLUS = ("+", 100, 200, (("CDS", 100, 150), ("CDS", 181, 200),
                            ("exon", 100, 150), ("exon", 181, 200)))
    MINUS = ("-", 100, 200, PLUS[3])

    def test_accepts_a_plus_strand_extension(self):
        after = ("+", 100, 203, (("CDS", 100, 150), ("CDS", 181, 203),
                                 ("exon", 100, 150), ("exon", 181, 203)))
        assert _is_terminal_stop_extension(self.PLUS, after) is True

    def test_accepts_a_minus_strand_extension(self):
        after = ("-", 97, 200, (("CDS", 97, 150), ("CDS", 181, 200),
                                ("exon", 97, 150), ("exon", 181, 200)))
        assert _is_terminal_stop_extension(self.MINUS, after) is True

    def test_accepts_a_single_exon_model(self):
        before = ("+", 10, 33, (("CDS", 10, 33), ("exon", 10, 33)))
        after = ("+", 10, 36, (("CDS", 10, 36), ("exon", 10, 36)))
        assert _is_terminal_stop_extension(before, after) is True

    def test_rejects_the_wrong_exon_moving(self):
        after = ("+", 100, 203, (("CDS", 100, 153), ("CDS", 181, 200),
                                 ("exon", 100, 153), ("exon", 181, 200)))
        assert _is_terminal_stop_extension(self.PLUS, after) is False

    def test_rejects_more_than_one_codon(self):
        after = ("+", 100, 206, (("CDS", 100, 150), ("CDS", 181, 206),
                                 ("exon", 100, 150), ("exon", 181, 206)))
        assert _is_terminal_stop_extension(self.PLUS, after) is False

    def test_rejects_a_cds_that_outgrew_its_exon(self):
        after = ("+", 100, 203, (("CDS", 100, 150), ("CDS", 181, 203),
                                 ("exon", 100, 150), ("exon", 181, 200)))
        assert _is_terminal_stop_extension(self.PLUS, after) is False

    def test_rejects_a_changed_strand(self):
        after = ("-", 100, 203, (("CDS", 100, 150), ("CDS", 181, 203),
                                 ("exon", 100, 150), ("exon", 181, 203)))
        assert _is_terminal_stop_extension(self.PLUS, after) is False

    def test_rejects_a_changed_child_count(self):
        after = ("+", 100, 203, (("CDS", 100, 150), ("CDS", 181, 203),
                                 ("exon", 100, 150)))
        assert _is_terminal_stop_extension(self.PLUS, after) is False
