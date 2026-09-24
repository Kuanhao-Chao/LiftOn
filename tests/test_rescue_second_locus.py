"""One reference gene, two target genes.

A whole-genome duplication gives the target two genes where the reference has
one. The miniprot-only rescue refuses the second, because Iteration 23 dedups
by reference gene: a reference gene already emitted is never rescued again.
That rule is what makes the pass strictly additive at an OCCUPIED locus, and it
is exactly wrong at a free one.

Measured against zebrafish's own GRCz11 annotation, it hides 2,051 real target
genes on human -> zebrafish -- 84.5 % of everything miniprot finds at a locus
nothing occupies (`notes/coortholog_recall_measurement_2026-09.md`).

These tests pin the three things that make relaxing it safe: it is inert when
off, it still refuses an occupied locus when on, and it is capped.
"""
from __future__ import annotations

import textwrap

import pytest

from lifton import lifton, miniprot_rescue
from tests.test_integration_pipeline import hermetic_pipeline  # noqa: F401


GENE1_EXON1 = "ATG" + "GCT" * 32
GENE1_EXON2 = "GCT" * 32 + "TAA"
GENE2 = "ATG" + "GCT" * 31 + "TAA"          # 99 nt, a clean ORF
# The second target copy carries an intron, so its genomic span is 178 nt
# against a 99 nt reference gene -- ratio 1.80, outside Step 8's tight
# 0.9-1.5 band and inside the rescue's wider 0.5-2.0 one. Without that the hit
# never reaches the rescue at all: Step 8 emits it directly.
GENE2_COPY_EXON1 = GENE2[:51]
GENE2_COPY_EXON2 = GENE2[51:]
COPY_START, COPY_EXON2_START = 1101, 1231
COPY_END = COPY_EXON2_START + len(GENE2_COPY_EXON2) - 1


def _wrap(sequence):
    return "\n".join(textwrap.wrap(sequence, 60)) + "\n"


def _chromosome(placements, length=1500):
    chromosome = ["A"] * length
    for start, sequence in placements:
        chromosome[start - 1:start - 1 + len(sequence)] = sequence
    return "".join(chromosome)


def _build_workspace(work, *, occupy_second_locus=False):
    """gene2 exists once in the reference at 601-699 and twice in the target:
    the DNA lift places it at 601-699, and miniprot also finds it at 1101-1199,
    where nothing else is -- unless ``occupy_second_locus``, which puts a
    lifted gene3 there so the rescue must refuse even when it is allowed."""
    work.mkdir(parents=True, exist_ok=True)
    shared = [(101, GENE1_EXON1), (301, GENE1_EXON2), (601, GENE2)]
    target = shared + [(COPY_START, GENE2_COPY_EXON1),
                       (COPY_EXON2_START, GENE2_COPY_EXON2)]
    (work / "ref.fa").write_text(">chr1\n" + _wrap(_chromosome(shared)))
    (work / "tgt.fa").write_text(">chr1\n" + _wrap(_chromosome(target)))

    gene1 = (
        "chr1\t{src}\tgene\t101\t399\t.\t+\t.\tID=gene1;gene_biotype=protein_coding\n"
        "chr1\t{src}\tmRNA\t101\t399\t.\t+\t.\tID=tx1;Parent=gene1\n"
        "chr1\t{src}\texon\t101\t199\t.\t+\t.\tID=exon1;Parent=tx1\n"
        "chr1\t{src}\texon\t301\t399\t.\t+\t.\tID=exon2;Parent=tx1\n"
        "chr1\t{src}\tCDS\t101\t199\t.\t+\t0\tID=cds1;Parent=tx1\n"
        "chr1\t{src}\tCDS\t301\t399\t.\t+\t0\tID=cds1;Parent=tx1\n")
    gene2 = (
        "chr1\t{src}\tgene\t601\t699\t.\t+\t.\tID=gene2;gene_biotype=protein_coding\n"
        "chr1\t{src}\tmRNA\t601\t699\t.\t+\t.\tID=tx2;Parent=gene2\n"
        "chr1\t{src}\texon\t601\t699\t.\t+\t.\tID=exon3;Parent=tx2\n"
        "chr1\t{src}\tCDS\t601\t699\t.\t+\t0\tID=cds2;Parent=tx2\n")
    # gene3 exists only to occupy the second locus in the OCCUPIED variant.
    gene3 = ((
        f"chr1\t{{src}}\tgene\t{COPY_START}\t{COPY_END}\t.\t+\t.\t"
        f"ID=gene3;gene_biotype=protein_coding\n"
        f"chr1\t{{src}}\tmRNA\t{COPY_START}\t{COPY_END}\t.\t+\t.\t"
        f"ID=tx3;Parent=gene3\n"
        f"chr1\t{{src}}\texon\t{COPY_START}\t{COPY_START + 50}\t.\t+\t.\t"
        f"ID=exon4;Parent=tx3\n"
        f"chr1\t{{src}}\texon\t{COPY_EXON2_START}\t{COPY_END}\t.\t+\t.\t"
        f"ID=exon5;Parent=tx3\n"
        f"chr1\t{{src}}\tCDS\t{COPY_START}\t{COPY_START + 50}\t.\t+\t0\t"
        f"ID=cds3;Parent=tx3\n"
        f"chr1\t{{src}}\tCDS\t{COPY_EXON2_START}\t{COPY_END}\t.\t+\t0\t"
        f"ID=cds3;Parent=tx3\n")
        if occupy_second_locus else "")

    (work / "ref.gff3").write_text(
        "##gff-version 3\n" + gene1.format(src="test") + gene2.format(src="test")
        + gene3.format(src="test"))
    (work / "liftoff.gff3").write_text(
        "##gff-version 3\n" + gene1.format(src="Liftoff")
        + gene2.format(src="Liftoff") + gene3.format(src="Liftoff"))
    (work / "miniprot.gff3").write_text(
        "##gff-version 3\n"
        "chr1\tminiprot\tmRNA\t101\t399\t.\t+\t.\tID=MP1;Rank=1;Identity=1.0000;"
        "Target=tx1 1 65\n"
        "chr1\tminiprot\tCDS\t101\t199\t.\t+\t0\tID=MP1.c1;Parent=MP1\n"
        "chr1\tminiprot\tCDS\t301\t399\t.\t+\t0\tID=MP1.c2;Parent=MP1\n"
        # The hit at the locus the DNA lift already holds.
        "chr1\tminiprot\tmRNA\t601\t699\t.\t+\t.\tID=MP2;Rank=1;Identity=1.0000;"
        "Target=tx2 1 32\n"
        "chr1\tminiprot\tCDS\t601\t699\t.\t+\t0\tID=MP2.c1;Parent=MP2\n"
        # The SECOND target copy of the same reference gene.
        f"chr1\tminiprot\tmRNA\t{COPY_START}\t{COPY_END}\t.\t+\t.\t"
        f"ID=MP3;Rank=1;Identity=1.0000;Target=tx2 1 32\n"
        f"chr1\tminiprot\tCDS\t{COPY_START}\t{COPY_START + 50}\t.\t+\t0\t"
        f"ID=MP3.c1;Parent=MP3\n"
        f"chr1\tminiprot\tCDS\t{COPY_EXON2_START}\t{COPY_END}\t.\t+\t0\t"
        f"ID=MP3.c2;Parent=MP3\n")
    (work / "out").mkdir()
    return work


def _run(work, *flags):
    output = work / "out" / "lifton.gff3"
    lifton.run_all_lifton_steps(lifton.parse_args([
        str(work / "tgt.fa"), str(work / "ref.fa"),
        "-g", str(work / "ref.gff3"), "-L", str(work / "liftoff.gff3"),
        "-M", str(work / "miniprot.gff3"), "-o", str(output),
        "-ad", "RefSeq", "--force", *flags]))
    return output.read_text()


def _genes(text):
    out = []
    for line in text.splitlines():
        fields = line.split("\t")
        if len(fields) > 8 and fields[2] == "gene":
            attrs = dict(kv.split("=", 1) for kv in fields[8].split(";")
                         if "=" in kv)
            out.append((attrs.get("ID"), int(fields[3]), int(fields[4])))
    return out


@pytest.fixture(autouse=True)
def _clean_env(monkeypatch):
    for name in ("LIFTON_RESCUE_SECOND_LOCUS", "LIFTON_RESCUE_SECOND_LOCUS_MAX",
                 "LIFTON_MINIPROT_RESCUE", "LIFTON_RESCUE_ISOFORMS"):
        monkeypatch.delenv(name, raising=False)


class TestSwitch:
    def _args(self, *flags):
        return lifton.parse_args(["t.fa", "r.fa", "-g", "r.gff3", *flags])

    def test_on_by_default(self):
        # PROMOTED: +690 real GRCz11 genes on human to zebrafish with 0 lost at
        # gene and transcript level, safety gate 8/8 across the ladder.
        assert miniprot_rescue._second_locus_on(self._args())

    def test_the_opt_out_turns_it_off_and_the_alias_is_a_no_op(self):
        assert not miniprot_rescue._second_locus_on(
            self._args("--no-rescue-second-locus"))
        assert miniprot_rescue._second_locus_on(
            self._args("--rescue-second-locus"))

    def test_the_cap_defaults_to_one(self):
        assert miniprot_rescue._second_locus_max(self._args()) == 1

    def test_the_environment_wins(self, monkeypatch):
        monkeypatch.setenv("LIFTON_RESCUE_SECOND_LOCUS", "0")
        assert not miniprot_rescue._second_locus_on(self._args())
        assert not miniprot_rescue._second_locus_on(
            self._args("--rescue-second-locus"))
        monkeypatch.setenv("LIFTON_RESCUE_SECOND_LOCUS", "1")
        assert miniprot_rescue._second_locus_on(
            self._args("--no-rescue-second-locus"))

    def test_cross_locus_enabled_by_environment_also_takes_precedence(self, monkeypatch):
        """Both passes answer the same question with opposite policies, so the
        explicitly requested cross-locus rescue wins -- whether it was asked for
        by flag or, until now ignored here, by LIFTON_CROSS_LOCUS_RESCUE."""
        monkeypatch.setenv("LIFTON_CROSS_LOCUS_RESCUE", "1")
        assert not miniprot_rescue._second_locus_on(self._args())
        monkeypatch.setenv("LIFTON_CROSS_LOCUS_RESCUE", "0")
        assert miniprot_rescue._second_locus_on(self._args())

    def test_the_cap_gates_each_reference_gene_separately(self):
        args = self._args("--rescue-second-locus")
        counts = {}
        assert miniprot_rescue._second_locus_allowed("gene-A", counts, args)
        counts["gene-A"] = 1
        assert not miniprot_rescue._second_locus_allowed("gene-A", counts, args)
        assert miniprot_rescue._second_locus_allowed("gene-B", counts, args)

    def test_a_cap_of_zero_is_the_same_as_off(self):
        args = self._args("--rescue-second-locus", "--rescue-second-locus-max", "0")
        assert not miniprot_rescue._second_locus_allowed("gene-A", {}, args)


class TestFreeLocus:
    def test_the_opt_out_emits_one_model_for_the_reference_gene(
            self, tmp_path, hermetic_pipeline):
        genes = _genes(_run(_build_workspace(tmp_path / "off"),
                            "--no-rescue-second-locus"))
        gene2 = [g for g in genes if g[0].startswith("gene2")]
        assert [(g[1], g[2]) for g in gene2] == [(601, 699)]

    def test_on_places_the_second_target_copy_and_tags_it(
            self, tmp_path, hermetic_pipeline):
        off = _run(_build_workspace(tmp_path / "off"), "--no-rescue-second-locus")
        on = _run(_build_workspace(tmp_path / "on"))
        off_genes, on_genes = _genes(off), _genes(on)

        # off subset of on: nothing the default emitted may move or disappear.
        assert set(off_genes) <= set(on_genes)
        added = set(on_genes) - set(off_genes)
        assert len(added) == 1
        added_id, start, end = added.pop()
        assert (start, end) == (COPY_START, COPY_END)
        assert added_id.startswith("gene2")
        # It must be distinguishable in the output.
        assert "lifton_rescue_second_locus=true" in on
        assert "lifton_rescue_second_locus" not in off

    def test_the_first_model_keeps_its_own_id(self, tmp_path, hermetic_pipeline):
        on = _genes(_run(_build_workspace(tmp_path / "on")))
        first = [g for g in on if (g[1], g[2]) == (601, 699)]
        assert [g[0] for g in first] == ["gene2"]

    def test_the_cap_bounds_how_many_extra_loci_a_gene_gets(
            self, tmp_path, hermetic_pipeline):
        capped = _genes(_run(_build_workspace(tmp_path / "capped"),
                             "--rescue-second-locus-max", "0"))
        assert [(g[1], g[2]) for g in capped if g[0].startswith("gene2")] \
            == [(601, 699)]


class TestOccupiedLocusIsStillRefused:
    def test_an_occupied_second_locus_is_never_taken(
            self, tmp_path, hermetic_pipeline):
        # This is the property that keeps the pass additive: a gene may be
        # placed again only where no emitted model reaches.
        workspace = _build_workspace(tmp_path / "occupied",
                                     occupy_second_locus=True)
        off = _genes(_run(workspace, "--no-rescue-second-locus"))
        workspace = _build_workspace(tmp_path / "occupied_on",
                                     occupy_second_locus=True)
        on = _genes(_run(workspace))
        assert set(on) == set(off)
        assert [(g[1], g[2]) for g in on if g[0].startswith("gene2")] \
            == [(601, 699)]


def test_a_hit_on_an_occupied_locus_is_refused_before_its_cds_query(monkeypatch):
    """The sub-pass queried every candidate's CDS rows before checking whether
    its locus was already taken -- most of its 40-55 s per whole genome. The
    placement loop only ever adds intervals, so a hit overlapping an emitted
    model when the sub-pass starts is refused either way; it is now refused
    before the query."""
    from types import SimpleNamespace
    from intervaltree import Interval, IntervalTree
    from lifton import miniprot_rescue

    queried = []

    class FeatureDB:
        def children(self, feature, featuretype=None):
            queried.append(feature.attributes["ID"][0])
            return iter(())

    hit = SimpleNamespace(seqid="chr1", start=100, end=400, strand="+", id="MP1",
                          attributes={"ID": ["MP1"]})
    tree = {"chr1": IntervalTree([Interval(90, 410, "gene-emitted")])}
    args = SimpleNamespace(overlap=0.1, miniprot_rescue_len=(0.5, 2.0),
                           rescue_second_locus=True, rescue_second_locus_max=1)
    added = miniprot_rescue._second_locus_subpass(
        [hit], 0.5, FeatureDB(), None, tree, None, {"tx1": "M"}, {"tx1": "A"},
        {}, {"MP1": "tx1"}, {"gene1": 300}, {"tx1": 1}, {"tx1": "gene1"},
        {"gene1"}, None, args)
    assert added == 0
    assert queried == []
