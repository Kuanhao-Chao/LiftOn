"""Protein-coverage rescue sub-pass (v1.0.12 A1).

The miniprot-only rescue's length band divides a hit's genomic span by the
reference gene's CDS span. Both include introns, so a gene lifted into a genome
with compact introns fails the band even when miniprot aligns the whole
protein. Sub-pass B reconsiders exactly those candidates, gated on protein
coverage. These tests pin its helpers, its flag resolution, and its end-to-end
behaviour on a fixture whose target intron is 150 bp where the reference
intron is 1,300 bp: span ratio 0.18, coverage 1.0.
"""
from __future__ import annotations

import textwrap
import types

import pytest

from lifton import lifton, miniprot_rescue
from tests.test_integration_pipeline import hermetic_pipeline  # noqa: F401


class _Feature:
    def __init__(self, attributes, seqid="chr1", start=1, end=100):
        self.attributes = attributes
        self.seqid, self.start, self.end = seqid, start, end


PROTEIN = "M" + "A" * 31 + "*"   # 32 residues plus the terminal stop


# --------------------------------------------------------------------------- #
# Coverage
# --------------------------------------------------------------------------- #
class TestProteinCoverage:
    def _coverage(self, target, proteins=None):
        attributes = {"ID": ["MP"]}
        if target is not None:
            attributes["Target"] = [target]
        return miniprot_rescue.miniprot_protein_coverage(
            _Feature(attributes),
            {"tx": PROTEIN} if proteins is None else proteins, "tx")

    def test_full_length_hit_ignores_the_stop(self):
        assert self._coverage("tx 1 32") == 1.0

    def test_partial_hit(self):
        assert self._coverage("tx 1 8") == pytest.approx(8 / 32)
        assert self._coverage("tx 17 32") == pytest.approx(16 / 32)

    def test_hit_through_the_stop_is_clamped(self):
        assert self._coverage("tx 1 33") == 1.0

    def test_protein_without_a_stop_uses_its_full_length(self):
        assert self._coverage("tx 1 16", {"tx": "M" * 32}) == 0.5

    @pytest.mark.parametrize("target", [None, "", "tx", "tx 1", "tx a 5",
                                        "tx 9 3"])
    def test_unusable_target_is_none(self, target):
        assert self._coverage(target) is None

    def test_missing_or_empty_protein_is_none(self):
        assert self._coverage("tx 1 5", {}) is None
        assert self._coverage("tx 1 5", {"tx": "*"}) is None


# --------------------------------------------------------------------------- #
# Candidate order
# --------------------------------------------------------------------------- #
class TestCandidateOrder:
    def _key(self, identifier, rank, identity, coverage, start):
        attributes = {"ID": [identifier]}
        if rank is not None:
            attributes["Rank"] = [str(rank)]
        if identity is not None:
            attributes["Identity"] = [str(identity)]
        return miniprot_rescue._candidate_quality_key(
            _Feature(attributes, start=start, end=start + 10), coverage)

    def test_best_hit_first_regardless_of_position(self):
        keys = {
            "late_rank1": self._key("late_rank1", 1, 0.60, 0.9, 900),
            "early_rank2": self._key("early_rank2", 2, 0.99, 1.0, 10),
            "rank1_higher_identity": self._key("rank1_higher_identity", 1, 0.90,
                                               0.8, 500),
            "rank1_same_identity_more_coverage": self._key(
                "rank1_same_identity_more_coverage", 1, 0.60, 1.0, 950),
            "no_rank": self._key("no_rank", None, None, 1.0, 1),
        }
        ordered = sorted(keys, key=keys.get)
        assert ordered == ["rank1_higher_identity",
                           "rank1_same_identity_more_coverage", "late_rank1",
                           "early_rank2", "no_rank"]

    def test_order_is_total(self):
        first = self._key("a", 1, 0.5, 1.0, 10)
        second = self._key("b", 1, 0.5, 1.0, 10)
        assert first != second and sorted([second, first]) == [first, second]


# --------------------------------------------------------------------------- #
# Flag and environment resolution
# --------------------------------------------------------------------------- #
class TestResolution:
    def _parse(self, *extra):
        return lifton.parse_args(["t.fa", "r.fa", "-g", "r.gff3", "-o", "o.gff3",
                                  *extra])

    @pytest.fixture(autouse=True)
    def _clean_env(self, monkeypatch):
        for key in ("LIFTON_RESCUE_COVERAGE_GATE", "LIFTON_RESCUE_COVERAGE_MIN"):
            monkeypatch.delenv(key, raising=False)

    def test_unset_flag_defers_to_the_module_default(self, monkeypatch):
        args = lifton.resolve_miniprot_rescue_args(self._parse())
        assert args.coverage_rescue_gate is None
        assert args.miniprot_rescue_coverage_min == 0.8
        for default in (False, True):
            monkeypatch.setattr(miniprot_rescue, "COVERAGE_GATE_DEFAULT", default)
            assert miniprot_rescue._coverage_gate_on(args) is default

    def test_explicit_flags(self):
        on = lifton.resolve_miniprot_rescue_args(
            self._parse("--coverage-rescue-gate"))
        off = lifton.resolve_miniprot_rescue_args(
            self._parse("--no-coverage-rescue-gate"))
        assert miniprot_rescue._coverage_gate_on(on) is True
        assert miniprot_rescue._coverage_gate_on(off) is False

    @pytest.mark.parametrize("value, expected", [("1", True), ("yes", True),
                                                 ("0", False), ("false", False)])
    def test_environment_wins_over_the_flag(self, monkeypatch, value, expected):
        monkeypatch.setenv("LIFTON_RESCUE_COVERAGE_GATE", value)
        flag = "--no-coverage-rescue-gate" if expected else "--coverage-rescue-gate"
        args = lifton.resolve_miniprot_rescue_args(self._parse(flag))
        assert miniprot_rescue._coverage_gate_on(args) is expected

    def test_coverage_threshold_from_environment(self, monkeypatch):
        monkeypatch.setenv("LIFTON_RESCUE_COVERAGE_MIN", "0.95")
        args = lifton.resolve_miniprot_rescue_args(self._parse())
        assert miniprot_rescue._coverage_min(args) == 0.95
        monkeypatch.setenv("LIFTON_RESCUE_COVERAGE_MIN", "garbage")
        args = lifton.resolve_miniprot_rescue_args(self._parse())
        assert miniprot_rescue._coverage_min(args) == 0.8
        assert miniprot_rescue._coverage_min(types.SimpleNamespace()) == 0.8


# --------------------------------------------------------------------------- #
# End to end (pre-baked Liftoff/miniprot; external tools never run)
# --------------------------------------------------------------------------- #
GENE1_EXON1 = "ATG" + "GCT" * 32
GENE1_EXON2 = "GCT" * 32 + "TAA"
GENE2_EXON1 = "ATG" + "GCT" * 16          # 51 nt
GENE2_EXON2 = "GCT" * 15 + "TAA"          # 48 nt -> M + 31 A + stop


def _wrap(sequence):
    return "\n".join(textwrap.wrap(sequence, 60)) + "\n"


def _chromosome(placements, length=2200):
    chromosome = ["A"] * length
    for start, sequence in placements:
        chromosome[start - 1:start - 1 + len(sequence)] = sequence
    return "".join(chromosome)


def _build_workspace(work, *, target_cds):
    """gene1 is lifted by Liftoff. gene2 has a 1,300 bp intron in the reference
    (CDS span 601-1998) and is compact in the target; miniprot places it at
    ``target_cds`` (two segments) with full coverage. Liftoff misses gene2."""
    work.mkdir(parents=True, exist_ok=True)
    (work / "ref.fa").write_text(">chr1\n" + _wrap(_chromosome([
        (101, GENE1_EXON1), (301, GENE1_EXON2),
        (601, GENE2_EXON1), (1951, GENE2_EXON2)])))
    (second_start, second_end) = target_cds[1]
    second = (GENE2_EXON2 if second_end - second_start + 1 == len(GENE2_EXON2)
              else "GCT" * ((second_end - second_start + 1 - 3) // 3) + "TAA")
    (work / "tgt.fa").write_text(">chr1\n" + _wrap(_chromosome([
        (101, GENE1_EXON1), (301, GENE1_EXON2),
        (target_cds[0][0], GENE2_EXON1), (second_start, second)])))
    gene1 = (
        "chr1\t{src}\tgene\t101\t399\t.\t+\t.\tID=gene1;gene_biotype=protein_coding\n"
        "chr1\t{src}\tmRNA\t101\t399\t.\t+\t.\tID=tx1;Parent=gene1\n"
        "chr1\t{src}\texon\t101\t199\t.\t+\t.\tID=exon1;Parent=tx1\n"
        "chr1\t{src}\texon\t301\t399\t.\t+\t.\tID=exon2;Parent=tx1\n"
        "chr1\t{src}\tCDS\t101\t199\t.\t+\t0\tID=cds1;Parent=tx1\n"
        "chr1\t{src}\tCDS\t301\t399\t.\t+\t0\tID=cds2;Parent=tx1\n")
    gene2 = (
        "chr1\ttest\tgene\t601\t1998\t.\t+\t.\tID=gene2;gene_biotype=protein_coding\n"
        "chr1\ttest\tmRNA\t601\t1998\t.\t+\t.\tID=tx2;Parent=gene2\n"
        "chr1\ttest\texon\t601\t651\t.\t+\t.\tID=exon3;Parent=tx2\n"
        "chr1\ttest\texon\t1951\t1998\t.\t+\t.\tID=exon4;Parent=tx2\n"
        "chr1\ttest\tCDS\t601\t651\t.\t+\t0\tID=cds3;Parent=tx2\n"
        "chr1\ttest\tCDS\t1951\t1998\t.\t+\t0\tID=cds3;Parent=tx2\n")
    (work / "ref.gff3").write_text("##gff-version 3\n" + gene1.format(src="test")
                                   + gene2)
    (work / "liftoff.gff3").write_text("##gff-version 3\n"
                                       + gene1.format(src="Liftoff"))
    (first_start, first_end) = target_cds[0]
    (work / "miniprot.gff3").write_text(
        "##gff-version 3\n"
        "chr1\tminiprot\tmRNA\t101\t399\t.\t+\t.\tID=MP1;Rank=1;Identity=1.0000;Target=tx1 1 65\n"
        "chr1\tminiprot\tCDS\t101\t199\t.\t+\t0\tID=MP1.c1;Parent=MP1\n"
        "chr1\tminiprot\tCDS\t301\t399\t.\t+\t0\tID=MP1.c2;Parent=MP1\n"
        f"chr1\tminiprot\tmRNA\t{first_start}\t{second_end}\t.\t+\t.\t"
        "ID=MP2;Rank=1;Identity=1.0000;Target=tx2 1 32\n"
        f"chr1\tminiprot\tCDS\t{first_start}\t{first_end}\t.\t+\t0\tID=MP2.c1;Parent=MP2\n"
        f"chr1\tminiprot\tCDS\t{second_start}\t{second_end}\t.\t+\t0\tID=MP2.c2;Parent=MP2\n")
    (work / "out").mkdir()
    return work


def _run(work, *flags):
    output = work / "out" / "lifton.gff3"
    argv = [str(work / "tgt.fa"), str(work / "ref.fa"), "-g", str(work / "ref.gff3"),
            "-L", str(work / "liftoff.gff3"), "-M", str(work / "miniprot.gff3"),
            "-o", str(output), "-ad", "RefSeq", "--force", *flags]
    lifton.run_all_lifton_steps(lifton.parse_args(argv))
    return output.read_text()


def _mrna_rows(body):
    return [row for row in body.splitlines()
            if not row.startswith("#") and row.split("\t")[2:3] == ["mRNA"]]


COMPACT = ((601, 651), (801, 848))       # target span 248 / ref CDS span 1398


@pytest.fixture
def _clean_rescue_env(monkeypatch):
    for key in ("LIFTON_RESCUE_COVERAGE_GATE", "LIFTON_RESCUE_COVERAGE_MIN",
                "LIFTON_MINIPROT_RESCUE", "LIFTON_MINIPROT_RESCUE_MIN_ID",
                "LIFTON_MINIPROT_RESCUE_LEN", "LIFTON_RESCUE_ADAPTIVE_FLOOR"):
        monkeypatch.delenv(key, raising=False)


@pytest.mark.usefixtures("_clean_rescue_env")
class TestCoverageGatePipeline:
    def test_off_leaves_the_compact_gene_unlifted(self, tmp_path, hermetic_pipeline):
        body = _run(_build_workspace(tmp_path / "w", target_cds=COMPACT),
                    "--no-coverage-rescue-gate")
        assert "ID=tx2" not in body
        assert "rescue_gate=protein_coverage" not in body

    def test_on_rescues_it_and_only_appends(self, tmp_path, hermetic_pipeline):
        off = _run(_build_workspace(tmp_path / "off", target_cds=COMPACT),
                   "--no-coverage-rescue-gate")
        on = _run(_build_workspace(tmp_path / "on", target_cds=COMPACT),
                  "--coverage-rescue-gate")
        rows = [row for row in _mrna_rows(on) if "ID=tx2" in row]
        assert len(rows) == 1
        attributes = rows[0].split("\t")[8]
        for expected in ("lifton_rescue=miniprot_only",
                         "rescue_gate=protein_coverage",
                         "miniprot_protein_coverage=1.000",
                         "miniprot_annotation_ratio=0.177",
                         "protein_identity=1.000"):
            assert expected in attributes
        # Sub-pass B runs last and only adds: the OFF output is a prefix.
        assert on.startswith(off) and len(on) > len(off)

    def test_environment_enables_it(self, tmp_path, hermetic_pipeline, monkeypatch):
        monkeypatch.setenv("LIFTON_RESCUE_COVERAGE_GATE", "1")
        body = _run(_build_workspace(tmp_path / "w", target_cds=COMPACT))
        assert "rescue_gate=protein_coverage" in body

    def test_partial_hit_is_not_rescued(self, tmp_path, hermetic_pipeline):
        work = _build_workspace(tmp_path / "w", target_cds=COMPACT)
        text = (work / "miniprot.gff3").read_text().replace("Target=tx2 1 32",
                                                            "Target=tx2 1 12")
        (work / "miniprot.gff3").write_text(text)
        assert "ID=tx2" not in _run(work, "--coverage-rescue-gate")

    @pytest.mark.parametrize("second_cds, rescued", [
        ((801, 890), True),     # 51 + 90 = 141 nt <= 1.5 x 96
        ((801, 950), False),    # 51 + 150 = 201 nt > 1.5 x 96
    ])
    def test_cds_length_bound(self, tmp_path, hermetic_pipeline, monkeypatch,
                              second_cds, rescued):
        # Drop the identity floor so only the length bound can decide.
        monkeypatch.setenv("LIFTON_MINIPROT_RESCUE_MIN_ID", "0")
        monkeypatch.setenv("LIFTON_RESCUE_ADAPTIVE_FLOOR", "0")
        work = _build_workspace(tmp_path / "w",
                                target_cds=((601, 651), second_cds))
        assert ("ID=tx2" in _run(work, "--coverage-rescue-gate")) is rescued
