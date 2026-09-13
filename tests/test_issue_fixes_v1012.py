"""Regressions for two reported failures fixed in v1.0.12.

GH #37: a flat annotation -- a prokaryotic bakta GFF, a miniprot GFF -- has
top-level ``CDS`` rows and no ``gene``, so the gene-like auto-detection selected
nothing and the run died several steps later inside vendored Liftoff with a bare
"Use -f to provide a list of other feature types to lift over".

GH #14: every run writing beside the same output directory shared one
``lifton_output/``, so concurrent jobs overwrote each other's intermediate
files.
"""
from __future__ import annotations

import os
import textwrap

import pytest

from lifton import annotation, lifton, lifton_utils
from lifton.exceptions import LiftOnInputError
from tests.test_integration_pipeline import hermetic_pipeline  # noqa: F401


ORF = "ATG" + "GCT" * 31 + "TAA"


def _wrap(sequence):
    return "\n".join(textwrap.wrap(sequence, 60)) + "\n"


def _chromosome(placements, length=800):
    chromosome = ["A"] * length
    for start, sequence in placements:
        chromosome[start - 1:start - 1 + len(sequence)] = sequence
    return "".join(chromosome)


def _genomes(work):
    work.mkdir(parents=True, exist_ok=True)
    body = ">chr1\n" + _wrap(_chromosome([(101, ORF)]))
    (work / "ref.fa").write_text(body)
    (work / "tgt.fa").write_text(body)


FLAT_GFF = (
    "##gff-version 3\n"
    "chr1\tbakta\tregion\t1\t800\t.\t+\t.\tID=region1\n"
    "chr1\tbakta\tCDS\t101\t199\t.\t+\t0\tID=cds1;product=hypothetical\n"
)


class TestFlatAnnotation:
    """GH #37 -- a top-level-CDS annotation selects those CDS rows."""

    def _ref_db(self, tmp_path, text):
        _genomes(tmp_path)
        path = tmp_path / "ref.gff3"
        path.write_text(text)
        return annotation.Annotation(str(path), False, False, force=True)

    def test_top_level_cds_is_detected(self, tmp_path):
        ref_db = self._ref_db(tmp_path, FLAT_GFF)
        assert lifton_utils.get_gene_like_feature_types(ref_db) == ["CDS"]

    @pytest.mark.parametrize("rows", [
        # landmarks
        "chr1\tt\tregion\t1\t800\t.\t+\t.\tID=r1\n"
        "chr1\tt\tchromosome\t1\t800\t.\t+\t.\tID=c1\n",
        # childless regulatory features are not lift targets either
        "chr1\tt\tenhancer\t10\t20\t.\t+\t.\tID=e1\n"
        "chr1\tt\tregion\t1\t800\t.\t+\t.\tID=r1\n",
    ])
    def test_nothing_liftable_still_falls_back_to_gene(self, tmp_path, rows):
        ref_db = self._ref_db(tmp_path, "##gff-version 3\n" + rows)
        assert lifton_utils.get_gene_like_feature_types(ref_db) == ["gene"]

    def test_a_gene_bearing_annotation_is_unchanged(self, tmp_path):
        ref_db = self._ref_db(
            tmp_path,
            "##gff-version 3\n"
            "chr1\tt\tregion\t1\t800\t.\t+\t.\tID=r1\n"
            "chr1\tt\tgene\t101\t199\t.\t+\t.\tID=g1\n"
            "chr1\tt\tmRNA\t101\t199\t.\t+\t.\tID=t1;Parent=g1\n"
            "chr1\tt\texon\t101\t199\t.\t+\t.\tID=e1;Parent=t1\n"
            "chr1\tt\tCDS\t101\t199\t.\t+\t0\tID=c1;Parent=t1\n")
        assert lifton_utils.get_gene_like_feature_types(ref_db) == ["gene"]

    def test_empty_selection_names_the_types_present(self, tmp_path,
                                                     hermetic_pipeline):
        _genomes(tmp_path)
        (tmp_path / "ref.gff3").write_text(
            "##gff-version 3\n"
            "chr1\tt\tregion\t1\t800\t.\t+\t.\tID=r1\n")
        (tmp_path / "features.txt").write_text("gene\n")
        (tmp_path / "out").mkdir()
        argv = [str(tmp_path / "tgt.fa"), str(tmp_path / "ref.fa"),
                "-g", str(tmp_path / "ref.gff3"),
                "-f", str(tmp_path / "features.txt"),
                "-o", str(tmp_path / "out" / "lifton.gff3"), "--force"]
        with pytest.raises(LiftOnInputError) as raised:
            lifton.run_all_lifton_steps(lifton.parse_args(argv))
        message = str(raised.value)
        assert "No features to lift" in message
        assert "gene" in message
        assert "region" in message          # what the annotation actually has
        assert "-f/--features" in message


class TestIntermediateDirectory:
    """GH #14 -- --intermediate-dir keeps concurrent runs independent."""

    def test_default_is_beside_the_output(self):
        assert lifton._run_outdirs("/data/run/out.gff3") == (
            "/data/run", "/data/run/lifton_output")

    def test_flag_relocates_the_artifact_directory(self):
        outdir, artifacts = lifton._run_outdirs("/data/run/out.gff3",
                                                "/scratch/job7")
        assert (outdir, artifacts) == ("/data/run", "/scratch/job7")

    def test_flag_is_resolved_and_expanded(self):
        _, artifacts = lifton._run_outdirs("out.gff3", "~/jobs/a")
        assert artifacts == os.path.join(os.path.expanduser("~"), "jobs", "a")
        _, relative = lifton._run_outdirs("out.gff3", "rel")
        assert os.path.isabs(relative)

    def test_two_runs_sharing_an_output_directory_get_separate_trees(
            self, tmp_path, hermetic_pipeline):
        _genomes(tmp_path)
        (tmp_path / "ref.gff3").write_text(
            "##gff-version 3\n"
            "chr1\ttest\tgene\t101\t199\t.\t+\t.\tID=g1;gene_biotype=protein_coding\n"
            "chr1\ttest\tmRNA\t101\t199\t.\t+\t.\tID=t1;Parent=g1\n"
            "chr1\ttest\texon\t101\t199\t.\t+\t.\tID=e1;Parent=t1\n"
            "chr1\ttest\tCDS\t101\t199\t.\t+\t0\tID=c1;Parent=t1\n")
        (tmp_path / "liftoff.gff3").write_text(
            "##gff-version 3\n"
            "chr1\tLiftoff\tgene\t101\t199\t.\t+\t.\tID=g1;gene_biotype=protein_coding\n"
            "chr1\tLiftoff\tmRNA\t101\t199\t.\t+\t.\tID=t1;Parent=g1\n"
            "chr1\tLiftoff\texon\t101\t199\t.\t+\t.\tID=e1;Parent=t1\n"
            "chr1\tLiftoff\tCDS\t101\t199\t.\t+\t0\tID=c1;Parent=t1\n")
        (tmp_path / "miniprot.gff3").write_text(
            "##gff-version 3\n"
            "chr1\tminiprot\tmRNA\t101\t199\t.\t+\t.\tID=MP1;Rank=1;Identity=1.0000;Target=t1 1 33\n"
            "chr1\tminiprot\tCDS\t101\t199\t.\t+\t0\tID=MP1.c1;Parent=MP1\n")
        shared = tmp_path / "shared"
        shared.mkdir()
        for name in ("a", "b"):
            argv = [str(tmp_path / "tgt.fa"), str(tmp_path / "ref.fa"),
                    "-g", str(tmp_path / "ref.gff3"),
                    "-L", str(tmp_path / "liftoff.gff3"),
                    "-M", str(tmp_path / "miniprot.gff3"),
                    "-o", str(shared / f"{name}.gff3"),
                    "--intermediate-dir", str(tmp_path / f"work_{name}"),
                    "-ad", "RefSeq", "--force"]
            lifton.run_all_lifton_steps(lifton.parse_args(argv))
        for name in ("a", "b"):
            work = tmp_path / f"work_{name}"
            assert (work / "intermediate_files").is_dir()
            assert (work / "stats").is_dir()
            assert (work / "score.txt").is_file()
        assert not (shared / "lifton_output").exists()
        assert (shared / "a.gff3").is_file() and (shared / "b.gff3").is_file()


def test_an_input_error_exits_cleanly_rather_than_tracebacking(tmp_path, capsys):
    """An unreadable input is the user's to fix: the message says what is wrong,
    and a traceback only buries it."""
    _genomes(tmp_path)
    (tmp_path / "ref.gff3").write_text(
        "##gff-version 3\nchr1\tt\tregion\t1\t800\t.\t+\t.\tID=r1\n")
    (tmp_path / "features.txt").write_text("gene\n")
    (tmp_path / "out").mkdir()
    argv = [str(tmp_path / "tgt.fa"), str(tmp_path / "ref.fa"),
            "-g", str(tmp_path / "ref.gff3"),
            "-f", str(tmp_path / "features.txt"),
            "-o", str(tmp_path / "out" / "lifton.gff3"), "--force"]
    with pytest.raises(SystemExit) as raised:
        lifton.main(argv)
    assert raised.value.code == 2
    assert "No features to lift" in capsys.readouterr().err


def test_a_meta_type_bearing_a_hierarchy_is_still_detected(tmp_path):
    """The meta-type filter applies to the fallback only. A type that genuinely
    bears a hierarchy must keep being lifted, or annotations the old code
    handled would silently lose features."""
    _genomes(tmp_path)
    path = tmp_path / "ref.gff3"
    path.write_text(
        "##gff-version 3\n"
        "chr1\tt\tregion\t1\t800\t.\t+\t.\tID=r1\n"
        "chr1\tt\tmatch\t101\t199\t.\t+\t.\tID=m1\n"
        "chr1\tt\tmatch_part\t101\t199\t.\t+\t.\tID=mp1;Parent=m1\n")
    ref_db = annotation.Annotation(str(path), False, False, force=True)
    assert lifton_utils.get_gene_like_feature_types(ref_db) == ["match"]
