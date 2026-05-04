"""End-to-end "Golden Path" integration test for run_all_lifton_steps.

The vendored Liftoff and the miniprot subprocess are bypassed by:
  (a) supplying pre-built `args.liftoff` / `args.miniprot` GFF files so
      exec_liftoff / exec_miniprot short-circuit on the os.path.exists
      check, and
  (b) monkey-patching lifton_utils.check_miniprot_installed so the
      miniprot binary need not be on PATH.

The aim is structural validation of the pipeline, not bit-exact GFF
output. We assert that:
  - intermediate directories are created,
  - lifton.gff3 exists, is non-empty, and contains gene + mRNA + exon + CDS
    feature lines all marked source=LiftOn,
  - the score / unmapped / extra-copy side-files are written.
"""

from __future__ import annotations

import os
import textwrap
from pathlib import Path

import pytest


# ---------------------------------------------------------------------------
# Fixtures specific to the integration test
# ---------------------------------------------------------------------------

def _wrap(seq: str) -> str:
    return "\n".join(textwrap.wrap(seq, 60)) + "\n"


@pytest.fixture
def integration_workspace(tmp_path):
    """Build a self-contained workspace with reference + target FASTAs,
    a reference GFF, and pre-baked liftoff / miniprot GFF outputs."""
    work = tmp_path / "work"
    work.mkdir()

    # Reference and target genomes — identical for this test so that the
    # transferred coordinates match exactly. (Liftoff would normally find
    # a near-identical mapping; we just hand-author the mapping here.)
    chrom = ["A"] * 600
    exon1 = "ATG" + "GCT" * 32          # 99 nt
    exon2 = ("GCT" * 32) + "TAA"        # 99 nt
    for i, ch in enumerate(exon1):
        chrom[100 + i] = ch
    for i, ch in enumerate(exon2):
        chrom[300 + i] = ch
    seq = "".join(chrom)

    ref_fa = work / "ref.fa"
    ref_fa.write_text(">chr1\n" + _wrap(seq))
    tgt_fa = work / "tgt.fa"
    tgt_fa.write_text(">chr1\n" + _wrap(seq))

    ref_gff = work / "ref.gff3"
    ref_gff.write_text(
        "##gff-version 3\n"
        "chr1\ttest\tgene\t101\t399\t.\t+\t.\tID=gene1;gene_biotype=protein_coding\n"
        "chr1\ttest\tmRNA\t101\t399\t.\t+\t.\tID=tx1;Parent=gene1\n"
        "chr1\ttest\texon\t101\t199\t.\t+\t.\tID=exon1;Parent=tx1\n"
        "chr1\ttest\texon\t301\t399\t.\t+\t.\tID=exon2;Parent=tx1\n"
        "chr1\ttest\tCDS\t101\t199\t.\t+\t0\tID=cds1;Parent=tx1\n"
        "chr1\ttest\tCDS\t301\t399\t.\t+\t0\tID=cds2;Parent=tx1\n"
    )

    # Pre-baked Liftoff output (identity transfer)
    liftoff_gff = work / "liftoff.gff3"
    liftoff_gff.write_text(
        "##gff-version 3\n"
        "chr1\tLiftoff\tgene\t101\t399\t.\t+\t.\tID=gene1;gene_biotype=protein_coding\n"
        "chr1\tLiftoff\tmRNA\t101\t399\t.\t+\t.\tID=tx1;Parent=gene1\n"
        "chr1\tLiftoff\texon\t101\t199\t.\t+\t.\tID=exon1;Parent=tx1\n"
        "chr1\tLiftoff\texon\t301\t399\t.\t+\t.\tID=exon2;Parent=tx1\n"
        "chr1\tLiftoff\tCDS\t101\t199\t.\t+\t0\tID=cds1;Parent=tx1\n"
        "chr1\tLiftoff\tCDS\t301\t399\t.\t+\t0\tID=cds2;Parent=tx1\n"
    )

    # Pre-baked Miniprot output — one mRNA whose Target attribute points
    # back at the reference transcript id, so miniprot_id_mapping links
    # it. Coordinates overlap the gene, single CDS to keep it simple.
    miniprot_gff = work / "miniprot.gff3"
    miniprot_gff.write_text(
        "##gff-version 3\n"
        "chr1\tminiprot\tmRNA\t101\t399\t.\t+\t.\tID=MP1;Target=tx1 1 66\n"
        "chr1\tminiprot\tCDS\t101\t199\t.\t+\t0\tID=MP1.cds1;Parent=MP1\n"
        "chr1\tminiprot\tCDS\t301\t399\t.\t+\t0\tID=MP1.cds2;Parent=MP1\n"
    )

    out_dir = work / "out"
    out_dir.mkdir()

    return {
        "ref_fa": ref_fa, "tgt_fa": tgt_fa, "ref_gff": ref_gff,
        "liftoff": liftoff_gff, "miniprot": miniprot_gff,
        "out": out_dir, "work": work,
    }


@pytest.fixture
def hermetic_pipeline(monkeypatch):
    """Patch the miniprot installation guard so the test does not require
    the miniprot binary on PATH. Also patch the run_liftoff fall-through
    so that even if args.liftoff was missing, we wouldn't spawn anything."""
    from lifton import lifton_utils, run_liftoff, run_miniprot
    monkeypatch.setattr(lifton_utils, "check_miniprot_installed",
                        lambda: None)
    monkeypatch.setattr(run_miniprot, "check_miniprot_installed",
                        lambda: True)

    def _fail(*args, **kwargs):
        raise AssertionError(
            "External liftoff/miniprot must NOT be invoked in tests"
        )

    monkeypatch.setattr(run_liftoff, "run_liftoff", _fail)
    monkeypatch.setattr(run_miniprot, "run_miniprot", _fail)


# ---------------------------------------------------------------------------
# The golden-path integration test
# ---------------------------------------------------------------------------

def test_run_all_lifton_steps_golden_path(integration_workspace,
                                          hermetic_pipeline, monkeypatch):
    from lifton import lifton as lifton_main

    out_gff = integration_workspace["out"] / "lifton.gff3"
    argv = [
        str(integration_workspace["tgt_fa"]),
        str(integration_workspace["ref_fa"]),
        "-g", str(integration_workspace["ref_gff"]),
        "-L", str(integration_workspace["liftoff"]),
        "-M", str(integration_workspace["miniprot"]),
        "-o", str(out_gff),
        "-ad", "RefSeq",
        "--force",
    ]
    args = lifton_main.parse_args(argv)
    # parse_args requires args.proteins / args.transcripts to be None so
    # extract_features runs from the gffutils DB. Default already None.
    assert args.proteins is None and args.transcripts is None

    lifton_main.run_all_lifton_steps(args)

    # ---- structural assertions on the final GFF -----------------------
    assert out_gff.exists() and out_gff.stat().st_size > 0
    body = out_gff.read_text()
    feat_types = [line.split("\t")[2] for line in body.splitlines()
                  if line.strip() and not line.startswith("#")]
    assert "gene" in feat_types
    assert "mRNA" in feat_types
    assert feat_types.count("exon") >= 2
    assert feat_types.count("CDS") >= 2
    # All output rows must be sourced from LiftOn
    sources = {line.split("\t")[1] for line in body.splitlines()
               if line.strip() and not line.startswith("#")}
    assert sources == {"LiftOn"}

    # ---- side-files ---------------------------------------------------
    out_dir = integration_workspace["out"]
    assert (out_dir / "lifton_output" / "score.txt").exists()
    assert (out_dir / "lifton_output" / "stats" /
            "unmapped_features.txt").exists()
    assert (out_dir / "lifton_output" / "stats" /
            "extra_copy_features.txt").exists()
    assert (out_dir / "lifton_output" / "stats" /
            "mapped_feature.txt").exists()


def test_pipeline_does_not_call_external_tools(integration_workspace,
                                               hermetic_pipeline):
    """If the test above passes, the _fail patches inside
    hermetic_pipeline guarantee no subprocess invocation. This second
    test pins the contract explicitly: providing -L and -M must skip
    both external runners."""
    from lifton import lifton as lifton_main

    out_gff = integration_workspace["out"] / "lifton.gff3"
    argv = [
        str(integration_workspace["tgt_fa"]),
        str(integration_workspace["ref_fa"]),
        "-g", str(integration_workspace["ref_gff"]),
        "-L", str(integration_workspace["liftoff"]),
        "-M", str(integration_workspace["miniprot"]),
        "-o", str(out_gff),
        "-ad", "RefSeq",
        "--force",
    ]
    args = lifton_main.parse_args(argv)
    # If exec_liftoff / exec_miniprot tried to fall through to the real
    # runners, hermetic_pipeline's _fail patches would raise. Reaching
    # this assertion proves the short-circuit fired.
    lifton_main.run_all_lifton_steps(args)
    assert out_gff.exists()
