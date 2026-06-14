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
def merge_firing_workspace(tmp_path):
    """A workspace where the Liftoff/miniprot **merge actually fires**.

    ``integration_workspace`` uses identical ref/tgt genomes with a perfect
    ORF, so ``liftoff_aln.identity == 1`` and the protein-maximization merge
    branch in ``run_liftoff.process_liftoff_with_protein`` is never entered —
    any test driven off it is vacuous w.r.t. the merge.

    Here the **target** carries a lesion in the lifted CDS region (101-199):
    three ``GCT`` (Ala) codons are mutated to ``GAT`` (Asp), so Liftoff's
    lifted protein is imperfect (identity < 1) and the merge branch is taken.
    A pre-baked miniprot annotation supplies a CDS over the **correct** ORF
    copy the target also carries at 301-399; the miniprot mRNA *span*
    (101-399) overlaps the Liftoff locus so the miniprot-overlap gate
    (``lifton_utils.LiftOn_miniprot_alignment`` Check 1) passes and the chunk
    is admitted, so ``has_valid_miniprot`` is True and the chaining algorithm
    runs. The emitted transcript carries ``status=LiftOn_chaining_algorithm``.

    This is the basis for ``TestMergePromotion`` — it exercises the promoted
    best-of-outcome default end-to-end (the full per-candidate compare branch,
    since the merge protein identity here is < 1.0)."""
    work = tmp_path / "work"
    work.mkdir()

    correct_cds = "ATG" + "GCT" * 31 + "TAA"                       # M + A*31 + stop
    mutated_cds = "ATG" + "GCT" * 9 + "GAT" * 3 + "GCT" * 19 + "TAA"  # 3x A->D
    assert len(correct_cds) == 99 and len(mutated_cds) == 99

    def build_chrom(region_101, region_301):
        chrom = ["A"] * 600
        for i, ch in enumerate(region_101):
            chrom[100 + i] = ch
        for i, ch in enumerate(region_301):
            chrom[300 + i] = ch
        return "".join(chrom)

    ref_seq = build_chrom(correct_cds, "A" * 99)            # reference: correct ORF at 101-199
    tgt_seq = build_chrom(mutated_cds, correct_cds)         # target: lesion at 101-199, correct at 301-399

    ref_fa = work / "ref.fa"
    ref_fa.write_text(">chr1\n" + _wrap(ref_seq))
    tgt_fa = work / "tgt.fa"
    tgt_fa.write_text(">chr1\n" + _wrap(tgt_seq))

    ref_gff = work / "ref.gff3"
    ref_gff.write_text(
        "##gff-version 3\n"
        "chr1\ttest\tgene\t101\t199\t.\t+\t.\tID=gene1;gene_biotype=protein_coding\n"
        "chr1\ttest\tmRNA\t101\t199\t.\t+\t.\tID=tx1;Parent=gene1\n"
        "chr1\ttest\texon\t101\t199\t.\t+\t.\tID=exon1;Parent=tx1\n"
        "chr1\ttest\tCDS\t101\t199\t.\t+\t0\tID=cds1;Parent=tx1\n"
    )

    # Liftoff: CDS at 101-199 (the lesioned region) -> imperfect protein.
    liftoff_gff = work / "liftoff.gff3"
    liftoff_gff.write_text(
        "##gff-version 3\n"
        "chr1\tLiftoff\tgene\t101\t199\t.\t+\t.\tID=gene1;gene_biotype=protein_coding\n"
        "chr1\tLiftoff\tmRNA\t101\t199\t.\t+\t.\tID=tx1;Parent=gene1\n"
        "chr1\tLiftoff\texon\t101\t199\t.\t+\t.\tID=exon1;Parent=tx1\n"
        "chr1\tLiftoff\tCDS\t101\t199\t.\t+\t0\tID=cds1;Parent=tx1\n"
    )

    # Miniprot: mRNA SPAN 101-399 overlaps the Liftoff locus (101-199) so the
    # overlap gate passes; its CDS sits at 301-399 where the correct ORF lives,
    # so the chained chunk wins on protein identity and the merge fires.
    miniprot_gff = work / "miniprot.gff3"
    miniprot_gff.write_text(
        "##gff-version 3\n"
        "chr1\tminiprot\tmRNA\t101\t399\t.\t+\t.\tID=MP1;Target=tx1 1 33\n"
        "chr1\tminiprot\tCDS\t301\t399\t.\t+\t0\tID=MP1.cds1;Parent=MP1\n"
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


# ---------------------------------------------------------------------------
# Iteration 5: --lift-gene-like (lift gene-like features beyond `gene`)
# ---------------------------------------------------------------------------

@pytest.fixture
def gene_like_workspace(tmp_path):
    """A protein-coding gene PLUS a gene-like pseudogene (has an exon child) at
    401-499. The pseudogene is present in both the ref and the pre-baked Liftoff
    GFF, so it is lifted ONLY when --lift-gene-like expands the feature list
    beyond `gene`. (Kept separate from `integration_workspace` so the byte-
    identity matrix fixture is untouched.)"""
    work = tmp_path / "work"
    work.mkdir()
    chrom = ["A"] * 600
    exon1 = "ATG" + "GCT" * 32          # 99 nt coding
    exon2 = ("GCT" * 32) + "TAA"        # 99 nt coding
    for i, ch in enumerate(exon1):
        chrom[100 + i] = ch
    for i, ch in enumerate(exon2):
        chrom[300 + i] = ch
    for i, ch in enumerate("GCT" * 33):  # 99 nt pseudogene body (no ATG/stop needed)
        chrom[400 + i] = ch
    seq = "".join(chrom)

    ref_fa = work / "ref.fa"
    ref_fa.write_text(">chr1\n" + _wrap(seq))
    tgt_fa = work / "tgt.fa"
    tgt_fa.write_text(">chr1\n" + _wrap(seq))

    gene_block = (
        "chr1\t{src}\tgene\t101\t399\t.\t+\t.\tID=gene1;gene_biotype=protein_coding\n"
        "chr1\t{src}\tmRNA\t101\t399\t.\t+\t.\tID=tx1;Parent=gene1\n"
        "chr1\t{src}\texon\t101\t199\t.\t+\t.\tID=exon1;Parent=tx1\n"
        "chr1\t{src}\texon\t301\t399\t.\t+\t.\tID=exon2;Parent=tx1\n"
        "chr1\t{src}\tCDS\t101\t199\t.\t+\t0\tID=cds1;Parent=tx1\n"
        "chr1\t{src}\tCDS\t301\t399\t.\t+\t0\tID=cds2;Parent=tx1\n"
        "chr1\t{src}\tpseudogene\t401\t499\t.\t+\t.\tID=pg1;gene_biotype=pseudogene\n"
        "chr1\t{src}\texon\t401\t499\t.\t+\t.\tID=pgexon1;Parent=pg1\n"
    )
    ref_gff = work / "ref.gff3"
    ref_gff.write_text("##gff-version 3\n" + gene_block.format(src="test"))
    liftoff_gff = work / "liftoff.gff3"
    liftoff_gff.write_text("##gff-version 3\n" + gene_block.format(src="Liftoff"))
    miniprot_gff = work / "miniprot.gff3"
    miniprot_gff.write_text(
        "##gff-version 3\n"
        "chr1\tminiprot\tmRNA\t101\t399\t.\t+\t.\tID=MP1;Target=tx1 1 66\n"
        "chr1\tminiprot\tCDS\t101\t199\t.\t+\t0\tID=MP1.cds1;Parent=MP1\n"
        "chr1\tminiprot\tCDS\t301\t399\t.\t+\t0\tID=MP1.cds2;Parent=MP1\n"
    )
    out_dir = work / "out"
    out_dir.mkdir()
    return {"ref_fa": ref_fa, "tgt_fa": tgt_fa, "ref_gff": ref_gff,
            "liftoff": liftoff_gff, "miniprot": miniprot_gff, "out": out_dir}


def _run_gene_like(ws, lift_gene_like):
    from lifton import lifton as lifton_main
    out_gff = ws["out"] / "lifton.gff3"
    argv = [str(ws["tgt_fa"]), str(ws["ref_fa"]),
            "-g", str(ws["ref_gff"]), "-L", str(ws["liftoff"]),
            "-M", str(ws["miniprot"]), "-o", str(out_gff),
            "-ad", "RefSeq", "--force"]
    if lift_gene_like:
        argv.append("--lift-gene-like")
    lifton_main.run_all_lifton_steps(lifton_main.parse_args(argv))
    body = out_gff.read_text()
    return [ln.split("\t")[2] for ln in body.splitlines()
            if ln.strip() and not ln.startswith("#")]


class TestLiftGeneLike:
    def test_flag_on_emits_pseudogene(self, gene_like_workspace, hermetic_pipeline):
        types = _run_gene_like(gene_like_workspace, lift_gene_like=True)
        assert "pseudogene" in types          # gene-like pseudogene captured
        assert "gene" in types                # the coding gene still lifted

    def test_default_omits_pseudogene(self, gene_like_workspace, hermetic_pipeline):
        types = _run_gene_like(gene_like_workspace, lift_gene_like=False)
        assert "pseudogene" not in types      # default lifts only `gene`
        assert "gene" in types
