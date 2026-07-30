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
import json
import textwrap
from pathlib import Path
from types import SimpleNamespace

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

    # Every completed run publishes an auditable, content-addressed manifest.
    manifest_path = out_dir / "lifton_output" / "run_manifest.json"
    manifest = json.loads(manifest_path.read_text())
    assert manifest["run"]["status"] == "success"
    assert manifest["inputs"]["target_genome"]["sha256"]
    assert manifest["inputs"]["reference_annotation"]["sha256"]
    assert manifest["counts"]["processed_features"] >= 1
    assert manifest["counts"]["miniprot_candidates_processed"] == 1
    assert manifest["counts"]["miniprot_genes_emitted"] == 0
    assert manifest["counts"]["step8_max_inflight_observed"] == 1
    assert manifest["counts"]["miniprot_candidates_pre_overlap_elided"] == 0
    assert manifest["counts"]["miniprot_candidates_submitted_for_analysis"] == 1
    assert manifest["counts"]["step8_child_batch_calls"] == 0
    assert manifest["counts"]["step8_child_batch_fallbacks"] == 0
    assert manifest["counts"]["step8_child_scalar_materializations"] == 0
    assert manifest["run"]["backend"]["step8_dispatch"] == "serial"
    assert manifest["run"]["backend"]["step8_overlap_filter"] == "in_process"
    assert manifest["validation"]["gff3_structural"]["passed"] is True
    assert manifest["run"]["cache"]["reference_annotation"] == "rebuilt"


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


def test_invalid_validation_preserves_previous_output(
        integration_workspace, hermetic_pipeline, monkeypatch):
    from lifton import gff3_validator
    from lifton import lifton as lifton_main
    from lifton.exceptions import LiftOnValidationError

    out_gff = integration_workspace["out"] / "lifton.gff3"
    out_gff.write_text("known-good-output\n")
    monkeypatch.setattr(
        gff3_validator, "validate_gff3_file",
        lambda **kwargs: SimpleNamespace(
            is_valid=False, errors=[object()], warnings=[],
        ),
    )
    monkeypatch.setattr(
        gff3_validator, "print_validation_report", lambda *args, **kwargs: None,
    )
    args = lifton_main.parse_args([
        str(integration_workspace["tgt_fa"]),
        str(integration_workspace["ref_fa"]),
        "-g", str(integration_workspace["ref_gff"]),
        "-L", str(integration_workspace["liftoff"]),
        "-M", str(integration_workspace["miniprot"]),
        "-o", str(out_gff),
        "-ad", "RefSeq",
        "--validate-output",
        "--force",
    ])

    with pytest.raises(LiftOnValidationError):
        lifton_main.run_all_lifton_steps(args)

    assert out_gff.read_text() == "known-good-output\n"
    partial = out_gff.with_name("lifton.partial.gff3")
    assert partial.exists() and "##gff-version 3" in partial.read_text()
    manifest = json.loads(
        (integration_workspace["out"] / "lifton_output" /
         "run_manifest.json").read_text()
    )
    assert manifest["run"]["status"] == "failed"
    assert manifest["validation"]["gff3"]["passed"] is False


def test_mandatory_structural_validation_cannot_be_bypassed(
        integration_workspace, hermetic_pipeline, monkeypatch):
    from lifton import gff3_validator
    from lifton import lifton as lifton_main
    from lifton.exceptions import LiftOnValidationError

    out_gff = integration_workspace["out"] / "lifton.gff3"
    out_gff.write_text("known-good-output\n")
    monkeypatch.setattr(
        gff3_validator,
        "validate_gff3_structure",
        lambda *args, **kwargs: SimpleNamespace(
            is_valid=False,
            errors=[object()],
            warnings=[],
            # `result.issues` is capped per check, so the pipeline reports the
            # uncapped severity total; a ValidationResult always carries it.
            severity_totals={},
            issue_totals={},
            stats={},
            file_path=str(out_gff),
            total_lines=1,
            data_lines=1,
            comment_lines=0,
        ),
    )
    monkeypatch.setattr(
        gff3_validator, "print_validation_report", lambda *args, **kwargs: None,
    )
    args = lifton_main.parse_args([
        str(integration_workspace["tgt_fa"]),
        str(integration_workspace["ref_fa"]),
        "-g", str(integration_workspace["ref_gff"]),
        "-L", str(integration_workspace["liftoff"]),
        "-M", str(integration_workspace["miniprot"]),
        "-o", str(out_gff),
        "-ad", "RefSeq",
        "--allow-partial-output",
        "--force",
    ])

    with pytest.raises(LiftOnValidationError):
        lifton_main.run_all_lifton_steps(args)

    assert out_gff.read_text() == "known-good-output\n"
    assert out_gff.with_name("lifton.partial.gff3").exists()
    manifest = json.loads(
        (integration_workspace["out"] / "lifton_output" /
         "run_manifest.json").read_text()
    )
    assert manifest["run"]["status"] == "failed"
    assert manifest["validation"]["gff3_structural"]["passed"] is False


@pytest.mark.parametrize(
    "extra_flags,publishes",
    [
        ([], True),                                              # default: publish + report
        (["--strict-completeness"], False),                      # opt back in to blocking
        (["--strict-completeness", "--allow-partial-output"], True),
    ],
)
def test_skipped_locus_publishes_unless_strict_completeness(
        integration_workspace, hermetic_pipeline, monkeypatch,
        extra_flags, publishes):
    """A skipped locus is a COMPLETENESS loss, not corruption.

    The staged GFF3 has already passed the mandatory structural validation, so it is
    well-formed; withholding a whole genome because one locus failed is not useful.
    LiftOn therefore publishes by default and reports `partial_success`, recording every
    skipped locus in the manifest. `--strict-completeness` restores the pre-v1.0.10
    behaviour of blocking, and `--allow-partial-output` overrides that in turn.
    """
    from lifton import parallel
    from lifton import lifton as lifton_main
    from lifton.exceptions import LiftOnPartialOutputError

    original_step7 = parallel.parallel_step7

    def step7_with_reported_failure(*args, **kwargs):
        processed = original_step7(*args, **kwargs)
        args[2].failure_records.append({
            "kind": "serialization_error",
            "locus_id": "synthetic-locus",
            "detail": "synthetic serialization failure",
        })
        return processed

    monkeypatch.setattr(parallel, "parallel_step7", step7_with_reported_failure)
    out_gff = integration_workspace["out"] / "lifton.gff3"
    out_gff.write_text("known-good-output\n")
    argv = [
        str(integration_workspace["tgt_fa"]),
        str(integration_workspace["ref_fa"]),
        "-g", str(integration_workspace["ref_gff"]),
        "-L", str(integration_workspace["liftoff"]),
        "-M", str(integration_workspace["miniprot"]),
        "-o", str(out_gff),
        "-ad", "RefSeq",
        "--force",
        *extra_flags,
    ]
    args = lifton_main.parse_args(argv)

    if publishes:
        lifton_main.run_all_lifton_steps(args)
        assert out_gff.read_text().startswith("##gff-version 3\n")
        expected_status = "partial_success"
    else:
        with pytest.raises(LiftOnPartialOutputError):
            lifton_main.run_all_lifton_steps(args)
        assert out_gff.read_text() == "known-good-output\n"
        assert out_gff.with_name("lifton.partial.gff3").exists()
        expected_status = "failed"

    manifest = json.loads(
        (integration_workspace["out"] / "lifton_output" /
         "run_manifest.json").read_text()
    )
    assert manifest["run"]["status"] == expected_status
    # The skipped locus is always recorded, published or not.
    assert manifest["counts"]["pipeline_failures"] >= 1


def test_stdout_mode_emits_only_gff3(
        integration_workspace, hermetic_pipeline, monkeypatch, capsys):
    from lifton import lifton as lifton_main

    monkeypatch.chdir(integration_workspace["out"])
    args = lifton_main.parse_args([
        str(integration_workspace["tgt_fa"]),
        str(integration_workspace["ref_fa"]),
        "-g", str(integration_workspace["ref_gff"]),
        "-L", str(integration_workspace["liftoff"]),
        "-M", str(integration_workspace["miniprot"]),
        "-o", "stdout",
        "-ad", "RefSeq",
        "--force",
    ])

    lifton_main.run_all_lifton_steps(args)

    stdout = capsys.readouterr().out
    assert stdout.startswith("##gff-version 3\n")
    assert "LiftOn processed" not in stdout
    assert "Total features" not in stdout
    for line in stdout.splitlines():
        if line and not line.startswith("#"):
            assert len(line.split("\t")) == 9


def test_direct_annotation_database_runs_without_text_validation(
        integration_workspace, hermetic_pipeline):
    import gffutils

    from lifton import lifton as lifton_main

    direct_db = integration_workspace["work"] / "reference.sqlite"
    built = gffutils.create_db(
        str(integration_workspace["ref_gff"]),
        dbfn=str(direct_db),
        force=True,
        keep_order=True,
        merge_strategy="create_unique",
        disable_infer_genes=True,
        disable_infer_transcripts=True,
    )
    built.conn.close()
    out_gff = integration_workspace["out"] / "direct-db.gff3"
    args = lifton_main.parse_args([
        str(integration_workspace["tgt_fa"]),
        str(integration_workspace["ref_fa"]),
        "-g", str(direct_db),
        "-L", str(integration_workspace["liftoff"]),
        "-M", str(integration_workspace["miniprot"]),
        "-o", str(out_gff),
        "-ad", "RefSeq",
    ])

    lifton_main.run_all_lifton_steps(args)

    assert "\tgene\t" in out_gff.read_text()
    manifest = json.loads(
        (integration_workspace["out"] / "lifton_output" /
         "run_manifest.json").read_text()
    )
    assert manifest["run"]["backend"][
        "reference_annotation_input"
    ] == "database"
    assert manifest["run"]["cache"][
        "reference_annotation"
    ] == "direct_database"


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


def _run_gene_like(ws, *, gene_only=False, lift_gene_like=False):
    from lifton import lifton as lifton_main
    out_gff = ws["out"] / "lifton.gff3"
    argv = [str(ws["tgt_fa"]), str(ws["ref_fa"]),
            "-g", str(ws["ref_gff"]), "-L", str(ws["liftoff"]),
            "-M", str(ws["miniprot"]), "-o", str(out_gff),
            "-ad", "RefSeq", "--force"]
    if gene_only:
        argv.append("--gene-only")
    if lift_gene_like:
        argv.append("--lift-gene-like")
    lifton_main.run_all_lifton_steps(lifton_main.parse_args(argv))
    body = out_gff.read_bytes()
    types = [ln.split("\t")[2] for ln in body.decode().splitlines()
             if ln.strip() and not ln.startswith("#")]
    return types, body


class TestLiftGeneLike:
    def test_default_emits_pseudogene(self, gene_like_workspace, hermetic_pipeline):
        # Iteration 12: the gene-like lift is now the DEFAULT.
        types, _ = _run_gene_like(gene_like_workspace)
        assert "pseudogene" in types          # gene-like pseudogene captured by default
        assert "gene" in types                # the coding gene still lifted

    def test_gene_only_omits_pseudogene(self, gene_like_workspace, hermetic_pipeline):
        # --gene-only is the opt-out: restores the pre-Iteration-12 gene-only lift.
        types, _ = _run_gene_like(gene_like_workspace, gene_only=True)
        assert "pseudogene" not in types      # opt-out lifts only `gene`
        assert "gene" in types

    def test_lift_gene_like_alias_is_noop(self, gene_like_workspace, hermetic_pipeline):
        # --lift-gene-like is now a kept no-op alias == default (byte-identical).
        _, default_bytes = _run_gene_like(gene_like_workspace)
        alias_types, alias_bytes = _run_gene_like(gene_like_workspace, lift_gene_like=True)
        assert alias_bytes == default_bytes
        assert "pseudogene" in alias_types


# ---------------------------------------------------------------------------
# Iteration 23: clean separate-pass miniprot-only rescue (--miniprot-rescue)
# ---------------------------------------------------------------------------

def _build_rescue_workspace(work, *, lift_gene2):
    """A reference with TWO coding genes. gene1 (101-399) is always lifted by the
    pre-baked Liftoff GFF. gene2 (single-CDS ORF at 601-699) is lifted ONLY when
    ``lift_gene2`` is True; otherwise the DNA lift "misses" it. The miniprot GFF
    covers gene2 with a UTR-padded mRNA (601-770 -> span/CDS-span = 170/99 = 1.72,
    OUTSIDE the default -min/-max_miniprot band [0.9,1.5] but INSIDE the rescue
    band (0.5,2.0)) and a full-coverage CDS (601-699 -> protein identity 1.0).
    So default Step 8 drops the gene2 candidate (out of band) and only the
    separate rescue pass can recover it. MP1 over gene1 is overlap-suppressed."""
    work.mkdir(parents=True, exist_ok=True)
    chrom = ["A"] * 900
    exon1 = "ATG" + "GCT" * 32          # 99 nt (gene1 exon1)
    exon2 = ("GCT" * 32) + "TAA"        # 99 nt (gene1 exon2)
    g2cds = "ATG" + "GCT" * 31 + "TAA"  # 99 nt (gene2 single-CDS ORF: M + 31A + stop)
    for i, ch in enumerate(exon1):
        chrom[100 + i] = ch
    for i, ch in enumerate(exon2):
        chrom[300 + i] = ch
    for i, ch in enumerate(g2cds):
        chrom[600 + i] = ch             # gene2 CDS at 601-699
    seq = "".join(chrom)

    ref_fa = work / "ref.fa"; ref_fa.write_text(">chr1\n" + _wrap(seq))
    tgt_fa = work / "tgt.fa"; tgt_fa.write_text(">chr1\n" + _wrap(seq))

    gene1_block = (
        "chr1\t{src}\tgene\t101\t399\t.\t+\t.\tID=gene1;gene_biotype=protein_coding\n"
        "chr1\t{src}\tmRNA\t101\t399\t.\t+\t.\tID=tx1;Parent=gene1\n"
        "chr1\t{src}\texon\t101\t199\t.\t+\t.\tID=exon1;Parent=tx1\n"
        "chr1\t{src}\texon\t301\t399\t.\t+\t.\tID=exon2;Parent=tx1\n"
        "chr1\t{src}\tCDS\t101\t199\t.\t+\t0\tID=cds1;Parent=tx1\n"
        "chr1\t{src}\tCDS\t301\t399\t.\t+\t0\tID=cds2;Parent=tx1\n"
    )
    gene2_block = (
        "chr1\t{src}\tgene\t601\t699\t.\t+\t.\tID=gene2;gene_biotype=protein_coding\n"
        "chr1\t{src}\tmRNA\t601\t699\t.\t+\t.\tID=tx2;Parent=gene2\n"
        "chr1\t{src}\texon\t601\t699\t.\t+\t.\tID=exon3;Parent=tx2\n"
        "chr1\t{src}\tCDS\t601\t699\t.\t+\t0\tID=cds3;Parent=tx2\n"
    )
    ref_gff = work / "ref.gff3"
    ref_gff.write_text("##gff-version 3\n"
                       + gene1_block.format(src="test")
                       + gene2_block.format(src="test"))

    liftoff_body = gene1_block.format(src="Liftoff")
    if lift_gene2:
        liftoff_body += gene2_block.format(src="Liftoff")
    (work / "liftoff.gff3").write_text("##gff-version 3\n" + liftoff_body)

    (work / "miniprot.gff3").write_text(
        "##gff-version 3\n"
        "chr1\tminiprot\tmRNA\t101\t399\t.\t+\t.\tID=MP1;Target=tx1 1 66\n"
        "chr1\tminiprot\tCDS\t101\t199\t.\t+\t0\tID=MP1.cds1;Parent=MP1\n"
        "chr1\tminiprot\tCDS\t301\t399\t.\t+\t0\tID=MP1.cds2;Parent=MP1\n"
        "chr1\tminiprot\tmRNA\t601\t770\t.\t+\t.\tID=MP2;Target=tx2 1 32\n"
        "chr1\tminiprot\tCDS\t601\t699\t.\t+\t0\tID=MP2.cds1;Parent=MP2\n"
    )
    out_dir = work / "out"; out_dir.mkdir()
    return {"ref_fa": ref_fa, "tgt_fa": tgt_fa, "ref_gff": ref_gff,
            "liftoff": work / "liftoff.gff3", "miniprot": work / "miniprot.gff3",
            "out": out_dir}


def _run_rescue(ws, *, rescue=False):
    from lifton import lifton as lifton_main
    out_gff = ws["out"] / "lifton.gff3"
    argv = [str(ws["tgt_fa"]), str(ws["ref_fa"]),
            "-g", str(ws["ref_gff"]), "-L", str(ws["liftoff"]),
            "-M", str(ws["miniprot"]), "-o", str(out_gff),
            "-ad", "RefSeq", "--force"]
    # Iteration-23 promotion: rescue is default-ON, so the OFF case must opt out
    # explicitly with --no-miniprot-rescue.
    argv.append("--miniprot-rescue" if rescue else "--no-miniprot-rescue")
    lifton_main.run_all_lifton_steps(lifton_main.parse_args(argv))
    return out_gff.read_text()


def _gene_ids(body):
    ids = []
    for ln in body.splitlines():
        if not ln.strip() or ln.startswith("#"):
            continue
        cols = ln.split("\t")
        if len(cols) >= 9 and cols[2] == "gene":
            for kv in cols[8].split(";"):
                if kv.startswith("ID="):
                    ids.append(kv[3:])
    return ids


class TestMiniprotRescuePass:
    def test_rescue_off_omits_missing_gene(self, tmp_path, hermetic_pipeline):
        # Default (flag OFF): the DNA lift missed gene2 and default Step 8 drops
        # the out-of-band miniprot candidate -> gene2 absent.
        ws = _build_rescue_workspace(tmp_path / "w", lift_gene2=False)
        body = _run_rescue(ws, rescue=False)
        ids = _gene_ids(body)
        assert "gene1" in ids
        assert "gene2" not in ids
        assert "lifton_rescue=miniprot_only" not in body

    def test_rescue_on_emits_missing_gene(self, tmp_path, hermetic_pipeline):
        # Flag ON: the separate post-pass recovers gene2 exactly once, tagged.
        ws = _build_rescue_workspace(tmp_path / "w", lift_gene2=False)
        body = _run_rescue(ws, rescue=True)
        ids = _gene_ids(body)
        assert "gene1" in ids
        assert ids.count("gene2") == 1
        assert "lifton_rescue=miniprot_only" in body
        # the tag rides the rescued mRNA row, and every row is source=LiftOn
        mrna_rows = [ln for ln in body.splitlines() if "\tmRNA\t" in ln]
        assert any("lifton_rescue=miniprot_only" in ln for ln in mrna_rows)
        sources = {ln.split("\t")[1] for ln in body.splitlines()
                   if ln.strip() and not ln.startswith("#")}
        assert sources == {"LiftOn"}

    def test_rescue_dedup_skips_already_lifted_gene(self, tmp_path, hermetic_pipeline):
        # gene2 IS lifted by the DNA lift -> the rescue must NOT add a duplicate
        # (ref-id dedup), and the surviving gene2 is the DNA-lifted one (untagged).
        ws = _build_rescue_workspace(tmp_path / "w", lift_gene2=True)
        body = _run_rescue(ws, rescue=True)
        ids = _gene_ids(body)
        assert ids.count("gene2") == 1
        assert "lifton_rescue=miniprot_only" not in body

    def test_rescue_on_superset_of_off(self, tmp_path, hermetic_pipeline):
        # off ⊆ on: every gene id present with the flag OFF is still present ON.
        ws_off = _build_rescue_workspace(tmp_path / "off", lift_gene2=False)
        off_ids = set(_gene_ids(_run_rescue(ws_off, rescue=False)))
        ws_on = _build_rescue_workspace(tmp_path / "on", lift_gene2=False)
        on_ids = set(_gene_ids(_run_rescue(ws_on, rescue=True)))
        assert off_ids <= on_ids
        assert "gene2" in (on_ids - off_ids)


def _build_xlocus_workspace(work):
    """A reference with ONE coding gene, gene1 (spliced tx1, exon1 101-199 +
    exon2 301-399 -> the same clean 65-aa ORF as ``_build_rescue_workspace``).
    Liftoff lifts gene1 to a WRONG target locus (401-499, pure 'A' filler ->
    no matching ORF exists there -> catastrophically weak protein identity,
    the exact "syntenic stub" failure mode from the Figure-4 investigation).
    miniprot ALSO has a model for tx1 (Target=tx1) but at a THIRD, entirely
    DIFFERENT locus (601-699 + 750-848, two CDS segments mirroring the
    reference's own exon1/exon2 split -- a single contiguous CDS would trip
    the processed-pseudogene filter, since the reference tx1 has >1 exon) that
    carries an exact copy of the true ORF -> protein identity ~1.0 -- the
    cross-locus rescue candidate. Its mRNA span (248) vs gene1's reference
    length (399-101+1=299) gives ratio 0.83, BELOW the default -min_miniprot
    band (0.9), so default Step 8 declines it (mirrors the out-of-band
    condition `_build_rescue_workspace` uses for gene2) -- only the
    cross-locus-rescue pass can recover it."""
    work.mkdir(parents=True, exist_ok=True)
    ref_chrom = ["A"] * 900
    exon1 = "ATG" + "GCT" * 32          # 99 nt
    exon2 = ("GCT" * 32) + "TAA"        # 99 nt
    for i, ch in enumerate(exon1):
        ref_chrom[100 + i] = ch          # 101-199
    for i, ch in enumerate(exon2):
        ref_chrom[300 + i] = ch          # 301-399
    ref_seq = "".join(ref_chrom)

    tgt_chrom = ["A"] * 900             # positions 401-499 stay pure filler
    for i, ch in enumerate(exon1):
        tgt_chrom[600 + i] = ch          # 601-699: cross-locus copy of exon1
    for i, ch in enumerate(exon2):
        tgt_chrom[749 + i] = ch          # 750-848: cross-locus copy of exon2
    tgt_seq = "".join(tgt_chrom)

    ref_fa = work / "ref.fa"; ref_fa.write_text(">chr1\n" + _wrap(ref_seq))
    tgt_fa = work / "tgt.fa"; tgt_fa.write_text(">chr1\n" + _wrap(tgt_seq))

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
    (work / "liftoff.gff3").write_text(
        "##gff-version 3\n"
        "chr1\tLiftoff\tgene\t401\t499\t.\t+\t.\tID=gene1;gene_biotype=protein_coding\n"
        "chr1\tLiftoff\tmRNA\t401\t499\t.\t+\t.\tID=tx1;Parent=gene1\n"
        "chr1\tLiftoff\texon\t401\t499\t.\t+\t.\tID=exon1w;Parent=tx1\n"
        "chr1\tLiftoff\tCDS\t401\t499\t.\t+\t0\tID=cds1w;Parent=tx1\n"
    )
    (work / "miniprot.gff3").write_text(
        "##gff-version 3\n"
        "chr1\tminiprot\tmRNA\t601\t848\t.\t+\t.\tID=MPX;Target=tx1 1 65\n"
        "chr1\tminiprot\tCDS\t601\t699\t.\t+\t0\tID=MPX.cds1;Parent=MPX\n"
        "chr1\tminiprot\tCDS\t750\t848\t.\t+\t0\tID=MPX.cds2;Parent=MPX\n"
    )
    out_dir = work / "out"; out_dir.mkdir()
    return {"ref_fa": ref_fa, "tgt_fa": tgt_fa, "ref_gff": ref_gff,
            "liftoff": work / "liftoff.gff3", "miniprot": work / "miniprot.gff3",
            "out": out_dir}


def _run_xlocus(ws, *, rescue):
    from lifton import lifton as lifton_main
    out_gff = ws["out"] / "lifton.gff3"
    argv = [str(ws["tgt_fa"]), str(ws["ref_fa"]),
            "-g", str(ws["ref_gff"]), "-L", str(ws["liftoff"]),
            "-M", str(ws["miniprot"]), "-o", str(out_gff),
            "-ad", "RefSeq", "--force", "--no-miniprot-rescue"]
    if rescue:
        argv.append("--miniprot-cross-locus-rescue")
    lifton_main.run_all_lifton_steps(lifton_main.parse_args(argv))
    return out_gff.read_text()


def _protein_identity_by_gene(body):
    """gene id -> protein_identity parsed off its mRNA row's attributes."""
    out = {}
    gene_of_parent = {}
    for ln in body.splitlines():
        if not ln.strip() or ln.startswith("#"):
            continue
        cols = ln.split("\t")
        if len(cols) < 9:
            continue
        attrs = dict(kv.split("=", 1) for kv in cols[8].split(";") if "=" in kv)
        if cols[2] == "mRNA" and "protein_identity" in attrs:
            gid = attrs.get("Parent", "")
            out[gid] = float(attrs["protein_identity"])
    return out


class TestCrossLocusRescuePass:
    """Full-pipeline regression coverage for ``lifton.cross_locus_rescue`` — the
    mocked unit tests in tests/test_cross_locus_rescue.py stub out
    run_miniprot.lifton_miniprot_with_ref_protein entirely, so they could not
    catch a bug in the REAL build path. This exercises the real
    Lifton_GENE/Lifton_TRANS/write_entry machinery end-to-end (the same class
    of gap that let candidate-3's update_cds_list corruption slip past its own
    mocked tests earlier this investigation) — it caught a real
    KeyError('coding') in _build_replacement_text (a bare {} passed to
    write_entry's coding/non-coding stats-counter branch instead of the
    {'coding': {}, 'non-coding': {}, 'other': {}} shape lifton.py:757 uses)."""

    def test_rescue_off_keeps_weak_locus(self, tmp_path, hermetic_pipeline):
        ws = _build_xlocus_workspace(tmp_path / "w")
        body = _run_xlocus(ws, rescue=False)
        ids = _gene_ids(body)
        assert ids.count("gene1") == 1
        assert "lifton_rescue=cross_locus" not in body
        pi = _protein_identity_by_gene(body)
        assert pi["gene1"] < 0.5, f"expected the weak Liftoff stub, got {pi}"

    def test_rescue_on_replaces_with_cross_locus_model(self, tmp_path, hermetic_pipeline):
        ws = _build_xlocus_workspace(tmp_path / "w")
        body = _run_xlocus(ws, rescue=True)
        ids = _gene_ids(body)
        # REPLACE, not ADD: gene1 appears EXACTLY once (duplicate-safe).
        assert ids.count("gene1") == 1
        assert "lifton_rescue=cross_locus" in body
        pi = _protein_identity_by_gene(body)
        assert pi["gene1"] > 0.9, f"expected the rescued high-identity model, got {pi}"
        # the surviving gene1 sits at the cross-locus target, not the weak one.
        gene_lines = [ln for ln in body.splitlines()
                     if "\tgene\t" in ln and "ID=gene1" in ln]
        assert len(gene_lines) == 1
        assert "\t601\t798\t" in gene_lines[0] or "\t601\t" in gene_lines[0]

    def test_rescue_flag_off_byte_identical_gene_set(self, tmp_path, hermetic_pipeline):
        # Same gene-id SET both ways (replace, not add) -- unlike Iter-23's
        # off-subset-of-on ADD semantics, Phase D's off/on gene-id sets are
        # equal size; only gene1's CONTENT differs.
        ws_off = _build_xlocus_workspace(tmp_path / "off")
        off_ids = _gene_ids(_run_xlocus(ws_off, rescue=False))
        ws_on = _build_xlocus_workspace(tmp_path / "on")
        on_ids = _gene_ids(_run_xlocus(ws_on, rescue=True))
        assert sorted(off_ids) == sorted(on_ids)
