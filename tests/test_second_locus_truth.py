"""Independent target truth separates new coding loci from existing coverage."""

from pathlib import Path

import pytest
from pyfaidx import Fasta

from benchmarks.compare.second_locus_truth import evaluate


def _row(kind, start, end, attrs, strand="+"):
    phase = "0" if kind == "CDS" else "."
    return f"chr1\ttest\t{kind}\t{start}\t{end}\t.\t{strand}\t{phase}\t{attrs}\n"


def _gene(gene, transcript, start, end, *, second=False):
    tag = ";lifton_rescue_second_locus=true" if second else ""
    return (
        _row("gene", start, end, f"ID={gene};gene_biotype=protein_coding")
        + _row("mRNA", start, end, f"ID={transcript};Parent={gene}{tag}")
        + _row("CDS", start, end, f"ID=cds-{transcript};Parent={transcript}")
    )


def test_new_target_gene_is_supported_but_unannotated_addition_is_not(tmp_path):
    fasta = tmp_path / "target.fa"
    fasta.write_text(">chr1\n" + "A" * 1000 + "\n")
    Fasta(str(fasta)).close()
    truth = tmp_path / "truth.gff3"
    truth.write_text(
        "##gff-version 3\n#!genome-build-accession NCBI_Assembly:GCF_TEST.1\n"
        "##sequence-region chr1 1 1000\n"
        + _gene("truth-old", "truth-old-t", 20, 80)
        + _gene("truth-new", "truth-new-t", 320, 380)
    )
    off = tmp_path / "off.gff3"
    off.write_text("##gff-version 3\n" + _gene("old", "old-t", 20, 80))
    on = tmp_path / "on.gff3"
    on.write_text(off.read_text() + _gene("new", "new-t", 320, 380, second=True)
                  + _gene("novel", "novel-t", 600, 660, second=True))

    report = evaluate(off, on, truth, fasta, "GCF_TEST.1", replicates=0)
    assert report["target_covered_by_off"] == 1
    assert report["added_models"] == 2
    assert report["independently_supported_new_loci"] == 1
    assert report["placements"]["new-t"]["assigned_target_gene"] == "truth-new"
    assert "assigned_target_gene" not in report["placements"]["novel-t"]
    with pytest.raises(ValueError, match="expected GCF_WRONG"):
        evaluate(off, on, truth, fasta, "GCF_WRONG", replicates=0)
