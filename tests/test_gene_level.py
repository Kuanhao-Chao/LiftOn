"""Gene-level, primary-assembly and GeneID-collapsed recall (benchmarks E1)."""

from __future__ import annotations

import json
from pathlib import Path

from benchmarks.compare import evaluator, gene_level


REFSEQ_GFF = """##gff-version 3
chr1\tRefSeq\tgene\t1\t900\t.\t+\t.\tID=gene-A;Dbxref=GeneID:11;Name=A
chr1\tRefSeq\tmRNA\t1\t900\t.\t+\t.\tID=rna-A1;Parent=gene-A;Dbxref=GeneID:11,GenBank:NM_1
chr1\tRefSeq\texon\t1\t900\t.\t+\t.\tID=exon-A1-1;Parent=rna-A1
chr1\tRefSeq\tmRNA\t1\t800\t.\t+\t.\tID=rna-A2;Parent=gene-A
chr6_GL000251v2_alt\tRefSeq\tgene\t1\t900\t.\t+\t.\tID=gene-A-2;Dbxref=GeneID:11;Name=A
chr6_GL000251v2_alt\tRefSeq\tmRNA\t1\t900\t.\t+\t.\tID=rna-A1-2;Parent=gene-A-2;Dbxref=GeneID:11
chr9_KN196479v1_fix\tRefSeq\tgene\t1\t900\t.\t+\t.\tID=gene-B;Dbxref=GeneID:22%2C;Name=B
chr9_KN196479v1_fix\tRefSeq\tmRNA\t1\t900\t.\t+\t.\tID=rna-B1;Parent=gene-B;product=x%3By
chr1_KI270706v1_random\tEnsembl\tgene\t1\t900\t.\t+\t.\tID=ENSG1
chr1_KI270706v1_random\tEnsembl\ttranscript\t1\t900\t.\t+\t.\tID=ENST1;Parent=ENSG1
"""


def _write(tmp_path: Path, name: str, text: str) -> Path:
    path = tmp_path / name
    path.write_text(text)
    return path


def test_parse_attributes_decodes_and_splits():
    attributes = gene_level.parse_attributes(
        "ID=x;Dbxref=GeneID:1,HGNC:HGNC:2;product=a%3Bb%2Cc;empty")
    assert attributes["Dbxref"] == ["GeneID:1", "HGNC:HGNC:2"]
    assert attributes["product"] == ["a;b,c"]
    assert "empty" not in attributes


def test_geneid_from_dbxref():
    assert gene_level.geneid_from_dbxref(["GenBank:NM_1", "GeneID:7"]) == "GeneID:7"
    assert gene_level.geneid_from_dbxref(["GeneID:"]) is None
    assert gene_level.geneid_from_dbxref(None) is None


def test_seqid_class_keeps_unplaced_scaffolds_primary():
    assert gene_level.seqid_class("chr1") == "primary"
    assert gene_level.seqid_class("chr6_GL000251v2_alt") == "alt"
    assert gene_level.seqid_class("chr9_KN196479v1_fix") == "fix"
    assert gene_level.seqid_class("chr1_KI270706v1_random") == "primary"
    assert gene_level.seqid_class("chrUn_KI270302v1") == "primary"
    assert gene_level.seqid_class("NT_187633.1") == "primary"
    assert gene_level.seqid_class(
        "NT_187633.1", frozenset({"NT_187633.1"})) == "alt"


def test_load_transcript_genes_groups_haplotype_copies(tmp_path):
    genes = gene_level.load_transcript_genes(_write(tmp_path, "r.gff", REFSEQ_GFF))
    # The transcript's own GeneID wins; a transcript without one inherits its
    # gene's; alt copies share the primary gene's key.
    assert genes["rna-A1"] == gene_level.TranscriptGene(
        "gene-A", "GeneID:11", "Dbxref:GeneID", "primary")
    assert genes["rna-A2"].gene_key == "GeneID:11"
    assert genes["rna-A1-2"] == gene_level.TranscriptGene(
        "gene-A-2", "GeneID:11", "Dbxref:GeneID", "alt")
    # Percent-encoded comma inside the value is part of the value, not a split.
    assert genes["rna-B1"].gene_key == "GeneID:22,"
    assert genes["rna-B1"].seqid_class == "fix"
    # No cross-reference: fall back to the parent ID without collapsing.
    assert genes["ENST1"] == gene_level.TranscriptGene(
        "ENSG1", "ENSG1", "gene_id", "primary")
    # Sub-feature rows are not transcript containers.
    assert "exon-A1-1" not in genes


def _rows(recovered):
    return [
        {"ref_mrna_id": tx, "recovered": flag, "is_coding": "1"}
        for tx, flag in recovered.items()
    ]


def test_gene_level_summary_separates_the_three_views(tmp_path):
    genes = gene_level.load_transcript_genes(_write(tmp_path, "r.gff", REFSEQ_GFF))
    rows = gene_level.annotate_rows(_rows({
        "rna-A1": "0", "rna-A2": "1",   # primary gene A, one isoform recovered
        "rna-A1-2": "0",                # its alt copy: the target has no locus
        "rna-B1": "0",                  # fix-patch gene, missed
        "ENST1": "True",                # unplaced scaffold gene, recovered
    }), genes)
    summary = gene_level.gene_level_summary(rows)

    assert summary["n_coding_transcripts"] == 5
    assert summary["transcript_recall_coding"] == 0.4
    assert summary["n_coding_genes"] == 4
    assert summary["gene_recall_coding"] == 0.5
    assert summary["n_coding_genes_by_seqid_class"] == {
        "alt": 1, "fix": 1, "primary": 2}
    assert summary["primary"] == {
        "n_coding_transcripts": 3, "transcript_recall_coding": 0.66667,
        "n_coding_genes": 2, "n_coding_genes_recovered": 2,
        "gene_recall_coding": 1.0,
    }
    # GeneID:11 (primary + alt copy) collapses to one recovered group.
    assert summary["geneid_collapsed"] == {
        "n_groups": 3, "n_groups_recovered": 2, "recall": 0.66667}


def test_noncoding_rows_are_ignored_and_unknown_transcripts_still_count():
    rows = gene_level.annotate_rows([
        {"ref_mrna_id": "nc", "recovered": 1, "is_coding": 0},
        {"ref_mrna_id": "orphan", "recovered": 1, "is_coding": 1},
    ], {})
    summary = gene_level.gene_level_summary(rows)
    assert summary["n_coding_genes"] == 1
    assert summary["gene_recall_coding"] == 1.0
    assert rows[1]["ref_gene_id"] == rows[1]["ref_gene_key"] == "orphan"


def test_empty_input_has_null_ratios():
    summary = gene_level.gene_level_summary([])
    assert summary["gene_recall_coding"] is None
    assert summary["primary"]["gene_recall_coding"] is None
    assert summary["geneid_collapsed"]["recall"] is None


def test_cli_scores_existing_tsvs(tmp_path, capsys):
    ref = _write(tmp_path, "r.gff", REFSEQ_GFF)
    tsv = _write(tmp_path, "t.tsv", "ref_mrna_id\trecovered\tis_coding\n"
                 "rna-A1\t1\t1\nrna-A1-2\t0\t1\n")
    out = tmp_path / "out.json"
    assert gene_level.main([
        "--ref-gff", str(ref), "--tsv", f"tool={tsv}", "--json", str(out)]) == 0
    summaries = json.loads(out.read_text())
    assert summaries["tool"]["gene_recall_coding"] == 0.5
    assert summaries["tool"]["primary"]["gene_recall_coding"] == 1.0
    assert summaries["tool"]["geneid_collapsed"]["recall"] == 1.0
    assert "tool" in capsys.readouterr().out


def test_alt_seqids_file(tmp_path):
    path = _write(tmp_path, "alts.txt", "# alternate loci\nNT_1\n\nNT_2  # patch\n")
    assert gene_level.load_alt_seqids(path) == frozenset({"NT_1", "NT_2"})
    assert gene_level.load_alt_seqids(None) == frozenset()


# --- evaluator hook: opt-in only -------------------------------------------

GENOME = "N" * 10 + "ATGAAACCCGGGTTTTAA" + "N" * 12
EVAL_GFF = """##gff-version 3
chr1\tt\tgene\t11\t28\t.\t+\t.\tID=g1;Dbxref=GeneID:5
chr1\tt\tmRNA\t11\t28\t.\t+\t.\tID=t1;Parent=g1
chr1\tt\texon\t11\t28\t.\t+\t.\tID=e1;Parent=t1
chr1\tt\tCDS\t11\t28\t.\t+\t0\tID=c1;Parent=t1
"""


def _evaluate(tmp_path, gene_map):
    fasta = _write(tmp_path, "g.fa", ">chr1\n" + GENOME + "\n")
    ref_gff = _write(tmp_path, "ref.gff3", EVAL_GFF)
    tool_gff = _write(tmp_path, "tool.gff3", EVAL_GFF)
    ref, ref_index = evaluator.build_reference(str(ref_gff), str(fasta),
                                               log=lambda *_: None)
    out_dir = tmp_path / ("with" if gene_map is not None else "without")
    summary = evaluator.evaluate_tool(
        "lifton", str(tool_gff), str(fasta), ref,
        {"id": "x", "species": "s", "cross_species": False},
        out_dir, None, log=lambda *_: None, ref_index=ref_index,
        gene_map=gene_map)
    return summary, out_dir


def test_evaluate_tool_is_unchanged_without_gene_map(tmp_path):
    summary, out_dir = _evaluate(tmp_path, None)
    assert "gene_level" not in summary
    assert "gene_level" not in json.loads(
        (out_dir / "lifton.summary.json").read_text())
    header = (out_dir / "lifton.transcripts.tsv").read_text().splitlines()[0]
    assert header.split("\t")[-1] == "orf_valid"


def test_evaluate_tool_adds_gene_level_block_on_request(tmp_path):
    ref_gff = _write(tmp_path, "ref_for_map.gff3", EVAL_GFF)
    genes = gene_level.load_transcript_genes(ref_gff)
    summary, out_dir = _evaluate(tmp_path, genes)
    block = summary["gene_level"]
    assert block["n_coding_genes"] == 1
    assert block["gene_recall_coding"] == 1.0
    assert block["geneid_collapsed"]["n_groups"] == 1
    written = json.loads((out_dir / "lifton.summary.json").read_text())
    assert written["gene_level"] == block
