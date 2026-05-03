"""Unit tests for lifton.intervals.initialize_interval_tree."""

from __future__ import annotations

import pytest

from lifton import annotation, intervals


def test_interval_tree_built_from_gene(gff_standard):
    db = annotation.Annotation(
        str(gff_standard), False, False, "create_unique", None, True, False,
    )
    tree_dict = intervals.initialize_interval_tree(db.db_connection, ["gene"])
    assert "chr1" in tree_dict
    overlaps = tree_dict["chr1"].overlap(150, 160)
    ids = {iv.data for iv in overlaps}
    assert "gene1" in ids


def test_empty_db_yields_empty_dict(tmp_path, gff_standard):
    # Use a feature type that does not appear -> empty dict
    db = annotation.Annotation(
        str(gff_standard), False, False, "create_unique", None, True, False,
    )
    tree_dict = intervals.initialize_interval_tree(db.db_connection, ["pseudogene"])
    assert tree_dict == {}


def test_multiple_genes_share_chromosome(tmp_path):
    fp = tmp_path / "two.gff3"
    fp.write_text(
        "##gff-version 3\n"
        "chr1\tt\tgene\t100\t200\t.\t+\t.\tID=gA;gene_biotype=protein_coding\n"
        "chr1\tt\tgene\t300\t400\t.\t+\t.\tID=gB;gene_biotype=protein_coding\n"
    )
    db = annotation.Annotation(
        str(fp), False, False, "create_unique", None, True, False,
    )
    tree_dict = intervals.initialize_interval_tree(db.db_connection, ["gene"])
    assert {iv.data for iv in tree_dict["chr1"]} == {"gA", "gB"}
