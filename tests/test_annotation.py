"""Unit tests for lifton.annotation.Annotation."""

from __future__ import annotations

import pytest

from lifton import annotation


# ---------------------------------------------------------------------------
# Database build / cache
# ---------------------------------------------------------------------------

class TestAnnotationDB:
    def test_build_creates_sidecar_db_file(self, gff_standard):
        ann = annotation.Annotation(
            str(gff_standard), False, False, "create_unique", None, True, False,
        )
        # gffutils writes <gff>_db
        sidecar = str(gff_standard) + "_db"
        import os
        assert os.path.exists(sidecar)
        assert ann.db_connection is not None

    def test_lookup_by_id(self, gff_standard):
        ann = annotation.Annotation(
            str(gff_standard), False, False, "create_unique", None, True, False,
        )
        gene = ann.db_connection["gene1"]
        assert gene.start == 101 and gene.end == 399
        assert gene.strand == "+"


# ---------------------------------------------------------------------------
# Coding / non-coding partitioning
# ---------------------------------------------------------------------------

class TestProteinCodingPartition:
    def test_protein_coding_features_returns_gene_id(self, gff_standard):
        ann = annotation.Annotation(
            str(gff_standard), False, False, "create_unique", None, True, False,
        )
        coding = ann.get_protein_coding_features(["gene"])
        assert "gene1" in coding

    def test_noncoding_features_returns_gene_id(self, gff_noncoding):
        ann = annotation.Annotation(
            str(gff_noncoding), False, False, "create_unique", None, True, False,
        )
        nc = ann.get_noncoding_features(["gene"])
        assert "ncg" in nc

    def test_protein_coding_excludes_noncoding(self, gff_noncoding):
        ann = annotation.Annotation(
            str(gff_noncoding), False, False, "create_unique", None, True, False,
        )
        coding = ann.get_protein_coding_features(["gene"])
        assert coding == []


# ---------------------------------------------------------------------------
# get_features_of_type / get_feature_dict
# ---------------------------------------------------------------------------

class TestGetFeatures:
    def test_get_features_of_type_returns_list(self, gff_standard):
        ann = annotation.Annotation(
            str(gff_standard), False, False, "create_unique", None, True, False,
        )
        exons = ann.get_features_of_type(["exon"])
        assert {e.id for e in exons} == {"exon1", "exon2"}

    def test_get_feature_dict_keyed_by_id(self, gff_standard):
        ann = annotation.Annotation(
            str(gff_standard), False, False, "create_unique", None, True, False,
        )
        d = ann.get_feature_dict(["CDS"])
        assert set(d.keys()) == {"cds1", "cds2"}
        assert d["cds1"].featuretype == "CDS"


# ---------------------------------------------------------------------------
# Hierarchy helpers
# ---------------------------------------------------------------------------

class TestHierarchyHelpers:
    def test_is_highest_parent_for_gene(self, gff_standard):
        ann = annotation.Annotation(
            str(gff_standard), False, False, "create_unique", None, True, False,
        )
        assert ann.is_highest_parent("gene1") is True

    def test_is_highest_parent_false_for_child(self, gff_standard):
        ann = annotation.Annotation(
            str(gff_standard), False, False, "create_unique", None, True, False,
        )
        assert ann.is_highest_parent("tx1") is False

    def test_is_lowest_child_true_for_cds(self, gff_standard):
        ann = annotation.Annotation(
            str(gff_standard), False, False, "create_unique", None, True, False,
        )
        assert ann.is_lowest_child("cds1") is True

    def test_get_num_levels_three(self, gff_standard):
        ann = annotation.Annotation(
            str(gff_standard), False, False, "create_unique", None, True, False,
        )
        # gene -> mRNA -> exon/CDS  ⇒  the loop returns the *first* level with
        # zero children, which is 3 (after seeing children at lvl 1 + lvl 2).
        assert ann.get_num_levels("gene1") == 3


# ---------------------------------------------------------------------------
# merge_children_intervals (module-level helper, not the class method)
# ---------------------------------------------------------------------------

class TestMergeChildrenIntervalsModule:
    def test_handles_empty(self):
        assert annotation.merge_children_intervals([]) == []

    def test_overlapping(self):
        from types import SimpleNamespace
        c1 = SimpleNamespace(start=1, end=20)
        c2 = SimpleNamespace(start=10, end=30)
        merged = annotation.merge_children_intervals([c1, c2])
        assert merged == [[1, 30]]
