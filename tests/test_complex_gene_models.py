"""Phase 4.5 Step 3 — complex gene-model topology tests.

Covers nested genes, isoform-shared exons, multi-Parent attributes,
opposite-strand overlaps, gene → CDS prokaryote pattern, and genes with
no transcript children.
"""

from __future__ import annotations

from types import SimpleNamespace

import pytest
from pyfaidx import Fasta

from lifton import annotation, extract_sequence, intervals, lifton_utils


# ---------------------------------------------------------------------------
# Nested genes (gene-in-intron)
# ---------------------------------------------------------------------------

class TestNestedGenes:
    def test_both_genes_present_in_interval_tree(self, gff_nested_genes):
        ann = annotation.Annotation(
            str(gff_nested_genes), False, False, "create_unique", None,
            True, False,
        )
        tree_dict = intervals.initialize_interval_tree(
            ann.db_connection, ["gene"]
        )
        # Nested gene B (220..280) must be findable inside outer gene A's
        # genomic span (101..399).
        ovps = tree_dict["chr1"].overlap(220, 280)
        ids = {iv.data for iv in ovps}
        assert {"gA", "gB"}.issubset(ids)

    def test_inner_gene_does_not_contain_outer_exons(self, gff_nested_genes):
        ann = annotation.Annotation(
            str(gff_nested_genes), False, False, "create_unique", None,
            True, False,
        )
        # Outer gene's two exons are at 101..199 and 301..399. Inner gene's
        # span 220..280 must not overlap either exon.
        ex_a1 = ann.db_connection["exA1"]
        ex_a2 = ann.db_connection["exA2"]
        gB = ann.db_connection["gB"]
        _, ovp1 = lifton_utils.segments_overlap_length(
            (ex_a1.start, ex_a1.end), (gB.start, gB.end)
        )
        _, ovp2 = lifton_utils.segments_overlap_length(
            (ex_a2.start, ex_a2.end), (gB.start, gB.end)
        )
        assert ovp1 is False and ovp2 is False


# ---------------------------------------------------------------------------
# Isoform-shared exons (Parent=tx1,tx2)
# ---------------------------------------------------------------------------

class TestSharedExonAcrossIsoforms:
    def test_shared_exon_resolved_under_both_transcripts(
            self, gff_shared_exon_isoforms):
        ann = annotation.Annotation(
            str(gff_shared_exon_isoforms), False, False, "create_unique",
            None, True, False,
        )
        # gffutils should expose the shared exon as a child of both
        # transcripts because Parent=tx1,tx2 lists both.
        children_tx1 = list(ann.db_connection.children(
            "tx1", featuretype="exon",
        ))
        children_tx2 = list(ann.db_connection.children(
            "tx2", featuretype="exon",
        ))
        ids_tx1 = {c.id for c in children_tx1}
        ids_tx2 = {c.id for c in children_tx2}
        assert "exShared" in ids_tx1
        assert "exShared" in ids_tx2

    def test_multivalue_parent_round_trip(self, gff_shared_exon_isoforms):
        ann = annotation.Annotation(
            str(gff_shared_exon_isoforms), False, False, "create_unique",
            None, True, False,
        )
        shared = ann.db_connection["exShared"]
        # Per NCBI multi-value attr §, Parent must be a list of two ids.
        assert shared.attributes["Parent"] == ["tx1", "tx2"]

    def test_extract_features_handles_shared_exon(
            self, gff_shared_exon_isoforms, fasta_standard):
        ann = annotation.Annotation(
            str(gff_shared_exon_isoforms), False, False, "create_unique",
            None, True, False,
        )
        fa = Fasta(str(fasta_standard))
        ref_trans, ref_proteins = extract_sequence.extract_features(
            ann, ["gene"], fa,
        )
        # Both transcripts should have produced an entry in ref_trans
        assert "tx1" in ref_trans
        assert "tx2" in ref_trans


# ---------------------------------------------------------------------------
# Opposite-strand overlap
# ---------------------------------------------------------------------------

class TestOppositeStrandOverlap:
    def test_both_genes_partitioned_into_tree(
            self, gff_overlapping_opposite_strand):
        ann = annotation.Annotation(
            str(gff_overlapping_opposite_strand), False, False,
            "create_unique", None, True, False,
        )
        tree_dict = intervals.initialize_interval_tree(
            ann.db_connection, ["gene"]
        )
        # Both genes share genomic span; intervaltree uses half-open
        # ranges so query a strict subset of the overlap region.
        ovps = tree_dict["chr1"].overlap(180, 220)
        ids = {iv.data for iv in ovps}
        assert ids == {"gFwd", "gRev"}

    def test_strand_attribute_preserved(
            self, gff_overlapping_opposite_strand):
        ann = annotation.Annotation(
            str(gff_overlapping_opposite_strand), False, False,
            "create_unique", None, True, False,
        )
        assert ann.db_connection["gFwd"].strand == "+"
        assert ann.db_connection["gRev"].strand == "-"


# ---------------------------------------------------------------------------
# Gene with no transcripts
# ---------------------------------------------------------------------------

class TestGeneWithoutChildren:
    def test_gene_only_loads(self, gff_gene_no_transcripts):
        ann = annotation.Annotation(
            str(gff_gene_no_transcripts), False, False, "create_unique",
            None, True, False,
        )
        gene = ann.db_connection["gOnly"]
        assert gene.start == 100 and gene.end == 200
        # No children
        children = list(ann.db_connection.children("gOnly"))
        assert children == []

    def test_extract_features_does_not_crash_for_childless_gene(
            self, gff_gene_no_transcripts, fasta_standard):
        ann = annotation.Annotation(
            str(gff_gene_no_transcripts), False, False, "create_unique",
            None, True, False,
        )
        fa = Fasta(str(fasta_standard))
        ref_trans, ref_proteins = extract_sequence.extract_features(
            ann, ["gene"], fa,
        )
        # No exon/CDS children -> no entries populated
        assert "gOnly" not in ref_trans
        assert "gOnly" not in ref_proteins


# ---------------------------------------------------------------------------
# Prokaryote gene → CDS pattern (NCBI NOTE 2)
# ---------------------------------------------------------------------------

class TestProkaryoteGeneToCDS:
    def test_recurses_into_cds_children(
            self, gff_prokaryote_gene_to_cds, fasta_standard):
        ann = annotation.Annotation(
            str(gff_prokaryote_gene_to_cds), False, False, "create_unique",
            None, True, False,
        )
        fa = Fasta(str(fasta_standard))
        ref_trans, ref_proteins = extract_sequence.extract_features(
            ann, ["gene"], fa,
        )
        # Direct gene → CDS: legacy code path may attribute the protein
        # to the gene id. Pin observed behaviour.
        assert "gP" in ref_proteins or "cdsP" in ref_proteins


# ---------------------------------------------------------------------------
# Trans-spliced / non-contiguous exons
# ---------------------------------------------------------------------------

class TestNonContiguousExons:
    def test_distant_exons_not_merged(self):
        # Two intervals far apart -> merge_children_intervals should keep
        # them separate.
        from types import SimpleNamespace
        a = SimpleNamespace(start=100, end=200)
        b = SimpleNamespace(start=10000, end=10100)
        merged = extract_sequence.merge_children_intervals([a, b])
        assert merged == [[100, 200], [10000, 10100]]

    def test_adjacent_intervals_merged(self):
        from types import SimpleNamespace
        a = SimpleNamespace(start=100, end=200)
        b = SimpleNamespace(start=200, end=300)
        merged = extract_sequence.merge_children_intervals([a, b])
        # 200 == 200, so they touch and current logic merges them
        assert merged == [[100, 300]]
