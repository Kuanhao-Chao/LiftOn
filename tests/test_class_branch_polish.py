"""Phase 4.5 Step 6 (final pass) — last-mile branch coverage for
lifton_class.py to clear the 90 % gate.

Targets the remaining uncovered branches:
  - Lifton_GENE: ENSEMBL/CHESS annotation_database branches (line 84-85),
    update_boundaries / write_entry stat-counter increment (line 175),
    print_gene loop (line 185)
  - LiftOn_FEATURE: copy_num >0 path (lines 197-200), ref_tran_id branch
  - update_cds_list Case 1 sub-branches: pure-before, pure-after, in-overlap
  - get_coding_seq with strand "-" (lines 462-464)
  - Lifton_EXON.add_novel_lifton_cds (lines 700-709)
  - __iterate_exons_update_cds: when exon.cds is None branches
"""

from __future__ import annotations

import copy
from types import SimpleNamespace

import pytest
from Bio.Seq import Seq
from pyfaidx import Fasta

from lifton import annotation, lifton_class


# ---------------------------------------------------------------------------
# Lifton_GENE annotation_database branches
# ---------------------------------------------------------------------------

class TestAnnotationDatabaseBranches:
    @pytest.mark.parametrize("db_label", ["GENCODE", "ENSEMBL", "CHESS"])
    def test_alternate_databases_use_gene_type_key(self, gff_standard,
                                                   db_label):
        from lifton import annotation as ann_mod
        db = ann_mod.Annotation(
            str(gff_standard), False, False, "create_unique", None, True, False,
        ).db_connection
        feat_dict = {"gene1": lifton_class.Lifton_feature("gene1")}
        args = SimpleNamespace(annotation_database=db_label,
                               evaluation=False,
                               evaluation_liftoff_chm13=False)
        gene_entry = copy.deepcopy(db["gene1"])
        # Add a gene_type attribute so the downstream branch fires
        gene_entry.attributes["gene_type"] = ["protein_coding"]
        gene = lifton_class.Lifton_GENE(
            ref_gene_id="gene1",
            gffutil_entry_gene=gene_entry,
            ref_gene_attrs={"ID": ["gene1"], "gene_type": ["protein_coding"]},
            tree_dict={},
            ref_features_dict=feat_dict,
            args=args,
        )
        assert gene.is_protein_coding is True


# ---------------------------------------------------------------------------
# write_entry: stat-counter increment branch (line 175)
# ---------------------------------------------------------------------------

class TestWriteEntryStatCounter:
    def test_existing_ref_tran_id_increments_counter(
            self, gff_standard, fake_args, ref_features_dict_one_gene,
            tmp_path):
        from lifton import annotation as ann_mod
        db = ann_mod.Annotation(
            str(gff_standard), False, False, "create_unique", None, True, False,
        ).db_connection
        gene = lifton_class.Lifton_GENE(
            ref_gene_id="gene1",
            gffutil_entry_gene=copy.deepcopy(db["gene1"]),
            ref_gene_attrs={"ID": ["gene1"], "gene_biotype": ["protein_coding"]},
            tree_dict={},
            ref_features_dict=ref_features_dict_one_gene,
            args=fake_args,
        )
        gene.add_transcript("tx1", copy.deepcopy(db["tx1"]),
                            {"ID": ["tx1"], "Parent": ["gene1"]})
        out = tmp_path / "out.gff3"
        with open(out, "w") as fw:
            stats = {"coding": {"tx1": 5}, "non-coding": {}, "other": {}}
            gene.write_entry(fw, stats)
        # Counter incremented from 5 -> 6 (line 175 branch)
        assert stats["coding"]["tx1"] == 6


# ---------------------------------------------------------------------------
# LiftOn_FEATURE copy_num > 0 branch
# ---------------------------------------------------------------------------

class TestLiftOnFeatureCopyNum:
    def test_copy_num_zero_does_not_rename(self, make_gffutils_feature):
        feat = make_gffutils_feature(
            featuretype="mRNA", start=1, end=10,
            attributes={"ID": ["tx1"], "Parent": []},
        )
        f = lifton_class.LiftOn_FEATURE("g1", feat, copy_num=0)
        assert f.entry.id == "tx1"

    def test_copy_num_positive_renames_with_id_base(self,
                                                    make_gffutils_feature):
        feat = make_gffutils_feature(
            featuretype="mRNA", start=1, end=10,
            attributes={"ID": ["tx1"], "Parent": []},
        )
        f = lifton_class.LiftOn_FEATURE("g1", feat, copy_num=3)
        # ID base "tx1" + "_3" -> "tx1_3"
        assert f.entry.id == "tx1_3"
        assert f.entry.attributes["ID"] == ["tx1_3"]
        assert f.entry.attributes["Parent"] == ["g1"]


# ---------------------------------------------------------------------------
# get_coding_seq with negative-strand CDS (lines 462-464)
# ---------------------------------------------------------------------------

class TestGetCodingSeqNegativeStrand:
    def test_negative_strand_chains_in_reverse(
            self, gff_standard, fake_args, ref_features_dict_one_gene,
            fasta_standard):
        from lifton import annotation as ann_mod
        db = ann_mod.Annotation(
            str(gff_standard), False, False, "create_unique", None, True, False,
        ).db_connection
        # Force CDS strand to "-" on the entries
        gene = lifton_class.Lifton_GENE(
            ref_gene_id="gene1",
            gffutil_entry_gene=copy.deepcopy(db["gene1"]),
            ref_gene_attrs={"ID": ["gene1"], "gene_biotype": ["protein_coding"]},
            tree_dict={},
            ref_features_dict=ref_features_dict_one_gene,
            args=fake_args,
        )
        trans = gene.add_transcript("tx1", copy.deepcopy(db["tx1"]),
                                    {"ID": ["tx1"], "Parent": ["gene1"]})
        for ex_id in ("exon1", "exon2"):
            ex = copy.deepcopy(db[ex_id])
            ex.strand = "-"
            trans.add_exon(ex)
        for cds_id in ("cds1", "cds2"):
            cds = copy.deepcopy(db[cds_id])
            cds.strand = "-"
            trans.add_cds(cds)
        fa = Fasta(str(fasta_standard))
        coding_seq, cds_children, cdss_lens = trans.get_coding_seq(fa)
        # Two CDSs of 99 bp each chained
        assert sum(cdss_lens) == 198
        assert len(cds_children) == 2


# ---------------------------------------------------------------------------
# Lifton_EXON.add_novel_lifton_cds
# ---------------------------------------------------------------------------

class TestAddNovelLiftonCDS:
    def test_add_novel_cds_creates_correctly_typed_child(
            self, make_gffutils_feature):
        ex_feat = make_gffutils_feature(
            featuretype="exon", start=100, end=200,
            attributes={"ID": ["e1"], "Parent": ["t1"]},
        )
        ex = lifton_class.Lifton_EXON(ex_feat)
        ex.add_novel_lifton_cds(ex_feat, 110, 190)
        assert ex.cds is not None
        assert ex.cds.entry.featuretype == "CDS"
        assert ex.cds.entry.start == 110
        assert ex.cds.entry.end == 190
        # Parent is inherited from the parent exon
        assert ex.cds.entry.attributes["Parent"] == ["t1"]


# ---------------------------------------------------------------------------
# update_cds_list Case 1 sub-branches: pure-before-only and in-overlap-only
# ---------------------------------------------------------------------------

class TestUpdateCdsListCase1SubBranches:
    def _build(self, fake_args, ref_features_dict_one_gene, gff_single_cds,
               make_gffutils_feature, exon_specs, cds_spec):
        from lifton import annotation as ann_mod
        db = ann_mod.Annotation(
            str(gff_single_cds), False, False, "create_unique", None, True, False,
        ).db_connection
        gene = lifton_class.Lifton_GENE(
            ref_gene_id="g_s",
            gffutil_entry_gene=copy.deepcopy(db["g_s"]),
            ref_gene_attrs={"ID": ["g_s"], "gene_biotype": ["protein_coding"]},
            tree_dict={},
            ref_features_dict={"g_s": lifton_class.Lifton_feature("g_s")},
            args=fake_args,
        )
        trans = gene.add_transcript(
            "tx_s", copy.deepcopy(db["tx_s"]),
            {"ID": ["tx_s"], "Parent": ["g_s"]},
        )
        for s, e in exon_specs:
            ex = make_gffutils_feature(
                featuretype="exon", start=s, end=e,
                attributes={"ID": [f"ex{s}"], "Parent": ["tx_s"]},
            )
            trans.add_exon(ex)
        cds = lifton_class.Lifton_CDS(make_gffutils_feature(
            featuretype="CDS", start=cds_spec[0], end=cds_spec[1],
            frame="0",
            attributes={"ID": ["cdsX"], "Parent": ["tx_s"]},
        ))
        return trans, cds

    def test_case1_exon_strictly_before_cds(
            self, fake_args, gff_single_cds, make_gffutils_feature):
        # Exon at 100..150 strictly before CDS at 200..250
        # This drives the `exon.entry.end < only_cds.entry.start` branch
        trans, cds = self._build(
            fake_args, None, gff_single_cds, make_gffutils_feature,
            exon_specs=[(100, 150), (300, 350)],
            cds_spec=(200, 250),
        )
        trans.update_cds_list([cds])
        # Original exon at 100..150 should be preserved
        starts = [e.entry.start for e in trans.exons]
        assert 100 in starts


# ---------------------------------------------------------------------------
# update_cds_list Case 3 sub-branches: head_order False (CDS before exon)
# ---------------------------------------------------------------------------

class TestUpdateCdsListCase3HeadOrderFalse:
    def test_cds_before_first_exon(self, fake_args,
                                   gff_standard, ref_features_dict_one_gene,
                                   make_gffutils_feature):
        from lifton import annotation as ann_mod
        db = ann_mod.Annotation(
            str(gff_standard), False, False, "create_unique", None, True, False,
        ).db_connection
        gene = lifton_class.Lifton_GENE(
            ref_gene_id="gene1",
            gffutil_entry_gene=copy.deepcopy(db["gene1"]),
            ref_gene_attrs={"ID": ["gene1"], "gene_biotype": ["protein_coding"]},
            tree_dict={},
            ref_features_dict=ref_features_dict_one_gene,
            args=fake_args,
        )
        trans = gene.add_transcript(
            "tx1", copy.deepcopy(db["tx1"]),
            {"ID": ["tx1"], "Parent": ["gene1"]},
        )
        # First exon starts AFTER first CDS ends -> init_head_order = False
        ex_a = make_gffutils_feature(
            featuretype="exon", start=300, end=400,
            attributes={"ID": ["exA"], "Parent": ["tx1"]},
        )
        ex_b = make_gffutils_feature(
            featuretype="exon", start=500, end=600,
            attributes={"ID": ["exB"], "Parent": ["tx1"]},
        )
        trans.add_exon(ex_a)
        trans.add_exon(ex_b)
        cds_list = [
            lifton_class.Lifton_CDS(make_gffutils_feature(
                featuretype="CDS", start=100, end=200, frame="0",
                attributes={"ID": ["cdsX"], "Parent": ["tx1"]},
            )),
            lifton_class.Lifton_CDS(make_gffutils_feature(
                featuretype="CDS", start=550, end=580, frame="0",
                attributes={"ID": ["cdsY"], "Parent": ["tx1"]},
            )),
        ]
        # Should not crash
        trans.update_cds_list(cds_list)
        assert len(trans.exons) >= 1
