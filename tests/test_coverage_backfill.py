"""Phase 4.5 Step 6 — coverage backfill for lifton_class.py and
lifton_utils.py to push both modules above 90 %.

Targets the uncovered branches identified by coverage.report:
  - lifton_class.py: Lifton_Alignment.write_alignment, all Lifton_GENE
    helper methods, update_cds_list Case 3 multi-exon × multi-CDS,
    __find_orfs, __update_cds_boundary on both strands, LiftOn_FEATURE.
  - lifton_utils.py: miniprot_id_mapping multi-Target, get_ref_ids_*,
    write_lifton_status / eval_status / chains, get_ref_liffover_features
    branches (RefSeq vs GENCODE), LiftOn_eval/liftoff_alignment.
"""

from __future__ import annotations

import copy
import io
import os
from pathlib import Path
from types import SimpleNamespace

import pytest
from Bio.Seq import Seq
from pyfaidx import Fasta

from lifton import (
    annotation,
    intervals,
    lifton_class,
    lifton_utils,
)


# ===========================================================================
# Lifton_Alignment.write_alignment
# ===========================================================================

class TestLiftonAlignmentWrite:
    def test_write_alignment_creates_subdir_and_file(self, tmp_path):
        aln = lifton_class.Lifton_Alignment(
            extracted_identity=0.95,
            cds_children=None,
            alignment_query="MAGT*",
            alignment_comp="|||||",
            alignment_ref="MAGT*",
            cdss_protein_boundary=None,
            cdss_protein_aln_boundary=None,
            extracted_seq="MAGT*",
            reference_seq="MAGT*",
            db_entry=None,
        )
        aln.write_alignment(str(tmp_path), "miniprot", "frameshift", "tx1")
        out = tmp_path / "miniprot" / "frameshift" / "tx1.fa"
        assert out.exists()
        body = out.read_text()
        assert "> Reference" in body
        assert "> Target" in body
        assert "MAGT*" in body


# ===========================================================================
# LiftOn_FEATURE
# ===========================================================================

class TestLiftOnFeature:
    def test_construct_and_add_child_feature(self, make_gffutils_feature,
                                             tmp_path):
        parent_feat = make_gffutils_feature(
            featuretype="gene", start=100, end=200,
            attributes={"ID": ["g1"], "Parent": []},
        )
        feature = lifton_class.LiftOn_FEATURE("g1", parent_feat, copy_num=0)
        child_feat = make_gffutils_feature(
            featuretype="mRNA", start=120, end=180,
            attributes={"ID": ["tx1"], "Parent": ["g1"]},
        )
        added = feature.add_feature(child_feat)
        assert added.entry.id == "tx1"
        assert "tx1" in feature.features
        # write_entry chain
        out = tmp_path / "feat.gff3"
        with open(out, "w") as fw:
            feature.write_entry(fw)
        body = out.read_text()
        assert "g1" in body
        assert "tx1" in body

    def test_copy_num_appends_suffix_to_feature_id(self,
                                                   make_gffutils_feature):
        # Post-merge: LiftOn_FEATURE calls get_ID_base WITHOUT a
        # ref_features_dict, which is now conservative. The "_5" suffix
        # is preserved, so the renamed id is "tx_5_2" (not "tx_2").
        feat = make_gffutils_feature(
            featuretype="mRNA", start=1, end=10,
            attributes={"ID": ["tx_5"], "Parent": []},
        )
        f = lifton_class.LiftOn_FEATURE("g1", feat, copy_num=2)
        assert f.entry.id == "tx_5_2"

    def test_print_feature_does_not_crash(self, make_gffutils_feature,
                                          capsys):
        feat = make_gffutils_feature(
            featuretype="mRNA", start=1, end=10,
            attributes={"ID": ["tx1"], "Parent": []},
        )
        f = lifton_class.LiftOn_FEATURE("g1", feat, copy_num=0)
        f.print_feature()
        out = capsys.readouterr().out
        assert "tx1" in out


# ===========================================================================
# Lifton_GENE helper methods (update_*, add_*, miniprot_transcript)
# ===========================================================================

def _build_gene(gff_standard, fake_args, ref_features_dict_one_gene):
    db = annotation.Annotation(
        str(gff_standard), False, False, "create_unique", None, True, False,
    ).db_connection
    gene_entry = copy.deepcopy(db["gene1"])
    gene = lifton_class.Lifton_GENE(
        ref_gene_id="gene1",
        gffutil_entry_gene=gene_entry,
        ref_gene_attrs={"ID": ["gene1"], "gene_biotype": ["protein_coding"]},
        tree_dict={},
        ref_features_dict=ref_features_dict_one_gene,
        args=fake_args,
    )
    return db, gene


class TestLiftonGeneHelpers:
    def test_update_gene_info_mutates_seqid_start_end(
            self, gff_standard, fake_args, ref_features_dict_one_gene):
        _, gene = _build_gene(gff_standard, fake_args,
                              ref_features_dict_one_gene)
        gene.update_gene_info("chr2", 5000, 6000)
        assert gene.entry.seqid == "chr2"
        assert gene.entry.start == 5000
        assert gene.entry.end == 6000

    def test_add_miniprot_transcript_routes_to_lifton_trans(
            self, gff_standard, fake_args, ref_features_dict_one_gene):
        db, gene = _build_gene(gff_standard, fake_args,
                               ref_features_dict_one_gene)
        trans = gene.add_miniprot_transcript(
            ref_trans_id="tx1",
            gffutil_entry_trans=copy.deepcopy(db["tx1"]),
            ref_trans_attrs={"ID": ["tx1"], "Parent": ["gene1"]},
            ref_features_dict=ref_features_dict_one_gene,
        )
        assert trans.entry.id == "tx1"
        assert "tx1" in gene.transcripts

    def test_update_trans_info_mutates(self, gff_standard, fake_args,
                                       ref_features_dict_one_gene):
        db, gene = _build_gene(gff_standard, fake_args,
                               ref_features_dict_one_gene)
        gene.add_transcript("tx1", copy.deepcopy(db["tx1"]),
                            {"ID": ["tx1"], "Parent": ["gene1"]})
        gene.update_trans_info("tx1", "chr3", 1, 999)
        assert gene.transcripts["tx1"].entry.seqid == "chr3"

    def test_add_feature_via_LiftOn_FEATURE(
            self, gff_standard, fake_args, ref_features_dict_one_gene):
        db, gene = _build_gene(gff_standard, fake_args,
                               ref_features_dict_one_gene)
        feat = gene.add_feature(copy.deepcopy(db["tx1"]))
        assert feat.entry.id in gene.transcripts

    def test_add_exon_and_add_cds_dispatch_to_transcript(
            self, gff_standard, fake_args, ref_features_dict_one_gene):
        db, gene = _build_gene(gff_standard, fake_args,
                               ref_features_dict_one_gene)
        gene.add_transcript("tx1", copy.deepcopy(db["tx1"]),
                            {"ID": ["tx1"], "Parent": ["gene1"]})
        gene.add_exon("tx1", copy.deepcopy(db["exon1"]))
        gene.add_exon("tx1", copy.deepcopy(db["exon2"]))
        gene.add_cds("tx1", copy.deepcopy(db["cds1"]))
        gene.add_cds("tx1", copy.deepcopy(db["cds2"]))
        assert len(gene.transcripts["tx1"].exons) == 2

    def test_add_lifton_gene_status_attrs(
            self, gff_standard, fake_args, ref_features_dict_one_gene):
        _, gene = _build_gene(gff_standard, fake_args,
                              ref_features_dict_one_gene)
        gene.add_lifton_gene_status_attrs("liftoff")
        assert gene.entry.attributes["source"] == ["liftoff"]

    def test_add_lifton_trans_status_attrs(
            self, gff_standard, fake_args, ref_features_dict_one_gene):
        db, gene = _build_gene(gff_standard, fake_args,
                               ref_features_dict_one_gene)
        gene.add_transcript("tx1", copy.deepcopy(db["tx1"]),
                            {"ID": ["tx1"], "Parent": ["gene1"]})
        status = lifton_class.Lifton_Status()
        status.lifton_aa = 0.987
        status.lifton_dna = 0.876
        status.annotation = "lifton"
        gene.add_lifton_trans_status_attrs("tx1", status)
        attrs = gene.transcripts["tx1"].entry.attributes
        assert attrs["protein_identity"] == ["0.987"]
        assert attrs["dna_identity"] == ["0.876"]
        assert attrs["status"] == ["lifton"]

    def test_print_gene_does_not_crash(
            self, gff_standard, fake_args, ref_features_dict_one_gene,
            capsys):
        _, gene = _build_gene(gff_standard, fake_args,
                              ref_features_dict_one_gene)
        gene.print_gene()
        # Just confirm output appeared
        assert len(capsys.readouterr().out) > 0

    def test_update_cds_list_delegates_to_transcript(
            self, gff_standard, fake_args, ref_features_dict_one_gene):
        db, gene = _build_gene(gff_standard, fake_args,
                               ref_features_dict_one_gene)
        gene.add_transcript("tx1", copy.deepcopy(db["tx1"]),
                            {"ID": ["tx1"], "Parent": ["gene1"]})
        gene.add_exon("tx1", copy.deepcopy(db["exon1"]))
        gene.add_exon("tx1", copy.deepcopy(db["exon2"]))
        cds_list = [
            lifton_class.Lifton_CDS(copy.deepcopy(db["cds1"])),
            lifton_class.Lifton_CDS(copy.deepcopy(db["cds2"])),
        ]
        gene.update_cds_list("tx1", cds_list)
        # Boundaries should reflect children
        assert gene.entry.start <= 101 and gene.entry.end >= 399


# ===========================================================================
# Lifton_GENE.orf_search_protein dispatch (with eval_liftoff_chm13)
# ===========================================================================

class TestGeneOrfSearchDispatch:
    def test_orf_search_protein_returns_none_when_trans_id_missing(
            self, gff_standard, fake_args, ref_features_dict_one_gene,
            fasta_standard):
        _, gene = _build_gene(gff_standard, fake_args,
                              ref_features_dict_one_gene)
        fa = Fasta(str(fasta_standard))
        status = lifton_class.Lifton_Status()
        result = gene.orf_search_protein(
            "missing_tx", "ref_tx", fa, {}, {}, status,
        )
        assert result == (None, False)

    def test_orf_search_protein_with_eval_liftoff_chm13_prefix(
            self, gff_standard, fake_args, ref_features_dict_one_gene,
            fasta_standard):
        db, gene = _build_gene(gff_standard, fake_args,
                               ref_features_dict_one_gene)
        gene.add_transcript("tx1", copy.deepcopy(db["tx1"]),
                            {"ID": ["tx1"], "Parent": ["gene1"]})
        for ex_id in ("exon1", "exon2"):
            gene.add_exon("tx1", copy.deepcopy(db[ex_id]))
        for cds_id in ("cds1", "cds2"):
            gene.add_cds("tx1", copy.deepcopy(db[cds_id]))

        fa = Fasta(str(fasta_standard))
        # eval_liftoff_chm13=True triggers the "rna-" prefix lookup path
        ref_proteins = {"rna-tx1": "MAAAAA*"}
        ref_trans = {"rna-tx1": "ATG" + "GCT" * 5 + "TAA"}
        status = lifton_class.Lifton_Status()
        result = gene.orf_search_protein(
            "tx1", "tx1", fa, ref_proteins, ref_trans, status,
            eval_only=True, eval_liftoff_chm13=True,
        )
        assert result is not None


# ===========================================================================
# update_cds_list — Case 3 (multi-CDS × multi-exon)
# ===========================================================================

class TestUpdateCdsListCase3:
    def test_multi_cds_multi_exon_head_order_True(
            self, gff_standard, fake_args, ref_features_dict_one_gene,
            make_gffutils_feature):
        db, gene = _build_gene(gff_standard, fake_args,
                               ref_features_dict_one_gene)
        trans = gene.add_transcript(
            "tx1", copy.deepcopy(db["tx1"]),
            {"ID": ["tx1"], "Parent": ["gene1"]},
        )
        # Exons span 100..200 and 300..400 — exon comes first (head order True)
        ex_a = make_gffutils_feature(featuretype="exon", start=100, end=200,
                                     attributes={"ID": ["exA"], "Parent": ["tx1"]})
        ex_b = make_gffutils_feature(featuretype="exon", start=300, end=400,
                                     attributes={"ID": ["exB"], "Parent": ["tx1"]})
        trans.add_exon(ex_a)
        trans.add_exon(ex_b)
        # Two CDSs that together overlap both exons
        cds_a = lifton_class.Lifton_CDS(make_gffutils_feature(
            featuretype="CDS", start=150, end=200, frame="0",
            attributes={"ID": ["cdsA"], "Parent": ["tx1"]},
        ))
        cds_b = lifton_class.Lifton_CDS(make_gffutils_feature(
            featuretype="CDS", start=300, end=380, frame="0",
            attributes={"ID": ["cdsB"], "Parent": ["tx1"]},
        ))
        trans.update_cds_list([cds_a, cds_b])
        # Assert: at least the original two exon-shaped boundaries are
        # preserved as separate exons (Case 3 does not collapse them)
        assert len(trans.exons) >= 2

    def test_multi_cds_multi_exon_negative_strand(
            self, gff_standard, fake_args, ref_features_dict_one_gene,
            make_gffutils_feature):
        db, gene = _build_gene(gff_standard, fake_args,
                               ref_features_dict_one_gene)
        trans = gene.add_transcript(
            "tx1", copy.deepcopy(db["tx1"]),
            {"ID": ["tx1"], "Parent": ["gene1"]},
        )
        # Force the transcript onto the - strand to drive cds_list.reverse()
        trans.entry.strand = "-"
        ex_a = make_gffutils_feature(featuretype="exon", start=100, end=200,
                                     strand="-",
                                     attributes={"ID": ["exA"], "Parent": ["tx1"]})
        ex_b = make_gffutils_feature(featuretype="exon", start=300, end=400,
                                     strand="-",
                                     attributes={"ID": ["exB"], "Parent": ["tx1"]})
        trans.add_exon(ex_a)
        trans.add_exon(ex_b)
        cds_a = lifton_class.Lifton_CDS(make_gffutils_feature(
            featuretype="CDS", start=150, end=200, strand="-", frame="0",
            attributes={"ID": ["cdsA"], "Parent": ["tx1"]},
        ))
        cds_b = lifton_class.Lifton_CDS(make_gffutils_feature(
            featuretype="CDS", start=300, end=380, strand="-", frame="0",
            attributes={"ID": ["cdsB"], "Parent": ["tx1"]},
        ))
        trans.update_cds_list([cds_a, cds_b])
        # Did not crash; exons preserved
        assert len(trans.exons) >= 2


# ===========================================================================
# Frame computation (private helper exposed indirectly)
# ===========================================================================

class TestCDSFrameAccumulator:
    def test_frame_after_exact_codon_boundary(
            self, gff_standard, fake_args, ref_features_dict_one_gene,
            fasta_standard):
        db, gene = _build_gene(gff_standard, fake_args,
                               ref_features_dict_one_gene)
        trans = gene.add_transcript(
            "tx1", copy.deepcopy(db["tx1"]),
            {"ID": ["tx1"], "Parent": ["gene1"]},
        )
        for ex_id in ("exon1", "exon2"):
            trans.add_exon(copy.deepcopy(db[ex_id]))
        for cds_id in ("cds1", "cds2"):
            trans.add_cds(copy.deepcopy(db[cds_id]))
        # Both CDSs are 99 bp = multiple of 3 → frame on second exon == 0
        fa = Fasta(str(fasta_standard))
        coding_seq, trans_seq = trans.get_coding_trans_seq(fa)
        # Frame stored on second CDS should be "0"
        assert trans.exons[1].cds.entry.frame == "0"


# ===========================================================================
# lifton_utils — get_ref_ids_liftoff branches
# ===========================================================================

class TestGetRefIdsLiftoff:
    def _make_dict(self):
        d = {
            "gene1": lifton_class.Lifton_feature("gene1"),
            "tx1": lifton_class.Lifton_feature("tx1"),
        }
        return d

    def test_gene_only_path(self):
        d = self._make_dict()
        ref_gene, ref_trans = lifton_utils.get_ref_ids_liftoff(
            d, "gene1", None,
        )
        assert ref_gene == "gene1"
        assert ref_trans is None

    def test_trans_only_path(self):
        d = self._make_dict()
        ref_gene, ref_trans = lifton_utils.get_ref_ids_liftoff(
            d, None, "tx1",
        )
        assert ref_gene is None
        assert ref_trans == "tx1"

    def test_trans_only_with_extra_copy_suffix(self):
        d = self._make_dict()
        ref_gene, ref_trans = lifton_utils.get_ref_ids_liftoff(
            d, None, "tx1_2",
        )
        assert ref_trans == "tx1"

    def test_both_known(self):
        # Post-merge: with both ids supplied, get_ref_ids_liftoff calls
        # get_ID_base(trans_id, None) on the trans path (line 489 in
        # lifton_utils.py), which is now conservative and preserves the
        # _<int> suffix.
        d = self._make_dict()
        ref_gene, ref_trans = lifton_utils.get_ref_ids_liftoff(
            d, "gene1", "tx1_2",
        )
        assert ref_gene == "gene1"
        assert ref_trans == "tx1_2"

    def test_both_unknown_returns_none_pair(self):
        d = self._make_dict()
        ref_gene, ref_trans = lifton_utils.get_ref_ids_liftoff(
            d, "ghost", "phantom",
        )
        assert (ref_gene, ref_trans) == (None, None)


class TestGetRefIdsMiniprot:
    def test_known_miniprot_id(self):
        m_id_2_ref = {"MP1": "tx1"}
        rev = {"tx1": "gene1"}
        assert lifton_utils.get_ref_ids_miniprot(rev, "MP1", m_id_2_ref) == \
               ("gene1", "tx1")

    def test_unknown_miniprot_id_returns_pair_of_none(self):
        m_id_2_ref = {"MP1": "tx1"}
        assert lifton_utils.get_ref_ids_miniprot({}, "MP_unknown", m_id_2_ref) == \
               (None, None)

    def test_known_miniprot_unknown_in_reverse_dict(self):
        m_id_2_ref = {"MP1": "tx_orphan"}
        ref_gene, ref_trans = lifton_utils.get_ref_ids_miniprot(
            {}, "MP1", m_id_2_ref,
        )
        assert ref_gene is None
        assert ref_trans == "tx_orphan"


# ===========================================================================
# miniprot_id_mapping
# ===========================================================================

class TestMiniprotIdMapping:
    def _build_db(self, tmp_path, mrnas):
        fp = tmp_path / "mp.gff3"
        body = "##gff-version 3\n"
        for mp_id, target in mrnas:
            body += (f"chr1\tminiprot\tmRNA\t100\t200\t.\t+\t.\t"
                     f"ID={mp_id};Target={target} 1 50\n")
        fp.write_text(body)
        ann = annotation.Annotation(
            str(fp), False, False, "create_unique", None, True, False,
        )
        return ann.db_connection

    def test_single_mapping(self, tmp_path):
        db = self._build_db(tmp_path, [("MP1", "tx1")])
        ref_to_m, m_to_ref = lifton_utils.miniprot_id_mapping(db)
        assert ref_to_m == {"tx1": ["MP1"]}
        assert m_to_ref == {"MP1": "tx1"}

    def test_multiple_miniprots_share_target(self, tmp_path):
        db = self._build_db(tmp_path, [
            ("MP1", "tx1"),
            ("MP2", "tx1"),
            ("MP3", "tx2"),
        ])
        ref_to_m, m_to_ref = lifton_utils.miniprot_id_mapping(db)
        assert set(ref_to_m["tx1"]) == {"MP1", "MP2"}
        assert ref_to_m["tx2"] == ["MP3"]


# ===========================================================================
# write_lifton_status / write_lifton_eval_status / write_lifton_chains
# ===========================================================================

class TestWriteLiftonStatus:
    def test_write_lifton_status_format(self, tmp_path):
        status = lifton_class.Lifton_Status()
        status.liftoff = 0.95
        status.miniprot = 0.92
        status.lifton_dna = 0.99
        status.lifton_aa = 0.98
        status.annotation = "lifton"
        status.status = ["identical"]
        trans = SimpleNamespace(seqid="chr1", start=100, end=200)
        out = tmp_path / "score.txt"
        with open(out, "w") as fw:
            lifton_utils.write_lifton_status(fw, "tx1", trans, status)
        body = out.read_text()
        # Tab-separated, contains all numeric fields and the location
        assert "tx1" in body
        assert "chr1:100-200" in body
        assert "identical" in body

    def test_write_lifton_eval_status_format(self, tmp_path):
        status = lifton_class.Lifton_Status()
        status.eval_dna = 0.9
        status.eval_aa = 0.8
        status.annotation = "eval"
        status.status = ["nonsynonymous"]
        trans = SimpleNamespace(seqid="chr2", start=10, end=20)
        out = tmp_path / "eval.txt"
        with open(out, "w") as fw:
            lifton_utils.write_lifton_eval_status(fw, "tx2", trans, status)
        body = out.read_text()
        assert "tx2" in body
        assert "chr2:10-20" in body

    def test_write_lifton_chains_format(self, tmp_path):
        out = tmp_path / "chains.txt"
        with open(out, "w") as fw:
            lifton_utils.write_lifton_chains(fw, "tx1", ["L1", "M1", "L2"])
        body = out.read_text()
        assert body == "tx1\tL1;M1;L2\n"

    def test_print_lifton_status_emits_to_stderr(self, capsys):
        status = lifton_class.Lifton_Status()
        status.status = ["identical"]
        status.annotation = "lifton"
        trans = SimpleNamespace(seqid="chr1", start=1, end=10)
        lifton_utils.print_lifton_status("tx1", trans, status, DEBUG=True)
        err = capsys.readouterr().err
        assert "tx1" in err
        assert "chr1:1-10" in err


# ===========================================================================
# get_ref_liffover_features — RefSeq + GENCODE branches
# ===========================================================================

@pytest.fixture
def gff_with_gene_biotype_lncrna(tmp_path):
    fp = tmp_path / "lnc.gff3"
    fp.write_text(
        "##gff-version 3\n"
        "chr1\tt\tgene\t100\t200\t.\t+\t.\tID=lncgene;gene_biotype=lncRNA\n"
        "chr1\tt\tlnc_RNA\t100\t200\t.\t+\t.\tID=lnctx;Parent=lncgene\n"
        "chr1\tt\texon\t100\t200\t.\t+\t.\tID=lncex;Parent=lnctx\n"
    )
    return fp


@pytest.fixture
def gff_gencode_protein_coding(tmp_path):
    fp = tmp_path / "gencode.gff3"
    fp.write_text(
        "##gff-version 3\n"
        "chr1\tt\tgene\t100\t300\t.\t+\t.\tID=g1;gene_type=protein_coding\n"
        "chr1\tt\tmRNA\t100\t300\t.\t+\t.\tID=tx1;Parent=g1\n"
        "chr1\tt\texon\t100\t199\t.\t+\t.\tID=ex1;Parent=tx1\n"
        "chr1\tt\texon\t201\t300\t.\t+\t.\tID=ex2;Parent=tx1\n"
        "chr1\tt\tCDS\t100\t199\t.\t+\t0\tID=cds1;Parent=tx1\n"
        "chr1\tt\tCDS\t201\t300\t.\t+\t0\tID=cds2;Parent=tx1\n"
    )
    return fp


class TestGetRefLiftoverFeatures:
    def _args(self, db_label="RefSeq"):
        return SimpleNamespace(
            annotation_database=db_label,
            evaluation=False,
            evaluation_liftoff_chm13=False,
        )

    def test_refseq_protein_coding_branch(self, gff_standard, tmp_path):
        ann = annotation.Annotation(
            str(gff_standard), False, False, "create_unique", None, True, False,
        )
        d, len_d, rev_d, exon_d = lifton_utils.get_ref_liffover_features(
            ["gene"], ann, str(tmp_path), self._args("RefSeq"),
        )
        assert "gene1" in d
        assert d["gene1"].is_protein_coding is True
        assert "tx1" in rev_d
        assert exon_d["tx1"] == 2  # two CDS children

    def test_refseq_lncrna_branch(self, gff_with_gene_biotype_lncrna,
                                  tmp_path):
        ann = annotation.Annotation(
            str(gff_with_gene_biotype_lncrna), False, False, "create_unique",
            None, True, False,
        )
        d, _, _, _ = lifton_utils.get_ref_liffover_features(
            ["gene"], ann, str(tmp_path), self._args("RefSeq"),
        )
        assert d["lncgene"].is_non_coding is True

    def test_gencode_protein_coding_branch(self, gff_gencode_protein_coding,
                                           tmp_path):
        ann = annotation.Annotation(
            str(gff_gencode_protein_coding), False, False, "create_unique",
            None, True, False,
        )
        d, _, _, _ = lifton_utils.get_ref_liffover_features(
            ["gene"], ann, str(tmp_path), self._args("GENCODE"),
        )
        assert d["g1"].is_protein_coding is True


# ===========================================================================
# LiftOn_eval_alignment / LiftOn_liftoff_alignment
# ===========================================================================

class TestEvalAndLiftoffAlignment:
    def test_eval_alignment_returns_none_when_id_missing(self):
        from types import SimpleNamespace
        # Stub a Lifton_TRANS so lifton_parasail_align bails on the
        # ref_proteins lookup
        class _StubTrans:
            entry = SimpleNamespace(strand="+", start=1, end=10)
            exons = []
        status = lifton_class.Lifton_Status()
        result = lifton_utils.LiftOn_eval_alignment(
            _StubTrans(), SimpleNamespace(start=1, end=10),
            None, {}, "missing_id", status,
        )
        assert result is None
        # status not updated
        assert status.lifton_aa == 0

    def test_liftoff_alignment_returns_none_when_id_missing(self):
        from types import SimpleNamespace
        class _StubTrans:
            entry = SimpleNamespace(strand="+", start=1, end=10)
            exons = []
        status = lifton_class.Lifton_Status()
        result = lifton_utils.LiftOn_liftoff_alignment(
            _StubTrans(), SimpleNamespace(start=1, end=10),
            None, {}, "missing_id", status,
        )
        assert result is None
        assert status.liftoff == 0


# ===========================================================================
# check_miniprot_installed (simulate not-installed branch via monkeypatch)
# ===========================================================================

class TestCheckMiniprotInstalled:
    def test_exits_when_miniprot_absent(self, monkeypatch):
        from lifton import run_miniprot
        monkeypatch.setattr(run_miniprot, "check_miniprot_installed",
                            lambda: False)
        with pytest.raises(SystemExit):
            lifton_utils.check_miniprot_installed()
