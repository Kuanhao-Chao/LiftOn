"""Unit tests for lifton.lifton_class — Lifton_GENE / Lifton_TRANS / Lifton_EXON / Lifton_CDS.

The constructors mutate their gffutil_entry arguments. We exercise the
legacy code path on real gffutils.Feature objects pulled from a small
in-memory FeatureDB so attribute handling is identical to production.
"""

from __future__ import annotations

import copy
from types import SimpleNamespace

import pytest
from pyfaidx import Fasta

from lifton import annotation, lifton_class, lifton_utils


# ---------------------------------------------------------------------------
# Construction helpers
# ---------------------------------------------------------------------------

def _gene_attrs(extra_copy=None, gene_biotype="protein_coding"):
    attrs = {"ID": ["gene1"], "gene_biotype": [gene_biotype]}
    if extra_copy is not None:
        attrs["extra_copy_number"] = [str(extra_copy)]
    return attrs


def _trans_attrs(parent="gene1", trans_id="tx1"):
    return {"ID": [trans_id], "Parent": [parent]}


def _build_lifton_gene(db, fake_args, ref_features_dict, tree_dict=None):
    gene_entry = copy.deepcopy(db["gene1"])
    return lifton_class.Lifton_GENE(
        ref_gene_id="gene1",
        gffutil_entry_gene=gene_entry,
        ref_gene_attrs=_gene_attrs(),
        tree_dict=tree_dict if tree_dict is not None else {},
        ref_features_dict=ref_features_dict,
        args=fake_args,
    )


# ---------------------------------------------------------------------------
# Lifton_Status / Lifton_ORF / Lifton_feature
# ---------------------------------------------------------------------------

class TestStatusAndOrfDataclasses:
    def test_status_default_zero(self):
        s = lifton_class.Lifton_Status()
        assert s.lifton_aa == 0 and s.lifton_dna == 0
        assert s.status == [] and s.annotation is None

    def test_orf_round_trips_coords(self):
        o = lifton_class.Lifton_ORF(10, 50)
        assert (o.start, o.end) == (10, 50)

    def test_lifton_feature_has_empty_children(self):
        f = lifton_class.Lifton_feature("xyz")
        assert f.id == "xyz"
        assert f.copy_num == 0
        assert f.is_protein_coding is False
        assert f.children == set()


# ---------------------------------------------------------------------------
# Lifton_GENE constructor
# ---------------------------------------------------------------------------

class TestLiftonGene:
    def test_protein_coding_flag_set_from_gene_biotype(
            self, gff_standard, fake_args, ref_features_dict_one_gene):
        db = annotation.Annotation(
            str(gff_standard), False, False, "create_unique", None, True, False,
        ).db_connection
        gene = _build_lifton_gene(db, fake_args, ref_features_dict_one_gene)
        assert gene.is_protein_coding is True
        assert gene.is_non_coding is False

    def test_noncoding_flag_set_from_lncRNA_biotype(
            self, gff_noncoding, fake_args):
        db = annotation.Annotation(
            str(gff_noncoding), False, False, "create_unique", None, True, False,
        ).db_connection
        ref_features_dict = {"ncg": lifton_class.Lifton_feature("ncg")}
        gene_entry = copy.deepcopy(db["ncg"])
        gene = lifton_class.Lifton_GENE(
            ref_gene_id="ncg",
            gffutil_entry_gene=gene_entry,
            ref_gene_attrs={"ID": ["ncg"], "gene_biotype": ["lncRNA"]},
            tree_dict={},
            ref_features_dict=ref_features_dict,
            args=fake_args,
        )
        assert gene.is_non_coding is True
        assert gene.is_protein_coding is False

    def test_extra_copy_number_appends_suffix(
            self, gff_standard, fake_args, ref_features_dict_one_gene):
        # NOTE: __get_gene_copy reads `extra_copy_number` from the *original*
        # gffutils entry attributes (executed BEFORE ref_gene_attrs is
        # assigned at line 70). So the marker must live on the gff entry.
        db = annotation.Annotation(
            str(gff_standard), False, False, "create_unique", None, True, False,
        ).db_connection
        gene_entry = copy.deepcopy(db["gene1"])
        gene_entry.attributes["extra_copy_number"] = ["2"]
        gene = lifton_class.Lifton_GENE(
            ref_gene_id="gene1",
            gffutil_entry_gene=gene_entry,
            ref_gene_attrs=_gene_attrs(extra_copy=2),
            tree_dict={},
            ref_features_dict=ref_features_dict_one_gene,
            args=fake_args,
        )
        assert gene.copy_num == 2
        # Phase 5 bug fix #1 verified: entry.id now correctly carries the
        # full gene id with the copy-number suffix.
        assert gene.entry.id == "gene1_2"
        assert gene.entry.attributes["ID"] == ["gene1_2"]

    def test_constructor_seeds_tree_dict(
            self, gff_standard, fake_args, ref_features_dict_one_gene):
        db = annotation.Annotation(
            str(gff_standard), False, False, "create_unique", None, True, False,
        ).db_connection
        tree_dict = {}
        _build_lifton_gene(db, fake_args, ref_features_dict_one_gene, tree_dict)
        assert "chr1" in tree_dict
        ivs = list(tree_dict["chr1"])
        assert len(ivs) == 1
        # Phase 5 bug fix #1 verified: tree data is now the full gene id.
        assert ivs[0].data == "gene1"
        assert ivs[0].begin == 101 and ivs[0].end == 399


# ---------------------------------------------------------------------------
# Lifton_TRANS construction + add_exon + add_cds
# ---------------------------------------------------------------------------

class TestLiftonTransBasics:
    def _gene_with_trans(self, gff_standard, fake_args, ref_features_dict_one_gene):
        db = annotation.Annotation(
            str(gff_standard), False, False, "create_unique", None, True, False,
        ).db_connection
        gene = _build_lifton_gene(db, fake_args, ref_features_dict_one_gene)
        trans_entry = copy.deepcopy(db["tx1"])
        trans = gene.add_transcript("tx1", trans_entry, _trans_attrs())
        return db, gene, trans

    def test_transcript_id_assigned(
            self, gff_standard, fake_args, ref_features_dict_one_gene):
        _, _, trans = self._gene_with_trans(
            gff_standard, fake_args, ref_features_dict_one_gene
        )
        assert trans.entry.id == "tx1"
        # Phase 5 bug fix #2 verified: Parent now correctly inherits the
        # full gene id from Lifton_GENE.entry.id.
        assert trans.entry.attributes["Parent"] == ["gene1"]

    def test_add_exon_inserts_in_sorted_order(
            self, gff_standard, fake_args, ref_features_dict_one_gene):
        db, _, trans = self._gene_with_trans(
            gff_standard, fake_args, ref_features_dict_one_gene
        )
        exon2 = copy.deepcopy(db["exon2"])
        exon1 = copy.deepcopy(db["exon1"])
        # Insert in deliberately reversed order, expect sorted by entry.end
        trans.add_exon(exon2)
        trans.add_exon(exon1)
        assert [e.entry.end for e in trans.exons] == [199, 399]

    def test_add_cds_attaches_to_overlapping_exon(
            self, gff_standard, fake_args, ref_features_dict_one_gene):
        db, _, trans = self._gene_with_trans(
            gff_standard, fake_args, ref_features_dict_one_gene
        )
        trans.add_exon(copy.deepcopy(db["exon1"]))
        trans.add_exon(copy.deepcopy(db["exon2"]))
        trans.add_cds(copy.deepcopy(db["cds1"]))
        trans.add_cds(copy.deepcopy(db["cds2"]))
        assert trans.exons[0].cds is not None
        assert trans.exons[1].cds is not None
        assert trans.exons[0].cds.entry.start == 101


# ---------------------------------------------------------------------------
# update_cds_list — Case 1 (single CDS), Case 2 (single exon)
# ---------------------------------------------------------------------------

class TestUpdateCdsListSingleCds:
    def test_case1_single_cds_single_exon_collapses(
            self, gff_single_cds, fake_args):
        db = annotation.Annotation(
            str(gff_single_cds), False, False, "create_unique", None, True, False,
        ).db_connection
        ref_features_dict = {"g_s": lifton_class.Lifton_feature("g_s")}
        gene_entry = copy.deepcopy(db["g_s"])
        gene = lifton_class.Lifton_GENE(
            ref_gene_id="g_s",
            gffutil_entry_gene=gene_entry,
            ref_gene_attrs={"ID": ["g_s"], "gene_biotype": ["protein_coding"]},
            tree_dict={},
            ref_features_dict=ref_features_dict,
            args=fake_args,
        )
        trans = gene.add_transcript(
            "tx_s", copy.deepcopy(db["tx_s"]),
            {"ID": ["tx_s"], "Parent": ["g_s"]},
        )
        trans.add_exon(copy.deepcopy(db["ex_s"]))
        cds = lifton_class.Lifton_CDS(copy.deepcopy(db["cds_s"]))
        trans.update_cds_list([cds])
        # After Case 1 reconciliation we expect a merged exon spanning the
        # original exon (101..199) with a CDS attached.
        assert len(trans.exons) == 1
        merged = trans.exons[0]
        assert (merged.entry.start, merged.entry.end) == (101, 199)
        assert merged.cds is not None
        assert (merged.cds.entry.start, merged.cds.entry.end) == (101, 199)


class TestUpdateCdsListSingleExonMultipleCds:
    def test_case2_single_exon_multiple_cds_splits_exon(
            self, gff_single_cds, fake_args, make_gffutils_feature):
        db = annotation.Annotation(
            str(gff_single_cds), False, False, "create_unique", None, True, False,
        ).db_connection
        ref_features_dict = {"g_s": lifton_class.Lifton_feature("g_s")}
        gene = lifton_class.Lifton_GENE(
            ref_gene_id="g_s",
            gffutil_entry_gene=copy.deepcopy(db["g_s"]),
            ref_gene_attrs={"ID": ["g_s"], "gene_biotype": ["protein_coding"]},
            tree_dict={},
            ref_features_dict=ref_features_dict,
            args=fake_args,
        )
        trans = gene.add_transcript(
            "tx_s", copy.deepcopy(db["tx_s"]),
            {"ID": ["tx_s"], "Parent": ["g_s"]},
        )
        # The fixture has exactly 1 exon -> Case 2 path with two synthetic CDSs
        trans.add_exon(copy.deepcopy(db["ex_s"]))
        cds_a = lifton_class.Lifton_CDS(make_gffutils_feature(
            featuretype="CDS", start=110, end=140, frame="0",
            attributes={"ID": ["cdsA"], "Parent": ["tx_s"]},
        ))
        cds_b = lifton_class.Lifton_CDS(make_gffutils_feature(
            featuretype="CDS", start=160, end=190, frame="0",
            attributes={"ID": ["cdsB"], "Parent": ["tx_s"]},
        ))
        trans.update_cds_list([cds_a, cds_b])
        # Case 2 produces one exon per CDS (len == len(cds_list))
        assert len(trans.exons) == 2
        for exon in trans.exons:
            assert exon.cds is not None


# ---------------------------------------------------------------------------
# update_boundaries
# ---------------------------------------------------------------------------

class TestUpdateBoundaries:
    def test_boundaries_reflect_first_and_last_exon(
            self, gff_standard, fake_args, ref_features_dict_one_gene):
        db = annotation.Annotation(
            str(gff_standard), False, False, "create_unique", None, True, False,
        ).db_connection
        gene = _build_lifton_gene(db, fake_args, ref_features_dict_one_gene)
        trans = gene.add_transcript(
            "tx1", copy.deepcopy(db["tx1"]), _trans_attrs(),
        )
        trans.add_exon(copy.deepcopy(db["exon1"]))
        trans.add_exon(copy.deepcopy(db["exon2"]))
        trans.update_boundaries()
        assert trans.entry.start == 101
        assert trans.entry.end == 399


# ---------------------------------------------------------------------------
# Sequence assembly + translation
# ---------------------------------------------------------------------------

class TestSequenceAssembly:
    def _assemble(self, gff_standard, fasta_standard, fake_args,
                  ref_features_dict_one_gene):
        db = annotation.Annotation(
            str(gff_standard), False, False, "create_unique", None, True, False,
        ).db_connection
        gene = _build_lifton_gene(db, fake_args, ref_features_dict_one_gene)
        trans = gene.add_transcript(
            "tx1", copy.deepcopy(db["tx1"]), _trans_attrs(),
        )
        trans.add_exon(copy.deepcopy(db["exon1"]))
        trans.add_exon(copy.deepcopy(db["exon2"]))
        trans.add_cds(copy.deepcopy(db["cds1"]))
        trans.add_cds(copy.deepcopy(db["cds2"]))
        return trans, Fasta(str(fasta_standard))

    def test_get_coding_seq_starts_with_atg(
            self, gff_standard, fasta_standard, fake_args,
            ref_features_dict_one_gene):
        trans, fa = self._assemble(
            gff_standard, fasta_standard, fake_args, ref_features_dict_one_gene
        )
        coding_seq, cds_children, cdss_lens = trans.get_coding_seq(fa)
        assert coding_seq.startswith("ATG")
        assert sum(cdss_lens) == len(coding_seq)
        assert len(cds_children) == 2

    def test_translate_coding_seq_gives_valid_protein(
            self, gff_standard, fasta_standard, fake_args,
            ref_features_dict_one_gene):
        trans, fa = self._assemble(
            gff_standard, fasta_standard, fake_args, ref_features_dict_one_gene
        )
        coding_seq, _ = trans.get_coding_trans_seq(fa)
        protein = trans.translate_coding_seq(coding_seq)
        assert protein.startswith("M")
        assert protein.endswith("*")
        assert protein.count("*") == 1

    def test_translate_empty_returns_none(
            self, gff_standard, fasta_standard, fake_args,
            ref_features_dict_one_gene):
        trans, _ = self._assemble(
            gff_standard, fasta_standard, fake_args, ref_features_dict_one_gene
        )
        assert trans.translate_coding_seq("") is None


# ---------------------------------------------------------------------------
# ORF rescue: orf_search_protein on a clean transcript should NOT mutate
# (no mutations detected -> no ORF search triggered)
# ---------------------------------------------------------------------------

class TestOrfRescue:
    def test_clean_transcript_no_mutations(
            self, gff_standard, fasta_standard, fake_args,
            ref_features_dict_one_gene):
        db = annotation.Annotation(
            str(gff_standard), False, False, "create_unique", None, True, False,
        ).db_connection
        gene = _build_lifton_gene(db, fake_args, ref_features_dict_one_gene)
        trans = gene.add_transcript(
            "tx1", copy.deepcopy(db["tx1"]), _trans_attrs(),
        )
        trans.add_exon(copy.deepcopy(db["exon1"]))
        trans.add_exon(copy.deepcopy(db["exon2"]))
        trans.add_cds(copy.deepcopy(db["cds1"]))
        trans.add_cds(copy.deepcopy(db["cds2"]))

        fa = Fasta(str(fasta_standard))
        coding_seq, _ = trans.get_coding_trans_seq(fa)
        ref_protein = trans.translate_coding_seq(coding_seq)
        # Ref protein equals translated protein -> identity status
        status = lifton_class.Lifton_Status()
        lifton_tran_aln, lifton_aa_aln = trans.orf_search_protein(
            fa, ref_protein, coding_seq, status, is_non_coding=False,
        )
        assert lifton_aa_aln is not None
        assert lifton_aa_aln.identity == pytest.approx(1.0)
        # No mutation tags pushed onto attributes when sequences match
        assert "mutation" not in trans.entry.attributes


# ---------------------------------------------------------------------------
# Serialization (write_entry chain)
# ---------------------------------------------------------------------------

class TestWriteEntry:
    def test_write_entry_emits_gene_then_trans_then_exon_then_cds(
            self, gff_standard, fake_args, ref_features_dict_one_gene, tmp_path):
        db = annotation.Annotation(
            str(gff_standard), False, False, "create_unique", None, True, False,
        ).db_connection
        gene = _build_lifton_gene(db, fake_args, ref_features_dict_one_gene)
        trans = gene.add_transcript(
            "tx1", copy.deepcopy(db["tx1"]), _trans_attrs(),
        )
        trans.add_exon(copy.deepcopy(db["exon1"]))
        trans.add_exon(copy.deepcopy(db["exon2"]))
        trans.add_cds(copy.deepcopy(db["cds1"]))
        trans.add_cds(copy.deepcopy(db["cds2"]))

        out = tmp_path / "out.gff3"
        with open(out, "w") as fw:
            gene.write_entry(fw, {"coding": {}, "non-coding": {}, "other": {}})
        body = out.read_text()
        # Order: gene → mRNA → 2x exon → 2x CDS
        feature_types = [line.split("\t")[2]
                         for line in body.splitlines() if line.strip()]
        assert feature_types[0] == "gene"
        assert feature_types[1] == "mRNA"
        assert feature_types.count("exon") == 2
        assert feature_types.count("CDS") == 2
        # Source rewritten to "LiftOn"
        assert all(line.split("\t")[1] == "LiftOn"
                   for line in body.splitlines() if line.strip())

    def test_tmp_gene_skips_gene_line(
            self, gff_standard, fake_args, ref_features_dict_one_gene, tmp_path):
        db = annotation.Annotation(
            str(gff_standard), False, False, "create_unique", None, True, False,
        ).db_connection
        gene_entry = copy.deepcopy(db["gene1"])
        gene = lifton_class.Lifton_GENE(
            ref_gene_id="gene1",
            gffutil_entry_gene=gene_entry,
            ref_gene_attrs=_gene_attrs(),
            tree_dict={},
            ref_features_dict=ref_features_dict_one_gene,
            args=fake_args,
            tmp=True,
        )
        gene.add_transcript(
            "tx1", copy.deepcopy(db["tx1"]), _trans_attrs(),
        )
        out = tmp_path / "tmp.gff3"
        with open(out, "w") as fw:
            gene.write_entry(fw, {"coding": {}, "non-coding": {}, "other": {}})
        body = out.read_text()
        feature_types = [line.split("\t")[2]
                         for line in body.splitlines() if line.strip()]
        assert "gene" not in feature_types
        assert "mRNA" in feature_types


# ---------------------------------------------------------------------------
# Lifton_EXON / Lifton_CDS leaf classes
# ---------------------------------------------------------------------------

class TestLeafClasses:
    def test_exon_starts_with_no_cds(self, make_gffutils_feature):
        feat = make_gffutils_feature(
            featuretype="exon", start=100, end=200,
            attributes={"ID": ["e1"], "Parent": ["t1"]},
        )
        ex = lifton_class.Lifton_EXON(feat)
        assert ex.cds is None
        assert ex.entry.source == "LiftOn"

    def test_exon_update_clears_cds(self, make_gffutils_feature):
        feat = make_gffutils_feature(
            featuretype="exon", start=100, end=200,
            attributes={"ID": ["e1"], "Parent": ["t1"]},
        )
        ex = lifton_class.Lifton_EXON(feat)
        ex.cds = "fake-cds"
        ex.update_exon_info(150, 250)
        assert ex.cds is None
        assert ex.entry.start == 150 and ex.entry.end == 250

    def test_extra_copy_number_attribute_stripped(self, make_gffutils_feature):
        feat = make_gffutils_feature(
            featuretype="exon", start=1, end=10,
            attributes={"ID": ["e1"], "Parent": ["t1"],
                        "extra_copy_number": ["3"]},
        )
        ex = lifton_class.Lifton_EXON(feat)
        assert "extra_copy_number" not in ex.entry.attributes
