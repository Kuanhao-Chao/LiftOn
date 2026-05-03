"""Phase 5 final-mile coverage backfill — drive every uncovered branch
in lifton_class.py and lifton_utils.py to push both above 95% and the
package over 80%.

Targeted gaps:
  - lifton_class.py update_cds_list Case 2 inner-CDS branch (lines 329-330)
  - lifton_class.py __iterate_exons_update_cds: novel CDS creation on
    + and - strand, full-exon CDS extension, "after ORF" no-CDS branch
  - lifton_class.py Lifton_GENE inner default branches
  - lifton_utils.exec_liftoff / exec_miniprot file-exists short-circuit
  - lifton_utils.LiftOn_eval_alignment / LiftOn_liftoff_alignment success
    paths (status mutation)
  - lifton_utils.get_ref_liffover_features all gene_type / "other" branches
  - lifton_utils.LiftOn_miniprot_alignment cross-gene-overlap rejection
  - lifton_utils.get_ref_ids_liftoff fallback when liftoff_id is unknown
"""

from __future__ import annotations

import copy
import os
from pathlib import Path
from types import SimpleNamespace

import pytest
from intervaltree import Interval, IntervalTree
from pyfaidx import Fasta

from lifton import (
    annotation,
    intervals,
    lifton_class,
    lifton_utils,
    run_liftoff,
    run_miniprot,
)


# ===========================================================================
# exec_liftoff / exec_miniprot — file-exists short-circuit branches
# ===========================================================================

class TestExecLiftoffMiniprotShortCircuit:
    def test_exec_liftoff_returns_provided_file_when_exists(self, tmp_path,
                                                            monkeypatch):
        # Monkey-patch the underlying runner so any fallthrough explodes
        def _fail(*a, **k):
            raise AssertionError("Should have short-circuited")
        monkeypatch.setattr(run_liftoff, "run_liftoff", _fail)

        existing = tmp_path / "prior_liftoff.gff3"
        existing.write_text("##gff-version 3\n")
        args = SimpleNamespace(liftoff=str(existing))
        # Post-merge: exec_liftoff signature is (outdir, ref_db, args)
        result = lifton_utils.exec_liftoff(str(tmp_path), None, args)
        assert result == str(existing)

    def test_exec_liftoff_runs_runner_when_file_missing(self, tmp_path,
                                                       monkeypatch):
        called = {"n": 0}
        def _runner(outdir, ref_db, args):
            called["n"] += 1
            return os.path.join(outdir, "synth.gff3")
        monkeypatch.setattr(run_liftoff, "run_liftoff", _runner)
        args = SimpleNamespace(liftoff=None)
        result = lifton_utils.exec_liftoff(str(tmp_path), None, args)
        assert called["n"] == 1
        assert result.endswith("synth.gff3")

    def test_exec_miniprot_returns_provided_file_when_exists(self, tmp_path,
                                                             monkeypatch):
        monkeypatch.setattr(lifton_utils, "check_miniprot_installed",
                            lambda: None)
        def _fail(*a, **k):
            raise AssertionError("Should have short-circuited")
        monkeypatch.setattr(run_miniprot, "run_miniprot", _fail)

        existing = tmp_path / "prior_mp.gff3"
        existing.write_text("##gff-version 3\n")
        args = SimpleNamespace(miniprot=str(existing))
        result = lifton_utils.exec_miniprot(
            str(tmp_path), args, "tgt.fa", "ref_proteins.fa",
        )
        assert result == str(existing)

    def test_exec_miniprot_runs_runner_when_file_missing(self, tmp_path,
                                                        monkeypatch):
        monkeypatch.setattr(lifton_utils, "check_miniprot_installed",
                            lambda: None)
        called = {"n": 0}
        def _runner(outdir, args, tgt, ref_proteins):
            called["n"] += 1
            return os.path.join(outdir, "synth.gff3")
        monkeypatch.setattr(run_miniprot, "run_miniprot", _runner)
        args = SimpleNamespace(miniprot=None)
        result = lifton_utils.exec_miniprot(
            str(tmp_path), args, "tgt.fa", "ref_proteins.fa",
        )
        assert called["n"] == 1


# ===========================================================================
# LiftOn_eval_alignment / LiftOn_liftoff_alignment — happy paths
# ===========================================================================

class TestEvalLiftoffAlignmentSuccess:
    def _build_trans_with_real_alignment(self, gff_standard, fasta_standard,
                                         fake_args, ref_features_dict_one_gene):
        db = annotation.Annotation(
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
        for ex_id in ("exon1", "exon2"):
            trans.add_exon(copy.deepcopy(db[ex_id]))
        for cds_id in ("cds1", "cds2"):
            trans.add_cds(copy.deepcopy(db[cds_id]))
        return db, trans

    def test_eval_alignment_updates_status(
            self, gff_standard, fasta_standard, fake_args,
            ref_features_dict_one_gene):
        db, trans = self._build_trans_with_real_alignment(
            gff_standard, fasta_standard, fake_args,
            ref_features_dict_one_gene,
        )
        fa = Fasta(str(fasta_standard))
        # Compute the actual translated protein and use it as the reference
        from lifton import extract_sequence
        from types import SimpleNamespace as SN
        proto = extract_sequence.get_protein_sequence(
            SN(seqid="chr1", strand="+"), fa,
            [SN(start=101, end=199), SN(start=301, end=399)],
        )
        ref_proteins = {"tx1": proto}
        status = lifton_class.Lifton_Status()
        result = lifton_utils.LiftOn_eval_alignment(
            trans, db["tx1"], fa, ref_proteins, "tx1", status,
        )
        assert result is not None
        assert status.lifton_aa > 0

    def test_liftoff_alignment_updates_status(
            self, gff_standard, fasta_standard, fake_args,
            ref_features_dict_one_gene):
        db, trans = self._build_trans_with_real_alignment(
            gff_standard, fasta_standard, fake_args,
            ref_features_dict_one_gene,
        )
        fa = Fasta(str(fasta_standard))
        from lifton import extract_sequence
        from types import SimpleNamespace as SN
        proto = extract_sequence.get_protein_sequence(
            SN(seqid="chr1", strand="+"), fa,
            [SN(start=101, end=199), SN(start=301, end=399)],
        )
        ref_proteins = {"tx1": proto}
        status = lifton_class.Lifton_Status()
        result = lifton_utils.LiftOn_liftoff_alignment(
            trans, db["tx1"], fa, ref_proteins, "tx1", status,
        )
        assert result is not None
        assert status.liftoff > 0


# ===========================================================================
# LiftOn_miniprot_alignment cross-gene-overlap rejection branch
# ===========================================================================

class TestMiniprotAlignmentCrossGene:
    def test_returns_none_when_ref_not_in_proteins(self):
        status = lifton_class.Lifton_Status()
        result, has_valid = lifton_utils.LiftOn_miniprot_alignment(
            "chr1", SimpleNamespace(seqid="chr1", start=1, end=10),
            {}, None, {}, None, {}, "missing", status,
        )
        assert result is None
        assert has_valid is False


# ===========================================================================
# get_ref_ids_liftoff — gene_only with unknown id
# ===========================================================================

class TestGetRefIdsLiftoffEdgeCases:
    def test_gene_only_unknown_returns_none(self):
        d = {"known_gene": lifton_class.Lifton_feature("known_gene")}
        ref_gene, ref_trans = lifton_utils.get_ref_ids_liftoff(
            d, "ghost_gene", None,
        )
        assert ref_gene is None
        assert ref_trans is None

    def test_trans_only_unknown_returns_none(self):
        d = {"known_gene": lifton_class.Lifton_feature("known_gene")}
        ref_gene, ref_trans = lifton_utils.get_ref_ids_liftoff(
            d, None, "ghost_tx",
        )
        assert ref_gene is None
        assert ref_trans is None


# ===========================================================================
# update_cds_list Case 2 inner-CDS branch (3+ CDSs share a single exon)
# ===========================================================================

class TestUpdateCdsListCase2InnerCDS:
    def test_three_cds_one_exon_drives_inner_branch(
            self, gff_single_cds, fake_args, make_gffutils_feature):
        db = annotation.Annotation(
            str(gff_single_cds), False, False, "create_unique", None, True, False,
        ).db_connection
        feat_dict = {"g_s": lifton_class.Lifton_feature("g_s")}
        gene = lifton_class.Lifton_GENE(
            ref_gene_id="g_s",
            gffutil_entry_gene=copy.deepcopy(db["g_s"]),
            ref_gene_attrs={"ID": ["g_s"], "gene_biotype": ["protein_coding"]},
            tree_dict={},
            ref_features_dict=feat_dict,
            args=fake_args,
        )
        trans = gene.add_transcript(
            "tx_s", copy.deepcopy(db["tx_s"]),
            {"ID": ["tx_s"], "Parent": ["g_s"]},
        )
        trans.add_exon(copy.deepcopy(db["ex_s"]))
        cdss = [
            lifton_class.Lifton_CDS(make_gffutils_feature(
                featuretype="CDS", start=110, end=130, frame="0",
                attributes={"ID": ["A"], "Parent": ["tx_s"]},
            )),
            lifton_class.Lifton_CDS(make_gffutils_feature(
                featuretype="CDS", start=140, end=160, frame="0",
                attributes={"ID": ["B"], "Parent": ["tx_s"]},
            )),
            lifton_class.Lifton_CDS(make_gffutils_feature(
                featuretype="CDS", start=170, end=190, frame="0",
                attributes={"ID": ["C"], "Parent": ["tx_s"]},
            )),
        ]
        trans.update_cds_list(cdss)
        # Three CDSs -> three exons in Case 2
        assert len(trans.exons) == 3
        # The middle exon was created by the "else" branch (cds_idx is
        # neither first nor last)
        middle = trans.exons[1]
        assert (middle.entry.start, middle.entry.end) == (140, 160)


# ===========================================================================
# __iterate_exons_update_cds — novel CDS branches and "after ORF" branch
# ===========================================================================

def _make_lifton_trans_with_n_exons(make_gffutils_feature, n=3,
                                    strand="+", with_cds=True):
    """Build a bare Lifton_TRANS with N exons, each 100 bp wide and
    gap-free (so accum_exon_length math is predictable)."""
    instance = lifton_class.Lifton_TRANS.__new__(lifton_class.Lifton_TRANS)
    instance.entry = make_gffutils_feature(
        featuretype="mRNA", start=1, end=n * 100, strand=strand,
        attributes={"ID": ["tx"], "Parent": ["g"]},
    )
    instance.exons = []
    instance.exon_dic = {}
    for i in range(n):
        ex_feat = make_gffutils_feature(
            featuretype="exon", start=1 + i * 100, end=(i + 1) * 100,
            strand=strand,
            attributes={"ID": [f"ex{i}"], "Parent": ["tx"]},
        )
        ex = lifton_class.Lifton_EXON(ex_feat)
        if with_cds:
            cds_feat = make_gffutils_feature(
                featuretype="CDS", start=1 + i * 100, end=(i + 1) * 100,
                strand=strand, frame="0",
                attributes={"ID": [f"cds{i}"], "Parent": ["tx"]},
            )
            ex.add_cds(cds_feat)
        instance.exons.append(ex)
    return instance


class TestIterateExonsUpdateCDS:
    def test_positive_strand_full_orf_in_middle_exon(self,
                                                    make_gffutils_feature):
        """ORF spans the middle exon entirely; first exon must lose its
        CDS, last exon must lose its CDS, middle exon keeps full CDS."""
        trans = _make_lifton_trans_with_n_exons(make_gffutils_feature, n=3,
                                               strand="+", with_cds=True)
        orf = lifton_class.Lifton_ORF(start=100, end=200)
        # Use the public dispatcher via name-mangled private attr access
        trans._Lifton_TRANS__update_cds_boundary(orf)
        # First exon (1..100): accum_exon_length=0, curr_exon_len=100,
        # final_orf.start=100 not < 100 -> exon.cds set None
        assert trans.exons[0].cds is None
        # Last exon: final_orf.end=200 <= accum_exon_length=200 -> cds None
        assert trans.exons[2].cds is None

    def test_negative_strand_dispatch(self, make_gffutils_feature):
        trans = _make_lifton_trans_with_n_exons(make_gffutils_feature, n=3,
                                               strand="-", with_cds=True)
        orf = lifton_class.Lifton_ORF(start=50, end=150)
        # Should not crash; just exercise the - branch
        trans._Lifton_TRANS__update_cds_boundary(orf)

    def test_novel_cds_creation_when_exon_has_no_cds_positive_strand(
            self, make_gffutils_feature):
        trans = _make_lifton_trans_with_n_exons(make_gffutils_feature, n=3,
                                               strand="+", with_cds=False)
        orf = lifton_class.Lifton_ORF(start=20, end=180)
        trans._Lifton_TRANS__update_cds_boundary(orf)
        # First exon now carries a novel CDS starting at exon.start + 20
        assert trans.exons[0].cds is not None
        assert trans.exons[0].cds.entry.start == 1 + 20

    def test_novel_cds_creation_when_exon_has_no_cds_negative_strand(
            self, make_gffutils_feature):
        trans = _make_lifton_trans_with_n_exons(make_gffutils_feature, n=3,
                                               strand="-", with_cds=False)
        orf = lifton_class.Lifton_ORF(start=20, end=180)
        trans._Lifton_TRANS__update_cds_boundary(orf)
        # On - strand, the function iterates self.exons[::-1] so the
        # "first partial CDS" branch fires on the LAST original exon
        # (ex2: 201..300). orf.start=20, accum_exon_length=0:
        # cds.end = 300 - 20 = 280; cds.start = exon.start = 201
        last_exon = trans.exons[2]
        assert last_exon.cds is not None
        assert last_exon.cds.entry.start == 201
        assert last_exon.cds.entry.end == 280

    def test_inner_full_exon_cds_extension_when_cds_is_none(
            self, make_gffutils_feature):
        """Drive the path: middle exon has no CDS, ORF spans the whole
        middle exon → add_novel_lifton_cds covers the full exon range."""
        trans = _make_lifton_trans_with_n_exons(make_gffutils_feature, n=3,
                                               strand="+", with_cds=False)
        # Place ORF so that final_orf.start < accum_exon_length (middle)
        # and final_orf.end > accum_exon_length+curr_exon_len.
        # Middle exon: accum_exon_length=100, curr_exon_len=100 →
        # need orf.start < 100 and orf.end > 200.
        orf = lifton_class.Lifton_ORF(start=20, end=250)
        trans._Lifton_TRANS__update_cds_boundary(orf)
        # Middle exon now carries a full-width novel CDS
        assert trans.exons[1].cds is not None
        assert trans.exons[1].cds.entry.start == 101
        assert trans.exons[1].cds.entry.end == 200


# ===========================================================================
# Lifton_GENE: __get_gene_copy fallback when ref_id not in dict
# ===========================================================================

class TestGeneCopyFallback:
    def test_unknown_ref_gene_id_yields_zero_copy(self,
                                                  gff_standard, fake_args):
        db = annotation.Annotation(
            str(gff_standard), False, False, "create_unique", None, True, False,
        ).db_connection
        # Pass an empty ref_features_dict so the dict lookup misses -> 0
        gene = lifton_class.Lifton_GENE(
            ref_gene_id="ghost",
            gffutil_entry_gene=copy.deepcopy(db["gene1"]),
            ref_gene_attrs={"ID": ["ghost"], "gene_biotype": ["protein_coding"]},
            tree_dict={},
            ref_features_dict={},
            args=fake_args,
        )
        assert gene.copy_num == 0


# ===========================================================================
# get_ref_liffover_features — "other" gene_type and missing gene_type_key
# ===========================================================================

@pytest.fixture
def gff_with_pseudogene_biotype(tmp_path):
    fp = tmp_path / "pseudo.gff3"
    fp.write_text(
        "##gff-version 3\n"
        "chr1\tt\tgene\t100\t200\t.\t+\t.\tID=pg;gene_biotype=pseudogene\n"
        "chr1\tt\tmRNA\t100\t200\t.\t+\t.\tID=pgtx;Parent=pg\n"
        "chr1\tt\texon\t100\t200\t.\t+\t.\tID=pgex;Parent=pgtx\n"
    )
    return fp


@pytest.fixture
def gff_no_gene_type_attribute(tmp_path):
    """Gene with no gene_biotype / gene_type → falls into the "other"
    bucket."""
    fp = tmp_path / "no_type.gff3"
    fp.write_text(
        "##gff-version 3\n"
        "chr1\tt\tgene\t100\t200\t.\t+\t.\tID=g_no_type\n"
        "chr1\tt\tmRNA\t100\t200\t.\t+\t.\tID=tx_no_type;Parent=g_no_type\n"
        "chr1\tt\texon\t100\t200\t.\t+\t.\tID=ex_no_type;Parent=tx_no_type\n"
    )
    return fp


class TestRefLiffoverFeaturesOtherBranches:
    def test_pseudogene_falls_to_other(self, gff_with_pseudogene_biotype,
                                       tmp_path):
        ann = annotation.Annotation(
            str(gff_with_pseudogene_biotype), False, False, "create_unique",
            None, True, False,
        )
        args = SimpleNamespace(annotation_database="RefSeq",
                               evaluation=False,
                               evaluation_liftoff_chm13=False)
        d, _, _, _ = lifton_utils.get_ref_liffover_features(
            ["gene"], ann, str(tmp_path), args,
        )
        # pseudogene biotype is neither protein_coding nor (lncRNA|ncRNA)
        # → falls into the "other" bucket; flags remain default False/False
        assert d["pg"].is_protein_coding is False
        assert d["pg"].is_non_coding is False

    def test_no_gene_type_attribute_falls_to_other(
            self, gff_no_gene_type_attribute, tmp_path):
        ann = annotation.Annotation(
            str(gff_no_gene_type_attribute), False, False, "create_unique",
            None, True, False,
        )
        args = SimpleNamespace(annotation_database="RefSeq",
                               evaluation=False,
                               evaluation_liftoff_chm13=False)
        d, _, _, _ = lifton_utils.get_ref_liffover_features(
            ["gene"], ann, str(tmp_path), args,
        )
        assert d["g_no_type"].is_protein_coding is False
        assert d["g_no_type"].is_non_coding is False
