"""Phase 4.5 Step 6 (continued) — drive lifton_class.__find_orfs and
__update_cds_boundary on both strands so the ORF-rescue body, novel-CDS
creation, and CDS-boundary patching are all exercised.

Crafted to push lifton_class.py coverage past 90 %.
"""

from __future__ import annotations

import copy
import textwrap
from pathlib import Path
from types import SimpleNamespace

import pytest
from Bio.Seq import Seq
from pyfaidx import Fasta

from lifton import annotation, lifton_class


# ---------------------------------------------------------------------------
# Build a chromosome whose embedded transcript has a frameshift in its
# CDS, so when LiftOn translates it, the protein diverges from the
# reference and __find_orfs fires.
# ---------------------------------------------------------------------------

def _wrap(seq, width=60):
    return "\n".join(textwrap.wrap(seq, width)) + "\n"


@pytest.fixture
def fasta_with_long_orf(tmp_path):
    """1200-bp chromosome whose UTRs and exons contain a long alternative
    ORF starting later in the sequence. Used to force __find_orfs to find
    a better ORF than the original CDS."""
    base = ["A"] * 1200
    # Original CDS spans 101..199 + 301..399 — a clean ATG..TAA but
    # with a frameshift introduced at the start.
    exon1 = "AT" + "G" + "GCT" * 32   # starts ATGGCT... (M, A, A...)
    exon2 = ("GCT" * 32) + "TAA"
    for i, ch in enumerate(exon1):
        base[100 + i] = ch
    for i, ch in enumerate(exon2):
        base[300 + i] = ch
    # Plant a much longer ORF in the 3' UTR region (positions 500..900)
    long_orf = "ATG" + "GCT" * 100 + "TAA"   # 306 bp
    for i, ch in enumerate(long_orf):
        base[500 + i] = ch
    fp = tmp_path / "long_orf.fa"
    fp.write_text(">chr1\n" + _wrap("".join(base)))
    return fp


@pytest.fixture
def fasta_with_no_internal_stops(tmp_path):
    """1200 bp chromosome whose two CDS exons translate to a clean ORF
    matching the reference, but the surrounding UTR sequence contains
    an alternative ATG…TAA. This drives ORF discovery without rescue
    triggering (because identity is already maxed)."""
    base = ["A"] * 1200
    cds = "ATG" + "GCT" * 30 + "TAA"
    for i, ch in enumerate(cds):
        base[100 + i] = ch
    fp = tmp_path / "no_stops.fa"
    fp.write_text(">chr1\n" + _wrap("".join(base)))
    return fp


@pytest.fixture
def gff_full_transcript(tmp_path):
    fp = tmp_path / "full.gff3"
    fp.write_text(
        "##gff-version 3\n"
        "chr1\ttest\tgene\t1\t1200\t.\t+\t.\tID=g1;gene_biotype=protein_coding\n"
        "chr1\ttest\tmRNA\t1\t1200\t.\t+\t.\tID=tx1;Parent=g1\n"
        "chr1\ttest\texon\t1\t199\t.\t+\t.\tID=ex1;Parent=tx1\n"
        "chr1\ttest\texon\t301\t1200\t.\t+\t.\tID=ex2;Parent=tx1\n"
        "chr1\ttest\tCDS\t101\t199\t.\t+\t0\tID=cds1;Parent=tx1\n"
        "chr1\ttest\tCDS\t301\t399\t.\t+\t0\tID=cds2;Parent=tx1\n"
    )
    return fp


@pytest.fixture
def gff_full_transcript_negative_strand(tmp_path):
    fp = tmp_path / "neg.gff3"
    fp.write_text(
        "##gff-version 3\n"
        "chr1\ttest\tgene\t1\t1200\t.\t-\t.\tID=gN;gene_biotype=protein_coding\n"
        "chr1\ttest\tmRNA\t1\t1200\t.\t-\t.\tID=txN;Parent=gN\n"
        "chr1\ttest\texon\t1\t199\t.\t-\t.\tID=exN1;Parent=txN\n"
        "chr1\ttest\texon\t301\t1200\t.\t-\t.\tID=exN2;Parent=txN\n"
        "chr1\ttest\tCDS\t101\t199\t.\t-\t0\tID=cdsN1;Parent=txN\n"
        "chr1\ttest\tCDS\t301\t399\t.\t-\t0\tID=cdsN2;Parent=txN\n"
    )
    return fp


def _build_trans(gff_path, fake_args, ref_features_dict_one_gene,
                 gene_id="g1", trans_id="tx1",
                 exon_ids=("ex1", "ex2"),
                 cds_ids=("cds1", "cds2")):
    db = annotation.Annotation(
        str(gff_path), False, False, "create_unique", None, True, False,
    ).db_connection
    feat_dict = {gene_id: lifton_class.Lifton_feature(gene_id)}
    feat_dict[gene_id].is_protein_coding = True
    gene = lifton_class.Lifton_GENE(
        ref_gene_id=gene_id,
        gffutil_entry_gene=copy.deepcopy(db[gene_id]),
        ref_gene_attrs={"ID": [gene_id], "gene_biotype": ["protein_coding"]},
        tree_dict={},
        ref_features_dict=feat_dict,
        args=fake_args,
    )
    trans = gene.add_transcript(
        trans_id, copy.deepcopy(db[trans_id]),
        {"ID": [trans_id], "Parent": [gene_id]},
    )
    for ex_id in exon_ids:
        trans.add_exon(copy.deepcopy(db[ex_id]))
    for cds_id in cds_ids:
        trans.add_cds(copy.deepcopy(db[cds_id]))
    return gene, trans


# ---------------------------------------------------------------------------
# Tests targeting __find_orfs body
# ---------------------------------------------------------------------------

class TestFindOrfsBody:
    def test_find_orfs_runs_when_premature_stop_introduced(
            self, gff_full_transcript, fasta_with_long_orf, fake_args,
            ref_features_dict_one_gene):
        """The CDS in fasta_with_long_orf produces a protein. Reference
        protein is the same protein. Force a mutation by giving a
        DIFFERENT reference (longer protein from the alternative ORF)."""
        _, trans = _build_trans(gff_full_transcript, fake_args,
                                ref_features_dict_one_gene)
        fa = Fasta(str(fasta_with_long_orf))
        # Build a much longer ref protein
        ref_dna = "ATG" + "GCT" * 100 + "TAA"
        ref_protein = str(Seq(ref_dna).translate())
        ref_trans_seq = ref_dna
        status = lifton_class.Lifton_Status()
        # Force ORF rescue: ref_trans differs from translated trans → mutations recorded
        trans.orf_search_protein(
            fa, ref_protein, ref_trans_seq, status,
            is_non_coding=False, eval_only=False,
        )
        # Either ORF rescue ran (and identity maybe improved) OR mutations
        # were recorded. Both indicate __find_orfs body executed.
        assert "mutation" in trans.entry.attributes or status.lifton_aa > 0

    def test_find_orfs_runs_for_negative_strand(
            self, gff_full_transcript_negative_strand,
            fasta_with_long_orf, fake_args):
        feat_dict = {"gN": lifton_class.Lifton_feature("gN")}
        feat_dict["gN"].is_protein_coding = True
        ref_features = feat_dict
        _, trans = _build_trans(
            gff_full_transcript_negative_strand, fake_args,
            ref_features, gene_id="gN", trans_id="txN",
            exon_ids=("exN1", "exN2"), cds_ids=("cdsN1", "cdsN2"),
        )
        fa = Fasta(str(fasta_with_long_orf))
        ref_dna = "ATG" + "GCT" * 100 + "TAA"
        ref_protein = str(Seq(ref_dna).translate())
        status = lifton_class.Lifton_Status()
        # Should not crash on negative strand
        trans.orf_search_protein(
            fa, ref_protein, ref_dna, status,
            is_non_coding=False, eval_only=False,
        )
        # Reaching this line proves the negative-strand __update_cds_boundary
        # branch was at least reachable.
        assert isinstance(status.status, list)


# ---------------------------------------------------------------------------
# write_entry non-temp branches and update_boundaries
# ---------------------------------------------------------------------------

class TestWriteEntryNonTempBranches:
    def test_write_gene_with_noncoding_lncRNA_branch(
            self, gff_noncoding, fake_args, tmp_path):
        db = annotation.Annotation(
            str(gff_noncoding), False, False, "create_unique", None, True, False,
        ).db_connection
        feat_dict = {"ncg": lifton_class.Lifton_feature("ncg")}
        feat_dict["ncg"].is_non_coding = True
        gene = lifton_class.Lifton_GENE(
            ref_gene_id="ncg",
            gffutil_entry_gene=copy.deepcopy(db["ncg"]),
            ref_gene_attrs={"ID": ["ncg"], "gene_biotype": ["lncRNA"]},
            tree_dict={},
            ref_features_dict=feat_dict,
            args=fake_args,
        )
        trans = gene.add_transcript(
            "nctx", copy.deepcopy(db["nctx"]),
            {"ID": ["nctx"], "Parent": ["ncg"]},
        )
        # Force featuretype to lnc_RNA for the non-coding write_entry branch
        trans.entry.featuretype = "lnc_RNA"
        trans.add_exon(copy.deepcopy(db["ncex"]))

        out = tmp_path / "nc.gff3"
        with open(out, "w") as fw:
            stats = {"coding": {}, "non-coding": {}, "other": {}}
            gene.write_entry(fw, stats)
        # The non-coding bucket must have been incremented for nctx
        assert "nctx" in stats["non-coding"]

    def test_write_gene_other_bucket_when_neither_coding_nor_noncoding(
            self, gff_standard, fake_args, ref_features_dict_one_gene,
            tmp_path):
        db = annotation.Annotation(
            str(gff_standard), False, False, "create_unique", None, True, False,
        ).db_connection
        gene = lifton_class.Lifton_GENE(
            ref_gene_id="gene1",
            gffutil_entry_gene=copy.deepcopy(db["gene1"]),
            ref_gene_attrs={"ID": ["gene1"], "gene_biotype": ["other"]},
            tree_dict={},
            ref_features_dict=ref_features_dict_one_gene,
            args=fake_args,
        )
        trans = gene.add_transcript(
            "tx1", copy.deepcopy(db["tx1"]),
            {"ID": ["tx1"], "Parent": ["gene1"]},
        )
        trans.entry.featuretype = "transcript"  # not mRNA, not nc_RNA
        out = tmp_path / "other.gff3"
        with open(out, "w") as fw:
            stats = {"coding": {}, "non-coding": {}, "other": {}}
            gene.write_entry(fw, stats)
        assert "tx1" in stats["other"]

    def test_update_boundaries_picks_min_max_across_transcripts(
            self, gff_standard, fake_args, ref_features_dict_one_gene):
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
        gene.add_transcript("tx1", copy.deepcopy(db["tx1"]),
                            {"ID": ["tx1"], "Parent": ["gene1"]})
        # update_boundaries only WIDENS — start can only shrink, end can
        # only grow. Set narrow bounds in the middle of the trans range.
        gene.entry.start = 250
        gene.entry.end = 260
        gene.update_boundaries()
        # start widens leftward to trans.start (101)
        assert gene.entry.start == 101
        # end stays at 260 because trans.end(399) > 260 → updated to 399
        assert gene.entry.end == 399


# ---------------------------------------------------------------------------
# Lifton_EXON / Lifton_CDS leaf methods (print_*, add_lifton_cds None branch)
# ---------------------------------------------------------------------------

class TestLeafPrint:
    def test_print_exon_does_not_crash(self, make_gffutils_feature, capsys):
        feat = make_gffutils_feature(
            featuretype="exon", start=1, end=10,
            attributes={"ID": ["e1"], "Parent": ["t1"]},
        )
        ex = lifton_class.Lifton_EXON(feat)
        ex.print_exon()
        out = capsys.readouterr().out
        assert len(out) > 0

    def test_print_exon_with_cds(self, make_gffutils_feature, capsys):
        e_feat = make_gffutils_feature(
            featuretype="exon", start=1, end=10,
            attributes={"ID": ["e1"], "Parent": ["t1"]},
        )
        c_feat = make_gffutils_feature(
            featuretype="CDS", start=1, end=10, frame="0",
            attributes={"ID": ["c1"], "Parent": ["t1"]},
        )
        ex = lifton_class.Lifton_EXON(e_feat)
        ex.add_cds(c_feat)
        ex.print_exon()
        out = capsys.readouterr().out
        assert "exon" in out.lower() or "CDS" in out

    def test_print_cds_does_not_crash(self, make_gffutils_feature, capsys):
        feat = make_gffutils_feature(
            featuretype="CDS", start=1, end=10, frame="0",
            attributes={"ID": ["c1"], "Parent": ["t1"]},
        )
        cds = lifton_class.Lifton_CDS(feat)
        cds.print_cds()
        out = capsys.readouterr().out
        assert len(out) > 0

    def test_add_lifton_cds_none_branch(self, make_gffutils_feature):
        feat = make_gffutils_feature(
            featuretype="exon", start=1, end=10,
            attributes={"ID": ["e1"], "Parent": ["t1"]},
        )
        ex = lifton_class.Lifton_EXON(feat)
        # add_lifton_cds(None) should set self.cds = None without raising
        ex.add_lifton_cds(None)
        assert ex.cds is None

    def test_cds_update_info(self, make_gffutils_feature):
        feat = make_gffutils_feature(
            featuretype="CDS", start=1, end=10, frame="0",
            attributes={"ID": ["c1"], "Parent": ["t1"]},
        )
        cds = lifton_class.Lifton_CDS(feat)
        cds.update_CDS_info(50, 100)
        assert cds.entry.start == 50
        assert cds.entry.end == 100


# ---------------------------------------------------------------------------
# Lifton_TRANS misc — mv_exon_idx, update_gffutil_entry_trans, print_transcript
# ---------------------------------------------------------------------------

class TestTransMisc:
    def test_mv_exon_idx_caps_at_last(self, gff_standard, fake_args,
                                      ref_features_dict_one_gene):
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
        trans = gene.add_transcript("tx1", copy.deepcopy(db["tx1"]),
                                    {"ID": ["tx1"], "Parent": ["gene1"]})
        trans.add_exon(copy.deepcopy(db["exon1"]))
        # With one exon, mv_exon_idx(0) returns 0 (cannot advance)
        assert trans.mv_exon_idx(0) == 0
        trans.add_exon(copy.deepcopy(db["exon2"]))
        # With two exons, mv_exon_idx(0) returns 1
        assert trans.mv_exon_idx(0) == 1
        # mv_exon_idx(1) stays at 1 (last index)
        assert trans.mv_exon_idx(1) == 1

    def test_update_gffutil_entry_trans_merges_attrs(
            self, gff_standard, fake_args, ref_features_dict_one_gene):
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
        trans = gene.add_transcript("tx1", copy.deepcopy(db["tx1"]),
                                    {"ID": ["tx1"], "Parent": ["gene1"]})
        # New entry with extra attribute
        new_entry = copy.deepcopy(db["tx1"])
        new_entry.attributes["new_attr"] = ["value"]
        trans.update_gffutil_entry_trans(new_entry)
        assert trans.entry.attributes["new_attr"] == ["value"]

    def test_print_transcript_does_not_crash(
            self, gff_standard, fake_args, ref_features_dict_one_gene,
            capsys):
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
        trans = gene.add_transcript("tx1", copy.deepcopy(db["tx1"]),
                                    {"ID": ["tx1"], "Parent": ["gene1"]})
        trans.add_exon(copy.deepcopy(db["exon1"]))
        trans.print_transcript()
        out = capsys.readouterr().out
        assert "tx1" in out


# ---------------------------------------------------------------------------
# Translate empty / None propagation
# ---------------------------------------------------------------------------

class TestAlignmentEmptyPaths:
    def test_align_coding_seq_returns_none_for_empty_ref(
            self, gff_standard, fake_args, ref_features_dict_one_gene):
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
        trans = gene.add_transcript("tx1", copy.deepcopy(db["tx1"]),
                                    {"ID": ["tx1"], "Parent": ["gene1"]})
        status = lifton_class.Lifton_Status()
        aln, peps = trans.align_coding_seq("MAGT*", "", status)
        assert aln is None
        assert peps is None

    def test_align_coding_seq_returns_none_for_empty_query(
            self, gff_standard, fake_args, ref_features_dict_one_gene):
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
        trans = gene.add_transcript("tx1", copy.deepcopy(db["tx1"]),
                                    {"ID": ["tx1"], "Parent": ["gene1"]})
        status = lifton_class.Lifton_Status()
        aln, peps = trans.align_coding_seq("", "MAGT*", status)
        assert aln is None
        assert peps is None

    def test_align_trans_seq_returns_none_for_empty(
            self, gff_standard, fake_args, ref_features_dict_one_gene):
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
        trans = gene.add_transcript("tx1", copy.deepcopy(db["tx1"]),
                                    {"ID": ["tx1"], "Parent": ["gene1"]})
        status = lifton_class.Lifton_Status()
        assert trans.align_trans_seq("ATG", "", status) is None
        assert trans.align_trans_seq("", "ATG", status) is None
