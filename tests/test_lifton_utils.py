"""Unit tests for lifton.lifton_utils — pure helpers only.

Tests target the legacy semantics. If any assertion changes during
later refactors, the diff itself is the audit trail.
"""

from __future__ import annotations

import os
from types import SimpleNamespace

import pytest

from lifton import lifton_utils


# ---------------------------------------------------------------------------
# segments_overlap_length
# ---------------------------------------------------------------------------

class TestSegmentsOverlapLength:
    def test_clear_overlap(self):
        ovp_len, ovp = lifton_utils.segments_overlap_length((10, 20), (15, 25))
        assert ovp is True
        assert ovp_len == 6  # 20 - 15 + 1

    def test_touching_endpoints_count_as_overlap(self):
        ovp_len, ovp = lifton_utils.segments_overlap_length((10, 20), (20, 30))
        assert ovp is True
        assert ovp_len == 1

    def test_disjoint_segments(self):
        ovp_len, ovp = lifton_utils.segments_overlap_length((1, 5), (10, 15))
        assert ovp is False
        assert ovp_len <= 0

    def test_full_containment(self):
        # Phase 5 bug fix #4: overlap length is now the true intersection
        # width (min(end) - max(start) + 1), not the legacy quirk.
        # (20..30) is fully inside (10..50) → 11 bp overlap.
        ovp_len, ovp = lifton_utils.segments_overlap_length((10, 50), (20, 30))
        assert ovp is True
        assert ovp_len == 11

    def test_order_independence(self):
        a = lifton_utils.segments_overlap_length((30, 40), (35, 45))
        b = lifton_utils.segments_overlap_length((35, 45), (30, 40))
        assert a == b

    def test_invalid_tuple_raises(self):
        with pytest.raises(ValueError):
            lifton_utils.segments_overlap_length((1, 2, 3), (4, 5))


# ---------------------------------------------------------------------------
# custom_bisect_insert  (sorts by .entry.end)
# ---------------------------------------------------------------------------

class _BiscectStub:
    def __init__(self, end):
        self.entry = SimpleNamespace(end=end)
    def __repr__(self):
        return f"<S end={self.entry.end}>"


class TestCustomBisectInsert:
    def test_inserts_in_ascending_end_order(self):
        bucket = []
        for e in (50, 10, 30, 40, 20):
            lifton_utils.custom_bisect_insert(bucket, _BiscectStub(e))
        assert [x.entry.end for x in bucket] == [10, 20, 30, 40, 50]

    def test_insert_into_empty_list(self):
        bucket = []
        lifton_utils.custom_bisect_insert(bucket, _BiscectStub(7))
        assert len(bucket) == 1 and bucket[0].entry.end == 7

    def test_duplicates_are_kept(self):
        bucket = []
        for e in (10, 10, 10):
            lifton_utils.custom_bisect_insert(bucket, _BiscectStub(e))
        assert [x.entry.end for x in bucket] == [10, 10, 10]


# ---------------------------------------------------------------------------
# get_ID_base / get_ID
# ---------------------------------------------------------------------------

class TestGetIDBase:
    """Post-merge note: `get_ID_base` now accepts an optional
    `ref_features_dict`. Without it the function is conservative and
    refuses to strip any trailing _<int> suffix (to avoid mangling ids
    like FMUND_1 where the trailing integer is intrinsic). With the
    dict, it strips only when the stripped base appears in the dict."""

    def test_strips_trailing_integer_suffix_with_dict(self):
        d = {"ENST00000123": object()}
        assert lifton_utils.get_ID_base("ENST00000123_2", d) == "ENST00000123"

    def test_no_dict_keeps_suffix_conservatively(self):
        assert lifton_utils.get_ID_base("ENST00000123_2") == "ENST00000123_2"

    def test_no_suffix_returns_original(self):
        assert lifton_utils.get_ID_base("ENST00000123") == "ENST00000123"

    def test_non_integer_suffix_kept(self):
        assert lifton_utils.get_ID_base("ENST00000123_alt") == "ENST00000123_alt"

    def test_multiple_underscores_strips_only_last_int_with_dict(self):
        d = {"rna_NM_001": object()}
        assert lifton_utils.get_ID_base("rna_NM_001_3", d) == "rna_NM_001"

    def test_dict_lookup_miss_keeps_id(self):
        # base "FMUND" not in dict → don't strip
        assert lifton_utils.get_ID_base("FMUND_1", {"otherbase": None}) == "FMUND_1"

    def test_get_ID_pair_returns_conservative_base(self):
        feat = SimpleNamespace(id="rna_NM_001_3")
        full, base = lifton_utils.get_ID(feat)
        # get_ID calls get_ID_base without the dict → conservative no-strip
        assert full == "rna_NM_001_3"
        assert base == "rna_NM_001_3"


# ---------------------------------------------------------------------------
# check_protein_valid + get_truncated_protein
# ---------------------------------------------------------------------------

class TestProteinValidity:
    def test_valid_simple(self):
        assert lifton_utils.check_protein_valid("MAGT*") is True

    def test_empty_invalid(self):
        assert lifton_utils.check_protein_valid("") is False

    def test_no_start_methionine(self):
        assert lifton_utils.check_protein_valid("AAGT*") is False

    def test_no_stop(self):
        assert lifton_utils.check_protein_valid("MAGT") is False

    def test_internal_stop_invalid(self):
        # two * means internal stop
        assert lifton_utils.check_protein_valid("MA*GT*") is False


class TestGetTruncatedProtein:
    def test_filters_only_invalid(self):
        proteins = {
            "good": "MAGT*",
            "bad_no_stop": "MAGT",
            "bad_no_start": "AAGT*",
        }
        truncated = lifton_utils.get_truncated_protein(proteins)
        assert set(truncated.keys()) == {"bad_no_stop", "bad_no_start"}

    def test_empty_input(self):
        assert lifton_utils.get_truncated_protein({}) == {}


# ---------------------------------------------------------------------------
# write_seq_2_file
# ---------------------------------------------------------------------------

class TestWriteSeq2File:
    @pytest.mark.parametrize("kind,suffix", [
        ("proteins", "proteins.fa"),
        ("transcripts", "transcripts.fa"),
        ("truncated_proteins", "proteins_truncated.fa"),
    ])
    def test_writes_fasta_with_expected_suffix(self, tmp_path, kind, suffix):
        out = lifton_utils.write_seq_2_file(
            str(tmp_path), {"a": "MAGT*", "b": "MQQQ*"}, kind
        )
        assert out.endswith(suffix)
        body = open(out).read()
        assert ">a\nMAGT*\n" in body
        assert ">b\nMQQQ*\n" in body


# ---------------------------------------------------------------------------
# get_parent_features_to_lift
# ---------------------------------------------------------------------------

class TestGetParentFeaturesToLift:
    def test_default_when_no_file(self):
        assert lifton_utils.get_parent_features_to_lift(None) == ["gene"]

    def test_reads_one_per_line(self, tmp_path):
        f = tmp_path / "f.txt"
        f.write_text("gene\npseudogene\nlnc_RNA\n")
        assert lifton_utils.get_parent_features_to_lift(str(f)) == [
            "gene", "pseudogene", "lnc_RNA",
        ]


# ---------------------------------------------------------------------------
# get_gene_like_feature_types  (Iteration 5: --lift-gene-like auto-detect)
# ---------------------------------------------------------------------------

class TestGetGeneLikeFeatureTypes:
    """Auto-detect the top-level parent types that have a transcript/exon
    hierarchy, so the lift can go beyond the hardcoded `["gene"]`."""

    @staticmethod
    def _db(tmp_path, gff_text):
        from lifton import annotation
        fp = tmp_path / "ref.gff3"
        fp.write_text(gff_text)
        return annotation.Annotation(str(fp), None, None, "create_unique",
                                     "ID", False, False, True)

    def test_detects_gene_and_mixed_pseudogene_excludes_childless(self, tmp_path):
        # gene (has mRNA->exon), pseudogene #1 childless, pseudogene #2 with an
        # exon child, enhancer childless. Expect gene-like = gene + pseudogene;
        # enhancer excluded; child types (mRNA/exon) excluded (they carry Parent).
        gff = (
            "##gff-version 3\n"
            "chr1\ts\tgene\t100\t200\t.\t+\t.\tID=gene1;gene_biotype=protein_coding\n"
            "chr1\ts\tmRNA\t100\t200\t.\t+\t.\tID=rna1;Parent=gene1\n"
            "chr1\ts\texon\t100\t200\t.\t+\t.\tID=ex1;Parent=rna1\n"
            "chr1\ts\tpseudogene\t300\t400\t.\t+\t.\tID=pg_childless;gene_biotype=pseudogene\n"
            "chr1\ts\tpseudogene\t500\t600\t.\t+\t.\tID=pg_child;gene_biotype=pseudogene\n"
            "chr1\ts\texon\t500\t600\t.\t+\t.\tID=pgex;Parent=pg_child\n"
            "chr1\ts\tenhancer\t700\t800\t.\t+\t.\tID=enh1\n"
        )
        got = lifton_utils.get_gene_like_feature_types(self._db(tmp_path, gff))
        assert got == ["gene", "pseudogene"]
        assert "enhancer" not in got
        assert "exon" not in got and "mRNA" not in got

    def test_falls_back_to_gene_when_nothing_gene_like(self, tmp_path):
        # only childless meta/regulatory features -> fall back to ["gene"]
        gff = (
            "##gff-version 3\n"
            "chr1\ts\tenhancer\t10\t20\t.\t+\t.\tID=e1\n"
            "chr1\ts\tregion\t1\t1000\t.\t+\t.\tID=r1\n"
        )
        assert lifton_utils.get_gene_like_feature_types(self._db(tmp_path, gff)) == ["gene"]


# ---------------------------------------------------------------------------
# check_ovps_ratio
# ---------------------------------------------------------------------------

class TestCheckOvpsRatio:
    def test_unknown_chromosome_returns_false(self):
        mtrans = SimpleNamespace(seqid="chrUNKNOWN")
        assert lifton_utils.check_ovps_ratio(
            mtrans, (10, 100), 0.1, {"chr1": None}
        ) is False

    def test_overlap_ratio_triggers_true(self):
        """Phase 5 bug fix #3 verified: check_ovps_ratio now unpacks a
        tuple (or Interval) before calling IntervalTree.overlap()."""
        from intervaltree import Interval, IntervalTree
        tree = IntervalTree()
        tree.add(Interval(50, 150, "ref_gene"))
        tree_dict = {"chr1": tree}
        mtrans = SimpleNamespace(seqid="chr1")
        # tuple input
        assert lifton_utils.check_ovps_ratio(
            mtrans, (40, 200), 0.1, tree_dict
        ) is True
        # Interval input also works
        assert lifton_utils.check_ovps_ratio(
            mtrans, Interval(40, 200, "tx"), 0.1, tree_dict
        ) is True
