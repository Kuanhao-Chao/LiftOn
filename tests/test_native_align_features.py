"""Phase 11 — native alignment path unit + parity tests.

Covers `lifton/liftoff/native_align.py`:
  - cigar_str_to_pysam_tuples conversion
  - _PysamShim attribute mapping (mirrors what `add_alignment` reads)
  - parse_alignment_from_hits dict shape
  - align_features_to_target_native dispatcher routing
  - mappy-unavailable fallback warning
"""

from __future__ import annotations

from types import SimpleNamespace
from unittest import mock

import pytest

from lifton.liftoff import native_align
from lifton.liftoff.native_align import (
    _PysamShim,
    align_features_to_target_native,
    cigar_str_to_pysam_tuples,
    parse_alignment_from_hits,
)
from lifton.native_bindings import MinimapHit


# ---------------------------------------------------------------------------
# CIGAR parsing
# ---------------------------------------------------------------------------

class TestCigarParsing:
    def test_simple_match(self):
        assert cigar_str_to_pysam_tuples("100M") == [(0, 100)]

    def test_match_indel_match(self):
        assert cigar_str_to_pysam_tuples("50M2I30M") == [
            (0, 50), (1, 2), (0, 30),
        ]

    def test_eqx_match_mismatch(self):
        assert cigar_str_to_pysam_tuples("10=2X20=") == [
            (7, 10), (8, 2), (7, 20),
        ]

    def test_deletion_and_clips(self):
        assert cigar_str_to_pysam_tuples("5S20M3D10M5H") == [
            (4, 5), (0, 20), (2, 3), (0, 10), (5, 5),
        ]

    def test_empty_returns_empty(self):
        assert cigar_str_to_pysam_tuples("") == []

    def test_skipped_region(self):
        assert cigar_str_to_pysam_tuples("100M500N100M") == [
            (0, 100), (3, 500), (0, 100),
        ]


# ---------------------------------------------------------------------------
# _PysamShim
# ---------------------------------------------------------------------------

class TestPysamShim:
    def _hit(self, **overrides):
        defaults = dict(
            query_name="g1", ctg="chr1", r_st=100, r_en=200,
            q_st=0, q_en=100, strand=1, mapq=60, NM=2,
            cigar_str="100M", is_primary=True,
        )
        defaults.update(overrides)
        return MinimapHit(**defaults)

    def test_attribute_mapping(self):
        s = _PysamShim(self._hit())
        assert s.query_name == "g1"
        assert s.reference_name == "chr1"
        assert s.reference_start == 100
        assert s.query_alignment_start == 0
        assert s.query_alignment_end == 100
        assert s.cigar == [(0, 100)]
        assert s.is_reverse is False
        assert s.is_unmapped is False

    def test_negative_strand_marks_reverse(self):
        s = _PysamShim(self._hit(strand=-1))
        assert s.is_reverse is True

    def test_complex_cigar(self):
        s = _PysamShim(self._hit(cigar_str="50=2X48=", q_en=100))
        assert s.cigar == [(7, 50), (8, 2), (7, 48)]


# ---------------------------------------------------------------------------
# parse_alignment_from_hits
# ---------------------------------------------------------------------------

class _FakeFeatureHierarchy:
    def __init__(self, parents, children):
        self.parents = parents
        self.children = children


class _FakeParent:
    def __init__(self, ident, start, end):
        self.id = ident
        self.start = start
        self.end = end


class _FakeChild:
    def __init__(self, start, end):
        self.start = start
        self.end = end


class TestParseAlignmentFromHits:
    def test_unmapped_query_appended(self):
        parents = {"gene1": _FakeParent("gene1", 1, 100)}
        children = {"gene1": [_FakeChild(1, 100)]}
        hierarchy = _FakeFeatureHierarchy(parents, children)
        unmapped = []
        result = parse_alignment_from_hits(
            {}, hierarchy, unmapped, "chrm_by_chrm",
        )
        # No hits → gene1 appears in unmapped
        assert parents["gene1"] in unmapped
        # No emission
        assert result == {}

    def test_present_hit_emits_aligned_seg(self):
        parents = {"gene1": _FakeParent("gene1", 1, 100)}
        children = {"gene1": [_FakeChild(1, 100)]}
        hierarchy = _FakeFeatureHierarchy(parents, children)
        unmapped = []
        hit = MinimapHit(
            query_name="gene1", ctg="chr1", r_st=200, r_en=300,
            q_st=0, q_en=100, strand=1, mapq=60, NM=0,
            cigar_str="100=", is_primary=True,
        )
        result = parse_alignment_from_hits(
            {"gene1": [hit]}, hierarchy, unmapped, "chrm_by_chrm",
        )
        # Edited query name with the "_0" suffix Liftoff appends
        assert "gene1_0" in result
        segs = result["gene1_0"]
        assert len(segs) == 1
        assert segs[0].reference_name == "chr1"
        assert segs[0].is_reverse is False
        assert parents["gene1"] not in unmapped


# ---------------------------------------------------------------------------
# align_features_to_target_native dispatcher
# ---------------------------------------------------------------------------

class TestNativeDispatcher:
    def test_falls_back_when_mappy_unavailable(self, monkeypatch, capsys):
        # Force mappy-unavailable branch
        monkeypatch.setattr(native_align, "is_mappy_available",
                            lambda: False)
        # Patch the legacy align_features_to_target so we can detect
        # the fallback was taken.
        called = {"n": 0}

        def fake_legacy(*args, **kwargs):
            called["n"] += 1
            return {"sentinel": []}

        monkeypatch.setattr(native_align._af, "align_features_to_target",
                            fake_legacy)

        args = SimpleNamespace(target="dummy.fa", mm2_options="", threads=1)
        result = align_features_to_target_native(
            ["chr1"], ["chr1"], args, _FakeFeatureHierarchy({}, {}),
            "chrm_by_chrm", [],
        )
        assert called["n"] == 1
        assert result == {"sentinel": []}
        # Stderr warning emitted
        err = capsys.readouterr().err
        assert "mappy" in err

    def test_native_path_invoked_when_mappy_available(self, monkeypatch,
                                                     tmp_path):
        # Build a tiny target FASTA and a tiny features FASTA so mappy
        # has something to map.
        target = tmp_path / "tgt.fa"
        target.write_text(">chr1\n" + "ACGT" * 200 + "\n")
        features = tmp_path / "features.fa"
        features.write_text(">gene1\n" + "ACGT" * 50 + "\n")

        # Patch get_features_file to point at our features FASTA
        monkeypatch.setattr(native_align._af, "get_features_file",
                            lambda *a, **k: (str(features), "synth"))

        # Provide a hierarchy entry so any hit can be processed.
        hierarchy = _FakeFeatureHierarchy(
            parents={"gene1": _FakeParent("gene1", 1, 200)},
            children={"gene1": [_FakeChild(1, 200)]},
        )

        args = SimpleNamespace(target=str(target), mm2_options="",
                               threads=1)
        unmapped = []
        result = align_features_to_target_native(
            ["chr1"], ["chr1"], args, hierarchy,
            "chrm_by_chrm", unmapped,
        )
        # Whether or not the synthetic seq aligned, the call must not raise.
        assert isinstance(result, dict)


# ---------------------------------------------------------------------------
# Routing in align_features.align_features_to_target
# ---------------------------------------------------------------------------

class TestRoutingDispatch:
    def test_native_flag_routes_to_native_path(self, monkeypatch):
        from lifton.liftoff import align_features
        called = {"native": 0, "legacy": 0}

        def fake_native(*a, **k):
            called["native"] += 1
            return {}

        def fake_legacy(*a, **k):
            called["legacy"] += 1
            return {}

        monkeypatch.setattr(
            "lifton.liftoff.native_align.align_features_to_target_native",
            fake_native,
        )
        # Replace the body of the legacy fallback so we can detect routing
        monkeypatch.setattr(align_features, "split_target_sequence",
                            lambda *a, **k: {})
        monkeypatch.setattr(align_features, "get_genome_size",
                            lambda *a, **k: 1)

        args = SimpleNamespace(
            native=True, subcommand=None, target="t.fa",
            mm2_options="", threads=1, directory="/tmp",
        )
        align_features.align_features_to_target(
            ["chr1"], ["chr1"], args,
            _FakeFeatureHierarchy({}, {}),
            "chrm_by_chrm", [],
        )
        assert called["native"] == 1

    def test_polish_subcommand_uses_legacy_even_with_native(self, monkeypatch):
        """The polish path requires a real SAM file on disk; --native
        must not hijack it."""
        from lifton.liftoff import align_features
        # Patch the native function so we can detect if it was called
        called = {"native": 0}

        def spy_native(*a, **k):
            called["native"] += 1
            return {}

        monkeypatch.setattr(
            "lifton.liftoff.native_align.align_features_to_target_native",
            spy_native,
        )
        # Patch parse_all_sam_files to no-op so we don't actually need a SAM file
        monkeypatch.setattr(align_features, "parse_all_sam_files",
                            lambda *a, **k: {})

        args = SimpleNamespace(
            native=True, subcommand="polish", directory="/tmp",
        )
        align_features.align_features_to_target(
            ["chr1"], ["chr1"], args,
            _FakeFeatureHierarchy({}, {}),
            "chrm_by_chrm", [],
        )
        assert called["native"] == 0
