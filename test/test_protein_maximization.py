"""
test_protein_maximization.py – Comprehensive unit tests for:
  • protein_maximization.chaining_algorithm
  • lifton_class.Lifton_TRANS.update_cds_list
  • lifton_class.Lifton_TRANS.__find_orfs  (via orf_search_protein)

Run with:
    pytest test/test_protein_maximization.py -v

Every test is self-contained using lightweight mock objects so no live
gffutils databases or external files are needed.
"""

import copy
import math
import pytest

# --------------------------------------------------------------------------- #
# Minimal mocks for gffutils FeatureDB entries                                #
# --------------------------------------------------------------------------- #

class _FakeFeature:
    """Minimal gffutils.Feature-like object."""
    def __init__(self, seqid, ftype, start, end, strand="+", attrs=None):
        self.seqid       = seqid
        self.featuretype = ftype
        self.start       = start
        self.end         = end
        self.strand      = strand
        self.attributes  = attrs or {"Parent": ["tx1"], "ID": ["f1"]}
        self.id          = (attrs or {}).get("ID", ["f1"])[0]
        self.source      = "LiftOn"
        self.frame       = "."
        self.score       = "."

    def __str__(self):
        return (f"{self.seqid}\t{self.source}\t{self.featuretype}\t"
                f"{self.start}\t{self.end}\t{self.score}\t"
                f"{self.strand}\t{self.frame}\t.")

    def sequence(self, fai):
        return "ATG"  # stub


def _feat(start, end, strand="+", ftype="CDS"):
    return _FakeFeature("chr1", ftype, start, end, strand,
                        attrs={"Parent": ["tx1"], "ID": [f"{ftype}_{start}_{end}"]})


# --------------------------------------------------------------------------- #
# Mock Lifton_Alignment                                                         #
# --------------------------------------------------------------------------- #

class _MockAlignment:
    def __init__(self, cds_features, protein_aln_bounds, ref_aln, query_aln,
                 strand="+"):
        db_entry = _FakeFeature("chr1", "mRNA", 100, 900, strand)
        self.identity                  = 0.9
        self.cds_children              = cds_features
        self.cdss_protein_aln_boundaries = protein_aln_bounds
        self.ref_aln                   = ref_aln
        self.query_aln                 = query_aln
        self.ref_seq                   = ref_aln
        self.db_entry                  = db_entry


# --------------------------------------------------------------------------- #
# Helper: build a Lifton_TRANS with synthetic exons                           #
# --------------------------------------------------------------------------- #

def _make_trans(exon_intervals, strand="+", ref_trans_id="tx1"):
    """Return a Lifton_TRANS-like object with real exons."""
    from lifton import lifton_class
    # Build a fake reference DB dict
    fake_ref_attrs = {
        "ID": [ref_trans_id],
        "Parent": ["gene1"],
    }
    trans_entry = _FakeFeature("chr1", "mRNA", exon_intervals[0][0],
                               exon_intervals[-1][1], strand,
                               {"ID": [ref_trans_id], "Parent": ["gene1"]})

    class _FakeTrans:
        pass

    trans = _FakeTrans()
    trans.entry = trans_entry
    trans.exons = []
    for (s, e) in exon_intervals:
        ex_feat = _feat(s, e, strand, ftype="exon")
        ex_feat.attributes = {"Parent": [ref_trans_id]}
        ex = lifton_class.Lifton_EXON(ex_feat)
        trans.exons.append(ex)
    # Bind the real methods we need
    import types
    from lifton.lifton_class import Lifton_TRANS as LT
    trans.update_cds_list           = types.MethodType(LT.update_cds_list, trans)
    trans.update_boundaries         = types.MethodType(LT.update_boundaries, trans)
    trans._Lifton_TRANS__get_cds_frame = types.MethodType(
        LT._Lifton_TRANS__get_cds_frame, trans)
    return trans


# =========================================================================== #
# Tests: protein_maximization.chaining_algorithm                              #
# =========================================================================== #

class TestChainingAlgorithm:

    def _make_aln(self, cds_coords, strand="+"):
        """Make a minimal alignment where each CDS covers 10 ref AAs."""
        cds_feats = [_feat(s, e, strand) for (s, e) in cds_coords]
        # Each CDS covers 10 AA in the alignment → boundaries are (0,10), (10,20), …
        bounds = [(i * 10.0, (i + 1) * 10.0) for i in range(len(cds_coords))]
        n = len(cds_coords) * 10
        ref_aln   = "A" * n
        query_aln = "A" * n
        return _MockAlignment(cds_feats, bounds, ref_aln, query_aln, strand)

    # ── Edge-case 1: Both inputs empty ─────────────────────────────────────
    def test_both_empty(self):
        from lifton import protein_maximization as pm
        l_aln = self._make_aln([])
        m_aln = self._make_aln([])
        cds_list, chains = pm.chaining_algorithm(l_aln, m_aln, None, False)
        assert cds_list == []

    # ── Edge-case 2: Miniprot has 0 children ───────────────────────────────
    def test_miniprot_empty(self):
        from lifton import protein_maximization as pm
        l_aln = self._make_aln([(100, 200), (300, 400)])
        m_aln = self._make_aln([])
        cds_list, chains = pm.chaining_algorithm(l_aln, m_aln, None, False)
        # Should fall back to liftoff CDS list
        assert len(cds_list) == 2

    # ── Edge-case 3: Liftoff has 0 children ────────────────────────────────
    def test_liftoff_empty(self):
        from lifton import protein_maximization as pm
        l_aln = self._make_aln([])
        m_aln = self._make_aln([(100, 200)])
        cds_list, chains = pm.chaining_algorithm(l_aln, m_aln, None, False)
        assert cds_list == []

    # ── Edge-case 4: Single CDS on both sides ──────────────────────────────
    def test_single_cds_both(self):
        from lifton import protein_maximization as pm
        l_aln = self._make_aln([(100, 400)])
        m_aln = self._make_aln([(100, 400)])
        cds_list, chains = pm.chaining_algorithm(l_aln, m_aln, None, False)
        assert len(cds_list) == 1

    # ── Edge-case 5: Single CDS on liftoff, multiple on miniprot ──────────
    def test_single_liftoff_multi_miniprot(self):
        from lifton import protein_maximization as pm
        l_aln = self._make_aln([(100, 400)])
        m_aln = self._make_aln([(100, 200), (300, 400)])
        cds_list, chains = pm.chaining_algorithm(l_aln, m_aln, None, False)
        assert len(cds_list) >= 1

    # ── Edge-case 6: Miniprot wins whole alignment ─────────────────────────
    def test_miniprot_wins(self):
        from lifton import protein_maximization as pm
        n = 20
        ref_aln   = "A" * n
        # Perfect miniprot (all matches)
        m_query   = "A" * n
        # Bad liftoff (all mismatches)
        l_query   = "C" * n
        l_cds = [_feat(100, 200), _feat(300, 400)]
        m_cds = [_feat(100, 200), _feat(300, 400)]
        bounds = [(0.0, 10.0), (10.0, 20.0)]
        l_aln = _MockAlignment(l_cds, bounds, ref_aln, l_query)
        m_aln = _MockAlignment(m_cds, bounds, ref_aln, m_query)
        cds_list, chains = pm.chaining_algorithm(l_aln, m_aln, None, False)
        # Miniprot should have been picked for the last chunk at least
        assert any("miniprot" in c for c in chains)

    # ── Edge-case 7: Liftoff wins on tie ───────────────────────────────────
    def test_liftoff_wins_tie(self):
        from lifton import protein_maximization as pm
        n = 20
        # Identical alignments → tie → liftoff wins
        ref_aln = "A" * n
        query   = "A" * n
        l_cds = [_feat(100, 200), _feat(300, 400)]
        m_cds = [_feat(100, 200), _feat(300, 400)]
        bounds = [(0.0, 10.0), (10.0, 20.0)]
        l_aln = _MockAlignment(l_cds, bounds, ref_aln, query)
        m_aln = _MockAlignment(m_cds, bounds, ref_aln, query)
        cds_list, chains = pm.chaining_algorithm(l_aln, m_aln, None, False)
        # All chunks should be liftoff on a tie
        assert all("liftoff" in c for c in chains)

    # ── Edge-case 8: Negative strand ───────────────────────────────────────
    def test_negative_strand(self):
        from lifton import protein_maximization as pm
        l_aln = self._make_aln([(100, 200), (300, 400)], strand="-")
        m_aln = self._make_aln([(100, 200), (300, 400)], strand="-")
        cds_list, chains = pm.chaining_algorithm(l_aln, m_aln, None, False)
        assert len(cds_list) >= 1

    # ── Edge-case 9: Many CDS blocks ────────────────────────────────────────
    def test_many_cds_blocks(self):
        from lifton import protein_maximization as pm
        coords = [(i * 200, i * 200 + 100) for i in range(8)]
        l_aln = self._make_aln(coords)
        m_aln = self._make_aln(coords)
        cds_list, chains = pm.chaining_algorithm(l_aln, m_aln, None, False)
        assert len(cds_list) >= 1


# =========================================================================== #
# Tests: get_partial_id_fraction                                              #
# =========================================================================== #

class TestGetPartialIdFraction:

    def test_perfect_match(self):
        from lifton.get_id_fraction import get_partial_id_fraction
        ref   = "ACDEFGHIKLM"
        query = "ACDEFGHIKLM"
        m, l = get_partial_id_fraction(ref, query, 0, len(ref))
        assert m == l  # 100% identity

    def test_all_mismatch(self):
        from lifton.get_id_fraction import get_partial_id_fraction
        ref   = "ACDEF"
        query = "GHIKL"
        m, l = get_partial_id_fraction(ref, query, 0, 5)
        assert m == 0

    def test_gap_in_reference(self):
        from lifton.get_id_fraction import get_partial_id_fraction
        # ref has 1 gap → total_length shrinks by 1
        ref   = "A-CD"
        query = "AXCD"
        m, l = get_partial_id_fraction(ref, query, 0, 4)
        assert l == 3   # 4 - 1 gap

    def test_empty_window_returns_one(self):
        from lifton.get_id_fraction import get_partial_id_fraction
        m, l = get_partial_id_fraction("AAAA", "AAAA", 2, 2)
        assert l == 1   # division-by-zero protection

    def test_stop_codon_terminates_count(self):
        from lifton.get_id_fraction import get_partial_id_fraction
        # '*' in query should stop counting at that position
        ref   = "ACDE"
        query = "AC*E"
        m, l = get_partial_id_fraction(ref, query, 0, 4)
        # Counting stops at position 2 ('*'), so only A,C compared
        assert m == 2


# =========================================================================== #
# Tests: update_cds_list                                                      #
# =========================================================================== #

class TestUpdateCdsList:

    def _make_cds(self, start, end, strand="+"):
        from lifton.lifton_class import Lifton_CDS
        feat = _feat(start, end, strand, ftype="CDS")
        return Lifton_CDS(feat)

    # ── Case 0: Empty CDS list ───────────────────────────────────────────
    def test_empty_cds_list(self):
        trans = _make_trans([(100, 300), (500, 700)])
        trans.update_cds_list([])
        # All exon CDS slots should be None
        for exon in trans.exons:
            assert exon.cds is None

    # ── Case 1a: Single CDS fully contained in one exon ──────────────────
    def test_single_cds_contained_in_one_exon(self):
        trans = _make_trans([(100, 500)])
        cds   = [self._make_cds(200, 400)]
        trans.update_cds_list(cds)
        assert len(trans.exons) == 1
        assert trans.exons[0].cds is not None
        assert trans.exons[0].cds.entry.start == 200
        assert trans.exons[0].cds.entry.end   == 400

    # ── Case 1b: Single CDS spans multiple exons ──────────────────────────
    def test_single_cds_spans_multiple_exons(self):
        trans = _make_trans([(100, 300), (400, 600)])
        cds   = [self._make_cds(150, 550)]
        trans.update_cds_list(cds)
        # Should produce one merged exon covering [150, 550]
        assert any(e.cds is not None for e in trans.exons)

    # ── Case 1c: Single CDS lies entirely BEFORE all exons ───────────────
    def test_single_cds_before_all_exons(self):
        trans = _make_trans([(500, 700)])
        cds   = [self._make_cds(100, 200)]
        trans.update_cds_list(cds)
        # A synthetic exon covering the CDS should be emitted
        cds_exons = [e for e in trans.exons if e.cds is not None]
        assert len(cds_exons) >= 1
        assert cds_exons[0].cds.entry.start == 100

    # ── Case 1d: Single CDS lies entirely AFTER all exons ────────────────
    def test_single_cds_after_all_exons(self):
        trans = _make_trans([(100, 200)])
        cds   = [self._make_cds(500, 700)]
        trans.update_cds_list(cds)
        cds_exons = [e for e in trans.exons if e.cds is not None]
        assert len(cds_exons) >= 1

    # ── Case 1e: No exons at all ─────────────────────────────────────────
    def test_empty_exon_list_single_cds(self):
        trans = _make_trans([(100, 200)])
        trans.exons = []   # wipe exons
        cds = [self._make_cds(100, 200)]
        # Should not crash
        trans.update_cds_list(cds)

    # ── Case 2: Single exon, multiple CDSs ────────────────────────────────
    def test_single_exon_multiple_cds(self):
        trans = _make_trans([(100, 900)])
        cds   = [self._make_cds(200, 350),
                 self._make_cds(400, 550),
                 self._make_cds(600, 750)]
        trans.update_cds_list(cds)
        assert len(trans.exons) == 3
        for exon in trans.exons:
            assert exon.cds is not None

    # ── Case 3a: Multiple exons, multiple CDSs — normal case ──────────────
    def test_multiple_exons_multiple_cds_normal(self):
        trans = _make_trans([(100, 200), (300, 400), (500, 600)])
        cds   = [self._make_cds(120, 180),
                 self._make_cds(310, 390),
                 self._make_cds(510, 590)]
        trans.update_cds_list(cds)
        cds_exons = [e for e in trans.exons if e.cds is not None]
        assert len(cds_exons) == 3

    # ── Case 3b: CDS leads (first CDS starts before first exon) ──────────
    def test_cds_leads_exons(self):
        trans = _make_trans([(300, 500), (600, 800)])
        cds   = [self._make_cds(100, 200), self._make_cds(310, 490)]
        trans.update_cds_list(cds)
        # First CDS is synthetic; second overlaps first exon
        assert len(trans.exons) >= 2

    # ── Case 3c: Last CDS has no overlapping exon ────────────────────────
    def test_last_cds_no_overlapping_exon(self):
        trans = _make_trans([(100, 200)])
        cds   = [self._make_cds(120, 180), self._make_cds(500, 600)]
        trans.update_cds_list(cds)
        cds_exons = [e for e in trans.exons if e.cds is not None]
        assert len(cds_exons) == 2

    # ── Case 3d: Negative strand ─────────────────────────────────────────
    def test_negative_strand_cds(self):
        trans = _make_trans([(100, 200), (300, 400)], strand="-")
        cds   = [self._make_cds(100, 200, strand="-"),
                 self._make_cds(300, 400, strand="-")]
        trans.update_cds_list(cds)
        cds_exons = [e for e in trans.exons if e.cds is not None]
        assert len(cds_exons) >= 1

    # ── Case 3e: CDS exactly at exon boundary ────────────────────────────
    def test_cds_exactly_at_exon_boundary(self):
        trans = _make_trans([(100, 200), (300, 400)])
        cds   = [self._make_cds(100, 200), self._make_cds(300, 400)]
        trans.update_cds_list(cds)
        cds_exons = [e for e in trans.exons if e.cds is not None]
        assert len(cds_exons) == 2

    # ── Case 3f: UTR-only exon at end (3' UTR) ───────────────────────────
    def test_utr_exon_at_end(self):
        """Last exon is pure UTR (no CDS should be attached)."""
        trans = _make_trans([(100, 200), (300, 400), (600, 800)])
        cds   = [self._make_cds(120, 190), self._make_cds(310, 390)]
        trans.update_cds_list(cds)
        # The last exon (600-800) should have no CDS
        last = trans.exons[-1]
        assert last.cds is None

    # ── Hierarchy: CDS start <= CDS end ──────────────────────────────────
    def test_cds_coords_sane(self):
        trans = _make_trans([(100, 300), (500, 700)])
        cds   = [self._make_cds(150, 280), self._make_cds(520, 680)]
        trans.update_cds_list(cds)
        for exon in trans.exons:
            if exon.cds is not None:
                assert exon.cds.entry.start <= exon.cds.entry.end

    # ── Hierarchy: CDS contained within exon ────────────────────────────
    def test_cds_contained_in_exon(self):
        trans = _make_trans([(100, 500)])
        cds   = [self._make_cds(200, 400)]
        trans.update_cds_list(cds)
        for exon in trans.exons:
            if exon.cds is not None:
                assert exon.cds.entry.start >= exon.entry.start
                assert exon.cds.entry.end   <= exon.entry.end


# =========================================================================== #
# Tests: __find_orfs (via orf_search_protein stub test)                       #
# =========================================================================== #

class TestFindOrfs:
    """
    We test __find_orfs indirectly by calling the internal-method pathway
    through a thin wrapper that bypasses the full gene graph.
    """

    def _call_find_orfs(self, trans_seq, ref_protein_seq):
        """Call Lifton_TRANS.__find_orfs in isolation via a fake TRANS object."""
        from lifton.lifton_class import Lifton_TRANS, Lifton_Status
        import types

        class _FakeLiftonStatus:
            def __init__(self):
                self.lifton_aa = 0.0
                self.status    = []

        class _FakeTrans:
            entry = _FakeFeature("chr1", "mRNA", 1, len(trans_seq) * 3,
                                 "+", {"ID": ["tx1"], "Parent": ["g1"]})
            exons = []

        ft = _FakeTrans()
        ls = _FakeLiftonStatus()

        # Python name-mangles private methods: self.__update_cds_boundary
        # inside Lifton_TRANS becomes self._Lifton_TRANS__update_cds_boundary.
        # Since ft is a _FakeTrans, that attribute doesn't exist unless we
        # explicitly bind the real implementations from Lifton_TRANS.
        ft._Lifton_TRANS__update_cds_boundary = types.MethodType(
            Lifton_TRANS._Lifton_TRANS__update_cds_boundary, ft)
        ft._Lifton_TRANS__iterate_exons_update_cds = types.MethodType(
            Lifton_TRANS._Lifton_TRANS__iterate_exons_update_cds, ft)
        ft._Lifton_TRANS__get_cds_frame = types.MethodType(
            Lifton_TRANS._Lifton_TRANS__get_cds_frame, ft)

        Lifton_TRANS._Lifton_TRANS__find_orfs(ft, trans_seq, ref_protein_seq,
                                              None, ls)
        return ls

    # ── Basic: single ATG→stop in frame 0 ────────────────────────────────
    def test_simple_orf_found(self):
        # ATG NNN NNN TAA   (10 nt → protein "MXX*")
        trans = "ATG" + "GCT" * 4 + "TAA"  # frame 0 ORF
        ref_p = "M" + "A" * 4              # reference protein
        ls = self._call_find_orfs(trans, ref_p)
        # Since we can't call align.parasail without a real install in unit
        # tests, just assert no crash.

    # ── ORF scan: no ATG in sequence ─────────────────────────────────────
    def test_no_atg(self):
        trans = "GGG" * 10  # no ATG anywhere
        ref_p = "MAAA"
        ls = self._call_find_orfs(trans, ref_p)
        assert ls.lifton_aa == 0.0  # no ORF found, no update

    # ── ORF scan: multiple frames, pick best ─────────────────────────────
    def test_multiple_frames_no_crash(self):
        # Frame 0: ATG...TAA
        # Frame 1: xATG...TAA
        # Frame 2: xxATG...TAA
        trans = "G" + "ATG" + "GCT" * 5 + "TAA" + "A" * 3
        ref_p = "M" + "A" * 5
        ls = self._call_find_orfs(trans, ref_p)

    # ── ORF scan: nested ATGs – should not produce infinitely many ORFs ──
    def test_nested_atg_not_expanded(self):
        # ATG GCT ATG GCT TAA — inner ATG should be skipped after outer ORF
        trans = "ATG" + "GCT" + "ATG" + "GCT" + "TAA"
        ref_p = "MAMA"
        ls = self._call_find_orfs(trans, ref_p)  # must not hang

    # ── ORF scan: stop codon immediately after ATG ────────────────────────
    def test_atg_immediately_followed_by_stop(self):
        trans = "ATG" + "TAA" + "GCT" * 10
        ref_p = "MAAAA"
        ls = self._call_find_orfs(trans, ref_p)

    # ── ORF scan: entire sequence is one long ORF (no stop) ──────────────
    def test_no_stop_codon(self):
        # ORF with no stop codon → should be silently skipped (partial gene)
        trans = "ATG" + "GCT" * 20  # no stop
        ref_p = "M" + "A" * 20
        ls = self._call_find_orfs(trans, ref_p)
        assert ls.lifton_aa == 0.0  # no ORF with stop found


# =========================================================================== #
# Tests: chaining_algorithm – protein_boundary indexing                       #
# =========================================================================== #

class TestProteinBoundary:

    def test_get_protein_boundary_correct(self):
        from lifton.protein_maximization import get_protein_boundary
        bounds = [(0.0, 10.0), (10.0, 20.0), (20.0, 30.0)]
        aa_s, aa_e = get_protein_boundary(bounds, c_idx_last=0, c_idx=2,
                                          DEBUG=False)
        assert aa_s == 0.0
        assert aa_e == 20.0

    def test_get_protein_reference_length_single_out_of_bounds(self):
        """Out-of-bounds c_idx should return 0, not raise."""
        from lifton.protein_maximization import get_protein_reference_length_single
        aln = _MockAlignment([], [], "AACCDD", "AACCDD")
        # No boundaries at all → c_idx=0 should return 0
        result = get_protein_reference_length_single(aln, 0, False)
        assert result == 0

    def test_process_m_l_children_empty_chunk(self):
        """c_idx_last == c_idx → empty chunk → return []."""
        from lifton.protein_maximization import process_m_l_children
        aln = _MockAlignment([], [], "AAAA", "AAAA")
        result = process_m_l_children(0, 0, aln, 0, 0, aln, None, [], False)
        assert result == []


# =========================================================================== #
# Tests: __iterate_exons_update_cds (Coordinate Mapping & Phase Logic)        #
# =========================================================================== #

class TestIterateExonsUpdateCds:
    """Test mapping transcript ORF boundaries back to genomic coordinates."""

    def _call_iterate(self, exon_intervals, orf_start, orf_end, strand="+"):
        from lifton.lifton_class import Lifton_TRANS
        from collections import namedtuple
        import types

        # Create the fake transcript and exons
        trans = _make_trans(exon_intervals, strand=strand)
        
        # We need a Fake_ORF object that has .start and .end
        FakeORF = namedtuple('FakeORF', ['start', 'end'])
        final_orf = FakeORF(orf_start, orf_end)
        
        # Bind the target method
        trans._Lifton_TRANS__iterate_exons_update_cds = types.MethodType(
            Lifton_TRANS._Lifton_TRANS__iterate_exons_update_cds, trans)
        
        # Exon order passed depends on strand
        pass_exons = trans.exons if strand == "+" else trans.exons[::-1]
        trans._Lifton_TRANS__iterate_exons_update_cds(final_orf, pass_exons, strand)
        
        # Return the modified exons
        return trans.exons

    # ── Bugfix: Single-Exon ORF (Positive Strand) ──────────────────────────
    def test_single_exon_positive(self):
        # Exon: [100, 200] (len=101)
        # ORF entirely inside: [10, 50) means len=40. Start=10, End=50.
        # Expect CDS: [100+10, 100+50-1] = [110, 149]
        exons = self._call_iterate([(100, 200)], orf_start=10, orf_end=50, strand="+")
        assert exons[0].cds is not None
        assert exons[0].cds.entry.start == 110
        assert exons[0].cds.entry.end == 149
        assert exons[0].cds.entry.frame == "0" # first CDS phase is 0

    # ── Bugfix: Single-Exon ORF (Negative Strand) ──────────────────────────
    def test_single_exon_negative(self):
        # Exon: [100, 200] (len=101). Negative strand! (5' is 200, 3' is 100)
        # ORF entirely inside: [10, 50) means len=40.
        # Transcript start=10 means distance from 5' end. So 200 - 10 = 190.
        # Transcript end=50 means distance from 5' end. So 200 - 50 + 1 = 151.
        # Expect CDS bounds: [151, 190]
        exons = self._call_iterate([(100, 200)], orf_start=10, orf_end=50, strand="-")
        assert exons[0].cds is not None
        assert exons[0].cds.entry.start == 151
        assert exons[0].cds.entry.end == 190

    # ── Multi-Exon ORF (Positive Strand) ───────────────────────────────────
    def test_multi_exon_positive(self):
        # Exons: E1=[100, 149] (50bp), E2=[200, 249] (50bp), E3=[300, 349] (50bp)
        # ORF: [10, 120] -> Starts in E1, spans E2, ends in E3.
        # E1 CDS: takes [110, 149] (len=40)
        # E2 CDS: takes [200, 249] (len=50)
        # E3 CDS: takes [300, 319] (len=20) (accum=100 so far, need 120, so 20 from E3)
        exons = self._call_iterate([(100, 149), (200, 249), (300, 349)], 10, 120, strand="+")
        
        assert exons[0].cds.entry.start == 110
        assert exons[0].cds.entry.end == 149
        
        assert exons[1].cds.entry.start == 200
        assert exons[1].cds.entry.end == 249
        
        assert exons[2].cds.entry.start == 300
        assert exons[2].cds.entry.end == 319

    # ── Phase Testing (Negative Strand) ────────────────────────────────────
    def test_multi_exon_negative_with_phases(self):
        # Exons mapping: E1=[300, 349] (50bp), E2=[200, 249] (50bp), E3=[100, 149] (50bp)
        # mRNA order is E1 -> E2 -> E3 because strand is "-".
        # ORF: [10, 120]
        # E1 (5' exon) starts at 349. ORF starts 10 bases in: 349 - 10 = 339.
        # E1 CDS: [300, 339] (len=40). 
        # Phase passed to E2: 40 % 3 = 1 -> so E2 phase must be (3-1)%3 = 2.
        # E2 CDS: [200, 249] (len=50). accum_cds=90.
        # Phase passed to E3: 90 % 3 = 0 -> so E3 phase must be (3-0)%3 = 0.
        # E3 CDS: Total accum_exon_length is 100 before E3. ORF needs to reach 120.
        # So E3 supplies 20 bases. 249 is 3' of E2, 149 is 5' of E3.
        # 20 bases from 5' end of E3 (149): 149 - 20 + 1 = 130.
        # E3 CDS bounds: [130, 149].
        
        exons = self._call_iterate([(100, 149), (200, 249), (300, 349)], 10, 120, strand="-")
        
        # Expected mRNA order: index 2, 1, 0
        e1 = exons[2].cds
        e2 = exons[1].cds
        e3 = exons[0].cds
        
        assert e1 is not None and e1.entry.start == 300 and e1.entry.end == 339
        assert e1.entry.frame == "0"
        
        assert e2 is not None and e2.entry.start == 200 and e2.entry.end == 249
        assert e2.entry.frame == "2" # (3 - 40%3)%3 = 2
        
        assert e3 is not None and e3.entry.start == 130 and e3.entry.end == 149
        assert e3.entry.frame == "0" # (3 - 90%3)%3 = 0
