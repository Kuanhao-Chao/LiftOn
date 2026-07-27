"""Unit tests for the Phase-3 structural evaluation metrics
(benchmarks/compare/evaluator.py): coordinate-independent intron-chain
exactness / exon Sn·Sp (boundary-level) and ORF-validity.

These are EVAL-only (they score emitted GFF3s in the A/B harness); they do not
touch the LiftOn pipeline, so the 24-cell byte-identity matrix is unaffected.
The metrics are coordinate-independent by construction: the reference and the
lifted transcript live on different genomes, so every comparison is on spliced
5'->3' boundary POSITIONS (cumulative segment lengths), never absolute coords.
"""
from types import SimpleNamespace

from benchmarks.compare import evaluator as ev


def _seg(start, end):
    return SimpleNamespace(start=start, end=end)


# --------------------------- internal boundaries --------------------------- #

class TestInternalBoundaries:
    def test_single_segment_has_no_internal_junction(self):
        assert ev._internal_boundaries([_seg(1, 100)], "+") == []

    def test_empty_is_empty(self):
        assert ev._internal_boundaries([], "+") == []

    def test_two_segments_forward(self):
        # exon1 len 100 -> one internal boundary at cumulative length 100.
        assert ev._internal_boundaries([_seg(1, 100), _seg(201, 300)], "+") == [100]

    def test_three_segments_forward(self):
        # lens 100, 150, 50 -> boundaries at 100 and 250.
        segs = [_seg(1, 100), _seg(201, 350), _seg(401, 450)]
        assert ev._internal_boundaries(segs, "+") == [100, 250]

    def test_unsorted_input_is_sorted(self):
        segs = [_seg(401, 450), _seg(1, 100), _seg(201, 350)]
        assert ev._internal_boundaries(segs, "+") == [100, 250]

    def test_reverse_strand_counts_from_3prime_end(self):
        # Same genomic segments (lens 100,150,50) but on '-' -> 5'->3' order is
        # reversed (50,150,100) -> boundaries at 50 and 200.
        segs = [_seg(1, 100), _seg(201, 350), _seg(401, 450)]
        assert ev._internal_boundaries(segs, "-") == [50, 200]

    def test_coordinate_independence(self):
        # Same structure at a different genomic locus -> identical boundaries.
        a = [_seg(1, 100), _seg(201, 300)]
        b = [_seg(9001, 9100), _seg(9500, 9599)]
        assert ev._internal_boundaries(a, "+") == ev._internal_boundaries(b, "+")


# --------------------------- boundary Sn/Sp/exact -------------------------- #

class TestBoundarySnSp:
    def test_identical_chains(self):
        sn, sp, exact = ev._boundary_snsp([100, 250], [100, 250])
        assert (sn, sp, exact) == (1.0, 1.0, 1)

    def test_single_exon_both(self):
        sn, sp, exact = ev._boundary_snsp([], [])
        assert (sn, sp, exact) == (1.0, 1.0, 1)

    def test_missing_one_junction_lowers_sn(self):
        # ref has 2 junctions, lifted recovered 1 -> sn 0.5, sp 1.0, not exact.
        sn, sp, exact = ev._boundary_snsp([100, 250], [100])
        assert sn == 0.5 and sp == 1.0 and exact == 0

    def test_extra_junction_lowers_sp(self):
        sn, sp, exact = ev._boundary_snsp([100], [100, 250])
        assert sn == 1.0 and sp == 0.5 and exact == 0

    def test_shifted_boundary_is_not_exact(self):
        # a 1 bp indel shifts the junction -> no exact match (correctly strict).
        sn, sp, exact = ev._boundary_snsp([100, 250], [100, 251])
        assert exact == 0 and sn == 0.5

    def test_empty_ref_is_vacuously_recovered(self):
        sn, sp, exact = ev._boundary_snsp([], [100])
        assert sn == 1.0 and sp == 0.0 and exact == 0


# --------------------------- ORF validity ---------------------------------- #

class TestOrfValidity:
    def test_clean_orf(self):
        assert ev._orf_validity("MKLLW*") == (1, 1, 1, 1)

    def test_missing_start(self):
        assert ev._orf_validity("KLLW*") == (0, 1, 1, 0)

    def test_missing_stop(self):
        assert ev._orf_validity("MKLLW") == (1, 0, 1, 0)

    def test_internal_stop(self):
        # internal '*' before the terminal -> not valid.
        assert ev._orf_validity("MK*LW*") == (1, 1, 0, 0)

    def test_empty(self):
        assert ev._orf_validity("") == (0, 0, 0, 0)

    def test_start_stop_but_internal_stop_only(self):
        # starts M, ends *, but has an internal stop -> valid 0.
        s, e, ni, v = ev._orf_validity("M*A*")
        assert s == 1 and e == 1 and ni == 0 and v == 0
