"""Tests for the anchor-windowed aligner (lifton/windowed_align.py) and the
length-gated routing in lifton/align.py.

The windowed path replaces O(L²) full parasail DP for giant genes (titin et al.)
that otherwise cost tens of GB + tens of seconds. These tests pin:
  - losslessness (the gapped alignment reconstructs the inputs exactly),
  - identity equivalence to full-DP on near-diagonal sequences,
  - parasail-compatible CIGAR (.cigar.decode bytes; =/X/I/D convention),
  - the length gate: BELOW threshold the exact full-DP path is used unchanged
    (byte-identical), ABOVE threshold the windowed shim is used.
"""

import random

import parasail
import pytest

from lifton import align, windowed_align, get_id_fraction

MAA = parasail.Matrix("blosum62")
MDNA = parasail.matrix_create("ACGTN*", 1, -3)
AA = "ACDEFGHIKLMNPQRSTVWY"
NT = "ACGT"


def _rand(alphabet, n, seed):
    r = random.Random(seed)
    return "".join(r.choice(alphabet) for _ in range(n))


def _mut(s, alphabet, rate, seed):
    r = random.Random(seed)
    return "".join(c if r.random() > rate else r.choice(alphabet) for c in s)


def _aa_id(res):
    m, l = get_id_fraction.get_AA_id_fraction(res.traceback.ref, res.traceback.query)
    return m / l


def _dna_id(res):
    m, l = get_id_fraction.get_DNA_id_fraction(res.traceback.ref, res.traceback.query)
    return m / l


def _ungap(s):
    return s.replace("-", "")


# ---------------------------------------------------------------------------
# Losslessness + structural invariants
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("seed", [1, 2, 3])
def test_protein_alignment_is_lossless(seed):
    ref = _rand(AA, 3000, seed)
    q = _mut(ref, AA, 0.03, seed + 100)
    q = q[:1000] + q[1003:] + "KK"            # a deletion + insertion
    res = windowed_align.windowed_traceback(q, ref, 11, 1, MAA, is_dna=False)
    assert len(res.traceback.query) == len(res.traceback.ref)
    assert _ungap(res.traceback.query) == q
    assert _ungap(res.traceback.ref) == ref


@pytest.mark.parametrize("seed", [4, 5])
def test_dna_alignment_is_lossless(seed):
    ref = _rand(NT, 4000, seed)
    q = _mut(ref, NT, 0.03, seed + 100)
    q = q[:1500] + "ACG" + q[1500:]
    res = windowed_align.windowed_traceback(q, ref, 5, 2, MDNA, is_dna=True)
    assert len(res.traceback.query) == len(res.traceback.ref)
    assert _ungap(res.traceback.query) == q
    assert _ungap(res.traceback.ref) == ref


def test_identical_sequence_is_perfect():
    ref = _rand(AA, 2500, 9)
    res = windowed_align.windowed_traceback(ref, ref, 11, 1, MAA, is_dna=False)
    assert res.traceback.query == ref
    assert res.traceback.ref == ref
    assert _aa_id(res) == 1.0
    assert res.cigar.decode.decode() == f"{len(ref)}="


# ---------------------------------------------------------------------------
# Identity equivalence to full-DP (the core correctness claim)
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("seed", [11, 12, 13])
def test_protein_identity_matches_full_dp(seed):
    ref = _rand(AA, 3000, seed)
    q = _mut(ref, AA, 0.02, seed + 50)
    q = q[:800] + q[803:]                      # small deletion
    q = q[:1900] + "WW" + q[1900:]             # small insertion
    full = parasail.nw_trace_scan_sat(q, ref, 11, 1, MAA)
    win = windowed_align.windowed_traceback(q, ref, 11, 1, MAA, is_dna=False)
    assert abs(_aa_id(win) - _aa_id(full)) < 1e-3


@pytest.mark.parametrize("seed", [21, 22, 23])
def test_dna_identity_matches_full_dp(seed):
    ref = _rand(NT, 5000, seed)
    q = _mut(ref, NT, 0.02, seed + 50)
    q = q[:1200] + q[1206:]                     # 6 nt deletion
    full = parasail.nw_trace_scan_sat(q, ref, 5, 2, MDNA)
    win = windowed_align.windowed_traceback(q, ref, 5, 2, MDNA, is_dna=True)
    assert abs(_dna_id(win) - _dna_id(full)) < 1e-3


def test_no_anchor_divergent_sequences_still_valid():
    """Two unrelated sequences share few/no unique k-mers → blind-split
    fallback. Must still return a valid, lossless, equal-length alignment."""
    a = _rand(AA, 2000, 31)
    b = _rand(AA, 2000, 99)
    res = windowed_align.windowed_traceback(a, b, 11, 1, MAA, is_dna=False)
    assert len(res.traceback.query) == len(res.traceback.ref)
    assert _ungap(res.traceback.query) == a
    assert _ungap(res.traceback.ref) == b


# ---------------------------------------------------------------------------
# CIGAR shim (parasail-compatible)
# ---------------------------------------------------------------------------

def test_cigar_matches_parasail_convention():
    cases = [("AAAAWWCCCC", "AAAACCCC"),    # 4=2I4=  (gap in ref)
             ("AAAACCCC", "AAAAWWCCCC"),    # 4=2D4=  (gap in query)
             ("AAAAKAAAA", "AAAADAAAA")]    # 4=1X4=  (mismatch)
    for q, r in cases:
        res = parasail.nw_trace_scan_sat(q, r, 11, 1, MAA)
        derived = windowed_align._cigar_from_aln(
            res.traceback.query, res.traceback.ref).decode()
        assert derived == res.cigar.decode.decode()


def test_cigar_is_consistent_with_traceback():
    ref = _rand(AA, 2000, 7)
    q = _mut(ref, AA, 0.05, 77)
    q = q[:500] + "PP" + q[500:]
    res = windowed_align.windowed_traceback(q, ref, 11, 1, MAA, is_dna=False)
    cig = res.cigar.decode.decode()
    # op lengths sum to the alignment length
    import re
    total = sum(int(n) for n in re.findall(r"(\d+)[=XID]", cig))
    assert total == len(res.traceback.query) == len(res.traceback.ref)
    # re-deriving from the returned strings reproduces the same cigar
    assert windowed_align._cigar_from_aln(
        res.traceback.query, res.traceback.ref).decode() == cig


# ---------------------------------------------------------------------------
# Length gate in align.py (byte-identity below threshold)
# ---------------------------------------------------------------------------

def test_protein_base_below_threshold_is_exact_full_dp():
    ref = _rand(AA, 1500, 5)             # < the band-everything gate (2500 aa)
    q = _mut(ref, AA, 0.05, 6)
    got = align.parasail_align_protein_base(q, ref)
    exp = parasail.nw_trace_scan_sat(q, ref, 11, 1, MAA)
    # below threshold the real parasail kernel is used → byte-identical traceback
    assert got.traceback.query == exp.traceback.query
    assert got.traceback.ref == exp.traceback.ref
    assert not isinstance(got, windowed_align._ShimResult)


def test_protein_base_above_threshold_is_windowed():
    ref = _rand(AA, 9000, 5)             # > _ALIGN_WINDOW_AA (any mode)
    q = _mut(ref, AA, 0.02, 6)
    got = align.parasail_align_protein_base(q, ref)
    assert isinstance(got, windowed_align._ShimResult)
    assert got.cigar.decode.decode()    # cigar path works for downstream boundary code


def test_dna_base_above_threshold_is_windowed():
    ref = _rand(NT, 30000, 5)            # > _ALIGN_WINDOW_NT (any mode)
    q = _mut(ref, NT, 0.02, 6)
    got = align.parasail_align_DNA_base(q, ref)
    assert isinstance(got, windowed_align._ShimResult)
    assert len(got.traceback.query) == len(got.traceback.ref)


# ---------------------------------------------------------------------------
# Iteration 3: "band everything" is the DEFAULT — low gate (2500 aa / 8000 nt),
# fine band (cap 1500), exact-DP fallback raised to the giant boundary. The
# exact giant-only path is restored by configure_alignment(band=False)
# (--full-dp-align). These pin the env cap, the max_fulldp exact fallback, and
# the symmetric band/giant toggle.
# ---------------------------------------------------------------------------

# Disjoint AA alphabets → the two sequences share NO k-mer at any k → genuinely
# anchor-less, so the fallback branch is actually exercised (random sequences
# would still produce spurious low-k anchors and never reach it).
AA_LEFT = "ACDEFGHIKL"
AA_RIGHT = "MNPQRSTVWY"


@pytest.fixture(autouse=True)
def reset_align_state():
    """Reset the module-global alignment config to the band-everything DEFAULT
    before each test and restore the prior state after, so neither these tests
    nor pipeline tests in other files leak mutated thresholds/cap across the
    suite (the config lives in module globals)."""
    saved = (align._ALIGN_WINDOW_AA, align._ALIGN_WINDOW_NT,
             align._ALIGN_MAX_FULLDP_AA, align._ALIGN_MAX_FULLDP_NT,
             align._FAST_ALIGN_ACTIVE, windowed_align.WINDOW_CAP)
    align.configure_alignment(band=True)       # known default state
    try:
        yield
    finally:
        (align._ALIGN_WINDOW_AA, align._ALIGN_WINDOW_NT,
         align._ALIGN_MAX_FULLDP_AA, align._ALIGN_MAX_FULLDP_NT,
         align._FAST_ALIGN_ACTIVE, windowed_align.WINDOW_CAP) = saved


def test_default_mode_is_band_everything():
    """No configure call: the imported module state is band-everything, so a
    3000 aa homolog pair (> 2500 gate) routes through the windowed shim."""
    assert align._FAST_ALIGN_ACTIVE is True
    assert align._ALIGN_WINDOW_AA == align._BAND_AA == 2500
    assert windowed_align.WINDOW_CAP == align._BAND_CAP == 1500
    ref = _rand(AA, 3000, 71)
    q = _mut(ref, AA, 0.03, 72)
    got = align.parasail_align_protein_base(q, ref)
    assert isinstance(got, windowed_align._ShimResult)
    exp = parasail.nw_trace_scan_sat(q, ref, 11, 1, MAA)
    assert abs(_aa_id(got) - _aa_id(exp)) < 1e-3       # homolog → identity preserved


def test_band_everything_identity_matches_full_dp():
    """With a finely-banded cap, a HOMOLOGOUS mid-size pair (anchors abundant
    and correct) reproduces full-DP identity — the core band-everything claim."""
    windowed_align.set_window_cap(1500)
    ref = _rand(AA, 3000, 41)
    q = _mut(ref, AA, 0.05, 42)
    q = q[:900] + q[904:] + "KK"               # indels off the diagonal
    full = parasail.nw_trace_scan_sat(q, ref, 11, 1, MAA)
    win = windowed_align.windowed_traceback(q, ref, 11, 1, MAA, is_dna=False,
                                            max_fulldp=8000)
    assert _ungap(win.traceback.query) == q
    assert _ungap(win.traceback.ref) == ref
    assert abs(_aa_id(win) - _aa_id(full)) < 1e-3


def test_exact_fallback_anchorless_nongiant_matches_full_dp():
    """An anchor-less (disjoint-alphabet) NON-giant region with max_fulldp raised
    above its length must use exact full DP → byte-identical traceback to
    parasail (no approximate blind-split)."""
    windowed_align.set_window_cap(1500)
    q = _rand(AA_LEFT, 3000, 7)
    r = _rand(AA_RIGHT, 3000, 8)               # shares no k-mer with q
    full = parasail.nw_trace_scan_sat(q, r, 11, 1, MAA)
    win = windowed_align.windowed_traceback(q, r, 11, 1, MAA, is_dna=False,
                                            max_fulldp=8000)
    assert win.traceback.query == full.traceback.query
    assert win.traceback.ref == full.traceback.ref


def test_anchorless_giant_uses_bounded_blind_split():
    """An anchor-less region ABOVE max_fulldp falls back to the bounded
    blind-split — still lossless and equal-length (memory-safe approximation)."""
    windowed_align.set_window_cap(1500)
    q = _rand(AA_LEFT, 4000, 11)
    r = _rand(AA_RIGHT, 4000, 12)
    win = windowed_align.windowed_traceback(q, r, 11, 1, MAA, is_dna=False,
                                            max_fulldp=2000)   # 4000 > 2000
    assert len(win.traceback.query) == len(win.traceback.ref)
    assert _ungap(win.traceback.query) == q
    assert _ungap(win.traceback.ref) == r


def test_max_fulldp_defaults_to_window_cap():
    """max_fulldp=None reproduces the giant-only behaviour: the exact fallback
    is off, so an anchor-less region > WINDOW_CAP blind-splits rather than
    full-DP'ing (identity differs from exact for non-homologous input)."""
    windowed_align.set_window_cap(1500)
    q = _rand(AA_LEFT, 3000, 21)
    r = _rand(AA_RIGHT, 3000, 22)
    none_res = windowed_align.windowed_traceback(q, r, 11, 1, MAA, is_dna=False)
    full = parasail.nw_trace_scan_sat(q, r, 11, 1, MAA)
    # default (None→cap) blind-splits → not byte-identical to exact full DP,
    # but still a valid lossless alignment.
    assert none_res.traceback.query != full.traceback.query
    assert _ungap(none_res.traceback.query) == q


def test_set_window_cap_roundtrip():
    # autouse reset leaves the band-everything default cap in place
    assert windowed_align.WINDOW_CAP == align._BAND_CAP == 1500
    windowed_align.set_window_cap(1234)
    assert windowed_align.WINDOW_CAP == 1234


def test_configure_alignment_is_symmetric():
    align.configure_alignment(band=True)
    assert align._ALIGN_WINDOW_AA == align._BAND_AA == 2500
    assert align._ALIGN_WINDOW_NT == align._BAND_NT == 8000
    assert windowed_align.WINDOW_CAP == align._BAND_CAP == 1500
    assert align._ALIGN_MAX_FULLDP_AA == align._GIANT_AA == 8000
    assert align._ALIGN_MAX_FULLDP_NT == align._GIANT_NT == 25000
    assert align._FAST_ALIGN_ACTIVE is True

    align.configure_alignment(band=False)      # --full-dp-align / giant-only
    assert align._ALIGN_WINDOW_AA == align._GIANT_AA == 8000
    assert align._ALIGN_WINDOW_NT == align._GIANT_NT == 25000
    assert windowed_align.WINDOW_CAP == align._GIANT_CAP == 6000
    assert align._ALIGN_MAX_FULLDP_AA is None
    assert align._ALIGN_MAX_FULLDP_NT is None
    assert align._FAST_ALIGN_ACTIVE is False

    # legacy alias still works
    align.configure_fast_align(True)
    assert align._FAST_ALIGN_ACTIVE is True


def test_full_dp_align_keeps_midsize_on_full_dp():
    """--full-dp-align (band=False): a 3000 aa pair (< 8000 giant gate) stays
    exact full DP. Default band mode (> 2500 gate) routes through the shim."""
    ref = _rand(AA, 3000, 5)
    q = _mut(ref, AA, 0.03, 6)
    align.configure_alignment(band=False)
    off = align.parasail_align_protein_base(q, ref)
    assert not isinstance(off, windowed_align._ShimResult)
    exp = parasail.nw_trace_scan_sat(q, ref, 11, 1, MAA)
    assert off.traceback.query == exp.traceback.query

    align.configure_alignment(band=True)
    on = align.parasail_align_protein_base(q, ref)
    assert isinstance(on, windowed_align._ShimResult)
    assert abs(_aa_id(on) - _aa_id(exp)) < 1e-3       # homolog → identity preserved
