"""Phase 4.5 Step 4 — property-based tests via Hypothesis.

Asserts core algebraic invariants over the coordinate / overlap / ID
helpers so that no random valid input can produce negative lengths,
end < start, or off-by-one errors.
"""

from __future__ import annotations

from types import SimpleNamespace

import pytest
from hypothesis import HealthCheck, given, settings, strategies as st

from lifton import (
    extract_sequence,
    get_id_fraction,
    lifton_utils,
)


# Bound integers to a sensible genomic range — large enough to
# stress-test arithmetic, small enough to keep examples fast.
COORD = st.integers(min_value=1, max_value=10**9)


# ---------------------------------------------------------------------------
# segments_overlap_length invariants
# ---------------------------------------------------------------------------

@settings(max_examples=200, deadline=None,
          suppress_health_check=[HealthCheck.function_scoped_fixture])
@given(a=COORD, b=COORD, c=COORD, d=COORD)
def test_segments_overlap_length_symmetric(a, b, c, d):
    """Phase 5 bug fix #4 verified: f((s1,e1),(s2,e2)) ==
    f((s2,e2),(s1,e1)) for any valid segments, regardless of tied
    start endpoints."""
    s1, e1 = sorted([a, b])
    s2, e2 = sorted([c, d])
    forward = lifton_utils.segments_overlap_length((s1, e1), (s2, e2))
    reverse = lifton_utils.segments_overlap_length((s2, e2), (s1, e1))
    assert forward == reverse


@settings(max_examples=200, deadline=None)
@given(a=COORD, b=COORD, c=COORD, d=COORD)
def test_segments_overlap_length_consistency(a, b, c, d):
    """ovp is True iff ovp_len > 0."""
    s1, e1 = sorted([a, b])
    s2, e2 = sorted([c, d])
    ovp_len, ovp = lifton_utils.segments_overlap_length((s1, e1), (s2, e2))
    assert (ovp_len > 0) == bool(ovp)


@settings(max_examples=200, deadline=None)
@given(a=COORD, b=COORD)
def test_segment_with_itself_full_overlap(a, b):
    s, e = sorted([a, b])
    ovp_len, ovp = lifton_utils.segments_overlap_length((s, e), (s, e))
    assert ovp is True
    assert ovp_len == e - s + 1


# ---------------------------------------------------------------------------
# merge_children_intervals invariants
# ---------------------------------------------------------------------------

@st.composite
def interval_list(draw):
    n = draw(st.integers(min_value=0, max_value=15))
    out = []
    for _ in range(n):
        a = draw(COORD)
        b = draw(COORD)
        s, e = sorted([a, b])
        out.append(SimpleNamespace(start=s, end=e))
    return out


@settings(max_examples=200, deadline=None)
@given(intervals=interval_list())
def test_merge_children_intervals_no_negative_length(intervals):
    merged = extract_sequence.merge_children_intervals(intervals)
    for s, e in merged:
        assert s <= e, f"merged interval has start > end: ({s},{e})"


@settings(max_examples=200, deadline=None)
@given(intervals=interval_list())
def test_merge_children_intervals_monotonic_starts(intervals):
    merged = extract_sequence.merge_children_intervals(intervals)
    starts = [s for s, _ in merged]
    assert starts == sorted(starts)


@settings(max_examples=200, deadline=None)
@given(intervals=interval_list())
def test_merge_children_intervals_covers_all_inputs(intervals):
    """Every input interval must be contained in some merged output
    interval."""
    if not intervals:
        return
    merged = extract_sequence.merge_children_intervals(intervals)
    for src in intervals:
        assert any(m_s <= src.start and src.end <= m_e
                   for m_s, m_e in merged), (
            f"input interval ({src.start},{src.end}) not covered by any "
            f"merged interval {merged}"
        )


# ---------------------------------------------------------------------------
# custom_bisect_insert invariant
# ---------------------------------------------------------------------------

class _BisectStub:
    def __init__(self, end):
        self.entry = SimpleNamespace(end=end)


@settings(max_examples=200, deadline=None)
@given(ends=st.lists(st.integers(min_value=1, max_value=10**9),
                     min_size=0, max_size=20))
def test_custom_bisect_insert_keeps_sorted_order(ends):
    bucket = []
    for e in ends:
        lifton_utils.custom_bisect_insert(bucket, _BisectStub(e))
        ordered = [x.entry.end for x in bucket]
        assert ordered == sorted(ordered), (
            f"bucket lost sort order after inserting {e}: {ordered}"
        )


# ---------------------------------------------------------------------------
# get_padding_length invariant
# ---------------------------------------------------------------------------

@settings(max_examples=300, deadline=None)
@given(L=st.integers(min_value=0, max_value=10**6))
def test_padding_yields_codon_multiple(L):
    pad = extract_sequence.get_padding_length(L)
    assert pad in (0, 1, 2)
    assert (L + pad) % 3 == 0


# ---------------------------------------------------------------------------
# get_ID_base idempotence
# ---------------------------------------------------------------------------

@settings(max_examples=200, deadline=None)
@given(stem=st.text(alphabet=st.characters(min_codepoint=33, max_codepoint=126,
                                           blacklist_characters="_\t\n"),
                    min_size=1, max_size=20),
       n=st.integers(min_value=0, max_value=999))
def test_get_id_base_strips_trailing_int_once(stem, n):
    """Phase 5 bug fix #5 verified: get_ID_base is idempotent and never
    reduces a single-component id to the empty string."""
    assume_id = f"{stem}_{n}"
    base = lifton_utils.get_ID_base(assume_id)
    base2 = lifton_utils.get_ID_base(base)
    # The function is idempotent EXCEPT when the stem itself happens to
    # end in `_<int>` — then a second pass strips that too. That is the
    # documented legacy behaviour; we only assert idempotence on stems
    # without a trailing _<int>.
    if not _ends_in_underscore_int(base):
        assert base == base2
    # And it never returns the empty string for a non-empty input
    assert base != ""


@settings(max_examples=50, deadline=None)
@given(n=st.integers(min_value=0, max_value=999))
def test_get_id_base_numeric_only_id_preserved(n):
    """Phase 5 bug fix #5 spot-check: a numeric-only id (no underscore)
    must NOT be reduced to the empty string."""
    assert lifton_utils.get_ID_base(str(n)) == str(n)


def _ends_in_underscore_int(s: str) -> bool:
    parts = s.split("_")
    if len(parts) < 2:
        return False
    try:
        int(parts[-1])
        return True
    except ValueError:
        return False


# ---------------------------------------------------------------------------
# get_AA_id_fraction perfect-match property
# ---------------------------------------------------------------------------

PROTEIN_ALPHABET = st.text(
    alphabet="ACDEFGHIKLMNPQRSTVWYM",
    min_size=1, max_size=80,
)


@settings(max_examples=100, deadline=None)
@given(seq=PROTEIN_ALPHABET)
def test_aa_id_fraction_self_is_perfect(seq):
    matches, length = get_id_fraction.get_AA_id_fraction(seq, seq)
    # The function breaks early if it sees '*' in target; our alphabet
    # excludes '*' so the full string is compared.
    assert matches == len(seq)
    assert length == len(seq)


@settings(max_examples=100, deadline=None)
@given(seq=st.text(alphabet="ACGT", min_size=1, max_size=80))
def test_dna_id_fraction_self_is_perfect(seq):
    matches, length = get_id_fraction.get_DNA_id_fraction(seq, seq)
    assert matches == len(seq)
    assert length == len(seq)


# ---------------------------------------------------------------------------
# IntervalTree-style coordinate sanity
# ---------------------------------------------------------------------------

@settings(max_examples=200, deadline=None)
@given(a=COORD, b=COORD)
def test_overlap_self_for_normalised_segment(a, b):
    """Sorting two random integers always yields start <= end with
    non-negative length."""
    s, e = sorted([a, b])
    length = e - s + 1
    assert length >= 1
    # And the segment must overlap itself
    _, ovp = lifton_utils.segments_overlap_length((s, e), (s, e))
    assert ovp is True
