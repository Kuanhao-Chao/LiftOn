"""Tier-3 perf fix — vectorized identity counters must be EXACTLY equivalent.

`lifton/get_id_fraction.py` runs after every alignment (~6-12 scans per transcript with
the best-of-outcome merge + candidate-3 + ORF rescue) and was pure per-character Python,
measuring 8-25% of the alignment+scoring pair and holding the GIL, which also capped
`--threads N` scaling. The functions are now expressed with C-level primitives.

These tests fuzz the new implementations against LITERAL COPIES of the original loops --
including the awkward parts: the `*` scan stops AFTER the position it is found at and
only within the visited window; gaps are counted over that truncated window; but
`total_length` uses the FULL requested span.
"""
import random
import string

import pytest

from lifton import get_id_fraction


# ---------------------------------------------------------------------------
# Literal copies of the pre-vectorization implementations (the oracle).
# ---------------------------------------------------------------------------

def _orig_partial(reference, target, start, end):
    matches = 0
    reference = reference.upper()
    target = target.upper()
    gaps_in_ref = 0
    for i, letter in enumerate(reference[start:end]):
        if letter == '-':
            gaps_in_ref += 1
        if letter == target[i + start]:
            matches += 1
        if target[i + start] == "*":
            break
    total_length = (end - start) - gaps_in_ref
    if total_length == 0:
        return matches, 1
    return matches, total_length


def _orig_aa(reference, target):
    matches = 0
    reference = reference.upper()
    target = target.upper()
    gaps_in_ref = 0
    for i, letter in enumerate(reference):
        if letter == '-':
            gaps_in_ref += 1
        if letter == target[i]:
            matches += 1
        if target[i] == "*":
            break
    if max(len(reference), len(target)) == 0:
        return matches, 1
    total_length = max(len(reference), len(target)) - gaps_in_ref
    if total_length <= 0:
        return matches, 1
    return matches, total_length


def _orig_dna(reference, target):
    matches = 0
    reference = reference.upper()
    target = target.upper()
    if len(reference) != len(target):
        raise ValueError("length mismatch")
    for i, letter in enumerate(reference):
        if letter == target[i]:
            matches += 1
    if max(len(reference), len(target)) == 0:
        return matches, 1
    return matches, max(len(reference), len(target))


# ---------------------------------------------------------------------------
# Fuzz corpora — heavy on gaps, stops and case, which is where the edges live.
# ---------------------------------------------------------------------------

AA_ALPHABET = "ACDEFGHIKLMNPQRSTVWY" + "-" * 6 + "*" * 3 + "acdefg"
DNA_ALPHABET = "ACGTN" + "-" * 3 + "acgtn"


def _pair(rng, alphabet, length):
    return ("".join(rng.choice(alphabet) for _ in range(length)),
            "".join(rng.choice(alphabet) for _ in range(length)))


@pytest.mark.parametrize("seed", range(40))
def test_aa_matches_original(seed):
    rng = random.Random(seed)
    for length in (0, 1, 2, 5, 17, 64, 257):
        ref, tgt = _pair(rng, AA_ALPHABET, length)
        assert get_id_fraction.get_AA_id_fraction(ref, tgt) == _orig_aa(ref, tgt), \
            (ref, tgt)


@pytest.mark.parametrize("seed", range(40))
def test_dna_matches_original(seed):
    rng = random.Random(seed)
    for length in (0, 1, 2, 5, 17, 64, 257):
        ref, tgt = _pair(rng, DNA_ALPHABET, length)
        assert get_id_fraction.get_DNA_id_fraction(ref, tgt) == _orig_dna(ref, tgt), \
            (ref, tgt)


@pytest.mark.parametrize("seed", range(40))
def test_partial_matches_original(seed):
    rng = random.Random(seed)
    for length in (1, 2, 5, 17, 64, 129):
        ref, tgt = _pair(rng, AA_ALPHABET, length)
        for _ in range(4):
            start = rng.randrange(0, length)
            end = rng.randrange(start, length + 1)
            assert (get_id_fraction.get_partial_id_fraction(ref, tgt, start, end)
                    == _orig_partial(ref, tgt, start, end)), (ref, tgt, start, end)


# ---------------------------------------------------------------------------
# Explicit edge cases (documented semantics, not just fuzz luck).
# ---------------------------------------------------------------------------

def test_stop_codon_truncates_after_its_own_position():
    ref, tgt = "MKKKK", "MK*KK"
    assert get_id_fraction.get_AA_id_fraction(ref, tgt) == _orig_aa(ref, tgt)
    # positions 0,1 match and position 2 is processed then breaks -> 2 matches
    assert get_id_fraction.get_AA_id_fraction(ref, tgt)[0] == 2


def test_gaps_counted_only_up_to_the_stop():
    ref, tgt = "--M--", "AB*AB"
    assert get_id_fraction.get_AA_id_fraction(ref, tgt) == _orig_aa(ref, tgt)


def test_all_gap_reference_returns_bounded_denominator():
    matches, total = get_id_fraction.get_AA_id_fraction("-----", "AAAAA")
    assert (matches, total) == _orig_aa("-----", "AAAAA")
    assert total == 1


def test_partial_total_length_uses_the_requested_span_not_the_truncation():
    ref, tgt = "MKKKKKK", "MK*KKKK"
    assert (get_id_fraction.get_partial_id_fraction(ref, tgt, 0, 7)
            == _orig_partial(ref, tgt, 0, 7))


def test_case_insensitivity():
    assert (get_id_fraction.get_DNA_id_fraction("acgt", "ACGT")
            == _orig_dna("acgt", "ACGT") == (4, 4))


def test_dna_length_mismatch_still_raises():
    with pytest.raises(ValueError):
        get_id_fraction.get_DNA_id_fraction("ACGT", "ACG")
