"""Unit tests for lifton.get_id_fraction — alignment identity helpers."""

from __future__ import annotations

import pytest

from lifton import get_id_fraction


class TestGetAAIdFraction:
    def test_perfect_match(self):
        matches, length = get_id_fraction.get_AA_id_fraction("MAGT", "MAGT")
        assert matches == 4
        assert length == 4

    def test_lowercase_normalised(self):
        matches, length = get_id_fraction.get_AA_id_fraction("magt", "MAGT")
        assert matches == 4 and length == 4

    def test_gaps_in_reference_shorten_total_length(self):
        # one ref gap is excluded from the denominator
        matches, length = get_id_fraction.get_AA_id_fraction("MA-GT", "MAQGT")
        assert length == len("MAQGT") - 1
        assert matches == 4  # M, A, G, T match (ref '-' vs query 'Q' is mismatch)

    def test_empty_pair_returns_one_for_length(self):
        matches, length = get_id_fraction.get_AA_id_fraction("", "")
        assert matches == 0 and length == 1

    def test_target_stop_breaks_iteration(self):
        # iteration breaks when target sees '*'
        matches, length = get_id_fraction.get_AA_id_fraction("MAGT", "M*XX")
        # only 'M' was matched before break
        assert matches == 1


class TestGetPartialIdFraction:
    def test_simple_window(self):
        matches, length = get_id_fraction.get_partial_id_fraction(
            "MAGTQQQ", "MAGTAAA", 0, 4
        )
        assert matches == 4 and length == 4

    def test_zero_length_window(self):
        matches, length = get_id_fraction.get_partial_id_fraction(
            "MAGT", "MAGT", 2, 2
        )
        assert length == 1  # protected from /0


class TestGetDNAIdFraction:
    def test_perfect_match(self):
        m, n = get_id_fraction.get_DNA_id_fraction("ACGT", "ACGT")
        assert m == 4 and n == 4

    def test_partial_match(self):
        m, n = get_id_fraction.get_DNA_id_fraction("ACGT", "ACCT")
        assert m == 3 and n == 4

    def test_empty_pair(self):
        m, n = get_id_fraction.get_DNA_id_fraction("", "")
        assert m == 0 and n == 1
