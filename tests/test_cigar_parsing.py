"""`coreutils.parse_cigar_items` replaces the third-party `cigar` package.

`cigar` 0.1.3 is the only dependency LiftOn had that ships an sdist which cannot
be built in a modern environment: its `setup.py` calls
`ez_setup.use_setuptools()`, which downloads a setuptools egg from a URL that no
longer serves one, so a clean `pip install lifton` died with
`tarfile.ReadError: not a gzip file`. LiftOn used exactly one call from it.

These tests pin the replacement against the behaviour it replaces. Where the
`cigar` package is still importable they compare output directly; otherwise they
assert the same expectations the package's own doctests declare.
"""

from __future__ import annotations

import pytest

from lifton.coreutils import parse_cigar_items


CIGARS = [
    "100M",
    "20H20M20S",
    "10M20S10M",
    "1S1H1S5H1S5M10H",
    "50M2I30M",
    "10=2X20=",
    "5S20M3D10M5H",
    "100M500N100M",
    "1M",
    "999999M",
    "3=1X3=1X10=",
    "12M3I4D5N6S7H8P9=10X",
]


class TestMatchesTheReplacedPackage:
    @pytest.mark.parametrize("cigar_string", CIGARS)
    def test_identical_to_cigar_package(self, cigar_string):
        Cigar = pytest.importorskip("cigar").Cigar
        assert list(parse_cigar_items(cigar_string)) == list(
            Cigar(cigar_string).items(),
        )


class TestParsing:
    def test_single_operation(self):
        assert list(parse_cigar_items("100M")) == [(100, "M")]

    def test_multiple_operations_keep_order(self):
        assert list(parse_cigar_items("20H20M20S")) == [
            (20, "H"), (20, "M"), (20, "S"),
        ]

    def test_repeated_operation_is_not_merged(self):
        assert list(parse_cigar_items("10M20S10M")) == [
            (10, "M"), (20, "S"), (10, "M"),
        ]

    def test_eqx_operations(self):
        """parasail emits `=`/`X`, which is what LiftOn actually parses."""
        assert list(parse_cigar_items("10=2X20=")) == [
            (10, "="), (2, "X"), (20, "="),
        ]

    def test_multi_digit_lengths(self):
        assert list(parse_cigar_items("999999M")) == [(999999, "M")]

    def test_empty_string_yields_nothing(self):
        assert list(parse_cigar_items("")) == []

    def test_unavailable_alignment(self):
        """`*` is the SAM "no alignment" CIGAR; upstream yields (0, None)."""
        assert list(parse_cigar_items("*")) == [(0, None)]

    def test_star_does_not_raise(self):
        """Upstream ran `raise StopIteration` here.

        PEP 479 turns that into a RuntimeError inside a generator on Python
        >= 3.7, so the package LiftOn depended on was already broken for this
        input on every supported interpreter.
        """
        list(parse_cigar_items("*"))  # must not raise

    def test_is_a_generator_not_a_list(self):
        items = parse_cigar_items("10M10I")
        assert next(items) == (10, "M")
        assert next(items) == (10, "I")
        with pytest.raises(StopIteration):
            next(items)

    def test_total_length_is_preserved(self):
        assert sum(n for n, _ in parse_cigar_items("12M3I4D5N6S7H8P9=10X")) == 64


class TestCigarIsNoLongerImported:
    def test_align_module_does_not_import_cigar(self):
        from pathlib import Path

        source = (
            Path(__file__).resolve().parents[1] / "lifton" / "align.py"
        ).read_text(encoding="utf-8")
        assert "from cigar import" not in source
        assert "import cigar" not in source
