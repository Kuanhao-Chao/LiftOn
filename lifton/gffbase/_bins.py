# ---------------------------------------------------------------------------
# Author: Kuan-Hao Chao <kuanhao.chao@gmail.com>
# ---------------------------------------------------------------------------
"""UCSC genomic binning. Used only for the legacy SQLite export path; the
runtime query layer uses DuckDB's R-tree (or the seqstart B-tree fallback)
and never touches `bin`. Port of the relevant logic from `gffutils/bins.py`.
"""

from __future__ import annotations

_BINOFFSETS = (512 + 64 + 8 + 1, 64 + 8 + 1, 8 + 1, 1, 0)
_BINFIRSTSHIFT = 17
_BINNEXTSHIFT = 3


def bin_from_coords(start: int, end: int) -> int:
    """Return the smallest UCSC bin that contains [start, end] (1-based, inclusive).
    `start`/`end` here are converted to 0-based half-open as the original code
    does. Used only by the SQLite export to populate the legacy `bin` column.
    """
    start0 = start - 1
    end0 = end
    start_bin = start0 >> _BINFIRSTSHIFT
    end_bin = (end0 - 1) >> _BINFIRSTSHIFT
    for offset in _BINOFFSETS:
        if start_bin == end_bin:
            return offset + start_bin
        start_bin >>= _BINNEXTSHIFT
        end_bin >>= _BINNEXTSHIFT
    return 0
