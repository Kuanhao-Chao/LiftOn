# ---------------------------------------------------------------------------
# Author: Kuan-Hao Chao <kuanhao.chao@gmail.com>
# ---------------------------------------------------------------------------
"""Compatibility helpers."""

from __future__ import annotations

from pathlib import Path


def example_filename(fn: str) -> str:
    """Return the absolute path to a packaged example file under ``tests/data``.

    Mirrors the legacy ``gffutils.example_filename`` API. The legacy package
    shipped fixtures inside the package itself; we point at the test fixtures
    directory bundled with the source distribution.
    """
    here = Path(__file__).resolve().parent
    candidates = [
        here.parent.parent / "tests" / "data" / fn,
        here / "data" / fn,
    ]
    for c in candidates:
        if c.is_file():
            return str(c)
    raise FileNotFoundError(f"example file not found: {fn}")
