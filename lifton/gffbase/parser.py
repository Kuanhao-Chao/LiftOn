# ---------------------------------------------------------------------------
# Author: Kuan-Hao Chao <kuanhao.chao@gmail.com>
# ---------------------------------------------------------------------------
"""Public parser entry points. Dispatches to the Rust extension when available,
falls back to the pure-Python implementation otherwise. The two implementations
are required by tests to produce identical output.
"""

from __future__ import annotations

from typing import Iterator, List, Optional

from .feature import ParsedFeature
from ._pyfallback import parser as _pyparser

try:  # pragma: no cover - import availability is env-dependent
    from . import _native as _rust  # type: ignore[attr-defined]

    _NATIVE = True
except ImportError:
    _rust = None
    _NATIVE = False


def native_available() -> bool:
    """True if the compiled extension is importable."""
    return _NATIVE


class _Iterator:
    """Adapter: wraps either the Rust iterator (yielding 11-tuples) or the
    pure-Python iterator (yielding ParsedFeature) and always yields
    ParsedFeature.

    Phase 16: also exposes ``.warnings`` — a list of structured-error
    dicts (``line_no``, ``kind``, ``message``) collected when the
    iterator was created with ``strict=False``.
    """

    __slots__ = ("_inner", "_native")

    def __init__(self, inner, native: bool):
        self._inner = inner
        self._native = native

    def __iter__(self):
        return self

    def __next__(self) -> ParsedFeature:
        item = next(self._inner)
        if self._native:
            return ParsedFeature.from_tuple(item)
        return item

    def dialect(self) -> dict:
        return self._inner.dialect()

    def directives(self) -> list:
        return list(self._inner.directives())

    @property
    def warnings(self) -> List[dict]:
        """Errors that were demoted to warnings during a non-strict run.

        Each item is a dict with keys ``line_no``, ``kind``, and
        ``message``. Empty for strict runs (errors raise instead).
        """
        if hasattr(self._inner, "warnings"):
            w = self._inner.warnings
            return list(w() if callable(w) else w)
        return []


def _resolve_engine(engine: Optional[str]) -> str:
    if engine is None or engine == "auto":
        return "rust" if _NATIVE else "python"
    if engine == "rust" and not _NATIVE:
        raise RuntimeError(
            "Rust extension not built. Run `maturin develop --release` or pass engine='python'."
        )
    if engine not in {"rust", "python"}:
        raise ValueError(f"unknown engine '{engine}'")
    return engine


def parse_gff(
    path: str,
    *,
    checklines: int = 10,
    force_dialect_check: bool = False,
    force_gff: bool = False,
    strict: bool = True,
    engine: Optional[str] = "auto",
) -> _Iterator:
    """Parse a GFF3/GTF file (plain text or ``.gz``).

    Returns an iterator of ``ParsedFeature`` plus ``.dialect()``,
    ``.directives()``, and (Phase 16) ``.warnings`` accessors.

    Parameters
    ----------
    strict : bool
        When True (default), the iterator raises ``GFFFormatError`` on
        the first malformed line. When False, malformed lines are
        skipped silently and recorded in ``iterator.warnings``.
    """
    eng = _resolve_engine(engine)
    if eng == "rust":
        it = _rust.parse_file(  # type: ignore[union-attr]
            path,
            checklines=checklines,
            force_dialect_check=force_dialect_check,
            force_gff=force_gff,
            strict=strict,
        )
        return _Iterator(it, native=True)
    it = _pyparser.parse_file(
        path,
        checklines=checklines,
        force_dialect_check=force_dialect_check,
        force_gff=force_gff,
        strict=strict,
    )
    return _Iterator(it, native=False)


def parse_bytes(
    data: bytes,
    *,
    checklines: int = 10,
    force_dialect_check: bool = False,
    force_gff: bool = False,
    strict: bool = True,
    engine: Optional[str] = "auto",
) -> _Iterator:
    eng = _resolve_engine(engine)
    if eng == "rust":
        it = _rust.parse_bytes(  # type: ignore[union-attr]
            data,
            checklines=checklines,
            force_dialect_check=force_dialect_check,
            force_gff=force_gff,
            strict=strict,
        )
        return _Iterator(it, native=True)
    it = _pyparser.parse_bytes(
        data,
        checklines=checklines,
        force_dialect_check=force_dialect_check,
        force_gff=force_gff,
        strict=strict,
    )
    return _Iterator(it, native=False)


def detect_dialect(
    path: str, *, checklines: int = 10, engine: Optional[str] = "auto"
) -> dict:
    eng = _resolve_engine(engine)
    if eng == "rust":
        return _rust.detect_dialect(path, checklines=checklines)  # type: ignore[union-attr]
    return _pyparser.detect_dialect(path, checklines=checklines)
