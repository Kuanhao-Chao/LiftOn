# ---------------------------------------------------------------------------
# Author: Kuan-Hao Chao <kuanhao.chao@gmail.com>
# ---------------------------------------------------------------------------
"""Exception classes — verbatim port of legacy `gffutils.exceptions`.

Same names, same constructors, same attributes. Downstream code that catches
`gffutils.FeatureNotFoundError` works against `gffbase.FeatureNotFoundError`.
"""

from __future__ import annotations


class GFFFormatError(ValueError):
    """Raised when a GFF3 line violates the spec.

    Inherits from ``ValueError`` so legacy code that catches
    ``ValueError`` keeps working — but the dedicated subclass carries
    structured fields ``line_no``, ``kind``, and ``message`` that point
    callers straight to the offender in the input.

    Phase 16: when the Rust extension is loaded, the canonical class is
    ``gffbase._native.GFFFormatError`` (a PyO3 ``create_exception!``
    type). At import time, ``gffbase.__init__`` rebinds the public name
    to whichever class is actually live so ``isinstance`` and
    ``except`` clauses keep working regardless of which path raised.
    """

    def __init__(self, message: str = "", *, line_no: int = 0, kind: str = ""):
        super().__init__(message)
        self.message = message
        self.line_no = line_no
        self.kind = kind


class FeatureNotFoundError(Exception):
    """Raised by ``FeatureDB.__getitem__`` when the requested ID is absent."""

    def __init__(self, feature_id: str):
        Exception.__init__(self, f"feature not found: {feature_id}")
        self.feature_id = feature_id


class DuplicateIDError(Exception):
    """Raised during ingestion when a duplicate ID is encountered with
    ``merge_strategy='error'`` (Phase 6 will wire `merge_strategy`)."""


class AttributeStringError(Exception):
    """Raised on malformed col-9 attributes."""


class EmptyInputError(Exception):
    """Raised when an input file or iterable yields no features."""
