# ---------------------------------------------------------------------------
# Author: Kuan-Hao Chao <kuanhao.chao@gmail.com>
# ---------------------------------------------------------------------------
"""gffbase — modernized successor to gffutils.

Phase 5: full drop-in public API surface (FeatureDB, Feature, create_db,
DataIterator, GFFWriter, exceptions) on top of the Phase 4 DuckDB ingestion
engine.
"""

from .exceptions import (
    AttributeStringError,
    DuplicateIDError,
    EmptyInputError,
    FeatureNotFoundError,
    GFFFormatError as _PyGFFFormatError,
)

# Phase 16: prefer the Rust-defined exception class when the extension
# is loaded — that's the type Rust will actually raise. Fall back to
# the pure-Python definition otherwise. Both inherit from `ValueError`
# so legacy `pytest.raises(ValueError)` callers keep working.
try:  # pragma: no cover — import-time branch
    from ._native import GFFFormatError  # type: ignore[attr-defined]
except ImportError:
    GFFFormatError = _PyGFFFormatError  # type: ignore[assignment]
from .feature import Feature, ParsedFeature
from .parser import parse_gff, parse_bytes, detect_dialect, native_available
from . import ingest
from . import merge_criteria
from .helpers import example_filename
from .gffwriter import GFFWriter
from .iterators import DataIterator
from .interface import FeatureDB
from .create_db import create_db
from .sqlite_export import export_sqlite

__all__ = [
    # Drop-in legacy surface
    "create_db",
    "FeatureDB",
    "Feature",
    "DataIterator",
    "GFFWriter",
    "example_filename",
    "FeatureNotFoundError",
    "DuplicateIDError",
    "AttributeStringError",
    "EmptyInputError",
    "GFFFormatError",
    "merge_criteria",
    # gffbase extras
    "ParsedFeature",
    "parse_gff",
    "parse_bytes",
    "detect_dialect",
    "native_available",
    "ingest",
    "export_sqlite",
    "__version__",
]

__version__ = "0.1.0"
