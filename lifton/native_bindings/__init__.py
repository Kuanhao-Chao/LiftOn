# ---------------------------------------------------------------------------
# Phase 10 — native binding facade.
#
# Two backends, one shape:
#   * minimap2 → mappy (real PyO3 binding, on PyPI today)
#   * miniprot → pyminiprot (real native binding does NOT yet exist;
#     until it does, the facade transparently falls back to the
#     subprocess path and exposes the same Index/align surface a real
#     binding would, so call sites can be wired once and the binding
#     can be swapped in later.)
# ---------------------------------------------------------------------------

from .minimap_facade import (
    MinimapAligner,
    is_mappy_available,
)
from .miniprot_facade import (
    MiniprotIndex,
    is_pyminiprot_native_available,
)
from .types import GFF3Bundle, GFF3Hit, MinimapHit

__all__ = [
    "MinimapAligner",
    "MiniprotIndex",
    "GFF3Bundle",
    "GFF3Hit",
    "MinimapHit",
    "is_mappy_available",
    "is_pyminiprot_native_available",
]
