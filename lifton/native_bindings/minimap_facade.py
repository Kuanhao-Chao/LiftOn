# ---------------------------------------------------------------------------
# Phase 10 — `mappy` (Python bindings to minimap2) facade.
#
# Real PyPI binding wrapped in a thread-safe class. The Aligner's `map()`
# method is documented to release the GIL during the C-level alignment
# kernel, so multiple threads can drive a single Aligner concurrently.
#
# Hermetic CI: when mappy is not installed (e.g. in a stripped-down
# environment), `is_mappy_available()` returns False and consumers
# must fall back to the subprocess path. The facade itself never
# spawns a subprocess — that's the legacy align_features.py.
# ---------------------------------------------------------------------------

from __future__ import annotations

import threading
from typing import Iterator, List, Optional, Tuple

from .types import MinimapHit

# minimap2 feature flag for extended (`=`/`X`) CIGAR ops — the C constant
# `MM_F_EQX` that `--eqx` sets. mappy ORs `extra_flags` into the option set,
# so passing this makes `Aligner.map()` emit `=`/`X` instead of `M`.
MM_F_EQX = 0x4000000


def is_mappy_available() -> bool:
    try:
        import mappy  # noqa: F401
        return True
    except ImportError:
        return False


def _translate_mm2_options(mm2_options: str, best_n: int) -> Tuple[int, int]:
    """Translate the subset of Liftoff's ``mm2_options`` string that mappy's
    ``Aligner`` constructor can honour into ``(extra_flags, best_n)``.

    Liftoff's default is ``-a --end-bonus 5 --eqx -N 50 -p 0.5``. Critically,
    **``--eqx`` is mandatory**: the downstream Liftoff CIGAR parser
    (``align_features.get_cigar_operations`` → match=7/mismatch=8) only counts
    ``=``/``X`` ops, so without ``MM_F_EQX`` mappy's default ``M`` CIGAR parses
    to *zero* aligned bases and every feature is reported unmapped (the
    fresh-``--native`` "maps nothing" bug). ``-N`` maps to ``best_n``. ``-a``
    (SAM output) is irrelevant to mappy's object API. ``-p`` (secondary/primary
    score ratio) and ``--end-bonus`` are **not exposed** by mappy's constructor,
    so they cannot be replicated here — they only fine-tune copy/secondary
    sensitivity (``-copies``) and end alignment, not whether a feature maps;
    the small resulting difference from the subprocess path is acceptable
    (Liftoff is itself non-deterministic across runs).
    """
    extra_flags = 0
    toks = [t for t in mm2_options.split() if t]
    i = 0
    while i < len(toks):
        t = toks[i]
        if t == "--eqx":
            extra_flags |= MM_F_EQX
        elif t == "-N" and i + 1 < len(toks):
            try:
                best_n = int(toks[i + 1])
            except ValueError:
                pass
            i += 1
        i += 1
    return extra_flags, best_n


class MinimapAligner:
    """Wraps :class:`mappy.Aligner` with a thread-safe `map(query)` that
    yields :class:`MinimapHit` records (the subset of fields LiftOn
    consumes). Construction is parameter-equivalent to the
    ``minimap2`` CLI: ``preset`` matches ``-x splice`` etc., and
    ``mm2_options`` is the same string args the legacy code used to
    pass via ``args.mm2_options``."""

    def __init__(
        self,
        target_fa: str,
        *,
        preset: Optional[str] = None,
        threads: int = 1,
        mm2_options: str = "",
        best_n: int = 50,
    ):
        if not is_mappy_available():
            raise RuntimeError(
                "mappy is not installed; install it via "
                "`pip install mappy` or `conda install -c bioconda mappy`."
            )
        import mappy
        # Translate Liftoff's mm2_options into the mappy-honourable subset.
        # MM_F_EQX (`--eqx`) is load-bearing: without it mappy emits `M`
        # CIGAR and the downstream Liftoff parser counts zero aligned bases,
        # mapping nothing (the fresh-`--native` empty-output bug).
        extra_flags, best_n = _translate_mm2_options(mm2_options, best_n)
        # mappy.Aligner is thread-safe for `.map(seq)`; storing a
        # single instance per process is the recommended pattern.
        self._aligner = mappy.Aligner(
            fn_idx_in=target_fa,
            preset=preset,
            n_threads=threads,
            best_n=best_n,
            extra_flags=extra_flags,
        )
        if not self._aligner:
            raise RuntimeError(
                f"mappy.Aligner failed to open target FASTA: {target_fa}"
            )
        self._lock = threading.Lock()  # only used in tests / debug

    @property
    def n_seq(self) -> int:
        return self._aligner.n_seq

    def seq_names(self) -> List[str]:
        return list(self._aligner.seq_names)

    def map(self, query_name: str, query_seq: str) -> Iterator[MinimapHit]:
        """Yield one :class:`MinimapHit` per alignment. Wraps
        :meth:`mappy.Aligner.map`, which releases the GIL during the
        C-level alignment kernel and therefore scales well under
        :class:`concurrent.futures.ThreadPoolExecutor`."""
        # mappy.Aligner.map yields Alignment objects with attributes
        # ctg, r_st, r_en, q_st, q_en, strand, mapq, NM, cigar_str,
        # is_primary. We project to a stable dataclass.
        for hit in self._aligner.map(query_seq):
            yield MinimapHit(
                query_name=query_name,
                ctg=hit.ctg,
                r_st=int(hit.r_st),
                r_en=int(hit.r_en),
                q_st=int(hit.q_st),
                q_en=int(hit.q_en),
                strand=int(hit.strand),
                mapq=int(hit.mapq),
                NM=int(getattr(hit, "NM", 0)),
                cigar_str=str(hit.cigar_str),
                is_primary=bool(getattr(hit, "is_primary", True)),
            )

    def map_one_best(self, query_name: str, query_seq: str
                      ) -> Optional[MinimapHit]:
        """Convenience: return the highest-MAPQ primary hit, or None."""
        best: Optional[MinimapHit] = None
        for h in self.map(query_name, query_seq):
            if not h.is_primary:
                continue
            if best is None or h.mapq > best.mapq:
                best = h
        return best
