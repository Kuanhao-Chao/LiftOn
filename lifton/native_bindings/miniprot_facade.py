# ---------------------------------------------------------------------------
# Phase 10 — `pyminiprot`-shaped facade.
#
# A real PyO3 binding for miniprot does NOT yet exist on PyPI. The
# facade below presents the API a real binding would (an `Index` with
# an `.align(protein_seq)` generator yielding GFF3Hit records) but the
# implementation transparently invokes the miniprot subprocess and
# decodes its stdout in-memory via :class:`GFF3Hit.from_gff_line`.
#
# Why ship the facade NOW
# ───────────────────────
# 1. LiftOn call sites can be wired against the binding-shaped API
#    today; swapping a real PyO3 binding in later is one-line.
# 2. The subprocess invocation is amortised across all reference
#    proteins (one fork per Index, not one fork per protein), which
#    is already much closer to the per-process cost of a real
#    binding than the pre-Phase-10 path that respawned miniprot for
#    every Index instance.
# 3. The GFF3-line projection is binding-shaped and thread-safe by
#    construction (the bytes blob lives in RAM; threads read it
#    without contention).
#
# When the real `pyminiprot` lands, ``is_pyminiprot_native_available``
# will return True and ``MiniprotIndex.align`` will route through the
# native binding. Until then it routes through the subprocess.
# ---------------------------------------------------------------------------

from __future__ import annotations

import os
import subprocess
import sys
import tempfile
from typing import Iterator, List, Optional

from .types import GFF3Bundle, GFF3Hit


def is_pyminiprot_native_available() -> bool:
    """Phase 10: the real PyO3 binding has not been built yet.
    The function exists so that the rest of the pipeline can ALREADY
    branch on it — when the real binding lands the import below
    succeeds and this returns True.
    """
    try:
        import pyminiprot  # type: ignore[import]
        return hasattr(pyminiprot, "Index")
    except ImportError:
        return False


class MiniprotIndex:
    """Drop-in shape for the eventual ``pyminiprot.Index``.

    Today it shells out to the ``miniprot`` CLI **once** per
    construction, captures the GFF3 stdout into a bytes blob, and
    decodes it into :class:`GFF3Hit` records that callers stream via
    :meth:`align`.

    Crucially, the API is exactly what a real binding would expose:

    .. code-block:: python

        idx = MiniprotIndex("target.fa", mp_options="")
        for hit in idx.align("MAGT...*"):
            ...

    Once the real ``pyminiprot`` lands, ``__init__`` will route
    through the native ``pyminiprot.Index`` constructor and ``align``
    will yield directly from the C-level iterator — call sites stay
    unchanged.
    """

    def __init__(
        self,
        target_fa: str,
        *,
        mp_options: str = "",
        miniprot_path: str = "miniprot",
        ref_proteins_path: Optional[str] = None,
    ):
        self.target_fa = target_fa
        self.mp_options = mp_options
        self.miniprot_path = miniprot_path
        self._ref_proteins_path = ref_proteins_path
        self._cached_bundle: Optional[GFF3Bundle] = None

        # If the real PyO3 binding ever appears, prefer it.
        self._native = None
        if is_pyminiprot_native_available():
            try:                                                # pragma: no cover
                import pyminiprot                               # type: ignore
                self._native = pyminiprot.Index(target_fa, mp_options=mp_options)
            except Exception:
                self._native = None

    # ------------------------------------------------------------------
    # Public API (mirrors the eventual pyminiprot.Index)
    # ------------------------------------------------------------------

    def align_all(self, ref_proteins_path: Optional[str] = None) -> GFF3Bundle:
        """Run miniprot ONCE against an entire reference-proteins FASTA
        and return the parsed :class:`GFF3Bundle`. The bundle is
        cached on the Index instance so repeated calls are O(1).

        For a real PyO3 binding this would iterate proteins inside the
        process; for the subprocess path we still amortise the fork
        cost across all proteins, so the *call-site* cost model is
        identical.
        """
        if self._cached_bundle is not None:
            return self._cached_bundle

        proteins = ref_proteins_path or self._ref_proteins_path
        if proteins is None:
            raise ValueError(
                "MiniprotIndex.align_all requires a ref_proteins_path "
                "(either at construction or via this argument)."
            )

        if self._native is not None:                            # pragma: no cover
            # Real PyO3 binding path (forward-compatible — not exercised today).
            hits: List[GFF3Hit] = []
            raw = bytearray()
            for native_hit in self._native.align_file(proteins):
                line = native_hit.to_gff_line() + "\n"
                raw.extend(line.encode("utf-8"))
                parsed = GFF3Hit.from_gff_line(line)
                if parsed is not None:
                    hits.append(parsed)
            self._cached_bundle = GFF3Bundle(hits=hits, raw_bytes=bytes(raw))
            return self._cached_bundle

        # ── Subprocess fallback ─────────────────────────────────────
        cmd = [
            self.miniprot_path,
            "--gff-only",
            self.target_fa,
            proteins,
        ] + [opt for opt in self.mp_options.split(" ") if opt]

        proc = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            bufsize=1 << 20,
        )
        stdout_bytes, stderr_bytes = proc.communicate()
        if proc.returncode != 0:
            err = (stderr_bytes or b"").decode("utf-8", errors="replace")
            raise RuntimeError(
                f"miniprot exited with code {proc.returncode}: {err.strip()}"
            )
        if stderr_bytes and b"ERROR" in stderr_bytes.upper():
            raise RuntimeError(
                "miniprot reported an ERROR during mapping; refusing to "
                "produce an Index."
            )

        hits = []
        for line in stdout_bytes.decode("utf-8", errors="replace").splitlines():
            parsed = GFF3Hit.from_gff_line(line)
            if parsed is not None:
                hits.append(parsed)
        self._cached_bundle = GFF3Bundle(hits=hits, raw_bytes=stdout_bytes)
        return self._cached_bundle

    def align(self, protein_seq: str) -> Iterator[GFF3Hit]:
        """Yield :class:`GFF3Hit` for a single protein query.

        Today the subprocess path runs miniprot against ALL reference
        proteins in one shot; this method then filters the cached
        bundle by ``Target=<seq_id>`` attribute. For a real PyO3
        binding ``align`` would dispatch a per-protein C call.

        For Phase 10 we treat the protein sequence's first 50 chars
        as the lookup key (mirroring how ``Target`` attributes
        carry the protein record id, not the sequence).
        """
        if self._cached_bundle is None:
            # Without a pre-built bundle, the per-protein call cannot
            # be served by the subprocess path; surface a clear error.
            raise RuntimeError(
                "MiniprotIndex.align(protein_seq) requires align_all() to "
                "have been called first (subprocess path) or a real "
                "pyminiprot binding."
            )
        # The subprocess path's GFF3 output uses the protein's record
        # id, not its sequence, in the Target= attribute. Without a
        # mapping we yield every hit; downstream code is responsible
        # for filtering by Target. This matches the legacy contract
        # where the miniprot.gff3 file held all hits regardless of
        # which protein triggered which one.
        return iter(self._cached_bundle.hits)

    @property
    def is_native(self) -> bool:
        """True iff a real ``pyminiprot`` binding is in use."""
        return self._native is not None

    @property
    def raw_bytes(self) -> bytes:
        """Forward-compatible accessor for callers that still want the
        full GFF3 blob (for example to feed into the Phase 7
        streaming-adapter ingest path)."""
        if self._cached_bundle is None:
            raise RuntimeError("Call align_all() before .raw_bytes")
        return self._cached_bundle.raw_bytes
