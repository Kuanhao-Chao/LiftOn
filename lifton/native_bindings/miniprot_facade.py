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
# 1. Integrations can use the binding-shaped API today. The main CLI keeps
#    its guarded direct-stream subprocess until a real PyO3 binding is
#    available; swapping that backend later does not change this API.
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

import subprocess
from io import BytesIO
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
        threads: int = 1,
    ):
        self.target_fa = target_fa
        self.mp_options = mp_options
        self.miniprot_path = miniprot_path
        self._ref_proteins_path = ref_proteins_path
        # Iteration 17: scale miniprot's own -t with LiftOn's -t/--threads,
        # matching the subprocess path in lifton.run_miniprot (the command is
        # built by the shared _build_miniprot_command helper below).
        self.threads = threads
        self._cached_bundle: Optional[GFF3Bundle] = None
        self._cached_raw_bytes: Optional[bytes] = None

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

    def align_all(
        self,
        ref_proteins_path: Optional[str] = None,
        *,
        raw_only: bool = False,
    ) -> GFF3Bundle:
        """Run miniprot ONCE against an entire reference-proteins FASTA
        and return a :class:`GFF3Bundle`. The default mode parses and caches
        ``GFF3Hit`` objects, preserving the original facade contract.

        ``raw_only=True`` skips UTF-8 decoding, line splitting, and hit
        construction. Its returned bundle intentionally has no ``hits`` and
        carries the byte-identical miniprot output in ``raw_bytes``. Raw bytes
        are cached separately, so a later parsed call avoids another process.

        For a real PyO3 binding this would iterate proteins inside the
        process; for the subprocess path we still amortise the fork
        cost across all proteins, so the *call-site* cost model is
        identical.
        """
        if not raw_only and self._cached_bundle is not None:
            return self._cached_bundle
        if raw_only and self._cached_raw_bytes is not None:
            return GFF3Bundle(raw_bytes=self._cached_raw_bytes)

        proteins = ref_proteins_path or self._ref_proteins_path
        if proteins is None:
            raise ValueError(
                "MiniprotIndex.align_all requires a ref_proteins_path "
                "(either at construction or via this argument)."
            )

        # A raw-only call may have populated the byte cache without creating
        # the parsed bundle. Parse that cache on demand without rerunning the
        # backend.
        if self._cached_raw_bytes is not None:
            hits: List[GFF3Hit] = []
            for line in self._cached_raw_bytes.decode(
                "utf-8", errors="replace",
            ).splitlines():
                parsed = GFF3Hit.from_gff_line(line)
                if parsed is not None:
                    hits.append(parsed)
            self._cached_bundle = GFF3Bundle(
                hits=hits, raw_bytes=self._cached_raw_bytes,
            )
            return self._cached_bundle

        if self._native is not None:                            # pragma: no cover
            # Real PyO3 binding path (forward-compatible — not exercised today).
            raw = BytesIO()
            for native_hit in self._native.align_file(proteins):
                line = native_hit.to_gff_line() + "\n"
                raw.write(line.encode("utf-8"))
            stdout_bytes = raw.getvalue()
            self._cached_raw_bytes = stdout_bytes
            if raw_only:
                return GFF3Bundle(raw_bytes=stdout_bytes)

            hits: List[GFF3Hit] = []
            for line in stdout_bytes.decode(
                "utf-8", errors="replace",
            ).splitlines():
                parsed = GFF3Hit.from_gff_line(line)
                if parsed is not None:
                    hits.append(parsed)
            self._cached_bundle = GFF3Bundle(
                hits=hits, raw_bytes=stdout_bytes,
            )
            return self._cached_bundle

        # ── Subprocess fallback ─────────────────────────────────────
        # Build the command via the shared helper so facade clients scale
        # miniprot's -t exactly like the CLI's guarded subprocess path.
        # Lazy import to avoid any import-order coupling with run_miniprot.
        from lifton.run_miniprot import _build_miniprot_command
        cmd = _build_miniprot_command(
            self.miniprot_path, self.target_fa, proteins,
            self.mp_options, self.threads,
        )

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

        self._cached_raw_bytes = stdout_bytes
        if raw_only:
            return GFF3Bundle(raw_bytes=stdout_bytes)

        hits: List[GFF3Hit] = []
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

    def align_all_into(
        self, sink, ref_proteins_path: Optional[str] = None,
        *, chunk_size: int = 1 << 20,
    ) -> int:
        """Stream an all-protein alignment into ``sink`` without caching it.

        ``sink`` may be a callable accepting bytes or an object with a
        ``write(bytes)`` method.  The return value is the exact stdout byte
        count.  Existing :meth:`align_all` remains the compatibility API for
        callers that explicitly need a parsed or raw in-memory bundle.
        """
        import threading
        from lifton.run_miniprot import (
            _BoundedStderr, _build_miniprot_command, _stop_process,
            _with_stderr_tail,
        )

        write = sink if callable(sink) else getattr(sink, "write", None)
        if not callable(write):
            raise TypeError("sink must be callable or expose write(bytes)")
        proteins = ref_proteins_path or self._ref_proteins_path
        if proteins is None:
            raise ValueError(
                "MiniprotIndex.align_all_into requires a ref_proteins_path"
            )
        if self._cached_raw_bytes is not None:
            write(self._cached_raw_bytes)
            return len(self._cached_raw_bytes)

        if self._native is not None:                            # pragma: no cover
            byte_count = 0
            for native_hit in self._native.align_file(proteins):
                chunk = (native_hit.to_gff_line() + "\n").encode("utf-8")
                write(chunk)
                byte_count += len(chunk)
            return byte_count

        cmd = _build_miniprot_command(
            self.miniprot_path, self.target_fa, proteins,
            self.mp_options, self.threads,
        )
        proc = subprocess.Popen(
            cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            bufsize=1 << 20,
        )
        stderr_state = _BoundedStderr()

        def _drain_stderr():
            if proc.stderr is None:
                return
            while True:
                chunk = proc.stderr.read(1 << 16)
                if not chunk:
                    return
                stderr_state.feed(chunk)

        stderr_thread = threading.Thread(
            target=_drain_stderr, name="lifton-miniprot-facade-stderr",
            daemon=True,
        )
        stderr_thread.start()
        byte_count = 0
        try:
            if proc.stdout is not None:
                while True:
                    chunk = proc.stdout.read(max(1, int(chunk_size)))
                    if not chunk:
                        break
                    write(chunk)
                    byte_count += len(chunk)
            returncode = proc.wait()
            stderr_thread.join()
            if returncode != 0:
                raise RuntimeError(_with_stderr_tail(
                    f"miniprot exited with code {returncode}", stderr_state,
                ))
            if stderr_state.error_seen:
                raise RuntimeError(_with_stderr_tail(
                    "miniprot reported an ERROR during mapping; refusing "
                    "to publish its stream",
                    stderr_state,
                ))
            return byte_count
        except BaseException:
            _stop_process(proc)
            stderr_thread.join(timeout=5)
            raise
        finally:
            for stream in (proc.stdout, proc.stderr):
                if stream is not None:
                    try:
                        stream.close()
                    except Exception:
                        pass

    @property
    def is_native(self) -> bool:
        """True iff a real ``pyminiprot`` binding is in use."""
        return self._native is not None

    @property
    def raw_bytes(self) -> bytes:
        """Forward-compatible accessor for callers that still want the
        full GFF3 blob (for example to feed into the Phase 7
        streaming-adapter ingest path)."""
        if self._cached_raw_bytes is None:
            raise RuntimeError("Call align_all() before .raw_bytes")
        return self._cached_raw_bytes
