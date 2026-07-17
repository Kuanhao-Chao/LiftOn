"""Transactional publication of LiftOn GFF3 output.

The pipeline can write and validate a complete GFF3 document at
``working_path`` before making it visible at the requested destination.  A
successful transaction publishes a file with one atomic :func:`os.replace`;
an unsuccessful transaction preserves the staged bytes as a
``*.partial.gff3`` diagnostic artifact.

``stdout`` cannot be replaced atomically, so that mode is buffered on disk and
copied to the supplied stream only when :meth:`OutputTransaction.commit` is
called.
"""

from __future__ import annotations

import os
import shutil
import signal
import sys
import tempfile
import threading
from enum import Enum
from pathlib import Path
from typing import IO, TextIO


class OutputTransactionStateError(RuntimeError):
    """Raised when an operation conflicts with a terminal transaction state."""


class OutputTransactionState(str, Enum):
    """Observable lifecycle states for :class:`OutputTransaction`."""

    OPEN = "open"
    CLOSED = "closed"
    COMMITTED = "committed"
    ABORTED = "aborted"
    FAILED = "failed"


def _partial_name(destination: Path) -> str:
    """Return a predictable, human-readable partial-output filename."""

    lower_name = destination.name.lower()
    for suffix in (".gff3", ".gff"):
        if lower_name.endswith(suffix):
            return destination.name[: -len(suffix)] + ".partial.gff3"
    if destination.suffix:
        return destination.stem + ".partial" + destination.suffix
    return destination.name + ".partial"


def _fsync_directory(path: Path) -> None:
    """Persist a same-directory rename where the platform supports it."""

    flags = os.O_RDONLY | getattr(os, "O_DIRECTORY", 0)
    descriptor = os.open(path, flags)
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


class OutputTransaction:
    """Stage GFF3 text and publish it only after successful completion.

    Parameters
    ----------
    output
        Destination path, or the literal string ``"stdout"``.
    stdout
        Optional text stream used in ``stdout`` mode.  Defaults to the current
        :data:`sys.stdout`; supplying a stream is useful for embedding and
        tests.
    work_dir
        Staging directory for ``stdout`` mode.  File destinations always stage
        in their own directory so that commit uses an atomic same-filesystem
        replacement.  Defaults to the current directory.

    Notes
    -----
    ``commit()``, ``abort()``, and ``close()`` are safe to repeat.  Switching
    from one terminal outcome to the other is an error: committing an aborted
    transaction (or aborting a committed one) raises
    :class:`OutputTransactionStateError`.
    """

    def __init__(
        self,
        output: str | os.PathLike[str],
        *,
        stdout: TextIO | None = None,
        work_dir: str | os.PathLike[str] | None = None,
        encoding: str = "utf-8",
    ) -> None:
        output_text = os.fspath(output)
        self._stdout_mode = output_text == "stdout"
        self._destination = None if self._stdout_mode else Path(output_text)
        self._stdout = sys.stdout if stdout is None else stdout
        self._encoding = encoding

        if self._stdout_mode:
            staging_dir = Path.cwd() if work_dir is None else Path(work_dir)
            self._partial_path = staging_dir / "lifton.partial.gff3"
            prefix = ".lifton.stdout."
        else:
            assert self._destination is not None
            staging_dir = self._destination.parent
            self._partial_path = staging_dir / _partial_name(self._destination)
            prefix = f".{self._destination.name}."

        if not staging_dir.is_dir():
            raise FileNotFoundError(
                f"output staging directory does not exist: {staging_dir}"
            )

        stream = tempfile.NamedTemporaryFile(
            mode="w",
            encoding=encoding,
            newline="",
            prefix=prefix,
            suffix=".tmp",
            dir=staging_dir,
            delete=False,
        )
        self._stream: IO[str] = stream
        self._working_path = Path(stream.name)
        self._state = OutputTransactionState.OPEN
        self._signal_handlers: dict[int, object] = {}

    @property
    def state(self) -> OutputTransactionState:
        """Current transaction lifecycle state."""

        return self._state

    @property
    def working_path(self) -> Path:
        """Path containing staged output until commit or abort."""

        return self._working_path

    @property
    def destination_path(self) -> Path | None:
        """Requested file destination, or ``None`` in stdout mode."""

        return self._destination

    @property
    def partial_path(self) -> Path:
        """Path used to preserve staged output on abort."""

        return self._partial_path

    def open(self) -> IO[str]:
        """Return the staging stream while the transaction is open.

        Calling ``open()`` repeatedly returns the same stream.  A closed or
        terminal transaction cannot be reopened.
        """

        if self._state is OutputTransactionState.OPEN:
            return self._stream
        raise OutputTransactionStateError(
            f"cannot open an output transaction in state {self._state.value!r}"
        )

    def close(self) -> None:
        """Flush, fsync, and close staged output without publishing it."""

        if self._state in {
            OutputTransactionState.CLOSED,
            OutputTransactionState.COMMITTED,
            OutputTransactionState.ABORTED,
        }:
            return
        if self._state is OutputTransactionState.FAILED:
            raise OutputTransactionStateError(
                "cannot close a failed output transaction"
            )

        if not self._stream.closed:
            try:
                self._stream.flush()
                os.fsync(self._stream.fileno())
                self._stream.close()
            except BaseException:
                self._state = OutputTransactionState.FAILED
                try:
                    self._stream.close()
                except BaseException:
                    pass
                raise
        self._state = OutputTransactionState.CLOSED

    def commit(self) -> Path | None:
        """Publish staged output and return its path (``None`` for stdout)."""

        if self._state is OutputTransactionState.COMMITTED:
            return self._destination
        if self._state is OutputTransactionState.ABORTED:
            raise OutputTransactionStateError(
                "cannot commit an aborted output transaction"
            )
        if self._state is OutputTransactionState.FAILED:
            raise OutputTransactionStateError(
                "cannot commit a failed output transaction"
            )

        self.close()
        if self._stdout_mode:
            with self._working_path.open(
                "r", encoding=self._encoding, newline=""
            ) as staged:
                shutil.copyfileobj(staged, self._stdout)
            self._stdout.flush()
            self._working_path.unlink()
            self._state = OutputTransactionState.COMMITTED
            try:
                _fsync_directory(self._working_path.parent)
            finally:
                self.restore_signal_handlers()
        else:
            assert self._destination is not None
            os.replace(self._working_path, self._destination)
            self._state = OutputTransactionState.COMMITTED
            try:
                _fsync_directory(self._destination.parent)
            finally:
                self.restore_signal_handlers()

        self._state = OutputTransactionState.COMMITTED
        return self._destination

    def abort(self) -> Path:
        """Preserve staged output under ``partial_path`` and return that path."""

        if self._state is OutputTransactionState.ABORTED:
            return self._partial_path
        if self._state is OutputTransactionState.COMMITTED:
            raise OutputTransactionStateError(
                "cannot abort a committed output transaction"
            )

        if self._state is OutputTransactionState.OPEN:
            try:
                self.close()
            except BaseException:
                # Preserve whatever reached the staging file; callers still
                # receive an error if publication was attempted.
                pass
        elif self._state is OutputTransactionState.FAILED and not self._stream.closed:
            try:
                self._stream.close()
            except BaseException:
                pass
        try:
            os.replace(self._working_path, self._partial_path)
            self._state = OutputTransactionState.ABORTED
            _fsync_directory(self._partial_path.parent)
        finally:
            # Cleanup must not mask the publication/durability error that
            # brought us here. On success, restoration failures still surface.
            if sys.exc_info()[0] is None:
                self.restore_signal_handlers()
            else:
                try:
                    self.restore_signal_handlers()
                except BaseException:
                    pass
        return self._partial_path

    def install_signal_handlers(self) -> bool:
        """Preserve staged output on SIGTERM/SIGHUP in the main thread."""

        if threading.current_thread() is not threading.main_thread():
            return False
        if self._signal_handlers:
            return True
        for signum in (signal.SIGTERM, getattr(signal, "SIGHUP", None)):
            if signum is None:
                continue
            self._signal_handlers[int(signum)] = signal.getsignal(signum)
            signal.signal(signum, self._handle_signal)
        return True

    def restore_signal_handlers(self) -> None:
        if not self._signal_handlers:
            return
        if threading.current_thread() is not threading.main_thread():
            return
        handlers = self._signal_handlers
        self._signal_handlers = {}
        for signum, previous in handlers.items():
            signal.signal(signum, previous)

    def _handle_signal(self, signum, frame) -> None:
        previous = self._signal_handlers.get(int(signum), signal.SIG_DFL)
        try:
            if self._state not in {
                OutputTransactionState.COMMITTED,
                OutputTransactionState.ABORTED,
            }:
                self.abort()
        finally:
            self.restore_signal_handlers()
        if callable(previous):
            previous(signum, frame)
        elif previous != signal.SIG_IGN:
            raise SystemExit(128 + int(signum))

    def __enter__(self) -> IO[str]:
        return self.open()

    def __exit__(self, exc_type, exc_value, traceback) -> bool:
        if exc_type is None:
            self.commit()
        else:
            self.abort()
        return False
