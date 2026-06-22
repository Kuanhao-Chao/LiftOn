"""Phase 16 Tier 2 — regression test for run_liftoff.run_liftoff's
exception handler.

Before the Tier 2 fix, the broad `except Exception as e:` block at
`lifton/run_liftoff.py` interpolated only the exception's str() form
into the error log, which discarded the traceback. For a RecursionError
(which is what every benchmark dataset hit in the most recent Phase 16
run), the str() form is the unhelpful message
"maximum recursion depth exceeded while calling a Python object" — no
hint at all as to which file is recursing. The fix logs
``traceback.format_exc()`` so the deepest LiftOn / vendored-Liftoff
frame is visible.
"""
from __future__ import annotations

import sys
from types import SimpleNamespace
from unittest.mock import patch

import pytest

from lifton import run_liftoff as run_liftoff_module


def _make_recursion(*_args, **_kwargs):
    """Trigger a real RecursionError (matches the failure mode the
    Phase 16 benchmark hits in vendored Liftoff)."""
    def _deep(n):
        return _deep(n + 1)
    return _deep(0)


def test_recursion_error_logs_full_traceback(tmp_path, capsys):
    args = SimpleNamespace(
        polish=False,
        inmemory_liftoff=False,
        output=str(tmp_path / "out.gff3"),
    )
    with patch(
        "lifton.run_liftoff.liftoff_main.run_all_liftoff_steps",
        side_effect=_make_recursion,
    ):
        with pytest.raises(SystemExit) as excinfo:
            run_liftoff_module.run_liftoff(
                str(tmp_path) + "/", None, args,
            )
    assert excinfo.value.code == 1

    err = capsys.readouterr().err
    # The original two error messages must still be present (we only
    # wedged the traceback in between, did not replace the prologue or
    # epilogue).
    assert "Liftoff encountered a fatal error during native execution" in err
    assert "LiftOn cannot proceed without a valid Liftoff baseline" in err
    # The new behaviour: the full traceback is dumped.
    assert "Full Python traceback" in err
    assert "Traceback (most recent call last)" in err
    # The recursive frame must appear by name so a future debugger can
    # identify the recursion site from stderr alone.
    assert "_deep" in err
    # The exception's str() form must also be there so summary tools
    # that grep for it (like the benchmark harness) keep working.
    assert "maximum recursion depth exceeded" in err


def test_non_recursion_exception_also_logs_traceback(tmp_path, capsys):
    """Sanity check: the new traceback dump is not RecursionError-
    specific. Any uncaught exception inside the vendored Liftoff call
    surfaces with full context."""
    args = SimpleNamespace(
        polish=False,
        inmemory_liftoff=False,
        output=str(tmp_path / "out.gff3"),
    )

    def _raise_value_error(*_a, **_kw):
        raise ValueError("synthetic Liftoff failure for the test")

    with patch(
        "lifton.run_liftoff.liftoff_main.run_all_liftoff_steps",
        side_effect=_raise_value_error,
    ):
        with pytest.raises(SystemExit):
            run_liftoff_module.run_liftoff(
                str(tmp_path) + "/", None, args,
            )

    err = capsys.readouterr().err
    assert "synthetic Liftoff failure for the test" in err
    assert "Traceback (most recent call last)" in err
    assert "ValueError" in err


# ---------------------------------------------------------------------------
# Phase 16 Tier 3 — recursion-limit guard around vendored Liftoff
# ---------------------------------------------------------------------------

def test_recursion_limit_is_raised_during_vendored_call(tmp_path):
    """The wrapper must elevate sys.getrecursionlimit() to at least
    10000 before invoking liftoff_main.run_all_liftoff_steps so that
    deep-but-acyclic feature-hierarchy traversal does not trip on
    Python's default 1000-frame limit."""
    args = SimpleNamespace(
        polish=False,
        inmemory_liftoff=False,
        output=str(tmp_path / "out.gff3"),
    )
    captured: dict = {}

    def _capture_limit(*_a, **_kw):
        captured["limit"] = sys.getrecursionlimit()

    pre = sys.getrecursionlimit()
    with patch(
        "lifton.run_liftoff.liftoff_main.run_all_liftoff_steps",
        side_effect=_capture_limit,
    ):
        run_liftoff_module.run_liftoff(str(tmp_path) + "/", None, args)
    post = sys.getrecursionlimit()

    assert captured["limit"] >= 10000, (
        "recursion limit during vendored Liftoff call was "
        f"{captured['limit']}, expected ≥ 10000"
    )
    assert post == pre, (
        f"recursion limit not restored: pre={pre}, post={post}"
    )


# ---------------------------------------------------------------------------
# Phase 16 follow-up — bytes-blob log summary
# ---------------------------------------------------------------------------

def test_describe_annotation_source_path_passthrough():
    """Path-shaped values must pass through unchanged so the existing
    "Creating X annotation database : <path>" log lines keep working."""
    from lifton.lifton import _describe_annotation_source
    assert _describe_annotation_source("/tmp/liftoff.gff3") == "/tmp/liftoff.gff3"
    assert _describe_annotation_source(None) is None


def test_describe_annotation_source_bytes_summary():
    """Bytes blobs (from --inmemory-liftoff / --stream) must NOT be
    interpolated verbatim; they get summarised. Without this, a
    33 MB GFF3 blob renders as one 33 MB stderr line and drowns the
    log file."""
    from lifton.lifton import _describe_annotation_source
    blob = b"##gff-version 3\nchr1\tlifton\tgene\t1\t100\t.\t+\t.\tID=g1\n" * 1000
    rendered = _describe_annotation_source(blob)
    # Output must be short (so it can be logged without flooding).
    assert len(rendered) < 100, rendered
    # The byte length must appear so an operator can still tell the
    # blob is roughly the right size.
    assert str(len(blob)) in rendered.replace(",", "")
    # The blob's actual content must NOT leak into the rendered form.
    assert b"##gff-version" not in rendered.encode()


def test_describe_annotation_source_bytearray_summary():
    """bytearray (a mutable cousin of bytes) must hit the same path."""
    from lifton.lifton import _describe_annotation_source
    blob = bytearray(b"X" * 5_000_000)
    rendered = _describe_annotation_source(blob)
    assert "in-memory bytes" in rendered
    assert "5,000,000" in rendered


def test_recursion_limit_restored_on_exception(tmp_path):
    """The finally clause must restore the original recursion limit
    even when the vendored call raises (otherwise long-running
    processes leak the elevated limit forever)."""
    args = SimpleNamespace(
        polish=False,
        inmemory_liftoff=False,
        output=str(tmp_path / "out.gff3"),
    )
    pre = sys.getrecursionlimit()
    with patch(
        "lifton.run_liftoff.liftoff_main.run_all_liftoff_steps",
        side_effect=ValueError("synthetic"),
    ):
        with pytest.raises(SystemExit):
            run_liftoff_module.run_liftoff(str(tmp_path) + "/", None, args)
    post = sys.getrecursionlimit()
    assert post == pre, (
        f"recursion limit not restored after exception: "
        f"pre={pre}, post={post}"
    )
