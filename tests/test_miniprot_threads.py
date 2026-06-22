"""Iteration 17 — miniprot thread plumbing (un-cap from the binary default 4).

Before Iteration 17 miniprot was always launched with NO ``-t`` flag, so it ran
at its hard-coded binary default of 4 threads regardless of LiftOn's
``-t/--threads N`` — becoming the serial tail of the concurrent Step 4 on
cross-species data (where miniprot does heavy protein alignment, while minimap2
already scaled with ``-t``). Iteration 17 plumbs ``args.threads`` into miniprot's
own ``-t`` via :func:`lifton.run_miniprot._build_miniprot_command`.

The feasibility gate proved miniprot's GFF output is **byte-identical** across
``-t 1/4/8/16`` (input-order-preserving, Heng Li ``kt_pipeline`` lineage), so
this is a byte-neutral change: the default ``-t 1`` run emits NO ``-t`` (gate on
``> 1``) and is byte-identical to the pre-Iteration-17 command, and higher-``-t``
runs change only miniprot's thread count, not its output.

These tests pin the command-construction seam directly (no real miniprot, no
subprocess) — the integration/24-cell suites use cached ``-L``/``-M`` so
``exec_miniprot`` short-circuits and the builder never runs there.
"""

from __future__ import annotations

import subprocess
from types import SimpleNamespace
from unittest import mock

from lifton.run_miniprot import (
    _build_miniprot_command,
    _resolve_miniprot_threads,
)
from lifton.native_bindings.miniprot_facade import MiniprotIndex

BASE = ["miniprot", "--gff-only", "tgt.fa", "ref_proteins.fa"]


def _build(mp_options="", threads=1):
    return _build_miniprot_command("miniprot", "tgt.fa", "ref_proteins.fa",
                                   mp_options, threads)


# ---------------------------------------------------------------------------
# The >1 gate — the default path stays byte-identical to pre-Iteration-17.
# ---------------------------------------------------------------------------

class TestThreadGate:
    def test_threads_1_emits_no_t(self):
        """The default `-t 1` invocation must add NO `-t` (miniprot keeps its
        own default of 4 → command byte-identical to the old code)."""
        assert _build(threads=1) == BASE
        assert "-t" not in _build(threads=1)

    def test_threads_gt1_appends_t(self):
        assert _build(threads=8) == BASE + ["-t", "8"]
        assert _build(threads=16) == BASE + ["-t", "16"]

    def test_threads_2_appends_t2(self):
        # >1 gate fires at the smallest meaningful budget.
        assert _build(threads=2) == BASE + ["-t", "2"]


# ---------------------------------------------------------------------------
# User-supplied -t in mp_options is respected (not duplicated).
# ---------------------------------------------------------------------------

class TestUserOverride:
    def test_explicit_t_in_mp_options_preserved(self):
        cmd = _build(mp_options="-t 2", threads=8)
        # user's -t 2 kept; our -t 8 NOT appended
        assert cmd == BASE + ["-t", "2"]
        assert cmd.count("-t") == 1

    def test_glued_t_in_mp_options_preserved(self):
        # `-t8` (glued, getopt-style) is a user override too — must not dup.
        cmd = _build(mp_options="-t8", threads=8)
        assert cmd == BASE + ["-t8"]
        assert not any(tok == "-t" for tok in cmd)

    def test_trans_does_not_false_positive(self):
        """The bug-catcher: `'-t' in '--trans'` is True as a *substring*, so a
        naive substring check would wrongly suppress the thread flag. The
        exact/glued-token check (`opt.startswith('-t')`) must NOT match the
        double-dashed `--trans`, so `-t 8` is still injected."""
        cmd = _build(mp_options="--trans", threads=8)
        assert "--trans" in cmd
        assert cmd[-2:] == ["-t", "8"]


# ---------------------------------------------------------------------------
# LIFTON_MINIPROT_THREADS escape hatch.
# ---------------------------------------------------------------------------

class TestEnvHatch:
    def test_env_zero_disables(self, monkeypatch):
        """`=0` reproduces the pre-Iteration-17 fixed default (no -t) even at
        high -t — the clean A/B baseline / opt-out."""
        monkeypatch.setenv("LIFTON_MINIPROT_THREADS", "0")
        assert _build(threads=8) == BASE
        assert _build(threads=16) == BASE

    def test_env_forces_count(self, monkeypatch):
        monkeypatch.setenv("LIFTON_MINIPROT_THREADS", "3")
        assert _build(threads=8) == BASE + ["-t", "3"]
        # forced count is independent of the threads arg
        assert _build(threads=1) == BASE + ["-t", "3"]

    def test_env_one_bypasses_gate(self, monkeypatch):
        """An explicit `=1` is honoured (forces `-t 1`), bypassing the >1 gate
        that only applies to the implicit args.threads-derived path."""
        monkeypatch.setenv("LIFTON_MINIPROT_THREADS", "1")
        assert _build(threads=8) == BASE + ["-t", "1"]

    def test_env_garbage_falls_back(self, monkeypatch):
        """A non-integer env value is ignored (can't crash a run)."""
        monkeypatch.setenv("LIFTON_MINIPROT_THREADS", "lots")
        assert _build(threads=8) == BASE + ["-t", "8"]

    def test_env_empty_falls_back(self, monkeypatch):
        monkeypatch.setenv("LIFTON_MINIPROT_THREADS", "")
        assert _build(threads=8) == BASE + ["-t", "8"]

    def test_env_respects_user_override(self, monkeypatch):
        # Even a forced env count must not duplicate an explicit user -t.
        monkeypatch.setenv("LIFTON_MINIPROT_THREADS", "9")
        assert _build(mp_options="-t 2", threads=8) == BASE + ["-t", "2"]


# ---------------------------------------------------------------------------
# _resolve_miniprot_threads unit semantics.
# ---------------------------------------------------------------------------

class TestResolve:
    def test_gate_on_gt1(self, monkeypatch):
        monkeypatch.delenv("LIFTON_MINIPROT_THREADS", raising=False)
        assert _resolve_miniprot_threads(1) is None
        assert _resolve_miniprot_threads(2) == 2
        assert _resolve_miniprot_threads(8) == 8

    def test_non_int_threads_arg(self, monkeypatch):
        monkeypatch.delenv("LIFTON_MINIPROT_THREADS", raising=False)
        assert _resolve_miniprot_threads(None) is None
        assert _resolve_miniprot_threads("oops") is None

    def test_env_overrides(self, monkeypatch):
        monkeypatch.setenv("LIFTON_MINIPROT_THREADS", "0")
        assert _resolve_miniprot_threads(8) is None
        monkeypatch.setenv("LIFTON_MINIPROT_THREADS", "5")
        assert _resolve_miniprot_threads(8) == 5
        assert _resolve_miniprot_threads(1) == 5


# ---------------------------------------------------------------------------
# Facade parity: the native+streaming path threads the SAME command.
# ---------------------------------------------------------------------------

class TestFacadeParity:
    @staticmethod
    def _mock_proc(blob=b"##gff-version 3\n"):
        p = mock.MagicMock()
        p.communicate.return_value = (blob, b"")
        p.returncode = 0
        return p

    def test_facade_threads_stored(self):
        idx = MiniprotIndex("tgt.fa", ref_proteins_path="ref_proteins.fa",
                            threads=8)
        assert idx.threads == 8

    def test_facade_default_threads_no_t(self):
        # Default (threads=1) → command carries no -t (parity with the
        # streaming subprocess path at -t1, which the parity test pins).
        idx = MiniprotIndex("tgt.fa", ref_proteins_path="ref_proteins.fa")
        with mock.patch.object(subprocess, "Popen",
                               return_value=self._mock_proc()) as mp:
            idx.align_all()
        cmd = mp.call_args.args[0]
        assert "-t" not in cmd

    def test_facade_threads_in_command(self):
        idx = MiniprotIndex("tgt.fa", ref_proteins_path="ref_proteins.fa",
                            threads=8)
        with mock.patch.object(subprocess, "Popen",
                               return_value=self._mock_proc()) as mp:
            idx.align_all()
        cmd = mp.call_args.args[0]
        assert cmd[-2:] == ["-t", "8"]
        # the command shape matches the subprocess builder
        assert cmd == ["miniprot", "--gff-only", "tgt.fa", "ref_proteins.fa",
                       "-t", "8"]


# ---------------------------------------------------------------------------
# Call-site robustness: a SimpleNamespace lacking `.threads` must not crash
# (mirrors tests/test_native_bindings.py:345's args object).
# ---------------------------------------------------------------------------

class TestCallSiteRobustness:
    def test_run_miniprot_missing_threads_attr(self):
        from lifton import run_miniprot
        from io import BytesIO

        def _make_proc():
            p = mock.MagicMock()
            p.wait.return_value = 0
            p.returncode = 0
            p.stdout = BytesIO(b"##gff-version 3\n")
            p.stderr = BytesIO(b"")
            return p

        # No `threads` attribute → getattr(args, "threads", 1) → 1 → no -t,
        # and the run must not raise AttributeError.
        args = SimpleNamespace(mp_options="", stream=True, miniprot=None)
        with mock.patch.object(run_miniprot.subprocess, "Popen",
                               return_value=_make_proc()):
            out = run_miniprot.run_miniprot("/tmp/", args, "tgt.fa", "rp.fa")
        assert out == b"##gff-version 3\n"
