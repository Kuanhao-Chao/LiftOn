"""Iteration 6 — concurrent Step 4 byte-identity gate (PROMOTED to default).

Step 4 overlaps the two external aligner programs (Liftoff + miniprot) by
default: miniprot (an independent subprocess) runs on a background worker
thread while Liftoff runs on the main thread. `--serial-aligners` restores
the old sequential path; `--parallel-aligners` is a kept no-op alias. The
overlap is a pure scheduling change, so the default (concurrent) output must
be byte-identical to `--serial-aligners`.

This is a *dedicated* test rather than a 5th axis on the 24-cell matrix
(`tests/test_native_matrix.py`): the change is orthogonal to the
stream/inmemory/threads/native axes (Step-4 scheduling only).

The `hermetic_pipeline` fixture supplies pre-baked `-L`/`-M`, so
`exec_liftoff`/`exec_miniprot` hit their `os.path.exists` short-circuit and
just *load* the GFFs — exercising the ThreadPoolExecutor wrapper,
`.result()` plumbing, and thread placement without invoking real
minimap2/miniprot (both are monkey-patched to raise if reached).
"""

from __future__ import annotations

import sys

from tests.test_integration_pipeline import (  # noqa: F401
    integration_workspace,
    hermetic_pipeline,
)


def _drive(workspace, *, serial: bool = False, alias: bool = False,
           stream: bool = False, inmem: bool = False, threads: int = 1,
           native: bool = False, suffix: str) -> bytes:
    from lifton import lifton as lifton_main

    out_gff = workspace["out"] / f"lifton_{suffix}.gff3"
    argv = [
        str(workspace["tgt_fa"]),
        str(workspace["ref_fa"]),
        "-g", str(workspace["ref_gff"]),
        "-L", str(workspace["liftoff"]),
        "-M", str(workspace["miniprot"]),
        "-o", str(out_gff),
        "-ad", "RefSeq",
        "--force",
    ]
    if serial:
        argv.append("--serial-aligners")
    if alias:
        argv.append("--parallel-aligners")   # kept no-op alias
    if stream:
        argv.append("--stream")
    if inmem:
        argv.append("--inmemory-liftoff")
    if native:
        argv.append("--native")
    if threads > 1:
        argv += ["-t", str(threads), "--locus-pipeline"]
    args = lifton_main.parse_args(argv)
    lifton_main.run_all_lifton_steps(args)
    return out_gff.read_bytes()


class TestConcurrentStep4ByteIdentical:
    def test_default_concurrent_vs_serial_byte_identical(
            self, integration_workspace, hermetic_pipeline):
        """The core contract: the default (concurrent) Step 4 emits the same
        bytes as `--serial-aligners`."""
        concurrent = _drive(integration_workspace, suffix="default")
        serial = _drive(integration_workspace, serial=True, suffix="serial")
        assert len(serial) > 0
        assert concurrent == serial, (
            "default concurrent Step 4 diverged from --serial-aligners — "
            "the overlap must be a pure scheduling change"
        )

    def test_parallel_aligners_alias_is_noop(self, integration_workspace,
                                             hermetic_pipeline):
        """`--parallel-aligners` is a kept no-op alias (concurrent is now the
        default), so it is byte-identical to plain default."""
        default = _drive(integration_workspace, suffix="alias_default")
        alias = _drive(integration_workspace, alias=True, suffix="alias_on")
        assert default == alias

    def test_concurrent_with_stream_and_inmemory(self, integration_workspace,
                                                 hermetic_pipeline):
        """Concurrency must survive the in-memory return shapes: `--stream`
        returns a miniprot bytes blob and `--inmemory-liftoff` a Liftoff
        bytes blob — produced/loaded across the threaded path and assigned
        on the parent after `.result()`."""
        serial = _drive(integration_workspace, serial=True, suffix="si_serial")
        concurrent = _drive(integration_workspace, stream=True, inmem=True,
                            suffix="si_default")
        assert serial == concurrent

    def test_concurrent_with_native_and_locus_pipeline(
            self, integration_workspace, hermetic_pipeline):
        """Concurrency must not interact with the Step-7 thread pool: the
        default Step 4 + `--native` + `-t 4 --locus-pipeline` is
        byte-identical to `--serial-aligners`."""
        serial = _drive(integration_workspace, serial=True, suffix="nl_serial")
        concurrent = _drive(integration_workspace, native=True, threads=4,
                            suffix="nl_default")
        assert serial == concurrent


class TestConcurrentStep4ThreadAffinity:
    def test_default_runs_liftoff_on_main_miniprot_on_worker(
            self, integration_workspace, hermetic_pipeline, monkeypatch):
        """By default Liftoff must run on the main thread (its SQLite ref_db
        is thread-bound) while miniprot runs on a background worker. This
        pins the fix for the `sqlite3.ProgrammingError: SQLite objects
        created in a thread can only be used in that same thread` crash the
        fresh-Step-4 A/B surfaced (cached-`-L` byte tests can't reach it —
        `exec_liftoff` short-circuits before touching ref_db)."""
        import threading
        from lifton import lifton_utils

        main_id = threading.get_ident()
        seen = {}
        orig_lift = lifton_utils.exec_liftoff
        orig_mini = lifton_utils.exec_miniprot

        def _rec_lift(*a, **k):
            seen["liftoff"] = threading.get_ident()
            return orig_lift(*a, **k)

        def _rec_mini(*a, **k):
            seen["miniprot"] = threading.get_ident()
            return orig_mini(*a, **k)

        monkeypatch.setattr(lifton_utils, "exec_liftoff", _rec_lift)
        monkeypatch.setattr(lifton_utils, "exec_miniprot", _rec_mini)

        _drive(integration_workspace, suffix="threadaffinity")

        assert seen.get("liftoff") == main_id, (
            "Liftoff must run on the main thread (SQLite ref_db is thread-bound)"
        )
        assert seen.get("miniprot") is not None and seen["miniprot"] != main_id, (
            "miniprot must run on a background worker thread"
        )

    def test_serial_runs_both_on_main_thread(
            self, integration_workspace, hermetic_pipeline, monkeypatch):
        """`--serial-aligners` runs both aligners inline on the main thread."""
        import threading
        from lifton import lifton_utils

        main_id = threading.get_ident()
        seen = {}
        orig_lift = lifton_utils.exec_liftoff
        orig_mini = lifton_utils.exec_miniprot

        def _rec_lift(*a, **k):
            seen["liftoff"] = threading.get_ident()
            return orig_lift(*a, **k)

        def _rec_mini(*a, **k):
            seen["miniprot"] = threading.get_ident()
            return orig_mini(*a, **k)

        monkeypatch.setattr(lifton_utils, "exec_liftoff", _rec_lift)
        monkeypatch.setattr(lifton_utils, "exec_miniprot", _rec_mini)

        _drive(integration_workspace, serial=True, suffix="serial_affinity")

        assert seen.get("liftoff") == main_id
        assert seen.get("miniprot") == main_id


class TestConcurrentStep4RecursionLimit:
    def test_recursion_limit_not_leaked(self, integration_workspace,
                                        hermetic_pipeline):
        """The default concurrent Step-4 path must not leak an elevated
        recursion limit: Liftoff runs on the main thread and does its own
        setrecursionlimit save/restore (run_liftoff.py), so the limit is the
        same before and after a run."""
        original = sys.getrecursionlimit()
        try:
            _drive(integration_workspace, suffix="reclimit")
        finally:
            restored = sys.getrecursionlimit()
        assert restored == original, (
            f"recursion limit leaked: {original} -> {restored}"
        )
