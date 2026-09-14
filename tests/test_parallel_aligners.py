"""Step-4 scheduling byte identity and thread-affinity gates.

Step 4 overlaps the two external aligner programs (Liftoff + miniprot) by
default for targets up to 4 Gb and serializes them above that resource-policy
boundary. ``--serial-aligners`` and ``--parallel-aligners`` are explicit
overrides. Scheduling is byte-neutral.

This is a *dedicated* test rather than a 5th axis on the 24-cell matrix
(`tests/test_native_matrix.py`): the change is orthogonal to the
stream/inmemory/threads/native axes (Step-4 scheduling only).

The `hermetic_pipeline` fixture supplies pre-baked results for byte checks.
Tests of active scheduling omit `-L`/`-M` and replace only LiftOn's two aligner
wrappers with functions that return those same fixtures. This reaches the real
ThreadPoolExecutor and `.result()` plumbing without invoking native tools.
"""

from __future__ import annotations

import json
import sys

from tests.test_integration_pipeline import (  # noqa: F401
    integration_workspace,
    hermetic_pipeline,
)


def _drive(workspace, *, serial: bool = False, alias: bool = False,
           stream: bool = False, inmem: bool = False, threads: int = 1,
           native: bool = False, cached_aligners: bool = True,
           suffix: str) -> bytes:
    from lifton import lifton as lifton_main

    out_gff = workspace["out"] / f"lifton_{suffix}.gff3"
    argv = [
        str(workspace["tgt_fa"]),
        str(workspace["ref_fa"]),
        "-g", str(workspace["ref_gff"]),
        "-o", str(out_gff),
        "-ad", "RefSeq",
        "--force",
    ]
    if cached_aligners:
        argv += [
            "-L", str(workspace["liftoff"]),
            "-M", str(workspace["miniprot"]),
        ]
    if serial:
        argv.append("--serial-aligners")
    if alias:
        argv.append("--parallel-aligners")
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
            self, integration_workspace, hermetic_pipeline, monkeypatch):
        """The core contract: the default (concurrent) Step 4 emits the same
        bytes as `--serial-aligners`."""
        from lifton import lifton_utils

        monkeypatch.setattr(
            lifton_utils, "exec_liftoff",
            lambda *args, **kwargs: str(integration_workspace["liftoff"]),
        )
        monkeypatch.setattr(
            lifton_utils, "exec_miniprot",
            lambda *args, **kwargs: str(integration_workspace["miniprot"]),
        )
        concurrent = _drive(
            integration_workspace, cached_aligners=False, suffix="default",
        )
        concurrent_manifest = json.loads((
            integration_workspace["out"] / "lifton_output" /
            "run_manifest.json"
        ).read_text())
        serial = _drive(
            integration_workspace, serial=True, cached_aligners=False,
            suffix="serial",
        )
        serial_manifest = json.loads((
            integration_workspace["out"] / "lifton_output" /
            "run_manifest.json"
        ).read_text())
        assert len(serial) > 0
        assert concurrent == serial, (
            "default concurrent Step 4 diverged from --serial-aligners — "
            "the overlap must be a pure scheduling change"
        )
        assert concurrent_manifest["aligners"]["schedule"]["mode"] == "parallel"
        assert serial_manifest["aligners"]["schedule"]["mode"] == "serial"
        target_stats = concurrent_manifest["input_statistics"]["target_genome"]
        assert target_stats["sequence_count"] == 1
        assert target_stats["total_bases"] > 0

    def test_parallel_override_is_byte_neutral_with_cached_inputs(
            self, integration_workspace, hermetic_pipeline):
        """A schedule override cannot change annotation selection or bytes."""
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
    def test_large_target_default_runs_both_on_main_thread(
            self, integration_workspace, hermetic_pipeline, monkeypatch):
        """The CLI must apply the automatic large-target policy, not merely
        expose it through the pure schedule resolver."""
        import threading
        from lifton import lifton as lifton_main
        from lifton import lifton_utils

        main_id = threading.get_ident()
        seen = []
        monkeypatch.setattr(
            lifton_main, "target_fasta_statistics",
            lambda _fasta: {
                "sequence_count": 1,
                "total_bases": 4_000_000_001,
                "maximum_sequence_bases": 8,
            },
        )

        def _rec_lift(*_args, **_kwargs):
            seen.append(("liftoff", threading.get_ident()))
            return str(integration_workspace["liftoff"])

        def _rec_mini(*_args, **_kwargs):
            seen.append(("miniprot", threading.get_ident()))
            return str(integration_workspace["miniprot"])

        monkeypatch.setattr(lifton_utils, "exec_liftoff", _rec_lift)
        monkeypatch.setattr(lifton_utils, "exec_miniprot", _rec_mini)

        _drive(
            integration_workspace, cached_aligners=False,
            suffix="large_target_default",
        )

        manifest = json.loads((
            integration_workspace["out"] / "lifton_output" /
            "run_manifest.json"
        ).read_text())
        assert seen == [("liftoff", main_id), ("miniprot", main_id)]
        assert manifest["aligners"]["schedule"]["mode"] == "serial"
        assert manifest["aligners"]["schedule"]["reason"] == (
            "target_exceeds_concurrent_boundary"
        )

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
        def _rec_lift(*a, **k):
            seen["liftoff"] = threading.get_ident()
            return str(integration_workspace["liftoff"])

        def _rec_mini(*a, **k):
            seen["miniprot"] = threading.get_ident()
            return str(integration_workspace["miniprot"])

        monkeypatch.setattr(lifton_utils, "exec_liftoff", _rec_lift)
        monkeypatch.setattr(lifton_utils, "exec_miniprot", _rec_mini)

        _drive(
            integration_workspace, cached_aligners=False,
            suffix="threadaffinity",
        )

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
        def _rec_lift(*a, **k):
            seen["liftoff"] = threading.get_ident()
            return str(integration_workspace["liftoff"])

        def _rec_mini(*a, **k):
            seen["miniprot"] = threading.get_ident()
            return str(integration_workspace["miniprot"])

        monkeypatch.setattr(lifton_utils, "exec_liftoff", _rec_lift)
        monkeypatch.setattr(lifton_utils, "exec_miniprot", _rec_mini)

        _drive(
            integration_workspace, serial=True, cached_aligners=False,
            suffix="serial_affinity",
        )

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
