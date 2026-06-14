"""Iteration 8 — fresh parallel Step 7 WITHOUT --native.

Before Iteration 8, ``--threads N --locus-pipeline`` only fanned Step 7
out across worker threads when ``--native`` (or ``LIFTON_PARALLEL_FORCE``)
was set; on a plain default (gffutils SQLite) run it silently downgraded
to serial. Iteration 8 routes the parallel path through the
materialise + proxy-DB machinery unconditionally, so it works on any
backend without ``--native`` — and stays byte-identical.

These tests pin the new contract directly on the default backend:
  1. ``-t4 --locus-pipeline`` (no --native) is byte-identical to ``-t1``;
  2. a ThreadPoolExecutor is genuinely created (it is no longer serial);
  3. no-native parallel output == ``--native`` parallel output.

Hermeticity: the integration_workspace fixture supplies pre-baked
liftoff/miniprot GFF3s via -L / -M, so minimap2 / miniprot are never
invoked (Step 4 short-circuits to a load). Only the Step-7 dispatch is
exercised.
"""

from __future__ import annotations

import pytest

from tests.test_integration_pipeline import (  # noqa: F401
    integration_workspace,
    hermetic_pipeline,
)


def _drive(workspace, *, threads: int, native: bool, locus_pipeline: bool,
           suffix: str) -> bytes:
    from lifton import lifton as lifton_main

    out_gff = workspace["out"] / f"fresh_p7_{suffix}.gff3"
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
    if threads > 1:
        argv += ["-t", str(threads)]
    if locus_pipeline:
        argv.append("--locus-pipeline")
    if native:
        argv.append("--native")
    args = lifton_main.parse_args(argv)
    lifton_main.run_all_lifton_steps(args)
    return out_gff.read_bytes()


class TestFreshParallelStep7:
    def test_parallel_no_native_byte_identical_to_serial(
            self, integration_workspace, hermetic_pipeline):
        """`-t4 --locus-pipeline` with NO --native == serial `-t1`."""
        serial = _drive(integration_workspace, threads=1, native=False,
                        locus_pipeline=False, suffix="serial")
        parallel = _drive(integration_workspace, threads=4, native=False,
                          locus_pipeline=True, suffix="par_nonnative")
        assert len(serial) > 0
        assert parallel == serial

    def test_parallel_no_native_actually_pools(
            self, integration_workspace, hermetic_pipeline, monkeypatch):
        """The default-backend parallel run genuinely constructs a
        ThreadPoolExecutor — it is no longer a silent serial downgrade."""
        from lifton import parallel as _parallel
        original = _parallel.ThreadPoolExecutor
        constructed = {"n": 0}

        class TPESpy(original):
            def __init__(self, *a, **kw):
                constructed["n"] += 1
                super().__init__(*a, **kw)

        monkeypatch.setattr(_parallel, "ThreadPoolExecutor", TPESpy)
        out = _drive(integration_workspace, threads=4, native=False,
                     locus_pipeline=True, suffix="pool_probe")
        assert len(out) > 0
        assert constructed["n"] >= 1, (
            "parallel Step 7 on the default backend did not create a "
            "ThreadPoolExecutor without --native"
        )

    def test_no_native_parallel_equals_native_parallel(
            self, integration_workspace, hermetic_pipeline):
        """The win is now available without --native: the parallel
        output with and without --native is byte-for-byte identical."""
        no_native = _drive(integration_workspace, threads=4, native=False,
                           locus_pipeline=True, suffix="par_nonnative2")
        with_native = _drive(integration_workspace, threads=4, native=True,
                             locus_pipeline=True, suffix="par_native")
        assert len(no_native) > 0
        assert no_native == with_native

    @pytest.mark.parametrize("threads", [1, 2, 4])
    def test_thread_counts_byte_identical_no_native(
            self, integration_workspace, hermetic_pipeline, threads):
        """Determinism across thread counts on the default backend
        without --native."""
        baseline = _drive(integration_workspace, threads=1, native=False,
                          locus_pipeline=False, suffix="axis_baseline")
        result = _drive(integration_workspace, threads=threads, native=False,
                        locus_pipeline=(threads > 1),
                        suffix=f"axis_t{threads}")
        assert result == baseline
