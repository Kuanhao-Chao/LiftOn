"""Phase 9 — full 12-cell parallelism + I/O matrix.

The Phase 5 zero-bug contract says the output GFF3 must be byte-identical
across every supported flag combination. Phase 9 expands that gate to
include a third axis: ``--threads ∈ {1, 2, 4}``.

Cells: ``--stream={off,on}`` × ``--inmemory-liftoff={off,on}`` ×
``--threads ∈ {1,2,4}`` (with the locus-pipeline flag on for any
threads > 1) = **12 cells**, all expected byte-identical.
"""

from __future__ import annotations

import os
from pathlib import Path

import pytest

from tests.test_integration_pipeline import (  # noqa: F401
    integration_workspace,
    hermetic_pipeline,
)


def _drive(workspace, *, stream: bool, inmem: bool, threads: int,
           suffix: str) -> bytes:
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
    if stream:
        argv.append("--stream")
    if inmem:
        argv.append("--inmemory-liftoff")
    if threads > 1:
        argv += ["-t", str(threads), "--locus-pipeline"]
    args = lifton_main.parse_args(argv)
    lifton_main.run_all_lifton_steps(args)
    return out_gff.read_bytes()


# ---------------------------------------------------------------------------
# Threads-only axis (with stream + inmemory both off)
# ---------------------------------------------------------------------------

class TestThreadsAxis:
    @pytest.mark.parametrize("threads", [1, 2, 4])
    def test_thread_count_yields_valid_output(self, integration_workspace,
                                              hermetic_pipeline, threads):
        b = _drive(integration_workspace,
                   stream=False, inmem=False, threads=threads,
                   suffix=f"t{threads}")
        assert len(b) > 0

    def test_threads_1_2_4_all_byte_identical(self, integration_workspace,
                                              hermetic_pipeline):
        outputs = {
            t: _drive(integration_workspace,
                      stream=False, inmem=False, threads=t,
                      suffix=f"axis_t{t}")
            for t in (1, 2, 4)
        }
        assert outputs[1] == outputs[2] == outputs[4], (
            "Output diverged across thread counts — determinism gate failed"
        )


# ---------------------------------------------------------------------------
# Determinism with stream + inmemory both on
# ---------------------------------------------------------------------------

class TestThreadsAxisStreamingInmemory:
    @pytest.mark.parametrize("threads", [1, 2, 4])
    def test_under_streaming_path(self, integration_workspace,
                                  hermetic_pipeline, threads):
        b = _drive(integration_workspace,
                   stream=True, inmem=True, threads=threads,
                   suffix=f"si_t{threads}")
        assert len(b) > 0

    def test_streaming_threads_byte_identical(self, integration_workspace,
                                              hermetic_pipeline):
        outputs = {
            t: _drive(integration_workspace,
                      stream=True, inmem=True, threads=t,
                      suffix=f"si_axis_t{t}")
            for t in (1, 2, 4)
        }
        assert outputs[1] == outputs[2] == outputs[4]


# ---------------------------------------------------------------------------
# Full 12-cell matrix
# ---------------------------------------------------------------------------

class TestFull12CellMatrix:
    def test_all_twelve_combinations_byte_identical(
            self, integration_workspace, hermetic_pipeline):
        """The Phase 9 golden output gate: every cell of the
        2 × 2 × 3 matrix produces identical bytes."""
        outputs = {}
        for stream in (False, True):
            for inmem in (False, True):
                for threads in (1, 2, 4):
                    key = f"s{int(stream)}_i{int(inmem)}_t{threads}"
                    outputs[key] = _drive(
                        integration_workspace,
                        stream=stream, inmem=inmem, threads=threads,
                        suffix=key,
                    )
        baseline = outputs["s0_i0_t1"]
        diverged = [k for k, v in outputs.items() if v != baseline]
        assert not diverged, (
            f"Output diverged for: {diverged}. "
            f"Lengths: {[(k, len(outputs[k])) for k in diverged]}"
        )
        # Sanity check: every output non-trivial
        assert all(len(v) > 0 for v in outputs.values())


# ---------------------------------------------------------------------------
# Locus-pipeline path requires threads>1 to fan out
# ---------------------------------------------------------------------------

class TestLocusPipelineSemantics:
    def test_locus_pipeline_with_threads_1_uses_serial(
            self, integration_workspace, hermetic_pipeline, monkeypatch):
        """When --locus-pipeline is set but --threads is 1 (the
        default), the dispatcher should still take the serial path
        (no ThreadPoolExecutor created)."""
        from lifton import parallel
        original = parallel.ThreadPoolExecutor
        constructed = {"n": 0}

        class TPESpy(original):
            def __init__(self, *a, **kw):
                constructed["n"] += 1
                super().__init__(*a, **kw)

        monkeypatch.setattr(parallel, "ThreadPoolExecutor", TPESpy)

        from lifton import lifton as lifton_main
        out_gff = integration_workspace["out"] / "lifton_lp_t1.gff3"
        argv = [
            str(integration_workspace["tgt_fa"]),
            str(integration_workspace["ref_fa"]),
            "-g", str(integration_workspace["ref_gff"]),
            "-L", str(integration_workspace["liftoff"]),
            "-M", str(integration_workspace["miniprot"]),
            "-o", str(out_gff),
            "-ad", "RefSeq", "--force",
            "--locus-pipeline",
        ]
        args = lifton_main.parse_args(argv)
        lifton_main.run_all_lifton_steps(args)
        assert out_gff.exists()
        assert constructed["n"] == 0

    def test_locus_pipeline_with_threads_4_creates_pool(
            self, integration_workspace, hermetic_pipeline, monkeypatch):
        """With LIFTON_PARALLEL_FORCE + a gffbase backend the dispatcher
        creates at least one ThreadPoolExecutor.

        ``>= 1`` (not ``== 1``) keeps this robust to the dispatch shape.
        Iteration 10 fused the materialise + process phases into ONE pool
        on the on-disk path (so a viable run now builds a single fused pool,
        not a prefetcher pool + a worker pool); the exact-count contracts are
        pinned in ``tests/test_fuse_step7.py`` (fused → 1, ``LIFTON_FUSE_STEP7=0``
        two-phase → 2).
        """
        monkeypatch.setenv("LIFTON_USE_GFFBASE", "1")
        monkeypatch.setenv("LIFTON_PARALLEL_FORCE", "1")

        from lifton import parallel
        original = parallel.ThreadPoolExecutor
        constructed = {"n": 0}

        class TPESpy(original):
            def __init__(self, *a, **kw):
                constructed["n"] += 1
                super().__init__(*a, **kw)

        monkeypatch.setattr(parallel, "ThreadPoolExecutor", TPESpy)

        from lifton import lifton as lifton_main
        out_gff = integration_workspace["out"] / "lifton_lp_t4.gff3"
        argv = [
            str(integration_workspace["tgt_fa"]),
            str(integration_workspace["ref_fa"]),
            "-g", str(integration_workspace["ref_gff"]),
            "-L", str(integration_workspace["liftoff"]),
            "-M", str(integration_workspace["miniprot"]),
            "-o", str(out_gff),
            "-ad", "RefSeq", "--force",
            "--locus-pipeline", "-t", "4",
        ]
        args = lifton_main.parse_args(argv)
        lifton_main.run_all_lifton_steps(args)
        assert out_gff.exists()
        assert constructed["n"] >= 1

    def test_locus_pipeline_parallel_on_sqlite_by_default(
            self, integration_workspace, hermetic_pipeline, monkeypatch):
        """Iteration 8: --locus-pipeline + threads>1 on the DEFAULT
        gffutils (SQLite) backend now runs PARALLEL (via the
        materialised-payload + proxy-DB path) WITHOUT --native — it no
        longer silently downgrades to serial. A ThreadPoolExecutor is
        created and the output stays byte-identical to a serial run.
        """
        from lifton import parallel
        original = parallel.ThreadPoolExecutor
        constructed = {"n": 0}

        class TPESpy(original):
            def __init__(self, *a, **kw):
                constructed["n"] += 1
                super().__init__(*a, **kw)

        monkeypatch.setattr(parallel, "ThreadPoolExecutor", TPESpy)

        from lifton import lifton as lifton_main
        out_serial = integration_workspace["out"] / "default_serial.gff3"
        out_parallel = integration_workspace["out"] / "default_parallel.gff3"

        base = [
            str(integration_workspace["tgt_fa"]),
            str(integration_workspace["ref_fa"]),
            "-g", str(integration_workspace["ref_gff"]),
            "-L", str(integration_workspace["liftoff"]),
            "-M", str(integration_workspace["miniprot"]),
            "-ad", "RefSeq", "--force",
        ]
        # Serial reference (no pool expected)
        lifton_main.run_all_lifton_steps(
            lifton_main.parse_args(base + ["-o", str(out_serial)]))
        assert constructed["n"] == 0, "serial run unexpectedly created a pool"

        # Parallel on the default gffutils backend — a pool IS created.
        lifton_main.run_all_lifton_steps(lifton_main.parse_args(
            base + ["-o", str(out_parallel), "--locus-pipeline", "-t", "4"]))
        assert constructed["n"] >= 1, (
            "parallel Step 7 did not create a ThreadPoolExecutor on the "
            "default gffutils backend (silently fell back to serial?)"
        )
        # Byte-identical to serial — the whole point of the contract.
        assert out_serial.read_bytes() == out_parallel.read_bytes()

    def test_block_gffutils_restores_serial_fallback(
            self, integration_workspace, hermetic_pipeline, monkeypatch, capsys):
        """The LIFTON_PARALLEL_BLOCK_GFFUTILS opt-out restores the
        pre-Iteration-8 strict serial fallback on gffutils: no pool is
        created, a warning is emitted, and output stays byte-identical.
        """
        monkeypatch.setenv("LIFTON_PARALLEL_BLOCK_GFFUTILS", "1")

        from lifton import parallel
        original = parallel.ThreadPoolExecutor
        constructed = {"n": 0}

        class TPESpy(original):
            def __init__(self, *a, **kw):
                constructed["n"] += 1
                super().__init__(*a, **kw)

        monkeypatch.setattr(parallel, "ThreadPoolExecutor", TPESpy)

        from lifton import lifton as lifton_main
        out_serial = integration_workspace["out"] / "block_serial.gff3"
        out_attempted = integration_workspace["out"] / "block_attempted.gff3"

        base = [
            str(integration_workspace["tgt_fa"]),
            str(integration_workspace["ref_fa"]),
            "-g", str(integration_workspace["ref_gff"]),
            "-L", str(integration_workspace["liftoff"]),
            "-M", str(integration_workspace["miniprot"]),
            "-ad", "RefSeq", "--force",
        ]
        lifton_main.run_all_lifton_steps(
            lifton_main.parse_args(base + ["-o", str(out_serial)]))

        lifton_main.run_all_lifton_steps(lifton_main.parse_args(
            base + ["-o", str(out_attempted), "--locus-pipeline", "-t", "4"]))

        # Blocked → fell back to serial → no pool created, byte-identical.
        assert constructed["n"] == 0
        assert out_serial.read_bytes() == out_attempted.read_bytes()
        assert "LIFTON_PARALLEL_BLOCK_GFFUTILS" in capsys.readouterr().err
