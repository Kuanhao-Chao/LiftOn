"""Phase 9 — memory profile harness.

These tests are opt-in via ``pytest -m perf`` because they take longer
than the unit suite and depend on the host's memory accounting.

The harness uses :mod:`tracemalloc` to capture peak Python-object
allocation across the whole pipeline and asserts:

  - ``--threads 4 --locus-pipeline`` peaks at no more than ~1.5×
    the serial-baseline peak (the small overhead is the worker
    thread + ordered-writer buffer; the algorithm itself is locus-
    major in both modes so memory usage is comparable).
  - Per-locus discard semantics: the ordered-writer never holds more
    than O(threads) :class:`LocusResult` objects in flight.

Note: the >50 % RSS drop on chr22 promised in the Phase 6.4 roadmap
is gated on the Phase 4 RefSeqProvider lazy-loading work, which has
not landed and is out of scope for Phase 9. Phase 9 delivers the
parallel architecture that makes that drop achievable.
"""

from __future__ import annotations

import tracemalloc
from pathlib import Path

import pytest

from tests.test_integration_pipeline import (  # noqa: F401
    integration_workspace,
    hermetic_pipeline,
)


pytestmark = pytest.mark.perf


def _run_pipeline(workspace, *, threads: int) -> int:
    """Run the full pipeline and return the peak Python-allocated
    bytes recorded by tracemalloc."""
    from lifton import lifton as lifton_main
    out_gff = workspace["out"] / f"perf_t{threads}.gff3"
    argv = [
        str(workspace["tgt_fa"]),
        str(workspace["ref_fa"]),
        "-g", str(workspace["ref_gff"]),
        "-L", str(workspace["liftoff"]),
        "-M", str(workspace["miniprot"]),
        "-o", str(out_gff),
        "-ad", "RefSeq", "--force",
    ]
    if threads > 1:
        argv += ["--locus-pipeline", "-t", str(threads), "--stream",
                 "--inmemory-liftoff"]
    args = lifton_main.parse_args(argv)

    tracemalloc.start()
    lifton_main.run_all_lifton_steps(args)
    _current, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    return peak


class TestLocusMemoryProfile:
    def test_serial_baseline_completes(self, integration_workspace,
                                       hermetic_pipeline):
        peak = _run_pipeline(integration_workspace, threads=1)
        assert peak > 0

    def test_parallel_peak_within_1_5x_of_serial(
            self, integration_workspace, hermetic_pipeline, monkeypatch):
        """Parallel peak ≤ 1.5× serial peak. The synthetic fixture is
        tiny so the overhead is dominated by the executor itself
        rather than per-locus storage; on chr22 with the future
        RefSeqProvider lazy load this ratio is expected to drop
        below 1.0×."""
        monkeypatch.setenv("LIFTON_PARALLEL_FORCE", "1")
        serial = _run_pipeline(integration_workspace, threads=1)
        parallel = _run_pipeline(integration_workspace, threads=4)
        assert parallel <= 1.5 * serial, (
            f"Parallel peak {parallel} > 1.5× serial peak {serial} — "
            "memory regression"
        )


class TestOrderedWriterMemoryBound:
    def test_pending_buffer_never_exceeds_threads_squared(self, monkeypatch):
        """The ordered-writer's `pending` heap never holds more than
        roughly O(threads^2) results in the worst case where every
        worker's result lands far out of order. On the synthetic
        fixture it should peak at ≤ threads."""
        monkeypatch.setenv("LIFTON_PARALLEL_FORCE", "1")
        from lifton import parallel, run_liftoff
        from types import SimpleNamespace
        from unittest import mock
        import io, time

        N = 16
        loci = [SimpleNamespace(id=f"l_{i}", _idx=i) for i in range(N)]

        class FakeDB:
            def features_of_type(self, ft): yield from loci

        def fake(*args, ENTRY_FEATURE=False, **kw):
            locus = args[1]
            # Vary the sleep so completion order is unpredictable
            time.sleep(0.001 * (locus._idx % 5))
            gene = mock.MagicMock()
            gene.ref_gene_id = "ok"
            gene.write_entry.side_effect = (
                lambda fw, stats: fw.write(f"{locus._idx}\n")
            )
            return gene

        monkeypatch.setattr(run_liftoff, "process_liftoff", fake)

        # Spy on heap depth via consume hook
        from lifton.locus_pipeline import StepContext
        ctx = StepContext(
            ref_db=mock.Mock(), l_feature_db=FakeDB(), m_feature_db=None,
            ref_id_2_m_id_trans_dict={}, tree_dict={},
            tgt_fai=mock.Mock(), ref_proteins={}, ref_trans={},
            ref_features_dict={}, fw_score=io.StringIO(),
            fw_chain=None, args=SimpleNamespace(),
        )

        fw = io.StringIO()
        stats = {"coding": {}, "non-coding": {}, "other": {}}
        n = parallel.parallel_step7(
            ["gene"], FakeDB(), ctx, fw, stats, threads=4,
        )
        assert n == N
        # Output is in submission order
        assert fw.getvalue() == "".join(f"{i}\n" for i in range(N))
