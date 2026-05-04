"""Phase 9 — locus-major task wrapper + ordered-writer determinism.

Three test categories:

1. **Unit tests** for the data classes and the per-locus task.
2. **Ordered-writer tests** asserting deterministic emit order
   regardless of completion order (the byte-identity guarantee).
3. **Concurrency tests** asserting exception isolation, thread fan-out,
   and the parasail GIL-release behaviour.
"""

from __future__ import annotations

import io
import threading
import time
from concurrent.futures import ThreadPoolExecutor
from types import SimpleNamespace
from unittest import mock

import pytest

from lifton import locus_pipeline, parallel
from lifton.locus_pipeline import LocusResult, StepContext, consume, process_locus


# ---------------------------------------------------------------------------
# 1. Data-class semantics
# ---------------------------------------------------------------------------

class TestLocusResult:
    def test_index_round_trips(self):
        r = LocusResult(index=42, locus_id="g1")
        assert r.index == 42

    def test_emittable_false_when_gene_is_none(self):
        r = LocusResult(index=0, locus_id="g")
        assert r.emittable is False

    def test_emittable_false_when_ref_gene_id_missing(self):
        gene = SimpleNamespace(ref_gene_id=None)
        r = LocusResult(index=0, locus_id="g", lifton_gene=gene)
        assert r.emittable is False

    def test_emittable_false_when_error_set(self):
        gene = SimpleNamespace(ref_gene_id="ok")
        r = LocusResult(index=0, locus_id="g", lifton_gene=gene,
                        error=RuntimeError("bad"))
        assert r.emittable is False

    def test_emittable_true_for_valid_gene(self):
        gene = SimpleNamespace(ref_gene_id="ok")
        r = LocusResult(index=0, locus_id="g", lifton_gene=gene)
        assert r.emittable is True


# ---------------------------------------------------------------------------
# 2. process_locus()
# ---------------------------------------------------------------------------

def _fake_ctx() -> StepContext:
    """All-Mock StepContext; only used by tests that monkey-patch
    process_liftoff to a sentinel callable."""
    return StepContext(
        ref_db=mock.Mock(name="ref_db"),
        l_feature_db=mock.Mock(name="l_feature_db"),
        m_feature_db=None,
        ref_id_2_m_id_trans_dict={},
        tree_dict={},
        tgt_fai=mock.Mock(name="tgt_fai"),
        ref_proteins={}, ref_trans={}, ref_features_dict={},
        fw_score=io.StringIO(),
        fw_chain=None,
        args=SimpleNamespace(),
    )


class TestProcessLocus:
    def test_returns_locus_result_with_submission_index(self, monkeypatch):
        from lifton import run_liftoff
        sentinel = SimpleNamespace(ref_gene_id="g1")
        monkeypatch.setattr(run_liftoff, "process_liftoff",
                            lambda *a, **k: sentinel)
        locus = SimpleNamespace(id="locus_X")
        r = process_locus(7, locus, ctx=_fake_ctx())
        assert isinstance(r, LocusResult)
        assert r.index == 7
        assert r.locus_id == "locus_X"
        assert r.lifton_gene is sentinel
        assert r.error is None

    def test_packages_exception_instead_of_propagating(self, monkeypatch):
        from lifton import run_liftoff

        def boom(*a, **k):
            raise RuntimeError("synthetic failure")
        monkeypatch.setattr(run_liftoff, "process_liftoff", boom)

        r = process_locus(0, SimpleNamespace(id="bad"), ctx=_fake_ctx())
        assert r.lifton_gene is None
        assert isinstance(r.error, RuntimeError)
        assert "synthetic failure" in str(r.error)

    def test_handles_locus_without_id_attribute(self, monkeypatch):
        from lifton import run_liftoff
        monkeypatch.setattr(run_liftoff, "process_liftoff",
                            lambda *a, **k: None)
        # An object with no `.id` should fall back to the placeholder.
        r = process_locus(0, object(), ctx=_fake_ctx())
        assert r.locus_id == "<unknown>"

    def test_passes_ENTRY_FEATURE_True(self, monkeypatch):
        """The legacy serial loop calls process_liftoff with
        ``ENTRY_FEATURE=True``; the locus-pipeline must mirror that
        for byte-identity."""
        from lifton import run_liftoff
        captured = {}
        def fake(*args, ENTRY_FEATURE=False, **kwargs):
            captured["ENTRY_FEATURE"] = ENTRY_FEATURE
            return None
        monkeypatch.setattr(run_liftoff, "process_liftoff", fake)
        process_locus(0, SimpleNamespace(id="x"), ctx=_fake_ctx())
        assert captured["ENTRY_FEATURE"] is True


# ---------------------------------------------------------------------------
# 3. consume()
# ---------------------------------------------------------------------------

class TestConsume:
    def _stats(self):
        return {"coding": {}, "non-coding": {}, "other": {}}

    def test_skips_when_lifton_gene_is_none(self):
        fw = io.StringIO()
        stats = self._stats()
        r = LocusResult(index=0, locus_id="g")
        assert consume(r, fw, stats) is False
        assert fw.getvalue() == ""
        assert stats == {"coding": {}, "non-coding": {}, "other": {}}

    def test_skips_when_ref_gene_id_is_none(self):
        fw = io.StringIO()
        stats = self._stats()
        gene = mock.MagicMock()
        gene.ref_gene_id = None
        r = LocusResult(index=0, locus_id="g", lifton_gene=gene)
        assert consume(r, fw, stats) is False
        gene.write_entry.assert_not_called()

    def test_writes_when_emittable(self):
        fw = io.StringIO()
        stats = self._stats()
        gene = mock.MagicMock()
        gene.ref_gene_id = "g_ok"
        r = LocusResult(index=0, locus_id="g_ok", lifton_gene=gene)
        assert consume(r, fw, stats) is True
        gene.write_entry.assert_called_once_with(fw, stats)

    def test_logs_error_and_does_not_write(self, capsys):
        fw = io.StringIO()
        stats = self._stats()
        r = LocusResult(index=0, locus_id="bad",
                        error=RuntimeError("kaboom"))
        assert consume(r, fw, stats) is False
        captured = capsys.readouterr()
        # logger writes to stderr/stdout depending on level; just
        # assert the locus id appears somewhere.
        assert "bad" in (captured.err + captured.out)


# ---------------------------------------------------------------------------
# 4. _iter_loci ordering
# ---------------------------------------------------------------------------

class TestIterLoci:
    def test_visits_features_then_loci_in_db_order(self):
        # Build a fake l_feature_db that returns predictable loci
        loci_by_type = {
            "gene": [SimpleNamespace(id="g1"), SimpleNamespace(id="g2")],
            "pseudogene": [SimpleNamespace(id="pg1")],
        }

        class FakeDB:
            def features_of_type(self, ft):
                yield from loci_by_type[ft]

        ordered = list(parallel._iter_loci(["gene", "pseudogene"], FakeDB()))
        assert [f for (f, _l) in ordered] == ["gene", "gene", "pseudogene"]
        assert [l.id for (_f, l) in ordered] == ["g1", "g2", "pg1"]


# ---------------------------------------------------------------------------
# 5. Ordered-writer determinism
# ---------------------------------------------------------------------------

def _build_step7_inputs(n_loci=10, slow_indices=()):
    """Build a fake l_feature_db + a fake `process_liftoff` that
    sleeps for `slow_indices` so worker completion order differs from
    submission order. This is the kernel of the determinism gate."""
    loci = [SimpleNamespace(id=f"locus_{i}", _idx=i) for i in range(n_loci)]

    class FakeDB:
        def features_of_type(self, ft):
            yield from loci

    def fake_process(*args, ENTRY_FEATURE=False, **kw):
        # `args[1]` is the locus argument
        locus = args[1]
        if locus._idx in slow_indices:
            time.sleep(0.05)
        gene = mock.MagicMock()
        gene.ref_gene_id = "ok"
        # Capture the locus' submission index in the gene so consume()
        # can write a deterministic record.
        def _write(fw, stats):
            fw.write(f"{locus._idx}\n")
        gene.write_entry.side_effect = _write
        return gene

    return loci, FakeDB(), fake_process


class TestOrderedWriterDeterminism:
    def test_serial_path_emits_in_submission_order(self, monkeypatch):
        from lifton import run_liftoff
        _, db, fake = _build_step7_inputs(n_loci=8)
        monkeypatch.setattr(run_liftoff, "process_liftoff", fake)

        fw = io.StringIO()
        stats = {"coding": {}, "non-coding": {}, "other": {}}
        ctx = _fake_ctx()
        n = parallel.parallel_step7(
            ["gene"], db, ctx, fw, stats, threads=1,
        )
        assert n == 8
        assert fw.getvalue() == "".join(f"{i}\n" for i in range(8))

    @pytest.mark.parametrize("threads", [2, 4, 8])
    def test_parallel_path_emits_in_submission_order(self, monkeypatch, threads):
        monkeypatch.setenv("LIFTON_PARALLEL_FORCE", "1")
        from lifton import run_liftoff
        # Slow-down indices that DO complete out of order with threads.
        _, db, fake = _build_step7_inputs(n_loci=12, slow_indices={0, 3, 7})
        monkeypatch.setattr(run_liftoff, "process_liftoff", fake)

        fw = io.StringIO()
        stats = {"coding": {}, "non-coding": {}, "other": {}}
        ctx = _fake_ctx()
        n = parallel.parallel_step7(
            ["gene"], db, ctx, fw, stats, threads=threads,
        )
        assert n == 12
        # Output MUST still be in submission order regardless of
        # completion order.
        assert fw.getvalue() == "".join(f"{i}\n" for i in range(12))

    def test_parallel_byte_identical_to_serial(self, monkeypatch):
        monkeypatch.setenv("LIFTON_PARALLEL_FORCE", "1")
        from lifton import run_liftoff
        _, db, fake = _build_step7_inputs(n_loci=15, slow_indices={2, 5, 11})
        monkeypatch.setattr(run_liftoff, "process_liftoff", fake)

        # Serial baseline
        fw1 = io.StringIO()
        stats1 = {"coding": {}, "non-coding": {}, "other": {}}
        parallel.parallel_step7(["gene"], db, _fake_ctx(), fw1, stats1, threads=1)
        # Parallel
        fw4 = io.StringIO()
        stats4 = {"coding": {}, "non-coding": {}, "other": {}}
        parallel.parallel_step7(["gene"], db, _fake_ctx(), fw4, stats4, threads=4)
        assert fw1.getvalue() == fw4.getvalue()


# ---------------------------------------------------------------------------
# 6. Exception isolation
# ---------------------------------------------------------------------------

class TestExceptionIsolation:
    def test_one_bad_locus_does_not_abort_siblings(self, monkeypatch):
        monkeypatch.setenv("LIFTON_PARALLEL_FORCE", "1")
        from lifton import run_liftoff

        loci = [SimpleNamespace(id=f"l_{i}", _idx=i) for i in range(6)]

        class FakeDB:
            def features_of_type(self, ft): yield from loci

        def fake(*args, ENTRY_FEATURE=False, **kw):
            locus = args[1]
            if locus._idx == 3:
                raise RuntimeError("locus 3 is poison")
            gene = mock.MagicMock()
            gene.ref_gene_id = "ok"
            gene.write_entry.side_effect = (
                lambda fw, stats: fw.write(f"{locus._idx}\n")
            )
            return gene

        monkeypatch.setattr(run_liftoff, "process_liftoff", fake)
        fw = io.StringIO()
        stats = {"coding": {}, "non-coding": {}, "other": {}}
        n = parallel.parallel_step7(
            ["gene"], FakeDB(), _fake_ctx(), fw, stats, threads=4,
        )
        # processed_features still counts every locus (the legacy
        # counter increments unconditionally), AND the surviving
        # outputs are emitted in order with locus 3 silently skipped.
        assert n == 6
        emitted = fw.getvalue().splitlines()
        assert "3" not in emitted
        assert emitted == ["0", "1", "2", "4", "5"]


# ---------------------------------------------------------------------------
# 7. Thread fan-out — workers really do run in parallel
# ---------------------------------------------------------------------------

class TestThreadFanout:
    def test_workers_run_concurrently(self, monkeypatch):
        """Use a barrier sized N — if all N workers can pass it, they
        were running concurrently. Threads are required (the GIL is
        released by `time.sleep`).

        Override the Phase 9 thread-safety guard via the env-var
        escape hatch so the parallel path is taken even though the
        test uses Mock FeatureDBs (no real backend involvement)."""
        monkeypatch.setenv("LIFTON_PARALLEL_FORCE", "1")
        from lifton import run_liftoff
        N = 4
        loci = [SimpleNamespace(id=f"l_{i}", _idx=i) for i in range(N)]

        class FakeDB:
            def features_of_type(self, ft): yield from loci

        barrier = threading.Barrier(N, timeout=2.0)

        def fake(*args, ENTRY_FEATURE=False, **kw):
            barrier.wait()  # blocks until N workers arrive
            gene = mock.MagicMock()
            gene.ref_gene_id = "ok"
            gene.write_entry.side_effect = (
                lambda fw, stats: fw.write(f"{args[1]._idx}\n")
            )
            return gene

        monkeypatch.setattr(run_liftoff, "process_liftoff", fake)
        fw = io.StringIO()
        stats = {"coding": {}, "non-coding": {}, "other": {}}
        # If threads weren't running concurrently, the barrier would
        # time out and Brokenfit raise.
        parallel.parallel_step7(
            ["gene"], FakeDB(), _fake_ctx(), fw, stats, threads=N,
        )
        # And output remains deterministic.
        assert fw.getvalue() == "0\n1\n2\n3\n"


# ---------------------------------------------------------------------------
# 8. Empty input + threads=1 fallback
# ---------------------------------------------------------------------------

class TestEdgeCases:
    def test_empty_db_returns_zero(self, monkeypatch):
        monkeypatch.setenv("LIFTON_PARALLEL_FORCE", "1")
        from lifton import run_liftoff
        monkeypatch.setattr(run_liftoff, "process_liftoff",
                            lambda *a, **k: None)

        class EmptyDB:
            def features_of_type(self, ft): return iter([])

        fw = io.StringIO()
        stats = {"coding": {}, "non-coding": {}, "other": {}}
        n = parallel.parallel_step7(
            ["gene"], EmptyDB(), _fake_ctx(), fw, stats, threads=4,
        )
        assert n == 0
        assert fw.getvalue() == ""

    def test_threads_zero_falls_back_to_serial(self, monkeypatch):
        from lifton import run_liftoff
        _, db, fake = _build_step7_inputs(n_loci=3)
        monkeypatch.setattr(run_liftoff, "process_liftoff", fake)
        fw = io.StringIO()
        stats = {"coding": {}, "non-coding": {}, "other": {}}
        # threads=0 should be treated as serial (no executor created).
        n = parallel.parallel_step7(
            ["gene"], db, _fake_ctx(), fw, stats, threads=0,
        )
        assert n == 3
        assert fw.getvalue() == "0\n1\n2\n"

    def test_threads_none_falls_back_to_serial(self, monkeypatch):
        from lifton import run_liftoff
        _, db, fake = _build_step7_inputs(n_loci=3)
        monkeypatch.setattr(run_liftoff, "process_liftoff", fake)
        fw = io.StringIO()
        stats = {"coding": {}, "non-coding": {}, "other": {}}
        n = parallel.parallel_step7(
            ["gene"], db, _fake_ctx(), fw, stats, threads=None,
        )
        assert n == 3


# ---------------------------------------------------------------------------
# 9. Stats dict mutation
# ---------------------------------------------------------------------------

class TestStatsDictMutation:
    def test_stats_dict_updated_on_emit(self, monkeypatch):
        monkeypatch.setenv("LIFTON_PARALLEL_FORCE", "1")
        from lifton import run_liftoff

        def fake(*a, ENTRY_FEATURE=False, **k):
            gene = mock.MagicMock()
            gene.ref_gene_id = "ok"
            def _write(fw, stats):
                stats["coding"]["tx_X"] = stats["coding"].get("tx_X", 0) + 1
                fw.write("x\n")
            gene.write_entry.side_effect = _write
            return gene

        loci = [SimpleNamespace(id=f"l_{i}", _idx=i) for i in range(3)]

        class FakeDB:
            def features_of_type(self, ft): yield from loci

        monkeypatch.setattr(run_liftoff, "process_liftoff", fake)
        fw = io.StringIO()
        stats = {"coding": {}, "non-coding": {}, "other": {}}
        parallel.parallel_step7(
            ["gene"], FakeDB(), _fake_ctx(), fw, stats, threads=4,
        )
        assert stats["coding"]["tx_X"] == 3


# ---------------------------------------------------------------------------
# 10. CLI flag plumbing
# ---------------------------------------------------------------------------

class TestLocusPipelineCLI:
    def test_flag_default_is_false(self):
        from lifton import lifton as lifton_main
        argv = ["t.fa", "r.fa", "-g", "r.gff3"]
        args = lifton_main.parse_args(argv)
        assert args.locus_pipeline is False

    def test_flag_set_to_true(self):
        from lifton import lifton as lifton_main
        argv = ["t.fa", "r.fa", "-g", "r.gff3", "--locus-pipeline"]
        args = lifton_main.parse_args(argv)
        assert args.locus_pipeline is True

    def test_threads_default_is_one(self):
        from lifton import lifton as lifton_main
        argv = ["t.fa", "r.fa", "-g", "r.gff3"]
        args = lifton_main.parse_args(argv)
        assert int(args.threads) == 1

    def test_threads_can_be_overridden(self):
        from lifton import lifton as lifton_main
        argv = ["t.fa", "r.fa", "-g", "r.gff3", "-t", "4"]
        args = lifton_main.parse_args(argv)
        assert int(args.threads) == 4

    def test_help_text_mentions_locus_pipeline(self, capsys):
        from lifton import lifton as lifton_main
        with pytest.raises(SystemExit):
            lifton_main.parse_args(["--help"])
        out = capsys.readouterr().out
        assert "--locus-pipeline" in out


# ---------------------------------------------------------------------------
# 11. parasail GIL release smoke test
# ---------------------------------------------------------------------------

class TestParasailGILRelease:
    def test_concurrent_alignment_completes(self):
        """Best-effort proof that parasail runs from multiple threads
        without serialising. Not a strict timing test (CI variance);
        we just assert that 8 concurrent alignments complete without
        deadlock or exception."""
        import parasail
        matrix = parasail.Matrix("blosum62")

        def _align(_):
            return parasail.nw_trace_scan_sat(
                "MAGTAAA", "MAGTBBB", 11, 1, matrix,
            )

        with ThreadPoolExecutor(max_workers=8) as ex:
            results = list(ex.map(_align, range(32)))
        assert len(results) == 32
        for r in results:
            # Each result has a traceback — proves the kernel ran.
            assert hasattr(r, "traceback")
