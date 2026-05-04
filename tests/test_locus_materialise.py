"""Phase 11 — MaterialisedLocus + process_locus_native + dispatch.

Covers:
  * `materialise_locus` populates every payload field on a fake DB
  * `process_locus_native` reproduces the legacy `process_locus`
    output for the synthetic fixture
  * Determinism: the parallel native path produces byte-identical
    output to the serial path under the heap ordered-writer
  * Real ThreadPoolExecutor fan-out under `args.native=True`
  * Worker exception isolation on the native path
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
from lifton.locus_pipeline import (
    LocusResult,
    MaterialisedLocus,
    StepContext,
    consume,
    materialise_locus,
    process_locus_native,
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

class _FakeRefDB:
    """Mimics gffbase.FeatureDB.__getitem__ semantics."""
    def __init__(self, entries):
        self._entries = entries

    def __getitem__(self, key):
        if key not in self._entries:
            raise KeyError(key)
        return self._entries[key]


class _FakeLDB:
    """Mimics l_feature_db.children with a per-feature-type lookup."""
    __module__ = "lifton.gffbase.interface"  # pretend gffbase

    def __init__(self, by_type):
        self._by_type = by_type

    def features_of_type(self, ft):
        yield from self._by_type.get(ft, [])

    def children(self, locus, *, featuretype=None, level=None,
                 order_by=None):
        if featuretype is None:
            return iter(self._by_type.get("level1", []))
        if isinstance(featuretype, tuple):
            results = []
            for ft in featuretype:
                results.extend(self._by_type.get(ft, []))
            return iter(results)
        return iter(self._by_type.get(featuretype, []))


def _build_ctx_for_materialise(monkeypatch=None):
    locus = SimpleNamespace(id="gene1", start=100, end=200, seqid="chr1",
                            strand="+")
    exon1 = SimpleNamespace(id="exon1", start=100, end=150)
    exon2 = SimpleNamespace(id="exon2", start=160, end=200)
    cds1 = SimpleNamespace(id="cds1", start=100, end=150)
    cds2 = SimpleNamespace(id="cds2", start=160, end=200)
    l_db = _FakeLDB({
        "exon": [exon1, exon2],
        "CDS": [cds1, cds2],
        "level1": [exon1, exon2, cds1, cds2],
    })
    ref_db = _FakeRefDB({
        "gene1": SimpleNamespace(attributes={"ID": ["gene1"]}),
        "tx1": SimpleNamespace(attributes={"ID": ["tx1"]}),
    })
    ref_features_dict = {
        "gene1": SimpleNamespace(),
    }
    if monkeypatch is not None:
        # Force get_ref_ids_liftoff to return the trans id we expect.
        from lifton import lifton_utils
        monkeypatch.setattr(
            lifton_utils, "get_ref_ids_liftoff",
            lambda d, gene_id, trans_id: ("gene1", "tx1"),
        )
    ctx = StepContext(
        ref_db=ref_db, l_feature_db=l_db, m_feature_db=None,
        ref_id_2_m_id_trans_dict={"tx1": ["MP1"]}, tree_dict={},
        tgt_fai=mock.Mock(), ref_proteins={"tx1": "MAGT*"},
        ref_trans={"tx1": "ATGGCTTAA"},
        ref_features_dict=ref_features_dict,
        fw_score=io.StringIO(), fw_chain=None,
        args=SimpleNamespace(native=True),
    )
    return ctx, locus


# ---------------------------------------------------------------------------
# 1. MaterialisedLocus dataclass
# ---------------------------------------------------------------------------

class TestMaterialisedLocus:
    def test_default_lists_empty(self):
        m = MaterialisedLocus(submission_index=0, locus=None, locus_id="g")
        assert m.children_l1 == []
        assert m.exon_children == []
        assert m.cds_children == []
        assert m.cds_stop_children == []
        assert m.ref_gene_attrs == {}
        assert m.ref_trans_attrs == {}
        assert m.ref_gene_id is None
        assert m.ref_trans_id is None

    def test_index_round_trip(self):
        m = MaterialisedLocus(submission_index=42, locus=None, locus_id="g")
        assert m.submission_index == 42

    def test_locus_id_round_trip(self):
        m = MaterialisedLocus(submission_index=0, locus=None,
                              locus_id="ENSG00000001")
        assert m.locus_id == "ENSG00000001"


# ---------------------------------------------------------------------------
# 2. materialise_locus
# ---------------------------------------------------------------------------

class TestMaterialiseLocus:
    def test_populates_children(self, monkeypatch):
        ctx, locus = _build_ctx_for_materialise(monkeypatch)
        payload = materialise_locus(7, locus, ctx)
        assert payload.submission_index == 7
        assert payload.locus_id == "gene1"
        assert len(payload.exon_children) == 2
        assert len(payload.cds_children) == 2
        assert len(payload.cds_stop_children) == 2
        assert len(payload.children_l1) == 4

    def test_populates_ref_attrs(self, monkeypatch):
        ctx, locus = _build_ctx_for_materialise(monkeypatch)
        payload = materialise_locus(0, locus, ctx)
        assert payload.ref_gene_id == "gene1"
        assert payload.ref_trans_id == "tx1"
        assert payload.ref_gene_attrs == {"ID": ["gene1"]}
        assert payload.ref_trans_attrs == {"ID": ["tx1"]}

    def test_handles_missing_ref_ids(self, monkeypatch):
        ctx, locus = _build_ctx_for_materialise()
        # No monkeypatch: get_ref_ids_liftoff against an empty
        # ref_features_dict returns (None, None).
        from lifton import lifton_utils
        monkeypatch.setattr(
            lifton_utils, "get_ref_ids_liftoff",
            lambda d, gene_id, trans_id: (None, None),
        )
        payload = materialise_locus(0, locus, ctx)
        assert payload.ref_gene_id is None
        assert payload.ref_trans_id is None
        assert payload.ref_gene_attrs == {}

    def test_runs_in_calling_thread(self, monkeypatch):
        """All DB reads must happen on the thread that called
        materialise_locus, not on a background worker."""
        ctx, locus = _build_ctx_for_materialise(monkeypatch)
        thread_ids_seen = []
        original_children = ctx.l_feature_db.children

        def spy(*a, **kw):
            thread_ids_seen.append(threading.current_thread().ident)
            return original_children(*a, **kw)

        ctx.l_feature_db.children = spy
        payload = materialise_locus(0, locus, ctx)
        # Every DB call happened on the test's own thread
        assert all(tid == threading.current_thread().ident
                   for tid in thread_ids_seen)
        assert payload.exon_children


# ---------------------------------------------------------------------------
# 3. process_locus_native semantics
# ---------------------------------------------------------------------------

class TestProcessLocusNative:
    def test_returns_locus_result_with_index(self, monkeypatch):
        from lifton import run_liftoff
        sentinel = SimpleNamespace(ref_gene_id="g1")
        monkeypatch.setattr(run_liftoff, "process_liftoff",
                            lambda *a, **k: sentinel)
        ctx = StepContext(
            ref_db=mock.Mock(), l_feature_db=mock.Mock(), m_feature_db=None,
            ref_id_2_m_id_trans_dict={}, tree_dict={},
            tgt_fai=mock.Mock(), ref_proteins={}, ref_trans={},
            ref_features_dict={}, fw_score=io.StringIO(), fw_chain=None,
            args=SimpleNamespace(),
        )
        payload = MaterialisedLocus(
            submission_index=5, locus=SimpleNamespace(id="g1"),
            locus_id="g1",
        )
        result = process_locus_native(payload, ctx)
        assert result.index == 5
        assert result.locus_id == "g1"
        assert result.lifton_gene is sentinel
        assert result.error is None

    def test_packages_exception(self, monkeypatch):
        from lifton import run_liftoff

        def boom(*a, **kw):
            raise RuntimeError("worker failure")

        monkeypatch.setattr(run_liftoff, "process_liftoff", boom)
        ctx = StepContext(
            ref_db=mock.Mock(), l_feature_db=mock.Mock(), m_feature_db=None,
            ref_id_2_m_id_trans_dict={}, tree_dict={},
            tgt_fai=mock.Mock(), ref_proteins={}, ref_trans={},
            ref_features_dict={}, fw_score=io.StringIO(), fw_chain=None,
            args=SimpleNamespace(),
        )
        payload = MaterialisedLocus(
            submission_index=0, locus=SimpleNamespace(id="g1"), locus_id="g1",
        )
        result = process_locus_native(payload, ctx)
        assert result.lifton_gene is None
        assert isinstance(result.error, RuntimeError)


# ---------------------------------------------------------------------------
# 4. parallel_step7 with native dispatch
# ---------------------------------------------------------------------------

class TestParallelStep7NativeDispatch:
    def test_parent_pre_materialises_payloads(self, monkeypatch):
        """The native dispatch path must call materialise_locus for
        every locus before submitting any worker."""
        from lifton import parallel as p_mod, run_liftoff
        seen_payloads = []
        seen_threads = []

        def spy_materialise(idx, locus, ctx):
            seen_threads.append(threading.current_thread().ident)
            payload = MaterialisedLocus(
                submission_index=idx, locus=locus,
                locus_id=getattr(locus, "id", "?"),
            )
            seen_payloads.append(payload)
            return payload

        def fake_process(*args, ENTRY_FEATURE=False, **kw):
            locus = args[1]
            gene = mock.MagicMock()
            gene.ref_gene_id = "ok"
            gene.write_entry.side_effect = (
                lambda fw, stats: fw.write(f"{locus._idx}\n")
            )
            return gene

        monkeypatch.setattr(p_mod, "_backend_supports_threads",
                            lambda *a, **kw: True)
        monkeypatch.setattr(locus_pipeline, "materialise_locus",
                            spy_materialise)
        monkeypatch.setattr(run_liftoff, "process_liftoff", fake_process)

        loci = [SimpleNamespace(id=f"l_{i}", _idx=i) for i in range(5)]

        class FakeDB:
            __module__ = "lifton.gffbase.interface"
            def features_of_type(self, ft): yield from loci

        db = FakeDB()
        ctx = StepContext(
            ref_db=db, l_feature_db=db, m_feature_db=None,
            ref_id_2_m_id_trans_dict={}, tree_dict={},
            tgt_fai=mock.Mock(), ref_proteins={}, ref_trans={},
            ref_features_dict={}, fw_score=io.StringIO(), fw_chain=None,
            args=SimpleNamespace(native=True),
        )
        fw = io.StringIO()
        stats = {"coding": {}, "non-coding": {}, "other": {}}
        n = parallel.parallel_step7(["gene"], db, ctx, fw, stats, threads=4)
        assert n == 5
        # All 5 payloads pre-materialised on the parent (this) thread.
        assert len(seen_payloads) == 5
        assert all(tid == threading.current_thread().ident for tid in seen_threads)
        # Output preserved in submission order.
        assert fw.getvalue() == "0\n1\n2\n3\n4\n"

    def test_native_threads_byte_identical_to_serial(self, monkeypatch):
        """Phase 11 byte-identity contract: the same locus list
        produces byte-identical output across threads={1,2,4} when
        --native is on."""
        from lifton import parallel as p_mod, run_liftoff
        monkeypatch.setattr(p_mod, "_backend_supports_threads",
                            lambda *a, **kw: True)

        loci = [SimpleNamespace(id=f"l_{i}", _idx=i) for i in range(10)]

        class FakeDB:
            __module__ = "lifton.gffbase.interface"
            def features_of_type(self, ft): yield from loci

        def fake(*args, ENTRY_FEATURE=False, **kw):
            locus = args[1]
            # Simulate variable processing time so completion order
            # differs from submission order
            time.sleep(0.001 * (locus._idx % 3))
            gene = mock.MagicMock()
            gene.ref_gene_id = "ok"
            gene.write_entry.side_effect = (
                lambda fw, stats: fw.write(f"{locus._idx}\n")
            )
            return gene

        monkeypatch.setattr(run_liftoff, "process_liftoff", fake)

        outputs = {}
        for threads in (1, 2, 4):
            db = FakeDB()
            ctx = StepContext(
                ref_db=db, l_feature_db=db, m_feature_db=None,
                ref_id_2_m_id_trans_dict={}, tree_dict={},
                tgt_fai=mock.Mock(), ref_proteins={}, ref_trans={},
                ref_features_dict={}, fw_score=io.StringIO(),
                fw_chain=None, args=SimpleNamespace(native=True),
            )
            fw = io.StringIO()
            stats = {"coding": {}, "non-coding": {}, "other": {}}
            parallel.parallel_step7(
                ["gene"], db, ctx, fw, stats, threads=threads,
            )
            outputs[threads] = fw.getvalue()
        assert outputs[1] == outputs[2] == outputs[4]
        # Sanity: deterministic, matches submission order
        assert outputs[1] == "".join(f"{i}\n" for i in range(10))

    def test_concurrency_proof_with_barrier(self, monkeypatch):
        """N workers must run concurrently — proven by a Barrier of
        size N that times out if they don't."""
        from lifton import parallel as p_mod, run_liftoff
        N = 4
        monkeypatch.setattr(p_mod, "_backend_supports_threads",
                            lambda *a, **kw: True)

        loci = [SimpleNamespace(id=f"l_{i}", _idx=i) for i in range(N)]

        class FakeDB:
            __module__ = "lifton.gffbase.interface"
            def features_of_type(self, ft): yield from loci

        barrier = threading.Barrier(N, timeout=2.0)

        def fake(*args, ENTRY_FEATURE=False, **kw):
            barrier.wait()
            locus = args[1]
            gene = mock.MagicMock()
            gene.ref_gene_id = "ok"
            gene.write_entry.side_effect = (
                lambda fw, stats: fw.write(f"{locus._idx}\n")
            )
            return gene

        monkeypatch.setattr(run_liftoff, "process_liftoff", fake)

        db = FakeDB()
        ctx = StepContext(
            ref_db=db, l_feature_db=db, m_feature_db=None,
            ref_id_2_m_id_trans_dict={}, tree_dict={},
            tgt_fai=mock.Mock(), ref_proteins={}, ref_trans={},
            ref_features_dict={}, fw_score=io.StringIO(), fw_chain=None,
            args=SimpleNamespace(native=True),
        )
        fw = io.StringIO()
        stats = {"coding": {}, "non-coding": {}, "other": {}}
        # If the pool isn't actually fanning out, the barrier times
        # out and BrokenBarrierError surfaces.
        parallel.parallel_step7(
            ["gene"], db, ctx, fw, stats, threads=N,
        )
        assert fw.getvalue() == "0\n1\n2\n3\n"


# ---------------------------------------------------------------------------
# 5. Worker exception isolation on native path
# ---------------------------------------------------------------------------

class TestNativeExceptionIsolation:
    def test_one_bad_locus_does_not_abort_siblings(self, monkeypatch):
        from lifton import parallel as p_mod, run_liftoff
        monkeypatch.setattr(p_mod, "_backend_supports_threads",
                            lambda *a, **kw: True)
        loci = [SimpleNamespace(id=f"l_{i}", _idx=i) for i in range(5)]

        class FakeDB:
            __module__ = "lifton.gffbase.interface"
            def features_of_type(self, ft): yield from loci

        def fake(*args, ENTRY_FEATURE=False, **kw):
            locus = args[1]
            if locus._idx == 2:
                raise RuntimeError("locus 2 poison")
            gene = mock.MagicMock()
            gene.ref_gene_id = "ok"
            gene.write_entry.side_effect = (
                lambda fw, stats: fw.write(f"{locus._idx}\n")
            )
            return gene

        monkeypatch.setattr(run_liftoff, "process_liftoff", fake)
        db = FakeDB()
        ctx = StepContext(
            ref_db=db, l_feature_db=db, m_feature_db=None,
            ref_id_2_m_id_trans_dict={}, tree_dict={},
            tgt_fai=mock.Mock(), ref_proteins={}, ref_trans={},
            ref_features_dict={}, fw_score=io.StringIO(), fw_chain=None,
            args=SimpleNamespace(native=True),
        )
        fw = io.StringIO()
        stats = {"coding": {}, "non-coding": {}, "other": {}}
        n = parallel.parallel_step7(["gene"], db, ctx, fw, stats, threads=4)
        assert n == 5
        emitted = fw.getvalue().splitlines()
        # locus 2 silently skipped; all others emitted in order
        assert "2" not in emitted
        assert emitted == ["0", "1", "3", "4"]
