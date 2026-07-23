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

import gffutils
import pytest

from lifton import locus_pipeline, parallel
from lifton.locus_pipeline import (
    LocusResult,
    HierarchyBatchLoader,
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


class _RecordingBatchedLDB:
    """Small hierarchy backend that records only set-based DB reads."""
    __module__ = "lifton.gffbase.interface"

    def __init__(self, direct_children):
        self._direct_children = direct_children
        self.batch_calls = []

    @staticmethod
    def _id(feature):
        return getattr(feature, "id", feature)

    def _select(self, anchor, *, featuretype=None, level=None,
                order_by=None):
        anchor_id = self._id(anchor)
        direct = list(self._direct_children.get(anchor_id, ()))
        if level == 1:
            candidates = direct
        else:
            candidates = []
            pending = direct
            visited = set()
            while pending:
                child = pending.pop(0)
                child_id = self._id(child)
                if child_id in visited:
                    continue
                visited.add(child_id)
                candidates.append(child)
                pending.extend(self._direct_children.get(child_id, ()))
        if featuretype is not None:
            accepted = (
                set(featuretype) if isinstance(featuretype, tuple)
                else {featuretype}
            )
            candidates = [
                child for child in candidates
                if child.featuretype in accepted
            ]
        if order_by == "start":
            candidates.sort(key=lambda child: (child.start, child.id))
        return tuple(candidates)

    def children(self, anchor, *, featuretype=None, level=None,
                 order_by=None):
        return iter(self._select(
            anchor, featuretype=featuretype, level=level,
            order_by=order_by,
        ))

    def children_batched_features(self, anchors, *, featuretype=None,
                                  level=None, order_by=None):
        anchor_ids = tuple(self._id(anchor) for anchor in anchors)
        self.batch_calls.append(
            (anchor_ids, featuretype, level, order_by)
        )
        return {
            anchor_id: self._select(
                anchor_id, featuretype=featuretype, level=level,
                order_by=order_by,
            )
            for anchor_id in anchor_ids
        }


def _build_ctx_for_materialise(monkeypatch=None):
    locus = SimpleNamespace(id="gene1", start=100, end=200, seqid="chr1",
                            strand="+")
    exon1 = SimpleNamespace(
        id="exon1", start=100, end=150, featuretype="exon",
    )
    exon2 = SimpleNamespace(
        id="exon2", start=160, end=200, featuretype="exon",
    )
    cds1 = SimpleNamespace(
        id="cds1", start=100, end=150, featuretype="CDS",
    )
    cds2 = SimpleNamespace(
        id="cds2", start=160, end=200, featuretype="CDS",
    )
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

    def test_batched_walk_queries_only_runtime_reachable_branches(self):
        gene = SimpleNamespace(
            id="gene1", start=1, end=500, featuretype="gene",
        )
        tx1 = SimpleNamespace(
            id="tx1", start=10, end=200, featuretype="mRNA",
        )
        tx2 = SimpleNamespace(
            id="tx2", start=250, end=490, featuretype="mRNA",
        )
        exon1 = SimpleNamespace(
            id="exon1", start=10, end=90, featuretype="exon",
        )
        cds1 = SimpleNamespace(
            id="cds1", start=20, end=80, featuretype="CDS",
        )
        stop1 = SimpleNamespace(
            id="stop1", start=81, end=83, featuretype="stop_codon",
        )
        exon2 = SimpleNamespace(
            id="exon2", start=250, end=490, featuretype="exon",
        )
        cds2 = SimpleNamespace(
            id="cds2", start=260, end=480, featuretype="CDS",
        )
        db = _RecordingBatchedLDB({
            "gene1": (tx1, tx2),
            "tx1": (exon1, cds1, stop1),
            "tx2": (exon2, cds2),
        })
        ctx = StepContext(
            ref_db=_FakeRefDB({}), l_feature_db=db, m_feature_db=None,
            ref_id_2_m_id_trans_dict={}, tree_dict={},
            tgt_fai=mock.Mock(), ref_proteins={}, ref_trans={},
            ref_features_dict={}, fw_score=io.StringIO(), fw_chain=None,
            args=SimpleNamespace(native=True, threads=8),
        )

        payload = materialise_locus(0, gene, ctx)

        # One exon probe plus one branch-specific query at each hierarchy
        # depth. No all-exon fallback, gene-level CDS query, or leaf probes.
        assert db.batch_calls == [
            (("gene1",), "exon", 1, "start"),
            (("gene1",), None, 1, None),
            (("tx1", "tx2"), "exon", 1, "start"),
            (
                ("tx1", "tx2"),
                ("CDS", "stop_codon"),
                None,
                "start",
            ),
        ]
        assert list(payload.feature_cache) == ["gene1", "tx1", "tx2"]
        assert [child.id for child in payload.cds_children] == [
            "cds1", "cds2",
        ]
        assert [child.id for child in payload.cds_stop_children] == [
            "cds1", "stop1", "cds2",
        ]

        def traversal_bytes(hierarchy_db):
            rows = []

            def walk(feature):
                direct_exons = list(hierarchy_db.children(
                    feature, featuretype="exon", level=1,
                    order_by="start",
                ))
                if direct_exons:
                    all_exons = hierarchy_db.children(
                        feature, featuretype="exon", order_by="start",
                    )
                    cds_stop = hierarchy_db.children(
                        feature, featuretype=("CDS", "stop_codon"),
                        order_by="start",
                    )
                    rows.append(
                        f"{feature.id}|exons="
                        f"{','.join(child.id for child in all_exons)}|cds="
                        f"{','.join(child.id for child in cds_stop)}\n"
                    )
                    return
                rows.append(f"{feature.id}|container\n")
                for child in hierarchy_db.children(feature, level=1):
                    walk(child)

            walk(gene)
            return "".join(rows).encode()

        proxy = locus_pipeline._LFeatureDbProxy(payload.feature_cache)
        assert traversal_bytes(proxy) == traversal_bytes(db)


class TestHierarchyBatchLoader:
    def test_missing_reference_keeps_valid_sibling_transcripts(self):
        valid = {
            "tx1": SimpleNamespace(id="tx1"),
            "tx2": SimpleNamespace(id="tx2"),
        }

        class ScalarDB:
            def __getitem__(self, feature_id):
                try:
                    return valid[feature_id]
                except KeyError:
                    raise gffutils.exceptions.FeatureNotFoundError(
                        feature_id
                    )

        loaded = HierarchyBatchLoader(ScalarDB()).features_many(
            ["tx1", "missing", "tx2"]
        )

        assert list(loaded) == ["tx1", "tx2"]
        assert loaded["tx1"] is valid["tx1"]
        assert loaded["tx2"] is valid["tx2"]

    def test_unexpected_lookup_error_is_not_suppressed(self):
        class BrokenDB:
            def __getitem__(self, _feature_id):
                raise IndexError("corrupt backend row")

        with pytest.raises(IndexError, match="corrupt backend row"):
            HierarchyBatchLoader(BrokenDB()).features_many(["tx1"])


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
    def test_reopenable_dispatch_materialises_lazily_on_workers(
            self, monkeypatch):
        """Safe connection-less fixtures take the fused lazy worker path."""
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
        # Every payload is materialised, but no whole-root parent snapshot is
        # required before the executor starts.
        assert len(seen_payloads) == 5
        assert all(tid != threading.current_thread().ident for tid in seen_threads)
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
