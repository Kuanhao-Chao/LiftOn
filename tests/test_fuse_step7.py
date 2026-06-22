"""Iteration 10 — fuse Step-7 materialise + process into one pipelined pool.

Pre-Iteration-10, the parallel Step-7 path ran TWO sequential phases with a
hard barrier: a 4-thread prefetcher pool built ALL `MaterialisedLocus`
payloads (SQLite-bound), then a separate N-thread worker pool ran ALL the
parasail processing. Iteration 10 fuses them — each worker materialises its
OWN locus (thread-local DB) then immediately processes it, so the SQLite I/O
overlaps the GIL-released parasail.

These tests pin the new contract on the on-disk default backend (where
`_ThreadLocalCtxFactory.viable` is True → the fused path is taken):
  1. the fused run constructs exactly ONE ThreadPoolExecutor (prefetcher gone),
     and is byte-identical to serial;
  2. `LIFTON_FUSE_STEP7=0` restores the two-phase path (TWO pools) and is
     byte-identical to the fused default and to serial;
  3. a failure in the materialise HALF becomes a per-locus error LocusResult,
     not a pool crash — siblings still emit in submission order.

Hermeticity: the integration_workspace fixture supplies pre-baked
liftoff/miniprot GFF3s via -L / -M, so minimap2 / miniprot are never invoked.
"""

from __future__ import annotations

import io
from types import SimpleNamespace
from unittest import mock

import pytest

from tests.test_integration_pipeline import (  # noqa: F401
    integration_workspace,
    hermetic_pipeline,
)


def _drive(workspace, *, threads: int, locus_pipeline: bool, suffix: str) -> bytes:
    """Run the full pipeline on the default (gffutils) backend and return the
    output bytes. On-disk DBs → factory.viable → fused Step-7 path."""
    from lifton import lifton as lifton_main

    out_gff = workspace["out"] / f"fuse_{suffix}.gff3"
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
    args = lifton_main.parse_args(argv)
    lifton_main.run_all_lifton_steps(args)
    return out_gff.read_bytes()


def _pool_spy(monkeypatch):
    """Install a ThreadPoolExecutor spy on lifton.parallel; return the counter."""
    from lifton import parallel as _parallel
    original = _parallel.ThreadPoolExecutor
    constructed = {"n": 0}

    class TPESpy(original):
        def __init__(self, *a, **kw):
            constructed["n"] += 1
            super().__init__(*a, **kw)

    monkeypatch.setattr(_parallel, "ThreadPoolExecutor", TPESpy)
    return constructed


class TestFusedDispatch:
    def test_fused_single_pool_and_byte_identical(
            self, integration_workspace, hermetic_pipeline, monkeypatch):
        """Fused path (default) builds exactly ONE pool (the prefetcher pool
        is gone) and is byte-identical to serial."""
        serial = _drive(integration_workspace, threads=1,
                        locus_pipeline=False, suffix="serial")
        constructed = _pool_spy(monkeypatch)
        fused = _drive(integration_workspace, threads=4,
                       locus_pipeline=True, suffix="fused")
        assert len(serial) > 0
        assert fused == serial
        assert constructed["n"] == 1, (
            "fused Step 7 must create exactly one ThreadPoolExecutor "
            f"(no prefetcher pool); saw {constructed['n']}"
        )

    def test_envgate_restores_two_phase_byte_identical(
            self, integration_workspace, hermetic_pipeline, monkeypatch):
        """LIFTON_FUSE_STEP7=0 restores the two-phase path, byte-identical to
        the fused default and to serial."""
        serial = _drive(integration_workspace, threads=1,
                        locus_pipeline=False, suffix="serial2")
        fused = _drive(integration_workspace, threads=4,
                       locus_pipeline=True, suffix="fused2")
        monkeypatch.setenv("LIFTON_FUSE_STEP7", "0")
        two_phase = _drive(integration_workspace, threads=4,
                           locus_pipeline=True, suffix="twophase")
        assert two_phase == fused == serial

    def test_envgate_restores_two_phase_two_pools(
            self, integration_workspace, hermetic_pipeline, monkeypatch):
        """With LIFTON_FUSE_STEP7=0 the on-disk path builds TWO pools
        (prefetcher + worker)."""
        monkeypatch.setenv("LIFTON_FUSE_STEP7", "0")
        constructed = _pool_spy(monkeypatch)
        out = _drive(integration_workspace, threads=4,
                     locus_pipeline=True, suffix="twophase_pools")
        assert len(out) > 0
        assert constructed["n"] == 2, (
            "two-phase path must create a prefetcher pool + a worker pool; "
            f"saw {constructed['n']}"
        )


class TestFusedMaterialiseFailureIsolation:
    def test_one_materialise_failure_does_not_abort_siblings(
            self, monkeypatch, tmp_path):
        """A raise in the materialise HALF of the fused task becomes a
        per-locus error LocusResult (siblings still emit in order), not a
        pool crash. Uses a fake on-disk-looking DB (dbfn set → factory.viable
        → fused path) with materialise_locus_with_factory spied to raise on
        one locus."""
        from lifton import parallel, run_liftoff, locus_pipeline
        from lifton import parallel as p_mod
        from lifton.locus_pipeline import StepContext, MaterialisedLocus

        # Force the thread-safety guard True (mirrors the existing harness).
        monkeypatch.setattr(p_mod, "_backend_supports_threads",
                            lambda *a, **kw: True)

        loci = [SimpleNamespace(id=f"l_{i}", _idx=i) for i in range(5)]
        dbfile = tmp_path / "fake.db"
        dbfile.write_text("x")  # exists → _extract_dbfn returns a path

        class FakeDB:
            __module__ = "lifton.gffbase.interface"
            dbfn = str(dbfile)      # → factory.viable True → fused path

            def features_of_type(self, ft):
                yield from loci

        # Spy the materialise half: raise on l_2, minimal payload for others.
        def spy_mat(idx, locus, factory):
            if locus.id == "l_2":
                raise ValueError("boom in materialise")
            return MaterialisedLocus(submission_index=idx, locus=locus,
                                     locus_id=locus.id)

        monkeypatch.setattr(locus_pipeline, "materialise_locus_with_factory",
                            spy_mat)

        def fake_process(*args, ENTRY_FEATURE=False, **kw):
            locus = args[1]
            gene = mock.MagicMock()
            gene.ref_gene_id = "ok"
            gene.write_entry.side_effect = (
                lambda fw, stats: fw.write(f"{locus._idx}\n"))
            return gene

        monkeypatch.setattr(run_liftoff, "process_liftoff", fake_process)

        db = FakeDB()
        ctx = StepContext(
            ref_db=db, l_feature_db=db, m_feature_db=None,
            ref_id_2_m_id_trans_dict={}, tree_dict={},
            tgt_fai=mock.Mock(), ref_proteins={}, ref_trans={},
            ref_features_dict={}, fw_score=io.StringIO(), fw_chain=None,
            args=SimpleNamespace(native=False),
        )
        fw = io.StringIO()
        stats = {"coding": {}, "non-coding": {}, "other": {}}
        # Must NOT raise despite the materialise failure on l_2.
        parallel.parallel_step7(["gene"], db, ctx, fw, stats, threads=4)
        # l_2 absent (its materialise raised → error LocusResult, skipped);
        # every sibling emitted in submission order.
        assert fw.getvalue() == "0\n1\n3\n4\n"
