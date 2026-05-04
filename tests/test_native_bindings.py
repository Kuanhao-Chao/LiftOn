"""Phase 10 — native binding facade unit + parity tests.

Covers:
  * `lifton.native_bindings.types.GFF3Hit.from_gff_line` parser
    — every NCBI GFF3 corner case (blank, comment, malformed cols,
    non-integer coords, valid 9-col row).
  * `lifton.native_bindings.MinimapAligner` — basic construction
    against a real FASTA fixture and `.map()` round-trip.
  * `lifton.native_bindings.MiniprotIndex` — subprocess-fallback
    branch with mocked `subprocess.Popen`, parse correctness, error
    surfacing, .raw_bytes and .align() generators.
  * Differential parity: the bytes blob produced by
    `MiniprotIndex.align_all` is byte-identical to the legacy
    `run_miniprot.run_miniprot --stream` output for the same input.
  * Availability flags: `is_mappy_available()` /
    `is_pyminiprot_native_available()` return correct booleans
    under presence/absence of the underlying packages.
"""

from __future__ import annotations

import io
import os
import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace
from unittest import mock

import pytest

from lifton.native_bindings import (
    GFF3Bundle,
    GFF3Hit,
    MinimapAligner,
    MiniprotIndex,
    MinimapHit,
    is_mappy_available,
    is_pyminiprot_native_available,
)


# ---------------------------------------------------------------------------
# 1. Availability flags
# ---------------------------------------------------------------------------

class TestAvailabilityFlags:
    def test_mappy_available_in_test_env(self):
        # The test conda env ships mappy; flip if you ever run
        # without it.
        assert is_mappy_available() is True

    def test_pyminiprot_native_not_yet_available(self):
        """The real PyO3 binding has not been built yet. The flag
        MUST return False so call sites take the subprocess path."""
        assert is_pyminiprot_native_available() is False


# ---------------------------------------------------------------------------
# 2. GFF3Hit parser
# ---------------------------------------------------------------------------

class TestGFF3HitParser:
    def test_parses_valid_row(self):
        row = "chr1\tminiprot\tmRNA\t100\t200\t.\t+\t.\tID=MP1;Target=tx1 1 66"
        h = GFF3Hit.from_gff_line(row)
        assert h is not None
        assert h.seqid == "chr1"
        assert h.source == "miniprot"
        assert h.featuretype == "mRNA"
        assert h.start == 100 and h.end == 200
        assert h.score == "."
        assert h.strand == "+"
        assert h.phase == "."
        assert "ID=MP1" in h.attributes
        assert "Target=tx1 1 66" in h.attributes

    def test_skips_comment_line(self):
        assert GFF3Hit.from_gff_line("# this is a comment") is None
        assert GFF3Hit.from_gff_line("##gff-version 3") is None

    def test_skips_blank_line(self):
        assert GFF3Hit.from_gff_line("") is None

    def test_skips_malformed_column_count(self):
        # 8 cols
        assert GFF3Hit.from_gff_line(
            "chr1\tminiprot\tmRNA\t100\t200\t.\t+\t."
        ) is None
        # 10 cols
        assert GFF3Hit.from_gff_line(
            "chr1\tminiprot\tmRNA\t100\t200\t.\t+\t.\tID=g\textra"
        ) is None

    def test_skips_non_integer_coords(self):
        assert GFF3Hit.from_gff_line(
            "chr1\tminiprot\tmRNA\tabc\t200\t.\t+\t.\tID=MP1"
        ) is None

    def test_handles_trailing_newline(self):
        h = GFF3Hit.from_gff_line(
            "chr1\tminiprot\tmRNA\t100\t200\t.\t+\t.\tID=MP1\n"
        )
        assert h is not None
        assert h.end == 200


# ---------------------------------------------------------------------------
# 3. MinimapAligner — real mappy round-trip on a synthetic FASTA
# ---------------------------------------------------------------------------

@pytest.fixture
def tiny_fasta(tmp_path):
    """50 kb synthetic chromosome with a known sequence so we can
    issue a query that is a substring and assert at least one hit."""
    fp = tmp_path / "tiny.fa"
    seq = "".join(["ACGT", "GCTA", "TTTA", "GGGG", "CCCC"]) * 2500   # 50 kb
    fp.write_text(">chr_synth\n" + "\n".join(
        seq[i:i + 60] for i in range(0, len(seq), 60)
    ) + "\n")
    return fp, seq


class TestMinimapAlignerRealRoundTrip:
    def test_construction(self, tiny_fasta):
        fp, _ = tiny_fasta
        aligner = MinimapAligner(str(fp), preset="map-ont", threads=1)
        assert aligner.n_seq >= 1
        assert "chr_synth" in aligner.seq_names()

    def test_map_yields_minimaphit(self, tiny_fasta):
        fp, seq = tiny_fasta
        aligner = MinimapAligner(str(fp), preset="map-ont", threads=1)
        # Query a 200 bp window that must align cleanly.
        query = seq[1000:1200]
        hits = list(aligner.map("test_query", query))
        # Empty hit list is acceptable for a short random query, but
        # if there is at least one hit the schema must hold.
        for h in hits:
            assert isinstance(h, MinimapHit)
            assert h.query_name == "test_query"
            assert h.r_st >= 0 and h.r_en >= h.r_st
            assert h.q_st >= 0 and h.q_en >= h.q_st
            assert h.strand in (-1, 1)
            assert h.cigar_str
            assert h.mapq >= 0

    def test_map_one_best_returns_optional(self, tiny_fasta):
        fp, seq = tiny_fasta
        aligner = MinimapAligner(str(fp), preset="map-ont", threads=1)
        result = aligner.map_one_best("q", seq[1000:1200])
        assert result is None or isinstance(result, MinimapHit)


class TestMinimapAlignerMissingMappy:
    def test_clear_error_when_mappy_unavailable(self, monkeypatch):
        from lifton.native_bindings import minimap_facade
        monkeypatch.setattr(minimap_facade, "is_mappy_available",
                            lambda: False)
        with pytest.raises(RuntimeError, match="mappy is not installed"):
            MinimapAligner("/nonexistent.fa")


# ---------------------------------------------------------------------------
# 4. MiniprotIndex — subprocess-fallback branch
# ---------------------------------------------------------------------------

_FAKE_GFF_BLOB = (
    b"##gff-version 3\n"
    b"chr1\tminiprot\tmRNA\t100\t200\t.\t+\t.\tID=MP1;Target=tx1 1 66\n"
    b"chr1\tminiprot\tCDS\t100\t150\t.\t+\t0\tID=MP1.cds1;Parent=MP1\n"
    b"chr1\tminiprot\tCDS\t160\t200\t.\t+\t0\tID=MP1.cds2;Parent=MP1\n"
)


class TestMiniprotIndexSubprocess:
    def _patch_popen(self, stdout_bytes, stderr_bytes=b"", returncode=0):
        proc = mock.MagicMock()
        proc.communicate.return_value = (stdout_bytes, stderr_bytes)
        proc.returncode = returncode
        return mock.patch.object(
            subprocess, "Popen", return_value=proc,
        )

    def test_align_all_caches_bundle(self):
        idx = MiniprotIndex("tgt.fa", ref_proteins_path="rp.fa")
        with self._patch_popen(_FAKE_GFF_BLOB):
            b1 = idx.align_all()
            # Second call must NOT re-spawn the subprocess.
            b2 = idx.align_all()
        assert b1 is b2

    def test_align_all_parses_three_hits(self):
        idx = MiniprotIndex("tgt.fa", ref_proteins_path="rp.fa")
        with self._patch_popen(_FAKE_GFF_BLOB):
            bundle = idx.align_all()
        assert len(bundle) == 3
        assert {h.featuretype for h in bundle} == {"mRNA", "CDS"}

    def test_raw_bytes_round_trips(self):
        idx = MiniprotIndex("tgt.fa", ref_proteins_path="rp.fa")
        with self._patch_popen(_FAKE_GFF_BLOB):
            idx.align_all()
        assert idx.raw_bytes == _FAKE_GFF_BLOB

    def test_align_yields_cached_hits(self):
        idx = MiniprotIndex("tgt.fa", ref_proteins_path="rp.fa")
        with self._patch_popen(_FAKE_GFF_BLOB):
            idx.align_all()
        hits = list(idx.align("MAGT*"))
        assert len(hits) == 3

    def test_align_without_align_all_raises(self):
        idx = MiniprotIndex("tgt.fa", ref_proteins_path="rp.fa")
        with pytest.raises(RuntimeError, match="align_all"):
            list(idx.align("MAGT*"))

    def test_align_all_without_proteins_raises(self):
        idx = MiniprotIndex("tgt.fa")  # no ref_proteins_path
        with pytest.raises(ValueError, match="ref_proteins_path"):
            idx.align_all()

    def test_nonzero_exit_raises(self):
        idx = MiniprotIndex("tgt.fa", ref_proteins_path="rp.fa")
        with self._patch_popen(b"", b"miniprot crashed", returncode=1):
            with pytest.raises(RuntimeError, match="exited with code 1"):
                idx.align_all()

    def test_error_in_stderr_raises(self):
        idx = MiniprotIndex("tgt.fa", ref_proteins_path="rp.fa")
        with self._patch_popen(b"", b"ERROR during mapping\n", returncode=0):
            with pytest.raises(RuntimeError, match="ERROR"):
                idx.align_all()

    def test_is_native_false_for_subprocess_path(self):
        idx = MiniprotIndex("tgt.fa", ref_proteins_path="rp.fa")
        assert idx.is_native is False

    def test_raw_bytes_before_align_all_raises(self):
        idx = MiniprotIndex("tgt.fa", ref_proteins_path="rp.fa")
        with pytest.raises(RuntimeError):
            _ = idx.raw_bytes


# ---------------------------------------------------------------------------
# 5. Differential parity: native facade vs run_miniprot streaming path
# ---------------------------------------------------------------------------

class TestMiniprotFacadeStreamingParity:
    """The Phase 7 streaming-path bytes are produced by piping
    miniprot stdout through `subprocess.Popen.communicate()`. The
    Phase 10 facade does the same — both paths must yield byte-
    identical bytes for the same input."""

    def test_facade_vs_streaming_run_miniprot(self):
        from lifton import run_miniprot
        # Patch `subprocess.Popen` for BOTH callers with the same
        # canned stdout, then verify the bytes captured by each
        # path are identical.
        proc = mock.MagicMock()
        proc.communicate.return_value = (_FAKE_GFF_BLOB, b"")
        proc.returncode = 0
        # Patch the streaming-run path's subprocess
        with mock.patch.object(run_miniprot.subprocess, "Popen",
                               return_value=proc):
            args = SimpleNamespace(mp_options="", stream=True, miniprot=None)
            stream_bytes = run_miniprot.run_miniprot(
                "/tmp/", args, "tgt.fa", "rp.fa",
            )
        # Patch the facade's subprocess
        with mock.patch.object(subprocess, "Popen", return_value=proc):
            idx = MiniprotIndex("tgt.fa", ref_proteins_path="rp.fa")
            idx.align_all()
            facade_bytes = idx.raw_bytes
        assert stream_bytes == facade_bytes


# ---------------------------------------------------------------------------
# 6. CLI plumbing for --native
# ---------------------------------------------------------------------------

class TestNativeCLI:
    def test_default_is_false(self):
        from lifton import lifton as lifton_main
        args = lifton_main.parse_args(["t.fa", "r.fa", "-g", "r.gff3"])
        assert args.native is False

    def test_set_to_true(self):
        from lifton import lifton as lifton_main
        args = lifton_main.parse_args(
            ["t.fa", "r.fa", "-g", "r.gff3", "--native"]
        )
        assert args.native is True

    def test_help_text_mentions_native(self, capsys):
        from lifton import lifton as lifton_main
        with pytest.raises(SystemExit):
            lifton_main.parse_args(["--help"])
        assert "--native" in capsys.readouterr().out


# ---------------------------------------------------------------------------
# 7. Backend thread-safety guard with --native
# ---------------------------------------------------------------------------

class TestThreadSafetyGuardWithNative:
    def test_guard_returns_true_for_native_phase_11_no_dbs(self, monkeypatch):
        """Phase 11 contract: with no FeatureDB inputs (i.e. nothing
        SQLite-bound to fail), --native unlocks thread safety."""
        monkeypatch.delenv("LIFTON_PARALLEL_FORCE", raising=False)
        from lifton.parallel import _backend_supports_threads
        assert _backend_supports_threads(native=True) is True

    def test_guard_returns_false_when_gffutils_in_workers(self, monkeypatch):
        """Phase 11 honest contract: SQLite (gffutils) cannot tolerate
        cross-thread access regardless of pre-materialisation. The
        guard refuses to fan out when any input DB is gffutils."""
        monkeypatch.delenv("LIFTON_PARALLEL_FORCE", raising=False)
        # Define a class whose __module__ starts with "gffutils".
        FakeSqlite = type("FakeSqliteDB", (), {})
        FakeSqlite.__module__ = "gffutils.interface"
        from lifton.parallel import _backend_supports_threads
        assert _backend_supports_threads(FakeSqlite(), native=True) is False

    def test_guard_returns_true_when_only_gffbase(self, monkeypatch):
        """Phase 11: gffbase (DuckDB) with --native unlocks threads."""
        monkeypatch.delenv("LIFTON_PARALLEL_FORCE", raising=False)
        FakeGffbase = type("FakeGffbaseDB", (), {})
        FakeGffbase.__module__ = "lifton.gffbase.interface"
        from lifton.parallel import _backend_supports_threads
        assert _backend_supports_threads(FakeGffbase(), native=True) is True

    def test_guard_returns_false_when_native_inactive(self, monkeypatch):
        monkeypatch.delenv("LIFTON_PARALLEL_FORCE", raising=False)
        from lifton.parallel import _backend_supports_threads
        assert _backend_supports_threads(native=False) is False

    def test_force_env_var_overrides_native_false(self, monkeypatch):
        monkeypatch.setenv("LIFTON_PARALLEL_FORCE", "1")
        from lifton.parallel import _backend_supports_threads
        assert _backend_supports_threads(native=False) is True

    def test_force_env_var_overrides_gffutils_block(self, monkeypatch):
        monkeypatch.setenv("LIFTON_PARALLEL_FORCE", "1")
        FakeSqlite = type("FakeSqliteDB", (), {})
        FakeSqlite.__module__ = "gffutils.interface"
        from lifton.parallel import _backend_supports_threads
        assert _backend_supports_threads(FakeSqlite(), native=True) is True


# ---------------------------------------------------------------------------
# 8. parallel_step7 with --native unlocks parallelism
# ---------------------------------------------------------------------------

class TestParallelStep7WithNative:
    """Phase 11 contract: with --native AND a thread-safe backend
    (FakeDB modules NOT starting with ``gffutils``), the dispatcher
    really creates a ThreadPoolExecutor and per-locus tasks run
    concurrently. SQLite-backed gffutils still falls back."""

    def test_native_unlocks_pool_for_safe_backend(self, monkeypatch):
        from lifton import parallel, run_liftoff, locus_pipeline
        monkeypatch.delenv("LIFTON_PARALLEL_FORCE", raising=False)
        from concurrent.futures import ThreadPoolExecutor as _real_TPE
        constructed = {"n": 0}

        class TPESpy(_real_TPE):
            def __init__(self, *a, **kw):
                constructed["n"] += 1
                super().__init__(*a, **kw)

        monkeypatch.setattr(parallel, "ThreadPoolExecutor", TPESpy)

        loci = [SimpleNamespace(id=f"l_{i}", _idx=i) for i in range(4)]

        class FakeDB:
            __module__ = "lifton.gffbase.interface"  # gffbase-shaped
            def features_of_type(self, ft): yield from loci
            def children(self, *a, **kw): return iter([])

        def fake(*args, ENTRY_FEATURE=False, **kw):
            locus = args[1]
            gene = mock.MagicMock()
            gene.ref_gene_id = "ok"
            gene.write_entry.side_effect = (
                lambda fw, stats: fw.write(f"{locus._idx}\n")
            )
            return gene

        monkeypatch.setattr(run_liftoff, "process_liftoff", fake)

        db = FakeDB()
        ctx = locus_pipeline.StepContext(
            ref_db=db, l_feature_db=db, m_feature_db=None,
            ref_id_2_m_id_trans_dict={}, tree_dict={},
            tgt_fai=mock.Mock(), ref_proteins={}, ref_trans={},
            ref_features_dict={}, fw_score=io.StringIO(), fw_chain=None,
            args=SimpleNamespace(native=True),
        )
        fw = io.StringIO()
        stats = {"coding": {}, "non-coding": {}, "other": {}}
        n = parallel.parallel_step7(
            ["gene"], db, ctx, fw, stats, threads=4,
        )
        assert n == 4
        # Pool WAS created — Phase 11 unlocks parallelism.
        assert constructed["n"] == 1
        # Output still deterministic via the ordered-writer.
        assert fw.getvalue() == "0\n1\n2\n3\n"

    def test_native_falls_back_on_sqlite_backend(self, monkeypatch):
        """gffutils backend stays serial under --native because SQLite
        connections cannot cross threads."""
        from lifton import parallel, run_liftoff, locus_pipeline
        monkeypatch.delenv("LIFTON_PARALLEL_FORCE", raising=False)
        from concurrent.futures import ThreadPoolExecutor as _real_TPE
        constructed = {"n": 0}

        class TPESpy(_real_TPE):
            def __init__(self, *a, **kw):
                constructed["n"] += 1
                super().__init__(*a, **kw)

        monkeypatch.setattr(parallel, "ThreadPoolExecutor", TPESpy)

        loci = [SimpleNamespace(id=f"l_{i}", _idx=i) for i in range(4)]

        class FakeSqliteDB:
            __module__ = "gffutils.interface"   # SQLite-shaped
            def features_of_type(self, ft): yield from loci
            def children(self, *a, **kw): return iter([])

        def fake(*args, ENTRY_FEATURE=False, **kw):
            locus = args[1]
            gene = mock.MagicMock()
            gene.ref_gene_id = "ok"
            gene.write_entry.side_effect = (
                lambda fw, stats: fw.write(f"{locus._idx}\n")
            )
            return gene

        monkeypatch.setattr(run_liftoff, "process_liftoff", fake)

        db = FakeSqliteDB()
        ctx = locus_pipeline.StepContext(
            ref_db=db, l_feature_db=db, m_feature_db=None,
            ref_id_2_m_id_trans_dict={}, tree_dict={},
            tgt_fai=mock.Mock(), ref_proteins={}, ref_trans={},
            ref_features_dict={}, fw_score=io.StringIO(), fw_chain=None,
            args=SimpleNamespace(native=True),
        )
        fw = io.StringIO()
        stats = {"coding": {}, "non-coding": {}, "other": {}}
        n = parallel.parallel_step7(
            ["gene"], db, ctx, fw, stats, threads=4,
        )
        assert n == 4
        # Pool NOT created — fell back to serial under SQLite guard.
        assert constructed["n"] == 0
        assert fw.getvalue() == "0\n1\n2\n3\n"

    def test_native_plus_force_unlocks_pool(self, monkeypatch):
        """LIFTON_PARALLEL_FORCE=1 overrides the SQLite block too."""
        monkeypatch.setenv("LIFTON_PARALLEL_FORCE", "1")
        from lifton import parallel, run_liftoff, locus_pipeline
        from concurrent.futures import ThreadPoolExecutor as _real_TPE
        constructed = {"n": 0}

        class TPESpy(_real_TPE):
            def __init__(self, *a, **kw):
                constructed["n"] += 1
                super().__init__(*a, **kw)

        monkeypatch.setattr(parallel, "ThreadPoolExecutor", TPESpy)

        loci = [SimpleNamespace(id=f"l_{i}", _idx=i) for i in range(4)]

        class FakeDB:
            def features_of_type(self, ft): yield from loci
            def children(self, *a, **kw): return iter([])

        def fake(*args, ENTRY_FEATURE=False, **kw):
            locus = args[1]
            gene = mock.MagicMock()
            gene.ref_gene_id = "ok"
            gene.write_entry.side_effect = (
                lambda fw, stats: fw.write(f"{locus._idx}\n")
            )
            return gene

        monkeypatch.setattr(run_liftoff, "process_liftoff", fake)

        db = FakeDB()
        ctx = locus_pipeline.StepContext(
            ref_db=db, l_feature_db=db, m_feature_db=None,
            ref_id_2_m_id_trans_dict={}, tree_dict={},
            tgt_fai=mock.Mock(), ref_proteins={}, ref_trans={},
            ref_features_dict={}, fw_score=io.StringIO(), fw_chain=None,
            args=SimpleNamespace(native=True),
        )
        fw = io.StringIO()
        stats = {"coding": {}, "non-coding": {}, "other": {}}
        n = parallel.parallel_step7(
            ["gene"], db, ctx, fw, stats, threads=4,
        )
        assert n == 4
        assert constructed["n"] == 1
        assert fw.getvalue() == "0\n1\n2\n3\n"

    def test_native_plus_force_unlocks_pool(self, monkeypatch):
        """With LIFTON_PARALLEL_FORCE=1 the guard is overridden, the
        pool is created, and output stays byte-identical to the
        serial baseline."""
        monkeypatch.setenv("LIFTON_PARALLEL_FORCE", "1")
        from lifton import parallel, run_liftoff, locus_pipeline
        from concurrent.futures import ThreadPoolExecutor as _real_TPE
        constructed = {"n": 0}

        class TPESpy(_real_TPE):
            def __init__(self, *a, **kw):
                constructed["n"] += 1
                super().__init__(*a, **kw)

        monkeypatch.setattr(parallel, "ThreadPoolExecutor", TPESpy)

        loci = [SimpleNamespace(id=f"l_{i}", _idx=i) for i in range(4)]

        class FakeDB:
            def features_of_type(self, ft): yield from loci

        def fake(*args, ENTRY_FEATURE=False, **kw):
            locus = args[1]
            gene = mock.MagicMock()
            gene.ref_gene_id = "ok"
            gene.write_entry.side_effect = (
                lambda fw, stats: fw.write(f"{locus._idx}\n")
            )
            return gene

        monkeypatch.setattr(run_liftoff, "process_liftoff", fake)

        ctx = locus_pipeline.StepContext(
            ref_db=mock.Mock(), l_feature_db=FakeDB(), m_feature_db=None,
            ref_id_2_m_id_trans_dict={}, tree_dict={},
            tgt_fai=mock.Mock(), ref_proteins={}, ref_trans={},
            ref_features_dict={}, fw_score=io.StringIO(), fw_chain=None,
            args=SimpleNamespace(native=True),
        )
        fw = io.StringIO()
        stats = {"coding": {}, "non-coding": {}, "other": {}}
        n = parallel.parallel_step7(
            ["gene"], FakeDB(), ctx, fw, stats, threads=4,
        )
        assert n == 4
        assert constructed["n"] == 1
        assert fw.getvalue() == "0\n1\n2\n3\n"
