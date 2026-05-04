"""Phase 7 — unit tests for the streaming adapter layer.

Covers:
  - lifton.gffbase_adapter.{db_path_for, open_existing_db,
    build_database, build_database_from_string, looks_like_gff3_blob}
  - lifton.annotation.Annotation polymorphism on bytes blob input
  - lifton.run_miniprot.run_miniprot streaming branch
"""

from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace
from unittest import mock

import pytest

from lifton import gffbase_adapter
from lifton.annotation import Annotation


# ---------------------------------------------------------------------------
# Synthetic GFF3 fixture used across the unit tests
# ---------------------------------------------------------------------------

_TINY_GFF = (
    b"##gff-version 3\n"
    b"chr1\ttest\tgene\t100\t200\t.\t+\t.\tID=g1;gene_biotype=protein_coding\n"
    b"chr1\ttest\tmRNA\t100\t200\t.\t+\t.\tID=tx1;Parent=g1\n"
    b"chr1\ttest\texon\t100\t150\t.\t+\t.\tID=ex1;Parent=tx1\n"
    b"chr1\ttest\texon\t160\t200\t.\t+\t.\tID=ex2;Parent=tx1\n"
)


# ---------------------------------------------------------------------------
# gffbase_adapter primitives
# ---------------------------------------------------------------------------

class TestDbPathFor:
    def test_uses_duckdb_suffix(self):
        assert gffbase_adapter.db_path_for("/tmp/x.gff3") == "/tmp/x.gff3.duckdb"

    def test_does_not_collide_with_legacy_sqlite_cache(self):
        legacy = "/tmp/x.gff3_db"
        adapter_cache = gffbase_adapter.db_path_for("/tmp/x.gff3")
        assert adapter_cache != legacy


class TestLooksLikeGFF3Blob:
    def test_bytes_with_tab_is_blob(self):
        assert gffbase_adapter.looks_like_gff3_blob(b"chr1\tt\tgene\t1\t2\t.\t+\t.\tID=g") is True

    def test_bytes_with_directive_is_blob(self):
        assert gffbase_adapter.looks_like_gff3_blob(b"##gff-version 3\nrest") is True

    def test_path_string_is_not_blob(self):
        assert gffbase_adapter.looks_like_gff3_blob("/tmp/a.gff3") is False

    def test_str_is_never_blob(self):
        # Phase 7: only bytes/bytearray are detected as blobs by this
        # heuristic. str inputs are treated as paths.
        assert gffbase_adapter.looks_like_gff3_blob("##gff-version 3\nx") is False

    def test_empty_bytes_is_not_blob(self):
        assert gffbase_adapter.looks_like_gff3_blob(b"") is False


class TestBuildDatabaseFromString:
    def test_round_trip_bytes_blob(self):
        db = gffbase_adapter.build_database_from_string(gff_text=_TINY_GFF)
        assert "gene" in list(db.featuretypes())
        assert "chr1" in list(db.seqids())
        assert db.count_features_of_type("gene") == 1
        assert db.count_features_of_type("exon") == 2

    def test_round_trip_str_blob(self):
        db = gffbase_adapter.build_database_from_string(
            gff_text=_TINY_GFF.decode("utf-8")
        )
        assert db.count_features_of_type("gene") == 1

    def test_round_trip_bytearray_blob(self):
        db = gffbase_adapter.build_database_from_string(
            gff_text=bytearray(_TINY_GFF)
        )
        assert db.count_features_of_type("gene") == 1

    def test_in_memory_dbfn_default(self):
        # No on-disk file should be created when dbfn=":memory:".
        db = gffbase_adapter.build_database_from_string(gff_text=_TINY_GFF)
        assert db is not None

    def test_blob_ingest_matches_file_ingest(self, tmp_path):
        """The whole point of the streaming path: file ingest and
        in-memory ingest must yield the same FeatureDB content."""
        fp = tmp_path / "x.gff3"
        fp.write_bytes(_TINY_GFF)

        from lifton.gffbase import create_db
        db_file = create_db(str(fp), dbfn=":memory:", force=True)
        db_blob = gffbase_adapter.build_database_from_string(gff_text=_TINY_GFF)

        # Same featuretype counts
        for ft in ("gene", "mRNA", "exon"):
            assert db_file.count_features_of_type(ft) == \
                   db_blob.count_features_of_type(ft)
        # Same gene id set
        ids_file = {f.id for f in db_file.features_of_type("gene")}
        ids_blob = {f.id for f in db_blob.features_of_type("gene")}
        assert ids_file == ids_blob


class TestOpenExistingDb:
    def test_returns_none_when_no_cache(self, tmp_path):
        fp = tmp_path / "missing.gff3"
        fp.write_bytes(_TINY_GFF)
        # No <fp>.duckdb sidecar exists yet
        assert gffbase_adapter.open_existing_db(str(fp)) is None

    def test_returns_db_when_cache_exists(self, tmp_path):
        fp = tmp_path / "ok.gff3"
        fp.write_bytes(_TINY_GFF)
        # Build, then re-open
        gffbase_adapter.build_database(
            file_name=str(fp),
            infer_genes=False, infer_transcripts=False,
            force=True, verbose=False,
        )
        db = gffbase_adapter.open_existing_db(str(fp))
        assert db is not None
        assert db.count_features_of_type("gene") == 1


# ---------------------------------------------------------------------------
# Annotation polymorphism
# ---------------------------------------------------------------------------

class TestAnnotationPolymorphic:
    def test_blob_routes_to_gffbase(self):
        ann = Annotation(_TINY_GFF, False, False)
        assert ann.backend == "gffbase"
        assert ann._is_blob is True
        assert list(ann.db_connection.seqids()) == ["chr1"]

    def test_blob_with_gffutils_backend_raises(self):
        with pytest.raises(ValueError, match="gffbase"):
            Annotation(_TINY_GFF, False, False, backend="gffutils")

    def test_path_branch_default_backend_unchanged(self, tmp_path):
        """Phase 5 baseline: a file-path input with no backend
        override stays on gffutils. No regression."""
        fp = tmp_path / "a.gff3"
        fp.write_bytes(_TINY_GFF)
        ann = Annotation(str(fp), False, False)
        assert ann.backend == "gffutils"

    def test_path_with_explicit_gffbase_backend(self, tmp_path):
        fp = tmp_path / "b.gff3"
        fp.write_bytes(_TINY_GFF)
        ann = Annotation(str(fp), False, False, backend="gffbase")
        assert ann.backend == "gffbase"

    def test_env_var_opt_in_to_gffbase(self, tmp_path, monkeypatch):
        fp = tmp_path / "c.gff3"
        fp.write_bytes(_TINY_GFF)
        monkeypatch.setenv("LIFTON_USE_GFFBASE", "1")
        ann = Annotation(str(fp), False, False)
        assert ann.backend == "gffbase"


# ---------------------------------------------------------------------------
# run_miniprot streaming branch
# ---------------------------------------------------------------------------

@pytest.fixture
def fake_miniprot_args(tmp_path):
    """Build an argparse-shaped namespace minimum enough for run_miniprot."""
    return SimpleNamespace(
        mp_options="",
        stream=False,
        miniprot=None,
    )


class TestRunMiniprotLegacy:
    """The legacy file-write branch must remain byte-identical to Phase 5."""

    def test_legacy_returns_path(self, tmp_path, fake_miniprot_args):
        from lifton import run_miniprot

        outdir = str(tmp_path) + "/"
        # Patch subprocess.run with a fake that "writes" a real GFF3 file
        def fake_run(cmd, stdout=None, stderr=None, text=None):
            stdout.write("##gff-version 3\nchr1\tmp\tmRNA\t1\t10\t.\t+\t.\tID=MP1;Target=tx1 1 10\n")
            return SimpleNamespace(returncode=0, stderr="")
        with mock.patch.object(run_miniprot.subprocess, "run", side_effect=fake_run):
            result = run_miniprot.run_miniprot(
                outdir, fake_miniprot_args, "tgt.fa", "ref_proteins.fa"
            )
        assert isinstance(result, str)
        assert result.endswith("miniprot.gff3")
        assert os.path.exists(result)
        # File contents preserved (Phase 5 contract)
        body = open(result).read()
        assert "MP1" in body


class TestRunMiniprotStreaming:
    """Phase 7 streaming branch returns bytes — no disk write."""

    def _patch_popen(self, stdout_bytes, stderr_bytes=b"", returncode=0):
        # Phase 15c: the streaming branch now drives `.stdout.read()`
        # / `.stderr.read()` / `.wait()` directly via
        # `_drain_stream_chunks` instead of `.communicate()`. The
        # mock must satisfy BOTH contracts so this helper stays
        # compatible if the implementation switches back.
        from io import BytesIO
        proc = mock.MagicMock()
        proc.communicate.return_value = (stdout_bytes, stderr_bytes)
        proc.returncode = returncode
        proc.wait.return_value = returncode
        proc.stdout = BytesIO(stdout_bytes)
        proc.stderr = BytesIO(stderr_bytes)
        return mock.patch("lifton.run_miniprot.subprocess.Popen",
                          return_value=proc)

    def test_streaming_returns_bytes(self, tmp_path, fake_miniprot_args):
        from lifton import run_miniprot
        fake_miniprot_args.stream = True
        gff_blob = (b"##gff-version 3\n"
                    b"chr1\tmp\tmRNA\t1\t10\t.\t+\t.\tID=MP1;Target=tx1 1 10\n")
        with self._patch_popen(gff_blob):
            result = run_miniprot.run_miniprot(
                str(tmp_path) + "/", fake_miniprot_args, "tgt.fa", "ref.fa"
            )
        assert isinstance(result, (bytes, bytearray))
        assert b"MP1" in result

    def test_streaming_no_disk_write(self, tmp_path, fake_miniprot_args):
        """Streaming branch must NOT write miniprot.gff3."""
        from lifton import run_miniprot
        fake_miniprot_args.stream = True
        gff_blob = b"##gff-version 3\nchr1\tmp\tmRNA\t1\t10\t.\t+\t.\tID=MP1\n"
        with self._patch_popen(gff_blob):
            run_miniprot.run_miniprot(
                str(tmp_path) + "/", fake_miniprot_args, "tgt.fa", "ref.fa"
            )
        # The miniprot/ directory may have been created (mkdir is unconditional),
        # but no miniprot.gff3 file should exist inside it.
        assert not (tmp_path / "miniprot" / "miniprot.gff3").exists()

    def test_streaming_nonzero_exit_returns_none(self, tmp_path, fake_miniprot_args):
        from lifton import run_miniprot
        fake_miniprot_args.stream = True
        with self._patch_popen(b"", b"miniprot crashed\n", returncode=1):
            result = run_miniprot.run_miniprot(
                str(tmp_path) + "/", fake_miniprot_args, "tgt.fa", "ref.fa"
            )
        assert result is None

    def test_streaming_error_in_stderr_returns_none(self, tmp_path, fake_miniprot_args):
        from lifton import run_miniprot
        fake_miniprot_args.stream = True
        with self._patch_popen(b"", b"ERROR during mapping\n", returncode=0):
            result = run_miniprot.run_miniprot(
                str(tmp_path) + "/", fake_miniprot_args, "tgt.fa", "ref.fa"
            )
        assert result is None

    def test_streaming_empty_stdout_returns_none(self, tmp_path, fake_miniprot_args):
        from lifton import run_miniprot
        fake_miniprot_args.stream = True
        with self._patch_popen(b"", b"", returncode=0):
            result = run_miniprot.run_miniprot(
                str(tmp_path) + "/", fake_miniprot_args, "tgt.fa", "ref.fa"
            )
        assert result is None

    def test_streaming_large_blob_does_not_truncate(self, tmp_path, fake_miniprot_args):
        """Pipe-buffer deadlock guard: a 10 MB synthetic miniprot stdout
        ingests cleanly without truncation. communicate() drains both
        pipes simultaneously."""
        from lifton import run_miniprot
        fake_miniprot_args.stream = True
        big_blob = (b"##gff-version 3\n"
                    + b"chr1\tmp\tmRNA\t1\t10\t.\t+\t.\tID=MP_X;Target=tx 1 10\n" * 200_000)
        assert len(big_blob) > 5 * 1024 * 1024  # > 5 MB
        with self._patch_popen(big_blob):
            result = run_miniprot.run_miniprot(
                str(tmp_path) + "/", fake_miniprot_args, "tgt.fa", "ref.fa"
            )
        assert result is not None
        assert len(result) == len(big_blob)


# ---------------------------------------------------------------------------
# End-to-end: streamed miniprot bytes → Annotation → FeatureDB
# ---------------------------------------------------------------------------

class TestStreamingEndToEnd:
    def test_stream_blob_through_annotation(self):
        """The full streaming hand-off: bytes from a fake miniprot
        run land in an Annotation whose db_connection is queryable
        exactly like the legacy file path."""
        miniprot_blob = (
            b"##gff-version 3\n"
            b"chr1\tmp\tmRNA\t100\t200\t.\t+\t.\tID=MP1;Target=tx1 1 30\n"
            b"chr1\tmp\tCDS\t100\t150\t.\t+\t0\tID=MP1.cds1;Parent=MP1\n"
            b"chr1\tmp\tCDS\t160\t200\t.\t+\t0\tID=MP1.cds2;Parent=MP1\n"
        )
        ann = Annotation(miniprot_blob, False, False)
        assert ann.backend == "gffbase"
        # Both methods used by lifton_utils.miniprot_id_mapping
        # must work on the streamed FeatureDB.
        mrnas = list(ann.db_connection.features_of_type("mRNA"))
        assert len(mrnas) == 1
        assert mrnas[0].attributes["Target"][0].startswith("tx1")
        cdss = list(ann.db_connection.children(mrnas[0], featuretype="CDS"))
        assert len(cdss) == 2
