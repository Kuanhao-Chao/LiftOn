"""Regression tests for content-aware annotation database caches."""

from __future__ import annotations

import json
from pathlib import Path

import gffutils
import pytest

from lifton import __version__, annotation, annotation_cache, gffbase_adapter


def _close(db) -> None:
    connection = getattr(db, "conn", None)
    if connection is not None:
        connection.close()


def _manifest(db_path: str) -> dict:
    return json.loads(Path(annotation_cache.manifest_path_for(db_path)).read_text())


class TestGffutilsCacheManifest:
    def test_manifest_records_source_versions_and_settings(self, gff_standard):
        ann = annotation.Annotation(str(gff_standard), False, False, force=True)
        db_path = str(gff_standard) + "_db"
        manifest = _manifest(db_path)

        assert manifest["source"]["size"] == gff_standard.stat().st_size
        assert len(manifest["source"]["sha256"]) == 64
        assert manifest["backend"]["name"] == "gffutils"
        assert manifest["backend"]["version"] == gffutils.__version__
        assert len(manifest["backend"]["schema_version"]) == 64
        assert manifest["tool"] == {"name": "LiftOn", "version": __version__}
        assert manifest["settings"]["inference"] == {
            "genes": False,
            "transcripts": False,
        }
        assert manifest["settings"]["merge_strategy"] == "create_unique"
        assert manifest["settings"]["id_spec"] is None
        _close(ann.db_connection)

    def test_matching_manifest_reuses_database(self, gff_standard, monkeypatch):
        first = annotation.Annotation(str(gff_standard), False, False, force=True)
        _close(first.db_connection)

        def unexpected_build(*args, **kwargs):
            raise AssertionError("matching cache should not be rebuilt")

        monkeypatch.setattr(annotation.gffutils, "create_db", unexpected_build)
        second = annotation.Annotation(str(gff_standard), False, False, force=False)
        assert second.db_connection["gene1"].id == "gene1"
        _close(second.db_connection)

    def test_same_size_source_change_invalidates_by_sha256(self, gff_standard, monkeypatch):
        first = annotation.Annotation(str(gff_standard), False, False, force=True)
        _close(first.db_connection)
        original_manifest = _manifest(str(gff_standard) + "_db")

        changed = gff_standard.read_text().replace("gene1", "gene9")
        assert len(changed.encode()) == gff_standard.stat().st_size
        gff_standard.write_text(changed)

        real_create_db = annotation.gffutils.create_db
        calls = []

        def recording_build(*args, **kwargs):
            calls.append(args[1])
            return real_create_db(*args, **kwargs)

        monkeypatch.setattr(annotation.gffutils, "create_db", recording_build)
        second = annotation.Annotation(str(gff_standard), False, False, force=False)
        assert calls
        assert second.db_connection["gene9"].id == "gene9"
        assert _manifest(str(gff_standard) + "_db")["source"]["sha256"] != \
            original_manifest["source"]["sha256"]
        _close(second.db_connection)

    def test_parser_settings_each_invalidate_cache(self, gff_standard, monkeypatch):
        first = annotation.Annotation(str(gff_standard), False, False, force=True)
        _close(first.db_connection)

        real_create_db = annotation.gffutils.create_db
        calls = []

        def recording_build(*args, **kwargs):
            calls.append(kwargs)
            return real_create_db(*args, **kwargs)

        monkeypatch.setattr(annotation.gffutils, "create_db", recording_build)

        changed_merge = annotation.Annotation(
            str(gff_standard), False, False, merge_strategy="merge", force=False,
        )
        _close(changed_merge.db_connection)
        assert len(calls) == 1

        changed_id = annotation.Annotation(
            str(gff_standard), False, False, merge_strategy="merge", id_spec="ID", force=False,
        )
        _close(changed_id.db_connection)
        assert len(calls) == 2

        changed_inference = annotation.Annotation(
            str(gff_standard), False, True, merge_strategy="merge", id_spec="ID", force=False,
        )
        _close(changed_inference.db_connection)
        assert len(calls) == 3

    def test_force_rebuilds_matching_cache(self, gff_standard, monkeypatch):
        first = annotation.Annotation(str(gff_standard), False, False, force=True)
        _close(first.db_connection)
        real_create_db = annotation.gffutils.create_db
        calls = []

        def recording_build(*args, **kwargs):
            calls.append(args[1])
            return real_create_db(*args, **kwargs)

        monkeypatch.setattr(annotation.gffutils, "create_db", recording_build)
        second = annotation.Annotation(str(gff_standard), False, False, force=True)
        assert len(calls) == 1
        _close(second.db_connection)

    @pytest.mark.parametrize("manifest_change", ["missing", "corrupt", "tool-version"])
    def test_untrusted_manifest_is_never_reused(
        self, gff_standard, monkeypatch, manifest_change,
    ):
        first = annotation.Annotation(str(gff_standard), False, False, force=True)
        _close(first.db_connection)
        db_path = str(gff_standard) + "_db"
        manifest_path = Path(annotation_cache.manifest_path_for(db_path))
        if manifest_change == "missing":
            manifest_path.unlink()
        elif manifest_change == "corrupt":
            manifest_path.write_text("not-json")
        else:
            value = json.loads(manifest_path.read_text())
            value["tool"]["version"] = "stale"
            manifest_path.write_text(json.dumps(value))

        real_create_db = annotation.gffutils.create_db
        calls = []

        def recording_build(*args, **kwargs):
            calls.append(args[1])
            return real_create_db(*args, **kwargs)

        monkeypatch.setattr(annotation.gffutils, "create_db", recording_build)
        rebuilt = annotation.Annotation(str(gff_standard), False, False, force=False)
        assert len(calls) == 1
        _close(rebuilt.db_connection)

    def test_failed_rebuild_keeps_published_database(self, gff_standard, monkeypatch):
        first = annotation.Annotation(str(gff_standard), False, False, force=True)
        _close(first.db_connection)
        db_path = Path(str(gff_standard) + "_db")
        original_db = db_path.read_bytes()
        original_manifest = Path(annotation_cache.manifest_path_for(str(db_path))).read_bytes()
        gff_standard.write_text(gff_standard.read_text().replace("gene1", "gene9"))

        def fail_build(*args, **kwargs):
            raise RuntimeError("injected build failure")

        monkeypatch.setattr(annotation.gffutils, "create_db", fail_build)
        with pytest.raises(SystemExit):
            annotation.Annotation(str(gff_standard), False, False, force=False)

        assert db_path.read_bytes() == original_db
        assert Path(annotation_cache.manifest_path_for(str(db_path))).read_bytes() == original_manifest
        assert not list(db_path.parent.glob(f".{db_path.name}.tmp.*"))


class TestGffbaseCacheManifest:
    def test_content_settings_and_force_control_reuse(self, gff_standard, monkeypatch):
        monkeypatch.setenv("LIFTON_DISABLE_RTREE", "1")
        first = annotation.Annotation(
            str(gff_standard), False, False, force=True, backend="gffbase",
        )
        _close(first.db_connection)
        db_path = str(gff_standard) + ".duckdb"
        manifest = _manifest(db_path)
        assert manifest["backend"]["name"] == "gffbase"
        assert manifest["source"]["size"] == gff_standard.stat().st_size

        real_create_db = gffbase_adapter._gffbase.create_db
        calls = []

        def recording_build(*args, **kwargs):
            calls.append(kwargs["dbfn"])
            return real_create_db(*args, **kwargs)

        monkeypatch.setattr(gffbase_adapter._gffbase, "create_db", recording_build)
        reused = annotation.Annotation(
            str(gff_standard), False, False, force=False, backend="gffbase",
        )
        assert not calls
        _close(reused.db_connection)

        changed = gff_standard.read_text().replace("gene1", "gene9")
        assert len(changed.encode()) == manifest["source"]["size"]
        gff_standard.write_text(changed)
        invalidated = annotation.Annotation(
            str(gff_standard), False, False, force=False, backend="gffbase",
        )
        assert len(calls) == 1
        assert invalidated.db_connection["gene9"].id == "gene9"
        _close(invalidated.db_connection)

        forced = annotation.Annotation(
            str(gff_standard), False, False, force=True, backend="gffbase",
        )
        assert len(calls) == 2
        _close(forced.db_connection)


class TestDirectDatabaseInput:
    def test_gffutils_database_bypasses_text_validation(self, gff_standard, tmp_path, monkeypatch):
        direct_path = tmp_path / "direct.sqlite"
        direct = gffutils.create_db(str(gff_standard), str(direct_path), force=True)
        _close(direct)

        def unexpected_validation(*args, **kwargs):
            raise AssertionError("direct database must not be text-validated")

        monkeypatch.setattr(annotation, "validate_annotation_file", unexpected_validation)
        ann = annotation.Annotation(str(direct_path), False, False)
        assert ann.backend == "gffutils"
        assert ann.db_connection["gene1"].id == "gene1"
        _close(ann.db_connection)

    def test_gffbase_database_bypasses_text_validation(self, gff_standard, monkeypatch):
        monkeypatch.setenv("LIFTON_DISABLE_RTREE", "1")
        built = gffbase_adapter.build_database(
            file_name=str(gff_standard),
            infer_genes=False,
            infer_transcripts=False,
            force=True,
        )
        direct_path = gffbase_adapter.db_path_for(str(gff_standard))
        _close(built)

        def unexpected_validation(*args, **kwargs):
            raise AssertionError("direct database must not be text-validated")

        monkeypatch.setattr(annotation, "validate_annotation_file", unexpected_validation)
        ann = annotation.Annotation(direct_path, False, False, backend="gffbase")
        assert ann.backend == "gffbase"
        assert ann.db_connection["gene1"].id == "gene1"
        _close(ann.db_connection)

    def test_renamed_gffbase_database_is_detected_by_signature(
            self, gff_standard, tmp_path, monkeypatch):
        monkeypatch.setenv("LIFTON_DISABLE_RTREE", "1")
        built = gffbase_adapter.build_database(
            file_name=str(gff_standard),
            infer_genes=False,
            infer_transcripts=False,
            force=True,
        )
        original_path = Path(gffbase_adapter.db_path_for(str(gff_standard)))
        _close(built)
        direct_path = tmp_path / "reference.annotation-db"
        original_path.replace(direct_path)

        def unexpected_validation(*args, **kwargs):
            raise AssertionError("direct database must not be text-validated")

        monkeypatch.setattr(annotation, "validate_annotation_file", unexpected_validation)
        ann = annotation.Annotation(str(direct_path), False, False)
        assert ann.backend == "gffbase"
        assert ann.cache_status == "direct_database"
        assert ann.db_connection["gene1"].id == "gene1"
        _close(ann.db_connection)
