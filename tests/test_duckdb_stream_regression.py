"""Regression coverage for GitHub #56's large ``--stream`` ingest failure.

The crash is in DuckDB's native-GEOMETRY row-group append path, not in
miniprot attribute parsing.  Keep these tests separate from the pipeline
fixtures: those pass ``-M`` and therefore never hand an in-memory blob to
gffbase.
"""

from __future__ import annotations

import importlib
import io
from pathlib import Path
from types import SimpleNamespace

import duckdb
import pytest

from lifton import gffbase_adapter
from lifton.gffbase import ingest


@pytest.mark.parametrize(
    ("version", "expected"),
    [
        ("1.5.2", False),
        ("1.5.2.dev9", False),
        ("1.5.3", True),
        ("1.5.3.dev12", True),
        ("v1.5.3-custom", True),
        ("1.5.4", True),
        ("1.5.4+cpu", True),
        ("1.5.5", False),
        ("1.5.5.dev1", False),
        ("1.6.0", False),
        ("not-a-version", False),
    ],
)
def test_geometry_shredding_workaround_version_boundaries(version, expected):
    assert ingest._duckdb_needs_geometry_shredding_workaround(version) is expected


@pytest.mark.parametrize(
    "env_name", ["LIFTON_DISABLE_RTREE", "GFFBASE_TEST_DISABLE_RTREE"]
)
@pytest.mark.parametrize("enabled_value", ["1", "true", "YES", "on"])
def test_rtree_disable_supported_and_legacy_env_names(
        monkeypatch, env_name, enabled_value):
    monkeypatch.delenv("LIFTON_DISABLE_RTREE", raising=False)
    monkeypatch.delenv("GFFBASE_TEST_DISABLE_RTREE", raising=False)
    monkeypatch.setenv(env_name, enabled_value)

    assert ingest._rtree_disabled_by_env() is True


def test_rtree_disable_env_false_values_do_not_disable(monkeypatch):
    monkeypatch.setenv("LIFTON_DISABLE_RTREE", "0")
    monkeypatch.setenv("GFFBASE_TEST_DISABLE_RTREE", "false")

    assert ingest._rtree_disabled_by_env() is False


@pytest.mark.parametrize(
    ("message", "expected"),
    [
        (
            "INTERNAL Error: Expected unified vector format of type LIST, "
            "but found type VARCHAR",
            True,
        ),
        ("Expected unified vector format of type STRUCT, but found type VARCHAR", False),
        ("Expected unified vector format of type LIST, but found type BLOB", False),
        ("ordinary DuckDB failure", False),
    ],
)
def test_geometry_append_error_signature_is_exact(message, expected):
    exc = duckdb.InternalException(message)
    assert ingest.is_geometry_append_error(exc) is expected


def test_geometry_append_error_rejects_non_duckdb_exception():
    message = "Expected unified vector format of type LIST, but found type VARCHAR"
    assert ingest.is_geometry_append_error(RuntimeError(message)) is False


@pytest.mark.parametrize("build_rtree", [True, False])
def test_create_db_propagates_build_rtree(monkeypatch, build_rtree):
    """The compatibility wrapper must not swallow this keyword in ``**kwargs``."""
    create_db_module = importlib.import_module("lifton.gffbase.create_db")
    seen = {}
    sentinel_connection = object()
    sentinel_stats = SimpleNamespace()

    def fake_from_file(path, **kwargs):
        seen.update(kwargs)
        return sentinel_connection, sentinel_stats

    monkeypatch.setattr(create_db_module._ingest, "from_file", fake_from_file)
    monkeypatch.setattr(
        create_db_module,
        "FeatureDB",
        lambda payload, **kwargs: payload,
    )

    result = create_db_module.create_db(
        "input.gff3", ":memory:", build_rtree=build_rtree
    )

    assert seen["build_rtree"] is build_rtree
    assert result == (sentinel_connection, sentinel_stats)


def test_create_db_retries_stream_blob_on_exact_geometry_error(monkeypatch):
    """Retry the same materialized stream on a fresh, geometry-free ingest."""
    create_db_module = importlib.import_module("lifton.gffbase.create_db")
    calls = []
    sentinel_connection = object()
    sentinel_stats = SimpleNamespace()

    def fake_from_file(path, **kwargs):
        calls.append(
            {
                "path": path,
                "contents": Path(path).read_text(),
                **kwargs,
            }
        )
        if len(calls) == 1:
            raise duckdb.InternalException(
                "Expected unified vector format of type LIST, "
                "but found type VARCHAR"
            )
        return sentinel_connection, sentinel_stats

    monkeypatch.setattr(create_db_module._ingest, "from_file", fake_from_file)
    monkeypatch.setattr(
        create_db_module,
        "FeatureDB",
        lambda payload, **kwargs: payload,
    )

    gff_text = "chr1\tminiprot\tmRNA\t1\t9\t.\t+\t.\tID=MP1\n"
    with pytest.warns(RuntimeWarning, match="fresh connection"):
        result = create_db_module.create_db(
            gff_text,
            ":memory:",
            from_string=True,
            build_rtree=True,
        )

    assert result == (sentinel_connection, sentinel_stats)
    assert [call["build_rtree"] for call in calls] == [True, False]
    assert [call["force"] for call in calls] == [False, True]
    assert calls[0]["path"] == calls[1]["path"]
    assert calls[0]["contents"] == calls[1]["contents"] == gff_text
    assert not Path(calls[0]["path"]).exists()


def test_create_db_does_not_retry_unrelated_duckdb_error(monkeypatch):
    create_db_module = importlib.import_module("lifton.gffbase.create_db")
    paths = []

    def fake_from_file(path, **kwargs):
        paths.append(path)
        raise duckdb.InternalException("unrelated internal failure")

    monkeypatch.setattr(create_db_module._ingest, "from_file", fake_from_file)

    with pytest.raises(duckdb.InternalException, match="unrelated"):
        create_db_module.create_db(
            "chr1\tminiprot\tmRNA\t1\t9\t.\t+\t.\tID=MP1\n",
            ":memory:",
            from_string=True,
            build_rtree=True,
        )

    assert len(paths) == 1
    assert not Path(paths[0]).exists()


def test_create_db_propagates_retry_failure_without_third_attempt(monkeypatch):
    create_db_module = importlib.import_module("lifton.gffbase.create_db")
    calls = []

    def fake_from_file(path, **kwargs):
        calls.append((path, kwargs["build_rtree"]))
        raise duckdb.InternalException(
            "Expected unified vector format of type LIST, "
            "but found type VARCHAR"
        )

    monkeypatch.setattr(create_db_module._ingest, "from_file", fake_from_file)

    with pytest.warns(RuntimeWarning, match="fresh connection"):
        with pytest.raises(duckdb.InternalException, match="unified vector"):
            create_db_module.create_db(
                "chr1\tminiprot\tmRNA\t1\t9\t.\t+\t.\tID=MP1\n",
                ":memory:",
                from_string=True,
                build_rtree=True,
            )

    assert [build_rtree for _path, build_rtree in calls] == [True, False]
    assert calls[0][0] == calls[1][0]
    assert not Path(calls[0][0]).exists()


@pytest.mark.parametrize(
    ("from_string", "build_rtree"),
    [
        (False, True),
        (True, False),
    ],
)
def test_create_db_does_not_retry_geometry_error_outside_stream_rtree(
        monkeypatch, from_string, build_rtree):
    """Only an R-tree failure while consuming a stream is retry-safe here."""
    create_db_module = importlib.import_module("lifton.gffbase.create_db")
    paths = []

    def fake_from_file(path, **kwargs):
        paths.append(path)
        raise duckdb.InternalException(
            "Expected unified vector format of type LIST, "
            "but found type VARCHAR"
        )

    monkeypatch.setattr(create_db_module._ingest, "from_file", fake_from_file)
    data = (
        "chr1\tminiprot\tmRNA\t1\t9\t.\t+\t.\tID=MP1\n"
        if from_string
        else "input.gff3"
    )

    with pytest.raises(duckdb.InternalException, match="unified vector"):
        create_db_module.create_db(
            data,
            ":memory:",
            from_string=from_string,
            build_rtree=build_rtree,
        )

    assert len(paths) == 1
    if from_string:
        assert not Path(paths[0]).exists()
    else:
        assert paths == ["input.gff3"]


@pytest.mark.parametrize("build_rtree", [True, False])
def test_blob_adapter_propagates_build_rtree(monkeypatch, build_rtree):
    seen = {}
    sentinel = object()

    def fake_create_db(data, dbfn, **kwargs):
        seen.update(kwargs)
        return sentinel

    monkeypatch.setattr(gffbase_adapter._gffbase, "create_db", fake_create_db)

    result = gffbase_adapter.build_database_from_string(
        gff_text=b"chr1\tminiprot\tmRNA\t1\t9\t.\t+\t.\tID=MP1\n",
        build_rtree=build_rtree,
    )

    assert seen["build_rtree"] is build_rtree
    assert result is sentinel


@pytest.mark.parametrize("build_rtree", [True, False])
def test_file_adapter_propagates_build_rtree(monkeypatch, build_rtree):
    seen = {}
    sentinel = object()

    def fake_create_db(data, dbfn, **kwargs):
        seen.update(kwargs)
        return sentinel

    monkeypatch.setattr(gffbase_adapter._gffbase, "create_db", fake_create_db)

    result = gffbase_adapter.build_database(
        file_name="input.gff3",
        infer_genes=False,
        infer_transcripts=False,
        build_rtree=build_rtree,
    )

    assert seen["build_rtree"] is build_rtree
    assert result is sentinel


def test_failed_spatial_batch_closes_invalidated_connection(monkeypatch):
    """Never issue unregister/SQL calls after DuckDB invalidates a connection."""
    builder = ingest._ArrowBatchBuilder({}, has_spatial=True)
    builder.f_id.append("MP1")
    monkeypatch.setattr(builder, "features_table", lambda: object())
    monkeypatch.setattr(builder, "attributes_table", lambda: object())

    class FailingConnection:
        def __init__(self):
            self.closed = False
            self.unregister_called = False

        def register(self, name, table):
            return None

        def execute(self, sql):
            raise duckdb.InternalException(
                "Expected unified vector format of type LIST, "
                "but found type VARCHAR"
            )

        def unregister(self, name):
            self.unregister_called = True

        def close(self):
            self.closed = True

    con = FailingConnection()
    with pytest.raises(duckdb.InternalException, match="unified vector"):
        builder.flush_into(con)

    assert con.closed is True
    assert con.unregister_called is False


def test_python_parser_closes_stream_after_exhaustion(monkeypatch):
    stream = io.StringIO(
        "chr1\tminiprot\tmRNA\t1\t9\t.\t+\t.\tID=MP1\n"
    )
    monkeypatch.setattr(ingest._parser._pyparser, "_open", lambda path: stream)

    records = list(ingest._parser.parse_gff("unused.gff3", engine="python"))

    assert len(records) == 1
    assert stream.closed is True


def _large_miniprot_blob(n_rows: int = 100_001) -> bytes:
    """Cross gffbase's two complete 50,000-row batches by exactly one row."""
    blob = bytearray(b"##gff-version 3\n")
    for i in range(n_rows):
        blob.extend(
            f"chr1\tminiprot\tmRNA\t{i + 1}\t{i + 9}\t.\t+\t.\t"
            f"ID=MP{i:07d};Rank=1;Identity=0.99;"
            f"Target=XP_{i:09d}.1 1 30\n".encode()
        )
    return bytes(blob)


def test_affected_duckdb_ingests_third_miniprot_batch_without_rtree():
    """DuckDB 1.5.3/1.5.4 used to fail on row 100,001 with an R-tree.

    Skip other DuckDB releases because they do not exercise the affected
    storage implementation. The version guard must route affected releases
    to the B-tree path before optional spatial-extension setup begins.
    """
    if not ingest._duckdb_needs_geometry_shredding_workaround():
        pytest.skip(f"DuckDB {duckdb.__version__} is outside the affected range")

    with pytest.warns(RuntimeWarning, match="known GEOMETRY append bug"):
        db = gffbase_adapter.build_database_from_string(
            gff_text=_large_miniprot_blob(),
            build_rtree=True,
        )

    try:
        assert db.count_features_of_type("mRNA") == 100_001
        assert db._rtree_built is False
    finally:
        db.conn.close()
