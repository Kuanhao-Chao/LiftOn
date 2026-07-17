# ---------------------------------------------------------------------------
# Phase 7 (Streaming Adapter Layer) — gffbase shim.
#
# Translates the few LiftOn-specific gffutils call patterns into gffbase
# semantics. Pure translation; no LiftOn algorithmic logic lives here.
#
# Three public functions:
#   - db_path_for(file_name)              — where the DuckDB cache goes
#   - open_existing_db(file_name)         — attach an existing DuckDB cache
#   - build_database(file_name, ...)      — ingest a path into a fresh DB
#   - build_database_from_string(blob)    — ingest a bytes blob held in RAM
# ---------------------------------------------------------------------------

from __future__ import annotations

import os
from typing import Callable, Optional, Union

import duckdb

from lifton import __version__ as _lifton_version
from lifton import annotation_cache as _cache
from lifton import gffbase as _gffbase
from lifton.gffbase.schema import SCHEMA_VERSION


def db_path_for(file_name: str) -> str:
    """LiftOn's gffutils backend caches its SQLite DB at ``<file>_db``.
    The gffbase backend produces a DuckDB file; use a different suffix
    so the two backends can coexist on disk without clobbering each
    other."""
    return file_name + ".duckdb"


def _cache_settings(
    *,
    infer_genes: bool,
    infer_transcripts: bool,
    merge_strategy: str,
    id_spec,
    transform: Optional[Callable],
    build_rtree: bool,
    parser_settings: Optional[dict],
) -> dict:
    """Return every setting that can affect gffbase ingestion."""
    return {
        "parser_contract_version": _cache.PARSER_CONTRACT_VERSION,
        "parser": {
            "native_available": _gffbase.native_available(),
            "settings": parser_settings or {},
        },
        "inference": {
            "genes": bool(infer_genes),
            "transcripts": bool(infer_transcripts),
        },
        "merge_strategy": merge_strategy,
        "id_spec": id_spec,
        "transform": transform,
        "build_rtree": bool(build_rtree),
        "rtree_disabled": bool(
            os.environ.get("LIFTON_DISABLE_RTREE")
            or os.environ.get("GFFBASE_TEST_DISABLE_RTREE")
        ),
        "duckdb_version": duckdb.__version__,
    }


def _expected_manifest(
    file_name: str,
    *,
    infer_genes: bool,
    infer_transcripts: bool,
    merge_strategy: str,
    id_spec,
    transform: Optional[Callable],
    build_rtree: bool,
    parser_settings: Optional[dict],
) -> dict:
    return _cache.build_manifest(
        file_name,
        backend="gffbase",
        backend_version=_gffbase.__version__,
        schema_version=SCHEMA_VERSION,
        tool_version=_lifton_version,
        settings=_cache_settings(
            infer_genes=infer_genes,
            infer_transcripts=infer_transcripts,
            merge_strategy=merge_strategy,
            id_spec=id_spec,
            transform=transform,
            build_rtree=build_rtree,
            parser_settings=parser_settings,
        ),
    )


def open_database_path(db_path: str) -> Optional["_gffbase.FeatureDB"]:
    """Open a direct gffbase database input after validating its schema."""
    if not os.path.isfile(db_path):
        return None
    db = None
    try:
        db = _gffbase.FeatureDB(os.fspath(db_path))
        row = db.conn.execute(
            "SELECT value FROM meta WHERE key = 'schema_version'"
        ).fetchone()
        if row is None or str(row[0]) != SCHEMA_VERSION:
            db.conn.close()
            return None
        db.conn.execute("SELECT 1 FROM features LIMIT 1").fetchone()
        return db
    except Exception:
        if db is not None:
            try:
                db.conn.close()
            except Exception:
                pass
        return None


def open_existing_db(
    file_name: str,
    *,
    infer_genes: bool = False,
    infer_transcripts: bool = False,
    merge_strategy: str = "create_unique",
    id_spec=None,
    transform: Optional[Callable] = None,
    build_rtree: bool = True,
    parser_settings: Optional[dict] = None,
) -> Optional["_gffbase.FeatureDB"]:
    """Return a :class:`gffbase.FeatureDB` attached to an existing
    DuckDB cache for ``file_name`` if one exists, else ``None``.

    Mirrors the legacy ``try: gffutils.FeatureDB(<file>_db)`` short-
    circuit at :mod:`lifton.annotation`.
    """
    db_path = db_path_for(file_name)
    try:
        expected = _expected_manifest(
            file_name,
            infer_genes=infer_genes,
            infer_transcripts=infer_transcripts,
            merge_strategy=merge_strategy,
            id_spec=id_spec,
            transform=transform,
            build_rtree=build_rtree,
            parser_settings=parser_settings,
        )
    except OSError:
        return None
    if not _cache.cache_matches(db_path, expected):
        return None
    return open_database_path(db_path)


def build_database(
    *,
    file_name: str,
    infer_genes: bool,
    infer_transcripts: bool,
    merge_strategy: str = "create_unique",
    id_spec: Optional[str] = None,
    force: bool = True,
    verbose: bool = False,
    transform: Optional[Callable] = None,
    build_rtree: bool = True,
    parser_settings: Optional[dict] = None,
) -> "_gffbase.FeatureDB":
    """Build a DuckDB-backed FeatureDB from a GFF3 file on disk.

    The Phase 6.2 design notes that gffbase's ingest already
    deduplicates IDs internally, so the legacy three-strategy retry
    in ``Annotation._build_database`` collapses to a single call.
    """
    db_path = db_path_for(file_name)
    expected = _expected_manifest(
        file_name,
        infer_genes=infer_genes,
        infer_transcripts=infer_transcripts,
        merge_strategy=merge_strategy,
        id_spec=id_spec,
        transform=transform,
        build_rtree=build_rtree,
        parser_settings=parser_settings,
    )
    if not force and _cache.cache_matches(db_path, expected):
        existing = open_database_path(db_path)
        if existing is not None:
            return existing

    temporary = _cache.temporary_db_path(db_path)
    built = None
    try:
        built = _gffbase.create_db(
            file_name,
            dbfn=temporary,
            force=True,
            verbose=verbose,
            merge_strategy=merge_strategy,
            id_spec=id_spec,
            transform=transform,
            disable_infer_genes=not infer_genes,
            disable_infer_transcripts=not infer_transcripts,
            build_rtree=build_rtree,
        )
        built.conn.close()
        built = None
        if _cache.source_fingerprint(file_name) != expected["source"]:
            raise OSError(f"Annotation changed while its database was being built: {file_name}")
        os.replace(temporary, db_path)
        _cache.write_manifest_atomic(db_path, expected)
        reopened = open_database_path(db_path)
        if reopened is None:
            raise OSError(f"Published gffbase database could not be reopened: {db_path}")
        return reopened
    finally:
        if built is not None:
            try:
                built.conn.close()
            except Exception:
                pass
        try:
            os.unlink(temporary)
        except FileNotFoundError:
            pass


def build_database_from_string(
    *,
    gff_text: Union[str, bytes, bytearray],
    dbfn: str = ":memory:",
    infer_genes: bool = False,
    infer_transcripts: bool = False,
    merge_strategy: str = "create_unique",
    id_spec: Optional[str] = None,
    force: bool = True,
    verbose: bool = False,
    transform: Optional[Callable] = None,
    build_rtree: bool = True,
) -> "_gffbase.FeatureDB":
    """Ingest a GFF3 blob held in RAM (no intermediate file).

    Phase 7 streaming-adapter entry point. Used by ``run_miniprot``
    when ``--stream`` is on so the miniprot stdout never touches
    disk.

    ``gff_text`` may be ``str`` or ``bytes``; gffbase's
    :func:`create_db` ``from_string=True`` materialises a tempfile
    internally, so callers do not have to.
    """
    if isinstance(gff_text, (bytes, bytearray)):
        gff_text = bytes(gff_text).decode("utf-8")
    return _gffbase.create_db(
        gff_text,
        dbfn=dbfn,
        from_string=True,
        force=force,
        verbose=verbose,
        merge_strategy=merge_strategy,
        id_spec=id_spec,
        transform=transform,
        disable_infer_genes=not infer_genes,
        disable_infer_transcripts=not infer_transcripts,
        build_rtree=build_rtree,
    )


def looks_like_gff3_blob(value) -> bool:
    """Cheap heuristic used by ``Annotation.__init__`` to decide
    between ``bytes``/``str`` blobs vs filesystem paths. A real GFF3
    blob always contains a tab character or a directive marker;
    plain filenames almost never do."""
    if isinstance(value, (bytes, bytearray)):
        sample = bytes(value)[:4096]
        return b"\t" in sample or sample.startswith(b"#")
    return False
