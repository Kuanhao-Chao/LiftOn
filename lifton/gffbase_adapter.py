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

from lifton import gffbase as _gffbase


def db_path_for(file_name: str) -> str:
    """LiftOn's gffutils backend caches its SQLite DB at ``<file>_db``.
    The gffbase backend produces a DuckDB file; use a different suffix
    so the two backends can coexist on disk without clobbering each
    other."""
    return file_name + ".duckdb"


def open_existing_db(file_name: str) -> Optional["_gffbase.FeatureDB"]:
    """Return a :class:`gffbase.FeatureDB` attached to an existing
    DuckDB cache for ``file_name`` if one exists, else ``None``.

    Mirrors the legacy ``try: gffutils.FeatureDB(<file>_db)`` short-
    circuit at :mod:`lifton.annotation`.
    """
    p = db_path_for(file_name)
    if not os.path.exists(p):
        return None
    try:
        return _gffbase.FeatureDB(p)
    except Exception:
        return None


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
) -> "_gffbase.FeatureDB":
    """Build a DuckDB-backed FeatureDB from a GFF3 file on disk.

    The Phase 6.2 design notes that gffbase's ingest already
    deduplicates IDs internally, so the legacy three-strategy retry
    in ``Annotation._build_database`` collapses to a single call.
    """
    return _gffbase.create_db(
        file_name,
        dbfn=db_path_for(file_name),
        force=force,
        verbose=verbose,
        merge_strategy=merge_strategy,
        id_spec=id_spec,
        transform=transform,
        disable_infer_genes=not infer_genes,
        disable_infer_transcripts=not infer_transcripts,
    )


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
