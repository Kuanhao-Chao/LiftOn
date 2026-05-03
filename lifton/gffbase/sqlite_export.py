# ---------------------------------------------------------------------------
# Author: Kuan-Hao Chao <kuanhao.chao@gmail.com>
# ---------------------------------------------------------------------------
"""``export_sqlite`` — write a legacy gffutils-format SQLite database from a
gffbase DuckDB connection.

Documented break: the ``attributes`` column in the exported file is the raw
col-9 bytes (UTF-8 string), not legacy-style JSON. Most downstream code reads
attributes via ``Feature.attributes`` and continues to work; raw-SQL queries
that depend on JSON-decoded attributes need migration.
"""

from __future__ import annotations

import json
import os
import sqlite3
from typing import Optional

import duckdb

from ._bins import bin_from_coords


_LEGACY_SCHEMA = """
CREATE TABLE features (
    id text,
    seqid text,
    source text,
    featuretype text,
    start int,
    end int,
    score text,
    strand text,
    frame text,
    attributes text,
    extra text,
    bin int,
    primary key (id)
);
CREATE TABLE relations (
    parent text,
    child text,
    level int,
    primary key (parent, child, level)
);
CREATE TABLE meta (
    dialect text,
    version text
);
CREATE TABLE directives (
    directive text
);
CREATE TABLE autoincrements (
    base text,
    n int,
    primary key (base)
);
CREATE TABLE duplicates (
    idspecid text,
    newid text,
    primary key (newid)
);
CREATE INDEX featuretype  ON features (featuretype);
CREATE INDEX seqidstartend ON features (seqid, start, end);
CREATE INDEX relationsparent ON relations (parent);
CREATE INDEX relationschild  ON relations (child);
CREATE INDEX binindex ON features (bin);
"""


def export_sqlite(con: duckdb.DuckDBPyConnection, path: str,
                  force: bool = False) -> str:
    """Write a legacy SQLite ``.db`` from the given DuckDB connection.

    Returns the absolute path on success.
    """
    if os.path.exists(path):
        if not force:
            raise ValueError(f"{path} already exists; pass force=True to overwrite")
        os.unlink(path)

    sqlite_con = sqlite3.connect(path)
    try:
        sqlite_con.executescript(_LEGACY_SCHEMA)

        # Stream features in file order.
        rows = con.execute(
            """
            SELECT id, seqid, source, featuretype, start, "end",
                   score, strand, frame,
                   CAST(attributes_blob AS VARCHAR) AS attributes,
                   CAST(extra_blob      AS VARCHAR) AS extra
            FROM features
            ORDER BY file_order NULLS LAST, id
            """
        ).fetchall()

        # Compute UCSC bin in Python (DuckDB has no native equivalent and
        # legacy SQLite users rely on this for `region()` queries).
        export_rows = []
        for r in rows:
            (fid, seqid, source, featuretype, start, end, score, strand,
             frame, attributes, extra) = r
            ucsc_bin = bin_from_coords(start, end) if start and end else None
            export_rows.append(
                (fid, seqid, source, featuretype, start, end, score, strand,
                 frame, attributes or "", extra or "", ucsc_bin)
            )
        sqlite_con.executemany(
            "INSERT INTO features VALUES (?,?,?,?,?,?,?,?,?,?,?,?)",
            export_rows,
        )

        # Closure → relations(parent, child, level=depth).
        rels = con.execute(
            "SELECT ancestor, descendant, depth FROM closure"
        ).fetchall()
        sqlite_con.executemany("INSERT INTO relations VALUES (?,?,?)", rels)

        # Meta — write the dialect (JSON) + version.
        meta = dict(con.execute("SELECT key, value FROM meta").fetchall())
        sqlite_con.execute(
            "INSERT INTO meta VALUES (?, ?)",
            (meta.get("dialect", json.dumps({"fmt": "gff3"})),
             "gffbase-export"),
        )

        # Directives.
        dirs = con.execute(
            "SELECT directive FROM directives ORDER BY seq"
        ).fetchall()
        sqlite_con.executemany("INSERT INTO directives VALUES (?)", dirs)

        # Autoincrements (typically empty in Phase 5).
        try:
            ai = con.execute("SELECT base, n FROM autoincrements").fetchall()
            if ai:
                sqlite_con.executemany("INSERT INTO autoincrements VALUES (?, ?)", ai)
        except duckdb.Error:
            pass

        sqlite_con.commit()
    finally:
        sqlite_con.close()
    return os.path.abspath(path)
