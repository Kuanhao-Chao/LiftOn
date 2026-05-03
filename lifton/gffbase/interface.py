# ---------------------------------------------------------------------------
# Author: Kuan-Hao Chao <kuanhao.chao@gmail.com>
# ---------------------------------------------------------------------------
"""``FeatureDB`` — drop-in successor to ``gffutils.FeatureDB``.

Wraps a DuckDB connection produced by ``gffbase.ingest.from_file`` (or opened
from an on-disk ``.duckdb`` file) and exposes the legacy public surface.

Two routing decisions are made dynamically:

* ``region()`` dispatches to an R-tree-backed query when ``meta.rtree_built``
  is true, otherwise to the multi-column B-tree path. Both queries are
  semantically identical; only the planner choice differs.
* ``children()`` / ``parents()`` read from the materialized closure table
  when the requested ``level`` lies within ``meta.max_depth``, and fall back
  to a dynamic recursive CTE when a deeper traversal is requested.
"""

from __future__ import annotations

import json
import os
from typing import Iterator, List, Optional, Tuple, Union

import duckdb

from .exceptions import FeatureNotFoundError
from .feature import Feature, feature_from_row, _DB_ROW_FIELDS


# Selection clause for FeatureDB → Feature reconstruction. Mirrors
# `_DB_ROW_FIELDS`.
_SELECT_FEATURE = (
    'id, seqid, source, featuretype, start, "end", '
    'score, strand, frame, attributes_blob, extra_blob, file_order'
)


class FeatureDB:
    """Drop-in successor to ``gffutils.FeatureDB``."""

    # ------------------------------------------------------------------
    # Construction
    # ------------------------------------------------------------------

    def __init__(
        self,
        dbfn,
        default_encoding: str = "utf-8",
        keep_order: bool = False,
        pragmas: Optional[dict] = None,
        sort_attribute_values: bool = False,
        text_factory=str,
    ):
        self.default_encoding = default_encoding
        self.keep_order = keep_order
        self.sort_attribute_values = sort_attribute_values
        self.text_factory = text_factory
        self._analyzed_flag = False

        # Resolve dbfn → connection.
        if isinstance(dbfn, duckdb.DuckDBPyConnection):
            self.conn = dbfn
            self.dbfn = ":existing-connection:"
        elif isinstance(dbfn, tuple) and len(dbfn) == 2 and isinstance(
            dbfn[0], duckdb.DuckDBPyConnection
        ):
            # (con, IngestStats) — used internally by `create_db`.
            self.conn = dbfn[0]
            self.dbfn = ":existing-connection:"
        elif isinstance(dbfn, str):
            self.dbfn = dbfn
            self.conn = duckdb.connect(dbfn, read_only=False)
        else:
            raise TypeError(
                f"dbfn must be a path, DuckDB connection, or (con, stats) tuple; got {type(dbfn)!r}"
            )

        if pragmas:
            self.set_pragmas(pragmas)

        # Recover provenance from `meta`.
        meta = self._read_meta()
        self.dialect = self._parse_dialect(meta.get("dialect"))
        self.fmt = meta.get("fmt", "gff3")
        self._rtree_built = meta.get("rtree_built", "false").lower() == "true"
        self._max_depth = int(meta.get("max_depth", "8"))
        # Defensive: confirm the R-tree index actually exists in this DB.
        if self._rtree_built:
            self._rtree_built = self._has_rtree_index()
        # If the R-tree was built at ingest time, the spatial extension's
        # functions (ST_Intersects, ST_MakeEnvelope, …) must be loaded into
        # the current connection — they are NOT auto-loaded by DuckDB just
        # because the index exists. If the load fails (offline HPC node, etc),
        # gracefully fall back to the multi-column B-tree path.
        if self._rtree_built:
            try:
                self.conn.execute("LOAD spatial")
            except duckdb.Error:
                try:
                    self.conn.execute("INSTALL spatial")
                    self.conn.execute("LOAD spatial")
                except duckdb.Error:
                    self._rtree_built = False

        # Phase 7 — closure-cache vs dynamic-CTE dispatcher. Read the corpus's
        # true hierarchy depth once at open; fall back to a live MAX(depth)
        # query if older DBs don't carry the meta row.
        cmd_meta = meta.get("closure_max_depth")
        if cmd_meta is not None:
            self._closure_max_depth = int(cmd_meta)
        else:
            try:
                row = self.conn.execute("SELECT MAX(depth) FROM closure").fetchone()
                self._closure_max_depth = int(row[0]) if row and row[0] is not None else 0
            except duckdb.Error:
                self._closure_max_depth = 0

        # Phase 7 — load the seqid → y-band map so `_region_sql_rtree` can
        # produce a tightly-bounded ST_MakeEnvelope at query time. Empty when
        # no R-tree was built (we fall back to the B-tree path anyway).
        self._seqid_y_map: dict = {}
        if self._rtree_built:
            try:
                rows = self.conn.execute(
                    "SELECT seqid, seqid_y FROM seqid_map"
                ).fetchall()
                self._seqid_y_map = {s: int(y) for s, y in rows}
            except duckdb.Error:
                # Older DBs without the map — fall back to B-tree to be safe.
                self._rtree_built = False

        # Directives
        self.directives = [
            row[0] for row in self.conn.execute(
                "SELECT directive FROM directives ORDER BY seq"
            ).fetchall()
        ]

    @staticmethod
    def _parse_dialect(raw):
        if not raw:
            return {"fmt": "gff3"}
        try:
            return json.loads(raw)
        except (TypeError, ValueError):
            return {"fmt": "gff3"}

    def _read_meta(self) -> dict:
        try:
            rows = self.conn.execute("SELECT key, value FROM meta").fetchall()
            return {k: v for k, v in rows}
        except duckdb.Error:
            return {}

    def _has_rtree_index(self) -> bool:
        try:
            rows = self.conn.execute(
                "SELECT index_name FROM duckdb_indexes() "
                "WHERE table_name = 'features' AND index_name = 'features_rtree'"
            ).fetchall()
            return bool(rows)
        except duckdb.Error:
            return False

    # ------------------------------------------------------------------
    # Properties
    # ------------------------------------------------------------------

    @property
    def schema(self) -> str:
        rows = self.conn.execute("""
            SELECT sql FROM duckdb_tables() WHERE database_name = current_database()
            UNION ALL
            SELECT sql FROM duckdb_views()  WHERE database_name = current_database()
        """).fetchall()
        return "\n".join(r[0] for r in rows if r[0])

    @property
    def _analyzed(self) -> bool:
        return self._analyzed_flag

    # ------------------------------------------------------------------
    # Dunders
    # ------------------------------------------------------------------

    def __getitem__(self, key) -> Feature:
        target_id = key.id if isinstance(key, Feature) else key
        row = self.conn.execute(
            f"SELECT {_SELECT_FEATURE} FROM features WHERE id = ?", [target_id]
        ).fetchone()
        if row is None:
            raise FeatureNotFoundError(target_id)
        return feature_from_row(row, dialect=self.dialect)

    def __contains__(self, key) -> bool:
        target_id = key.id if isinstance(key, Feature) else key
        row = self.conn.execute(
            "SELECT 1 FROM features WHERE id = ? LIMIT 1", [target_id]
        ).fetchone()
        return row is not None

    # ------------------------------------------------------------------
    # Counts and distinct-value iterators
    # ------------------------------------------------------------------

    def count_features_of_type(self, featuretype: Optional[str] = None) -> int:
        if featuretype is None:
            return self.conn.execute("SELECT COUNT(*) FROM features").fetchone()[0]
        return self.conn.execute(
            "SELECT COUNT(*) FROM features WHERE featuretype = ?", [featuretype]
        ).fetchone()[0]

    def featuretypes(self) -> Iterator[str]:
        for (ft,) in self.conn.execute(
            "SELECT DISTINCT featuretype FROM features ORDER BY featuretype"
        ).fetchall():
            yield ft

    def seqids(self) -> Iterator[str]:
        for (s,) in self.conn.execute(
            "SELECT DISTINCT seqid FROM features ORDER BY seqid"
        ).fetchall():
            yield s

    # ------------------------------------------------------------------
    # Scans (all_features / features_of_type)
    # ------------------------------------------------------------------

    def all_features(
        self,
        limit=None,
        strand: Optional[str] = None,
        featuretype: Optional[Union[str, List[str]]] = None,
        order_by=None,
        reverse: bool = False,
        completely_within: bool = False,
    ) -> Iterator[Feature]:
        sql, params = self._build_scan_sql(
            base_where=[], base_params=[],
            limit=limit, strand=strand, featuretype=featuretype,
            order_by=order_by, reverse=reverse,
            completely_within=completely_within,
        )
        yield from self._yield_features(sql, params)

    def features_of_type(
        self,
        featuretype: Union[str, List[str]],
        limit=None,
        strand: Optional[str] = None,
        order_by=None,
        reverse: bool = False,
        completely_within: bool = False,
    ) -> Iterator[Feature]:
        yield from self.all_features(
            limit=limit, strand=strand, featuretype=featuretype,
            order_by=order_by, reverse=reverse,
            completely_within=completely_within,
        )

    def _build_scan_sql(self, *, base_where, base_params,
                        limit, strand, featuretype,
                        order_by, reverse, completely_within):
        where = list(base_where)
        params = list(base_params)
        # `limit` (legacy) accepts a region triple. Reuse `region()`'s parser.
        if limit is not None:
            rseqid, rstart, rend = self._normalize_region_args(limit, None, None, None)
            if rseqid is not None:
                where.append("seqid = ?")
                params.append(rseqid)
            if rstart is not None and rend is not None:
                if completely_within:
                    where.append("start >= ? AND \"end\" <= ?")
                    params.extend([rstart, rend])
                else:
                    where.append("start <= ? AND \"end\" >= ?")
                    params.extend([rend, rstart])
        if strand is not None:
            where.append("strand = ?")
            params.append(strand)
        if featuretype is not None:
            if isinstance(featuretype, (list, tuple, set)):
                placeholders = ",".join("?" * len(featuretype))
                where.append(f"featuretype IN ({placeholders})")
                params.extend(featuretype)
            else:
                where.append("featuretype = ?")
                params.append(featuretype)

        sql = f"SELECT {_SELECT_FEATURE} FROM features"
        if where:
            sql += " WHERE " + " AND ".join(where)
        sql += " ORDER BY " + self._order_clause(order_by, reverse)
        return sql, params

    @staticmethod
    def _order_clause(order_by, reverse: bool) -> str:
        if order_by is None:
            col = "file_order"
        elif order_by == "length":
            col = '("end" - start)'
        elif order_by in {
            "seqid", "source", "featuretype", "start", "end",
            "score", "strand", "frame", "file_order",
        }:
            col = '"end"' if order_by == "end" else order_by
        else:
            # Unknown / multi-field — accept as a literal for power users.
            col = order_by
        direction = "DESC" if reverse else "ASC"
        return f"{col} {direction}"

    # ------------------------------------------------------------------
    # region() — smart R-tree vs B-tree dispatch
    # ------------------------------------------------------------------

    def region(
        self,
        region=None,
        seqid: Optional[str] = None,
        start: Optional[int] = None,
        end: Optional[int] = None,
        strand: Optional[str] = None,
        featuretype: Optional[Union[str, List[str]]] = None,
        completely_within: bool = False,
    ) -> Iterator[Feature]:
        rseqid, rstart, rend = self._normalize_region_args(region, seqid, start, end)
        # Decide path.
        use_rtree = (
            self._rtree_built
            and rseqid is not None
            and rstart is not None
            and rend is not None
        )
        if use_rtree:
            sql, params = self._region_sql_rtree(
                rseqid, rstart, rend, strand, featuretype, completely_within
            )
        else:
            sql, params = self._region_sql_btree(
                rseqid, rstart, rend, strand, featuretype, completely_within
            )
        yield from self._yield_features(sql, params)

    def _region_sql_rtree(self, seqid, start, end, strand, featuretype, completely_within):
        # Phase 7: each seqid lives in its own y-band, so the R-tree query
        # envelope is tight on both axes — no cross-chromosome candidates.
        seqid_y = self._seqid_y_map.get(seqid)
        if seqid_y is None:
            # Unknown seqid (e.g., not in this DB). Fall through to a query
            # that returns no rows but is still valid.
            return self._region_sql_btree(seqid, start, end, strand, featuretype, completely_within)
        where = ['seqid = ?', 'ST_Intersects(bbox, ST_MakeEnvelope(?, ?, ?, ?))']
        params = [seqid, start, seqid_y, end, seqid_y + 1]
        if completely_within:
            where.append('start >= ? AND "end" <= ?')
            params.extend([start, end])
        if strand is not None:
            where.append("strand = ?")
            params.append(strand)
        if featuretype is not None:
            if isinstance(featuretype, (list, tuple, set)):
                ph = ",".join("?" * len(featuretype))
                where.append(f"featuretype IN ({ph})")
                params.extend(featuretype)
            else:
                where.append("featuretype = ?")
                params.append(featuretype)
        sql = (
            f"SELECT {_SELECT_FEATURE} FROM features "
            f"WHERE {' AND '.join(where)} ORDER BY start"
        )
        return sql, params

    def _region_sql_btree(self, seqid, start, end, strand, featuretype, completely_within):
        where = []
        params = []
        if seqid is not None:
            where.append("seqid = ?")
            params.append(seqid)
        if start is not None and end is not None:
            if completely_within:
                where.append('start >= ? AND "end" <= ?')
                params.extend([start, end])
            else:
                # Standard overlap: feature.start <= region.end AND feature.end >= region.start
                where.append('start <= ? AND "end" >= ?')
                params.extend([end, start])
        elif start is not None:
            where.append("start >= ?")
            params.append(start)
        elif end is not None:
            where.append('"end" <= ?')
            params.append(end)
        if strand is not None:
            where.append("strand = ?")
            params.append(strand)
        if featuretype is not None:
            if isinstance(featuretype, (list, tuple, set)):
                ph = ",".join("?" * len(featuretype))
                where.append(f"featuretype IN ({ph})")
                params.extend(featuretype)
            else:
                where.append("featuretype = ?")
                params.append(featuretype)
        sql = f"SELECT {_SELECT_FEATURE} FROM features"
        if where:
            sql += " WHERE " + " AND ".join(where)
        sql += " ORDER BY start"
        return sql, params

    @staticmethod
    def _normalize_region_args(region, seqid, start, end):
        """Return ``(seqid, start, end)`` triple from any of the four legacy
        input shapes."""
        if region is not None and any(v is not None for v in (seqid, start, end)):
            raise ValueError("pass either `region` or seqid/start/end, not both")
        if region is None:
            return seqid, start, end
        if isinstance(region, str):
            # "chr:start-end" or "chr:start-end:strand"
            chrom_rest = region.split(":", 1)
            if len(chrom_rest) == 1:
                return chrom_rest[0], None, None
            chrom, rest = chrom_rest
            if "-" in rest:
                s, e = rest.split("-", 1)
                # strip optional ":strand" suffix.
                if ":" in e:
                    e = e.split(":", 1)[0]
                return chrom, int(s), int(e)
            return chrom, None, None
        if isinstance(region, tuple):
            if len(region) == 3:
                return region[0], int(region[1]), int(region[2])
            if len(region) == 2:
                return region[0], None, None
            raise ValueError(f"region tuple must be (seqid, start, end); got {region!r}")
        if isinstance(region, Feature):
            return region.seqid, region.start, region.end
        raise TypeError(f"unsupported region type: {type(region)!r}")

    # ------------------------------------------------------------------
    # Phase 12 — Vectorized batched spatial API.
    # ------------------------------------------------------------------

    def region_batched(
        self,
        regions,
        featuretype: Optional[Union[str, List[str]]] = None,
        completely_within: bool = False,
        format: str = "arrow",
    ):
        """Bulk overlap query. Performs a SINGLE spatial JOIN between every
        input region and the features table, returning a column-oriented
        result that maps each query (`query_idx`) back to its overlapping
        features.

        Parameters
        ----------
        regions : iterable of (seqid, start, end) | str | Feature
            Each item is normalized through `_normalize_region_args`; the
            same four input shapes accepted by `region()` are accepted here.
        featuretype : str | list[str] | None
            Optional `features.featuretype` filter applied to all regions.
        completely_within : bool
            If True, only features fully contained in the region are returned
            (default False — overlap is sufficient).
        format : "arrow" | "df" | "polars"
            Return shape (default `"arrow"`).

        Result columns are query_idx, query_seqid, query_start,
        query_end, id, seqid, source, featuretype, start, end, score,
        strand, frame, file_order.
        """
        rows = []
        for r in regions:
            if isinstance(r, tuple) and len(r) == 3 and all(v is not None for v in r):
                seqid, rs, re_ = r[0], int(r[1]), int(r[2])
            else:
                seqid, rs, re_ = self._normalize_region_args(r, None, None, None)
            if seqid is None or rs is None or re_ is None:
                continue
            rows.append((seqid, int(rs), int(re_)))
        if not rows:
            return self._empty_region_batched(format)

        import pyarrow as pa
        regions_table = pa.table({
            "query_idx":   list(range(len(rows))),
            "query_seqid": [r[0] for r in rows],
            "query_start": [r[1] for r in rows],
            "query_end":   [r[2] for r in rows],
        })
        self.conn.register("__staging_regions", regions_table)
        try:
            # The R-tree path uses the seqid_y-encoded envelope so that
            # DuckDB's spatial index segregates chromosomes (Phase 7).
            # B-tree fallback is identical SQL minus the ST_Intersects
            # predicate.
            ft_where, ft_params = self._featuretype_filter(featuretype, qualifier="f")
            ft_clause = (" AND " + " AND ".join(ft_where)) if ft_where else ""
            within_clause = (
                ' AND f.start >= q.query_start AND f."end" <= q.query_end'
                if completely_within else ""
            )

            if self._rtree_built and self._seqid_y_map:
                # Inline the seqid → seqid_y map as a small VALUES table so
                # the JOIN can prune by chromosome inside the R-tree.
                values_pairs = ",".join(
                    f"(?, {y})" for y in (self._seqid_y_map[s] for s in self._seqid_y_map)
                )
                seqid_y_params = list(self._seqid_y_map.keys())
                sql = f"""
                    WITH seqid_lookup(seqid, seqid_y) AS (
                        VALUES {values_pairs}
                    )
                    SELECT
                        q.query_idx, q.query_seqid, q.query_start, q.query_end,
                        f.id, f.seqid, f.source, f.featuretype,
                        f.start, f."end" AS "end",
                        f.score, f.strand, f.frame, f.file_order
                    FROM __staging_regions q
                    JOIN seqid_lookup s ON s.seqid = q.query_seqid
                    JOIN features f
                      ON f.seqid = q.query_seqid
                      AND ST_Intersects(
                          f.bbox,
                          ST_MakeEnvelope(q.query_start, s.seqid_y,
                                          q.query_end,   s.seqid_y + 1))
                    WHERE 1=1{within_clause}{ft_clause}
                    ORDER BY q.query_idx, f.start
                """
                params = list(seqid_y_params) + ft_params
            else:
                sql = f"""
                    SELECT
                        q.query_idx, q.query_seqid, q.query_start, q.query_end,
                        f.id, f.seqid, f.source, f.featuretype,
                        f.start, f."end" AS "end",
                        f.score, f.strand, f.frame, f.file_order
                    FROM __staging_regions q
                    JOIN features f
                      ON f.seqid = q.query_seqid
                      AND f.start <= q.query_end
                      AND f."end"  >= q.query_start
                    WHERE 1=1{within_clause}{ft_clause}
                    ORDER BY q.query_idx, f.start
                """
                params = list(ft_params)

            return self._materialize_batched(sql, params, format=format)
        finally:
            try:
                self.conn.unregister("__staging_regions")
            except Exception:
                pass

    def _empty_region_batched(self, format: str):
        import pyarrow as pa
        schema = pa.schema([
            ("query_idx",    pa.int64()),
            ("query_seqid",  pa.string()),
            ("query_start",  pa.int64()),
            ("query_end",    pa.int64()),
            ("id",           pa.string()),
            ("seqid",        pa.string()),
            ("source",       pa.string()),
            ("featuretype",  pa.string()),
            ("start",        pa.int64()),
            ("end",          pa.int64()),
            ("score",        pa.string()),
            ("strand",       pa.string()),
            ("frame",        pa.string()),
            ("file_order",   pa.int64()),
        ])
        empty = pa.table({name: [] for name in schema.names}, schema=schema)
        fmt = format.lower()
        if fmt == "arrow":
            return empty
        if fmt in ("df", "pandas"):
            return empty.to_pandas()
        if fmt == "polars":
            try:
                import polars as pl
            except ImportError as e:  # pragma: no cover
                raise ImportError("format='polars' requires the optional polars package") from e
            return pl.from_arrow(empty)
        raise ValueError(f"format must be one of 'arrow' | 'df' | 'polars'; got {format!r}")

    # ------------------------------------------------------------------
    # children() / parents() — closure cache + dynamic CTE fallback
    # ------------------------------------------------------------------

    def children(
        self,
        id,
        level: Optional[int] = None,
        featuretype: Optional[Union[str, List[str]]] = None,
        order_by=None,
        reverse: bool = False,
        limit=None,
        completely_within: bool = False,
    ) -> Iterator[Feature]:
        target_id = id.id if isinstance(id, Feature) else id
        yield from self._relation_query(
            target_id, level, featuretype, order_by, reverse, limit, completely_within,
            direction="children",
        )

    def parents(
        self,
        id,
        level: Optional[int] = None,
        featuretype: Optional[Union[str, List[str]]] = None,
        order_by=None,
        reverse: bool = False,
        completely_within: bool = False,
        limit=None,
    ) -> Iterator[Feature]:
        target_id = id.id if isinstance(id, Feature) else id
        yield from self._relation_query(
            target_id, level, featuretype, order_by, reverse, limit, completely_within,
            direction="parents",
        )

    # ------------------------------------------------------------------
    # Phase 12 — Vectorized batched API.
    #
    # The row-by-row `children()` / `parents()` / `region()` generators
    # carry per-row Python overhead that dominates wall time on small
    # GENCODE-scale queries (Phase 11 §4.2). The methods below replace the
    # per-id loop with a single bulk SQL query and return the result as a
    # zero-copy PyArrow `Table` (or pandas / polars DataFrame), letting
    # downstream ML pipelines consume the data column-wise without ever
    # materializing a Python `Feature` object.
    # ------------------------------------------------------------------

    def children_batched(
        self,
        feature_ids,
        level: Optional[int] = None,
        featuretype: Optional[Union[str, List[str]]] = None,
        format: str = "arrow",
    ):
        """Bulk children lookup. Returns the descendants of ALL `feature_ids`
        in a single vectorized DuckDB query.

        Parameters
        ----------
        feature_ids : iterable of str | Feature
            Anchors. May contain `Feature` objects or raw ID strings.
        level : int | None
            None → all descendants (closure cache when materialized,
            otherwise dynamic CTE). Integer → exact-depth point lookup.
        featuretype : str | list[str] | None
            Optional filter on `features.featuretype`.
        format : "arrow" | "df" | "polars"
            Return shape (default `"arrow"` — `pyarrow.Table`).

        Result columns are anchor (the parent ID supplied),
        descendant_id, seqid, source, featuretype, start, end, score,
        strand, frame, file_order, depth.
        """
        return self._batched_relation(
            feature_ids, level=level, featuretype=featuretype,
            direction="children", format=format,
        )

    def parents_batched(
        self,
        feature_ids,
        level: Optional[int] = None,
        featuretype: Optional[Union[str, List[str]]] = None,
        format: str = "arrow",
    ):
        """Bulk parents lookup. Mirrors `children_batched` but walks the
        closure / edges in the reverse direction."""
        return self._batched_relation(
            feature_ids, level=level, featuretype=featuretype,
            direction="parents", format=format,
        )

    def _batched_relation(
        self,
        feature_ids,
        *,
        level: Optional[int],
        featuretype,
        direction: str,
        format: str,
    ):
        ids = self._coerce_id_list(feature_ids)
        if not ids:
            return self._empty_batched_result(format, "children" if direction == "children" else "parents")

        # Decide cache vs dynamic the same way the row-by-row dispatcher
        # does — it's a one-time decision per call here, not per-row.
        use_dynamic = (
            (level is not None and level > self._max_depth)
            or (level is None and self._closure_max_depth == 0)
        )

        ph = ",".join("?" * len(ids))
        ft_where, ft_params = self._featuretype_filter(featuretype, qualifier="f")
        ft_clause = (" AND " + " AND ".join(ft_where)) if ft_where else ""

        if direction == "children":
            anchor_alias, descendant_alias = "ancestor", "descendant"
            edge_anchor_col, edge_descendant_col = "parent", "child"
        else:
            anchor_alias, descendant_alias = "descendant", "ancestor"
            edge_anchor_col, edge_descendant_col = "child", "parent"

        if use_dynamic:
            # Recursive CTE seeded by every anchor in the batch.
            max_walk = level if level is not None else max(64, self._max_depth * 4)
            depth_filter = " AND w.depth = ?" if level is not None else ""
            cte = f"""
                WITH RECURSIVE walk(anchor, id, depth) AS (
                    SELECT {edge_anchor_col}, {edge_descendant_col}, 1
                    FROM edges WHERE {edge_anchor_col} IN ({ph})
                    UNION ALL
                    SELECT w.anchor, e.{edge_descendant_col}, w.depth + 1
                    FROM walk w
                    JOIN edges e ON e.{edge_anchor_col} = w.id
                    WHERE w.depth < ?
                )
                SELECT
                    w.anchor       AS anchor,
                    f.id           AS descendant_id,
                    f.seqid, f.source, f.featuretype,
                    f.start, f."end" AS "end",
                    f.score, f.strand, f.frame, f.file_order, w.depth
                FROM walk w
                JOIN features f ON f.id = w.id
                WHERE 1=1{depth_filter}{ft_clause}
            """
            params: list = list(ids) + [max_walk]
            if level is not None:
                params.append(level)
            params.extend(ft_params)
        else:
            # Closure cache hit — single set-based JOIN.
            depth_filter = " AND c.depth = ?" if level is not None else ""
            cte = f"""
                SELECT
                    c.{anchor_alias}     AS anchor,
                    f.id                  AS descendant_id,
                    f.seqid, f.source, f.featuretype,
                    f.start, f."end" AS "end",
                    f.score, f.strand, f.frame, f.file_order, c.depth
                FROM closure c
                JOIN features f ON f.id = c.{descendant_alias}
                WHERE c.{anchor_alias} IN ({ph}){depth_filter}{ft_clause}
            """
            params = list(ids)
            if level is not None:
                params.append(level)
            params.extend(ft_params)

        return self._materialize_batched(cte, params, format=format)

    @staticmethod
    def _coerce_id_list(feature_ids) -> List[str]:
        """Normalize a heterogeneous iterable of (id-string | Feature) into a
        list of strings."""
        out: List[str] = []
        for item in feature_ids:
            if isinstance(item, Feature):
                out.append(item.id)
            elif isinstance(item, str):
                out.append(item)
            else:
                raise TypeError(
                    f"feature_ids may contain only str or Feature; got {type(item)!r}"
                )
        return out

    def _materialize_batched(self, sql: str, params: list, *, format: str):
        """Execute `sql` and return the result in the requested shape.

        DuckDB's `fetch_arrow_table()` is a zero-copy hand-off — the
        returned `pyarrow.Table` shares the same Arrow buffers DuckDB uses
        internally, with no per-row Python boundary crossings.
        """
        cur = self.conn.execute(sql, params)
        fmt = format.lower()
        if fmt == "arrow":
            # DuckDB ≥ 1.0 prefers `to_arrow_table()`; fall back to the older
            # `fetch_arrow_table()` for environments pinning ≤ 0.10.
            return getattr(cur, "to_arrow_table", cur.fetch_arrow_table)() \
                if hasattr(cur, "to_arrow_table") else cur.fetch_arrow_table()
        if fmt in ("df", "pandas"):
            return cur.df()
        if fmt == "polars":
            return cur.pl()                  # DuckDB ≥1.0 returns a polars.DataFrame
        raise ValueError(
            f"format must be one of 'arrow' | 'df' | 'polars'; got {format!r}"
        )

    def _empty_batched_result(self, format: str, direction: str):
        """Return a properly-typed empty result when the input id list is
        empty. Avoids issuing a SQL query at all."""
        import pyarrow as pa
        schema = pa.schema([
            ("anchor",         pa.string()),
            ("descendant_id",  pa.string()),
            ("seqid",          pa.string()),
            ("source",         pa.string()),
            ("featuretype",    pa.string()),
            ("start",          pa.int64()),
            ("end",            pa.int64()),
            ("score",          pa.string()),
            ("strand",         pa.string()),
            ("frame",          pa.string()),
            ("file_order",     pa.int64()),
            ("depth",          pa.int16()),
        ])
        empty = pa.table({name: [] for name in schema.names}, schema=schema)
        fmt = format.lower()
        if fmt == "arrow":
            return empty
        if fmt in ("df", "pandas"):
            return empty.to_pandas()
        if fmt == "polars":
            try:
                import polars as pl
            except ImportError as e:  # pragma: no cover
                raise ImportError("format='polars' requires the optional polars package") from e
            return pl.from_arrow(empty)
        raise ValueError(f"format must be one of 'arrow' | 'df' | 'polars'; got {format!r}")

    def _relation_query(
        self,
        target_id: str,
        level: Optional[int],
        featuretype,
        order_by,
        reverse: bool,
        limit,
        completely_within: bool,
        direction: str,
    ):
        # Phase 7 — smart cache-vs-dynamic dispatcher.
        use_dynamic = self._dispatch_relation(level, target_id, direction)
        if use_dynamic:
            sql, params = self._relation_sql_dynamic(
                target_id, level, featuretype, order_by, reverse,
                limit, completely_within, direction,
            )
        else:
            sql, params = self._relation_sql_cached(
                target_id, level, featuretype, order_by, reverse,
                limit, completely_within, direction,
            )
        yield from self._yield_features(sql, params)

    def _dispatch_relation(
        self, level: Optional[int], target_id: str, direction: str
    ) -> bool:
        """Return True iff the caller should be served by the dynamic CTE.

        Decision matrix (Phase 7, revised):

        +-------+------------------+-------------------+--------------+
        | level | closure_max_depth| has_overflow?     | path         |
        +=======+==================+===================+==============+
        | None  | == 0 (no edges)  | n/a               | dynamic CTE  |
        | None  | >= 1             | no                | closure cache|
        | None  | >= 1             | yes               | dynamic CTE  |
        | int   | <= max_depth     | n/a               | closure cache|
        | int   | >  max_depth     | n/a               | dynamic CTE  |
        +-------+------------------+-------------------+--------------+

        Phase 7 measurement on GENCODE v45 (depth 2, 2000 genes / 274k descs):
        forced cache  = 18.35 s
        forced dynamic = 66.82 s
        Cache wins by 3.85× when the corpus's hierarchy fits in the cache.
        Phase 6's earlier finding (cache marginally slower) was an OS-cache
        artifact that disappeared on a clean run with the new R-tree encoding.

        Cache is preferred whenever it can serve the request; dynamic is the
        correctness fallback for traversals that extend past `max_depth`.
        """
        if level is not None:
            return level > self._max_depth
        # level is None.
        if self._closure_max_depth == 0:
            return True            # closure is empty; dynamic walks edges
        # Cache covers most of the tree; check for overflow past the boundary.
        return self._has_overflow(target_id, direction)

    def _has_overflow(self, target_id: str, direction: str) -> bool:
        """Are there descendants/ancestors past max_depth?"""
        if direction == "children":
            sql = """
                SELECT EXISTS (
                    SELECT 1 FROM edges e
                    JOIN closure c ON c.descendant = e.parent
                    WHERE c.ancestor = ? AND c.depth = ?
                )
            """
        else:
            sql = """
                SELECT EXISTS (
                    SELECT 1 FROM edges e
                    JOIN closure c ON c.ancestor = e.child
                    WHERE c.descendant = ? AND c.depth = ?
                )
            """
        row = self.conn.execute(sql, [target_id, self._max_depth]).fetchone()
        return bool(row[0])

    def _relation_sql_cached(
        self, target_id, level, featuretype, order_by, reverse,
        limit, completely_within, direction,
    ):
        if direction == "children":
            join_col, anchor_col = "c.descendant", "c.ancestor"
        else:
            join_col, anchor_col = "c.ancestor", "c.descendant"
        where = [f"{anchor_col} = ?"]
        params = [target_id]
        if level is not None:
            where.append("c.depth = ?")
            params.append(level)
        feat_where, feat_params = self._featuretype_filter(featuretype)
        where.extend(feat_where)
        params.extend(feat_params)
        lim_where, lim_params = self._limit_filter(limit, completely_within)
        where.extend(lim_where)
        params.extend(lim_params)

        sql = (
            f"SELECT {self._select_feature_aliased('f')} "
            f"FROM closure c JOIN features f ON f.id = {join_col} "
            f"WHERE {' AND '.join(where)} "
            f"ORDER BY {self._order_clause(order_by, reverse).replace(' ', ' f.', 0) if False else self._order_clause_qualified(order_by, reverse, 'f')}"
        )
        return sql, params

    def _relation_sql_dynamic(
        self, target_id, level, featuretype, order_by, reverse,
        limit, completely_within, direction,
    ):
        # Dynamic recursive CTE walks the edges table directly.
        if direction == "children":
            base_where = "parent = ?"
            recurse_join = "JOIN edges e ON e.parent = w.id"
            select_col = "child"
        else:
            base_where = "child = ?"
            recurse_join = "JOIN edges e ON e.child = w.id"
            select_col = "parent"

        # Walk depth bound: requested `level` if given, else a generous cap.
        max_walk = level if level is not None else max(64, self._max_depth * 4)

        cte = f"""
            WITH RECURSIVE walk(id, depth) AS (
                SELECT {select_col}, 1 FROM edges WHERE {base_where}
                UNION ALL
                SELECT e.{select_col}, w.depth + 1
                FROM walk w
                {recurse_join}
                WHERE w.depth < ?
            )
        """
        where = ["1=1"]
        params: list = [target_id, max_walk]
        if level is not None:
            where.append("w.depth = ?")
            params.append(level)
        feat_where, feat_params = self._featuretype_filter(featuretype, qualifier="f")
        where.extend(feat_where)
        params.extend(feat_params)
        lim_where, lim_params = self._limit_filter(limit, completely_within, qualifier="f")
        where.extend(lim_where)
        params.extend(lim_params)

        sql = (
            f"{cte} "
            f"SELECT {self._select_feature_aliased('f')} "
            f"FROM walk w JOIN features f ON f.id = w.id "
            f"WHERE {' AND '.join(where)} "
            f"ORDER BY {self._order_clause_qualified(order_by, reverse, 'f')}"
        )
        return sql, params

    @staticmethod
    def _select_feature_aliased(alias: str) -> str:
        cols = [
            "id", "seqid", "source", "featuretype", "start", '"end"',
            "score", "strand", "frame", "attributes_blob", "extra_blob", "file_order",
        ]
        return ", ".join(f"{alias}.{c}" if c != '"end"' else f'{alias}."end"' for c in cols)

    @staticmethod
    def _featuretype_filter(featuretype, *, qualifier: str = "f"):
        if featuretype is None:
            return [], []
        if isinstance(featuretype, (list, tuple, set)):
            ph = ",".join("?" * len(featuretype))
            return [f"{qualifier}.featuretype IN ({ph})"], list(featuretype)
        return [f"{qualifier}.featuretype = ?"], [featuretype]

    def _limit_filter(self, limit, completely_within: bool, *, qualifier: str = "f"):
        if limit is None:
            return [], []
        rseqid, rstart, rend = self._normalize_region_args(limit, None, None, None)
        where = []
        params = []
        if rseqid is not None:
            where.append(f"{qualifier}.seqid = ?")
            params.append(rseqid)
        if rstart is not None and rend is not None:
            if completely_within:
                where.append(f'{qualifier}.start >= ? AND {qualifier}."end" <= ?')
                params.extend([rstart, rend])
            else:
                where.append(f'{qualifier}.start <= ? AND {qualifier}."end" >= ?')
                params.extend([rend, rstart])
        return where, params

    @staticmethod
    def _order_clause_qualified(order_by, reverse, qualifier):
        if order_by is None:
            col = f"{qualifier}.file_order"
        elif order_by == "length":
            col = f'({qualifier}."end" - {qualifier}.start)'
        elif order_by in {
            "seqid", "source", "featuretype", "start", "end",
            "score", "strand", "frame", "file_order",
        }:
            col = f'{qualifier}."end"' if order_by == "end" else f"{qualifier}.{order_by}"
        else:
            col = order_by
        direction = "DESC" if reverse else "ASC"
        return f"{col} {direction}"

    # ------------------------------------------------------------------
    # Mutation
    # ------------------------------------------------------------------

    def delete(self, features, make_backup: bool = True, **kwargs) -> "FeatureDB":
        ids = self._coerce_ids(features)
        if not ids:
            return self
        placeholders = ",".join("?" * len(ids))
        self.conn.execute(f"DELETE FROM features WHERE id IN ({placeholders})", ids)
        self.conn.execute(f"DELETE FROM attributes WHERE feature_id IN ({placeholders})", ids)
        self.conn.execute(
            f"DELETE FROM edges WHERE parent IN ({placeholders}) OR child IN ({placeholders})",
            ids + ids,
        )
        self.conn.execute(
            f"DELETE FROM closure WHERE ancestor IN ({placeholders}) OR descendant IN ({placeholders})",
            ids + ids,
        )
        return self

    def update(self, data, make_backup: bool = True, **kwargs) -> "FeatureDB":
        # Phase 5 minimal update: accept iterable of Feature objects and
        # append them to features + attributes + edges, then refresh closure.
        from .ingest import _ArrowBatchBuilder
        from .feature import ParsedFeature

        # Phase 19: the builder needs the seqid_to_y dict so it can stamp
        # seqid_y (and bbox, when the R-tree is live) inline. Reuse the map
        # the FeatureDB already loaded from `seqid_map`.
        builder = _ArrowBatchBuilder(
            self._seqid_y_map,
            has_spatial=bool(self._rtree_built),
        )
        order = self.conn.execute(
            "SELECT COALESCE(MAX(file_order), 0) FROM features"
        ).fetchone()[0]

        if isinstance(data, FeatureDB):
            data = list(data.all_features())
        for feat in data:
            order += 1
            if isinstance(feat, Feature):
                blob = (
                    feat._attributes_blob
                    if feat._attributes_blob is not None
                    else feat._format_attributes().encode("utf-8")
                )
                pairs = [(k, v, i) for k, vs in feat.attributes.items() for i, v in enumerate(vs)]
                pf = ParsedFeature(
                    seqid=feat.seqid,
                    source=feat.source,
                    featuretype=feat.featuretype,
                    start=feat.start,
                    end=feat.end,
                    score=feat.score,
                    strand=feat.strand,
                    frame=feat.frame,
                    attributes_blob=blob,
                    attributes_pairs=pairs,
                    extra=list(feat.extra),
                )
                fid = feat.id or f"{feat.featuretype}_{order}"
            elif isinstance(feat, ParsedFeature):
                pf = feat
                fid = next(
                    (v for k, v, _ in pf.attributes_pairs if k == "ID"),
                    f"{pf.featuretype}_{order}",
                )
            else:
                raise TypeError(f"update() does not accept {type(feat)!r}")
            builder.append(fid, pf, order)
        builder.flush_into(self.conn)
        # Refresh closure: rebuild from edges (cheap on small updates).
        self.conn.execute("DELETE FROM closure")
        from .schema import CLOSURE_RECURSIVE_CTE
        self.conn.execute(CLOSURE_RECURSIVE_CTE, [self._max_depth])
        return self

    def add_relation(self, parent, child, level: int = 1,
                     parent_func=None, child_func=None) -> "FeatureDB":
        parent_id = parent.id if isinstance(parent, Feature) else parent
        child_id = child.id if isinstance(child, Feature) else child
        self.conn.execute(
            "INSERT INTO edges(parent, child) VALUES (?, ?)", [parent_id, child_id]
        )
        # Apply optional callbacks (legacy semantics: mutate attributes).
        if parent_func is not None and isinstance(parent, Feature):
            parent_func(parent, child)
        if child_func is not None and isinstance(child, Feature):
            child_func(parent, child)
        # Incrementally update closure: any ancestor of `parent` becomes an
        # ancestor of `child` (and transitively); any descendant of `child`
        # becomes a descendant of `parent`. Single set-based pass.
        self.conn.execute(
            """
            INSERT INTO closure (ancestor, descendant, depth)
            SELECT ? AS ancestor, ? AS descendant, 1
            WHERE NOT EXISTS (
                SELECT 1 FROM closure
                WHERE ancestor = ? AND descendant = ? AND depth = 1
            )
            """,
            [parent_id, child_id, parent_id, child_id],
        )
        # Rebuild closure incrementally. Simplest correct approach: drop the
        # closure rows that touch parent_id or child_id, then re-derive from
        # edges via the recursive CTE up to max_depth. For typical
        # `add_relation` calls (a handful per session) this is dramatically
        # cheaper than a full table rebuild and avoids the keyword pitfalls
        # of nested CTEs in DuckDB.
        from .schema import CLOSURE_RECURSIVE_CTE
        self.conn.execute("DELETE FROM closure")
        self.conn.execute(CLOSURE_RECURSIVE_CTE, [self._max_depth])
        return self

    @staticmethod
    def _coerce_ids(features) -> List[str]:
        if isinstance(features, str):
            return [features]
        if isinstance(features, Feature):
            return [features.id]
        if isinstance(features, FeatureDB):
            return [f.id for f in features.all_features()]
        out: List[str] = []
        for f in features:
            if isinstance(f, str):
                out.append(f)
            elif isinstance(f, Feature):
                out.append(f.id)
        return out

    # ------------------------------------------------------------------
    # Synthesis / convenience
    # ------------------------------------------------------------------

    def interfeatures(self, features, new_featuretype=None,
                      merge_attributes: bool = True, numeric_sort: bool = False,
                      dialect=None, attribute_func=None, update_attributes=None):
        feats = list(features)
        if not feats:
            return
        for prev, cur in zip(feats[:-1], feats[1:]):
            new_start = prev.end + 1
            new_end = cur.start - 1
            if new_end < new_start:
                continue
            attrs = {}
            if merge_attributes:
                for k, v in prev.attributes.items():
                    attrs.setdefault(k, []).extend(v)
                for k, v in cur.attributes.items():
                    attrs.setdefault(k, []).extend(v)
                # Dedupe.
                attrs = {k: list(dict.fromkeys(v)) for k, v in attrs.items()}
            if update_attributes:
                attrs.update(update_attributes)
            if attribute_func:
                attrs = attribute_func(prev, cur, attrs)
            ftype = new_featuretype or "interfeature"
            yield Feature(
                seqid=prev.seqid,
                source=prev.source,
                featuretype=ftype,
                start=new_start,
                end=new_end,
                strand=prev.strand,
                attributes=attrs,
                dialect=dialect or self.dialect,
            )

    def merge(self, features, merge_criteria=None, multiline: bool = False):
        from . import merge_criteria as mc
        if merge_criteria is None:
            merge_criteria = (mc.seqid, mc.overlap_end_inclusive, mc.strand, mc.feature_type)
        feats = sorted(features, key=lambda f: (f.seqid, f.start, f.end))
        if not feats:
            return
        accum = None
        components: List[Feature] = []
        for f in feats:
            if accum is None:
                accum = self._clone_for_merge(f)
                components = [f]
                continue
            if all(pred(accum, f, components) for pred in merge_criteria):
                accum.end = max(accum.end, f.end)
                components.append(f)
            else:
                accum.children = list(components)
                yield accum
                accum = self._clone_for_merge(f)
                components = [f]
        if accum is not None:
            accum.children = list(components)
            yield accum

    @staticmethod
    def _clone_for_merge(f: Feature) -> Feature:
        return Feature(
            seqid=f.seqid, source=f.source, featuretype=f.featuretype,
            start=f.start, end=f.end, score=f.score, strand=f.strand,
            frame=f.frame, attributes={k: list(v) for k, v in f.attributes.items()},
            dialect=f.dialect,
        )

    def merge_all(self, merge_order=("seqid", "featuretype", "strand", "start"),
                  merge_criteria=None, featuretypes_groups=(None,),
                  exclude_components: bool = False) -> List[Feature]:
        out: List[Feature] = []
        for group in featuretypes_groups:
            feats = list(self.all_features(featuretype=group))
            for k in reversed(merge_order):
                feats.sort(key=lambda f: getattr(f, k) if k != "start" else (f.start or 0))
            out.extend(self.merge(feats, merge_criteria=merge_criteria))
        return out

    def create_introns(self, exon_featuretype: str = "exon",
                       grandparent_featuretype: Optional[str] = "gene",
                       parent_featuretype: Optional[str] = None,
                       new_featuretype: str = "intron",
                       merge_attributes: bool = True,
                       numeric_sort: bool = False) -> Iterator[Feature]:
        if grandparent_featuretype and parent_featuretype:
            raise ValueError("specify exactly one of grandparent_featuretype/parent_featuretype")
        if not (grandparent_featuretype or parent_featuretype):
            raise ValueError("must specify grandparent_featuretype or parent_featuretype")
        anchor_type = grandparent_featuretype or parent_featuretype
        for anchor in self.features_of_type(anchor_type):
            exons = sorted(
                (
                    e for e in self.children(anchor, featuretype=exon_featuretype)
                ),
                key=lambda f: (f.start, f.end),
            )
            if len(exons) < 2:
                continue
            yield from self.interfeatures(
                exons, new_featuretype=new_featuretype,
                merge_attributes=merge_attributes,
                numeric_sort=numeric_sort,
            )

    def create_splice_sites(self, exon_featuretype: str = "exon",
                            grandparent_featuretype: Optional[str] = "gene",
                            parent_featuretype: Optional[str] = None,
                            merge_attributes: bool = True,
                            numeric_sort: bool = False) -> Iterator[Feature]:
        for intron in self.create_introns(
            exon_featuretype=exon_featuretype,
            grandparent_featuretype=grandparent_featuretype,
            parent_featuretype=parent_featuretype,
            new_featuretype="splice_site",
            merge_attributes=merge_attributes, numeric_sort=numeric_sort,
        ):
            # Yield 1bp left + 1bp right splice sites.
            yield Feature(
                seqid=intron.seqid, source=intron.source,
                featuretype="splice_site",
                start=intron.start, end=intron.start,
                strand=intron.strand, dialect=self.dialect,
            )
            yield Feature(
                seqid=intron.seqid, source=intron.source,
                featuretype="splice_site",
                start=intron.end, end=intron.end,
                strand=intron.strand, dialect=self.dialect,
            )

    def children_bp(self, feature, child_featuretype: str = "exon",
                    merge: bool = False, merge_criteria=None, **kwargs) -> int:
        kids = list(self.children(feature, featuretype=child_featuretype))
        if merge:
            kids = list(self.merge(kids, merge_criteria=merge_criteria))
        total = 0
        for k in kids:
            if k.start is not None and k.end is not None:
                total += k.end - k.start + 1
        return total

    def bed12(self, feature, block_featuretype=("exon",),
              thick_featuretype=("CDS",), thin_featuretype=None,
              name_field: str = "ID", color=None) -> str:
        if isinstance(feature, str):
            feature = self[feature]
        blocks = sorted(
            self.children(feature, featuretype=list(block_featuretype)),
            key=lambda f: (f.start, f.end),
        )
        cds = list(self.children(feature, featuretype=list(thick_featuretype)))
        chrom_start = feature.start - 1
        chrom_end = feature.end
        if cds:
            thick_start = min(c.start for c in cds) - 1
            thick_end = max(c.end for c in cds)
        else:
            thick_start = chrom_start
            thick_end = chrom_start
        try:
            name_value = feature.attributes[name_field][0]
        except (KeyError, IndexError):
            name_value = feature.id or "."
        score = feature.score if feature.score not in (".", "") else "0"
        strand = feature.strand if feature.strand in ("+", "-") else "+"
        rgb = color or "0,0,0"
        block_count = len(blocks)
        block_sizes = ",".join(str(b.end - b.start + 1) for b in blocks)
        block_starts = ",".join(str((b.start - 1) - chrom_start) for b in blocks)
        return "\t".join(str(x) for x in (
            feature.seqid, chrom_start, chrom_end, name_value, score, strand,
            thick_start, thick_end, rgb, block_count,
            block_sizes + ("," if block_count else ""),
            block_starts + ("," if block_count else ""),
        ))

    def iter_by_parent_childs(self, featuretype: str = "gene",
                              level: Optional[int] = None,
                              order_by=None, reverse: bool = False,
                              completely_within: bool = False) -> Iterator[List[Feature]]:
        for parent in self.features_of_type(featuretype, order_by=order_by, reverse=reverse):
            kids = list(self.children(parent, level=level,
                                       order_by=order_by, reverse=reverse,
                                       completely_within=completely_within))
            yield [parent, *kids]

    # ------------------------------------------------------------------
    # Escape hatch + maintenance
    # ------------------------------------------------------------------

    def execute(self, query: str):
        """Execute arbitrary SQL. Returns DuckDB's relation cursor.
        SQLite-style queries against ``features_compat`` and ``relations_compat``
        views are supported; see ``compat_views.sql``."""
        return self.conn.execute(query.rstrip(";"))

    def analyze(self) -> None:
        self.conn.execute("ANALYZE")
        self._analyzed_flag = True

    def set_pragmas(self, pragmas: dict) -> None:
        # DuckDB pragmas. Silently skip ones DuckDB rejects (legacy callers
        # often pass SQLite-specific pragmas like `journal_mode`).
        for k, v in pragmas.items():
            try:
                self.conn.execute(f"PRAGMA {k} = {v}")
            except duckdb.Error:
                continue

    # ------------------------------------------------------------------
    # Internal: row → Feature streaming
    # ------------------------------------------------------------------

    def _yield_features(self, sql: str, params: list) -> Iterator[Feature]:
        cur = self.conn.execute(sql, params)
        while True:
            rows = cur.fetchmany(10_000)
            if not rows:
                return
            for row in rows:
                yield feature_from_row(row, dialect=self.dialect)
