# ---------------------------------------------------------------------------
# Author: Kuan-Hao Chao <kuanhao.chao@gmail.com>
# ---------------------------------------------------------------------------
"""Phase 4 ingestion engine.

Streams the Rust parser's output through PyArrow record batches into DuckDB,
then runs a fixed sequence of set-based SQL passes for normalization, GTF
synthesis, transitive closure, and indexing. No per-feature Python loops in
the hot path — everything that scales with feature count goes through Arrow
or pure SQL.
"""

from __future__ import annotations

import os
import re
import warnings
from dataclasses import dataclass
from typing import Iterable, Iterator, List, Optional, Tuple

import duckdb
import pyarrow as pa

from . import parser as _parser
from .feature import ParsedFeature
from .schema import (
    DDL,
    POST_LOAD_INDEXES,
    EDGES_FROM_PARENT,
    EDGES_FROM_GTF,
    GTF_SYNTHESIZE_TRANSCRIPTS,
    GTF_SYNTHESIZE_TRANSCRIPT_ATTRS,
    GTF_PROPAGATE_GENE_ID,
    GTF_SYNTHESIZE_GENES,
    CLOSURE_RECURSIVE_CTE,
    COMPAT_VIEWS_SQL,
    SCHEMA_VERSION,
)


DEFAULT_BATCH_SIZE = 50_000
DEFAULT_MAX_DEPTH = 8

# DuckDB 1.5.3 and 1.5.4 can corrupt the unified vector representation while
# appending a shredded GEOMETRY column across row groups (DuckDB #23737).
# LiftOn's bounding boxes are an optional query accelerator, so affected
# versions use the existing B-tree path until the upstream fix is available.
_GEOMETRY_APPEND_BUG_VERSIONS = {(1, 5, 3), (1, 5, 4)}
_GEOMETRY_APPEND_ERROR = (
    "Expected unified vector format of type LIST, but found type VARCHAR"
)


@dataclass
class IngestStats:
    """Reported back to the caller for benchmarking and tests."""
    n_features_raw: int = 0
    n_features_synthetic_transcripts: int = 0
    n_features_synthetic_genes: int = 0
    n_attributes: int = 0
    n_edges: int = 0
    n_closure_rows: int = 0
    rtree_built: bool = False
    fmt: str = "gff3"
    dialect: dict = None  # type: ignore[assignment]
    directives: List[str] = None  # type: ignore[assignment]


class GFF3ChunkDecoder:
    """Incrementally decode byte chunks into validated GFF3 records.

    Only an incomplete trailing line is retained between calls.  Oversized
    input chunks are internally split at ``max_chunk_bytes`` (1 MiB by
    default), and pathological unbroken lines fail before growing without a
    bound.
    """

    def __init__(
        self, *, max_chunk_bytes: int = 1 << 20,
        max_line_bytes: int = 8 << 20,
    ):
        self.max_chunk_bytes = max(1, int(max_chunk_bytes))
        self.max_line_bytes = max(self.max_chunk_bytes, int(max_line_bytes))
        self._pending = bytearray()
        self._line_no = 0
        self.byte_count = 0
        self.feature_count = 0
        self.max_pending_bytes = 0
        self.directives: List[str] = []

    def feed(self, chunk: bytes):
        if not isinstance(chunk, (bytes, bytearray, memoryview)):
            raise TypeError("GFF3ChunkDecoder.feed expects a bytes-like chunk")
        data = chunk if isinstance(chunk, bytes) else bytes(chunk)
        self.byte_count += len(data)
        for offset in range(0, len(data), self.max_chunk_bytes):
            piece_end = min(len(data), offset + self.max_chunk_bytes)
            cursor = offset
            while True:
                newline = data.find(b"\n", cursor, piece_end)
                if newline < 0:
                    self._pending.extend(data[cursor:piece_end])
                    self.max_pending_bytes = max(
                        self.max_pending_bytes, len(self._pending),
                    )
                    break
                if self._pending:
                    self._pending.extend(data[cursor:newline])
                    raw = bytes(self._pending)
                    self._pending.clear()
                else:
                    raw = data[cursor:newline]
                record = self._decode_line(raw.rstrip(b"\r"))
                if record is not None:
                    yield record
                cursor = newline + 1
            if len(self._pending) > self.max_line_bytes:
                raise ValueError(
                    "GFF3 record exceeds the configured maximum line size "
                    f"({self.max_line_bytes} bytes)"
                )

    def finish(self):
        if self._pending:
            raw = bytes(self._pending).rstrip(b"\r")
            self._pending.clear()
            record = self._decode_line(raw)
            if record is not None:
                yield record

    @property
    def pending_bytes(self) -> int:
        return len(self._pending)

    def _decode_line(self, raw: bytes):
        self._line_no += 1
        if not raw:
            return None
        line = raw.decode("utf-8", errors="strict")
        if line.startswith("#"):
            if line.startswith("##") or line.startswith("#!"):
                self.directives.append(line)
            return None
        from ._pyfallback.parser import _parse_line_into_feature
        record = _parse_line_into_feature(line, self._line_no)
        self.feature_count += 1
        return record


# ---------------------------------------------------------------------------
# Arrow batch builder. Accumulates ParsedFeature output until a fixed row
# count, then yields one PyArrow Table per category (features, attributes).
# Edges and directives are derived later by SQL.
# ---------------------------------------------------------------------------

class _ArrowBatchBuilder:
    """Accumulates parsed features into column-oriented Python lists, then
    produces PyArrow tables on flush. We deliberately keep the schema explicit
    so DuckDB sees the right column types (no INFER passes).

    Phase 19: the builder also stamps each row's ``seqid_y`` value during
    ``append()`` using a shared ``seqid_to_y`` dict (lazy band assignment in
    encounter order). This eliminates two full-table ``UPDATE`` passes that
    used to dominate ingest wall time on real GFF3 corpora.
    """

    FEATURES_SCHEMA = pa.schema([
        ("id",              pa.string()),
        ("seqid",           pa.string()),
        ("source",          pa.string()),
        ("featuretype",     pa.string()),
        ("start",           pa.int64()),
        ("end",             pa.int64()),
        ("score",           pa.string()),
        ("strand",          pa.string()),
        ("frame",           pa.string()),
        ("attributes_blob", pa.binary()),
        ("extra_blob",      pa.binary()),
        ("file_order",      pa.int64()),
        ("is_synthetic",    pa.bool_()),
        ("seqid_y",         pa.int64()),
    ])

    ATTRIBUTES_SCHEMA = pa.schema([
        ("feature_id", pa.string()),
        ("key",        pa.string()),
        ("value",      pa.string()),
        ("idx",        pa.int16()),
    ])

    def __init__(self, seqid_to_y: dict, has_spatial: bool = False):
        # Shared across all batches so seqid_y assignment is stable for the
        # whole file. Caller owns the dict; we mutate it in place.
        self._seqid_to_y = seqid_to_y
        self._has_spatial = has_spatial
        self._reset()

    def _reset(self):
        # Feature columns
        self.f_id: list = []
        self.f_seqid: list = []
        self.f_source: list = []
        self.f_type: list = []
        self.f_start: list = []
        self.f_end: list = []
        self.f_score: list = []
        self.f_strand: list = []
        self.f_frame: list = []
        self.f_blob: list = []
        self.f_extra: list = []
        self.f_order: list = []
        self.f_synth: list = []
        self.f_seqid_y: list = []
        # Attribute columns
        self.a_fid: list = []
        self.a_key: list = []
        self.a_val: list = []
        self.a_idx: list = []

    def append(self, feat_id: str, feat: ParsedFeature, file_order: int):
        seqid = feat.seqid
        # Lazy y-band assignment: each new seqid gets the next slot.
        # `dict.get` + assignment is faster than `setdefault` here because
        # we hit the cache path on >99 % of rows in real data.
        y = self._seqid_to_y.get(seqid)
        if y is None:
            y = len(self._seqid_to_y) * SEQID_Y_BAND
            self._seqid_to_y[seqid] = y
        self.f_id.append(feat_id)
        self.f_seqid.append(seqid)
        self.f_source.append(feat.source)
        self.f_type.append(feat.featuretype)
        self.f_start.append(feat.start if feat.start is not None else 0)
        self.f_end.append(feat.end if feat.end is not None else 0)
        self.f_score.append(feat.score)
        self.f_strand.append(feat.strand)
        self.f_frame.append(feat.frame)
        self.f_blob.append(feat.attributes_blob)
        self.f_extra.append(("\t".join(feat.extra)).encode("utf-8") if feat.extra else b"")
        self.f_order.append(file_order)
        self.f_synth.append(False)
        self.f_seqid_y.append(y)
        for k, v, idx in feat.attributes_pairs:
            self.a_fid.append(feat_id)
            self.a_key.append(k)
            self.a_val.append(v)
            self.a_idx.append(idx)

    def __len__(self) -> int:
        return len(self.f_id)

    def features_table(self) -> pa.Table:
        return pa.table({
            "id": self.f_id,
            "seqid": self.f_seqid,
            "source": self.f_source,
            "featuretype": self.f_type,
            "start": self.f_start,
            "end": self.f_end,
            "score": self.f_score,
            "strand": self.f_strand,
            "frame": self.f_frame,
            "attributes_blob": self.f_blob,
            "extra_blob": self.f_extra,
            "file_order": self.f_order,
            "is_synthetic": self.f_synth,
            "seqid_y": self.f_seqid_y,
        }, schema=self.FEATURES_SCHEMA)

    def attributes_table(self) -> pa.Table:
        return pa.table({
            "feature_id": self.a_fid,
            "key": self.a_key,
            "value": self.a_val,
            "idx": self.a_idx,
        }, schema=self.ATTRIBUTES_SCHEMA)

    def flush_into(self, con: duckdb.DuckDBPyConnection):
        if not self.f_id:
            return
        feats = self.features_table()
        attrs = self.attributes_table()
        # Register and INSERT ... SELECT — DuckDB's fastest Arrow path. When
        # the spatial extension is loaded we ALSO compute `bbox` inline so
        # the R-tree build at the end of ingest is a single CREATE INDEX
        # (no UPDATE pass over the table).
        registered = []
        try:
            con.register("__staging_features", feats)
            registered.append("__staging_features")
            con.register("__staging_attributes", attrs)
            registered.append("__staging_attributes")
            if self._has_spatial:
                con.execute(
                    "INSERT INTO features ("
                    "id, seqid, source, featuretype, start, \"end\", "
                    "score, strand, frame, attributes_blob, extra_blob, "
                    "file_order, is_synthetic, seqid_y, bbox"
                    ") SELECT id, seqid, source, featuretype, start, \"end\", "
                    "score, strand, frame, attributes_blob, extra_blob, "
                    "file_order, is_synthetic, seqid_y, "
                    "ST_MakeEnvelope(start, seqid_y, \"end\", seqid_y + 1) "
                    "FROM __staging_features"
                )
            else:
                con.execute(
                    "INSERT INTO features ("
                    "id, seqid, source, featuretype, start, \"end\", "
                    "score, strand, frame, attributes_blob, extra_blob, "
                    "file_order, is_synthetic, seqid_y"
                    ") SELECT * FROM __staging_features"
                )
            con.execute(
                "INSERT INTO attributes SELECT * FROM __staging_attributes"
            )
        except Exception:
            # DuckDB InternalException invalidates the connection. Release
            # the whole connection before create_db retries on a fresh,
            # geometry-free database. Do not issue more SQL against an
            # invalidated connection; close() releases Arrow registrations.
            try:
                con.close()
            except Exception:
                pass
            raise
        for name in registered:
            con.unregister(name)
        self._reset()


# ---------------------------------------------------------------------------
# ID resolution. Pulled out so the bulk loop has zero branches that hit Python
# attribute parsing twice.
# ---------------------------------------------------------------------------

def _derive_id(feat: ParsedFeature, dialect_fmt: str, autoincrement: dict) -> str:
    """Compute the row's primary key. GFF3 prefers `ID=`; GTF synthesizes
    `<featuretype>_<n>` so leaf rows still get unique IDs (gene/transcript
    IDs are filled in by the synthesis pass)."""
    if dialect_fmt == "gff3":
        for k, v, _idx in feat.attributes_pairs:
            if k == "ID":
                return v
    # Fall back: synthesize an auto-id by featuretype.
    n = autoincrement.get(feat.featuretype, 0) + 1
    autoincrement[feat.featuretype] = n
    return f"{feat.featuretype}_{n}"


# ---------------------------------------------------------------------------
# Public entry point.
# ---------------------------------------------------------------------------

def _from_record_stream_impl(
    records: Iterable[ParsedFeature],
    dbfn: str = ":memory:",
    *,
    dialect: Optional[dict] = None,
    directives: Optional[List[str]] = None,
    force: bool = False,
    batch_size: int = DEFAULT_BATCH_SIZE,
    max_depth: int = DEFAULT_MAX_DEPTH,
    disable_infer_genes: bool = False,
    disable_infer_transcripts: bool = False,
    gtf_subfeature: str = "exon",
    engine: Optional[str] = "auto",
    build_rtree: bool = True,
    _connection_out: Optional[List[duckdb.DuckDBPyConnection]] = None,
) -> Tuple[duckdb.DuckDBPyConnection, IngestStats]:
    """Ingest an incremental ``ParsedFeature`` stream into DuckDB.

    The iterable may expose parser-style ``dialect()`` and ``directives()``
    methods; explicit arguments take precedence.  Only the active Arrow batch
    (default 50,000 records) is retained in Python memory.
    """
    if dbfn != ":memory:":
        if os.path.exists(dbfn) and not force:
            raise ValueError(
                f"{dbfn} already exists. Pass force=True to overwrite."
            )
        if os.path.exists(dbfn) and force:
            os.unlink(dbfn)

    con = duckdb.connect(dbfn)
    if _connection_out is not None:
        _connection_out.append(con)
    _apply_pragmas(con)
    con.execute(DDL)

    # Phase 19: load the spatial extension UPFRONT (was: lazy after bulk
    # load). This lets us widen the `features` schema to include `bbox`
    # and stamp the R-tree envelope inline during the Arrow batch INSERT,
    # eliminating two full-table UPDATE passes that used to dominate
    # ingest wall time.
    spatial_requested = build_rtree and not _rtree_disabled_by_env()
    affected_duckdb = _duckdb_needs_geometry_shredding_workaround()
    if spatial_requested and affected_duckdb:
        warnings.warn(
            f"DuckDB {duckdb.__version__} has a known GEOMETRY append bug; "
            "LiftOn is disabling its optional R-tree for this ingest. Use "
            "DuckDB 1.5.2 or a fixed release to restore R-tree acceleration.",
            RuntimeWarning,
            stacklevel=2,
        )
    has_spatial = (
        spatial_requested
        and not affected_duckdb
        and _try_load_spatial(con)
    )
    if has_spatial:
        con.execute("ALTER TABLE features ADD COLUMN IF NOT EXISTS bbox GEOMETRY")

    seqid_to_y: dict = {}
    builder = _ArrowBatchBuilder(seqid_to_y, has_spatial=has_spatial)
    autoinc: dict = {}
    # NCBI RefSeq emits multiple GFF3 rows that share an `ID=cds-…` (a CDS
    # is "split" across exon-segments). Our schema has `id` as a primary
    # key, so we mimic gffutils' `merge_strategy="create_unique"`: track
    # how many times we've seen each base id and append `__N` (N >= 2)
    # when needed. The first occurrence keeps the bare id.
    id_counts: dict = {}
    duplicate_pairs: list = []  # (base_id, new_id)
    file_order = 0
    n_raw = 0
    # Resolve the dialect format ONCE — calling `it.dialect()` per record
    # is a Rust↔Python boundary cost we shouldn't pay 5 M times. The
    # parser commits to a dialect during the peek phase, so the value is
    # stable from the first yielded record onward.
    _fmt_cache: Optional[str] = None

    for feat in records:
        file_order += 1
        n_raw += 1
        if _fmt_cache is None:
            _fmt_cache = (
                (dialect or {}).get("fmt", "gff3")
                if dialect is not None else _dialect_fmt_safe(records)
            )
        fid = _derive_id(feat, _fmt_cache, autoinc)
        seen = id_counts.get(fid, 0)
        if seen:
            new_fid = f"{fid}__{seen + 1}"
            duplicate_pairs.append((fid, new_fid))
            id_counts[fid] = seen + 1
            fid = new_fid
        else:
            id_counts[fid] = 1
        builder.append(fid, feat, file_order)
        if len(builder) >= batch_size:
            builder.flush_into(con)
    builder.flush_into(con)

    # Record duplicate-id remappings (informational; the schema already has
    # this table — Phase 5).
    if duplicate_pairs:
        dup_tbl = pa.table({
            "original_id": [b for b, _ in duplicate_pairs],
            "new_id":      [n for _, n in duplicate_pairs],
        })
        con.register("__staging_dups", dup_tbl)
        con.execute(
            "INSERT INTO duplicates (original_id, new_id) "
            "SELECT original_id, new_id FROM __staging_dups"
        )
        con.unregister("__staging_dups")

    resolved_dialect = dialect
    if resolved_dialect is None:
        getter = getattr(records, "dialect", None)
        resolved_dialect = getter() if callable(getter) else {"fmt": "gff3"}
    resolved_directives = directives
    if resolved_directives is None:
        getter = getattr(records, "directives", None)
        resolved_directives = list(getter()) if callable(getter) else []
    else:
        resolved_directives = list(resolved_directives)
    fmt = (resolved_dialect or {}).get("fmt", "gff3")

    # Persist directives in one shot (set-based).
    if resolved_directives:
        dir_table = pa.table({"directive": resolved_directives})
        con.register("__staging_directives", dir_table)
        con.execute("INSERT INTO directives (directive) SELECT directive FROM __staging_directives")
        con.unregister("__staging_directives")

    # Set-based normalization passes.
    n_synth_t = 0
    n_synth_g = 0

    if fmt == "gtf":
        if not disable_infer_transcripts:
            n_synth_t = _synthesize_transcripts(con, gtf_subfeature)
        if not disable_infer_genes:
            n_synth_g = _synthesize_genes(con, gtf_subfeature)
        con.execute(EDGES_FROM_GTF)
        # GTF synthesis inserts new rows without seqid_y / bbox set. Patch
        # them up in a single targeted UPDATE (touches only synthesized
        # rows; ~4-9 % of features at GENCODE scale).
        if has_spatial:
            con.execute(
                "UPDATE features "
                "SET seqid_y = m.seqid_y, "
                "    bbox = ST_MakeEnvelope("
                "        features.start, m.seqid_y, "
                "        features.\"end\", m.seqid_y + 1) "
                "FROM seqid_map m "
                "WHERE features.seqid = m.seqid AND features.seqid_y IS NULL"
            )
        else:
            con.execute(
                "UPDATE features SET seqid_y = m.seqid_y "
                "FROM seqid_map m "
                "WHERE features.seqid = m.seqid AND features.seqid_y IS NULL"
            )
    else:
        con.execute(EDGES_FROM_PARENT)

    # Closure via recursive CTE.
    con.execute(CLOSURE_RECURSIVE_CTE, [max_depth])

    # Indexes — only after all data is materialized.
    con.execute(POST_LOAD_INDEXES)

    # Optional R-tree. Phase 19: when spatial is loaded, this is now a
    # single CREATE INDEX over the bbox column we already populated
    # inline during the Arrow batch INSERTs (no UPDATE pass).
    rtree_built = False
    if has_spatial:
        rtree_built = _finalize_rtree(con, seqid_to_y)

    # SQLite-compat views (must run after closure has been populated).
    con.execute(COMPAT_VIEWS_SQL)

    # Stats.
    n_attributes = con.execute("SELECT COUNT(*) FROM attributes").fetchone()[0]
    n_edges = con.execute("SELECT COUNT(*) FROM edges").fetchone()[0]
    n_closure = con.execute("SELECT COUNT(*) FROM closure").fetchone()[0]

    # Meta — record dialect, fmt, and the rtree availability so a re-opened
    # DB can route queries correctly without probing.
    _write_meta(
        con, resolved_dialect, fmt,
        rtree_built=rtree_built, max_depth=max_depth,
    )

    return con, IngestStats(
        n_features_raw=n_raw,
        n_features_synthetic_transcripts=n_synth_t,
        n_features_synthetic_genes=n_synth_g,
        n_attributes=n_attributes,
        n_edges=n_edges,
        n_closure_rows=n_closure,
        rtree_built=rtree_built,
        fmt=fmt,
        dialect=resolved_dialect,
        directives=resolved_directives,
    )


def from_record_stream(
    records: Iterable[ParsedFeature],
    dbfn: str = ":memory:",
    *,
    dialect: Optional[dict] = None,
    directives: Optional[List[str]] = None,
    force: bool = False,
    batch_size: int = DEFAULT_BATCH_SIZE,
    max_depth: int = DEFAULT_MAX_DEPTH,
    disable_infer_genes: bool = False,
    disable_infer_transcripts: bool = False,
    gtf_subfeature: str = "exon",
    engine: Optional[str] = "auto",
    build_rtree: bool = True,
) -> Tuple[duckdb.DuckDBPyConnection, IngestStats]:
    """Ingest records and explicitly close the database on every failure."""

    connection: List[duckdb.DuckDBPyConnection] = []
    try:
        return _from_record_stream_impl(
            records,
            dbfn=dbfn,
            dialect=dialect,
            directives=directives,
            force=force,
            batch_size=batch_size,
            max_depth=max_depth,
            disable_infer_genes=disable_infer_genes,
            disable_infer_transcripts=disable_infer_transcripts,
            gtf_subfeature=gtf_subfeature,
            engine=engine,
            build_rtree=build_rtree,
            _connection_out=connection,
        )
    except BaseException:
        if connection:
            try:
                connection[0].close()
            except Exception:
                pass
        raise


def from_file(
    path: str,
    dbfn: str = ":memory:",
    *,
    force: bool = False,
    batch_size: int = DEFAULT_BATCH_SIZE,
    max_depth: int = DEFAULT_MAX_DEPTH,
    disable_infer_genes: bool = False,
    disable_infer_transcripts: bool = False,
    gtf_subfeature: str = "exon",
    engine: Optional[str] = "auto",
    build_rtree: bool = True,
) -> Tuple[duckdb.DuckDBPyConnection, IngestStats]:
    """Parse a file and delegate to the incremental record ingest path."""
    records = _parser.parse_gff(path, engine=engine)
    try:
        return from_record_stream(
            records, dbfn=dbfn, force=force, batch_size=batch_size,
            max_depth=max_depth,
            disable_infer_genes=disable_infer_genes,
            disable_infer_transcripts=disable_infer_transcripts,
            gtf_subfeature=gtf_subfeature, build_rtree=build_rtree,
        )
    finally:
        records.close()


# ---------------------------------------------------------------------------
# Helpers.
# ---------------------------------------------------------------------------

def _dialect_fmt_safe(it) -> str:
    """The Rust iterator commits to a dialect after the peek phase. The
    Python fallback only sets it after the first record yields. Both are
    populated by the time we land in this function on the first feature."""
    try:
        d = it.dialect()
        if d:
            return d.get("fmt", "gff3")
    except Exception:
        pass
    return "gff3"


def _apply_pragmas(con: duckdb.DuckDBPyConnection):
    # DuckDB's defaults are excellent; we only nudge memory & threads. Anything
    # missing here is left to the caller via DUCKDB_THREADS env var.
    threads = os.environ.get("GFFUTILS2_THREADS")
    if threads:
        con.execute(f"PRAGMA threads = {int(threads)}")
    # Suppress the interactive progress bar — it floods stderr in batch and
    # subprocess scenarios and offers no value for benchmarking or scripting.
    try:
        con.execute("PRAGMA disable_progress_bar")
    except duckdb.Error:
        pass


def _synthesize_transcripts(con, subfeature: str) -> int:
    """Run the GROUP BY transcript synthesis. Returns rows inserted."""
    before = con.execute(
        "SELECT COUNT(*) FROM features WHERE featuretype = 'transcript'"
    ).fetchone()[0]
    con.execute(GTF_SYNTHESIZE_TRANSCRIPTS, [subfeature])
    after = con.execute(
        "SELECT COUNT(*) FROM features WHERE featuretype = 'transcript'"
    ).fetchone()[0]
    n = after - before
    # Make the synthesized transcripts visible to subsequent passes by giving
    # them a self-referential transcript_id attribute, then propagate the
    # gene_id from the children. The propagation joins attributes directly,
    # so no temporary edge inserts are needed (duplicate-edge-free).
    con.execute(GTF_SYNTHESIZE_TRANSCRIPT_ATTRS)
    con.execute(GTF_PROPAGATE_GENE_ID)
    return n


def _synthesize_genes(con, subfeature: str) -> int:
    before = con.execute(
        "SELECT COUNT(*) FROM features WHERE featuretype = 'gene'"
    ).fetchone()[0]
    con.execute(GTF_SYNTHESIZE_GENES, [subfeature])
    after = con.execute(
        "SELECT COUNT(*) FROM features WHERE featuretype = 'gene'"
    ).fetchone()[0]
    n = after - before
    # Mirror the synthesized gene_id into attributes so downstream queries
    # treat synthesized genes the same as authored ones.
    con.execute(
        """
        INSERT INTO attributes (feature_id, key, value, idx)
        SELECT f.id, 'gene_id', f.id, 0
        FROM features f
        WHERE f.featuretype = 'gene' AND f.is_synthetic = TRUE
        """
    )
    return n


SEQID_Y_BAND = 1_000_000  # gap between adjacent seqids' y-bands.


def _duckdb_needs_geometry_shredding_workaround(
    version: Optional[str] = None,
) -> bool:
    """Return whether ``version`` has DuckDB's row-group GEOMETRY bug.

    Development and vendor suffixes are intentionally accepted, for example
    ``1.5.4.dev12`` and ``v1.5.3-custom``.
    """
    raw_version = duckdb.__version__ if version is None else version
    match = re.search(r"(\d+)\.(\d+)\.(\d+)", raw_version)
    if match is None:
        return False
    parsed = tuple(int(part) for part in match.groups())
    return parsed in _GEOMETRY_APPEND_BUG_VERSIONS


def is_geometry_append_error(exc: BaseException) -> bool:
    """Recognize the upstream InternalException safe to retry without bbox."""
    return (
        isinstance(exc, duckdb.Error)
        and _GEOMETRY_APPEND_ERROR in str(exc)
    )


def _rtree_disabled_by_env() -> bool:
    """Return whether the process explicitly requested the B-tree path.

    ``LIFTON_DISABLE_RTREE`` is the supported user control. The older
    ``GFFBASE_TEST_DISABLE_RTREE`` name remains available for CI and callers
    that used it before the public alias was added.
    """
    truthy = ("1", "true", "yes", "on")
    return any(
        os.environ.get(name, "").strip().lower() in truthy
        for name in (
            "LIFTON_DISABLE_RTREE",
            "GFFBASE_TEST_DISABLE_RTREE",
        )
    )


def _try_load_spatial(con: duckdb.DuckDBPyConnection) -> bool:
    """Attempt to install + load the DuckDB spatial extension. Returns
    True iff it's now usable on this connection."""
    try:
        con.execute("INSTALL spatial")
        con.execute("LOAD spatial")
        return True
    except duckdb.Error:
        return False


def _finalize_rtree(con: duckdb.DuckDBPyConnection, seqid_to_y: dict) -> bool:
    """Persist the seqid_to_y dict into the `seqid_map` side table and
    create the R-tree index over the (already-populated) `bbox` column.

    Phase 19: this is the entire R-tree build — no UPDATE passes. The
    `bbox` column was filled in inline by ``_ArrowBatchBuilder.flush_into``
    using the per-row seqid_y stamped by the builder.
    """
    try:
        if seqid_to_y:
            seqid_rows = list(seqid_to_y.items())
            # Stable ordering (encounter order in the file) — preserves the
            # invariant that the first seqid sees seqid_y == 0.
            con.execute("DELETE FROM seqid_map")
            con.executemany(
                "INSERT INTO seqid_map(seqid, seqid_y) VALUES (?, ?)",
                seqid_rows,
            )
        con.execute(
            "CREATE INDEX IF NOT EXISTS features_rtree "
            "ON features USING RTREE (bbox)"
        )
        return True
    except duckdb.Error:
        return False


def _write_meta(con: duckdb.DuckDBPyConnection, dialect: dict, fmt: str,
                *, rtree_built: bool = False, max_depth: int = DEFAULT_MAX_DEPTH):
    import json
    # Closure max depth: used by FeatureDB's relational dispatcher (Phase 7) to
    # pick the cache vs. dynamic CTE without a per-call query.
    row = con.execute("SELECT MAX(depth) FROM closure").fetchone()
    closure_max_depth = int(row[0]) if row and row[0] is not None else 0
    rows = [
        ("schema_version", SCHEMA_VERSION),
        ("dialect", json.dumps(dialect or {})),
        ("fmt", fmt),
        ("rtree_built", "true" if rtree_built else "false"),
        ("max_depth", str(int(max_depth))),
        ("closure_max_depth", str(closure_max_depth)),
    ]
    con.executemany("INSERT OR REPLACE INTO meta(key, value) VALUES (?, ?)", rows)
