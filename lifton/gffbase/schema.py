# ---------------------------------------------------------------------------
# Author: Kuan-Hao Chao <kuanhao.chao@gmail.com>
# ---------------------------------------------------------------------------
"""DDL for the gffbase DuckDB schema (Phase 2 §4.1) and post-load index/synthesis SQL.

Centralized so the ingestion engine and tests share one source of truth.
"""

from __future__ import annotations


# Seven user-facing tables plus the bbox helper view if the spatial extension
# is available. Schema version is recorded in `meta`.
SCHEMA_VERSION = "1"

DDL = """
CREATE TABLE IF NOT EXISTS features (
    id              VARCHAR PRIMARY KEY,
    seqid           VARCHAR NOT NULL,
    source          VARCHAR,
    featuretype     VARCHAR NOT NULL,
    start           BIGINT  NOT NULL,
    "end"           BIGINT  NOT NULL,
    score           VARCHAR,
    strand          VARCHAR,
    frame           VARCHAR,
    attributes_blob BLOB,
    extra_blob      BLOB,
    file_order      BIGINT,
    is_synthetic    BOOLEAN DEFAULT FALSE,
    -- y-band populated post-load: each distinct seqid gets a unique band so
    -- the R-tree split heuristics segregate chromosomes (Phase 7 fix).
    seqid_y         BIGINT
);

-- One row per distinct seqid mapping to its R-tree y-band. Populated during
-- the R-tree build step; consulted at query time to translate `region()`
-- bounds into the same band space used by the index.
CREATE TABLE IF NOT EXISTS seqid_map (
    seqid  VARCHAR PRIMARY KEY,
    seqid_y BIGINT NOT NULL
);

CREATE TABLE IF NOT EXISTS attributes (
    feature_id      VARCHAR NOT NULL,
    key             VARCHAR NOT NULL,
    value           VARCHAR NOT NULL,
    idx             SMALLINT NOT NULL DEFAULT 0
);

CREATE TABLE IF NOT EXISTS edges (
    parent          VARCHAR NOT NULL,
    child           VARCHAR NOT NULL
);

CREATE TABLE IF NOT EXISTS closure (
    ancestor        VARCHAR NOT NULL,
    descendant      VARCHAR NOT NULL,
    depth           SMALLINT NOT NULL
);

CREATE TABLE IF NOT EXISTS meta (
    key             VARCHAR PRIMARY KEY,
    value           VARCHAR
);

CREATE SEQUENCE IF NOT EXISTS directive_seq START 1;

CREATE TABLE IF NOT EXISTS directives (
    seq             BIGINT PRIMARY KEY DEFAULT nextval('directive_seq'),
    directive       VARCHAR NOT NULL
);

CREATE TABLE IF NOT EXISTS autoincrements (
    base            VARCHAR PRIMARY KEY,
    n               BIGINT
);

CREATE TABLE IF NOT EXISTS duplicates (
    original_id     VARCHAR NOT NULL,
    new_id          VARCHAR PRIMARY KEY
);
"""


# Post-load index DDL. Built only after bulk insert + closure materialization,
# per the Phase 2 plan. Note: B-tree on `(seqid, start, end)` is the universal
# fallback; the R-tree (if the spatial extension loads) is added separately in
# `ingest.py::_build_rtree`.
POST_LOAD_INDEXES = """
CREATE INDEX IF NOT EXISTS features_type     ON features(featuretype);
CREATE INDEX IF NOT EXISTS features_seqstart ON features(seqid, start, "end");
CREATE INDEX IF NOT EXISTS attributes_kv     ON attributes(key, value);
CREATE INDEX IF NOT EXISTS attributes_fid    ON attributes(feature_id);
CREATE INDEX IF NOT EXISTS edges_parent      ON edges(parent);
CREATE INDEX IF NOT EXISTS edges_child       ON edges(child);
CREATE INDEX IF NOT EXISTS closure_ancestor  ON closure(ancestor, depth);
CREATE INDEX IF NOT EXISTS closure_descend   ON closure(descendant, depth);
"""
# Phase 19: dropped redundant `features_seqid` — every (seqid)
# predicate is satisfied by the leading prefix of `features_seqstart`.


# ---------------------------------------------------------------------------
# Set-based normalization SQL (Phase 2 §3).
# ---------------------------------------------------------------------------

# 1. Edges from GFF3 `Parent=` attributes.
EDGES_FROM_PARENT = """
INSERT INTO edges (parent, child)
SELECT a.value AS parent, a.feature_id AS child
FROM attributes a
WHERE a.key = 'Parent';
"""

# 2. Edges from GTF gene_id / transcript_id (after gene/transcript rows have
#    been synthesized). The transcript->child edge is for any feature with a
#    transcript_id attribute. The gene->transcript edge is for any feature with
#    featuretype='transcript'.
EDGES_FROM_GTF = """
INSERT INTO edges (parent, child)
SELECT a.value AS parent, a.feature_id AS child
FROM attributes a
JOIN features f ON f.id = a.feature_id
WHERE a.key = 'transcript_id'
  AND f.featuretype <> 'transcript'
  AND a.value <> f.id        -- guard: don't create self-loops
  AND EXISTS (SELECT 1 FROM features p WHERE p.id = a.value)
UNION ALL
SELECT a.value AS parent, a.feature_id AS child
FROM attributes a
JOIN features f ON f.id = a.feature_id
WHERE a.key = 'gene_id'
  AND f.featuretype = 'transcript'
  AND a.value <> f.id
  AND EXISTS (SELECT 1 FROM features p WHERE p.id = a.value);
"""


# 3. GTF transcript synthesis. One scan with GROUP BY transcript_id over the
#    subfeature rows (typically 'exon'). Replaces 250k+ MIN/MAX queries from
#    the legacy gffutils with a single vectorized aggregation.
GTF_SYNTHESIZE_TRANSCRIPTS = """
INSERT INTO features (id, seqid, source, featuretype, start, "end",
                      score, strand, frame,
                      attributes_blob, extra_blob, file_order, is_synthetic)
SELECT
    a.value                    AS id,
    ANY_VALUE(f.seqid)         AS seqid,
    'gffbase_derived'        AS source,
    'transcript'               AS featuretype,
    MIN(f.start)               AS start,
    MAX(f."end")               AS "end",
    '.'                        AS score,
    ANY_VALUE(f.strand)        AS strand,
    '.'                        AS frame,
    NULL                       AS attributes_blob,
    NULL                       AS extra_blob,
    NULL                       AS file_order,
    TRUE                       AS is_synthetic
FROM features f
JOIN attributes a ON a.feature_id = f.id AND a.key = 'transcript_id'
WHERE f.featuretype = ?                          -- subfeature, e.g. 'exon'
  AND a.value NOT IN (SELECT id FROM features)
GROUP BY a.value;
"""

# 3b. Mirror the synthesized transcript_id back into the attributes table so
#     subsequent gene synthesis sees them.
GTF_SYNTHESIZE_TRANSCRIPT_ATTRS = """
INSERT INTO attributes (feature_id, key, value, idx)
SELECT f.id, 'transcript_id', f.id, 0
FROM features f
WHERE f.featuretype = 'transcript' AND f.is_synthetic = TRUE;
"""

# 3c. Carry gene_id forward onto each synthesized transcript by joining
#     attributes directly: each subfeature row carries (transcript_id, gene_id)
#     pairs in the attributes table; we group by transcript and pick the most
#     common gene_id (handles inconsistent GTFs gracefully). Crucially this
#     does NOT depend on the edges table, so we can defer edge population to
#     a single pass after all synthesis completes — eliminating duplicate rows.
GTF_PROPAGATE_GENE_ID = """
INSERT INTO attributes (feature_id, key, value, idx)
WITH pairs AS (
    SELECT
        tid.value AS transcript_id,
        gid.value AS gene_id,
        COUNT(*) AS cnt
    FROM attributes tid
    JOIN attributes gid ON gid.feature_id = tid.feature_id
    WHERE tid.key = 'transcript_id'
      AND gid.key = 'gene_id'
    GROUP BY tid.value, gid.value
),
ranked AS (
    SELECT
        transcript_id, gene_id, cnt,
        ROW_NUMBER() OVER (PARTITION BY transcript_id ORDER BY cnt DESC, gene_id) AS rn
    FROM pairs
)
SELECT t.id, 'gene_id', r.gene_id, 0
FROM features t
JOIN ranked r ON r.transcript_id = t.id AND r.rn = 1
WHERE t.featuretype = 'transcript' AND t.is_synthetic = TRUE;
"""

# 4. GTF gene synthesis. Same shape as transcript synthesis, but groups over
#    all rows whose featuretype IN ('transcript', subfeature) carrying a
#    gene_id attribute.
GTF_SYNTHESIZE_GENES = """
INSERT INTO features (id, seqid, source, featuretype, start, "end",
                      score, strand, frame,
                      attributes_blob, extra_blob, file_order, is_synthetic)
SELECT
    a.value                    AS id,
    ANY_VALUE(f.seqid)         AS seqid,
    'gffbase_derived'        AS source,
    'gene'                     AS featuretype,
    MIN(f.start)               AS start,
    MAX(f."end")               AS "end",
    '.'                        AS score,
    ANY_VALUE(f.strand)        AS strand,
    '.'                        AS frame,
    NULL                       AS attributes_blob,
    NULL                       AS extra_blob,
    NULL                       AS file_order,
    TRUE                       AS is_synthetic
FROM features f
JOIN attributes a ON a.feature_id = f.id AND a.key = 'gene_id'
WHERE f.featuretype IN ('transcript', ?)        -- transcript + subfeature
  AND a.value NOT IN (SELECT id FROM features)
GROUP BY a.value;
"""


# ---------------------------------------------------------------------------
# SQLite-compat views (Phase 2 §6). Make `FeatureDB.execute()` accept queries
# written against the legacy gffutils schema names. The `attributes` column
# is the raw col-9 bytes (NOT JSON) — documented break.
# ---------------------------------------------------------------------------
COMPAT_VIEWS_SQL = """
CREATE OR REPLACE VIEW features_compat AS
    SELECT id, seqid, source, featuretype, start, "end",
           score, strand, frame,
           CAST(attributes_blob AS VARCHAR) AS attributes,
           CAST(extra_blob      AS VARCHAR) AS extra,
           0 AS bin
    FROM features;

CREATE OR REPLACE VIEW relations_compat AS
    SELECT ancestor AS parent, descendant AS child, depth AS level
    FROM closure;
"""


# ---------------------------------------------------------------------------
# Recursive CTE closure (Phase 2 §3.1) — replaces the N+1 grandchild loop.
# Single SQL statement; vectorized executor.
# ---------------------------------------------------------------------------
CLOSURE_RECURSIVE_CTE = """
INSERT INTO closure (ancestor, descendant, depth)
WITH RECURSIVE walk(ancestor, descendant, depth) AS (
    SELECT parent, child, 1 AS depth FROM edges
    UNION ALL
    SELECT w.ancestor, e.child, w.depth + 1
    FROM walk w
    JOIN edges e ON e.parent = w.descendant
    WHERE w.depth < ?
)
SELECT ancestor, descendant, depth FROM walk;
"""
