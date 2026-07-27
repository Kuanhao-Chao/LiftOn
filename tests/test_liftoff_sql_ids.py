"""Regression test for GitHub issue #35 — vendored Liftoff double-quoted SQL IDs.

``extract_features.{add_parents,add_children,add_intermediates}`` built their
``id IN (...)`` clauses as ``', '.join('"{0}"'.format(w) for w in ids)``. A
double-quoted ``"x"`` is an SQL *identifier* (a column name), not a string
literal: modern SQLite raises ``no such column: x`` (older SQLite's
double-quoted-string-literal misfeature silently accepted it — hence the
reporter's "works on one machine, not another"). The fix binds the IDs as query
parameters (``?``), which is correct on every SQLite build and also handles IDs
containing SQL-significant characters.

The second test is a *machine-independent* guard: a gene ID containing a double
quote (``gene"Q``) makes the OLD ``'"{0}"'.format(w)`` interpolation malformed
SQL (``near "Q": syntax error``) on ANY SQLite version, so a revert to the old
form fails here regardless of the local SQLite's DQS setting.
"""
import sqlite3

import gffutils

from lifton.liftoff import extract_features


def _db(tmp_path, gff_text):
    p = tmp_path / "ref.gff3"
    p.write_text(gff_text)
    return gffutils.create_db(
        str(p), str(tmp_path / "ref.db"),
        merge_strategy="create_unique", force=True, keep_order=True)


_GFF = """##gff-version 3
chr1\tsrc\tgene\t101\t399\t.\t+\t.\tID={gid}
chr1\tsrc\tmRNA\t101\t399\t.\t+\t.\tID=tx1;Parent={gid}
chr1\tsrc\texon\t101\t199\t.\t+\t.\tID=ex1;Parent=tx1
chr1\tsrc\texon\t301\t399\t.\t+\t.\tID=ex2;Parent=tx1
chr1\tsrc\tCDS\t101\t199\t.\t+\t0\tID=cds1;Parent=tx1
"""


def test_seperate_parents_and_children_plain_id(tmp_path):
    """The three parameterized IN-clause queries extract the full hierarchy:
    the top gene as a parent, the mRNA as an intermediate, exon+CDS as leaves."""
    db = _db(tmp_path, _GFF.format(gid="gene1"))
    fh, order = extract_features.seperate_parents_and_children(db, {"gene"})
    assert "gene1" in fh.parents
    assert "tx1" in fh.intermediates
    leaf_types = {c.featuretype for c in fh.children["gene1"]}
    assert {"exon", "CDS"} <= leaf_types
    assert [p[0] for p in order] == ["gene1"]


def test_id_with_double_quote_binds_as_parameter(tmp_path):
    """A gene ID containing `"` would make the OLD `'"{0}"'.format(w)` clause
    malformed SQL on EVERY sqlite build; the parameterized fix binds it cleanly."""
    db = _db(tmp_path, _GFF.format(gid='gene"Q'))
    fh, order = extract_features.seperate_parents_and_children(db, {"gene"})
    assert 'gene"Q' in fh.parents
    assert any(c.featuretype == "CDS" for c in fh.children['gene"Q'])


class _VariableLimitedCursor:
    def __init__(self, cursor, connection):
        self._cursor = cursor
        self.connection = connection

    def execute(self, query, parameters=()):
        parameter_count = len(parameters)
        self.connection.parameter_counts.append(parameter_count)
        if parameter_count > self.connection.variable_limit:
            raise sqlite3.OperationalError("too many SQL variables")
        return self._cursor.execute(query, parameters)


class _VariableLimitedConnection:
    def __init__(self, connection, variable_limit):
        self._connection = connection
        self.variable_limit = variable_limit
        self.parameter_counts = []

    def cursor(self):
        return _VariableLimitedCursor(self._connection.cursor(), self)

    def getlimit(self, category):
        return self.variable_limit


def test_large_leaf_set_is_batched_below_sqlite_variable_limit(tmp_path):
    """A leaf set above the bind limit is streamed without changing order."""
    exon_ids = [f"ex{24 - index:02d}" for index in range(25)]
    exons = [
        (
            f"chr1\tsrc\texon\t{position}\t{position}\t.\t+\t.\t"
            f"ID={exon_id};Parent=tx1"
        )
        for position, exon_id in zip(range(101, 151, 2), exon_ids)
    ]
    gff = "\n".join([
        "##gff-version 3",
        "chr1\tsrc\tgene\t101\t151\t.\t+\t.\tID=gene1",
        "chr1\tsrc\tmRNA\t101\t151\t.\t+\t.\tID=tx1;Parent=gene1",
        *exons,
        "",
    ])
    db = _db(tmp_path, gff)
    expected_hierarchy, expected_order = (
        extract_features.seperate_parents_and_children(db, {"gene"}))
    limited_connection = _VariableLimitedConnection(db.conn, variable_limit=8)
    db.conn = limited_connection

    hierarchy, order = extract_features.seperate_parents_and_children(
        db, {"gene"})

    assert [parent[0] for parent in order] == [
        parent[0] for parent in expected_order
    ]
    assert list(hierarchy.parents) == list(expected_hierarchy.parents)
    assert list(hierarchy.intermediates) == list(expected_hierarchy.intermediates)
    expected_children = [(exon_id, "exon") for exon_id in exon_ids]
    assert [
        (child.id, child.featuretype)
        for child in expected_hierarchy.children["gene1"]
    ] == expected_children
    assert [
        (child.id, child.featuretype)
        for child in hierarchy.children["gene1"]
    ] == expected_children
    assert max(limited_connection.parameter_counts) <= 8
    assert 0 in limited_connection.parameter_counts


def test_batched_hierarchy_preserves_database_order(tmp_path):
    """Shuffled IDs must retain parent, intermediate, and child row order."""
    gene_ids = [
        "z09", "a00", "m05", "b01", "y08", "c02",
        "x07", "d03", "w06", "e04", "v11", "f10",
    ]
    lines = ["##gff-version 3"]
    for index, gene_id in enumerate(gene_ids):
        start = 101 + index * 10
        transcript_id = f"tx-{gene_id}"
        lines.extend([
            f"chr1\tsrc\tgene\t{start}\t{start + 8}\t.\t+\t.\tID={gene_id}",
            (
                f"chr1\tsrc\tmRNA\t{start}\t{start + 8}\t.\t+\t.\t"
                f"ID={transcript_id};Parent={gene_id}"
            ),
            (
                f"chr1\tsrc\texon\t{start}\t{start + 8}\t.\t+\t.\t"
                f"ID=ex-{gene_id};Parent={transcript_id}"
            ),
        ])
    db = _db(tmp_path, "\n".join([*lines, ""]))
    expected_hierarchy, expected_order = (
        extract_features.seperate_parents_and_children(db, {"gene"}))
    limited_connection = _VariableLimitedConnection(db.conn, variable_limit=3)
    db.conn = limited_connection

    hierarchy, order = extract_features.seperate_parents_and_children(
        db, {"gene"})

    transcript_ids = [f"tx-{gene_id}" for gene_id in gene_ids]
    assert list(expected_hierarchy.parents) == gene_ids
    assert list(hierarchy.parents) == gene_ids
    assert list(expected_hierarchy.intermediates) == transcript_ids
    assert list(hierarchy.intermediates) == transcript_ids
    assert [parent[0] for parent in expected_order] == gene_ids
    assert [parent[0] for parent in order] == gene_ids
    for gene_id in gene_ids:
        assert [
            (child.id, child.featuretype)
            for child in hierarchy.children[gene_id]
        ] == [
            (child.id, child.featuretype)
            for child in expected_hierarchy.children[gene_id]
        ]
    assert max(limited_connection.parameter_counts) <= 3
    assert limited_connection.parameter_counts.count(3) == 8
