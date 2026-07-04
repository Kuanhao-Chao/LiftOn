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
