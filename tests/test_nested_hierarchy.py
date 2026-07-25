"""Tier-2 audit fix — a 3-level reference hierarchy must not silently lose transcripts.

`process_liftoff` assumed a 2-level hierarchy. When a gene's level-1 child has no direct
exons (RefSeq `gene -> primary_transcript -> miRNA -> exon`), that child becomes a
`LiftOn_FEATURE` and the recursion passes it in the `lifton_gene` slot. Two things then
went wrong:

1. The reference ids were resolved from the INTERMEDIATE feature's id rather than a gene
   id, so the lookup missed, `ref_db[None]` raised, and the transcript plus all of its
   exons and CDS were dropped -- with no log line and no failure record.
2. If the id did happen to resolve, `lifton_add_trans_exon_cds` called
   `add_transcript` on the `LiftOn_FEATURE`, which had no such method -> AttributeError
   -> `process_locus` dropped the whole locus.

There was no 3-level fixture anywhere in the suite.
"""
import gffutils
import pytest

from lifton import lifton_class


def _ent(ftype, start, end, fid, parent=None):
    attrs = {"ID": [fid]}
    if parent:
        attrs["Parent"] = [parent]
    f = gffutils.Feature(seqid="chr1", source="ref", featuretype=ftype, start=start,
                         end=end, strand="+", frame=".", attributes=attrs)
    f.id = fid
    return f


def _feature():
    """An intermediate LiftOn_FEATURE, as built by Lifton_GENE.add_feature."""
    return lifton_class.LiftOn_FEATURE(
        "gene1", _ent("primary_transcript", 100, 400, "pt1", "gene1"), 0)


def test_feature_exposes_the_transcript_bearing_api():
    """The delegates `process_liftoff` needs at depth >= 2."""
    for name in ("add_transcript", "add_exon", "add_cds", "orf_search_protein",
                 "add_lifton_gene_status_attrs", "add_lifton_trans_status_attrs"):
        assert hasattr(lifton_class.LiftOn_FEATURE, name), name


def test_nested_transcript_with_exons_and_cds_is_retained():
    feature = _feature()

    # Lifton_TRANS keys itself by the REFERENCE transcript id, exactly as
    # Lifton_GENE.add_transcript does; the caller uses the returned object's id.
    trans = feature.add_transcript(
        "rna-ref", _ent("miRNA", 100, 400, "mir1", "pt1"), {"Parent": ["pt1"]})
    assert trans is not None
    tid = trans.entry.id
    feature.add_exon(tid, _ent("exon", 100, 199, "e1", tid))
    feature.add_exon(tid, _ent("exon", 300, 400, "e2", tid))
    feature.add_cds(tid, _ent("CDS", 100, 199, "c1", tid))

    # The transcript is reachable from the intermediate feature...
    assert tid in feature.features
    child = feature.features[tid]
    assert isinstance(child, lifton_class.Lifton_TRANS)
    # ...and carries its exons and the CDS.
    assert [(e.entry.start, e.entry.end) for e in child.exons] == [(100, 199), (300, 400)]
    assert child.exons[0].cds is not None


def test_nested_transcript_is_written_under_its_parent(tmp_path):
    """The nested transcript must reach the output, not be dropped."""
    import io

    feature = _feature()
    trans = feature.add_transcript(
        "rna-ref", _ent("miRNA", 100, 400, "mir1", "pt1"), {"Parent": ["pt1"]})
    feature.add_exon(trans.entry.id, _ent("exon", 100, 199, "e1", trans.entry.id))

    buf = io.StringIO()
    feature.write_entry(buf)
    text = buf.getvalue()

    types = [ln.split("\t")[2] for ln in text.splitlines() if ln.strip()]
    assert "primary_transcript" in types
    assert "miRNA" in types, f"nested transcript was dropped: {types}"
    assert "exon" in types


def test_status_delegates_are_safe_for_plain_feature_children():
    """A LiftOn_FEATURE child (not a transcript) must not break the status calls."""
    feature = _feature()
    feature.add_feature(_ent("misc_RNA", 120, 380, "sub1", "pt1"))

    feature.add_lifton_gene_status_attrs("Liftoff")
    feature.add_lifton_trans_status_attrs("sub1", lifton_class.Lifton_Status())
    assert feature.entry.attributes["source"] == ["Liftoff"]

    # An unknown id is a no-op rather than a KeyError.
    feature.add_lifton_trans_status_attrs("nope", lifton_class.Lifton_Status())
    assert feature.orf_search_protein("nope", "r", None, {}, {}, None) == (None, False)


def test_update_boundaries_extends_to_cover_children():
    feature = _feature()
    feature.add_transcript(
        "rna-ref", _ent("miRNA", 50, 500, "mir1", "pt1"), {"Parent": ["pt1"]})
    feature.update_boundaries()
    assert (feature.entry.start, feature.entry.end) == (50, 500)
