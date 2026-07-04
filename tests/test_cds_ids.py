"""GitHub issue #32 / #8 — CDS features must carry an ``ID=``.

``Lifton_EXON.add_lifton_cds`` / ``add_novel_lifton_cds`` reset each CDS's
attributes to ``{Parent}`` only, so emitted CDS lines had no ``ID=``.
``Lifton_TRANS.normalize_containment`` (the write funnel, Iter-24) now assigns a
shared ``cds-<trans>`` to a transcript's CDS segments — the valid GFF3
discontinuous-CDS form (same ID + same Parent + same type), which
``gff3_validator.DISCONTINUOUS_FEATURE_TYPES = {"CDS"}`` already anticipates.
Gated by the same ``LIFTON_NO_CONTAINMENT_NORMALIZE`` escape hatch as the rest of
the write-funnel normalization, so it reproduces the pre-fix bytes when disabled.
"""
import gffutils

from lifton.lifton_class import Lifton_TRANS, Lifton_EXON, Lifton_CDS


def _ent(ftype, s, e, fid, parent):
    f = gffutils.Feature(seqid="chr1", source="LiftOn", featuretype=ftype, start=s,
                         end=e, strand="+", frame=".",
                         attributes={"ID": [fid], "Parent": [parent]})
    f.id = fid                       # the gffutils DB populates .id; set it here
    return f


def _trans(trans_id, cds_spans):
    trans = Lifton_TRANS.__new__(Lifton_TRANS)
    trans.entry = _ent("mRNA", 101, 399, trans_id, "gene9")
    trans.exons = []
    for i, (cs, ce) in enumerate(cds_spans, 1):
        ex = Lifton_EXON(_ent("exon", cs, ce, f"exon-A-{i}", trans_id))
        ex.cds = Lifton_CDS(_ent("CDS", cs, ce, "stale", trans_id))
        ex.cds.entry.attributes = {"Parent": [trans_id]}   # add_lifton_cds resets (no ID)
        trans.exons.append(ex)
    return trans


def test_cds_segments_share_one_cds_id():
    trans = _trans("rna-NM_9", [(101, 199), (301, 399)])
    trans.normalize_containment()
    ids = [e.cds.entry.attributes["ID"] for e in trans.exons]
    assert ids == [["cds-NM_9"], ["cds-NM_9"]]        # shared across segments, "rna-" stripped


def test_cds_id_without_rna_prefix():
    trans = _trans("mytx", [(101, 199)])
    trans.normalize_containment()
    assert trans.exons[0].cds.entry.attributes["ID"] == ["cds-mytx"]


def test_cds_id_reverts_under_escape_hatch(monkeypatch):
    monkeypatch.setenv("LIFTON_NO_CONTAINMENT_NORMALIZE", "1")
    trans = _trans("rna-NM_9", [(101, 199), (301, 399)])
    trans.normalize_containment()
    # normalize is a no-op -> CDS keep their ID-less {Parent} attributes.
    assert "ID" not in trans.exons[0].cds.entry.attributes
