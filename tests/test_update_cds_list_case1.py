"""Tier-2 audit fix — `update_cds_list` Case 1 must not duplicate the CDS-bearing exon.

Case 1 (a single CDS) walks the exon list and, once it has emitted the merged
CDS-bearing exon, must emit every FURTHER downstream exon as a plain 3'-UTR exon. That
guard was gated on `optimize`, so with ``optimize=False`` -- the ``--legacy-merge`` path
and the default of ``Lifton_GENE.update_cds_list`` -- every downstream exon re-entered
the else-branch and re-emitted ANOTHER copy of the CDS-bearing exon: duplicate exon IDs,
a doubled CDS, and a wrong protein on any multi-exon transcript whose CDS covers only
the leading exon(s). The bug scaled with exon count.

No existing test caught it because the merge tests compare identity SCORES, never the
emitted exon list — so these assert the structure directly.
"""
import pytest
import gffutils

from lifton.lifton_class import Lifton_TRANS, Lifton_EXON, Lifton_CDS


def _ent(ftype, start, end, fid, parent, strand="+"):
    f = gffutils.Feature(seqid="chr1", source="LiftOn", featuretype=ftype,
                         start=start, end=end, strand=strand, frame=".",
                         attributes={"ID": [fid], "Parent": [parent]})
    f.id = fid
    return f


def _trans(exon_spans, strand="+"):
    trans = Lifton_TRANS.__new__(Lifton_TRANS)
    trans.entry = _ent("mRNA", exon_spans[0][0], exon_spans[-1][1],
                       "rna-T1", "gene1", strand)
    trans.exons = [
        Lifton_EXON(_ent("exon", s, e, f"exon-T1-{i}", "rna-T1", strand))
        for i, (s, e) in enumerate(exon_spans, 1)
    ]
    return trans


def _structure(trans):
    return [((ex.entry.start, ex.entry.end),
             None if ex.cds is None else (ex.cds.entry.start, ex.cds.entry.end))
            for ex in trans.exons]


# CDS lies wholly inside exon 1; exons 2 and 3 are pure 3' UTR.
EXONS = [(100, 200), (300, 400), (500, 600)]
CDS_SPAN = (120, 180)


@pytest.mark.parametrize("optimize", [False, True])
def test_single_cds_emits_each_downstream_exon_once(optimize):
    trans = _trans(EXONS)
    cds = Lifton_CDS(_ent("CDS", *CDS_SPAN, "cds-T1", "rna-T1"))

    trans.update_cds_list([cds], optimize=optimize)

    structure = _structure(trans)
    # Exactly three exons: the CDS-bearing one plus the two untouched UTR exons.
    assert len(structure) == 3, structure
    coding = [s for s in structure if s[1] is not None]
    assert len(coding) == 1, f"CDS-bearing exon duplicated: {structure}"
    # The downstream UTR exons keep their own coordinates and carry no CDS.
    assert [s[0] for s in structure[1:]] == [(300, 400), (500, 600)]
    assert all(s[1] is None for s in structure[1:])


@pytest.mark.parametrize("optimize", [False, True])
def test_emitted_exon_ids_are_unique(optimize):
    """End-to-end ID contract, taken through the write funnel.

    Case 1 builds its merged exon as a deepcopy of the *downstream* exon, so the
    merged exon inherits that exon's ID -- a separate, pre-existing quirk that
    predates this fix and is masked because `normalize_containment` (the write
    funnel) renumbers exon IDs whenever they collide. Assert the guarantee that
    actually reaches the output rather than the intermediate state.
    """
    trans = _trans(EXONS)
    cds = Lifton_CDS(_ent("CDS", *CDS_SPAN, "cds-T1", "rna-T1"))

    trans.update_cds_list([cds], optimize=optimize)
    trans.normalize_containment()

    ids = [ex.entry.attributes["ID"][0] for ex in trans.exons]
    assert len(ids) == len(set(ids)), f"duplicate exon IDs: {ids}"


def test_optimize_false_matches_optimize_true():
    """The two paths must agree — the guard is correct for both."""
    a, b = _trans(EXONS), _trans(EXONS)
    a.update_cds_list([Lifton_CDS(_ent("CDS", *CDS_SPAN, "c", "rna-T1"))],
                      optimize=False)
    b.update_cds_list([Lifton_CDS(_ent("CDS", *CDS_SPAN, "c", "rna-T1"))],
                      optimize=True)
    assert _structure(a) == _structure(b)


def test_scales_with_exon_count():
    """The old bug re-emitted one duplicate per downstream exon, so a longer
    transcript degraded further; pin that it does not."""
    spans = [(100, 200)] + [(300 + i * 200, 400 + i * 200) for i in range(6)]
    trans = _trans(spans)
    trans.update_cds_list(
        [Lifton_CDS(_ent("CDS", *CDS_SPAN, "c", "rna-T1"))], optimize=False)
    structure = _structure(trans)
    assert len(structure) == len(spans), structure
    assert sum(1 for s in structure if s[1] is not None) == 1, structure
