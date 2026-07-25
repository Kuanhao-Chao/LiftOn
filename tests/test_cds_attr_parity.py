"""A rebuilt CDS must describe itself exactly like the reference row it came from.

The chaining merge (``Lifton_TRANS.update_cds_list``) and the ORF-rescue boundary
patch (``__update_cds_boundary``) both REBUILD a transcript's CDS list, and both
used to reset the emitted attributes to ``{Parent}``. The output was internally
inconsistent: measured on the full drosophila lift, 7,127 of 7,142 pure-Liftoff
CDS rows carried ``Dbxref``/``product``/``protein_id``/``gene``/``locus_tag``/…
while only 21 of 151,277 merged rows did — and every merged mRNA row kept those
same attributes, so the loss was confined to the CDS level.

No existing test caught it because the merge tests assert identity SCORES and CDS
IDs, never the descriptive attributes. These assert them directly, on both rebuild
paths, plus the two invariants that make the carry safe:
``ID`` stays owned by the write-funnel allocator, and sibling CDS never share one
attribute dict (``create_lifton_entries`` hands them all the same object).
"""
import gffutils
import pytest

from lifton.lifton_class import Lifton_TRANS, Lifton_EXON, Lifton_CDS


REF_CDS_ATTRS = {
    "ID": ["cds-NP_001.1"],
    "Parent": ["rna-NM_001.1"],
    "Dbxref": ["GeneID:1234", "GenBank:NP_001.1"],
    "Name": ["NP_001.1"],
    "gbkey": ["CDS"],
    "gene": ["actin"],
    "locus_tag": ["Dmel_CG1234"],
    "product": ["actin isoform A"],
    "protein_id": ["NP_001.1"],
}

DESCRIPTIVE_KEYS = {
    "Dbxref", "Name", "gbkey", "gene", "locus_tag", "product", "protein_id",
}


def _ent(ftype, start, end, fid, parent, strand="+", attributes=None):
    attrs = dict(attributes) if attributes else {"ID": [fid], "Parent": [parent]}
    f = gffutils.Feature(seqid="chr1", source="LiftOn", featuretype=ftype,
                         start=start, end=end, strand=strand, frame=".",
                         attributes=attrs)
    f.id = fid
    return f


def _trans(exon_spans, strand="+"):
    """A transcript with exons but no reference CDS attached yet."""
    trans = Lifton_TRANS.__new__(Lifton_TRANS)
    trans.entry = _ent("mRNA", exon_spans[0][0], exon_spans[-1][1],
                       "rna-NM_001.1", "gene-actin", strand)
    trans.exons = [
        Lifton_EXON(_ent("exon", s, e, f"exon-NM_001.1-{i}", "rna-NM_001.1",
                         strand))
        for i, (s, e) in enumerate(exon_spans, 1)
    ]
    return trans


def _ref_cds(start, end):
    return _ent("CDS", start, end, "cds-NP_001.1", "rna-NM_001.1",
                attributes=dict(REF_CDS_ATTRS))


def _emitted(trans):
    return [ex.cds.entry.attributes for ex in trans.exons if ex.cds is not None]


# ── the template ────────────────────────────────────────────────────────────

def test_add_cds_captures_the_transcript_attribute_template():
    trans = _trans([(100, 200)])
    trans.add_cds(_ref_cds(120, 180))

    assert set(trans._cds_attr_template) == DESCRIPTIVE_KEYS
    assert trans._cds_attr_template["product"] == ["actin isoform A"]
    # Structural keys are never part of the template.
    for key in ("ID", "Parent", "extra_copy_number"):
        assert key not in trans._cds_attr_template


def test_template_ignores_an_attributeless_reference_cds():
    trans = _trans([(100, 200)])
    trans.add_cds(_ent("CDS", 120, 180, "cds-x", "rna-NM_001.1"))
    assert trans._cds_attr_template is None


def test_template_is_a_copy_not_a_live_view_of_the_reference_row():
    trans = _trans([(100, 200)])
    ref = _ref_cds(120, 180)
    trans.add_cds(ref)
    ref.attributes["product"] = ["MUTATED"]
    assert trans._cds_attr_template["product"] == ["actin isoform A"]


# ── the chaining merge path ─────────────────────────────────────────────────

@pytest.mark.parametrize("optimize", [False, True])
def test_merged_cds_keeps_the_reference_descriptive_attributes(optimize):
    trans = _trans([(100, 200), (300, 400)])
    trans.add_cds(_ref_cds(120, 180))

    # The merge rebuilds the CDS list from scratch, as chaining does.
    trans.update_cds_list([Lifton_CDS(_ref_cds(120, 190))], optimize=optimize)

    emitted = _emitted(trans)
    assert emitted, "the rebuilt transcript emitted no CDS"
    for attrs in emitted:
        assert DESCRIPTIVE_KEYS <= set(attrs)
        assert attrs["product"] == ["actin isoform A"]
        assert attrs["protein_id"] == ["NP_001.1"]
        assert attrs["Dbxref"] == ["GeneID:1234", "GenBank:NP_001.1"]


def test_merged_cds_reparents_and_leaves_id_to_the_allocator():
    trans = _trans([(100, 200), (300, 400)])
    trans.add_cds(_ref_cds(120, 180))
    trans.update_cds_list([Lifton_CDS(_ref_cds(120, 190))])

    for attrs in _emitted(trans):
        assert attrs["Parent"] == ["rna-NM_001.1"]
        # `_normalize_cds_ids` owns the ID at write time; a stale reference ID
        # carried here would collide across segments and copies.
        assert "ID" not in attrs


def test_sibling_cds_do_not_share_one_attribute_dict():
    """`create_lifton_entries` gives every chained CDS the SAME attribute dict.

    If that object survived into the emitted rows, `_normalize_cds_ids`'
    per-segment ``attributes["ID"] = [...]`` would overwrite every sibling.
    """
    trans = _trans([(100, 200), (300, 400), (500, 600)])
    trans.add_cds(_ref_cds(120, 180))

    shared = dict(REF_CDS_ATTRS)
    cds_list = []
    for start, end in ((120, 200), (300, 400), (500, 560)):
        entry = _ent("CDS", start, end, "cds-NP_001.1", "rna-NM_001.1",
                     attributes=shared)
        entry.attributes = shared          # the aliasing the merge really does
        cds_list.append(Lifton_CDS(entry))

    trans.update_cds_list(cds_list)

    emitted = _emitted(trans)
    assert len(emitted) >= 2
    assert len({id(attrs) for attrs in emitted}) == len(emitted)

    # Mutating one row must not disturb the others.
    emitted[0]["ID"] = ["cds-first"]
    assert all("ID" not in attrs for attrs in emitted[1:])


def test_miniprot_sourced_cds_falls_back_to_its_own_attributes():
    """No reference CDS -> no template, so the source row's own set is used."""
    trans = _trans([(100, 200)])
    entry = _ent("CDS", 120, 180, "mp-cds", "rna-NM_001.1",
                 attributes={"ID": ["mp-cds"], "Parent": ["x"],
                             "Target": ["NP_001.1 1 60"], "Rank": ["1"]})
    trans.update_cds_list([Lifton_CDS(entry)])

    attrs = _emitted(trans)[0]
    assert attrs["Target"] == ["NP_001.1 1 60"]
    assert attrs["Parent"] == ["rna-NM_001.1"]
    assert "ID" not in attrs


# ── the ORF-rescue path ─────────────────────────────────────────────────────

def test_novel_cds_uses_the_template_not_the_exon_attributes():
    """An exon row carries ``gbkey=mRNA`` and no ``protein_id``.

    Copying it onto a CDS would mislabel the feature, so the ORF rescue takes
    the transcript's CDS template instead.
    """
    trans = _trans([(100, 200)])
    trans.add_cds(_ref_cds(120, 180))
    exon = trans.exons[0]
    exon.entry.attributes["gbkey"] = ["mRNA"]
    exon.entry.attributes["transcript_id"] = ["NM_001.1"]

    exon.add_novel_lifton_cds(exon.entry, 130, 170,
                              attr_template=trans._cds_attr_template)

    attrs = exon.cds.entry.attributes
    assert attrs["gbkey"] == ["CDS"]
    assert attrs["protein_id"] == ["NP_001.1"]
    assert "transcript_id" not in attrs
    assert attrs["Parent"] == ["rna-NM_001.1"]
    # A novel ORF boundary is not the reference CDS: the allocator synthesizes.
    assert exon.cds._source_id is None
    assert "ID" not in attrs


def test_novel_cds_without_a_template_stays_parent_only():
    trans = _trans([(100, 200)])
    exon = trans.exons[0]
    exon.add_novel_lifton_cds(exon.entry, 130, 170, attr_template=None)
    assert exon.cds.entry.attributes == {"Parent": ["rna-NM_001.1"]}


# ── the escape hatch ────────────────────────────────────────────────────────

def test_env_gate_restores_the_parent_only_shape(monkeypatch):
    monkeypatch.setenv("LIFTON_NO_CDS_ATTR_CARRY", "1")

    trans = _trans([(100, 200), (300, 400)])
    trans.add_cds(_ref_cds(120, 180))
    trans.update_cds_list([Lifton_CDS(_ref_cds(120, 190))])

    for attrs in _emitted(trans):
        assert attrs == {"Parent": ["rna-NM_001.1"]}


def test_env_gate_also_covers_the_orf_rescue_path(monkeypatch):
    monkeypatch.setenv("LIFTON_NO_CDS_ATTR_CARRY", "1")

    trans = _trans([(100, 200)])
    trans.add_cds(_ref_cds(120, 180))
    exon = trans.exons[0]
    exon.add_novel_lifton_cds(exon.entry, 130, 170,
                              attr_template=trans._cds_attr_template)

    assert exon.cds.entry.attributes == {"Parent": ["rna-NM_001.1"]}
