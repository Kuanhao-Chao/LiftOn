"""Every segment of one chained CDS must agree about its own identity.

`create_lifton_entries` handed every chained CDS of a transcript the SAME
attribute dict, and `Lifton_CDS.__init__` reads `extra_copy_number` and then
POPS it. On a copy gene the first segment therefore saw the key and had its
copy suffix stripped from `_source_id_base`, while every sibling -- reading the
now-emptied dict -- kept the suffixed form:

    ['cds-NP_9.1', 'cds-NP_9.1_1', 'cds-NP_9.1_1']

`_normalize_cds_ids` then sees two distinct "pre-existing" IDs where the
reference has one discontinuous CDS, takes its multi-ID branch, and emits a
transcript whose CDS segments carry inconsistent identifiers.

Found while building the CDS-attribute carry and deliberately left alone there
so that commit stayed about attributes.
"""
import gffutils
import pytest

from lifton import protein_maximization as pm
from lifton.lifton_class import Lifton_TRANS, Lifton_EXON, Lifton_CDS


def _cds(start, end, *, copy_number="1", cds_id="cds-NP_9.1_1"):
    attrs = {"ID": [cds_id], "Parent": ["rna-NM_9_1"], "product": ["thing"],
             "protein_id": ["NP_9.1"]}
    if copy_number is not None:
        attrs["extra_copy_number"] = [copy_number]
    feature = gffutils.Feature(
        seqid="chr1", source="Liftoff", featuretype="CDS", start=start,
        end=end, strand="+", frame="0", attributes=attrs)
    feature.id = cds_id
    return feature


class _Aln:
    def __init__(self, children, strand="+"):
        self.cds_children = children
        self.db_entry = type("Entry", (), {"strand": strand})()


def _chain(children, strand="+"):
    aln = _Aln(children, strand)
    return pm.create_lifton_entries(
        len(children), 0, aln, len(children), 0, aln, None, False)


def test_every_segment_of_a_copy_gene_shares_one_source_id_base():
    chained = _chain([_cds(100, 199), _cds(300, 399), _cds(500, 599)])
    bases = [cds._source_id_base for cds in chained]
    assert bases == ["cds-NP_9.1"] * 3, bases


def test_every_segment_reports_the_copy_number():
    chained = _chain([_cds(100, 199), _cds(300, 399)])
    assert [c._source_copy_number for c in chained] == ["1", "1"]


def test_segments_do_not_share_one_attribute_dict():
    chained = _chain([_cds(100, 199), _cds(300, 399), _cds(500, 599)])
    assert len({id(c.entry.attributes) for c in chained}) == 3


@pytest.mark.parametrize("strand", ["+", "-"])
def test_non_copy_genes_are_unaffected(strand):
    children = [_cds(100, 199, copy_number=None, cds_id="cds-NP_9.1"),
                _cds(300, 399, copy_number=None, cds_id="cds-NP_9.1")]
    chained = _chain(children, strand)
    assert [c._source_id_base for c in chained] == ["cds-NP_9.1"] * 2
    assert all(c._source_copy_number is None for c in chained)


def test_extra_copy_number_never_reaches_the_output():
    for cds in _chain([_cds(100, 199), _cds(300, 399)]):
        assert "extra_copy_number" not in cds.entry.attributes


def test_the_shared_base_collapses_to_one_emitted_id(make_gffutils_feature):
    """End-to-end: the write funnel must emit ONE discontinuous-CDS ID.

    With the segments disagreeing, `_normalize_cds_ids` took its multi-ID
    branch and produced inconsistent identifiers for one CDS.
    """
    trans = Lifton_TRANS.__new__(Lifton_TRANS)
    trans.entry = make_gffutils_feature(
        featuretype="mRNA", start=100, end=599, feature_id="rna-NM_9_1",
        attributes={"ID": ["rna-NM_9_1"], "Parent": ["gene-A_1"]})
    trans.entry.id = "rna-NM_9_1"
    trans.ref_tran_id = "rna-NM_9"
    trans.exons = []
    for index, (start, end) in enumerate(((100, 199), (300, 399), (500, 599)), 1):
        trans.exons.append(Lifton_EXON(make_gffutils_feature(
            featuretype="exon", start=start, end=end,
            feature_id=f"exon-NM_9_1-{index}",
            attributes={"ID": [f"exon-NM_9_1-{index}"],
                        "Parent": ["rna-NM_9_1"]})))

    for exon, cds in zip(trans.exons, _chain(
            [_cds(100, 199), _cds(300, 399), _cds(500, 599)])):
        exon.cds = cds
        cds.entry.attributes["Parent"] = ["rna-NM_9_1"]

    trans.normalize_containment()

    emitted = [e.cds.entry.attributes["ID"] for e in trans.exons]
    assert len({tuple(i) for i in emitted}) == 1, emitted
    # The copy suffix belongs on the emitted ID, applied once from the
    # transcript -- not baked into some segments and not others.
    assert emitted[0] == ["cds-NP_9.1_1"], emitted[0]


def test_lifton_cds_still_derives_provenance_when_built_alone():
    """The stripping itself is correct; only its per-segment consistency broke."""
    cds = Lifton_CDS(_cds(100, 199))
    assert cds._source_id == "cds-NP_9.1_1"
    assert cds._source_id_base == "cds-NP_9.1"
    assert cds._source_copy_number == "1"
