"""Two exons of one transcript must not overlap.

Measured on four whole genomes -- CHM13 27 transcripts (0.019%), human ->
zebrafish 96 (0.144%), rice 30, bee 13 -- against reference annotations that
have **zero**. This is the open half of issue #26.

All 27 CHM13 cases carry `status=LiftOn_chaining_algorithm`. They are born in
`Lifton_TRANS.update_cds_list`: a rebuilt exon takes its end from a chained CDS
end that can lie inside the *next* Liftoff exon, and the "append any remaining
3' UTR exons" step then appends that exon verbatim without ever comparing it to
the exon just emitted.

`normalize_containment` does not merge them, and its collision renumbering
gives the pair distinct IDs -- which is why the mandatory `duplicate_id` gate
passes and the overlap ships.

The canonical case is POLR2A / rna-NM_000937.5, whose two reference CDS blocks
(756 bp + 401 bp) chain into one contiguous 1,157 bp block at 7417085-7418241
while Liftoff's second exon, 7417841-7418678, is appended untouched.
"""
import copy

import pytest

from lifton import lifton_class


def _exon(make_feature, start, end, cds=None):
    exon = lifton_class.Lifton_EXON(make_feature(
        featuretype="exon", start=start, end=end,
        attributes={"Parent": ["tx1"]}))
    if cds is not None:
        exon.cds = lifton_class.Lifton_CDS(make_feature(
            featuretype="CDS", start=cds[0], end=cds[1], frame="0"))
    return exon


@pytest.fixture
def make_feature(make_gffutils_feature):
    return make_gffutils_feature


def _spans(exons):
    return [(e.entry.start, e.entry.end) for e in exons]


class TestReconcileOverlappingExons:
    def test_non_overlapping_exons_are_returned_untouched(self, make_feature):
        """Byte-identity: the common case must not even be rebuilt."""
        exons = [_exon(make_feature, 101, 199), _exon(make_feature, 301, 399)]
        out = lifton_class.reconcile_overlapping_exons(exons)
        assert out is exons

    def test_the_polr2a_shape_becomes_one_exon(self, make_feature):
        """A coding exon whose CDS ran into the next exon, plus that exon."""
        exons = [
            _exon(make_feature, 7417085, 7418241, cds=(7417085, 7418241)),
            _exon(make_feature, 7417841, 7418678),
        ]
        out = lifton_class.reconcile_overlapping_exons(exons)
        assert _spans(out) == [(7417085, 7418678)], (
            "the overlapping pair must collapse to the contiguous exon they "
            "describe, keeping the 3' UTR tail"
        )
        assert (out[0].cds.entry.start, out[0].cds.entry.end) == (7417085, 7418241), (
            "the CDS must not move -- the protein cannot change"
        )

    def test_a_contained_exon_is_absorbed(self, make_feature):
        exons = [_exon(make_feature, 1353467, 1353692),
                 _exon(make_feature, 1353604, 1353678)]
        assert _spans(lifton_class.reconcile_overlapping_exons(exons)) == [
            (1353467, 1353692)]

    def test_two_overlapping_coding_exons_merge_their_cds(self, make_feature):
        """FAM118A / rna-XM_024452255.2: both sides carry a CDS, and the two
        CDS themselves overlap by 27 bp -- a protein that double-counts those
        bases. One exon, one CDS spanning the union."""
        exons = [
            _exon(make_feature, 45813076, 45813324, cds=(45813076, 45813324)),
            _exon(make_feature, 45813298, 45813405, cds=(45813298, 45813384)),
        ]
        out = lifton_class.reconcile_overlapping_exons(exons)
        assert _spans(out) == [(45813076, 45813405)]
        assert (out[0].cds.entry.start, out[0].cds.entry.end) == (45813076, 45813384)

    def test_a_disjoint_cds_pair_is_left_alone(self, make_feature):
        """The unrepresentable case, reported rather than silently mangled.

        An exon holds at most one CDS, so two coding blocks separated by a real
        gap cannot merge into one exon without inventing coding sequence. Leave
        the structure as it is and let the validator flag it.
        """
        exons = [
            _exon(make_feature, 100, 300, cds=(100, 150)),
            _exon(make_feature, 250, 400, cds=(300, 400)),
        ]
        out = lifton_class.reconcile_overlapping_exons(exons)
        assert len(out) == 2, "must not fabricate a merged coding block"

    def test_it_is_idempotent(self, make_feature):
        exons = [_exon(make_feature, 10, 100, cds=(10, 100)),
                 _exon(make_feature, 80, 200)]
        once = lifton_class.reconcile_overlapping_exons(exons)
        assert _spans(lifton_class.reconcile_overlapping_exons(once)) == _spans(once)


class TestUpdateCdsListEmitsNoOverlap:
    def test_the_rebuild_reconciles_before_it_installs(self, make_feature,
                                                       monkeypatch):
        """The birth site, not the write funnel: downstream consumers (the
        best-of-outcome compare, the ORF rescue) read `self.exons` directly."""
        trans = lifton_class.Lifton_TRANS.__new__(lifton_class.Lifton_TRANS)
        trans.entry = make_feature(featuretype="mRNA", start=100, end=700)
        trans.exons = [
            _exon(make_feature, 100, 300),
            _exon(make_feature, 500, 700),
        ]
        trans.entry.attributes["Parent"] = ["gene1"]
        trans._cds_attr_template = None
        trans.__dict__.setdefault("transl_table_value", 1)

        # A chained CDS whose end reaches into the second Liftoff exon.
        cds_list = [
            lifton_class.Lifton_CDS(make_feature(
                featuretype="CDS", start=100, end=250, frame="0")),
            lifton_class.Lifton_CDS(make_feature(
                featuretype="CDS", start=260, end=600, frame="0")),
        ]
        trans.update_cds_list(cds_list)

        spans = _spans(trans.exons)
        overlaps = [(spans[i], spans[i + 1]) for i in range(len(spans) - 1)
                    if spans[i + 1][0] <= spans[i][1]]
        assert not overlaps, f"update_cds_list emitted overlapping exons: {overlaps}"


class TestTransSplicedPairIsNotMerged:
    """Exons of one transcript share a strand; a pair that does not is
    trans-spliced, and an overlap there is not evidence they are one exon.

    Across CHM13, human -> zebrafish, rice, bee and drosophila, not one
    overlapping pair crossed a strand or a seqid -- so this is the case the
    reconciler must refuse rather than the case it must handle.
    """

    def test_exons_on_opposite_strands_stay_separate(self, make_feature):
        plus = lifton_class.Lifton_EXON(make_feature(
            featuretype="exon", start=100, end=300, strand="+",
            attributes={"Parent": ["tx1"]}))
        minus = lifton_class.Lifton_EXON(make_feature(
            featuretype="exon", start=250, end=400, strand="-",
            attributes={"Parent": ["tx1"]}))
        out = lifton_class.reconcile_overlapping_exons([plus, minus])
        assert len(out) == 2, "a trans-spliced pair must not be collapsed"

    def test_the_same_pair_on_one_strand_does_merge(self, make_feature):
        a = lifton_class.Lifton_EXON(make_feature(
            featuretype="exon", start=100, end=300, strand="+",
            attributes={"Parent": ["tx1"]}))
        b = lifton_class.Lifton_EXON(make_feature(
            featuretype="exon", start=250, end=400, strand="+",
            attributes={"Parent": ["tx1"]}))
        out = lifton_class.reconcile_overlapping_exons([a, b])
        assert [(e.entry.start, e.entry.end) for e in out] == [(100, 400)]

    def test_exons_on_different_sequences_stay_separate(self, make_feature):
        """Rice's `nad5` really is written this way: `exception=trans-splicing`,
        with the transcript's rows on a different seqid from its gene."""
        a = lifton_class.Lifton_EXON(make_feature(
            seqid="CP132245.1", featuretype="exon", start=100, end=300,
            attributes={"Parent": ["tx1"]}))
        b = lifton_class.Lifton_EXON(make_feature(
            seqid="CP132246.1", featuretype="exon", start=250, end=400,
            attributes={"Parent": ["tx1"]}))
        out = lifton_class.reconcile_overlapping_exons([a, b])
        assert len(out) == 2, "coordinates on two sequences must not be merged"
