"""A CDS may belong to one exon, not to every exon it touches.

`Lifton_TRANS.add_cds` collected every exon a CDS overlapped and attached the
CDS to all of them, cloning the feature. The transcript then emitted the same
CDS two or more times -- the duplicate-CDS shape the overlapping-exon work of
the previous cycle existed to remove, arrived at from the other direction.

After that cycle this fires zero times across CHM13, human -> zebrafish,
drosophila, rice and bee (the "spans 2 exons" warning count went 8,546 / 4,468 /
3,614 / 2,632 -> 0), because the thing that used to produce a CDS straddling two
exons was miniprot's redundant `stop_codon` being ingested as a 3 bp exon. It
is still reachable on a reference CDS that genuinely spans an intron, so these
tests are constructed rather than drawn from the corpus.

The old warning also told the user "The reference model is malformed here",
which on the corpus was false -- the second exon was one LiftOn had just made.
"""
import copy

import pytest

from lifton import drop_ledger, lifton_class


@pytest.fixture
def make_feature(make_gffutils_feature):
    return make_gffutils_feature


def _trans(make_feature, exons):
    trans = lifton_class.Lifton_TRANS.__new__(lifton_class.Lifton_TRANS)
    trans.entry = make_feature(featuretype="mRNA", start=100, end=900,
                               attributes={"ID": ["tx1"], "Parent": ["gene1"]})
    trans.exons = [
        lifton_class.Lifton_EXON(make_feature(
            featuretype="exon", start=s, end=e,
            attributes={"Parent": ["tx1"]}))
        for s, e in exons
    ]
    trans._cds_attr_template = None
    return trans


def _cds(make_feature, start, end):
    return make_feature(featuretype="CDS", start=start, end=end, frame="0",
                        attributes={"ID": ["cds1"], "Parent": ["tx1"]})


class TestCdsSpanningTwoExons:
    def test_it_is_attached_once_not_to_every_exon(self, make_feature):
        """The whole point: one CDS row out, not two."""
        drop_ledger.reset()
        trans = _trans(make_feature, [(100, 300), (500, 700)])
        trans.add_cds(_cds(make_feature, 250, 550))
        carried = [e for e in trans.exons if e.cds is not None]
        assert len(carried) == 1, (
            f"the CDS was attached to {len(carried)} exons; emitting it more "
            f"than once duplicates coding sequence")

    def test_it_goes_to_the_exon_it_overlaps_most(self, make_feature):
        """250-550 shares 51 bp with exon 1 and 51 with exon 2 -- so use a
        clearly asymmetric case, where the answer is not a coin toss."""
        drop_ledger.reset()
        trans = _trans(make_feature, [(100, 300), (500, 700)])
        trans.add_cds(_cds(make_feature, 280, 690))   # 21 bp vs 191 bp
        carried = [e for e in trans.exons if e.cds is not None]
        assert [(e.entry.start, e.entry.end) for e in carried] == [(500, 700)]

    def test_the_discarded_attachment_is_counted(self, make_feature):
        drop_ledger.reset()
        trans = _trans(make_feature, [(100, 300), (500, 700)])
        trans.add_cds(_cds(make_feature, 250, 550))
        assert drop_ledger.counts().get("cds_spanning_exons"), (
            "a CDS that could not be placed in one exon must be counted, not "
            "silently duplicated or silently dropped")

    def test_a_cds_inside_one_exon_is_untouched(self, make_feature):
        """The overwhelmingly common case must not change at all."""
        drop_ledger.reset()
        trans = _trans(make_feature, [(100, 300), (500, 700)])
        entry = _cds(make_feature, 120, 280)
        trans.add_cds(entry)
        carried = [e for e in trans.exons if e.cds is not None]
        assert len(carried) == 1
        assert carried[0].cds.entry is entry, (
            "the single-exon path must attach the caller's own object, not a "
            "clone -- callers rely on that identity")
        assert drop_ledger.total() == 0

    def test_a_cds_touching_no_exon_attaches_nowhere(self, make_feature):
        drop_ledger.reset()
        trans = _trans(make_feature, [(100, 300), (500, 700)])
        trans.add_cds(_cds(make_feature, 350, 400))
        assert [e for e in trans.exons if e.cds is not None] == []
