"""An unusable miniprot CANDIDATE must be skipped, never take the gene with it.

`Lifton_TRANS.add_cds` now raises when a CDS spans several exons and cannot be
split at exon boundaries unambiguously -- the right call for a lifted model,
since the alternatives are a CDS written across an intron or coding bases
counted twice.

`lifton_utils.LiftOn_miniprot_alignment` builds a scaffold for each miniprot
candidate by calling `add_cds` with no guard. An exception there escapes to the
per-locus handler in Step 7, which drops the ENTIRE Liftoff gene -- even though
its DNA lift was fine and the problem was one optional miniprot hit. The same
function already documents that exact failure for a candidate whose alignment
cannot be built, and skips it; this is the same rule for the scaffold.

Measured on the corpus it never fires (no miniprot mRNA in rice, human ->
zebrafish, drosophila or dog -> cat has overlapping CDS rows), so the test is
constructed.
"""
from types import SimpleNamespace

import pytest
from intervaltree import IntervalTree

from lifton import drop_ledger, lifton_class, lifton_utils


@pytest.fixture
def make_feature(make_gffutils_feature):
    return make_gffutils_feature


class _MiniprotDb:
    """Just enough of a FeatureDB for one candidate."""

    def __init__(self, entry, cds_rows):
        self._entry = entry
        self._cds_rows = cds_rows

    def __getitem__(self, key):
        if key != self._entry.id:
            raise KeyError(key)
        return self._entry

    def children(self, feature, featuretype=None, order_by=None, **_):
        return iter(list(self._cds_rows))


def _candidate(make_feature, cds_spans):
    entry = make_feature(featuretype="mRNA", start=100, end=500,
                         attributes={"ID": ["MP1"]}, feature_id="MP1")
    rows = [make_feature(featuretype="CDS", start=s, end=e, frame="0",
                         attributes={"Parent": ["MP1"]})
            for s, e in cds_spans]
    return entry, rows


def _run(make_feature, cds_spans):
    entry, rows = _candidate(make_feature, cds_spans)
    transcript = make_feature(featuretype="mRNA", start=100, end=500,
                              attributes={"ID": ["tx1"]}, feature_id="tx1")
    tree = {"chr1": IntervalTree()}
    return lifton_utils.LiftOn_miniprot_alignment(
        "chr1", transcript, {"rna-ref": ["MP1"]}, _MiniprotDb(entry, rows),
        tree, fai=None, ref_proteins={"rna-ref": "M"}, ref_trans_id="rna-ref",
        lifton_status=lifton_class.Lifton_Status())


class TestAnAmbiguousCandidateIsSkipped:
    def test_it_does_not_raise_out_of_the_candidate_loop(self, make_feature):
        """Overlapping miniprot CDS rows make every split ambiguous."""
        drop_ledger.reset()
        aln, has_valid = _run(make_feature, [(100, 300), (250, 500)])
        assert aln is None
        assert has_valid is False, (
            "a discarded candidate must not report valid miniprot evidence")

    def test_the_rejection_is_counted(self, make_feature):
        drop_ledger.reset()
        _run(make_feature, [(100, 300), (250, 500)])
        assert drop_ledger.counts().get("cds_spanning_exons"), (
            "the rejected split must still reach the drop ledger")
