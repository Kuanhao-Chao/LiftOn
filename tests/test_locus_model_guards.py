"""Guards for two assumptions the locus model made but never checked.

* Lifton_TRANS.add_cds assumed "a CDS lies in exactly one exon". On a
  malformed reference the CDS was attached to every exon it overlapped -- and
  every one of them received the SAME Feature object, so a later edit to one
  copy silently rewrote the others.
* _MFeatureDbProxy.children ignored order_by and always served the one
  ordering its cache was built with, so a caller asking for file order (which
  is what the real backends return) would have been handed start-sorted rows
  instead. Every current caller asks for 'start', so this was latent.
"""
import pytest

from lifton.lifton_class import Lifton_TRANS, Lifton_EXON


# ── a CDS that spans several exons ──────────────────────────────────────────

class TestCdsSpanningSeveralExons:
    """`add_cds` assumed a CDS lies in exactly one exon and never checked it.

    On a malformed reference the CDS was attached to every exon it touched --
    and all of them received the SAME Feature object, so a later edit to one
    copy silently rewrote the others.
    """

    @staticmethod
    def _trans(make_feature, spans):
        trans = Lifton_TRANS.__new__(Lifton_TRANS)
        trans.entry = make_feature(
            featuretype="mRNA", start=spans[0][0], end=spans[-1][1],
            feature_id="rna-T", attributes={"ID": ["rna-T"]})
        trans.entry.id = "rna-T"
        trans.exons = [
            Lifton_EXON(make_feature(
                featuretype="exon", start=s, end=e, feature_id=f"exon-{i}",
                attributes={"ID": [f"exon-{i}"], "Parent": ["rna-T"]}))
            for i, (s, e) in enumerate(spans, 1)
        ]
        return trans

    def test_one_exon_is_the_normal_case(self, make_gffutils_feature):
        trans = self._trans(make_gffutils_feature, [(100, 200), (300, 400)])
        cds = make_gffutils_feature(featuretype="CDS", start=120, end=180,
                                    feature_id="cds-1",
                                    attributes={"ID": ["cds-1"]})
        trans.add_cds(cds)
        assert trans.exons[0].cds is not None
        assert trans.exons[1].cds is None
        # Unchanged for the single-exon case: the caller's object is used.
        assert trans.exons[0].cds.entry is cds

    def test_a_spanning_cds_gives_each_exon_its_own_object(self, make_gffutils_feature):
        trans = self._trans(make_gffutils_feature, [(100, 200), (300, 400)])
        cds = make_gffutils_feature(featuretype="CDS", start=150, end=350,
                                    feature_id="cds-1",
                                    attributes={"ID": ["cds-1"]})
        trans.add_cds(cds)

        first, second = trans.exons[0].cds, trans.exons[1].cds
        assert first is not None and second is not None
        assert first.entry is not second.entry

        first.entry.attributes["ID"] = ["rewritten"]
        assert second.entry.attributes["ID"] == ["cds-1"]

    def test_a_spanning_cds_is_reported(self, make_gffutils_feature, monkeypatch):
        from lifton import lifton_class
        warnings = []
        monkeypatch.setattr(lifton_class.logger, "log_warning", warnings.append)

        trans = self._trans(make_gffutils_feature, [(100, 200), (300, 400)])
        cds = make_gffutils_feature(featuretype="CDS", start=150, end=350,
                                    feature_id="cds-1",
                                    attributes={"ID": ["cds-1"]})
        trans.add_cds(cds)
        assert any("spans 2 exons" in message for message in warnings), warnings

    def test_the_normal_case_is_silent(self, make_gffutils_feature, monkeypatch):
        from lifton import lifton_class
        warnings = []
        monkeypatch.setattr(lifton_class.logger, "log_warning", warnings.append)

        trans = self._trans(make_gffutils_feature, [(100, 200), (300, 400)])
        trans.add_cds(make_gffutils_feature(
            featuretype="CDS", start=120, end=180, feature_id="cds-1",
            attributes={"ID": ["cds-1"]}))
        assert warnings == []


# ── the miniprot child proxy honours the requested ordering ─────────────────

class TestMiniprotProxyOrdering:
    @staticmethod
    def _children(make_feature):
        rows = []
        for index, (start, file_order) in enumerate(
                ((300, 2), (100, 1), (500, 3)), 1):
            feature = make_feature(featuretype="CDS", start=start,
                                   end=start + 50, feature_id=f"c{index}",
                                   attributes={"ID": [f"c{index}"]})
            feature.file_order = file_order
            rows.append(feature)
        return sorted(rows, key=lambda f: f.start)   # the cached order

    def test_start_order_is_served_from_the_cache(self, make_gffutils_feature):
        from lifton.locus_pipeline import _MFeatureDbProxy
        rows = self._children(make_gffutils_feature)
        assert _MFeatureDbProxy._ordered(rows, "start") is rows

    def test_file_order_is_re_sorted_not_silently_start_sorted(
            self, make_gffutils_feature):
        from lifton.locus_pipeline import _MFeatureDbProxy
        rows = self._children(make_gffutils_feature)
        ordered = _MFeatureDbProxy._ordered(rows, None)
        assert [f.file_order for f in ordered] == [1, 2, 3]
        assert [f.start for f in ordered] == [100, 300, 500]

    def test_rows_without_an_ingest_position_sort_last(self, make_gffutils_feature):
        from lifton.locus_pipeline import _MFeatureDbProxy
        rows = self._children(make_gffutils_feature)
        rows[0].file_order = None
        ordered = _MFeatureDbProxy._ordered(rows, None)
        assert ordered[-1] is rows[0]

    def test_an_unservable_ordering_is_loud(self, make_gffutils_feature):
        from lifton.locus_pipeline import _MFeatureDbProxy
        rows = self._children(make_gffutils_feature)
        with pytest.raises(NotImplementedError, match="order_by"):
            _MFeatureDbProxy._ordered(rows, "length")
