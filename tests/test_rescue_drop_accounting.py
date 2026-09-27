"""The miniprot-only rescue must count what it abandons for a missing lookup.

The rescue is default-ON and abandons candidates in roughly forty places. It
contained no `drop_ledger` call at all, so five whole genomes reported a drop
total of 0 -- which reads as "nothing was lost" and meant "nothing was counted".

About twenty-five of those sites are deliberate filters (identity floor, overlap
suppression, length band, dedup, processed-pseudogene) and must STAY uncounted:
they are decisions, not losses. The rest are lookups that failed.

After wiring, the five genomes still report 0 -- on this corpus every rescue
rejection really is a decision. These tests are what make that 0 mean something:
they show the counters fire when there IS something to count.
"""
import pytest

from lifton import drop_ledger, lifton_utils


class TestUnresolvedMiniprotHit:
    def test_a_hit_absent_from_the_id_map_is_counted(self):
        """`get_ref_ids_miniprot` returns (None, None) for this."""
        drop_ledger.reset()
        abandoned = lifton_utils.record_unresolved_miniprot_hit(
            None, None, "MP000123")
        assert abandoned is True
        assert drop_ledger.counts().get("miniprot_hit_unmapped") == 1
        assert drop_ledger.examples("miniprot_hit_unmapped") == ["MP000123"]

    def test_a_known_transcript_with_no_gene_is_counted_separately(self):
        """(None, ref_trans_id): the id map knew it, the reverse index did not.

        This is the case the organellar reference-index gap produced, and the
        two were previously collapsed into one `continue` by all five callers.
        """
        drop_ledger.reset()
        abandoned = lifton_utils.record_unresolved_miniprot_hit(
            None, "rna-OrsajCp058", "MP000123")
        assert abandoned is True
        assert drop_ledger.counts().get("miniprot_gene_unresolved") == 1
        assert drop_ledger.counts().get("miniprot_hit_unmapped") is None

    def test_a_resolved_hit_is_not_counted_and_not_abandoned(self):
        drop_ledger.reset()
        abandoned = lifton_utils.record_unresolved_miniprot_hit(
            "gene-X", "rna-X", "MP000123")
        assert abandoned is False
        assert drop_ledger.total() == 0


class TestReferenceLengthLookup:
    def test_a_missing_entry_is_counted(self):
        drop_ledger.reset()
        assert lifton_utils.reference_length_or_none({}, "gene-X") is None
        assert drop_ledger.counts().get("reference_feature_length_missing") == 1

    def test_a_genuine_zero_is_not_counted(self):
        """A reference gene with no CDS rows is a correct decision, not a loss.

        Both were consumed by `if not ref_len: continue`, so a lookup failure
        and a non-coding gene were indistinguishable.
        """
        drop_ledger.reset()
        assert lifton_utils.reference_length_or_none({"gene-X": 0}, "gene-X") == 0
        assert drop_ledger.total() == 0

    def test_a_real_length_is_returned_untouched(self):
        drop_ledger.reset()
        assert lifton_utils.reference_length_or_none(
            {"gene-X": 1234}, "gene-X") == 1234
        assert drop_ledger.total() == 0


class TestEveryClassIsExplainable:
    def test_the_rescue_classes_are_declared_with_a_sentence(self):
        """`drop_ledger.record` raises for an undeclared class on purpose: a
        counter nobody can interpret is barely better than no counter."""
        for kind in ("miniprot_hit_unmapped", "miniprot_gene_unresolved",
                     "reference_feature_length_missing",
                     "rescue_candidate_error", "cds_spanning_exons"):
            assert kind in drop_ledger.CLASSES
            assert len(drop_ledger.CLASSES[kind]) > 30

    def test_an_undeclared_class_is_refused(self):
        with pytest.raises(KeyError):
            drop_ledger.record("something_nobody_declared")


class TestTheWiringIsLive:
    """Guard against the failure this whole change is about: a counter that
    reports 0 because its call site is unreachable."""

    @pytest.mark.parametrize("module_name, expected", [
        ("lifton.miniprot_rescue", 3),
        ("lifton.run_miniprot", 1),
    ])
    def test_the_resolution_guard_is_called_where_candidates_are_abandoned(
            self, module_name, expected):
        import importlib
        import inspect
        source = inspect.getsource(importlib.import_module(module_name))
        assert source.count("record_unresolved_miniprot_hit") >= expected, (
            f"{module_name} abandons miniprot candidates in more places than "
            f"it counts them")

    def test_the_rescue_records_drops_at_all(self):
        import inspect
        from lifton import miniprot_rescue
        source = inspect.getsource(miniprot_rescue)
        assert "drop_ledger.record" in source, (
            "the rescue pass reports no losses because it counts none")
