"""A miniprot stop_codon must not be ingested as a second exon.

miniprot writes the terminal CDS of a hit *with* its stop codon, then repeats
those three bases as a nested `stop_codon` row. Ingesting both as exons gave
the transcript a 3 bp exon inside its terminal exon -- an overlapping-exon pair
the reference does not have -- and `add_cds` then overwrote the real terminal
CDS with the 3 bp one.

Measured on the regenerated CHM13 annotation: of 27 transcripts carrying an
overlapping exon pair, 11 had a 3 bp side; of 13 carrying an overlapping CDS
pair, 11 had a 3 bp side and 5 of those were an exact duplicate row. All 3,963
miniprot stop_codon rows in that run are fully inside a sibling CDS.
"""
from types import SimpleNamespace

from lifton import coreutils


def _row(ftype, start, end):
    return SimpleNamespace(featuretype=ftype, start=start, end=end,
                           id=f"{ftype}-{start}")


class TestDropRedundantStopCodons:
    def test_a_stop_codon_inside_its_cds_is_dropped(self):
        """miniprot's own shape: stop_codon 1135921-1135923 in CDS ...-1135923."""
        rows = [_row("CDS", 1135793, 1135923), _row("stop_codon", 1135921, 1135923)]
        kept = coreutils.drop_redundant_stop_codons(rows)
        assert [r.featuretype for r in kept] == ["CDS"]

    def test_a_stop_codon_at_the_cds_start_is_dropped(self):
        """The minus-strand form: stop_codon 823466-823468 in CDS 823466-823506."""
        rows = [_row("CDS", 823466, 823506), _row("stop_codon", 823466, 823468)]
        assert [r.featuretype for r in
                coreutils.drop_redundant_stop_codons(rows)] == ["CDS"]

    def test_a_stop_codon_outside_every_cds_is_kept(self):
        """The GTF convention, where the CDS stops short of the stop codon.

        That row carries three bases the CDS list does not have, so dropping it
        would shorten the coding sequence. This is why the filter tests
        containment instead of just removing the featuretype.
        """
        rows = [_row("CDS", 100, 199), _row("stop_codon", 200, 202)]
        assert [r.featuretype for r in
                coreutils.drop_redundant_stop_codons(rows)] == ["CDS", "stop_codon"]

    def test_a_partially_overlapping_stop_codon_is_kept(self):
        rows = [_row("CDS", 100, 200), _row("stop_codon", 199, 201)]
        assert len(coreutils.drop_redundant_stop_codons(rows)) == 2

    def test_order_and_identity_are_otherwise_untouched(self):
        rows = [_row("CDS", 10, 20), _row("CDS", 30, 40), _row("CDS", 50, 60)]
        kept = coreutils.drop_redundant_stop_codons(rows)
        assert kept == rows

    def test_a_stop_codon_with_no_cds_sibling_survives(self):
        rows = [_row("stop_codon", 10, 12)]
        assert len(coreutils.drop_redundant_stop_codons(rows)) == 1


class TestMiniprotCandidateScaffold:
    """The call site: the chaining's miniprot candidate must not gain a 3 bp exon."""

    def test_the_candidate_builder_filters_its_children(self):
        import inspect
        from lifton import lifton_utils
        source = inspect.getsource(lifton_utils.LiftOn_miniprot_alignment)
        assert "drop_redundant_stop_codons" in source, (
            "the miniprot candidate scaffold still ingests every "
            "('CDS','stop_codon') child as both an exon and a CDS"
        )
