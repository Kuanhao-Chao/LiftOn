"""A spanning CDS is split at exon boundaries, never written across an intron."""
import io

import pytest
from Bio.Seq import Seq

from lifton import drop_ledger, lifton_class
from lifton.gff3_validator import validate_gff3_file


@pytest.fixture
def make_feature(make_gffutils_feature):
    return make_gffutils_feature


def _trans(make_feature, exons, strand="+"):
    trans = lifton_class.Lifton_TRANS.__new__(lifton_class.Lifton_TRANS)
    trans.entry = make_feature(featuretype="mRNA", start=100, end=900, strand=strand,
                               attributes={"ID": ["tx1"], "Parent": ["gene1"]})
    trans.exons = [
        lifton_class.Lifton_EXON(make_feature(
            featuretype="exon", start=s, end=e, strand=strand,
            attributes={"ID": [f"exon{i}"], "Parent": ["tx1"]}))
        for i, (s, e) in enumerate(exons, 1)
    ]
    trans._cds_attr_template = None
    return trans


def _cds(make_feature, start, end, strand="+", frame="0"):
    return make_feature(featuretype="CDS", start=start, end=end, strand=strand, frame=frame,
                        attributes={"ID": ["cds1"], "Parent": ["tx1"]})


class TestCdsSpanningTwoExons:
    def test_split_serializes_only_exonic_bases_and_validates(self, make_feature, tmp_path):
        drop_ledger.reset()
        trans = _trans(make_feature, [(100, 300), (500, 700)])
        trans.add_cds(_cds(make_feature, 250, 550))
        assert [(e.cds.entry.start, e.cds.entry.end, e.cds.entry.frame)
                for e in trans.exons] == [(250, 300, "0"), (500, 550, "0")]
        assert drop_ledger.total() == 0

        trans.normalize_containment()
        assert [(e.entry.start, e.entry.end) for e in trans.exons] == [(100, 300), (500, 700)]
        buffer = io.StringIO()
        assert trans.write_entry(buffer)
        path = tmp_path / "split.gff3"
        path.write_text("##gff-version 3\nchr1\ttest\tgene\t100\t900\t.\t+\t.\tID=gene1\n" + buffer.getvalue())
        assert validate_gff3_file(path).errors == []
        cds_rows = [line.split("\t") for line in buffer.getvalue().splitlines()
                    if "\tCDS\t" in line]
        assert [(int(row[3]), int(row[4])) for row in cds_rows] == [(250, 300), (500, 550)]
        assert all("ID=cds1" in row[8] for row in cds_rows)

        # An independent synthetic target sequence makes intron inclusion or
        # double-counting observable at the amino-acid level.
        genome = ["C"] * 900
        first = "ATG" + "GCT" * 16
        second = "GCT" * 16 + "TAA"
        genome[249:300] = first
        genome[499:550] = second
        sequence = "".join(genome)
        coding = "".join(sequence[int(row[3])-1:int(row[4])] for row in cds_rows)
        assert str(Seq(coding).translate()) == "M" + "A" * 32 + "*"

    def test_minus_strand_phase_follows_transcript_order(self, make_feature):
        drop_ledger.reset()
        trans = _trans(make_feature, [(100, 300), (500, 700)], strand="-")
        trans.add_cds(_cds(make_feature, 250, 552, strand="-"))
        assert [(e.cds.entry.start, e.cds.entry.end, e.cds.entry.frame)
                for e in trans.exons] == [(250, 300, "1"), (500, 552, "0")]

    def test_ambiguous_overlapping_exons_are_counted_and_rejected(self, make_feature):
        drop_ledger.reset()
        trans = _trans(make_feature, [(100, 300), (290, 700)])
        with pytest.raises(ValueError, match="cannot be split unambiguously"):
            trans.add_cds(_cds(make_feature, 250, 550))
        assert drop_ledger.counts()["cds_spanning_exons"] == 1
        assert all(exon.cds is None for exon in trans.exons)

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
