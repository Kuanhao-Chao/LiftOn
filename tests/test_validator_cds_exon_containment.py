"""A CDS must be transcribed by an exon, not merely inside the mRNA envelope."""

import pytest

from lifton import gff3_validator


def _annotation(cds_rows, exon_rows=True):
    rows = [
        "##gff-version 3",
        "chr1\ttest\tgene\t100\t700\t.\t+\t.\tID=g1",
        "chr1\ttest\tmRNA\t100\t700\t.\t+\t.\tID=t1;Parent=g1",
    ]
    if exon_rows:
        rows += [
            "chr1\ttest\texon\t100\t300\t.\t+\t.\tID=e1;Parent=t1",
            "chr1\ttest\texon\t500\t700\t.\t+\t.\tID=e2;Parent=t1",
        ]
    rows += [f"chr1\ttest\tCDS\t{start}\t{end}\t.\t+\t{phase}\tID=c1;Parent=t1"
             for start, end, phase in cds_rows]
    return "\n".join(rows) + "\n"


@pytest.mark.parametrize("streaming", [True, False])
def test_cds_crossing_intron_is_an_error(tmp_path, monkeypatch, streaming):
    if not streaming:
        monkeypatch.setattr(gff3_validator, "_validate_streaming", lambda *args: None)
    path = tmp_path / "crossing.gff3"
    path.write_text(_annotation([(250, 550, 0)]))
    result = gff3_validator.validate_gff3_file(path)
    assert [issue.check for issue in result.errors] == ["cds_exon_containment"]
    assert not result.is_valid


@pytest.mark.parametrize("streaming", [True, False])
def test_discontinuous_exonic_cds_is_valid(tmp_path, monkeypatch, streaming):
    if not streaming:
        monkeypatch.setattr(gff3_validator, "_validate_streaming", lambda *args: None)
    path = tmp_path / "split.gff3"
    path.write_text(_annotation([(250, 300, 0), (500, 550, 0)]))
    result = gff3_validator.validate_gff3_file(path)
    assert result.errors == []


def test_sparse_cds_without_explicit_exon_is_allowed(tmp_path):
    path = tmp_path / "sparse.gff3"
    path.write_text(_annotation([(250, 300, 0)], exon_rows=False))
    result = gff3_validator.validate_gff3_file(path)
    assert not any(issue.check == "cds_exon_containment" for issue in result.errors)
