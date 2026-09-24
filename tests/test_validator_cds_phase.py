"""`cds_phase_consistency` must read a 5'-partial CDS's own first phase.

GFF3's phase is the number of bases to skip at the start of a segment to reach
the next codon. The check computed every later segment's expected phase as if
the first segment began on a codon (phase 0), so a model whose first CDS row
legitimately starts mid-codon -- `start_lost`, or any 5'-partial lift -- had
every later segment flagged. v1.0.14 keeps that initial phase (`9bdf8a1`), and
the check then reported 5,648 warnings on dog -> cat and 1,465 on human ->
gorilla against 26 and 0 for v1.0.13; read with the first phase, the same
v1.0.14 files have no inconsistent transcript at all.
"""
import textwrap

from lifton import gff3_validator


def _write(tmp_path, body):
    path = tmp_path / "out.gff3"
    path.write_text("##gff-version 3\n" + textwrap.dedent(body).lstrip())
    return str(path)


def _phase_issues(path):
    result = gff3_validator.validate_gff3_file(path)
    return [i for i in result.issues if i.check == "cds_phase_consistency"]


def _model(strand, cds_rows):
    rows = [f"chr1\tLiftOn\tgene\t100\t450\t.\t{strand}\t.\tID=g1\n",
            f"chr1\tLiftOn\tmRNA\t100\t450\t.\t{strand}\t.\tID=tx1;Parent=g1\n"]
    for n, (start, end, phase) in enumerate(cds_rows, 1):
        rows.append(f"chr1\tLiftOn\texon\t{start}\t{end}\t.\t{strand}\t.\tID=e{n};Parent=tx1\n")
    for start, end, phase in cds_rows:
        rows.append(f"chr1\tLiftOn\tCDS\t{start}\t{end}\t.\t{strand}\t{phase}\tID=c1;Parent=tx1\n")
    return "".join(rows)


class TestAFivePrimePartialModel:
    def test_plus_strand(self, tmp_path):
        # 51 bp in phase 2 leaves 49 coding bases: one over a codon, so the
        # next segment skips 2; 152 bp in, 150 coding bases, phase 0.
        path = _write(tmp_path, _model("+", [(100, 150, 2), (200, 300, 2), (400, 450, 0)]))
        assert _phase_issues(path) == []

    def test_minus_strand(self, tmp_path):
        # 5' is the highest coordinate: 400-450 first, phase 1 -> 50 coding
        # bases, two over, next phase 1; then 151 coding bases, next phase 2.
        path = _write(tmp_path, _model("-", [(100, 150, 2), (200, 300, 1), (400, 450, 1)]))
        assert _phase_issues(path) == []


class TestAnInconsistentPhaseStillWarns:
    def test_complete_model_with_a_wrong_phase(self, tmp_path):
        path = _write(tmp_path, _model("+", [(100, 150, 0), (200, 300, 1), (400, 450, 0)]))
        assert _phase_issues(path)

    def test_partial_model_with_a_wrong_phase(self, tmp_path):
        """Reading the first phase must not excuse a later segment."""
        path = _write(tmp_path, _model("+", [(100, 150, 2), (200, 300, 0), (400, 450, 0)]))
        assert _phase_issues(path)
