"""GitHub issue #50 — write the run summary as a file in the output directory.

``stats.print_report`` printed the summary block only to stderr; the per-category
side-car files (mapped_feature.txt, …) existed but not the summary itself. It now
also writes ``stats/summary.txt`` — derived from a side-car writer's path (no new
argument or caller change) and teed from the same stderr block, so the two can
never drift.
"""
import os

from lifton import stats


def _w(d, n):
    return open(os.path.join(str(d), n), "w")


def _run(tmp_path):
    fmf = _w(tmp_path, "mapped_feature.txt")
    fu = _w(tmp_path, "unmapped_features.txt")
    fe = _w(tmp_path, "extra_copy_features.txt")
    fmt = _w(tmp_path, "mapped_transcript.txt")
    fft = _w(tmp_path, "completeness_by_feature_type.txt")
    stats.print_report({}, {}, fu, fe, fmf, fmt, debug=False, fw_feature_type=fft)
    for fh in (fmf, fu, fe, fmt, fft):
        fh.close()


def test_print_report_writes_summary_txt(tmp_path):
    _run(tmp_path)
    summary = tmp_path / "summary.txt"
    assert summary.exists()
    text = summary.read_text()
    assert "Total features in reference" in text
    assert "Total features in reference\t\t: 0" in text
    assert "Completeness by feature type" in text
    assert text.rstrip().endswith("*********************************************")


def test_summary_file_mirrors_stderr(tmp_path, capsys):
    _run(tmp_path)
    err = capsys.readouterr().err
    file_text = (tmp_path / "summary.txt").read_text()
    # The block still reaches stderr (tee), and the file holds the same content.
    assert "Total features in reference" in err
    assert "Total features in reference" in file_text
