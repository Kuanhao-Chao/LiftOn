"""The deep validator streams one top-level block at a time, reporting identically.

`validate_gff3_file` materialised a `GFF3Record` per row and then walked that list
five times. Measured: **406 MB** for a 92,675-row output where the streaming
structural gate uses 36 MB, and **5.64 GB** for a full dog->cat lift (1,815,199
rows). It is reachable from `--validate-output` and the `gff3-validate` CLI, so a
mammalian genome could exhaust memory just to be checked.

Streaming is exact rather than approximate because LiftOn emits each top-level
feature immediately followed by its descendants -- verified on that full lift:
34,872 blocks, zero child rows whose Parent lies outside their block. The tests
below pin the three things that make the report survive being computed per block:
the per-check caps stay global, the issues come back in the whole-file order, and
a file that is NOT block-ordered falls back to the unchanged path.
"""
import ast
import pathlib
import re

import pytest

from lifton import gff3_validator as v


HEADER = "##gff-version 3\n"


def _write(tmp_path, name, rows):
    path = tmp_path / name
    path.write_text(HEADER + "".join(rows))
    return str(path)


def _row(ftype, start, end, attrs, seqid="chr1", strand="+", phase="."):
    return f"{seqid}\tLiftOn\t{ftype}\t{start}\t{end}\t.\t{strand}\t{phase}\t{attrs}\n"


def _gene_block(index, *, seqid="chr1", broken=None):
    """A complete gene -> mRNA -> exon + CDS block."""
    g, r = f"gene{index}", f"rna{index}"
    base = index * 1000
    rows = [
        _row("gene", base + 1, base + 900, f"ID={g}", seqid=seqid),
        _row("mRNA", base + 1, base + 900, f"ID={r};Parent={g}", seqid=seqid),
        _row("exon", base + 1, base + 300, f"ID=ex{index};Parent={r}", seqid=seqid),
    ]
    if broken != "no_cds":
        rows.append(_row("CDS", base + 1, base + 300,
                         f"ID=cds{index};Parent={r}", seqid=seqid, phase="0"))
    return rows


def _both_paths(path, **kwargs):
    """Return (streaming_result, whole_file_result) for the same input."""
    streamed = v.validate_gff3_file(path, **kwargs)
    original = v._validate_streaming
    try:
        v._validate_streaming = lambda *a, **k: None      # force the old path
        whole = v.validate_gff3_file(path, **kwargs)
    finally:
        v._validate_streaming = original
    return streamed, whole


def _fingerprint(result):
    return (
        [(str(i.severity), i.lineno, i.feature_id, i.check, i.message)
         for i in result.issues],
        dict(result.issue_totals),
        {str(k): val for k, val in result.severity_totals.items()},
        dict(result.stats),
        result.total_lines, result.data_lines, result.comment_lines,
        result.is_valid,
    )


# ── the two paths must agree, issue for issue ───────────────────────────────

class TestEquivalence:
    def test_clean_multi_gene_file(self, tmp_path):
        rows = [r for i in range(1, 6) for r in _gene_block(i)]
        path = _write(tmp_path, "clean.gff3", rows)
        streamed, whole = _both_paths(path)
        assert streamed.is_valid and whole.is_valid
        assert _fingerprint(streamed) == _fingerprint(whole)

    def test_file_exercising_several_rules_at_once(self, tmp_path):
        """Mixes rule groups so the per-block interleaving would show up."""
        rows = []
        rows += _gene_block(1)
        rows += _gene_block(2, broken="no_cds")          # mrna_has_cds
        rows += [_row("gene", 3000, 3900, "ID=gene3")]   # gene_has_transcripts
        rows += _gene_block(4, broken="no_cds")          # mrna_has_cds again
        rows += [_row("gene", 5000, 5900, "ID=gene5")]   # and again
        rows += _gene_block(6)
        path = _write(tmp_path, "mixed.gff3", rows)

        streamed, whole = _both_paths(path)
        assert len(streamed.issues) > 3, "fixture must actually produce issues"
        assert _fingerprint(streamed) == _fingerprint(whole)

    def test_issue_order_is_grouped_the_way_the_whole_file_path_groups_it(self, tmp_path):
        """Per-block emission interleaves rule groups; the sort undoes that."""
        rows = []
        for index in (1, 2, 3, 4):
            rows += _gene_block(index, broken="no_cds")
            rows += [_row("gene", index * 1000 + 500, index * 1000 + 600,
                          f"ID=lonely{index}")]
        path = _write(tmp_path, "grouped.gff3", rows)
        streamed, whole = _both_paths(path)
        assert [i.check for i in streamed.issues] == [i.check for i in whole.issues]

    def test_a_duplicate_id_falls_back_rather_than_reporting_differently(self, tmp_path):
        """A real duplicate ID cannot be emulated once blocks are discarded.

        The whole-file path keeps only the FIRST record under an id, so a later
        duplicate is invisible to every lookup and to the per-gene loop. Per
        block that record is long gone, and the duplicate would pick up an extra
        `gene_has_transcripts`. Invalid GFF3 either way, so the file goes back to
        the unchanged path instead of getting a different report.
        """
        from lifton.gff3_validator import ValidationResult
        rows = _gene_block(1) + _gene_block(2)
        rows.append(_row("gene", 9000, 9100, "ID=gene1"))   # collides with block 1
        path = _write(tmp_path, "dup.gff3", rows)

        assert v._validate_streaming(
            path, ValidationResult(file_path=path),
            True, True, True, True, 50) is None

        streamed, whole = _both_paths(path)
        assert any(i.check == "duplicate_id" for i in streamed.issues)
        assert _fingerprint(streamed) == _fingerprint(whole)

    def test_a_discontinuous_cds_is_not_a_duplicate(self, tmp_path):
        """Segments of one multi-exon CDS share an id legitimately, and must
        NOT cost the bounded path."""
        from lifton.gff3_validator import ValidationResult
        rows = [
            _row("gene", 1, 900, "ID=g"),
            _row("mRNA", 1, 900, "ID=r;Parent=g"),
            _row("exon", 1, 300, "ID=e1;Parent=r"),
            _row("exon", 601, 900, "ID=e2;Parent=r"),
            _row("CDS", 1, 300, "ID=c;Parent=r", phase="0"),
            _row("CDS", 601, 900, "ID=c;Parent=r", phase="0"),
        ]
        path = _write(tmp_path, "disc.gff3", rows)
        assert v._validate_streaming(
            path, ValidationResult(file_path=path),
            True, True, True, True, 50) is not None
        streamed, whole = _both_paths(path)
        assert _fingerprint(streamed) == _fingerprint(whole)

    @pytest.mark.parametrize("cap", [0, 1, 2, 50])
    def test_the_per_check_cap_stays_global(self, tmp_path, cap):
        """Capping per block would let a long file emit far more issues."""
        rows = [r for i in range(1, 9) for r in _gene_block(i, broken="no_cds")]
        path = _write(tmp_path, "capped.gff3", rows)
        streamed, whole = _both_paths(path, max_issues_per_check=cap)
        assert len(streamed.issues) == len(whole.issues)
        assert sum(1 for i in streamed.issues if i.check == "mrna_has_cds") <= cap
        assert _fingerprint(streamed) == _fingerprint(whole)

    def test_empty_and_directive_only_files(self, tmp_path):
        path = _write(tmp_path, "bare.gff3", [])
        streamed, whole = _both_paths(path)
        assert _fingerprint(streamed) == _fingerprint(whole)


# ── the fallback ────────────────────────────────────────────────────────────

class TestNonContiguousFallback:
    @staticmethod
    def _interleaved(tmp_path):
        """A child that refers back to an EARLIER block."""
        rows = _gene_block(1) + _gene_block(2)
        rows.append(_row("exon", 400, 500, "ID=late;Parent=rna1"))
        return _write(tmp_path, "interleaved.gff3", rows)

    def test_streaming_declines_and_returns_none(self, tmp_path):
        from lifton.gff3_validator import ValidationResult
        path = self._interleaved(tmp_path)
        assert v._validate_streaming(
            path, ValidationResult(file_path=path),
            True, True, True, True, 50) is None

    def test_the_result_is_the_whole_file_result(self, tmp_path):
        path = self._interleaved(tmp_path)
        streamed, whole = _both_paths(path)
        assert _fingerprint(streamed) == _fingerprint(whole)

    def test_contiguous_files_do_NOT_fall_back(self, tmp_path):
        """Otherwise the bounded path would silently never run."""
        from lifton.gff3_validator import ValidationResult
        path = _write(tmp_path, "ok.gff3",
                      [r for i in range(1, 4) for r in _gene_block(i)])
        assert v._validate_streaming(
            path, ValidationResult(file_path=path),
            True, True, True, True, 50) is not None


# ── the generator refactor ──────────────────────────────────────────────────

def test_parser_generator_matches_the_list_builder(tmp_path):
    path = _write(tmp_path, "parse.gff3",
                  [r for i in range(1, 4) for r in _gene_block(i)]
                  + [_row("CDS", 10, 20, "ID=x;Parent=nobody", phase="9")])
    records, issues, meta = v._parse_gff3(path, 50)

    parser = v._GFF3Parser(path, 50)
    streamed = list(parser)
    assert [r.lineno for r in streamed] == [r.lineno for r in records]
    assert [(i.check, i.lineno, i.message) for i in parser.issues] == \
           [(i.check, i.lineno, i.message) for i in issues]
    assert parser.meta == meta


# ── the grouping constant must track the source ─────────────────────────────

def test_hierarchy_groups_match_the_functions_loop_structure():
    """`_HIERARCHY_ISSUE_GROUPS` encodes which loop emits each check name.

    Re-derive it from the source so a rule added to a different loop fails here
    rather than silently reordering every report.
    """
    source = pathlib.Path(v.__file__).with_suffix(".py").read_text()
    lines = source.split("\n")
    node = next(n for n in ast.walk(ast.parse(source))
                if isinstance(n, ast.FunctionDef) and n.name == "_check_hierarchy")

    derived, group = {}, 0
    for stmt in node.body:
        if not isinstance(stmt, ast.For):
            continue
        segment = "\n".join(lines[stmt.lineno - 1:stmt.end_lineno])
        names = re.findall(
            r'GFF3Issue\(\s*Severity\.\w+,[^,]+,[^,]+,\s*"([a-z0-9_]+)"', segment)
        if not names:
            continue
        for name in dict.fromkeys(names):
            derived.setdefault(name, group)
        group += 1

    assert derived == v._HIERARCHY_ISSUE_GROUPS, (
        "the loop structure of _check_hierarchy changed; update "
        "_HIERARCHY_ISSUE_GROUPS or the streaming report order will drift"
    )
