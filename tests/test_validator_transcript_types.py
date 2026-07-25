"""Tier-1 audit fix — the mandatory structural validator must not reject LiftOn's
own structurally-valid output.

`validate_gff3_structure` runs unconditionally before publication and every issue it
raises is an ERROR that aborts the run. It recognised transcript levels via a closed
``TRANSCRIPT_TYPES`` allowlist, but the LIFT path auto-detects gene-like parent types
(`lifton_utils.get_gene_like_feature_types`) and preserves the reference's middle-level
type verbatim. Ensembl/GENCODE use ``pseudogene -> pseudogenic_transcript -> exon``, and
``pseudogenic_transcript`` is not in the allowlist — so LiftOn faithfully lifted the
feature, wrote valid GFF3, then rejected its own output and exited 2 with no annotation.

Recognition is now STRUCTURAL: a feature that itself sits under a gene-like root is a
transcript level whatever SO name the reference used.
"""
import pytest

from lifton.gff3_validator import validate_gff3_file, validate_gff3_structure


def _write(tmp_path, rows, name="out.gff3"):
    path = tmp_path / name
    path.write_text("##gff-version 3\n" + "\n".join(rows) + "\n", encoding="utf-8")
    return str(path)


def _row(ftype, start, end, fid, parent=None, phase="."):
    attrs = f"ID={fid}" + (f";Parent={parent}" if parent else "")
    return f"chr1\tLiftOn\t{ftype}\t{start}\t{end}\t.\t+\t{phase}\t{attrs}"


@pytest.mark.parametrize("transcript_type", [
    "pseudogenic_transcript",   # Ensembl/GENCODE pseudogene child (the reported case)
    "pseudogenic_tRNA",
    "misc_RNA",
    "pre_miRNA",
    "unconfirmed_transcript",
])
def test_exon_under_any_reference_transcript_type_is_valid(tmp_path, transcript_type):
    """Any middle level under a gene-like root is a valid exon parent."""
    parent_gene = "pseudogene" if transcript_type.startswith("pseudogenic") else "gene"
    path = _write(tmp_path, [
        _row(parent_gene, 100, 400, "g1"),
        _row(transcript_type, 100, 400, "t1", parent="g1"),
        _row("exon", 100, 199, "e1", parent="t1"),
        _row("exon", 300, 400, "e2", parent="t1"),
    ])
    result = validate_gff3_structure(path)
    assert result.is_valid, [str(e) for e in result.errors]


def test_cds_under_unlisted_transcript_type_is_valid(tmp_path):
    path = _write(tmp_path, [
        _row("gene", 100, 400, "g1"),
        _row("misc_RNA", 100, 400, "t1", parent="g1"),
        _row("exon", 100, 199, "e1", parent="t1"),
        _row("CDS", 100, 199, "c1", parent="t1", phase="0"),
    ])
    result = validate_gff3_structure(path)
    assert result.is_valid, [str(e) for e in result.errors]


def test_known_transcript_types_still_valid(tmp_path):
    """Control: the standard mRNA hierarchy is unaffected."""
    path = _write(tmp_path, [
        _row("gene", 100, 400, "g1"),
        _row("mRNA", 100, 400, "t1", parent="g1"),
        _row("exon", 100, 199, "e1", parent="t1"),
        _row("CDS", 100, 199, "c1", parent="t1", phase="0"),
    ])
    result = validate_gff3_structure(path)
    assert result.is_valid, [str(e) for e in result.errors]


def test_deep_validator_also_accepts_unlisted_transcript_types(tmp_path):
    """--validate-output uses validate_gff3_file, whose exon/CDS parent rules had the
    same closed allowlist -- and its failures abort the run too."""
    path = _write(tmp_path, [
        _row("pseudogene", 100, 400, "g1"),
        _row("pseudogenic_transcript", 100, 400, "t1", parent="g1"),
        _row("exon", 100, 199, "e1", parent="t1"),
    ])
    result = validate_gff3_file(path)
    parent_checks = {
        issue.check for issue in result.errors
        if issue.check in {"exon_parent_type", "cds_parent_type"}
    }
    assert not parent_checks, [str(e) for e in result.errors]


def test_exon_under_a_non_transcript_type_is_still_rejected(tmp_path):
    """The boundary: 'misc_feature' is outside the transcript namespace (it does not
    end in RNA or _transcript), so an exon beneath it is still a hierarchy error."""
    path = _write(tmp_path, [
        _row("gene", 100, 400, "g1"),
        _row("misc_feature", 100, 400, "p1", parent="g1"),
        _row("exon", 100, 199, "e1", parent="p1"),
    ])
    result = validate_gff3_structure(path)
    assert not result.is_valid
    assert any("exon_parent_type" in str(e) for e in result.errors), \
        [str(e) for e in result.errors]


def test_exon_under_a_genuinely_wrong_parent_is_still_rejected(tmp_path):
    """The check must keep its teeth: an exon parented to another EXON is invalid,
    because that parent is a leaf, not a middle level under a gene-like root."""
    path = _write(tmp_path, [
        _row("gene", 100, 400, "g1"),
        _row("mRNA", 100, 400, "t1", parent="g1"),
        _row("exon", 100, 199, "e1", parent="t1"),
        _row("exon", 150, 180, "e2", parent="e1"),      # exon under an exon
    ])
    result = validate_gff3_structure(path)
    assert not result.is_valid
    assert any("exon_parent_type" in str(e) or "exon_parent" in str(e)
               for e in result.errors), [str(e) for e in result.errors]
