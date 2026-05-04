"""Phase 5 — NCBI GFF3 strict validator (bug #6 fix).

Asserts that GFF3Validator catches every NCBI invariant the corruption
suite previously documented as silently accepted, and that strict-mode
pipeline invocation exits non-zero on errors.
"""

from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

import pytest

from lifton.io.gff3_validator import GFF3Validator, ValidationFinding
from lifton.io.ncbi_gff3_spec import (
    MULTI_VALUE_ATTRS,
    OFFICIAL_ATTRS,
    RESERVED_CHARS,
    VALID_PHASES,
    VALID_STRANDS,
)


# ---------------------------------------------------------------------------
# Spec constants — sanity checks
# ---------------------------------------------------------------------------

class TestSpecConstants:
    def test_official_attrs_includes_id_and_parent(self):
        assert "ID" in OFFICIAL_ATTRS
        assert "Parent" in OFFICIAL_ATTRS

    def test_multi_value_attrs_subset_of_official(self):
        assert MULTI_VALUE_ATTRS.issubset(OFFICIAL_ATTRS)

    def test_strands_match_gff3_spec(self):
        assert VALID_STRANDS == frozenset({"+", "-", ".", "?"})

    def test_phases_match_gff3_spec(self):
        assert VALID_PHASES == frozenset({"0", "1", "2", "."})

    def test_reserved_chars_includes_delimiters(self):
        for ch in (";", "=", ",", "&", "\t", "\n", "\r"):
            assert ch in RESERVED_CHARS


# ---------------------------------------------------------------------------
# Clean input — should produce zero findings
# ---------------------------------------------------------------------------

class TestCleanInput:
    def test_well_formed_gff_is_clean(self, gff_standard):
        findings = GFF3Validator().validate(gff_standard)
        # Allow warnings (e.g. unknown_seqid if no target supplied) but
        # zero errors.
        errors = [f for f in findings if f.severity == "error"]
        assert errors == [], f"unexpected errors: {errors}"

    def test_clean_with_target_seqids_yields_no_warning(self, gff_standard):
        findings = GFF3Validator(
            target_seqids={"chr1"},
        ).validate(gff_standard)
        assert findings == []


# ---------------------------------------------------------------------------
# Each NCBI invariant — one positive, one negative
# ---------------------------------------------------------------------------

class TestMissingDirective:
    def test_missing_gff_version_is_error(self, gff_missing_directive):
        findings = GFF3Validator().validate(gff_missing_directive)
        rules = {f.rule for f in findings if f.severity == "error"}
        assert "missing_gff_version" in rules


class TestColumnCount:
    def test_eight_column_row_is_error(self, gff_eight_columns):
        findings = GFF3Validator().validate(gff_eight_columns)
        rules = {f.rule for f in findings if f.severity == "error"}
        assert "bad_column_count" in rules

    def test_ten_column_row_is_error(self, gff_extra_tab_column):
        findings = GFF3Validator().validate(gff_extra_tab_column)
        rules = {f.rule for f in findings if f.severity == "error"}
        assert "bad_column_count" in rules


class TestCoordinates:
    def test_negative_start_is_error(self, gff_negative_coords):
        findings = GFF3Validator().validate(gff_negative_coords)
        rules = {f.rule for f in findings if f.severity == "error"}
        assert "negative_start" in rules

    def test_start_gt_end_is_error(self, gff_malformed_coords):
        findings = GFF3Validator().validate(gff_malformed_coords)
        rules = {f.rule for f in findings if f.severity == "error"}
        assert "start_gt_end" in rules

    def test_non_integer_coordinate_is_error(self, tmp_path):
        fp = tmp_path / "bad_int.gff3"
        fp.write_text(
            "##gff-version 3\n"
            "chr1\tt\tgene\tabc\t100\t.\t+\t.\tID=g\n"
        )
        findings = GFF3Validator().validate(fp)
        rules = {f.rule for f in findings if f.severity == "error"}
        assert "bad_coordinate" in rules


class TestStrand:
    def test_invalid_strand_is_error(self, tmp_path):
        fp = tmp_path / "bad_strand.gff3"
        fp.write_text(
            "##gff-version 3\n"
            "chr1\tt\tgene\t1\t10\t.\tX\t.\tID=g\n"
        )
        findings = GFF3Validator().validate(fp)
        rules = {f.rule for f in findings if f.severity == "error"}
        assert "bad_strand" in rules

    @pytest.mark.parametrize("strand", ["+", "-", ".", "?"])
    def test_valid_strand_is_clean(self, tmp_path, strand):
        fp = tmp_path / "ok.gff3"
        fp.write_text(
            "##gff-version 3\n"
            f"chr1\tt\tgene\t1\t10\t.\t{strand}\t.\tID=g\n"
        )
        errors = [f for f in GFF3Validator().validate(fp)
                  if f.severity == "error"]
        assert errors == []


class TestPhase:
    def test_invalid_phase_is_error(self, tmp_path):
        fp = tmp_path / "bad_phase.gff3"
        fp.write_text(
            "##gff-version 3\n"
            "chr1\tt\tCDS\t1\t9\t.\t+\t9\tID=c\n"
        )
        findings = GFF3Validator().validate(fp)
        rules = {f.rule for f in findings if f.severity == "error"}
        assert "bad_phase" in rules

    def test_cds_with_dot_phase_warns(self, tmp_path):
        fp = tmp_path / "cds_dot_phase.gff3"
        fp.write_text(
            "##gff-version 3\n"
            "chr1\tt\tCDS\t1\t9\t.\t+\t.\tID=c\n"
        )
        findings = GFF3Validator().validate(fp)
        rules = {f.rule for f in findings if f.severity == "warning"}
        assert "cds_missing_phase" in rules


class TestParentResolution:
    def test_dangling_parent_is_error(self, gff_dangling_parent):
        findings = GFF3Validator().validate(gff_dangling_parent)
        rules = {f.rule for f in findings if f.severity == "error"}
        assert "dangling_parent" in rules

    def test_resolved_parent_is_clean(self, gff_standard):
        findings = GFF3Validator().validate(gff_standard)
        rules = {f.rule for f in findings if f.severity == "error"}
        assert "dangling_parent" not in rules

    def test_multi_parent_resolves_both_ids(self, gff_shared_exon_isoforms):
        # Parent=tx1,tx2 should resolve cleanly because both tx1 and tx2
        # appear earlier in the file.
        findings = GFF3Validator().validate(gff_shared_exon_isoforms)
        rules = {f.rule for f in findings if f.severity == "error"}
        assert "dangling_parent" not in rules


class TestAttributeEncoding:
    def test_unencoded_semicolon_is_error(self, gff_unencoded_semicolon):
        findings = GFF3Validator().validate(gff_unencoded_semicolon)
        rules = {f.rule for f in findings if f.severity == "error"}
        # Note=alpha;beta is mis-split -> "beta" piece has no '='.
        assert "bad_attribute" in rules or "unencoded_reserved_char" in rules

    def test_percent_encoded_value_is_clean(self, tmp_path):
        # '%3B' is the percent-encoding of ';'
        fp = tmp_path / "encoded.gff3"
        fp.write_text(
            "##gff-version 3\n"
            "chr1\tt\tgene\t1\t10\t.\t+\t.\tID=g;Note=alpha%3Bbeta\n"
        )
        findings = GFF3Validator().validate(fp)
        rules = {f.rule for f in findings if f.severity == "error"}
        assert "unencoded_reserved_char" not in rules

    def test_attribute_without_equals_is_error(self, tmp_path):
        fp = tmp_path / "no_eq.gff3"
        fp.write_text(
            "##gff-version 3\n"
            "chr1\tt\tgene\t1\t10\t.\t+\t.\tID=g;loneFragment\n"
        )
        findings = GFF3Validator().validate(fp)
        rules = {f.rule for f in findings if f.severity == "error"}
        assert "bad_attribute" in rules


# ---------------------------------------------------------------------------
# BOM and CRLF
# ---------------------------------------------------------------------------

class TestEncodingQuirks:
    def test_bom_is_warning_not_error(self, gff_with_bom):
        findings = GFF3Validator().validate(gff_with_bom)
        warns = {f.rule for f in findings if f.severity == "warning"}
        assert "utf8_bom" in warns

    def test_crlf_line_endings_clean(self, gff_crlf_line_endings):
        # Both possible: zero errors. The validator strips \r before
        # tokenising and shouldn't flag anything.
        errors = [f for f in GFF3Validator().validate(gff_crlf_line_endings)
                  if f.severity == "error"]
        assert errors == []


# ---------------------------------------------------------------------------
# Unknown seqid warning
# ---------------------------------------------------------------------------

class TestUnknownSeqid:
    def test_unknown_seqid_warns_when_target_provided(self, gff_standard):
        findings = GFF3Validator(
            target_seqids={"chrZ"},
        ).validate(gff_standard)
        warns = {f.rule for f in findings if f.severity == "warning"}
        assert "unknown_seqid" in warns

    def test_no_unknown_seqid_warning_when_target_not_provided(self,
                                                               gff_standard):
        findings = GFF3Validator().validate(gff_standard)
        warns = {f.rule for f in findings if f.severity == "warning"}
        assert "unknown_seqid" not in warns


# ---------------------------------------------------------------------------
# has_errors / findings property / line-level entry point
# ---------------------------------------------------------------------------

class TestValidatorMisc:
    def test_has_errors_after_validate(self, gff_negative_coords):
        v = GFF3Validator()
        v.validate(gff_negative_coords)
        assert v.has_errors() is True

    def test_findings_property_returns_copy(self, gff_standard):
        v = GFF3Validator()
        v.validate(gff_standard)
        findings = v.findings
        findings.append(ValidationFinding("error", 1, "fake", "fake"))
        # Internal list is unchanged
        assert "fake" not in {f.rule for f in v.findings}

    def test_validate_line_returns_only_new_findings(self, gff_standard):
        v = GFF3Validator()
        v.validate(gff_standard)  # populate baseline
        baseline = len(v.findings)
        new = v.validate_line(99, "chr1\tt\tgene\t10\t1\t.\t+\t.\tID=ghost\n")
        assert any(f.rule == "start_gt_end" for f in new)
        # Baseline not duplicated
        assert len(v.findings) == baseline + len(new)

    def test_str_repr_of_finding(self):
        f = ValidationFinding("error", 7, "bad_strand", "bad strand value")
        s = str(f)
        assert "GFF3:error" in s
        assert "line 7" in s
        assert "bad_strand" in s


# ---------------------------------------------------------------------------
# CLI strict-mode integration
# ---------------------------------------------------------------------------

class TestStrictModeCLI:
    def test_strict_mode_exits_two_on_invalid_gff(self, tmp_path,
                                                  fasta_standard,
                                                  monkeypatch):
        """End-to-end: invoke run_all_lifton_steps with --strict-gff
        against a malformed GFF3 and assert SystemExit(2)."""
        from lifton import lifton as lifton_main, run_liftoff, run_miniprot
        from lifton import lifton_utils

        # Hermetic patches — never spawn liftoff/miniprot
        monkeypatch.setattr(lifton_utils, "check_miniprot_installed",
                            lambda: None)
        monkeypatch.setattr(run_miniprot, "check_miniprot_installed",
                            lambda: True)
        def _fail(*a, **k):
            raise AssertionError("Should not invoke external runners")
        monkeypatch.setattr(run_liftoff, "run_liftoff", _fail)
        monkeypatch.setattr(run_miniprot, "run_miniprot", _fail)

        bad_gff = tmp_path / "bad.gff3"
        bad_gff.write_text(
            # Missing ##gff-version directive AND start>end
            "chr1\tt\tgene\t200\t100\t.\t+\t.\tID=gbad\n"
        )
        out = tmp_path / "out.gff3"
        argv = [
            str(fasta_standard), str(fasta_standard),
            "-g", str(bad_gff),
            "-L", str(bad_gff),  # short-circuit liftoff exec
            "-M", str(bad_gff),  # short-circuit miniprot exec
            "-o", str(out),
            "-ad", "RefSeq",
            "--strict-gff",
            "--force",
        ]
        args = lifton_main.parse_args(argv)
        with pytest.raises(SystemExit) as exc:
            lifton_main.run_all_lifton_steps(args)
        assert exc.value.code == 2
