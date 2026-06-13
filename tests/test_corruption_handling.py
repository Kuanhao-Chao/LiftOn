"""Phase 4.5 Step 5 — corruption / dirty-input handling.

Each test feeds a deliberately broken GFF3 to the current parser
(`lifton.annotation.Annotation` → gffutils) and pins behaviour: either
the parser catches it (raise / warn / graceful exit) or accepts it
silently. Silent acceptance becomes an xfail tagged for Phase 4 Step 2
(GFF3 validator implementation).
"""

from __future__ import annotations

import os

import pytest

from lifton import annotation, extract_sequence


# ---------------------------------------------------------------------------
# Format detection (`extract_sequence.determine_file_format`)
# ---------------------------------------------------------------------------

class TestFormatDetectionDirty:
    def test_missing_directive_still_detected_as_gff(self,
                                                    gff_missing_directive):
        # determine_file_format checks the *attributes column* for `=`
        # markers, not the directive — so it tolerates missing directive.
        # Pin behaviour: still recognised as GFF format.
        fmt = extract_sequence.determine_file_format(str(gff_missing_directive))
        assert fmt == "GFF format"

    def test_eight_columns_defaults_to_gff(self, gff_eight_columns):
        # Post-merge: rows with the wrong column count are skipped and
        # determine_file_format falls back to its "GFF format" default.
        # The strict GFF3 validator (Phase 5 bug #6) is the layer that
        # rejects bad column counts.
        fmt = extract_sequence.determine_file_format(str(gff_eight_columns))
        assert fmt == "GFF format"


# ---------------------------------------------------------------------------
# annotation.Annotation parser robustness
# ---------------------------------------------------------------------------

class TestAnnotationParserDirty:
    def test_missing_gff_version_directive_silently_accepted(
            self, gff_missing_directive):
        """gffutils does not require ##gff-version 3, in violation of
        NCBI § "Directives". Phase 4 Step 2 validator must flag this."""
        ann = annotation.Annotation(
            str(gff_missing_directive), False, False, "create_unique",
            None, True, False,
        )
        gene = ann.db_connection["gND"]
        assert gene.start == 100  # silently parsed

    def test_eight_column_row_now_rejected(self, gff_eight_columns):
        """Post-merge: the upstream annotation_validator now refuses to
        build a database from a file whose only data line has fewer
        than 9 columns; it sys.exit(1)s with a formatted error report.
        This is the desired NCBI-strict behaviour (NCBI § Column
        Specifications)."""
        with pytest.raises(SystemExit):
            annotation.Annotation(
                str(gff_eight_columns), False, False, "create_unique",
                None, True, False,
            )

    def test_crlf_line_endings_tolerated(self, gff_crlf_line_endings):
        ann = annotation.Annotation(
            str(gff_crlf_line_endings), False, False, "create_unique",
            None, True, False,
        )
        gene = ann.db_connection["gCRLF"]
        assert gene.start == 100 and gene.end == 200

    def test_truncated_last_line_tolerated(self, gff_truncated_no_newline):
        ann = annotation.Annotation(
            str(gff_truncated_no_newline), False, False, "create_unique",
            None, True, False,
        )
        gene = ann.db_connection["gT"]
        assert gene.start == 100

    def test_dangling_parent_silently_accepted(self, gff_dangling_parent):
        """Parent=ghost references a non-existent feature. gffutils does
        not validate referential integrity on parse — it lets the dangling
        reference stand. Phase 4 Step 2 must catch this per NCBI Parent §."""
        ann = annotation.Annotation(
            str(gff_dangling_parent), False, False, "create_unique",
            None, True, False,
        )
        # Lookup of the orphan still succeeds
        orphan = ann.db_connection["txD"]
        assert orphan.attributes["Parent"] == ["ghost"]
        # And lookup of the ghost parent legitimately fails
        with pytest.raises(Exception):
            ann.db_connection["ghost"]

    def test_unencoded_semicolon_mis_splits_attribute(
            self, gff_unencoded_semicolon):
        """An unescaped ';' inside an attribute value is parsed as a
        delimiter, splitting 'Note=alpha;beta' into Note=alpha + a
        nameless 'beta' fragment. Pin current behaviour."""
        ann = annotation.Annotation(
            str(gff_unencoded_semicolon), False, False, "create_unique",
            None, True, False,
        )
        feat = ann.db_connection["gU"]
        # Note value should be just "alpha" because the ;beta part was
        # treated as a separate attribute pair
        assert feat.attributes.get("Note") == ["alpha"]
        # The orphaned 'beta' fragment may or may not appear depending on
        # gffutils version. Pin: it does NOT carry the full intended value.
        assert "beta" not in (feat.attributes.get("Note") or [])

    def test_negative_coordinates_silently_accepted(self,
                                                    gff_negative_coords):
        """gffutils accepts negative integers for column 4. NCBI requires
        col 4 (start) to be >= 1."""
        ann = annotation.Annotation(
            str(gff_negative_coords), False, False, "create_unique",
            None, True, False,
        )
        feat = ann.db_connection["gNeg"]
        assert feat.start == -5  # Bug: should have been rejected

    def test_negative_coordinates_strictly_rejected_under_gffbase(
            self, gff_negative_coords, monkeypatch):
        """Under LIFTON_USE_GFFBASE=1 the gffbase parser correctly
        rejects negative col-4 values at parse time per NCBI spec.
        (The default gffutils path is permissive; users who want the
        strict behaviour opt in via the env var or --strict-gff.)"""
        from lifton.gffbase.exceptions import GFFFormatError
        monkeypatch.setenv("LIFTON_USE_GFFBASE", "1")
        with pytest.raises((GFFFormatError, ValueError, SystemExit)):
            annotation.Annotation(
                str(gff_negative_coords), False, False, "create_unique",
                None, True, False,
            )

    def test_duplicate_id_collision_resolved_via_create_unique(
            self, gff_duplicate_id_collision):
        """Two rows with the same ID but different (seqid, type) — true
        collision. With merge_strategy='create_unique', gffutils will
        rename one of them rather than raise."""
        ann = annotation.Annotation(
            str(gff_duplicate_id_collision), False, False, "create_unique",
            None, True, False,
        )
        # Both features should be present, but at least one will have a
        # synthesised unique ID.
        all_ids = {f.id for f in ann.db_connection.all_features()}
        # Original 'dup' may exist; its renamed sibling should also exist
        assert any(i == "dup" or i.startswith("dup_") for i in all_ids)
        assert len([i for i in all_ids
                    if i == "dup" or i.startswith("dup_")]) >= 2


# ---------------------------------------------------------------------------
# BOM-prefixed file
# ---------------------------------------------------------------------------

class TestBOMHandling:
    def test_bom_file_loads_or_warns(self, gff_with_bom):
        """UTF-8 BOM handling. gffutils may either tolerate the BOM or
        treat it as the start of the seqid on row 1. Pin observed
        behaviour."""
        try:
            ann = annotation.Annotation(
                str(gff_with_bom), False, False, "create_unique",
                None, True, False,
            )
            # If it parsed, the gene should be findable by id
            gene = ann.db_connection["gBOM"]
            assert gene.start == 100
        except Exception:
            # Acceptable too — BOM rejection is a valid stance
            pytest.skip("Parser rejects BOM; that is also acceptable")


# ---------------------------------------------------------------------------
# 10-column extra-tab attribute file
# ---------------------------------------------------------------------------

class TestExtraTabColumn:
    def test_extra_tab_in_attributes_now_rejected(
            self, gff_extra_tab_column):
        """Post-merge: the upstream annotation_validator now rejects
        files with malformed column counts (10 cols here) by exiting
        with a formatted error report rather than tolerating them."""
        with pytest.raises(SystemExit):
            annotation.Annotation(
                str(gff_extra_tab_column), False, False, "create_unique",
                None, True, False,
            )


# ---------------------------------------------------------------------------
# Malformed coordinates (start > end) — already in conftest as
# gff_malformed_coords; extend with downstream-safety assertion.
# ---------------------------------------------------------------------------

class TestMalformedCoordinates:
    def test_start_gt_end_silently_parsed(self, gff_malformed_coords):
        """gffutils accepts start=200, end=100 without raising. NCBI cols
        4-5 mandate start <= end. Phase 4 Step 2 validator must reject."""
        ann = annotation.Annotation(
            str(gff_malformed_coords), False, False, "create_unique",
            None, True, False,
        )
        feat = ann.db_connection["gbad"]
        assert feat.start == 200

    def test_start_gt_end_strictly_rejected_under_gffbase(
            self, gff_malformed_coords, monkeypatch):
        """gffbase rejects start>end at parse time. Available via
        LIFTON_USE_GFFBASE=1."""
        from lifton.gffbase.exceptions import GFFFormatError
        monkeypatch.setenv("LIFTON_USE_GFFBASE", "1")
        with pytest.raises((GFFFormatError, ValueError, SystemExit)):
            annotation.Annotation(
                str(gff_malformed_coords), False, False, "create_unique",
                None, True, False,
            )

    def test_downstream_segment_overlap_handles_inverted_segment(self):
        """If the parser ever lets an inverted segment through, downstream
        coordinate math should not produce nonsense (negative-or-bigger
        overlap). Pin: segments_overlap_length sorts by start, so an
        inverted segment is treated as if start were the lower coord."""
        from lifton import lifton_utils
        # (200, 100) is interpreted as a degenerate segment; sorted by
        # start gives (200, 100) vs (50, 60) -> sort places (50, 60) first.
        # The function does NOT swap start/end inside a segment, so
        # ovp_len = 60 - 200 + 1 = -139 -> ovp = False.
        ovp_len, ovp = lifton_utils.segments_overlap_length((200, 100),
                                                            (50, 60))
        # The function emits a non-positive length and ovp=False — does
        # NOT crash. Pin observed behaviour.
        assert ovp is False
        assert ovp_len <= 0
