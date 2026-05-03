"""Streaming, single-pass NCBI GFF3 validator (Phase 5 bug #6 fix).

Reads a GFF3 file line-by-line and emits structured findings for
violations of the NCBI specification:
https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/file-formats/annotation-files/about-ncbi-gff3/

Severity model:
  - "error":   the input violates a hard NCBI invariant. Strict mode
               (`strict=True`) treats any error as fatal.
  - "warning": the input is technically allowed but suspicious or
               likely to mislead downstream tools.

The validator does NOT mutate the input file. Findings are returned
as a list; callers decide whether to log, fail, or filter.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

from lifton.io.ncbi_gff3_spec import (
    DIRECTIVE_PREFIX,
    GFF_VERSION_DIRECTIVE,
    MULTI_VALUE_ATTRS,
    NCBI_DIRECTIVE_PREFIX,
    OFFICIAL_ATTRS,
    RESERVED_CHARS,
    VALID_PHASES,
    VALID_STRANDS,
)


@dataclass(frozen=True)
class ValidationFinding:
    severity: str           # "error" or "warning"
    line_no: int            # 1-based; 0 means "file-level"
    rule: str               # short stable identifier
    message: str            # human-readable explanation

    def __str__(self) -> str:
        loc = f"line {self.line_no}" if self.line_no else "file"
        return f"[GFF3:{self.severity}] {loc}: {self.rule} — {self.message}"


# A pre-compiled regex for percent-encoded sequences in attribute values.
_PCT_RE = re.compile(r"%[0-9A-Fa-f]{2}")


def _has_unencoded_reserved(value: str) -> bool:
    """Return True iff the attribute value contains a reserved character
    that is NOT percent-encoded.

    We strip valid %XX sequences first, then scan for any reserved char
    in the remainder.
    """
    stripped = _PCT_RE.sub("", value)
    return any(ch in RESERVED_CHARS for ch in stripped)


class GFF3Validator:
    """Validate a GFF3 file against the NCBI specification.

    Parameters
    ----------
    target_seqids:
        Optional set of seqids known to be present in the target/reference
        FASTA. Rows whose seqid is missing produce a warning, not an
        error (NCBI § Col 1 recommends but does not require accession.version).
    strict:
        Affects only the caller — the validator always returns the full
        finding list. Callers use `has_errors()` + `strict` to decide.
    """

    def __init__(self, *, target_seqids: set[str] | None = None,
                 strict: bool = False) -> None:
        self.target_seqids = target_seqids
        self.strict = strict
        self._findings: list[ValidationFinding] = []
        self._declared_ids: set[str] = set()
        self._referenced_parents: list[tuple[int, str]] = []

    # ---------------------- public API -----------------------------

    def validate(self, path: str | Path) -> list[ValidationFinding]:
        self._findings = []
        self._declared_ids = set()
        self._referenced_parents = []
        path = Path(path)
        with open(path, "r", encoding="utf-8", newline="") as fh:
            for line_no, raw in enumerate(fh, start=1):
                self._validate_line(line_no, raw)
        # After the file is consumed, check that every Parent reference
        # resolves to a declared ID.
        for line_no, parent_id in self._referenced_parents:
            if parent_id not in self._declared_ids:
                self._add("error", line_no, "dangling_parent",
                          f"Parent={parent_id!r} references an ID that "
                          f"never appears in the file (NCBI § Parent).")
        return list(self._findings)

    def validate_line(self, line_no: int, raw: str) -> list[ValidationFinding]:
        """Public, line-level entry point. Useful for the post-write
        sanity check in lifton.io.gff_writer (Phase 5 Step 8)."""
        # Reset only the per-call buffer; do NOT touch declared_ids.
        before = len(self._findings)
        self._validate_line(line_no, raw)
        return list(self._findings[before:])

    def has_errors(self) -> bool:
        return any(f.severity == "error" for f in self._findings)

    @property
    def findings(self) -> list[ValidationFinding]:
        return list(self._findings)

    # ---------------------- internals ------------------------------

    def _add(self, severity: str, line_no: int, rule: str, message: str
             ) -> None:
        self._findings.append(
            ValidationFinding(severity, line_no, rule, message)
        )

    def _validate_line(self, line_no: int, raw: str) -> None:
        # Strip a UTF-8 BOM that may sneak in on the first line.
        if line_no == 1 and raw.startswith("﻿"):
            self._add("warning", 1, "utf8_bom",
                      "File begins with a UTF-8 BOM; some tools "
                      "mis-parse the first column.")
            raw = raw.lstrip("﻿")

        line = raw.rstrip("\r\n")
        if not line:
            return

        # Directives ----------------------------------------------------
        if line.startswith(DIRECTIVE_PREFIX) or \
                line.startswith(NCBI_DIRECTIVE_PREFIX):
            if line_no == 1 and not line.startswith(GFF_VERSION_DIRECTIVE):
                self._add("error", 1, "missing_gff_version",
                          "First line must be '##gff-version 3' "
                          "(NCBI § Directives).")
            return

        # Comment lines (single '#') are tolerated.
        if line.startswith("#"):
            return

        # If line 1 is a feature row, the gff-version directive is missing.
        if line_no == 1:
            self._add("error", 1, "missing_gff_version",
                      "First line must be '##gff-version 3' before any "
                      "feature rows (NCBI § Directives).")

        # Feature row ---------------------------------------------------
        cols = line.split("\t")
        if len(cols) != 9:
            self._add("error", line_no, "bad_column_count",
                      f"Expected 9 tab-separated columns, got {len(cols)} "
                      "(NCBI § Column Specifications).")
            return

        seqid, _src, _type, start_s, end_s, _score, strand, phase, attrs = cols

        # Col 1 ---------------------------------------------------------
        if not seqid:
            self._add("error", line_no, "empty_seqid",
                      "Column 1 (seqid) is empty.")
        elif self.target_seqids is not None and seqid not in self.target_seqids:
            self._add("warning", line_no, "unknown_seqid",
                      f"seqid {seqid!r} is not in the target FASTA index.")

        # Cols 4-5 ------------------------------------------------------
        try:
            start = int(start_s)
            end = int(end_s)
        except ValueError:
            self._add("error", line_no, "bad_coordinate",
                      f"Cols 4-5 must be integers (got {start_s!r}, "
                      f"{end_s!r}); NCBI § Cols 4-5.")
            return

        if start < 1:
            self._add("error", line_no, "negative_start",
                      f"Col 4 (start) must be >= 1 (got {start}); "
                      "NCBI § Cols 4-5 (1-based).")
        if end < 1:
            self._add("error", line_no, "negative_end",
                      f"Col 5 (end) must be >= 1 (got {end}); "
                      "NCBI § Cols 4-5 (1-based).")
        if start > end:
            self._add("error", line_no, "start_gt_end",
                      f"Col 4 (start={start}) must be <= col 5 (end={end}); "
                      "NCBI § Cols 4-5.")

        # Col 7 ---------------------------------------------------------
        if strand not in VALID_STRANDS:
            self._add("error", line_no, "bad_strand",
                      f"Col 7 (strand) must be in {sorted(VALID_STRANDS)}; "
                      f"got {strand!r}.")

        # Col 8 ---------------------------------------------------------
        if phase not in VALID_PHASES:
            self._add("error", line_no, "bad_phase",
                      f"Col 8 (phase) must be in {sorted(VALID_PHASES)}; "
                      f"got {phase!r}.")
        elif _type == "CDS" and phase == ".":
            # NCBI explicitly notes phase may be wrong on pseudogenes;
            # warn rather than reject.
            self._add("warning", line_no, "cds_missing_phase",
                      "CDS row has no phase (col 8 = '.'); NCBI § Col 8 "
                      "permits this only for pseudogenes / internal "
                      "frameshifts.")

        # Col 9 ---------------------------------------------------------
        attr_dict: dict[str, str] = {}
        for piece in attrs.split(";"):
            piece = piece.strip()
            if not piece:
                continue
            if "=" not in piece:
                self._add("error", line_no, "bad_attribute",
                          f"Attribute {piece!r} is not 'key=value'.")
                continue
            key, _, value = piece.partition("=")
            attr_dict[key] = value
            # Capitalisation: official attrs start with an uppercase
            # letter; misspellings (e.g. 'parent' instead of 'Parent')
            # are case-sensitive errors per the spec.
            if key not in OFFICIAL_ATTRS and key[:1].isupper() and \
                    key.lower() in {a.lower() for a in OFFICIAL_ATTRS}:
                self._add("warning", line_no, "miscapitalised_attribute",
                          f"Attribute {key!r} likely a miscapitalised "
                          f"official attribute.")
            if _has_unencoded_reserved(value):
                self._add("error", line_no, "unencoded_reserved_char",
                          f"Attribute {key!r} value contains an "
                          "unencoded reserved character (NCBI § "
                          "Attribute Specifications). Reserved chars "
                          "must be percent-encoded.")

        # Track ID for parent-resolution check.
        if "ID" in attr_dict:
            self._declared_ids.add(attr_dict["ID"])
        if "Parent" in attr_dict:
            for parent_id in attr_dict["Parent"].split(","):
                parent_id = parent_id.strip()
                if parent_id:
                    self._referenced_parents.append((line_no, parent_id))


def validate_path(path: str | Path, *,
                  target_seqids: set[str] | None = None,
                  strict: bool = False) -> list[ValidationFinding]:
    """Convenience function used by run_all_lifton_steps."""
    return GFF3Validator(target_seqids=target_seqids,
                         strict=strict).validate(path)
