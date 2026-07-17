"""
annotation_validator.py  —  Pre-flight validation of GFF3 / GTF annotation files.

Called before any gffutils database build to give the user an actionable error
message rather than a cryptic SQLite UNIQUE-constraint failure.
"""

from __future__ import annotations

import hashlib
import os
import sys
import re
from dataclasses import dataclass, field
from collections import defaultdict
from typing import List, Dict, Optional

from lifton.io.gff3_validator import GFF3Validator, ValidationFinding


# ─────────────────────────────────────────────────────────────────────────────
# Data classes
# ─────────────────────────────────────────────────────────────────────────────

@dataclass
class ValidationResult:
    """Result of validating a single annotation file."""
    file_path: str
    is_valid: bool = True            # False → fatal, abort DB build
    file_format: str = "unknown"     # "GFF3" | "GTF" | "unknown"
    total_lines: int = 0
    data_lines: int = 0
    comment_lines: int = 0
    duplicate_ids: Dict[str, List[int]] = field(default_factory=dict)
    orphan_parents: List[str] = field(default_factory=list)
    warnings: List[str] = field(default_factory=list)
    errors: List[str] = field(default_factory=list)
    fix_suggestions: List[str] = field(default_factory=list)


@dataclass(frozen=True)
class AnnotationScanResult:
    """Immutable metadata and preflight findings from one stable file pass."""

    file_path: str
    sha256: str = ""
    size: int = 0
    source_device: int = 0
    source_inode: int = 0
    source_mtime_ns: int = 0
    file_format: str = "unknown"
    total_lines: int = 0
    data_lines: int = 0
    comment_lines: int = 0
    directives: tuple[str, ...] = ()
    dialect_votes: tuple[tuple[str, int], ...] = ()
    feature_type_counts: tuple[tuple[str, int], ...] = ()
    seqids: tuple[str, ...] = ()
    missing_target_seqids: tuple[str, ...] = ()
    duplicate_ids: tuple[tuple[str, tuple[int, ...]], ...] = ()
    discontinuous_cds_ids: tuple[str, ...] = ()
    orphan_parents: tuple[str, ...] = ()
    bad_column_count_lines: tuple[int, ...] = ()
    ncbi_findings: tuple[ValidationFinding, ...] = ()
    ncbi_finding_counts: tuple[tuple[str, int], ...] = ()
    warnings: tuple[str, ...] = ()
    errors: tuple[str, ...] = ()

    @property
    def source_fingerprint(self) -> dict[str, object]:
        return {"sha256": self.sha256, "size": self.size}

    @property
    def duplicate_id_map(self) -> Dict[str, List[int]]:
        return {
            feature_id: list(lines)
            for feature_id, lines in self.duplicate_ids
        }


# ─────────────────────────────────────────────────────────────────────────────
# Public API
# ─────────────────────────────────────────────────────────────────────────────

def _attribute_values(attributes: str) -> Dict[str, tuple[str, ...]]:
    """Parse only hierarchy identity fields needed by the consolidated scan."""

    parsed: Dict[str, tuple[str, ...]] = {}
    for piece in attributes.split(";"):
        piece = piece.strip()
        if not piece:
            continue
        if piece.startswith("ID="):
            key, value = "ID", piece[3:]
        elif piece.startswith("Parent="):
            key, value = "Parent", piece[7:]
        else:
            key, separator, value = piece.partition("=")
            if not separator:
                continue
            key = key.strip()
            if key != "ID" and key != "Parent":
                continue
        values = []
        for item in value.split(","):
            item = item.strip()
            if item:
                values.append(item)
        parsed[key] = tuple(values)
    return parsed


def _cds_signature(columns: List[str], attrs: Dict[str, tuple[str, ...]]) -> tuple | None:
    parents = tuple(sorted(set(attrs.get("Parent", ()))))
    if columns[2] != "CDS" or not parents:
        return None
    return columns[2], columns[0], columns[1], columns[6], parents


def scan_annotation(
    file_path: str,
    target_seqids=None,
    max_examples: int = 20,
    max_findings_per_rule: int = 100,
) -> AnnotationScanResult:
    """Scan annotation bytes once and return stable, reusable metadata."""

    path = os.fspath(file_path)
    errors: List[str] = []
    warnings: List[str] = []
    try:
        before = os.stat(path)
    except OSError as exc:
        message = f"File not found: {path}" if not os.path.exists(path) else f"Cannot read file: {exc}"
        return AnnotationScanResult(file_path=path, errors=(message,))
    if not os.path.isfile(path):
        return AnnotationScanResult(
            file_path=path,
            errors=(f"Path is not a regular file: {path}",),
        )
    if before.st_size == 0:
        return AnnotationScanResult(
            file_path=path,
            errors=(f"File is empty: {path}",),
        )

    digest = hashlib.sha256()
    bytes_read = 0
    target_seqid_set = (
        set(target_seqids) if target_seqids is not None else None
    )
    total_lines = data_lines = comment_lines = 0
    directives: List[str] = []
    dialect_counts = {"GFF3": 0, "GTF": 0}
    feature_counts: Dict[str, int] = defaultdict(int)
    seqids: set[str] = set()
    missing_seqids: set[str] = set()
    bad_columns: List[int] = []
    all_ids: set[str] = set()
    parent_first_line: Dict[str, int] = {}
    id_state: Dict[str, dict] = {}
    ncbi_findings: List[ValidationFinding] = []
    ncbi_counts: Dict[str, int] = defaultdict(int)
    ncbi_validator = GFF3Validator(
        target_seqids=target_seqid_set,
        track_references=False,
    )
    gtf_key = re.compile(r'(?:^|;)\s*(?:gene_id|transcript_id)\s+"')

    def keep_finding(finding: ValidationFinding) -> None:
        ncbi_counts[finding.rule] += 1
        if ncbi_counts[finding.rule] <= max_findings_per_rule:
            ncbi_findings.append(finding)

    try:
        with open(path, "rb") as handle:
            for lineno, raw_line in enumerate(handle, start=1):
                digest.update(raw_line)
                bytes_read += len(raw_line)
                total_lines += 1
                line = raw_line.decode("utf-8", errors="replace").rstrip("\r\n")

                # Drain per-line findings to keep malformed inputs bounded.
                ncbi_validator._validate_line(lineno, line)
                for finding in ncbi_validator._findings:
                    keep_finding(finding)
                ncbi_validator._findings.clear()
                ncbi_validator._declared_ids.clear()
                ncbi_validator._referenced_parents.clear()

                if not line.strip() or line.startswith("#"):
                    comment_lines += 1
                    if line.startswith("##") or line.startswith("#!"):
                        directives.append(line)
                    continue

                columns = line.split("\t")
                if len(columns) != 9:
                    if len(bad_columns) < max_examples:
                        bad_columns.append(lineno)
                    continue

                data_lines += 1
                seqid, _source, feature_type, _start, _end, _score, _strand, _phase, attributes = columns
                feature_counts[feature_type] += 1
                seqids.add(seqid)
                if target_seqid_set is not None and seqid not in target_seqid_set:
                    missing_seqids.add(seqid)

                if any(key in attributes for key in ("ID=", "Parent=", "Name=")):
                    dialect_counts["GFF3"] += 1
                if gtf_key.search(attributes):
                    dialect_counts["GTF"] += 1

                attrs = _attribute_values(attributes)
                feature_values = attrs.get("ID", ())
                feature_id = feature_values[0] if feature_values else ""
                if feature_id:
                    all_ids.add(feature_id)
                    state = id_state.get(feature_id)
                    signature = _cds_signature(columns, attrs)
                    if state is None:
                        id_state[feature_id] = {
                            "first_line": lineno,
                            "lines": [lineno],
                            "signature": signature,
                            "count": 1,
                            "valid_cds": signature is not None,
                        }
                    else:
                        state["count"] += 1
                        if len(state["lines"]) < max_examples:
                            state["lines"].append(lineno)
                        if signature is None or signature != state["signature"]:
                            state["valid_cds"] = False

                for parent_id in attrs.get("Parent", ()):
                    parent_first_line.setdefault(parent_id, lineno)
    except OSError as exc:
        errors.append(f"Cannot read file: {exc}")

    try:
        after = os.stat(path)
    except OSError as exc:
        errors.append(f"Cannot stat file after scan: {exc}")
        after = before
    before_key = (before.st_dev, before.st_ino, before.st_size, before.st_mtime_ns)
    after_key = (after.st_dev, after.st_ino, after.st_size, after.st_mtime_ns)
    if before_key != after_key or bytes_read != after.st_size:
        errors.append(f"Annotation changed while it was being scanned: {path}")

    if dialect_counts["GFF3"] and dialect_counts["GFF3"] >= dialect_counts["GTF"]:
        file_format = "GFF3"
    elif dialect_counts["GTF"]:
        file_format = "GTF"
        warnings.append(
            "File appears to be in GTF format. LiftOn will attempt auto-conversion to GFF3."
        )
    else:
        file_format = "unknown"
        warnings.append(
            "Could not confidently determine file format (neither GFF3 nor GTF patterns detected). "
            "Will attempt to load as GFF3."
        )

    if data_lines == 0:
        errors.append(
            "File contains no valid 9-column data lines. "
            f"Total lines: {total_lines}, comment/blank lines: {comment_lines}."
        )
    if bad_columns:
        warnings.append(
            f"Lines with wrong number of columns (expected 9): {bad_columns}"
        )

    duplicate_items = []
    valid_discontinuous = []
    for feature_id, state in sorted(id_state.items()):
        if state["count"] <= 1:
            continue
        if state["valid_cds"]:
            valid_discontinuous.append(feature_id)
        else:
            duplicate_items.append((feature_id, tuple(state["lines"])))

    if duplicate_items:
        examples = "; ".join(
            f"ID={feature_id!r} (lines {list(lines[:3])}{'…' if len(lines) > 3 else ''})"
            for feature_id, lines in duplicate_items[:max_examples]
        )
        warnings.append(
            f"{len(duplicate_items)} colliding feature ID(s) detected. Examples: {examples}"
        )

    orphan_items = sorted(
        parent_id for parent_id in parent_first_line
        if parent_id and parent_id not in all_ids
    )
    for parent_id in orphan_items:
        keep_finding(ValidationFinding(
            "error",
            parent_first_line[parent_id],
            "dangling_parent",
            f"Parent={parent_id!r} references an ID that never appears in the file.",
        ))
    if orphan_items:
        warnings.append(
            f"{len(orphan_items)} Parent reference(s) point to unknown IDs. "
            f"Examples: {orphan_items[:5]}"
        )

    return AnnotationScanResult(
        file_path=path,
        sha256=digest.hexdigest(),
        size=after.st_size,
        source_device=after.st_dev,
        source_inode=after.st_ino,
        source_mtime_ns=after.st_mtime_ns,
        file_format=file_format,
        total_lines=total_lines,
        data_lines=data_lines,
        comment_lines=comment_lines,
        directives=tuple(directives),
        dialect_votes=tuple(sorted(dialect_counts.items())),
        feature_type_counts=tuple(sorted(feature_counts.items())),
        seqids=tuple(sorted(seqids)),
        missing_target_seqids=tuple(sorted(missing_seqids)),
        duplicate_ids=tuple(duplicate_items),
        discontinuous_cds_ids=tuple(valid_discontinuous),
        orphan_parents=tuple(orphan_items[:max_findings_per_rule]),
        bad_column_count_lines=tuple(bad_columns),
        ncbi_findings=tuple(ncbi_findings),
        ncbi_finding_counts=tuple(sorted(ncbi_counts.items())),
        warnings=tuple(warnings),
        errors=tuple(errors),
    )


def _validation_from_scan(
    scan: AnnotationScanResult,
    *,
    max_duplicate_examples: int,
    check_orphan_parents: bool,
) -> ValidationResult:
    duplicate_items = list(scan.duplicate_ids)[:max_duplicate_examples]
    result = ValidationResult(
        file_path=scan.file_path,
        is_valid=not scan.errors and not scan.duplicate_ids,
        file_format=scan.file_format,
        total_lines=scan.total_lines,
        data_lines=scan.data_lines,
        comment_lines=scan.comment_lines,
        duplicate_ids={key: list(lines) for key, lines in duplicate_items},
        orphan_parents=list(scan.orphan_parents) if check_orphan_parents else [],
        warnings=list(scan.warnings),
        errors=list(scan.errors),
    )
    if scan.file_format == "GTF":
        result.fix_suggestions.append(
            "For best results, convert GTF to GFF3 manually before running LiftOn:\n"
            f"  gffread -E {scan.file_path} -o {scan.file_path}.gff3\n"
            f"  or: agat_sp_gtf2gff.pl --gtf {scan.file_path} -o {scan.file_path}.gff3"
        )
    if scan.duplicate_ids:
        result.fix_suggestions.append(
            "Repeated IDs are accepted only for CDS segments with the same seqid, source, "
            "strand, and non-empty Parent set. Fix the reported collisions before reuse."
        )
    if check_orphan_parents and scan.orphan_parents:
        result.fix_suggestions.append(
            "Ensure every Parent value names an ID declared elsewhere in the annotation."
        )
    if scan.errors:
        result.fix_suggestions.append(
            "Check that the file is a stable, readable GFF3 or GTF annotation with nine columns."
        )
    return result


def validate_annotation_file(
    file_path: str,
    max_duplicate_examples: int = 20,
    check_orphan_parents: bool = True,
    *,
    scan_result: Optional[AnnotationScanResult] = None,
) -> ValidationResult:
    """
    Validate *file_path* before passing it to gffutils.

    Parameters
    ----------
    file_path : str
        Path to GFF3 or GTF annotation file.
    max_duplicate_examples : int
        Maximum number of duplicate-ID examples to include in the report.
    check_orphan_parents : bool
        Whether to verify that every Parent= reference resolves to a known ID.

    Returns
    -------
    ValidationResult
        Populated result. ``result.is_valid`` is False if the file cannot be
        used for a DB build without intervention.
    """
    scan = scan_result or scan_annotation(
        file_path,
        max_examples=max_duplicate_examples,
    )
    return _validation_from_scan(
        scan,
        max_duplicate_examples=max_duplicate_examples,
        check_orphan_parents=check_orphan_parents,
    )


# ─────────────────────────────────────────────────────────────────────────────
# Reporting helpers
# ─────────────────────────────────────────────────────────────────────────────

def print_validation_report(result: ValidationResult, always_show: bool = False) -> None:
    """
    Print a human-readable validation report to stderr.

    Parameters
    ----------
    result : ValidationResult
    always_show : bool
        If True, print even when the file is valid.
    """
    if result.is_valid and not always_show and not result.warnings:
        return

    width = 70

    def _box_line(text="", fill="─"):
        if text:
            pad = width - len(text) - 4
            right_pad = max(0, pad)
            return f"│ {text}{' ' * right_pad} │"
        return f"├{'─' * (width - 2)}┤"

    lines_out = []

    if result.errors:
        header = "╔" + "═" * (width - 2) + "╗"
        title = f"║  ❌  ANNOTATION FILE ERROR" + " " * (width - 28) + "║"
        footer = "╚" + "═" * (width - 2) + "╝"
    else:
        header = "╔" + "═" * (width - 2) + "╗"
        title = f"║  ⚠️   ANNOTATION FILE WARNING" + " " * (width - 31) + "║"
        footer = "╚" + "═" * (width - 2) + "╝"

    lines_out.append(header)
    lines_out.append(title)
    lines_out.append("╠" + "═" * (width - 2) + "╣")
    lines_out.append(_box_line(f"File   : {result.file_path}"))
    lines_out.append(_box_line(f"Format : {result.file_format}"))
    lines_out.append(_box_line(
        f"Lines  : {result.total_lines} total  "
        f"({result.data_lines} data, {result.comment_lines} comment/blank)"
    ))

    if result.errors:
        lines_out.append(_box_line())
        lines_out.append(_box_line("ERRORS:"))
        for err in result.errors:
            for chunk in _wrap(err, width - 6):
                lines_out.append(_box_line(f"  • {chunk}"))

    if result.warnings:
        lines_out.append(_box_line())
        lines_out.append(_box_line("WARNINGS:"))
        for warn in result.warnings:
            for chunk in _wrap(warn, width - 6):
                lines_out.append(_box_line(f"  ⚠ {chunk}"))

    if result.fix_suggestions:
        lines_out.append(_box_line())
        lines_out.append(_box_line("SUGGESTED FIXES:"))
        for sug in result.fix_suggestions:
            for chunk in _wrap(sug, width - 6):
                lines_out.append(_box_line(f"  → {chunk}"))

    lines_out.append(footer)

    for line in lines_out:
        # Trim to width to avoid misaligned boxes when text is very long
        print(line[:width + 4], file=sys.stderr)


def _wrap(text: str, max_width: int) -> List[str]:
    """Wrap *text* to *max_width* characters per line, preserving newlines."""
    result_lines: List[str] = []
    for paragraph in text.split("\n"):
        if len(paragraph) <= max_width:
            result_lines.append(paragraph)
        else:
            words = paragraph.split(" ")
            current = ""
            for word in words:
                if len(current) + len(word) + 1 <= max_width:
                    current = (current + " " + word).lstrip()
                else:
                    if current:
                        result_lines.append(current)
                    current = word
            if current:
                result_lines.append(current)
    return result_lines or [""]


def print_db_build_error(file_path: str, strategy: str, exception: Exception) -> None:
    """Print a boxed error when a gffutils DB build strategy fails."""
    width = 70
    exc_str = str(exception)
    print("╔" + "═" * (width - 2) + "╗", file=sys.stderr)
    print(f"║  ❌  gffutils DB BUILD FAILED (strategy={strategy!r})" +
          " " * max(0, width - 50) + "║", file=sys.stderr)
    print("╠" + "═" * (width - 2) + "╣", file=sys.stderr)
    print(f"│  File      : {file_path}"[:width] + " │", file=sys.stderr)
    exc_lines = _wrap(f"Error     : {exc_str}", width - 6)
    for ln in exc_lines:
        print(f"│  {ln:<{width - 6}} │", file=sys.stderr)
    print("╚" + "═" * (width - 2) + "╝", file=sys.stderr)


def print_db_build_success(file_path: str, strategy: str) -> None:
    """Print a one-liner when a gffutils DB build succeeds."""
    short = os.path.basename(file_path)
    print(
        f"  ✓ gffutils database built successfully "
        f"[file={short!r}, strategy={strategy!r}]",
        file=sys.stderr,
    )
