"""
annotation_validator.py  —  Pre-flight validation of GFF3 / GTF annotation files.

Called before any gffutils database build to give the user an actionable error
message rather than a cryptic SQLite UNIQUE-constraint failure.
"""

import os
import sys
import re
from dataclasses import dataclass, field
from collections import defaultdict
from typing import List, Dict, Optional


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


# ─────────────────────────────────────────────────────────────────────────────
# Public API
# ─────────────────────────────────────────────────────────────────────────────

def validate_annotation_file(
    file_path: str,
    max_duplicate_examples: int = 20,
    check_orphan_parents: bool = True,
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
    result = ValidationResult(file_path=file_path)

    # ── 1. File existence & readability ──────────────────────────────────────
    if not os.path.exists(file_path):
        result.is_valid = False
        result.errors.append(f"File not found: {file_path}")
        result.fix_suggestions.append(
            f"Check that the path is correct: {file_path}"
        )
        return result

    if not os.path.isfile(file_path):
        result.is_valid = False
        result.errors.append(f"Path is not a regular file: {file_path}")
        return result

    if os.path.getsize(file_path) == 0:
        result.is_valid = False
        result.errors.append(f"File is empty: {file_path}")
        return result

    if not os.access(file_path, os.R_OK):
        result.is_valid = False
        result.errors.append(f"File is not readable (check permissions): {file_path}")
        return result

    # ── 2. Parse the file (single streaming pass) ────────────────────────────
    id_linenos: Dict[str, List[int]] = defaultdict(list)   # id → [line numbers]
    all_ids: set = set()
    all_parent_refs: set = set()
    gff_count = 0
    gtf_count = 0
    data_lines = 0
    total_lines = 0
    comment_lines = 0
    bad_column_count_lines: List[int] = []

    _GFF_KEYS = {"ID=", "Parent=", "Name="}
    _GTF_KEYS_RE = re.compile(r'(gene_id|transcript_id)\s+"')
    _ID_RE = re.compile(r'(?:^|;)ID=([^;]+)')
    _PARENT_RE = re.compile(r'(?:^|;)Parent=([^;]+)')

    try:
        with open(file_path, "r", errors="replace") as fh:
            for lineno, raw_line in enumerate(fh, start=1):
                total_lines += 1
                line = raw_line.rstrip("\n\r")

                if line.startswith("#") or line.strip() == "":
                    comment_lines += 1
                    continue

                cols = line.split("\t")
                if len(cols) != 9:
                    if len(bad_column_count_lines) < 5:
                        bad_column_count_lines.append(lineno)
                    # Don't count as a data line
                    continue

                data_lines += 1
                attrs = cols[8]

                # Format detection
                if any(k in attrs for k in _GFF_KEYS) and "=" in attrs:
                    gff_count += 1
                elif _GTF_KEYS_RE.search(attrs):
                    gtf_count += 1

                # Extract IDs
                m_id = _ID_RE.search(attrs)
                if m_id:
                    feat_id = m_id.group(1).strip()
                    id_linenos[feat_id].append(lineno)
                    all_ids.add(feat_id)

                # Extract Parent references
                m_par = _PARENT_RE.search(attrs)
                if m_par:
                    # Parent can be a comma-separated list
                    for p in m_par.group(1).split(","):
                        all_parent_refs.add(p.strip())

    except OSError as exc:
        result.is_valid = False
        result.errors.append(f"Cannot read file: {exc}")
        return result

    result.total_lines = total_lines
    result.data_lines = data_lines
    result.comment_lines = comment_lines

    # ── 3. Format determination ───────────────────────────────────────────────
    if gff_count > 0 and gff_count >= gtf_count:
        result.file_format = "GFF3"
    elif gtf_count > 0:
        result.file_format = "GTF"
        result.warnings.append(
            "File appears to be in GTF format. "
            "LiftOn will attempt auto-conversion to GFF3."
        )
        result.fix_suggestions.append(
            "For best results, convert GTF to GFF3 manually before running LiftOn:\n"
            "  gffread -E {f} -o {f}.gff3\n"
            "  or: agat_sp_gtf2gff.pl --gtf {f} -o {f}.gff3".format(f=file_path)
        )
    else:
        result.file_format = "unknown"
        result.warnings.append(
            "Could not confidently determine file format (neither GFF3 nor GTF "
            "patterns detected). Will attempt to load as GFF3."
        )

    # ── 4. Structural sanity ─────────────────────────────────────────────────
    if data_lines == 0:
        result.is_valid = False
        result.errors.append(
            "File contains no valid 9-column data lines. "
            f"Total lines: {total_lines}, comment/blank lines: {comment_lines}."
        )
        result.fix_suggestions.append(
            "Check that the file is a valid GFF3 or GTF annotation file."
        )
        return result

    if bad_column_count_lines:
        result.warnings.append(
            f"Lines with wrong number of columns (expected 9): "
            f"{bad_column_count_lines[:5]}"
            + (" (showing first 5)" if len(bad_column_count_lines) > 5 else "")
        )

    # ── 5. Duplicate ID detection ─────────────────────────────────────────────
    duplicates = {
        feat_id: linenos
        for feat_id, linenos in id_linenos.items()
        if len(linenos) > 1
    }
    result.duplicate_ids = duplicates

    if duplicates:
        n_dup = len(duplicates)
        examples = list(duplicates.items())[:max_duplicate_examples]
        example_str = "; ".join(
            f"ID={eid!r} (lines {ls[:3]}{'…' if len(ls) > 3 else ''})"
            for eid, ls in examples
        )
        result.warnings.append(
            f"{n_dup} duplicate feature ID(s) detected. "
            f"Examples: {example_str}"
        )
        result.fix_suggestions.append(
            "Duplicate IDs are common in RefSeq GFF3 files where multiple CDS rows "
            "share the same ID. LiftOn will auto-retry with a unique-ID transformation.\n"
            "To inspect duplicates manually:\n"
            + "\n".join(
                f"  grep 'ID={eid}' {file_path} | head"
                for eid, _ in list(duplicates.items())[:5]
            )
        )
        # Duplicate IDs are a WARNING not a hard error — gffutils can handle
        # them if we use the right strategy.  Mark as not 100% valid so the
        # caller knows to use the fallback DB-build path.
        result.is_valid = False

    # ── 6. Orphan-parent check ────────────────────────────────────────────────
    if check_orphan_parents and result.file_format == "GFF3":
        orphans = [p for p in all_parent_refs if p and p not in all_ids]
        result.orphan_parents = orphans[:50]   # cap at 50
        if orphans:
            result.warnings.append(
                f"{len(orphans)} Parent reference(s) point to unknown IDs "
                f"(orphan parents). Examples: {orphans[:5]}"
            )
            result.fix_suggestions.append(
                "Some features reference Parent IDs that do not exist in the file. "
                "This can cause missing genes in the output.\n"
                "Check for typos or missing parent features with:\n"
                + "\n".join(
                    f"  grep 'ID={p}' {file_path}" for p in orphans[:3]
                )
            )

    return result


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
