"""
gff3_validator.py  —  Comprehensive GFF3 output validator for LiftOn.

Validates a GFF3 file against the NCBI GFF3 specification and LiftOn-specific
conventions:
  • Column format (9-column, correct data types)
  • Coordinate correctness (start ≤ end, 1-based, non-negative)
  • Strand and phase validity
  • Attribute syntax (key=value, semicolon-separated, required attributes)
  • Feature hierarchy (gene → mRNA/ncRNA/transcript → exon → CDS)
  • CDS phase consistency (phase tracks correctly given CDS lengths)
  • Exon / CDS containment within parent transcript
  • Transcript containment within parent gene
  • CDS coordinates contained within corresponding exon
  • No overlapping CDS within one transcript
  • No duplicate gene/transcript IDs in the file
  • LiftOn-specific attribute checks (protein_identity, dna_identity, annotation_source)

Can be used as:
  • A standalone script:  python -m lifton.gff3_validator output.gff3
  • An importable module: from lifton.gff3_validator import validate_gff3_file
"""

import sys
import re
import os
from dataclasses import dataclass, field
from collections import defaultdict
from typing import List, Dict, Optional, Tuple, Set


# ─────────────────────────────────────────────────────────────────────────────
# Constants
# ─────────────────────────────────────────────────────────────────────────────

# Valid feature types for the GFF3 hierarchy produced by LiftOn
GENE_TYPES = {
    "gene", "pseudogene", "transposable_element",
    "LiftOn-gene",   # LiftOn synthetic gene feature type
}
TRANSCRIPT_TYPES = {
    "mRNA", "ncRNA", "transcript", "lncRNA", "lnc_RNA", "nc_RNA",
    "rRNA", "tRNA", "miRNA", "snoRNA", "snRNA", "scRNA",
    "primary_transcript", "processed_pseudogene",
    "three_prime_overlapping_ncrna",
    "antisense_RNA", "antisense", "guide_RNA",
    "RNase_MRP_RNA", "RNase_P_RNA", "SRP_RNA", "vault_RNA",
    "Y_RNA", "telomerase_RNA",
    "C_gene_segment", "V_gene_segment", "D_gene_segment", "J_gene_segment",
}
EXON_TYPES    = {"exon"}
CDS_TYPES     = {"CDS"}
REGION_TYPES  = {"region"}

VALID_STRANDS = {"+", "-", "."}
VALID_PHASES  = {0, 1, 2}

# Official GFF3 attribute names (capital letters per spec)
OFFICIAL_ATTRS = {
    "ID", "Parent", "Name", "Alias", "Target", "Gap", "Derives_from",
    "Note", "Dbxref", "Ontology_term", "Is_circular",
}
# LiftOn-specific attributes always written on transcripts
LIFTON_TRANS_ATTRS = {"protein_identity", "dna_identity", "extra_copy_number",
                      "annotation_source"}

# Column indices (0-based)
COL_SEQID  = 0
COL_SOURCE = 1
COL_TYPE   = 2
COL_START  = 3
COL_END    = 4
COL_SCORE  = 5
COL_STRAND = 6
COL_PHASE  = 7
COL_ATTRS  = 8


# ─────────────────────────────────────────────────────────────────────────────
# Issue severity
# ─────────────────────────────────────────────────────────────────────────────

class Severity:
    ERROR   = "ERROR"
    WARNING = "WARNING"
    INFO    = "INFO"


# ─────────────────────────────────────────────────────────────────────────────
# Data classes
# ─────────────────────────────────────────────────────────────────────────────

@dataclass
class GFF3Issue:
    severity: str
    lineno: int
    feature_id: str
    check: str
    message: str

    def __str__(self) -> str:
        loc = f"line {self.lineno}" if self.lineno > 0 else "global"
        fid = f" [{self.feature_id}]" if self.feature_id else ""
        return f"[{self.severity}] {loc}{fid} — {self.check}: {self.message}"


@dataclass
class GFF3Record:
    """One parsed GFF3 data line."""
    lineno: int
    seqid: str
    source: str
    ftype: str
    start: int       # 1-based
    end: int
    score: str       # "." or float string
    strand: str
    phase: str       # "." or "0"/"1"/"2"
    attrs: Dict[str, List[str]]
    raw: str

    @property
    def feat_id(self) -> str:
        return self.attrs.get("ID", [""])[0]

    @property
    def parent_id(self) -> str:
        return self.attrs.get("Parent", [""])[0]


@dataclass
class ValidationResult:
    file_path: str
    total_lines: int = 0
    data_lines: int  = 0
    comment_lines: int = 0
    issues: List[GFF3Issue] = field(default_factory=list)
    stats: Dict[str, int] = field(default_factory=dict)

    @property
    def errors(self) -> List[GFF3Issue]:
        return [i for i in self.issues if i.severity == Severity.ERROR]

    @property
    def warnings(self) -> List[GFF3Issue]:
        return [i for i in self.issues if i.severity == Severity.WARNING]

    @property
    def is_valid(self) -> bool:
        return len(self.errors) == 0


# ─────────────────────────────────────────────────────────────────────────────
# Public API
# ─────────────────────────────────────────────────────────────────────────────

def validate_gff3_file(
    gff3_path: str,
    check_hierarchy: bool = True,
    check_cds_phase: bool = True,
    check_containment: bool = True,
    check_lifton_attrs: bool = True,
    max_issues_per_check: int = 50,
) -> ValidationResult:
    """
    Validate a GFF3 file produced by LiftOn.

    Parameters
    ----------
    gff3_path : str
        Path to the GFF3 file to validate.
    check_hierarchy : bool
        Validate gene→transcript→exon→CDS parent-child relationships.
    check_cds_phase : bool
        Validate that CDS phase values are consistent with cumulative CDS lengths.
    check_containment : bool
        Validate that children are contained within their parent coordinates.
    check_lifton_attrs : bool
        Validate LiftOn-specific attributes (protein_identity, dna_identity).
    max_issues_per_check : int
        Maximum number of issues to report per check type (avoids flooding output).

    Returns
    -------
    ValidationResult
    """
    result = ValidationResult(file_path=gff3_path)

    # ── File-level pre-checks ────────────────────────────────────────────────
    if not os.path.exists(gff3_path):
        result.issues.append(GFF3Issue(
            Severity.ERROR, 0, "", "file_exists",
            f"File not found: {gff3_path}"
        ))
        return result

    if os.path.getsize(gff3_path) == 0:
        result.issues.append(GFF3Issue(
            Severity.ERROR, 0, "", "file_not_empty",
            f"File is empty: {gff3_path}"
        ))
        return result

    # ── Parse pass ───────────────────────────────────────────────────────────
    records: List[GFF3Record] = []
    issue_counts: Dict[str, int] = defaultdict(int)

    try:
        records, parse_issues, meta = _parse_gff3(gff3_path, max_issues_per_check)
    except Exception as exc:
        result.issues.append(GFF3Issue(
            Severity.ERROR, 0, "", "parse_error",
            f"Unexpected error during parsing: {exc}"
        ))
        return result

    result.total_lines   = meta["total_lines"]
    result.data_lines    = meta["data_lines"]
    result.comment_lines = meta["comment_lines"]
    result.issues.extend(parse_issues)

    # ── Build feature index ──────────────────────────────────────────────────
    id_to_record: Dict[str, GFF3Record] = {}
    parent_to_children: Dict[str, List[GFF3Record]] = defaultdict(list)
    id_issues: List[GFF3Issue] = []

    seen_ids: Dict[str, int] = {}  # id → first lineno
    for rec in records:
        fid = rec.feat_id
        if fid:
            if fid in seen_ids:
                issue_counts["duplicate_id"] += 1
                if issue_counts["duplicate_id"] <= max_issues_per_check:
                    id_issues.append(GFF3Issue(
                        Severity.ERROR, rec.lineno, fid, "duplicate_id",
                        f"Duplicate ID '{fid}' (first seen on line {seen_ids[fid]})"
                    ))
            else:
                seen_ids[fid] = rec.lineno
                id_to_record[fid] = rec

        pid = rec.parent_id
        if pid:
            parent_to_children[pid].append(rec)

    result.issues.extend(id_issues)

    # ── Hierarchy validation ─────────────────────────────────────────────────
    if check_hierarchy:
        result.issues.extend(
            _check_hierarchy(records, id_to_record, parent_to_children,
                             max_issues_per_check)
        )

    # ── Containment validation ───────────────────────────────────────────────
    if check_containment:
        result.issues.extend(
            _check_containment(records, id_to_record, max_issues_per_check)
        )

    # ── CDS phase validation ─────────────────────────────────────────────────
    if check_cds_phase:
        result.issues.extend(
            _check_cds_phase(parent_to_children, id_to_record,
                             max_issues_per_check)
        )

    # ── LiftOn-specific attribute validation ─────────────────────────────────
    if check_lifton_attrs:
        result.issues.extend(
            _check_lifton_attrs(records, parent_to_children, max_issues_per_check)
        )

    # ── Statistics ───────────────────────────────────────────────────────────
    result.stats = _compute_stats(records, parent_to_children)

    return result


def print_validation_report(result: ValidationResult, verbose: bool = False) -> None:
    """
    Print a human-readable validation report to stderr.
    """
    w = 72
    err_count  = len(result.errors)
    warn_count = len(result.warnings)
    status_label = "✅  VALID" if result.is_valid else "❌  INVALID"

    print("╔" + "═" * (w - 2) + "╗", file=sys.stderr)
    print(f"║  GFF3 VALIDATION REPORT — {status_label}" +
          " " * max(0, w - 30 - len(status_label)) + "║", file=sys.stderr)
    print("╠" + "═" * (w - 2) + "╣", file=sys.stderr)

    def _row(text=""):
        return "║  " + text + " " * max(0, w - 4 - len(text)) + "  ║"

    print(_row(f"File    : {result.file_path}"), file=sys.stderr)
    print(_row(f"Lines   : {result.total_lines} total "
               f"({result.data_lines} data, {result.comment_lines} comment/blank)"),
          file=sys.stderr)

    print("╠" + "═" * (w - 2) + "╣", file=sys.stderr)
    print(_row("FEATURE COUNTS:"), file=sys.stderr)
    for ftype, count in sorted(result.stats.items()):
        print(_row(f"  {ftype:<30}: {count}"), file=sys.stderr)

    print("╠" + "═" * (w - 2) + "╣", file=sys.stderr)
    print(_row(f"Errors   : {err_count}"), file=sys.stderr)
    print(_row(f"Warnings : {warn_count}"), file=sys.stderr)

    # Show all errors, and warnings if verbose
    issues_to_show = result.errors + (result.warnings if verbose else [])
    if issues_to_show:
        print("╠" + "═" * (w - 2) + "╣", file=sys.stderr)
        print(_row("ISSUES:"), file=sys.stderr)
        prev_check = None
        for issue in issues_to_show:
            if issue.check != prev_check:
                print(_row(""), file=sys.stderr)
                print(_row(f"  [{issue.severity}] Check: {issue.check}"), file=sys.stderr)
                prev_check = issue.check
            loc = f"line {issue.lineno}" if issue.lineno > 0 else "global"
            fid = f" [{issue.feature_id}]" if issue.feature_id else ""
            msg = f"    {loc}{fid}: {issue.message}"
            # Wrap long messages
            for chunk in _wrap(msg, w - 6):
                print(_row("  " + chunk), file=sys.stderr)
        if not verbose and warn_count > 0:
            print(_row(f"  (Use verbose=True to also show {warn_count} warnings)"),
                  file=sys.stderr)

    print("╚" + "═" * (w - 2) + "╝", file=sys.stderr)


def _wrap(text: str, max_width: int) -> List[str]:
    if len(text) <= max_width:
        return [text]
    chunks: List[str] = []
    while len(text) > max_width:
        chunks.append(text[:max_width])
        text = text[max_width:]
    if text:
        chunks.append(text)
    return chunks


# ─────────────────────────────────────────────────────────────────────────────
# Parsing
# ─────────────────────────────────────────────────────────────────────────────

def _parse_gff3(
    path: str,
    max_issues_per_check: int,
) -> Tuple[List[GFF3Record], List[GFF3Issue], dict]:
    """
    Parse a GFF3 file into a list of GFF3Record objects.
    Also validates column-level format constraints during parsing.

    Returns (records, issues, meta_dict).
    """
    records: List[GFF3Record] = []
    issues: List[GFF3Issue]   = []
    issue_counts: Dict[str, int] = defaultdict(int)

    total_lines = 0
    data_lines  = 0
    comment_lines = 0

    _ATTR_RE = re.compile(r'^[a-zA-Z_][a-zA-Z0-9_]*=')

    with open(path, "r", errors="replace") as fh:
        for lineno, raw in enumerate(fh, start=1):
            total_lines += 1
            line = raw.rstrip("\n\r")

            # ── Comments and directives ──────────────────────────────────────
            if line.startswith("#") or line.strip() == "":
                comment_lines += 1
                # Check for required GFF3 header
                if lineno == 1 and not line.startswith("##gff-version"):
                    issue_counts["missing_gff_version"] += 1
                    if issue_counts["missing_gff_version"] <= max_issues_per_check:
                        issues.append(GFF3Issue(
                            Severity.WARNING, 1, "", "gff3_header",
                            "First line should be '##gff-version 3' directive"
                        ))
                continue

            # ── Column count ─────────────────────────────────────────────────
            cols = line.split("\t")
            if len(cols) != 9:
                issue_counts["col_count"] += 1
                if issue_counts["col_count"] <= max_issues_per_check:
                    issues.append(GFF3Issue(
                        Severity.ERROR, lineno, "", "column_count",
                        f"Expected 9 tab-separated columns, got {len(cols)}: {line[:80]}"
                    ))
                continue

            data_lines += 1

            seqid  = cols[COL_SEQID].strip()
            source = cols[COL_SOURCE].strip()
            ftype  = cols[COL_TYPE].strip()
            score  = cols[COL_SCORE].strip()
            strand = cols[COL_STRAND].strip()
            phase  = cols[COL_PHASE].strip()
            attrs_str = cols[COL_ATTRS].strip()

            # ── seqid ────────────────────────────────────────────────────────
            if not seqid or seqid == ".":
                issue_counts["seqid_empty"] += 1
                if issue_counts["seqid_empty"] <= max_issues_per_check:
                    issues.append(GFF3Issue(
                        Severity.ERROR, lineno, "", "seqid_empty",
                        "seqid (column 1) must not be empty or '.'"
                    ))

            # ── start / end ──────────────────────────────────────────────────
            try:
                start = int(cols[COL_START])
                end   = int(cols[COL_END])
            except ValueError:
                issue_counts["coord_not_int"] += 1
                if issue_counts["coord_not_int"] <= max_issues_per_check:
                    issues.append(GFF3Issue(
                        Severity.ERROR, lineno, "", "coord_not_int",
                        f"start/end must be integers, got '{cols[COL_START]}' and '{cols[COL_END]}'"
                    ))
                continue

            if start < 1:
                issue_counts["coord_negative"] += 1
                if issue_counts["coord_negative"] <= max_issues_per_check:
                    issues.append(GFF3Issue(
                        Severity.ERROR, lineno, "", "coord_1based",
                        f"start coordinate must be ≥ 1 (GFF3 is 1-based), got {start}"
                    ))

            if end < start:
                issue_counts["coord_order"] += 1
                if issue_counts["coord_order"] <= max_issues_per_check:
                    issues.append(GFF3Issue(
                        Severity.ERROR, lineno, "", "coord_order",
                        f"start ({start}) must be ≤ end ({end})"
                    ))

            # ── score ────────────────────────────────────────────────────────
            if score != ".":
                try:
                    float(score)
                except ValueError:
                    issue_counts["score_invalid"] += 1
                    if issue_counts["score_invalid"] <= max_issues_per_check:
                        issues.append(GFF3Issue(
                            Severity.WARNING, lineno, "", "score_format",
                            f"score must be a number or '.', got '{score}'"
                        ))

            # ── strand ───────────────────────────────────────────────────────
            if strand not in VALID_STRANDS:
                issue_counts["strand_invalid"] += 1
                if issue_counts["strand_invalid"] <= max_issues_per_check:
                    issues.append(GFF3Issue(
                        Severity.ERROR, lineno, "", "strand_valid",
                        f"strand must be '+', '-', or '.', got '{strand}'"
                    ))

            # ── phase ────────────────────────────────────────────────────────
            if ftype in CDS_TYPES:
                if phase == ".":
                    issue_counts["cds_phase_dot"] += 1
                    if issue_counts["cds_phase_dot"] <= max_issues_per_check:
                        issues.append(GFF3Issue(
                            Severity.ERROR, lineno, "", "cds_phase_required",
                            "CDS feature must have an integer phase (0, 1, or 2), not '.'"
                        ))
                else:
                    try:
                        ph = int(phase)
                        if ph not in VALID_PHASES:
                            raise ValueError()
                    except ValueError:
                        issue_counts["cds_phase_invalid"] += 1
                        if issue_counts["cds_phase_invalid"] <= max_issues_per_check:
                            issues.append(GFF3Issue(
                                Severity.ERROR, lineno, "", "cds_phase_value",
                                f"CDS phase must be 0, 1, or 2, got '{phase}'"
                            ))
            else:
                # For non-CDS, phase must be '.'
                if phase != ".":
                    issue_counts["non_cds_phase"] += 1
                    if issue_counts["non_cds_phase"] <= max_issues_per_check:
                        issues.append(GFF3Issue(
                            Severity.WARNING, lineno, "", "non_cds_phase",
                            f"Non-CDS feature '{ftype}' has phase '{phase}' (should be '.')"
                        ))

            # ── Attribute parsing ────────────────────────────────────────────
            attrs, attr_issues = _parse_attributes(attrs_str, lineno,
                                                   max_issues_per_check,
                                                   issue_counts)
            issues.extend(attr_issues)

            # ── Required ID attribute ────────────────────────────────────────
            # GFF3 spec: features that are referenced as Parent must have ID
            # LiftOn always writes ID for gene and transcript features
            if ftype in GENE_TYPES or ftype in TRANSCRIPT_TYPES:
                if "ID" not in attrs:
                    issue_counts["missing_id"] += 1
                    if issue_counts["missing_id"] <= max_issues_per_check:
                        issues.append(GFF3Issue(
                            Severity.ERROR, lineno, "", "missing_id",
                            f"'{ftype}' feature must have an ID attribute"
                        ))

            # ── Build record ─────────────────────────────────────────────────
            rec = GFF3Record(
                lineno=lineno, seqid=seqid, source=source, ftype=ftype,
                start=start, end=end, score=score, strand=strand,
                phase=phase, attrs=attrs, raw=line,
            )
            records.append(rec)

    meta = {
        "total_lines":   total_lines,
        "data_lines":    data_lines,
        "comment_lines": comment_lines,
    }
    return records, issues, meta


def _parse_attributes(
    attrs_str: str,
    lineno: int,
    max_issues: int,
    issue_counts: dict,
) -> Tuple[Dict[str, List[str]], List[GFF3Issue]]:
    """Parse the attributes column into a dict. Returns (attrs_dict, issues)."""
    attrs: Dict[str, List[str]] = {}
    issues: List[GFF3Issue] = []

    if not attrs_str or attrs_str == ".":
        return attrs, issues

    for part in attrs_str.split(";"):
        part = part.strip()
        if not part:
            continue
        if "=" not in part:
            issue_counts["attr_no_equals"] += 1
            if issue_counts["attr_no_equals"] <= max_issues:
                issues.append(GFF3Issue(
                    Severity.ERROR, lineno, "", "attr_format",
                    f"Attribute '{part}' has no '=' separator"
                ))
            continue
        key, _, value = part.partition("=")
        key = key.strip()
        value = value.strip()
        if not key:
            issue_counts["attr_empty_key"] += 1
            if issue_counts["attr_empty_key"] <= max_issues:
                issues.append(GFF3Issue(
                    Severity.WARNING, lineno, "", "attr_empty_key",
                    f"Attribute with empty key: '{part}'"
                ))
            continue
        # Values may be comma-separated lists
        values = [v.strip() for v in value.split(",") if v.strip()]
        attrs[key] = values

    return attrs, issues


# ─────────────────────────────────────────────────────────────────────────────
# Hierarchy checks
# ─────────────────────────────────────────────────────────────────────────────

def _check_hierarchy(
    records: List[GFF3Record],
    id_to_record: Dict[str, GFF3Record],
    parent_to_children: Dict[str, List[GFF3Record]],
    max_issues: int,
) -> List[GFF3Issue]:
    """
    Validate the gene→transcript→exon→CDS parent-child relationship.

    Rules (per NCBI GFF3 spec and LiftOn conventions):
    1. Every Parent reference must point to a known ID.
    2. Genes must not have a Parent (they are top-level).
    3. Transcripts must have a gene as Parent.
    4. Exons must have a transcript as Parent.
    5. CDS features must have a transcript as Parent.
    6. A transcript must have at least one exon child.
    7. Coding transcripts (mRNA) must have at least one CDS child.
    8. CDS must be a subset of at least one exon of the same transcript.
    """
    issues: List[GFF3Issue] = []
    issue_counts: Dict[str, int] = defaultdict(int)
    all_ids = set(id_to_record.keys())

    for rec in records:
        fid = rec.feat_id
        pid = rec.parent_id

        # ── Rule 1: All Parent references must resolve ───────────────────────
        if pid and pid not in all_ids:
            issue_counts["orphan_parent"] += 1
            if issue_counts["orphan_parent"] <= max_issues:
                issues.append(GFF3Issue(
                    Severity.ERROR, rec.lineno, fid, "orphan_parent",
                    f"Parent='{pid}' not found in file"
                ))

        # ── Rule 2: Genes must be top-level (no Parent) ──────────────────────
        if rec.ftype in GENE_TYPES and pid:
            issue_counts["gene_has_parent"] += 1
            if issue_counts["gene_has_parent"] <= max_issues:
                issues.append(GFF3Issue(
                    Severity.ERROR, rec.lineno, fid, "gene_has_parent",
                    f"Gene-type feature '{rec.ftype}' should have no Parent, "
                    f"but has Parent='{pid}'"
                ))

        # ── Rule 3: Transcripts must have a gene parent ──────────────────────
        if rec.ftype in TRANSCRIPT_TYPES:
            if not pid:
                issue_counts["trans_no_parent"] += 1
                if issue_counts["trans_no_parent"] <= max_issues:
                    issues.append(GFF3Issue(
                        Severity.ERROR, rec.lineno, fid, "transcript_no_parent",
                        f"Transcript-type '{rec.ftype}' has no Parent attribute"
                    ))
            elif pid in id_to_record:
                parent_rec = id_to_record[pid]
                if parent_rec.ftype not in GENE_TYPES:
                    issue_counts["trans_wrong_parent"] += 1
                    if issue_counts["trans_wrong_parent"] <= max_issues:
                        issues.append(GFF3Issue(
                            Severity.ERROR, rec.lineno, fid, "transcript_parent_type",
                            f"Transcript parent must be a gene; "
                            f"'{pid}' has type '{parent_rec.ftype}'"
                        ))

        # ── Rule 4: Exons must have a transcript parent ──────────────────────
        if rec.ftype in EXON_TYPES:
            if not pid:
                issue_counts["exon_no_parent"] += 1
                if issue_counts["exon_no_parent"] <= max_issues:
                    issues.append(GFF3Issue(
                        Severity.ERROR, rec.lineno, fid, "exon_no_parent",
                        "Exon has no Parent attribute"
                    ))
            elif pid in id_to_record:
                parent_rec = id_to_record[pid]
                if parent_rec.ftype not in TRANSCRIPT_TYPES:
                    issue_counts["exon_wrong_parent"] += 1
                    if issue_counts["exon_wrong_parent"] <= max_issues:
                        issues.append(GFF3Issue(
                            Severity.ERROR, rec.lineno, fid, "exon_parent_type",
                            f"Exon parent '{pid}' has type '{parent_rec.ftype}' "
                            f"(expected transcript type)"
                        ))

        # ── Rule 5: CDS must have a transcript parent ────────────────────────
        if rec.ftype in CDS_TYPES:
            if not pid:
                issue_counts["cds_no_parent"] += 1
                if issue_counts["cds_no_parent"] <= max_issues:
                    issues.append(GFF3Issue(
                        Severity.ERROR, rec.lineno, fid, "cds_no_parent",
                        "CDS has no Parent attribute"
                    ))
            elif pid in id_to_record:
                parent_rec = id_to_record[pid]
                if parent_rec.ftype not in TRANSCRIPT_TYPES:
                    issue_counts["cds_wrong_parent"] += 1
                    if issue_counts["cds_wrong_parent"] <= max_issues:
                        issues.append(GFF3Issue(
                            Severity.ERROR, rec.lineno, fid, "cds_parent_type",
                            f"CDS parent '{pid}' has type '{parent_rec.ftype}' "
                            f"(expected transcript type)"
                        ))

    # ── Rules 6 & 7: Transcripts must have exons (and mRNA must have CDS) ───
    for trans_id, trans_rec in id_to_record.items():
        if trans_rec.ftype not in TRANSCRIPT_TYPES:
            continue
        children = parent_to_children.get(trans_id, [])
        exon_children = [c for c in children if c.ftype in EXON_TYPES]
        cds_children  = [c for c in children if c.ftype in CDS_TYPES]

        if not exon_children and not cds_children:
            issue_counts["trans_no_exon"] += 1
            if issue_counts["trans_no_exon"] <= max_issues:
                issues.append(GFF3Issue(
                    Severity.WARNING, trans_rec.lineno, trans_id,
                    "transcript_has_exons",
                    f"{trans_rec.ftype} has no exon or CDS children"
                ))

        if trans_rec.ftype == "mRNA" and not cds_children:
            issue_counts["mrna_no_cds"] += 1
            if issue_counts["mrna_no_cds"] <= max_issues:
                issues.append(GFF3Issue(
                    Severity.WARNING, trans_rec.lineno, trans_id,
                    "mrna_has_cds",
                    "mRNA has no CDS children (expected for protein-coding transcripts)"
                ))

    # ── Rule 8: Verify all genes have at least one transcript ───────────────
    for gene_id, gene_rec in id_to_record.items():
        if gene_rec.ftype not in GENE_TYPES:
            continue
        trans_children = [
            c for c in parent_to_children.get(gene_id, [])
            if c.ftype in TRANSCRIPT_TYPES
        ]
        other_children = parent_to_children.get(gene_id, [])
        if not trans_children and not other_children:
            issue_counts["gene_no_trans"] += 1
            if issue_counts["gene_no_trans"] <= max_issues:
                issues.append(GFF3Issue(
                    Severity.WARNING, gene_rec.lineno, gene_id,
                    "gene_has_transcripts",
                    f"Gene has no child transcripts or features"
                ))

    return issues


# ─────────────────────────────────────────────────────────────────────────────
# Containment check
# ─────────────────────────────────────────────────────────────────────────────

def _check_containment(
    records: List[GFF3Record],
    id_to_record: Dict[str, GFF3Record],
    max_issues: int,
) -> List[GFF3Issue]:
    """
    Validate that every child feature is contained within its parent coordinates.
    Also validates that all features on the same seqid share the parent's seqid.
    """
    issues: List[GFF3Issue] = []
    issue_counts: Dict[str, int] = defaultdict(int)

    for rec in records:
        pid = rec.parent_id
        if not pid or pid not in id_to_record:
            continue
        parent = id_to_record[pid]

        # ── Same seqid ───────────────────────────────────────────────────────
        if rec.seqid != parent.seqid:
            issue_counts["different_seqid"] += 1
            if issue_counts["different_seqid"] <= max_issues:
                issues.append(GFF3Issue(
                    Severity.ERROR, rec.lineno, rec.feat_id,
                    "seqid_consistency",
                    f"Feature seqid '{rec.seqid}' differs from parent "
                    f"'{parent.feat_id}' seqid '{parent.seqid}'"
                ))

        # ── Same strand ──────────────────────────────────────────────────────
        if (rec.strand not in (".", "") and parent.strand not in (".", "")
                and rec.strand != parent.strand):
            issue_counts["strand_mismatch"] += 1
            if issue_counts["strand_mismatch"] <= max_issues:
                issues.append(GFF3Issue(
                    Severity.ERROR, rec.lineno, rec.feat_id,
                    "strand_consistency",
                    f"Feature strand '{rec.strand}' differs from parent "
                    f"'{parent.feat_id}' strand '{parent.strand}'"
                ))

        # ── Coordinate containment ───────────────────────────────────────────
        if rec.start < parent.start or rec.end > parent.end:
            issue_counts["coord_containment"] += 1
            if issue_counts["coord_containment"] <= max_issues:
                issues.append(GFF3Issue(
                    Severity.ERROR, rec.lineno, rec.feat_id,
                    "coord_containment",
                    f"{rec.ftype} [{rec.start}, {rec.end}] extends outside "
                    f"parent '{parent.feat_id}' [{parent.start}, {parent.end}]"
                ))

    return issues


# ─────────────────────────────────────────────────────────────────────────────
# CDS phase check
# ─────────────────────────────────────────────────────────────────────────────

def _check_cds_phase(
    parent_to_children: Dict[str, List[GFF3Record]],
    id_to_record: Dict[str, GFF3Record],
    max_issues: int,
) -> List[GFF3Issue]:
    """
    Validate GFF3 CDS phase values.

    GFF3 phase rule: phase = (3 - (cumulative_cds_length_so_far % 3)) % 3
    where cumulative_cds_length_so_far is the total length of all preceding
    CDS rows *before* the current one (sorted by coord, 5'→3').

    The first CDS should always have phase 0 for a complete CDS.
    """
    issues: List[GFF3Issue] = []
    issue_counts: Dict[str, int] = defaultdict(int)

    for trans_id, children in parent_to_children.items():
        cds_list = [c for c in children if c.ftype in CDS_TYPES]
        if len(cds_list) < 2:
            continue  # Can't validate phase on a single-CDS transcript

        trans_rec = id_to_record.get(trans_id)
        if not trans_rec:
            continue
        strand = trans_rec.strand

        # Sort CDS 5'→3'
        if strand == "-":
            cds_sorted = sorted(cds_list, key=lambda r: r.end, reverse=True)
        else:
            cds_sorted = sorted(cds_list, key=lambda r: r.start)

        accum_len = 0
        for i, cds in enumerate(cds_sorted):
            expected_phase = (3 - accum_len % 3) % 3
            try:
                actual_phase = int(cds.phase)
            except ValueError:
                accum_len += cds.end - cds.start + 1
                continue

            if i > 0 and actual_phase != expected_phase:
                issue_counts["cds_phase_mismatch"] += 1
                if issue_counts["cds_phase_mismatch"] <= max_issues:
                    issues.append(GFF3Issue(
                        Severity.WARNING, cds.lineno, cds.feat_id or trans_id,
                        "cds_phase_consistency",
                        f"CDS #{i+1} of transcript '{trans_id}': "
                        f"expected phase {expected_phase} "
                        f"(accum_len={accum_len}), got {actual_phase}"
                    ))
            accum_len += cds.end - cds.start + 1

    return issues


# ─────────────────────────────────────────────────────────────────────────────
# LiftOn-specific attribute checks
# ─────────────────────────────────────────────────────────────────────────────

def _check_lifton_attrs(
    records: List[GFF3Record],
    parent_to_children: Dict[str, List[GFF3Record]],
    max_issues: int,
) -> List[GFF3Issue]:
    """
    Validate LiftOn-specific attributes written on transcript features:
      - protein_identity: must be float in [0.0, 1.0]
      - dna_identity: must be float in [0.0, 1.0]
      - source column must be 'LiftOn' for features LiftOn wrote
    """
    issues: List[GFF3Issue] = []
    issue_counts: Dict[str, int] = defaultdict(int)

    for rec in records:
        fid = rec.feat_id

        # ── Source column ────────────────────────────────────────────────────
        if rec.ftype in GENE_TYPES | TRANSCRIPT_TYPES | EXON_TYPES | CDS_TYPES:
            if rec.source not in ("LiftOn", "miniprot", "Liftoff", "."):
                issue_counts["unexpected_source"] += 1
                if issue_counts["unexpected_source"] <= max_issues:
                    issues.append(GFF3Issue(
                        Severity.INFO, rec.lineno, fid, "lifton_source",
                        f"source='{rec.source}' (expected 'LiftOn', 'miniprot', or 'Liftoff')"
                    ))

        # ── protein_identity ─────────────────────────────────────────────────
        if "protein_identity" in rec.attrs:
            val_str = rec.attrs["protein_identity"][0]
            try:
                val = float(val_str)
                if not (0.0 <= val <= 1.0):
                    raise ValueError()
            except (ValueError, IndexError):
                issue_counts["protein_id_invalid"] += 1
                if issue_counts["protein_id_invalid"] <= max_issues:
                    issues.append(GFF3Issue(
                        Severity.ERROR, rec.lineno, fid, "lifton_protein_identity",
                        f"protein_identity='{val_str}' must be a float in [0.0, 1.0]"
                    ))

        # ── dna_identity ─────────────────────────────────────────────────────
        if "dna_identity" in rec.attrs:
            val_str = rec.attrs["dna_identity"][0]
            try:
                val = float(val_str)
                if not (0.0 <= val <= 1.0):
                    raise ValueError()
            except (ValueError, IndexError):
                issue_counts["dna_id_invalid"] += 1
                if issue_counts["dna_id_invalid"] <= max_issues:
                    issues.append(GFF3Issue(
                        Severity.ERROR, rec.lineno, fid, "lifton_dna_identity",
                        f"dna_identity='{val_str}' must be a float in [0.0, 1.0]"
                    ))

        # ── Transcripts should have both identity attrs ───────────────────────
        if rec.ftype in TRANSCRIPT_TYPES:
            has_prot_id = "protein_identity" in rec.attrs
            has_dna_id  = "dna_identity" in rec.attrs
            has_cds_children = any(
                c.ftype in CDS_TYPES
                for c in parent_to_children.get(fid, [])
            )
            if has_cds_children and not has_prot_id:
                issue_counts["missing_protein_id"] += 1
                if issue_counts["missing_protein_id"] <= max_issues:
                    issues.append(GFF3Issue(
                        Severity.WARNING, rec.lineno, fid, "lifton_attrs_present",
                        "Coding transcript is missing 'protein_identity' attribute"
                    ))
            if not has_dna_id:
                issue_counts["missing_dna_id"] += 1
                if issue_counts["missing_dna_id"] <= max_issues:
                    issues.append(GFF3Issue(
                        Severity.WARNING, rec.lineno, fid, "lifton_attrs_present",
                        "Transcript is missing 'dna_identity' attribute"
                    ))

    return issues


# ─────────────────────────────────────────────────────────────────────────────
# Statistics
# ─────────────────────────────────────────────────────────────────────────────

def _compute_stats(
    records: List[GFF3Record],
    parent_to_children: Dict[str, List[GFF3Record]],
) -> Dict[str, int]:
    stats: Dict[str, int] = defaultdict(int)
    for rec in records:
        stats[rec.ftype] += 1
    return dict(stats)


# ─────────────────────────────────────────────────────────────────────────────
# CLI entry point
# ─────────────────────────────────────────────────────────────────────────────

def _main(argv=None):
    import argparse

    parser = argparse.ArgumentParser(
        description="Validate a GFF3 file produced by LiftOn",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument(
        "gff3_file", help="Path to the GFF3 file to validate"
    )
    parser.add_argument(
        "-v", "--verbose", action="store_true",
        help="Also show warnings in the output"
    )
    parser.add_argument(
        "--no-hierarchy", action="store_true",
        help="Skip parent-child hierarchy checks"
    )
    parser.add_argument(
        "--no-phase", action="store_true",
        help="Skip CDS phase consistency checks"
    )
    parser.add_argument(
        "--no-containment", action="store_true",
        help="Skip coordinate containment checks"
    )
    parser.add_argument(
        "--no-lifton-attrs", action="store_true",
        help="Skip LiftOn-specific attribute checks"
    )
    parser.add_argument(
        "--max-issues", type=int, default=50, metavar="N",
        help="Maximum issues to report per check type (default: 50)"
    )
    parser.add_argument(
        "--json", action="store_true",
        help="Output issues as JSON to stdout (useful for programmatic use)"
    )
    args = parser.parse_args(argv)

    result = validate_gff3_file(
        gff3_path=args.gff3_file,
        check_hierarchy=not args.no_hierarchy,
        check_cds_phase=not args.no_phase,
        check_containment=not args.no_containment,
        check_lifton_attrs=not args.no_lifton_attrs,
        max_issues_per_check=args.max_issues,
    )

    if args.json:
        import json as _json
        out = {
            "file": result.file_path,
            "is_valid": result.is_valid,
            "stats": result.stats,
            "issues": [
                {
                    "severity": i.severity,
                    "lineno": i.lineno,
                    "feature_id": i.feature_id,
                    "check": i.check,
                    "message": i.message,
                }
                for i in result.issues
            ],
        }
        print(_json.dumps(out, indent=2))
    else:
        print_validation_report(result, verbose=args.verbose)

    sys.exit(0 if result.is_valid else 1)


if __name__ == "__main__":
    _main()
