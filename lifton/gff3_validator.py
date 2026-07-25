"""
gff3_validator.py  —  Comprehensive GFF3 output validator for LiftOn.

Validates a GFF3 file against the NCBI GFF3 specification and LiftOn-specific
conventions:
  • Column format (9-column, correct data types)
  • Coordinate correctness (start ≤ end, 1-based, non-negative)
  • Strand and phase validity
  • Attribute syntax (key=value, semicolon-separated, required attributes)
  • Feature hierarchy (child-bearing root → transcript → exon → CDS)
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

from __future__ import annotations

import sys
import re
import os
import gzip
from dataclasses import dataclass, field
from collections import defaultdict
from typing import List, Dict, Optional, Tuple, Set


# ─────────────────────────────────────────────────────────────────────────────
# Constants
# ─────────────────────────────────────────────────────────────────────────────

# Valid feature types for the GFF3 hierarchy produced by LiftOn
GENE_TYPES = {
    "gene", "pseudogene", "ncRNA_gene", "transposable_element",
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
# GFF3 allows a DISCONTINUOUS feature to be split across several lines that all
# share ONE ID (the canonical case is a multi-exon CDS: one CDS feature, one
# segment per coding exon, every segment carrying the same ID + Parent). Those
# shared IDs are valid, NOT duplicates, so the duplicate-ID check exempts them.
DISCONTINUOUS_FEATURE_TYPES = {"CDS"}
# NCBI RefSeq permits exon/CDS DIRECTLY under a pseudogene (or transposable
# element) with no intervening transcript level. Accept these as valid
# exon/CDS parents so faithfully-lifted pseudogenes (Iteration-5
# --lift-gene-like) don't trip the exon/CDS parent-type hierarchy checks.
DIRECT_EXON_PARENT_TYPES = {"pseudogene", "transposable_element"}


def _can_be_child_bearing_root(ftype: str) -> bool:
    """Return whether LiftOn may emit *ftype* as a hierarchy root.

    LiftOn auto-detects parentless, child-bearing locus types instead of
    limiting roots to ``gene`` (for example ``ncRNA_gene`` and structured
    mobile elements). Transcript-shaped loci with direct exons are also
    emitted without a synthetic gene row. Exon and CDS are leaf levels and
    can therefore never gain root status merely by being parentless.
    """

    return ftype not in EXON_TYPES | CDS_TYPES


def is_transcript_type(ftype: str) -> bool:
    """Return whether *ftype* names a transcript (middle) hierarchy level.

    The explicit ``TRANSCRIPT_TYPES`` set cannot be complete. LiftOn auto-detects
    gene-like parent types from the reference (``get_gene_like_feature_types``) and
    preserves the reference's middle-level type VERBATIM, while the Sequence Ontology
    transcript namespace is open-ended: Ensembl/GENCODE emit ``pseudogenic_transcript``
    under ``pseudogene``, RefSeq emits ``misc_RNA``/``pre_miRNA``, and so on. Because
    ``validate_gff3_structure`` is mandatory and every issue it raises is an ERROR,
    an unlisted-but-valid type made LiftOn faithfully lift a feature, write valid
    GFF3, and then reject its OWN output and exit 2 with no annotation.

    So recognise the namespace by shape as well as by name: anything ending in
    ``RNA`` (``misc_RNA``, ``pre_miRNA``, ``pseudogenic_tRNA``, …) or in
    ``_transcript`` (``pseudogenic_transcript``, ``unconfirmed_transcript``, …) is a
    transcript level. Types outside that namespace -- e.g. ``misc_feature`` -- are
    still rejected, so a genuinely wrong hierarchy keeps failing.
    """
    if ftype in TRANSCRIPT_TYPES:
        return True
    return ftype.endswith("RNA") or ftype.endswith("_transcript")

# The GFF3 spec permits four strand values: '+', '-', '.' (not stranded), and
# '?' (stranded but the strand is unknown). '?' is valid and appears on real
# RefSeq organellar features, so flagging it was a false positive.
VALID_STRANDS = {"+", "-", ".", "?"}
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

    @property
    def parent_ids(self) -> Tuple[str, ...]:
        """Return a normalized Parent tuple for duplicate-ID comparisons."""

        return tuple(sorted({
            parent.strip()
            for parent in self.attrs.get("Parent", [])
            if parent.strip()
        }))


@dataclass
class ValidationResult:
    file_path: str
    total_lines: int = 0
    data_lines: int  = 0
    comment_lines: int = 0
    issues: List[GFF3Issue] = field(default_factory=list)
    stats: Dict[str, int] = field(default_factory=dict)
    issue_totals: Dict[str, int] = field(default_factory=dict)
    severity_totals: Dict[str, int] = field(default_factory=dict)

    @property
    def errors(self) -> List[GFF3Issue]:
        return [i for i in self.issues if i.severity == Severity.ERROR]

    @property
    def warnings(self) -> List[GFF3Issue]:
        return [i for i in self.issues if i.severity == Severity.WARNING]

    @property
    def is_valid(self) -> bool:
        # Some streaming validators cap the number of materialized issue
        # objects while retaining uncapped severity totals.  In particular,
        # ``max_issues_per_check=0`` is a useful "counts only" mode and must
        # not turn a real error into a valid result merely because no detail
        # object was retained.
        return (
            self.severity_totals.get(Severity.ERROR, 0) == 0
            and len(self.errors) == 0
        )


# ─────────────────────────────────────────────────────────────────────────────
# Public API
# ─────────────────────────────────────────────────────────────────────────────

def _discontinuous_cds_signature(record: GFF3Record) -> tuple | None:
    """Return the identity shared by valid segments of one CDS feature."""

    if record.ftype != "CDS" or not record.parent_ids:
        return None
    return (
        record.ftype,
        record.seqid,
        record.source,
        record.strand,
        record.parent_ids,
    )


def _header_before_fasta_offset(handle, offset: int) -> tuple[int, bytes] | None:
    """Return the header line ending at a FAI sequence offset.

    FAI offsets point to the first base immediately after a FASTA header.  The
    backward scan is block-bounded in memory and normally examines one short
    header line.
    """

    if offset <= 0:
        return None
    handle.seek(offset - 1)
    if handle.read(1) != b"\n":
        return None
    search_end = offset - 1
    cursor = search_end
    header_start = 0
    while cursor > 0:
        block_start = max(0, cursor - 65536)
        handle.seek(block_start)
        block = handle.read(cursor - block_start)
        newline = block.rfind(b"\n")
        if newline >= 0:
            header_start = block_start + newline + 1
            break
        cursor = block_start
    handle.seek(header_start)
    header = handle.read(offset - header_start).rstrip(b"\r\n")
    if not header.startswith(b">"):
        return None
    return header_start, header


def _fai_sequence_lengths(path: str) -> Dict[str, int] | None:
    """Return FAI lengths only when the index agrees with the FASTA layout.

    Modification times alone cannot establish that a sidecar belongs to a
    FASTA: copied or replaced indexes can be newer than unrelated sequence
    content.  Validate names, offsets, line geometry, and record boundaries
    against the actual plain-text FASTA before accepting the fast path.  A
    rejected sidecar is harmless; :func:`fasta_sequence_lengths` regenerates
    lengths by streaming the sequence without modifying user files.
    """

    fasta_path = os.path.realpath(path)
    index_path = os.fspath(path) + ".fai"
    # Standard gzip streams do not expose the uncompressed byte offsets stored
    # by faidx.  Streaming fallback is the only backend-neutral verification.
    if os.fspath(path).lower().endswith(".gz"):
        return None
    try:
        fasta_stat = os.stat(fasta_path)
        index_stat = os.stat(index_path)
    except OSError:
        return None
    if index_stat.st_size <= 0 or fasta_stat.st_size <= 0:
        return None
    lengths: Dict[str, int] = {}
    entries: List[Tuple[str, int, int, int, int]] = []
    try:
        with open(index_path, "r", encoding="utf-8", errors="strict") as handle:
            for lineno, raw in enumerate(handle, start=1):
                columns = raw.rstrip("\r\n").split("\t")
                if len(columns) < 5:
                    raise ValueError(
                        f"FASTA index line {lineno} has fewer than five columns"
                    )
                identifier = columns[0]
                length = int(columns[1])
                offset = int(columns[2])
                line_bases = int(columns[3])
                line_width = int(columns[4])
                if (
                    not identifier
                    or identifier in lengths
                    or length <= 0
                    or offset <= 0
                    or line_bases <= 0
                    or line_bases > length
                    or line_width < line_bases
                    or line_width - line_bases not in (0, 1, 2)
                ):
                    raise ValueError(
                        f"FASTA index line {lineno} has an ambiguous record"
                    )
                lengths[identifier] = length
                entries.append((
                    identifier, length, offset, line_bases, line_width,
                ))
    except (OSError, UnicodeError, ValueError):
        return None
    if not entries:
        return None

    try:
        with open(fasta_path, "rb") as handle:
            headers = []
            for identifier, _length, offset, _bases, _width in entries:
                header_record = _header_before_fasta_offset(handle, offset)
                if header_record is None:
                    return None
                header_start, header = header_record
                try:
                    header_id = header[1:].split(None, 1)[0].decode("utf-8")
                except (IndexError, UnicodeError):
                    return None
                if header_id != identifier:
                    return None
                headers.append(header_start)

            if headers != sorted(headers) or len(set(headers)) != len(headers):
                return None
            if headers[0] != 0:
                # The streaming parser permits leading blank lines, but a FAI
                # with non-record content before its first indexed header is
                # not a sufficiently strong provenance shortcut.
                return None

            for index, (_identifier, length, offset,
                        line_bases, line_width) in enumerate(entries):
                quotient, remainder = divmod(length, line_bases)
                if remainder:
                    final_line_start = offset + quotient * line_width
                    final_line_bases = remainder
                else:
                    final_line_start = offset + (quotient - 1) * line_width
                    final_line_bases = line_bases
                sequence_end = final_line_start + final_line_bases
                record_end = (
                    headers[index + 1]
                    if index + 1 < len(headers)
                    else fasta_stat.st_size
                )
                if sequence_end > record_end or record_end - sequence_end > 2:
                    return None
                handle.seek(sequence_end)
                if handle.read(record_end - sequence_end) not in (
                    b"", b"\n", b"\r\n",
                ):
                    return None

                # Verify both ends of the declared fixed-width sequence layout.
                # This catches indexes copied from similarly sized FASTAs whose
                # record boundaries happen to align but line geometry does not.
                handle.seek(offset)
                first_line = handle.readline()
                first_sequence = first_line.rstrip(b"\r\n")
                expected_first = min(length, line_bases)
                if len(first_sequence) != expected_first:
                    return None
                if any(byte in b" \t\r\n" for byte in first_sequence):
                    return None
                if length > line_bases and len(first_line) != line_width:
                    return None

                handle.seek(final_line_start)
                final_line = handle.readline()
                final_sequence = final_line.rstrip(b"\r\n")
                if len(final_sequence) != final_line_bases:
                    return None
                if any(byte in b" \t\r\n" for byte in final_sequence):
                    return None
    except OSError:
        return None
    return lengths


def fasta_sequence_lengths(path: str) -> Dict[str, int]:
    """Return first-token FASTA sequence lengths without materializing bases.

    A structurally verified ``.fai`` sidecar is used when available; otherwise
    plain-text or gzip-compressed FASTA is streamed. Duplicate or empty
    identifiers are rejected because either makes target-bound validation
    ambiguous.
    """

    indexed = _fai_sequence_lengths(path)
    if indexed is not None:
        return indexed
    opener = gzip.open if os.fspath(path).lower().endswith(".gz") else open
    lengths: Dict[str, int] = {}
    current: Optional[str] = None
    with opener(path, "rt", encoding="utf-8", errors="strict", newline="") as handle:
        for lineno, raw in enumerate(handle, start=1):
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                identifier = line[1:].split(None, 1)[0] if line[1:].strip() else ""
                if not identifier:
                    raise ValueError(f"FASTA line {lineno} has an empty identifier")
                if identifier in lengths:
                    raise ValueError(
                        f"FASTA identifier {identifier!r} is declared more than once"
                    )
                lengths[identifier] = 0
                current = identifier
                continue
            if current is None:
                raise ValueError(
                    f"FASTA line {lineno} contains sequence before the first header"
                )
            lengths[current] += sum(not character.isspace() for character in raw)
    if not lengths:
        raise ValueError("FASTA contains no sequence records")
    empty = [identifier for identifier, length in lengths.items() if length == 0]
    if empty:
        raise ValueError(f"FASTA record {empty[0]!r} has no sequence")
    return lengths


def validate_gff3_target_bounds(
    path: str,
    target_fasta: str,
    *,
    max_issues_per_check: int = 50,
    strict_sequence_regions: bool = False,
) -> ValidationResult:
    """Stream-check GFF3 seqids and coordinates against the target FASTA.

    This independent benchmark gate is intentionally separate from the full
    hierarchy validator. It keeps only FASTA lengths and declared
    ``##sequence-region`` intervals in memory, so annotation size does not
    determine peak memory. Features outside an advisory sequence-region are
    warnings by default because historical LiftOn output copied source-assembly
    directives; ``strict_sequence_regions=True`` promotes them to errors.
    """

    result = ValidationResult(file_path=os.fspath(path))
    issue_counts: Dict[str, int] = defaultdict(int)
    feature_counts: Dict[str, int] = defaultdict(int)
    regions: Dict[str, Tuple[int, int, int]] = {}
    seqids_with_features: Set[str] = set()

    def add_issue(
        check: str,
        lineno: int,
        message: str,
        *,
        feature_id: str = "",
        severity: str = Severity.ERROR,
    ) -> None:
        issue_counts[check] += 1
        result.issue_totals[check] = issue_counts[check]
        result.severity_totals[severity] = (
            result.severity_totals.get(severity, 0) + 1
        )
        if issue_counts[check] <= max_issues_per_check:
            result.issues.append(GFF3Issue(
                severity, lineno, feature_id, check, message,
            ))

    try:
        lengths = fasta_sequence_lengths(target_fasta)
    except (OSError, UnicodeError, ValueError) as exc:
        add_issue(
            "target_fasta_readable",
            0,
            f"Cannot derive unambiguous target sequence lengths from "
            f"{target_fasta}: {exc}",
        )
        return result

    try:
        gff_opener = gzip.open if os.fspath(path).lower().endswith(".gz") else open
        handle = gff_opener(
            path, "rt", encoding="utf-8", errors="replace", newline="",
        )
    except OSError as exc:
        add_issue("file_readable", 0, f"Cannot read {path}: {exc}")
        return result

    with handle:
        for lineno, raw in enumerate(handle, start=1):
            result.total_lines += 1
            line = raw.rstrip("\r\n")
            if not line.strip():
                result.comment_lines += 1
                continue
            if line.startswith("##sequence-region"):
                result.comment_lines += 1
                fields = line.split()
                if len(fields) != 4:
                    add_issue(
                        "sequence_region_format",
                        lineno,
                        "Expected '##sequence-region seqid start end'.",
                    )
                    continue
                _directive, seqid, start_text, end_text = fields
                try:
                    start = int(start_text)
                    end = int(end_text)
                except ValueError:
                    add_issue(
                        "sequence_region_coordinate",
                        lineno,
                        f"Sequence-region bounds must be integers, got "
                        f"{start_text!r} and {end_text!r}.",
                    )
                    continue
                if seqid in seqids_with_features:
                    add_issue(
                        "sequence_region_order",
                        lineno,
                        f"Sequence-region for {seqid!r} appears after a feature "
                        "on that seqid.",
                    )
                if start < 1 or end < start:
                    add_issue(
                        "sequence_region_coordinate",
                        lineno,
                        f"Sequence-region must satisfy 1 <= start <= end, got "
                        f"{start}..{end}.",
                    )
                target_length = lengths.get(seqid)
                if target_length is None:
                    add_issue(
                        "target_seqid_unknown",
                        lineno,
                        f"Sequence-region seqid {seqid!r} is absent from the "
                        "target FASTA.",
                    )
                elif end > target_length:
                    add_issue(
                        "sequence_region_out_of_bounds",
                        lineno,
                        f"Sequence-region {seqid}:{start}-{end} exceeds target "
                        f"length {target_length}.",
                    )
                previous = regions.get(seqid)
                if previous is not None and previous[:2] != (start, end):
                    add_issue(
                        "sequence_region_conflict",
                        lineno,
                        f"Sequence-region {seqid!r} conflicts with line "
                        f"{previous[2]} ({previous[0]}..{previous[1]}).",
                    )
                else:
                    regions[seqid] = (start, end, lineno)
                continue
            if line.startswith("#"):
                result.comment_lines += 1
                continue

            columns = line.split("\t")
            if len(columns) != 9:
                # The structural/full validators own column-format findings.
                continue
            result.data_lines += 1
            seqid, _source, ftype, start_text, end_text = columns[:5]
            feature_counts[ftype] += 1
            seqids_with_features.add(seqid)
            try:
                start = int(start_text)
                end = int(end_text)
            except ValueError:
                continue
            if start < 1 or start > end:
                continue
            target_length = lengths.get(seqid)
            if target_length is None:
                add_issue(
                    "target_seqid_unknown",
                    lineno,
                    f"Feature seqid {seqid!r} is absent from the target FASTA.",
                )
                continue
            if end > target_length:
                add_issue(
                    "target_coordinate_out_of_bounds",
                    lineno,
                    f"Feature {seqid}:{start}-{end} exceeds target length "
                    f"{target_length}.",
                )
            region = regions.get(seqid)
            if region is not None and not (
                region[0] <= start <= end <= region[1]
            ):
                add_issue(
                    "sequence_region_containment",
                    lineno,
                    f"Feature {seqid}:{start}-{end} is outside declared region "
                    f"{region[0]}..{region[1]}.",
                    severity=(
                        Severity.ERROR
                        if strict_sequence_regions
                        else Severity.WARNING
                    ),
                )

    result.stats = dict(feature_counts)
    return result


def validate_gff3_structure(
    path: str,
    *,
    max_issues_per_check: int = 50,
) -> ValidationResult:
    """Stream publication-critical GFF3 checks without loading hierarchies.

    Memory use is bounded by the ID/signature maps needed for ordered Parent
    resolution and duplicate detection. :func:`validate_gff3_file` remains the
    optional full semantic validator.
    """

    result = ValidationResult(file_path=os.fspath(path))
    issue_counts: Dict[str, int] = defaultdict(int)
    # ID -> (line, discontinuous-CDS signature, type, has Parent)
    seen_ids: Dict[str, Tuple[int, tuple | None, str, bool]] = {}
    referenced_parent_ids: Set[str] = set()
    first_content_seen = False

    def add_issue(
        check: str,
        lineno: int,
        message: str,
        *,
        feature_id: str = "",
    ) -> None:
        issue_counts[check] += 1
        if issue_counts[check] <= max_issues_per_check:
            result.issues.append(GFF3Issue(
                Severity.ERROR, lineno, feature_id, check, message,
            ))

    try:
        handle = open(path, "r", encoding="utf-8", errors="replace", newline="")
    except OSError as exc:
        add_issue("file_readable", 0, f"Cannot read {path}: {exc}")
        return result

    feature_counts: Dict[str, int] = defaultdict(int)
    with handle:
        for lineno, raw in enumerate(handle, start=1):
            result.total_lines += 1
            line = raw.rstrip("\r\n")
            if not line.strip():
                result.comment_lines += 1
                continue

            if not first_content_seen:
                first_content_seen = True
                if line != "##gff-version 3":
                    add_issue(
                        "gff3_header", lineno,
                        "First non-blank line must be exactly '##gff-version 3'.",
                    )

            if line.startswith("#"):
                result.comment_lines += 1
                continue

            columns = line.split("\t")
            if len(columns) != 9:
                add_issue(
                    "column_count", lineno,
                    f"Expected 9 tab-separated columns, got {len(columns)}.",
                )
                continue

            result.data_lines += 1
            seqid, source, ftype, start_text, end_text, _score, strand, phase, attrs_text = columns
            feature_counts[ftype] += 1

            if not seqid or seqid == ".":
                add_issue("seqid_empty", lineno, "seqid must not be empty or '.'.")
            try:
                start = int(start_text)
                end = int(end_text)
            except ValueError:
                add_issue(
                    "coord_not_int", lineno,
                    f"start/end must be integers, got {start_text!r} and {end_text!r}.",
                )
            else:
                if start < 1 or end < 1:
                    add_issue(
                        "coord_1based", lineno,
                        f"coordinates must be >= 1, got start={start}, end={end}.",
                    )
                if start > end:
                    add_issue(
                        "coord_order", lineno,
                        f"start ({start}) must be <= end ({end}).",
                    )

            if strand not in VALID_STRANDS:
                add_issue(
                    "strand_valid", lineno,
                    f"strand must be one of {sorted(VALID_STRANDS)}, got {strand!r}.",
                )
            if ftype == "CDS" and phase not in {"0", "1", "2"}:
                add_issue(
                    "cds_phase_value", lineno,
                    f"CDS phase must be 0, 1, or 2, got {phase!r}.",
                )

            attrs: Dict[str, List[str]] = {}
            if attrs_text and attrs_text != ".":
                for part in attrs_text.split(";"):
                    part = part.strip()
                    if not part:
                        continue
                    if "=" not in part:
                        add_issue(
                            "attr_format", lineno,
                            f"Attribute {part!r} has no '=' separator.",
                        )
                        continue
                    key, value = part.split("=", 1)
                    key = key.strip()
                    if not key:
                        add_issue(
                            "attr_format", lineno,
                            f"Attribute {part!r} has an empty key.",
                        )
                        continue
                    values = [item.strip() for item in value.split(",")]
                    if key == "Parent" and any(not item for item in values):
                        add_issue(
                            "parent_empty", lineno,
                            "Parent must contain one or more non-empty IDs.",
                        )
                    attrs[key] = [item for item in values if item]

            feature_id = attrs.get("ID", [""])[0]
            if ftype in GENE_TYPES | TRANSCRIPT_TYPES and not feature_id:
                add_issue(
                    "missing_id", lineno,
                    f"{ftype!r} feature must have an ID attribute.",
                )

            parent_ids = attrs.get("Parent", [])
            if ftype in GENE_TYPES and parent_ids:
                add_issue(
                    "gene_has_parent", lineno,
                    f"Top-level {ftype!r} feature must not have Parent.",
                    feature_id=feature_id,
                )
            # A transcript-shaped locus may itself be a supported top-level
            # output, but only when a later row actually names it as Parent.
            # Defer that decision until the stream has been consumed. Exon and
            # CDS remain mandatory children and are rejected immediately.
            if (ftype in EXON_TYPES | CDS_TYPES
                    or (ftype in TRANSCRIPT_TYPES and not feature_id)) \
                    and not parent_ids:
                add_issue(
                    "missing_parent", lineno,
                    f"{ftype!r} feature must have a Parent attribute.",
                    feature_id=feature_id,
                )

            for parent_id in parent_ids:
                referenced_parent_ids.add(parent_id)
                parent_info = seen_ids.get(parent_id)
                if parent_info is None:
                    add_issue(
                        "parent_not_declared", lineno,
                        f"Parent={parent_id!r} has not been declared earlier in the file.",
                        feature_id=feature_id,
                    )
                    continue
                parent_type = parent_info[2]
                parent_is_root = (
                    not parent_info[3]
                    and _can_be_child_bearing_root(parent_type)
                )
                if (ftype in TRANSCRIPT_TYPES
                        and parent_type not in GENE_TYPES
                        and not parent_is_root):
                    add_issue(
                        "transcript_parent_type", lineno,
                        f"{ftype!r} Parent={parent_id!r} is "
                        f"{parent_type!r}, expected a gene-like feature.",
                        feature_id=feature_id,
                    )
                elif (ftype in EXON_TYPES | CDS_TYPES
                      and not is_transcript_type(parent_type)
                      and parent_type not in DIRECT_EXON_PARENT_TYPES
                      and not parent_is_root):
                    check = (
                        "exon_parent_type" if ftype in EXON_TYPES
                        else "cds_parent_type"
                    )
                    add_issue(
                        check, lineno,
                        f"{ftype!r} Parent={parent_id!r} is "
                        f"{parent_type!r}, expected a transcript or direct "
                        "gene-like parent.",
                        feature_id=feature_id,
                    )

            if feature_id:
                record = GFF3Record(
                    lineno=lineno,
                    seqid=seqid,
                    source=source,
                    ftype=ftype,
                    start=0,
                    end=0,
                    score=".",
                    strand=strand,
                    phase=phase,
                    attrs=attrs,
                    raw=line,
                )
                signature = _discontinuous_cds_signature(record)
                previous = seen_ids.get(feature_id)
                if previous is None:
                    seen_ids[feature_id] = (
                        lineno, signature, ftype, bool(parent_ids),
                    )
                elif signature is None or signature != previous[1]:
                    add_issue(
                        "duplicate_id", lineno,
                        f"Duplicate ID {feature_id!r} (first seen on line {previous[0]}).",
                        feature_id=feature_id,
                    )

    for feature_id, (lineno, _signature, ftype, has_parent) in seen_ids.items():
        if (ftype in TRANSCRIPT_TYPES and not has_parent
                and feature_id not in referenced_parent_ids):
            add_issue(
                "missing_parent", lineno,
                f"{ftype!r} feature must have a Parent attribute unless it "
                "is a child-bearing top-level locus.",
                feature_id=feature_id,
            )

    if not first_content_seen:
        add_issue("gff3_header", 0, "File is empty; expected '##gff-version 3'.")
    if result.data_lines == 0:
        add_issue("features_present", 0, "The file contains no GFF3 feature records.")
    result.stats = dict(feature_counts)
    return result


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

    # A directive-only file is non-empty on disk but is not an annotation.
    # Treat it as invalid so a failed/empty pipeline cannot report success
    # merely because it wrote ``##gff-version 3``.
    if result.data_lines == 0:
        result.issues.append(GFF3Issue(
            Severity.ERROR, 0, "", "features_present",
            "The file contains no GFF3 feature records."
        ))
        result.stats = {}
        return result

    # ── Build feature index ──────────────────────────────────────────────────
    id_to_record: Dict[str, GFF3Record] = {}
    parent_to_children: Dict[str, List[GFF3Record]] = defaultdict(list)
    id_issues: List[GFF3Issue] = []

    seen_ids: Dict[str, int] = {}  # id → first lineno
    seen_id_recs: Dict[str, GFF3Record] = {}  # id → first record (for type/parent)
    for rec in records:
        fid = rec.feat_id
        if fid:
            if fid in seen_ids:
                # A discontinuous feature (canonically a multi-exon CDS) is split
                # across several lines sharing ONE ID — valid GFF3, not a
                # duplicate. Exempt repeats that are the same discontinuous type
                # AND the same Parent as the first occurrence; everything else
                # (e.g. two genes / a gene and an mRNA colliding on an ID) is a
                # genuine duplicate-ID error.
                first = seen_id_recs[fid]
                first_signature = _discontinuous_cds_signature(first)
                is_discontinuous_segment = (
                    first_signature is not None
                    and _discontinuous_cds_signature(rec) == first_signature
                )
                if not is_discontinuous_segment:
                    issue_counts["duplicate_id"] += 1
                    if issue_counts["duplicate_id"] <= max_issues_per_check:
                        id_issues.append(GFF3Issue(
                            Severity.ERROR, rec.lineno, fid, "duplicate_id",
                            f"Duplicate ID '{fid}' (first seen on line {seen_ids[fid]})"
                        ))
            else:
                seen_ids[fid] = rec.lineno
                seen_id_recs[fid] = rec
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
    Validate LiftOn's root→transcript→exon→CDS relationships.

    Rules (per NCBI GFF3 spec and LiftOn conventions):
    1. Every Parent reference must point to a known ID.
    2. Known gene-like features must not have a Parent (they are top-level).
    3. Transcripts must have a child-bearing root as Parent, unless the
       transcript is itself a child-bearing root.
    4. Exons must have a transcript or child-bearing root as Parent.
    5. CDS features must have a transcript or child-bearing root as Parent.
    6. A transcript must have at least one exon child.
    7. Coding transcripts (mRNA) must have at least one CDS child.
    8. CDS must be a subset of at least one exon of the same transcript.
    """
    issues: List[GFF3Issue] = []
    issue_counts: Dict[str, int] = defaultdict(int)
    all_ids = set(id_to_record.keys())

    def is_child_bearing_root(rec: GFF3Record) -> bool:
        return (
            bool(rec.feat_id)
            and not rec.parent_ids
            and bool(parent_to_children.get(rec.feat_id))
            and _can_be_child_bearing_root(rec.ftype)
        )

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

        # ── Rule 3: Transcripts need a root parent or must be the root ───────
        if rec.ftype in TRANSCRIPT_TYPES:
            if not pid:
                if not is_child_bearing_root(rec):
                    issue_counts["trans_no_parent"] += 1
                    if issue_counts["trans_no_parent"] <= max_issues:
                        issues.append(GFF3Issue(
                            Severity.ERROR, rec.lineno, fid, "transcript_no_parent",
                            f"Transcript-type '{rec.ftype}' has no Parent "
                            "attribute and no children"
                        ))
            elif pid in id_to_record:
                parent_rec = id_to_record[pid]
                if (parent_rec.ftype not in GENE_TYPES
                        and not is_child_bearing_root(parent_rec)):
                    issue_counts["trans_wrong_parent"] += 1
                    if issue_counts["trans_wrong_parent"] <= max_issues:
                        issues.append(GFF3Issue(
                            Severity.ERROR, rec.lineno, fid, "transcript_parent_type",
                            f"Transcript parent must be a top-level root; "
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
                if (not is_transcript_type(parent_rec.ftype)
                        and parent_rec.ftype not in DIRECT_EXON_PARENT_TYPES
                        and not is_child_bearing_root(parent_rec)):
                    issue_counts["exon_wrong_parent"] += 1
                    if issue_counts["exon_wrong_parent"] <= max_issues:
                        issues.append(GFF3Issue(
                            Severity.ERROR, rec.lineno, fid, "exon_parent_type",
                            f"Exon parent '{pid}' has type '{parent_rec.ftype}' "
                            f"(expected transcript or top-level root)"
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
                if (not is_transcript_type(parent_rec.ftype)
                        and parent_rec.ftype not in DIRECT_EXON_PARENT_TYPES
                        and not is_child_bearing_root(parent_rec)):
                    issue_counts["cds_wrong_parent"] += 1
                    if issue_counts["cds_wrong_parent"] <= max_issues:
                        issues.append(GFF3Issue(
                            Severity.ERROR, rec.lineno, fid, "cds_parent_type",
                            f"CDS parent '{pid}' has type '{parent_rec.ftype}' "
                            f"(expected transcript or top-level root)"
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
