"""Deterministic synthetic target assemblies and exact-coordinate truth.

The generators implement the two canonical chr22 correctness cases. They are
also parameterised so small fixtures can exercise every transformation in
unit/property tests.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import os
import urllib.parse
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping, Sequence

from lifton.gff3_validator import GENE_TYPES, TRANSCRIPT_TYPES


SCHEMA_VERSION = 1
FRAGMENT_CUTS = (
    17_290_898,
    22_364_968,
    27_189_473,
    34_494_491,
    41_032_950,
    45_086_965,
    47_832_626,
    49_262_605,
)
SV_DELETION = (18_177_398, 18_516_337)
SV_INVERSION = (36_140_323, 38_318_084)
SV_DUPLICATION = (41_092_592, 44_331_714)
SV_INSERT_AFTER = 45_086_965
SV_INSERT_LENGTH = 100_000
SV_SEED = "lifton-canonical-v2-20260718"
_COMPLEMENT = str.maketrans(
    "ACGTRYMKBDHVNacgtrymkbdhvn",
    "TGCAYRKMVHDBNtgcayrkmvhdbn",
)


@dataclass(frozen=True)
class CoordinateSegment:
    source_start: int
    source_end: int
    target_seqid: str
    target_start: int
    target_end: int
    orientation: str = "+"
    copy_index: int | None = None

    def map_interval(self, start: int, end: int) -> tuple[int, int]:
        if not self.source_start <= start <= end <= self.source_end:
            raise ValueError(
                f"{start}..{end} is outside source segment "
                f"{self.source_start}..{self.source_end}"
            )
        if self.orientation == "+":
            mapped_start = self.target_start + start - self.source_start
            return mapped_start, mapped_start + end - start
        mapped_start = self.target_start + self.source_end - end
        return mapped_start, mapped_start + end - start


@dataclass(frozen=True)
class GffRow:
    seqid: str
    source: str
    feature_type: str
    start: int
    end: int
    score: str
    strand: str
    phase: str
    attributes: tuple[tuple[str, str], ...]
    lineno: int

    def values(self, key: str) -> tuple[str, ...]:
        for name, value in self.attributes:
            if name == key:
                return tuple(
                    urllib.parse.unquote(item, errors="strict")
                    for item in value.split(",")
                    if item
                )
        return ()

    @property
    def identifier(self) -> str | None:
        values = self.values("ID")
        return values[0] if values else None

    @property
    def parents(self) -> tuple[str, ...]:
        return self.values("Parent")


def _open_text(path: str | Path):
    path = Path(path)
    opener = gzip.open if path.name.lower().endswith(".gz") else open
    return opener(path, "rt", encoding="utf-8", errors="strict", newline="")


def _sha256_file(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_fasta_sequence(path: str | Path, seqid: str) -> str:
    """Read one first-token FASTA record and reject duplicates."""

    chunks: list[str] = []
    selected = False
    found = False
    with _open_text(path) as handle:
        for lineno, raw in enumerate(handle, start=1):
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                identifier = line[1:].split(None, 1)[0] if line[1:].strip() else ""
                if not identifier:
                    raise ValueError(f"{path}: empty FASTA ID on line {lineno}")
                selected = identifier == seqid
                if selected:
                    if found:
                        raise ValueError(f"{path}: duplicate FASTA ID {seqid!r}")
                    found = True
                continue
            if selected:
                chunks.append("".join(line.split()).upper())
    if not found:
        raise ValueError(f"{path}: FASTA ID {seqid!r} was not found")
    sequence = "".join(chunks)
    if not sequence:
        raise ValueError(f"{path}: FASTA ID {seqid!r} has no sequence")
    return sequence


def read_gff_rows(path: str | Path, seqid: str) -> list[GffRow]:
    rows = []
    with _open_text(path) as handle:
        for lineno, raw in enumerate(handle, start=1):
            if not raw.strip() or raw.startswith("#"):
                continue
            columns = raw.rstrip("\r\n").split("\t")
            if len(columns) != 9:
                raise ValueError(
                    f"{path}: line {lineno} has {len(columns)} columns, expected 9"
                )
            if columns[0] != seqid or columns[2] == "region":
                continue
            try:
                start, end = int(columns[3]), int(columns[4])
            except ValueError as exc:
                raise ValueError(
                    f"{path}: line {lineno} has non-integer coordinates"
                ) from exc
            if start < 1 or end < start:
                raise ValueError(
                    f"{path}: line {lineno} has invalid coordinates "
                    f"{start}..{end}"
                )
            attributes = []
            if columns[8] and columns[8] != ".":
                for item in columns[8].split(";"):
                    if not item:
                        continue
                    key, separator, value = item.partition("=")
                    if not separator:
                        raise ValueError(
                            f"{path}: line {lineno} has malformed attribute "
                            f"{item!r}"
                        )
                    attributes.append((key, value))
            rows.append(GffRow(
                seqid=columns[0],
                source=columns[1],
                feature_type=columns[2],
                start=start,
                end=end,
                score=columns[5],
                strand=columns[6],
                phase=columns[7],
                attributes=tuple(attributes),
                lineno=lineno,
            ))
    if not rows:
        raise ValueError(f"{path}: no feature rows found for {seqid!r}")
    return rows


def deterministic_dna(length: int, seed: str) -> str:
    """Generate platform-independent deterministic A/C/G/T sequence."""

    if length < 0:
        raise ValueError("deterministic sequence length must be non-negative")
    alphabet = "ACGT"
    sequence = []
    counter = 0
    while len(sequence) < length:
        digest = hashlib.sha256(
            seed.encode("utf-8") + b"\0" + counter.to_bytes(8, "big")
        ).digest()
        for byte in digest:
            sequence.extend((
                alphabet[(byte >> 6) & 3],
                alphabet[(byte >> 4) & 3],
                alphabet[(byte >> 2) & 3],
                alphabet[byte & 3],
            ))
            if len(sequence) >= length:
                break
        counter += 1
    return "".join(sequence[:length])


def _reverse_complement(sequence: str) -> str:
    return sequence.translate(_COMPLEMENT)[::-1]


def _write_fasta(path: Path, records: Sequence[tuple[str, str]]) -> None:
    with path.open("w", encoding="ascii", newline="\n") as handle:
        for seqid, sequence in records:
            handle.write(f">{seqid}\n")
            for offset in range(0, len(sequence), 60):
                handle.write(sequence[offset:offset + 60] + "\n")
        handle.flush()
        os.fsync(handle.fileno())


def _id_parent_index(
    rows: Sequence[GffRow],
) -> tuple[dict[str, tuple[str, ...]], dict[str, GffRow]]:
    parents: dict[str, tuple[str, ...]] = {}
    model_rows: dict[str, GffRow] = {}
    for row in rows:
        if row.identifier and row.identifier not in parents:
            parents[row.identifier] = row.parents
        if (
            row.identifier
            and row.feature_type in GENE_TYPES | TRANSCRIPT_TYPES
        ):
            if row.identifier in model_rows:
                raise ValueError(f"duplicate model ID {row.identifier!r}")
            model_rows[row.identifier] = row
    return parents, model_rows


def _roots_by_identifier(
    parents: Mapping[str, tuple[str, ...]],
) -> dict[str, tuple[str, ...]]:
    memo: dict[str, tuple[str, ...]] = {}

    def roots(identifier: str, stack: frozenset[str] = frozenset()) -> tuple[str, ...]:
        if identifier in memo:
            return memo[identifier]
        if identifier in stack:
            raise ValueError(f"cyclic GFF3 Parent hierarchy at {identifier!r}")
        parent_ids = parents.get(identifier, ())
        if not parent_ids:
            result = (identifier,)
        else:
            result = tuple(sorted({
                root
                for parent in parent_ids
                for root in roots(parent, stack | {identifier})
            }))
        memo[identifier] = result
        return result

    for identifier in parents:
        roots(identifier)
    return memo


def _row_roots(
    row: GffRow,
    roots_by_id: Mapping[str, tuple[str, ...]],
) -> tuple[str, ...]:
    identifiers = row.parents or ((row.identifier,) if row.identifier else ())
    return tuple(sorted({
        root
        for identifier in identifiers
        for root in roots_by_id.get(identifier, (identifier,))
    }))


def _crosses(start: int, end: int, boundary_after: int) -> bool:
    return start <= boundary_after < end


def _fragment_status(
    row: GffRow,
    cuts: Sequence[int],
) -> str:
    return (
        "breakpoint_crossing"
        if any(_crosses(row.start, row.end, cut) for cut in cuts)
        else "retained"
    )


def _sv_status(
    row: GffRow,
    *,
    deletion: tuple[int, int],
    inversion: tuple[int, int],
    duplication: tuple[int, int],
    insert_after: int,
) -> str:
    deletion_start, deletion_end = deletion
    if deletion_start <= row.start <= row.end <= deletion_end:
        return "deleted"
    if not (row.end < deletion_start or row.start > deletion_end):
        return "breakpoint_crossing"
    boundaries = (
        inversion[0] - 1,
        inversion[1],
        duplication[0] - 1,
        duplication[1],
        insert_after,
    )
    if any(_crosses(row.start, row.end, boundary) for boundary in boundaries):
        return "breakpoint_crossing"
    if duplication[0] <= row.start <= row.end <= duplication[1]:
        return "duplicated"
    return "retained"


def _root_statuses(
    rows: Sequence[GffRow],
    roots_by_id: Mapping[str, tuple[str, ...]],
    classifier,
) -> dict[str, str]:
    top_rows = {
        row.identifier: row
        for row in rows
        if row.identifier
        and not row.parents
    }
    roots = sorted({
        root
        for row in rows
        for root in _row_roots(row, roots_by_id)
    })
    statuses = {}
    for root in roots:
        row = top_rows.get(root)
        if row is None:
            members = [
                item for item in rows
                if root in _row_roots(item, roots_by_id)
            ]
            start = min(item.start for item in members)
            end = max(item.end for item in members)
            row = GffRow(
                members[0].seqid,
                members[0].source,
                "synthetic_root",
                start,
                end,
                ".",
                members[0].strand,
                ".",
                (("ID", root),),
                0,
            )
        statuses[root] = classifier(row)
    return statuses


def _fragment_records(
    sequence: str,
    source_seqid: str,
    cuts: Sequence[int],
) -> tuple[list[tuple[str, str]], list[CoordinateSegment]]:
    if tuple(sorted(set(cuts))) != tuple(cuts):
        raise ValueError("fragment cuts must be unique and strictly increasing")
    if any(cut < 1 or cut >= len(sequence) for cut in cuts):
        raise ValueError(
            f"fragment cuts must be within 1..{len(sequence) - 1}"
        )
    records = []
    segments = []
    source_start = 1
    for index, source_end in enumerate((*cuts, len(sequence)), start=1):
        target_seqid = f"{source_seqid}_frag_{index:03d}"
        fragment = sequence[source_start - 1:source_end]
        records.append((target_seqid, fragment))
        segments.append(CoordinateSegment(
            source_start,
            source_end,
            target_seqid,
            1,
            len(fragment),
        ))
        source_start = source_end + 1
    return records, segments


def _sv_record(
    sequence: str,
    source_seqid: str,
    *,
    deletion: tuple[int, int],
    inversion: tuple[int, int],
    duplication: tuple[int, int],
    insert_after: int,
    insertion: str,
) -> tuple[list[tuple[str, str]], list[CoordinateSegment]]:
    length = len(sequence)
    coordinates = (
        deletion[0], deletion[1], inversion[0], inversion[1],
        duplication[0], duplication[1], insert_after,
    )
    if not (
        1 <= deletion[0] <= deletion[1] < inversion[0] <= inversion[1]
        < duplication[0] <= duplication[1] <= insert_after < length
    ):
        raise ValueError(
            f"SV operations must be ordered, non-overlapping, and within "
            f"the source length {length}; got {coordinates}"
        )

    chunks: list[str] = []
    segments: list[CoordinateSegment] = []
    target_cursor = 1

    def append_source(
        source_start: int,
        source_end: int,
        *,
        orientation: str = "+",
        copy_index: int | None = None,
    ) -> None:
        nonlocal target_cursor
        chunk = sequence[source_start - 1:source_end]
        if orientation == "-":
            chunk = _reverse_complement(chunk)
        chunks.append(chunk)
        target_end = target_cursor + len(chunk) - 1
        segments.append(CoordinateSegment(
            source_start,
            source_end,
            source_seqid,
            target_cursor,
            target_end,
            orientation,
            copy_index,
        ))
        target_cursor = target_end + 1

    append_source(1, deletion[0] - 1)
    append_source(deletion[1] + 1, inversion[0] - 1)
    append_source(inversion[0], inversion[1], orientation="-")
    append_source(inversion[1] + 1, duplication[0] - 1)
    append_source(duplication[0], duplication[1], copy_index=1)
    append_source(duplication[0], duplication[1], copy_index=2)
    append_source(duplication[1] + 1, insert_after)
    chunks.append(insertion)
    target_cursor += len(insertion)
    append_source(insert_after + 1, length)
    return [(source_seqid, "".join(chunks))], segments


def _suffix(identifier: str, copy_index: int | None) -> str:
    return (
        identifier
        if copy_index is None
        else f"{identifier}__synthetic_copy{copy_index}"
    )


def _mapped_attributes(
    attributes: Sequence[tuple[str, str]],
    copy_index: int | None,
) -> str:
    mapped = []
    for key, value in attributes:
        if copy_index is not None and key in {"ID", "Parent", "Derives_from"}:
            value = ",".join(
                _suffix(identifier, copy_index)
                for identifier in value.split(",")
            )
        mapped.append(f"{key}={value}")
    return ";".join(mapped) if mapped else "."


def _find_segments(
    segments: Sequence[CoordinateSegment],
    start: int,
    end: int,
    copy_indices: Sequence[int | None],
) -> list[CoordinateSegment]:
    return [
        segment for segment in segments
        if segment.copy_index in copy_indices
        and segment.source_start <= start <= end <= segment.source_end
    ]


def _transform_rows(
    rows: Sequence[GffRow],
    segments: Sequence[CoordinateSegment],
    roots_by_id: Mapping[str, tuple[str, ...]],
    root_status: Mapping[str, str],
) -> list[str]:
    transformed = []
    for row in rows:
        roots = _row_roots(row, roots_by_id)
        statuses = {root_status[root] for root in roots}
        if not statuses or statuses & {"deleted", "breakpoint_crossing"}:
            continue
        if statuses == {"duplicated"}:
            copy_indices: tuple[int | None, ...] = (1, 2)
        elif statuses == {"retained"}:
            copy_indices = (None,)
        else:
            # A multi-parent row spanning roots with incompatible copy state is
            # not an exact synthetic truth model.
            continue
        mapped_segments = _find_segments(
            segments, row.start, row.end, copy_indices,
        )
        if len(mapped_segments) != len(copy_indices):
            raise ValueError(
                f"GFF3 line {row.lineno} does not map to exactly one target "
                "segment per expected copy"
            )
        for segment in sorted(
            mapped_segments,
            key=lambda item: (
                item.copy_index or 0, item.target_seqid, item.target_start,
            ),
        ):
            start, end = segment.map_interval(row.start, row.end)
            strand = row.strand
            if segment.orientation == "-":
                strand = {"+": "-", "-": "+"}.get(strand, strand)
            transformed.append("\t".join((
                segment.target_seqid,
                "LiftOn-synthetic-truth",
                row.feature_type,
                str(start),
                str(end),
                row.score,
                strand,
                row.phase,
                _mapped_attributes(row.attributes, segment.copy_index),
            )))
    return transformed


def _mapping_document(
    model_rows: Mapping[str, GffRow],
    roots_by_id: Mapping[str, tuple[str, ...]],
    root_status: Mapping[str, str],
) -> dict[str, Any]:
    mappings = []
    for identifier, row in sorted(model_rows.items()):
        roots = _row_roots(row, roots_by_id)
        statuses = {root_status[root] for root in roots}
        if statuses == {"duplicated"}:
            truth_ids = [_suffix(identifier, 1), _suffix(identifier, 2)]
            status = "retained"
        elif statuses == {"retained"}:
            truth_ids = [identifier]
            status = "retained"
        elif statuses == {"deleted"}:
            truth_ids = []
            status = "deleted"
        else:
            truth_ids = []
            status = "breakpoint_crossing"
        mappings.append({
            "source_id": identifier,
            "truth_ids": truth_ids,
            "feature_type": (
                "gene" if row.feature_type in GENE_TYPES else "transcript"
            ),
            "status": status,
        })
    return {
        "schema_version": SCHEMA_VERSION,
        "method": "deterministic-synthetic-coordinate-map-v1",
        "mappings": mappings,
    }


def _segment_document(segment: CoordinateSegment) -> dict[str, Any]:
    return {
        "source_start": segment.source_start,
        "source_end": segment.source_end,
        "target_seqid": segment.target_seqid,
        "target_start": segment.target_start,
        "target_end": segment.target_end,
        "orientation": segment.orientation,
        "copy_index": segment.copy_index,
    }


def _write_json(path: Path, document: Mapping[str, Any]) -> None:
    with path.open("w", encoding="utf-8") as handle:
        json.dump(document, handle, indent=2, sort_keys=True)
        handle.write("\n")
        handle.flush()
        os.fsync(handle.fileno())


def _publish_case(
    output_dir: Path,
    *,
    source_fasta: Path,
    source_gff: Path,
    source_seqid: str,
    records: Sequence[tuple[str, str]],
    segments: Sequence[CoordinateSegment],
    rows: Sequence[GffRow],
    roots_by_id: Mapping[str, tuple[str, ...]],
    model_rows: Mapping[str, GffRow],
    root_status: Mapping[str, str],
    transform: Mapping[str, Any],
    force: bool,
) -> dict[str, Path]:
    output_dir.mkdir(parents=True, exist_ok=True)
    destinations = {
        "target_fasta": output_dir / "target.fa",
        "truth_gff": output_dir / "target.truth.gff3",
        "ortholog_map": output_dir / "ortholog_map.json",
        "manifest": output_dir / "transform.manifest.json",
    }
    if not force:
        existing = [str(path) for path in destinations.values() if path.exists()]
        if existing:
            raise FileExistsError(
                "refusing to replace synthetic artifacts: " + ", ".join(existing)
            )
    temporary = {
        name: path.with_name(path.name + ".tmp")
        for name, path in destinations.items()
    }
    for path in temporary.values():
        if path.exists():
            path.unlink()

    _write_fasta(temporary["target_fasta"], records)
    transformed_rows = _transform_rows(
        rows, segments, roots_by_id, root_status,
    )
    with temporary["truth_gff"].open(
        "w", encoding="utf-8", newline="\n",
    ) as handle:
        handle.write("##gff-version 3\n")
        for seqid, sequence in records:
            handle.write(f"##sequence-region {seqid} 1 {len(sequence)}\n")
        for row in transformed_rows:
            handle.write(row + "\n")
        handle.flush()
        os.fsync(handle.fileno())

    mapping = _mapping_document(model_rows, roots_by_id, root_status)
    _write_json(temporary["ortholog_map"], mapping)
    status_counts = Counter(root_status.values())
    manifest = {
        "schema_version": SCHEMA_VERSION,
        "kind": "lifton-synthetic-target-truth",
        "source": {
            "seqid": source_seqid,
            "sequence_length": sum(
                segment.source_end - segment.source_start + 1
                for segment in segments
                if segment.copy_index in {None, 1}
            ) + (
                transform.get("deleted_length", 0)
            ),
            "fasta": {
                "name": source_fasta.name,
                "sha256": _sha256_file(source_fasta),
            },
            "gff": {
                "name": source_gff.name,
                "sha256": _sha256_file(source_gff),
            },
        },
        "transform": dict(transform),
        "coordinate_segments": [
            _segment_document(segment) for segment in segments
        ],
        "truth_policy": {
            "breakpoint_crossing": "excluded",
            "deleted": "expected_absent",
            "retained": "exact_coordinate_truth",
            "duplicated": "two_exact_expected_copies",
        },
        "root_status_counts": dict(sorted(status_counts.items())),
        "root_status": dict(sorted(root_status.items())),
        "outputs": {
            "target_fasta": {
                "name": destinations["target_fasta"].name,
                "sha256": _sha256_file(temporary["target_fasta"]),
                "sequence_lengths": {
                    seqid: len(sequence) for seqid, sequence in records
                },
            },
            "truth_gff": {
                "name": destinations["truth_gff"].name,
                "sha256": _sha256_file(temporary["truth_gff"]),
                "feature_records": len(transformed_rows),
            },
            "ortholog_map": {
                "name": destinations["ortholog_map"].name,
                "sha256": _sha256_file(temporary["ortholog_map"]),
                "entries": len(mapping["mappings"]),
            },
        },
    }
    _write_json(temporary["manifest"], manifest)
    for name in ("target_fasta", "truth_gff", "ortholog_map", "manifest"):
        os.replace(temporary[name], destinations[name])
    return destinations


def build_fragmented_case(
    source_fasta: str | Path,
    source_gff: str | Path,
    output_dir: str | Path,
    *,
    source_seqid: str = "chr22",
    cuts: Sequence[int] = FRAGMENT_CUTS,
    force: bool = False,
) -> dict[str, Path]:
    """Build the fixed-cut fragmented target and exact mapped annotation."""

    source_fasta = Path(source_fasta)
    source_gff = Path(source_gff)
    sequence = read_fasta_sequence(source_fasta, source_seqid)
    rows = read_gff_rows(source_gff, source_seqid)
    parents, model_rows = _id_parent_index(rows)
    roots_by_id = _roots_by_identifier(parents)
    records, segments = _fragment_records(sequence, source_seqid, tuple(cuts))
    root_status = _root_statuses(
        rows,
        roots_by_id,
        lambda row: _fragment_status(row, cuts),
    )
    return _publish_case(
        Path(output_dir),
        source_fasta=source_fasta,
        source_gff=source_gff,
        source_seqid=source_seqid,
        records=records,
        segments=segments,
        rows=rows,
        roots_by_id=roots_by_id,
        model_rows=model_rows,
        root_status=root_status,
        transform={
            "id": "v2_synth_chr22_fragmented",
            "kind": "fragmentation",
            "cuts_after_source_coordinate": list(cuts),
        },
        force=force,
    )


def build_sv_case(
    source_fasta: str | Path,
    source_gff: str | Path,
    output_dir: str | Path,
    *,
    source_seqid: str = "chr22",
    deletion: tuple[int, int] = SV_DELETION,
    inversion: tuple[int, int] = SV_INVERSION,
    duplication: tuple[int, int] = SV_DUPLICATION,
    insert_after: int = SV_INSERT_AFTER,
    insertion_length: int = SV_INSERT_LENGTH,
    seed: str = SV_SEED,
    force: bool = False,
) -> dict[str, Path]:
    """Build the canonical deletion/inversion/duplication/insertion target."""

    source_fasta = Path(source_fasta)
    source_gff = Path(source_gff)
    sequence = read_fasta_sequence(source_fasta, source_seqid)
    rows = read_gff_rows(source_gff, source_seqid)
    parents, model_rows = _id_parent_index(rows)
    roots_by_id = _roots_by_identifier(parents)
    insertion = deterministic_dna(insertion_length, seed)
    records, segments = _sv_record(
        sequence,
        source_seqid,
        deletion=deletion,
        inversion=inversion,
        duplication=duplication,
        insert_after=insert_after,
        insertion=insertion,
    )
    root_status = _root_statuses(
        rows,
        roots_by_id,
        lambda row: _sv_status(
            row,
            deletion=deletion,
            inversion=inversion,
            duplication=duplication,
            insert_after=insert_after,
        ),
    )
    return _publish_case(
        Path(output_dir),
        source_fasta=source_fasta,
        source_gff=source_gff,
        source_seqid=source_seqid,
        records=records,
        segments=segments,
        rows=rows,
        roots_by_id=roots_by_id,
        model_rows=model_rows,
        root_status=root_status,
        transform={
            "id": "v2_synth_chr22_sv",
            "kind": "structural_variation",
            "coordinate_basis": "1-based-inclusive-source-simultaneous",
            "deletion": list(deletion),
            "deleted_length": deletion[1] - deletion[0] + 1,
            "inversion": list(inversion),
            "tandem_duplication": list(duplication),
            "insert_after": insert_after,
            "insertion_length": insertion_length,
            "insertion_sha256": hashlib.sha256(
                insertion.encode("ascii")
            ).hexdigest(),
            "seed": seed,
        },
        force=force,
    )


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="kind", required=True)
    for kind in ("fragmented", "sv"):
        child = subparsers.add_parser(kind)
        child.add_argument("source_fasta", type=Path)
        child.add_argument("source_gff", type=Path)
        child.add_argument("output_dir", type=Path)
        child.add_argument("--source-seqid", default="chr22")
        child.add_argument("--force", action="store_true")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    arguments = _parser().parse_args(argv)
    builder = (
        build_fragmented_case
        if arguments.kind == "fragmented"
        else build_sv_case
    )
    outputs = builder(
        arguments.source_fasta,
        arguments.source_gff,
        arguments.output_dir,
        source_seqid=arguments.source_seqid,
        force=arguments.force,
    )
    print(json.dumps({
        name: str(path.resolve()) for name, path in outputs.items()
    }, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
