"""Neutral, ortholog-scoped scoring against an independent target annotation.

The target annotation is an evaluation-only input. Callers must run LiftOn
before invoking this module; no function here builds or executes a LiftOn
command. Metrics are based only on target coordinates and an optional,
independently curated source-to-target ortholog map.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import os
import sys
import urllib.parse
from collections import Counter, defaultdict
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

from lifton.gff3_validator import GENE_TYPES, TRANSCRIPT_TYPES


SCHEMA_VERSION = 1
METHOD = "ortholog-scoped-target-coordinate-v1"
FILTER_METHOD = "target-truth-target-seqid-filter-v1"
MAPPING_STATUSES = {
    "retained",
    "deleted",
    "excluded",
    "breakpoint_crossing",
    "unscored",
}


@dataclass(frozen=True)
class Model:
    identifier: str
    seqid: str
    start: int
    end: int
    strand: str
    parents: tuple[str, ...]
    feature_type: str


@dataclass(frozen=True)
class MappingEntry:
    source_id: str
    truth_ids: tuple[str, ...]
    feature_type: str | None = None
    status: str = "retained"

    @property
    def scorable(self) -> bool:
        return self.status in {"retained", "deleted"}


@dataclass
class Annotation:
    genes: dict[str, Model]
    transcripts: dict[str, Model]
    exons: dict[str, list[tuple[str, int, int, str]]]
    cds: dict[str, list[tuple[str, int, int, str]]]
    transcripts_by_gene: dict[str, list[str]]


@dataclass(frozen=True)
class _Group:
    key: str
    prediction_ids: tuple[str, ...]
    truth_ids: tuple[str, ...]


def _parse_attributes(text: str) -> dict[str, tuple[str, ...]]:
    attributes: dict[str, tuple[str, ...]] = {}
    if not text or text == ".":
        return attributes
    for item in text.split(";"):
        item = item.strip()
        if not item or "=" not in item:
            continue
        key, value = item.split("=", 1)
        attributes[key.strip()] = tuple(
            urllib.parse.unquote(part.strip(), errors="strict")
            for part in value.split(",")
            if part.strip()
        )
    return attributes


def parse_annotation(path: str | Path) -> Annotation:
    """Parse only the hierarchy and intervals needed by target-truth metrics."""

    genes: dict[str, Model] = {}
    transcripts: dict[str, Model] = {}
    exons: dict[str, list[tuple[str, int, int, str]]] = defaultdict(list)
    cds: dict[str, list[tuple[str, int, int, str]]] = defaultdict(list)
    transcripts_by_gene: dict[str, list[str]] = defaultdict(list)

    path = Path(path)
    opener = gzip.open if path.name.lower().endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8", errors="strict") as handle:
        for lineno, raw in enumerate(handle, start=1):
            if not raw.strip() or raw.startswith("#"):
                continue
            columns = raw.rstrip("\r\n").split("\t")
            if len(columns) != 9:
                raise ValueError(
                    f"{path}: line {lineno} has {len(columns)} columns, expected 9"
                )
            seqid, _source, feature_type, start_text, end_text, \
                _score, strand, _phase, attribute_text = columns
            try:
                start = int(start_text)
                end = int(end_text)
            except ValueError as exc:
                raise ValueError(
                    f"{path}: line {lineno} has non-integer coordinates"
                ) from exc
            if start < 1 or end < start:
                raise ValueError(
                    f"{path}: line {lineno} has invalid coordinates "
                    f"{start_text}..{end_text}"
                )
            attributes = _parse_attributes(attribute_text)
            feature_ids = attributes.get("ID", ())
            parents = attributes.get("Parent", ())
            if feature_type in GENE_TYPES | TRANSCRIPT_TYPES:
                if len(feature_ids) != 1:
                    raise ValueError(
                        f"{path}: line {lineno} {feature_type} must have one ID"
                    )
                model = Model(
                    feature_ids[0],
                    seqid,
                    start,
                    end,
                    strand,
                    parents,
                    feature_type,
                )
                destination = (
                    genes if feature_type in GENE_TYPES else transcripts
                )
                if model.identifier in destination:
                    raise ValueError(
                        f"{path}: duplicate {feature_type} ID "
                        f"{model.identifier!r}"
                    )
                destination[model.identifier] = model
                if feature_type in TRANSCRIPT_TYPES:
                    for parent in parents:
                        transcripts_by_gene[parent].append(model.identifier)
            elif feature_type in {"exon", "CDS"}:
                interval = (seqid, start, end, strand)
                destination = exons if feature_type == "exon" else cds
                for parent in parents:
                    destination[parent].append(interval)

    for values in (*exons.values(), *cds.values()):
        values.sort(key=lambda interval: (
            interval[0], interval[1], interval[2], interval[3],
        ))
    for values in transcripts_by_gene.values():
        values.sort()
    return Annotation(
        genes=genes,
        transcripts=transcripts,
        exons=dict(exons),
        cds=dict(cds),
        transcripts_by_gene=dict(transcripts_by_gene),
    )


def _fasta_sequence_lengths(path: str | Path) -> dict[str, int]:
    path = Path(path)
    opener = gzip.open if path.name.lower().endswith(".gz") else open
    lengths: dict[str, int] = {}
    identifier = None
    length = 0
    with opener(path, "rt", encoding="utf-8", errors="strict") as handle:
        for line_number, raw in enumerate(handle, start=1):
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if identifier is not None:
                    lengths[identifier] = length
                header = line[1:].strip()
                identifier = header.split()[0] if header else ""
                if not identifier:
                    raise ValueError(
                        f"{path}: FASTA line {line_number} has an empty ID"
                    )
                if identifier in lengths:
                    raise ValueError(
                        f"{path}: duplicate FASTA sequence ID {identifier!r}"
                    )
                length = 0
            elif identifier is None:
                raise ValueError(
                    f"{path}: FASTA sequence appears before the first header"
                )
            else:
                length += len(line)
    if identifier is not None:
        if identifier in lengths:
            raise ValueError(
                f"{path}: duplicate FASTA sequence ID {identifier!r}"
            )
        lengths[identifier] = length
    if not lengths or any(length <= 0 for length in lengths.values()):
        raise ValueError(f"{path}: FASTA has missing or empty sequences")
    return lengths


def filter_truth_to_target_fasta(
    truth_gff: str | Path,
    target_fasta: str | Path,
    output: str | Path,
    *,
    force: bool = False,
) -> dict[str, Any]:
    """Filter target truth to exactly the sequence IDs in a target FASTA.

    The operation is streaming and deterministic. It is intended for creating
    immutable subset-panel truth before a campaign, never during a timed cell.
    """

    truth_path = Path(truth_gff).expanduser().resolve()
    fasta_path = Path(target_fasta).expanduser().resolve()
    output_path = Path(output).expanduser().resolve()
    if output_path in {truth_path, fasta_path}:
        raise ValueError("filtered target truth output must be a new path")
    if output_path.exists() and not force:
        raise FileExistsError(
            f"refusing to replace filtered target truth: {output_path}"
        )
    sequence_lengths = _fasta_sequence_lengths(fasta_path)
    allowed = set(sequence_lengths)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    temporary = output_path.with_name(
        f".{output_path.name}.tmp-{os.getpid()}"
    )
    if temporary.exists() or temporary.is_symlink():
        raise FileExistsError(f"temporary output already exists: {temporary}")
    opener = gzip.open if truth_path.name.lower().endswith(".gz") else open
    kept_records = 0
    dropped_records = 0
    kept_models = 0
    kept_seqids: set[str] = set()
    try:
        with opener(
            truth_path, "rt", encoding="utf-8", errors="strict",
        ) as source, temporary.open(
            "x", encoding="utf-8", newline="\n",
        ) as destination:
            destination.write("##gff-version 3\n")
            for seqid, length in sequence_lengths.items():
                destination.write(
                    f"##sequence-region {seqid} 1 {length}\n"
                )
            in_embedded_fasta = False
            for line_number, raw in enumerate(source, start=1):
                if raw.startswith("##FASTA"):
                    in_embedded_fasta = True
                    continue
                if in_embedded_fasta or not raw.strip() or raw.startswith("#"):
                    continue
                columns = raw.rstrip("\r\n").split("\t")
                if len(columns) != 9:
                    raise ValueError(
                        f"{truth_path}: line {line_number} has "
                        f"{len(columns)} columns, expected 9"
                    )
                if columns[8] not in {"", "."} and "=" not in columns[8]:
                    raise ValueError(
                        f"{truth_path}: line {line_number} is not GFF3"
                    )
                if columns[0] not in allowed:
                    dropped_records += 1
                    continue
                destination.write("\t".join(columns) + "\n")
                kept_records += 1
                kept_seqids.add(columns[0])
                kept_models += columns[2] in GENE_TYPES | TRANSCRIPT_TYPES
            if kept_records == 0 or kept_models == 0:
                raise ValueError(
                    "target-seqid filtering retained no gene/transcript truth"
                )
            destination.flush()
            os.fsync(destination.fileno())
        os.replace(temporary, output_path)
    except Exception:
        if temporary.exists() or temporary.is_symlink():
            temporary.unlink()
        raise
    return {
        "schema_version": SCHEMA_VERSION,
        "method": FILTER_METHOD,
        "inputs": {
            "truth_gff": {
                "path": str(truth_path),
                "size": truth_path.stat().st_size,
                "sha256": _sha256_file(truth_path),
            },
            "target_fasta": {
                "path": str(fasta_path),
                "size": fasta_path.stat().st_size,
                "sha256": _sha256_file(fasta_path),
                "sequence_lengths": dict(sequence_lengths),
            },
        },
        "output": {
            "path": str(output_path),
            "size": output_path.stat().st_size,
            "sha256": _sha256_file(output_path),
            "feature_records": kept_records,
            "dropped_feature_records": dropped_records,
            "feature_seqids": sorted(kept_seqids),
        },
    }


def _normalise_mapping_item(item: Mapping[str, Any]) -> MappingEntry:
    if not isinstance(item, Mapping):
        raise ValueError("ortholog mapping entry is not an object")
    source_id = str(item.get("source_id", "")).strip()
    if not source_id:
        raise ValueError("ortholog mapping entry has no source_id")
    raw_truth_ids = item.get("truth_ids", item.get("truth_id", ()))
    if isinstance(raw_truth_ids, str):
        truth_ids = tuple(
            value.strip() for value in raw_truth_ids.split(",") if value.strip()
        )
    elif raw_truth_ids is None:
        truth_ids = ()
    elif isinstance(raw_truth_ids, Sequence):
        truth_ids = tuple(
            str(value).strip() for value in raw_truth_ids if str(value).strip()
        )
    else:
        raise ValueError(
            f"ortholog mapping {source_id!r} truth_ids must be a list or string"
        )
    if len(set(truth_ids)) != len(truth_ids):
        raise ValueError(
            f"ortholog mapping {source_id!r} repeats a target truth ID"
        )
    feature_type = item.get("feature_type")
    return MappingEntry(
        source_id=source_id,
        truth_ids=truth_ids,
        feature_type=(
            str(feature_type).strip().lower() if feature_type else None
        ),
        status=str(item.get("status", "retained")).strip().lower(),
    )


def _sha256_file(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _load_ortholog_map(
    path: str | Path | None,
) -> tuple[list[MappingEntry], dict[str, Any], Any]:
    """Load JSON or headered TSV ortholog mappings.

    JSON accepts ``{"mappings": [...]}``, a list of mapping objects, or the
    compact ``{"source": ["truth", ...]}`` form. TSV columns are
    ``source_id``, ``truth_id``/``truth_ids``, and optional ``feature_type`` and
    ``status``.
    """

    if path is None:
        return [], {"kind": "exact-id", "path": None}, None
    path = Path(path)
    raw_document = None
    if path.suffix.lower() == ".json":
        document = json.loads(path.read_text(encoding="utf-8"))
        raw_document = document
        raw = document.get("mappings", document) if isinstance(
            document, Mapping
        ) else document
        if isinstance(raw, Mapping):
            items = [
                {"source_id": source_id, "truth_ids": truth_ids}
                for source_id, truth_ids in raw.items()
                if source_id not in {"schema_version", "method", "metadata"}
            ]
        elif isinstance(raw, list):
            items = raw
        else:
            raise ValueError(f"{path}: unsupported ortholog JSON shape")
        entries = [_normalise_mapping_item(item) for item in items]
    else:
        with path.open("r", encoding="utf-8", newline="") as handle:
            entries = [
                _normalise_mapping_item(row)
                for row in csv.DictReader(handle, delimiter="\t")
            ]
    seen: set[tuple[str, str | None]] = set()
    for entry in entries:
        if entry.status not in MAPPING_STATUSES:
            raise ValueError(
                f"{path}: unsupported mapping status {entry.status!r}"
            )
        key = (entry.source_id, entry.feature_type)
        if key in seen:
            raise ValueError(
                f"{path}: duplicate mapping for source {entry.source_id!r} "
                f"at feature type {entry.feature_type!r}"
            )
        seen.add(key)
    return entries, {
        "kind": "ortholog-map",
        "path": str(path.resolve()),
        "entries": len(entries),
        "sha256": _sha256_file(path),
    }, raw_document


def load_ortholog_map(
    path: str | Path | None,
) -> tuple[list[MappingEntry], dict[str, Any]]:
    """Load a mapping while preserving the historical two-value API."""

    entries, summary, _document = _load_ortholog_map(path)
    return entries, summary


def _safe_ratio(numerator: int, denominator: int) -> float | None:
    return numerator / denominator if denominator else None


def _prf(true_positive: int, predicted: int, expected: int) -> dict[str, Any]:
    false_positive = predicted - true_positive
    false_negative = expected - true_positive
    precision = _safe_ratio(true_positive, predicted)
    recall = _safe_ratio(true_positive, expected)
    if precision is None or recall is None:
        f1 = None
    elif precision + recall:
        f1 = 2 * precision * recall / (precision + recall)
    else:
        f1 = 0.0
    return {
        "true_positive": true_positive,
        "false_positive": false_positive,
        "false_negative": false_negative,
        "predicted": predicted,
        "expected": expected,
        "precision": precision,
        "recall": recall,
        "f1": f1,
    }


def _source_id(identifier: str, source_ids: set[str]) -> str | None:
    if identifier in source_ids:
        return identifier
    candidate = identifier
    while "_" in candidate:
        prefix, suffix = candidate.rsplit("_", 1)
        if not suffix.isdigit():
            break
        if prefix in source_ids:
            return prefix
        candidate = prefix
    return None


def _reciprocal_overlap(left: Model, right: Model) -> tuple[float, float]:
    if left.seqid != right.seqid:
        return (0.0, 0.0)
    overlap = max(0, min(left.end, right.end) - max(left.start, right.start) + 1)
    return (
        overlap / (left.end - left.start + 1),
        overlap / (right.end - right.start + 1),
    )


def _match_models(
    predictions: Sequence[Model],
    truths: Sequence[Model],
    *,
    minimum_reciprocal_overlap: float,
) -> list[tuple[Model, Model]]:
    ordered_predictions = sorted(predictions, key=lambda model: model.identifier)
    ordered_truths = sorted(truths, key=lambda model: model.identifier)
    candidates: dict[int, list[tuple[int, tuple[float, int, int]]]] = defaultdict(list)
    for prediction_index, prediction in enumerate(ordered_predictions):
        for truth_index, truth in enumerate(ordered_truths):
            left, right = _reciprocal_overlap(prediction, truth)
            if min(left, right) < minimum_reciprocal_overlap:
                continue
            candidates[prediction_index].append((truth_index, (
                min(left, right),
                int(prediction.strand == truth.strand),
                -(abs(prediction.start - truth.start)
                  + abs(prediction.end - truth.end)),
            )))
    for values in candidates.values():
        values.sort(
            key=lambda item: (
                item[1],
                ordered_truths[item[0]].identifier,
            ),
            reverse=True,
        )

    # Ortholog copy groups are normally tiny. For up to 12 target models,
    # solve the assignment exactly: first maximize matched copies, then total
    # reciprocal overlap, strand agreement, and boundary proximity.
    if len(ordered_truths) <= 12 and len(ordered_predictions) <= 64:
        @lru_cache(maxsize=None)
        def solve(
            prediction_index: int,
            used_truths: int,
        ) -> tuple[tuple[int, float, int, int], tuple[tuple[int, int], ...]]:
            if prediction_index == len(ordered_predictions):
                return (0, 0.0, 0, 0), ()
            best_score, best_pairs = solve(prediction_index + 1, used_truths)
            for truth_index, weight in candidates.get(prediction_index, ()):
                bit = 1 << truth_index
                if used_truths & bit:
                    continue
                tail_score, tail_pairs = solve(
                    prediction_index + 1, used_truths | bit,
                )
                score = (
                    tail_score[0] + 1,
                    tail_score[1] + weight[0],
                    tail_score[2] + weight[1],
                    tail_score[3] + weight[2],
                )
                pairs = ((prediction_index, truth_index), *tail_pairs)
                if score > best_score or (
                    score == best_score and pairs < best_pairs
                ):
                    best_score, best_pairs = score, pairs
            return best_score, best_pairs

        _score, selected = solve(0, 0)
    else:
        # Avoid exponential state for unusually large paralog families while
        # still guaranteeing maximum-cardinality one-to-one matching.
        truth_to_prediction: dict[int, int] = {}

        def augment(prediction_index: int, seen: set[int]) -> bool:
            for truth_index, _weight in candidates.get(prediction_index, ()):
                if truth_index in seen:
                    continue
                seen.add(truth_index)
                previous = truth_to_prediction.get(truth_index)
                if previous is None or augment(previous, seen):
                    truth_to_prediction[truth_index] = prediction_index
                    return True
            return False

        for prediction_index in sorted(
            range(len(ordered_predictions)),
            key=lambda index: (
                len(candidates.get(index, ())),
                ordered_predictions[index].identifier,
            ),
        ):
            augment(prediction_index, set())
        selected = tuple(sorted(
            (prediction_index, truth_index)
            for truth_index, prediction_index in truth_to_prediction.items()
        ))
    return [
        (ordered_predictions[prediction_index], ordered_truths[truth_index])
        for prediction_index, truth_index in selected
    ]


def _mapping_entries(
    entries: Iterable[MappingEntry],
    level: str,
) -> list[MappingEntry]:
    aliases = {
        "gene": {"gene", "locus"},
        "transcript": {"transcript", "mrna", "rna"},
    }
    return [
        entry for entry in entries
        if entry.feature_type in aliases[level]
    ]


def _infer_mapping_types(
    entries: Sequence[MappingEntry],
    prediction: Annotation,
    truth: Annotation,
) -> list[MappingEntry]:
    return _infer_mapping_types_from_ids(
        entries,
        source_genes=set(prediction.genes),
        source_transcripts=set(prediction.transcripts),
        truth_genes=set(truth.genes),
        truth_transcripts=set(truth.transcripts),
    )


def _infer_mapping_types_from_ids(
    entries: Sequence[MappingEntry],
    *,
    source_genes: set[str],
    source_transcripts: set[str],
    truth_genes: set[str],
    truth_transcripts: set[str],
) -> list[MappingEntry]:
    inferred = []
    for entry in entries:
        if entry.feature_type is not None:
            aliases = {
                "gene": "gene",
                "locus": "gene",
                "transcript": "transcript",
                "mrna": "transcript",
                "rna": "transcript",
            }
            feature_type = aliases.get(entry.feature_type)
            if feature_type is None:
                raise ValueError(
                    f"unsupported ortholog mapping feature_type "
                    f"{entry.feature_type!r} for {entry.source_id!r}"
                )
            inferred.append(MappingEntry(
                source_id=entry.source_id,
                truth_ids=entry.truth_ids,
                feature_type=feature_type,
                status=entry.status,
            ))
            continue
        gene_evidence = (
            entry.source_id in source_genes
            or any(identifier in truth_genes for identifier in entry.truth_ids)
        )
        transcript_evidence = (
            entry.source_id in source_transcripts
            or any(
                identifier in truth_transcripts
                for identifier in entry.truth_ids
            )
        )
        if gene_evidence == transcript_evidence:
            raise ValueError(
                f"cannot infer feature_type for ortholog mapping "
                f"{entry.source_id!r}; specify gene or transcript"
            )
        inferred.append(MappingEntry(
            source_id=entry.source_id,
            truth_ids=entry.truth_ids,
            feature_type="gene" if gene_evidence else "transcript",
            status=entry.status,
        ))
    return inferred


def _validate_mapping_entries_against_scope(
    entries: Sequence[MappingEntry],
    *,
    source_genes: set[str],
    source_transcripts: set[str],
    truth_genes: set[str],
    truth_transcripts: set[str],
) -> dict[str, Any]:
    source_ids = {
        "gene": source_genes,
        "transcript": source_transcripts,
    }
    truth_ids = {
        "gene": truth_genes,
        "transcript": truth_transcripts,
    }
    retained_targets = {"gene": set(), "transcript": set()}
    counts = Counter()
    for entry in entries:
        feature_type = entry.feature_type
        if feature_type not in {"gene", "transcript"}:
            raise ValueError(
                f"ortholog mapping {entry.source_id!r} has no canonical type"
            )
        if entry.source_id not in source_ids[feature_type]:
            raise ValueError(
                f"ortholog mapping names unknown source {feature_type} ID "
                f"{entry.source_id!r}"
            )
        if entry.status == "retained":
            if not entry.truth_ids:
                raise ValueError(
                    f"retained {feature_type} mapping {entry.source_id!r} "
                    "has no target truth ID"
                )
            for truth_id in entry.truth_ids:
                if truth_id not in truth_ids[feature_type]:
                    raise ValueError(
                        f"{feature_type} ortholog mapping for "
                        f"{entry.source_id!r} names missing target truth ID "
                        f"{truth_id!r}"
                    )
                if truth_id in retained_targets[feature_type]:
                    raise ValueError(
                        f"{feature_type} target truth ID {truth_id!r} is "
                        "mapped to more than one source group"
                    )
                retained_targets[feature_type].add(truth_id)
        elif entry.truth_ids:
            raise ValueError(
                f"{entry.status} {feature_type} mapping "
                f"{entry.source_id!r} must not declare target truth IDs"
            )
        counts[(feature_type, entry.status)] += 1
    return {
        "source_scope_validated": True,
        "source_genes": len(source_genes),
        "source_transcripts": len(source_transcripts),
        "truth_genes": len(truth_genes),
        "truth_transcripts": len(truth_transcripts),
        "mapping_counts": {
            f"{feature_type}.{status}": count
            for (feature_type, status), count in sorted(counts.items())
        },
    }


def _mapped_groups(
    predictions: Mapping[str, Model],
    truths: Mapping[str, Model],
    entries: Sequence[MappingEntry],
    *,
    level: str,
) -> tuple[list[_Group], set[str], set[str]]:
    typed = _mapping_entries(entries, level)
    source_ids = {entry.source_id for entry in typed if entry.scorable}
    prediction_groups: dict[str, list[str]] = defaultdict(list)
    for identifier in predictions:
        source = _source_id(identifier, source_ids)
        if source is not None:
            prediction_groups[source].append(identifier)
    groups = []
    expected_truth_ids: set[str] = set()
    for entry in typed:
        if not entry.scorable:
            continue
        missing = [item for item in entry.truth_ids if item not in truths]
        if missing:
            raise ValueError(
                f"{level} ortholog mapping for {entry.source_id!r} names "
                f"missing target truth ID {missing[0]!r}"
            )
        overlap = expected_truth_ids.intersection(entry.truth_ids)
        if overlap:
            raise ValueError(
                f"{level} target truth ID {sorted(overlap)[0]!r} is mapped "
                "to more than one source group"
            )
        expected_truth_ids.update(entry.truth_ids)
        groups.append(_Group(
            key=entry.source_id,
            prediction_ids=tuple(sorted(prediction_groups[entry.source_id])),
            truth_ids=tuple(sorted(entry.truth_ids)),
        ))
    return groups, source_ids, expected_truth_ids


def _exact_id_groups(
    predictions: Mapping[str, Model],
    truths: Mapping[str, Model],
) -> list[_Group]:
    source_ids = set(truths)
    grouped_predictions: dict[str, list[str]] = defaultdict(list)
    for identifier in predictions:
        source = _source_id(identifier, source_ids) or identifier
        grouped_predictions[source].append(identifier)
        source_ids.add(source)
    return [
        _Group(
            key=identifier,
            prediction_ids=tuple(sorted(grouped_predictions.get(identifier, ()))),
            truth_ids=(identifier,) if identifier in truths else (),
        )
        for identifier in sorted(source_ids)
    ]


def _gene_inferred_transcript_groups(
    prediction: Annotation,
    truth: Annotation,
    entries: Sequence[MappingEntry],
    claimed_predictions: set[str],
    claimed_truths: set[str],
) -> list[_Group]:
    groups = []
    gene_entries = _mapping_entries(entries, "gene")
    source_gene_ids = {entry.source_id for entry in gene_entries if entry.scorable}
    prediction_gene_groups: dict[str, list[str]] = defaultdict(list)
    for gene_id, transcript_ids in prediction.transcripts_by_gene.items():
        source_gene = _source_id(gene_id, source_gene_ids)
        if source_gene is not None:
            prediction_gene_groups[source_gene].extend(transcript_ids)
    for entry in gene_entries:
        if not entry.scorable:
            continue
        prediction_ids = sorted({
            identifier
            for identifier in prediction_gene_groups.get(entry.source_id, ())
            if identifier not in claimed_predictions
        })
        truth_ids = sorted({
            identifier
            for truth_gene in entry.truth_ids
            for identifier in truth.transcripts_by_gene.get(truth_gene, ())
            if identifier not in claimed_truths
        })
        if prediction_ids or truth_ids:
            groups.append(_Group(
                key=f"gene:{entry.source_id}",
                prediction_ids=tuple(prediction_ids),
                truth_ids=tuple(truth_ids),
            ))
    return groups


def _groups_for_level(
    prediction: Annotation,
    truth: Annotation,
    entries: Sequence[MappingEntry],
    *,
    level: str,
) -> list[_Group]:
    predictions = prediction.genes if level == "gene" else prediction.transcripts
    truths = truth.genes if level == "gene" else truth.transcripts
    if not entries:
        return _exact_id_groups(predictions, truths)
    groups, _source_ids, _truth_ids = _mapped_groups(
        predictions, truths, entries, level=level,
    )
    if level == "transcript":
        claimed_predictions = {
            identifier for group in groups for identifier in group.prediction_ids
        }
        claimed_truths = {
            identifier for group in groups for identifier in group.truth_ids
        }
        groups.extend(_gene_inferred_transcript_groups(
            prediction,
            truth,
            entries,
            claimed_predictions,
            claimed_truths,
        ))
    # Once an explicit map is supplied, never widen its declared scope by
    # falling back to coincidentally identical IDs. An empty level is reported
    # as empty scope and therefore cannot produce a misleading perfect score.
    return groups


def _intervals(
    annotation: Annotation,
    transcript: Model,
    feature_type: str,
) -> list[tuple[str, int, int, str]]:
    source = annotation.exons if feature_type == "exon" else annotation.cds
    return list(source.get(transcript.identifier, ()))


def _introns(
    annotation: Annotation,
    transcript: Model,
) -> list[tuple[str, int, int, str]]:
    exons = sorted(
        annotation.exons.get(transcript.identifier, ()),
        key=lambda interval: (interval[0], interval[1], interval[2]),
    )
    introns = []
    for left, right in zip(exons, exons[1:]):
        if left[0] != right[0] or left[2] >= right[1] - 1:
            continue
        introns.append((left[0], left[2] + 1, right[1] - 1, transcript.strand))
    if transcript.strand == "-":
        introns.reverse()
    return introns


def _interval_metric(
    prediction: Annotation,
    truth: Annotation,
    matches: Sequence[tuple[Model, Model]],
    unmatched_predictions: Sequence[Model],
    unmatched_truths: Sequence[Model],
    feature_type: str,
) -> dict[str, Any]:
    true_positive = predicted = expected = 0
    getter = (
        _introns
        if feature_type == "intron"
        else lambda annotation, transcript: _intervals(
            annotation, transcript, feature_type,
        )
    )
    for prediction_model, truth_model in matches:
        prediction_intervals = getter(prediction, prediction_model)
        truth_intervals = getter(truth, truth_model)
        predicted += len(prediction_intervals)
        expected += len(truth_intervals)
        remaining = list(truth_intervals)
        for interval in prediction_intervals:
            if interval in remaining:
                true_positive += 1
                remaining.remove(interval)
    predicted += sum(
        len(getter(prediction, model)) for model in unmatched_predictions
    )
    expected += sum(len(getter(truth, model)) for model in unmatched_truths)
    return _prf(true_positive, predicted, expected)


def _score_level(
    prediction: Annotation,
    truth: Annotation,
    groups: Sequence[_Group],
    *,
    level: str,
    minimum_reciprocal_overlap: float,
) -> tuple[dict[str, Any], list[tuple[Model, Model]], list[Model], list[Model]]:
    predictions = prediction.genes if level == "gene" else prediction.transcripts
    truths = truth.genes if level == "gene" else truth.transcripts
    matches = []
    unmatched_predictions = []
    unmatched_truths = []
    copy_true_positive = 0
    copy_exact = 0
    for group in groups:
        prediction_models = [
            predictions[identifier]
            for identifier in group.prediction_ids
            if identifier in predictions
        ]
        truth_models = [
            truths[identifier]
            for identifier in group.truth_ids
            if identifier in truths
        ]
        group_matches = _match_models(
            prediction_models,
            truth_models,
            minimum_reciprocal_overlap=minimum_reciprocal_overlap,
        )
        matches.extend(group_matches)
        matched_predictions = {item[0].identifier for item in group_matches}
        matched_truths = {item[1].identifier for item in group_matches}
        unmatched_predictions.extend(
            model for model in prediction_models
            if model.identifier not in matched_predictions
        )
        unmatched_truths.extend(
            model for model in truth_models
            if model.identifier not in matched_truths
        )
        copy_true_positive += min(len(prediction_models), len(truth_models))
        copy_exact += int(len(prediction_models) == len(truth_models))

    n_predictions = sum(len(group.prediction_ids) for group in groups)
    n_truths = sum(len(group.truth_ids) for group in groups)
    strand_correct = sum(
        prediction_model.strand == truth_model.strand
        for prediction_model, truth_model in matches
    )
    metrics = {
        "locus": _prf(len(matches), n_predictions, n_truths),
        "strand": _prf(strand_correct, n_predictions, n_truths),
        "copy": _prf(copy_true_positive, n_predictions, n_truths),
        "copy_count_exact": {
            "groups_exact": copy_exact,
            "groups_total": len(groups),
            "rate": _safe_ratio(copy_exact, len(groups)),
        },
    }
    return metrics, matches, unmatched_predictions, unmatched_truths


def _scope_document(
    groups: Sequence[_Group],
    predictions: Mapping[str, Model],
    truths: Mapping[str, Model],
) -> dict[str, int]:
    prediction_ids = {
        identifier for group in groups for identifier in group.prediction_ids
    }
    truth_ids = {
        identifier for group in groups for identifier in group.truth_ids
    }
    return {
        "groups": len(groups),
        "predicted_scored": len(prediction_ids),
        "expected_scored": len(truth_ids),
        "prediction_models_total": len(predictions),
        "truth_models_total": len(truths),
        "prediction_models_ignored": len(set(predictions) - prediction_ids),
        "truth_models_ignored": len(set(truths) - truth_ids),
    }


def score_target_truth(
    prediction_gff: str | Path,
    truth_gff: str | Path,
    *,
    ortholog_map: str | Path | None = None,
    source_gff: str | Path | None = None,
    id_policy: str = "ortholog-map",
    minimum_reciprocal_overlap: float = 0.5,
) -> dict[str, Any]:
    """Score a prediction against independent target-coordinate truth."""

    if not 0 < minimum_reciprocal_overlap <= 1:
        raise ValueError("minimum_reciprocal_overlap must be in (0, 1]")
    if id_policy not in {"ortholog-map", "exact-id"}:
        raise ValueError("id_policy must be 'ortholog-map' or 'exact-id'")
    if id_policy == "ortholog-map" and ortholog_map is None:
        raise ValueError(
            "independent target-truth scoring requires an explicit non-empty "
            "ortholog map; use id_policy='exact-id' only for deliberately "
            "same-ID truth"
        )
    if id_policy == "exact-id" and ortholog_map is not None:
        raise ValueError("exact-id policy cannot be combined with ortholog_map")
    prediction = parse_annotation(prediction_gff)
    truth = parse_annotation(truth_gff)
    entries, mapping_document, raw_mapping_document = _load_ortholog_map(
        ortholog_map
    )
    if id_policy == "ortholog-map" and not entries:
        raise ValueError("ortholog map contains no mapping entries")
    mapping_validation = None
    if id_policy == "ortholog-map" and source_gff is not None:
        from . import ortholog_scope

        try:
            source_hierarchy = ortholog_scope.parse_hierarchy(source_gff)
            entries = _infer_mapping_types_from_ids(
                entries,
                source_genes=set(source_hierarchy.genes),
                source_transcripts=set(source_hierarchy.transcripts),
                truth_genes=set(truth.genes),
                truth_transcripts=set(truth.transcripts),
            )
            mapping_validation = _validate_mapping_entries_against_scope(
                entries,
                source_genes=set(source_hierarchy.genes),
                source_transcripts=set(source_hierarchy.transcripts),
                truth_genes=set(truth.genes),
                truth_transcripts=set(truth.transcripts),
            )
            is_protein_rbh_scope = (
                isinstance(raw_mapping_document, Mapping)
                and raw_mapping_document.get("schema_version")
                == ortholog_scope.SCHEMA_VERSION
                and raw_mapping_document.get("method")
                == ortholog_scope.METHOD
            )
            if is_protein_rbh_scope:
                ortholog_scope.validate_mapping_against_annotations(
                    raw_mapping_document,
                    source_annotation=source_gff,
                    target_annotation=truth_gff,
                )
                mapping_validation["producer_scope_validated"] = True
                mapping_validation["input_fingerprints_validated"] = True
        except ortholog_scope.ScopeBuildError as exc:
            raise ValueError(
                f"ortholog map does not match its annotation scope: {exc}"
            ) from exc
    else:
        entries = _infer_mapping_types(entries, prediction, truth)
        if (
            id_policy == "ortholog-map"
            and isinstance(raw_mapping_document, Mapping)
        ):
            from . import ortholog_scope

            if (
                raw_mapping_document.get("schema_version")
                == ortholog_scope.SCHEMA_VERSION
                and raw_mapping_document.get("method")
                == ortholog_scope.METHOD
            ):
                raise ValueError(
                    "protein-RBH ortholog scope validation requires source_gff"
                )
    if mapping_validation is not None:
        mapping_document["semantic_validation"] = mapping_validation
    gene_groups = _groups_for_level(
        prediction, truth, entries, level="gene",
    )
    transcript_groups = _groups_for_level(
        prediction, truth, entries, level="transcript",
    )
    gene_metrics, _gene_matches, _unmatched_genes, _missing_genes = _score_level(
        prediction,
        truth,
        gene_groups,
        level="gene",
        minimum_reciprocal_overlap=minimum_reciprocal_overlap,
    )
    transcript_metrics, transcript_matches, unmatched_predictions, \
        unmatched_truths = _score_level(
            prediction,
            truth,
            transcript_groups,
            level="transcript",
            minimum_reciprocal_overlap=minimum_reciprocal_overlap,
        )

    exact_intron_chains = sum(
        _introns(prediction, prediction_model)
        == _introns(truth, truth_model)
        for prediction_model, truth_model in transcript_matches
    )
    result = {
        "schema_version": SCHEMA_VERSION,
        "method": METHOD,
        "inputs": {
            "prediction_gff": str(Path(prediction_gff).resolve()),
            "truth_gff": str(Path(truth_gff).resolve()),
            "source_gff": (
                str(Path(source_gff).resolve())
                if source_gff is not None
                else None
            ),
            "mapping": mapping_document,
        },
        "parameters": {
            "minimum_reciprocal_overlap": minimum_reciprocal_overlap,
            "id_policy": id_policy,
            "mapping_required": id_policy == "ortholog-map",
            "mapping_requirement_satisfied": (
                id_policy == "exact-id" or bool(entries)
            ),
            "mapping_source_scope_validated": (
                mapping_validation is not None
            ),
        },
        "scope": {
            "gene_groups": len(gene_groups),
            "transcript_groups": len(transcript_groups),
            "mapping_entries": len(entries),
            "mapping_status_counts": dict(sorted(Counter(
                entry.status for entry in entries
            ).items())),
            "gene": _scope_document(
                gene_groups, prediction.genes, truth.genes,
            ),
            "transcript": _scope_document(
                transcript_groups,
                prediction.transcripts,
                truth.transcripts,
            ),
        },
        "gene": gene_metrics,
        "transcript": transcript_metrics,
        "structure": {
            "intron_chain": _prf(
                exact_intron_chains,
                sum(len(group.prediction_ids) for group in transcript_groups),
                sum(len(group.truth_ids) for group in transcript_groups),
            ),
            "intron": _interval_metric(
                prediction,
                truth,
                transcript_matches,
                unmatched_predictions,
                unmatched_truths,
                "intron",
            ),
            "exon": _interval_metric(
                prediction,
                truth,
                transcript_matches,
                unmatched_predictions,
                unmatched_truths,
                "exon",
            ),
            "CDS": _interval_metric(
                prediction,
                truth,
                transcript_matches,
                unmatched_predictions,
                unmatched_truths,
                "CDS",
            ),
        },
    }
    return result


def write_target_truth_metrics(path: str | Path, document: Mapping[str, Any]) -> Path:
    """Atomically publish a target-truth metrics document."""

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    with temporary.open("w", encoding="utf-8") as handle:
        json.dump(document, handle, indent=2, sort_keys=True)
        handle.write("\n")
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(temporary, path)
    return path


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("prediction_gff", type=Path)
    parser.add_argument("truth_gff", type=Path)
    parser.add_argument("--ortholog-map", type=Path)
    parser.add_argument(
        "--source-gff",
        type=Path,
        help="source annotation used to validate every mapped source ID",
    )
    parser.add_argument(
        "--exact-id",
        action="store_true",
        help="Explicitly score same-ID truth without an ortholog map.",
    )
    parser.add_argument(
        "--minimum-reciprocal-overlap",
        type=float,
        default=0.5,
    )
    parser.add_argument("--output", type=Path)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    raw_arguments = list(argv) if argv is not None else sys.argv[1:]
    if raw_arguments and raw_arguments[0] == "filter":
        parser = argparse.ArgumentParser(
            description="Filter target truth to target FASTA sequence IDs."
        )
        parser.add_argument("truth_gff", type=Path)
        parser.add_argument("target_fasta", type=Path)
        parser.add_argument("output", type=Path)
        parser.add_argument("--force", action="store_true")
        arguments = parser.parse_args(raw_arguments[1:])
        document = filter_truth_to_target_fasta(
            arguments.truth_gff,
            arguments.target_fasta,
            arguments.output,
            force=arguments.force,
        )
        print(json.dumps(document, indent=2, sort_keys=True))
        return 0
    arguments = _parser().parse_args(raw_arguments)
    document = score_target_truth(
        arguments.prediction_gff,
        arguments.truth_gff,
        ortholog_map=arguments.ortholog_map,
        source_gff=arguments.source_gff,
        id_policy="exact-id" if arguments.exact_id else "ortholog-map",
        minimum_reciprocal_overlap=arguments.minimum_reciprocal_overlap,
    )
    if arguments.output is not None:
        write_target_truth_metrics(arguments.output, document)
    else:
        print(json.dumps(document, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
