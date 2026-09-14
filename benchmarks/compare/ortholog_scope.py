#!/usr/bin/env python3
"""Build deterministic protein-RBH scopes for independent target truth.

The builder deliberately keeps ortholog discovery separate from LiftOn.  It
extracts transcript proteins with gffread, searches reciprocal best hits with
MMseqs2, resolves transcript and gene conflicts deterministically, and emits
the mapping format consumed by :mod:`benchmarks.compare.target_truth`.

No command in this module downloads data.  Output directories are published
atomically only after every input, tool, hit, mapping, and retained-group
check succeeds.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import math
import os
import re
import shutil
import subprocess
import sys
import tempfile
import urllib.parse
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable, Mapping, Sequence

from lifton.gff3_validator import GENE_TYPES


SCHEMA_VERSION = 1
METHOD = "protein-rbh-ortholog-scope-v1"
PREPARE_METHOD = "gffread-protein-preparation-v1"
MMSEQS_FORMAT = "query,target,fident,alnlen,qcov,tcov,evalue,bits"
DEFAULT_MINIMUM_RETAINED_GROUPS = 100
REGISTRY_SCHEMA_VERSION = 1
MAPPING_FEATURE_TYPES = ("gene", "transcript")
MAPPING_STATUSES = {"retained", "unscored"}
_GENE_TYPES_LOWER = {value.lower() for value in GENE_TYPES}
_PROTEIN_TRANSCRIPT_TYPES_LOWER = {"mrna", "transcript"}
_GTF_ATTRIBUTE_RE = re.compile(
    r"""^\s*([^\s;]+)\s+(?:"([^"]*)"|([^;\s]+))\s*$"""
)


class ScopeBuildError(RuntimeError):
    """Raised when an ortholog scope cannot be built safely."""


@dataclass(frozen=True)
class Hierarchy:
    """Annotation gene/transcript hierarchy used by protein and subset scopes."""

    genes: tuple[str, ...]
    transcript_to_gene: Mapping[str, str]
    format_counts: Mapping[str, int]
    excluded_genes: tuple[str, ...] = ()
    excluded_transcripts: tuple[str, ...] = ()

    @property
    def transcripts(self) -> tuple[str, ...]:
        return tuple(sorted(self.transcript_to_gene))

    def as_document(self) -> dict[str, Any]:
        return {
            "genes": list(self.genes),
            "transcript_to_gene": dict(sorted(self.transcript_to_gene.items())),
            "format_counts": dict(sorted(self.format_counts.items())),
            "excluded_ambiguous": {
                "genes": list(self.excluded_genes),
                "transcripts": list(self.excluded_transcripts),
            },
        }


@dataclass(frozen=True)
class Hit:
    source_id: str
    target_id: str
    identity: float
    alignment_length: int
    query_coverage: float
    target_coverage: float
    evalue: float
    bits: float
    line_number: int

    def metrics(self) -> dict[str, Any]:
        return {
            "identity": self.identity,
            "alignment_length": self.alignment_length,
            "query_coverage": self.query_coverage,
            "target_coverage": self.target_coverage,
            "evalue": self.evalue,
            "bits": self.bits,
        }


Runner = Callable[..., subprocess.CompletedProcess[str]]


def sha256_file(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def canonical_sha256(value: Any) -> str:
    payload = json.dumps(
        value,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=True,
        allow_nan=False,
    ).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def _json_bytes(document: Any) -> bytes:
    return (
        json.dumps(
            document,
            indent=2,
            sort_keys=True,
            ensure_ascii=True,
            allow_nan=False,
        )
        + "\n"
    ).encode("utf-8")


def _write_bytes(path: Path, payload: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("xb") as handle:
        handle.write(payload)
        handle.flush()
        os.fsync(handle.fileno())


def _write_json(path: Path, document: Any) -> None:
    _write_bytes(path, _json_bytes(document))


def _atomic_write_json(path: Path, document: Any) -> None:
    path = Path(path).resolve()
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists() or path.is_symlink():
        raise FileExistsError(f"refusing to replace ortholog artifact: {path}")
    temporary = path.with_name(f".{path.name}.tmp-{os.getpid()}")
    if temporary.exists() or temporary.is_symlink():
        raise FileExistsError(f"temporary output already exists: {temporary}")
    try:
        _write_json(temporary, document)
        os.replace(temporary, path)
    except Exception:
        if temporary.exists() or temporary.is_symlink():
            temporary.unlink()
        raise


def _publish_directory(
    output_dir: Path,
    builder: Callable[[Path], Any],
) -> Any:
    output_dir = Path(output_dir).resolve()
    output_dir.parent.mkdir(parents=True, exist_ok=True)
    if output_dir.exists() or output_dir.is_symlink():
        raise FileExistsError(
            f"refusing to replace ortholog output directory: {output_dir}"
        )
    temporary = Path(tempfile.mkdtemp(
        prefix=f".{output_dir.name}.tmp-",
        dir=output_dir.parent,
    ))
    try:
        result = builder(temporary)
        os.replace(temporary, output_dir)
        return result
    except Exception:
        shutil.rmtree(temporary, ignore_errors=True)
        raise


def _input_record(path: str | Path) -> dict[str, Any]:
    resolved = Path(path).expanduser().resolve()
    if not resolved.is_file():
        raise ScopeBuildError(f"input is missing or empty: {resolved}")
    before = resolved.stat()
    if before.st_size <= 0:
        raise ScopeBuildError(f"input is missing or empty: {resolved}")
    digest = sha256_file(resolved)
    after = resolved.stat()
    signature_fields = (
        "st_size",
        "st_mtime_ns",
        "st_ctime_ns",
        "st_dev",
        "st_ino",
    )
    if any(
        getattr(before, name) != getattr(after, name)
        for name in signature_fields
    ):
        raise ScopeBuildError(f"input changed while it was hashed: {resolved}")
    return {
        "path": str(resolved),
        "size": after.st_size,
        "sha256": digest,
        "mtime_ns": after.st_mtime_ns,
        "ctime_ns": after.st_ctime_ns,
        "st_dev": after.st_dev,
        "st_ino": after.st_ino,
    }


def _artifact_record(path: Path) -> dict[str, Any]:
    if not path.is_file() or path.stat().st_size <= 0:
        raise ScopeBuildError(f"artifact is missing or empty: {path}")
    return {
        "name": path.name,
        "size": path.stat().st_size,
        "sha256": sha256_file(path),
    }


def _assert_frozen_records(records: Mapping[str, Mapping[str, Any]]) -> None:
    identity_fields = ("path", "size", "sha256")
    for label, expected in records.items():
        observed = _input_record(expected["path"])
        if any(
            observed.get(field) != expected.get(field)
            for field in identity_fields
        ):
            changed = [
                field
                for field in sorted(set(expected) | set(observed))
                if expected.get(field) != observed.get(field)
            ]
            detail = "; ".join(
                f"{field}: expected {expected.get(field)!r}, "
                f"observed {observed.get(field)!r}"
                for field in changed
            )
            raise ScopeBuildError(
                f"frozen input changed during ortholog build: {label} "
                f"({detail})"
            )


def _assert_frozen_tools(tools: Mapping[str, Mapping[str, Any]]) -> None:
    for label, record in tools.items():
        path = Path(record["executable"])
        if sha256_file(path) != record["executable_sha256"]:
            raise ScopeBuildError(
                f"tool executable changed during ortholog build: {label}"
            )


def _resolve_executable(value: str | Path) -> Path:
    text = str(value)
    candidate = (
        shutil.which(text)
        if os.sep not in text and (os.altsep is None or os.altsep not in text)
        else text
    )
    if candidate is None:
        raise ScopeBuildError(f"required executable is not on PATH: {text}")
    path = Path(candidate).expanduser().resolve()
    if not path.is_file() or not os.access(path, os.X_OK):
        raise ScopeBuildError(f"required executable is not executable: {path}")
    return path


def _run_checked(
    argv: Sequence[str | Path],
    *,
    runner: Runner,
    cwd: Path | None = None,
    label: str,
) -> subprocess.CompletedProcess[str]:
    command = [str(value) for value in argv]
    try:
        result = runner(
            command,
            cwd=str(cwd) if cwd is not None else None,
            text=True,
            capture_output=True,
            check=False,
        )
    except OSError as exc:
        raise ScopeBuildError(f"{label} could not be launched: {exc}") from exc
    if not isinstance(result.returncode, int) or result.returncode != 0:
        stderr = str(result.stderr or "")[-2000:].strip()
        raise ScopeBuildError(
            f"{label} failed with exit code {result.returncode}: {stderr}"
        )
    return result


def _tool_record(
    executable: str | Path,
    *,
    version_arguments: Sequence[str],
    runner: Runner,
) -> dict[str, Any]:
    resolved = _resolve_executable(executable)
    result = _run_checked(
        [resolved, *version_arguments],
        runner=runner,
        label=f"{resolved.name} version probe",
    )
    version = "\n".join(
        part.strip()
        for part in (result.stdout, result.stderr)
        if isinstance(part, str) and part.strip()
    )
    if not version:
        raise ScopeBuildError(f"{resolved.name} returned an empty version")
    return {
        "executable": str(resolved),
        "executable_sha256": sha256_file(resolved),
        "version": version,
        "version_argv": [str(resolved), *version_arguments],
    }


def _gff3_attributes(text: str) -> dict[str, tuple[str, ...]]:
    result: dict[str, tuple[str, ...]] = {}
    for raw in text.split(";"):
        raw = raw.strip()
        if not raw:
            continue
        if "=" not in raw:
            raise ScopeBuildError(
                f"malformed GFF3 attribute without '=': {raw!r}"
            )
        key, value = raw.split("=", 1)
        key = key.strip()
        if not key or key in result:
            raise ScopeBuildError(f"duplicate or empty GFF3 attribute: {key!r}")
        result[key] = tuple(
            urllib.parse.unquote(item.strip(), errors="strict")
            for item in value.split(",")
            if item.strip()
        )
    return result


def _gtf_attributes(text: str) -> dict[str, tuple[str, ...]]:
    result: dict[str, tuple[str, ...]] = {}
    for raw in text.split(";"):
        raw = raw.strip()
        if not raw:
            continue
        match = _GTF_ATTRIBUTE_RE.fullmatch(raw)
        if match is None:
            raise ScopeBuildError(f"malformed GTF attribute: {raw!r}")
        key = match.group(1)
        value = match.group(2) if match.group(2) is not None else match.group(3)
        result[key] = (*result.get(key, ()), value)
    return result


def _one_attribute(
    attributes: Mapping[str, tuple[str, ...]],
    name: str,
    *,
    path: Path,
    line_number: int,
    required: bool,
) -> str | None:
    values = attributes.get(name, ())
    if not values and not required:
        return None
    if len(values) != 1 or not values[0]:
        raise ScopeBuildError(
            f"{path}: line {line_number} requires exactly one {name}"
        )
    value = values[0]
    if any(character.isspace() for character in value):
        raise ScopeBuildError(
            f"{path}: line {line_number} {name} contains whitespace"
        )
    return value


def parse_hierarchy(path: str | Path) -> Hierarchy:
    """Parse reproducible GFF3 or GTF gene/transcript identifiers."""

    path = Path(path).expanduser().resolve()
    _input_record(path)
    genes: set[str] = set()
    transcript_to_gene: dict[str, str] = {}
    model_records: set[tuple[str, str]] = set()
    duplicate_genes: set[str] = set()
    duplicate_transcripts: set[str] = set()
    format_counts = {"gff3": 0, "gtf": 0}
    opener = gzip.open if path.name.lower().endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8", errors="strict") as handle:
        for line_number, raw in enumerate(handle, start=1):
            if not raw.strip() or raw.startswith("#"):
                continue
            columns = raw.rstrip("\r\n").split("\t")
            if len(columns) != 9:
                raise ScopeBuildError(
                    f"{path}: line {line_number} has {len(columns)} columns"
                )
            feature_type = columns[2].strip()
            feature_lower = feature_type.lower()
            text = columns[8].strip()
            if text in {"", "."}:
                is_gff3 = ".gtf" not in path.name.lower()
                attributes = {}
            else:
                is_gff3 = re.match(r"^[^;\s=]+\s*=", text) is not None
                attributes = (
                    _gff3_attributes(text)
                    if is_gff3
                    else _gtf_attributes(text)
                )
            annotation_format = "gff3" if is_gff3 else "gtf"
            format_counts[annotation_format] += 1
            is_gene = feature_lower in _GENE_TYPES_LOWER
            is_transcript = feature_lower in _PROTEIN_TRANSCRIPT_TYPES_LOWER
            if is_gff3:
                if is_gene:
                    gene_id = _one_attribute(
                        attributes,
                        "ID",
                        path=path,
                        line_number=line_number,
                        required=True,
                    )
                    key = ("gene", gene_id)
                    if key in model_records:
                        duplicate_genes.add(gene_id)
                    model_records.add(key)
                    genes.add(gene_id)
                elif is_transcript:
                    transcript_id = _one_attribute(
                        attributes,
                        "ID",
                        path=path,
                        line_number=line_number,
                        required=True,
                    )
                    parents = attributes.get("Parent", ())
                    if len(parents) != 1:
                        raise ScopeBuildError(
                            f"{path}: line {line_number} transcript "
                            "requires exactly one Parent gene"
                        )
                    gene_id = parents[0]
                    key = ("transcript", transcript_id)
                    if key in model_records:
                        duplicate_transcripts.add(transcript_id)
                    model_records.add(key)
                    genes.add(gene_id)
                    previous = transcript_to_gene.get(transcript_id)
                    if previous is not None and previous != gene_id:
                        duplicate_transcripts.add(transcript_id)
                    else:
                        transcript_to_gene[transcript_id] = gene_id
                continue

            gene_id = _one_attribute(
                attributes,
                "gene_id",
                path=path,
                line_number=line_number,
                required=is_gene or is_transcript or feature_lower in {
                    "exon", "cds",
                },
            )
            transcript_id = _one_attribute(
                attributes,
                "transcript_id",
                path=path,
                line_number=line_number,
                required=is_transcript or feature_lower in {"exon", "cds"},
            )
            if gene_id is not None:
                genes.add(gene_id)
            if transcript_id is None:
                continue
            previous = transcript_to_gene.get(transcript_id)
            if previous is not None and previous != gene_id:
                duplicate_transcripts.add(transcript_id)
            else:
                transcript_to_gene[transcript_id] = str(gene_id)
            if is_transcript:
                key = ("transcript", transcript_id)
                if key in model_records:
                    duplicate_transcripts.add(transcript_id)
                model_records.add(key)
            if is_gene:
                key = ("gene", str(gene_id))
                if key in model_records:
                    duplicate_genes.add(str(gene_id))
                model_records.add(key)

    if not genes or not transcript_to_gene:
        raise ScopeBuildError(
            f"{path}: annotation has no usable gene/transcript hierarchy"
        )
    if all(format_counts.values()):
        if path.name.lower().endswith((".gtf", ".gtf.gz")):
            raise ScopeBuildError(f"{path}: GTF contains mixed GFF3 attributes")
        raise ScopeBuildError(f"{path}: GFF3 contains mixed GTF attributes")
    missing_parents = set(transcript_to_gene.values()) - genes
    if missing_parents:
        raise ScopeBuildError(
            f"{path}: transcript parent is absent: "
            f"{sorted(missing_parents)[0]!r}"
        )
    excluded_transcripts = duplicate_transcripts | {
        transcript_id
        for transcript_id, gene_id in transcript_to_gene.items()
        if gene_id in duplicate_genes
    }
    transcript_to_gene = {
        transcript_id: gene_id
        for transcript_id, gene_id in transcript_to_gene.items()
        if transcript_id not in excluded_transcripts
    }
    genes = set(transcript_to_gene.values()) - duplicate_genes
    return Hierarchy(
        genes=tuple(sorted(genes)),
        transcript_to_gene=dict(sorted(transcript_to_gene.items())),
        format_counts={
            key: value for key, value in sorted(format_counts.items()) if value
        },
        excluded_genes=tuple(sorted(duplicate_genes)),
        excluded_transcripts=tuple(sorted(excluded_transcripts)),
    )


def _read_fasta(path: Path) -> list[tuple[str, str]]:
    records = []
    identifier = None
    sequence: list[str] = []
    with path.open("r", encoding="utf-8", errors="strict") as handle:
        for line_number, raw in enumerate(handle, start=1):
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if identifier is not None:
                    records.append((identifier, "".join(sequence)))
                header = line[1:].strip()
                if not header:
                    raise ScopeBuildError(
                        f"{path}: FASTA line {line_number} has an empty header"
                    )
                identifier = header.split()[0]
                sequence = []
            elif identifier is None:
                raise ScopeBuildError(
                    f"{path}: FASTA sequence precedes the first header"
                )
            else:
                sequence.append("".join(line.split()))
    if identifier is not None:
        records.append((identifier, "".join(sequence)))
    return records


def normalize_protein_fasta(
    raw_path: str | Path,
    output_path: str | Path,
    hierarchy: Hierarchy,
) -> dict[str, Any]:
    """Sort gffread proteins by transcript ID and normalize sequence bytes."""

    raw_path = Path(raw_path)
    output_path = Path(output_path)
    allowed_ids = set(hierarchy.transcripts)
    excluded_ids = set(hierarchy.excluded_transcripts)
    sequences: dict[str, str] = {}
    ignored_ids = []
    outside_hierarchy_ids = []
    invalid_ids = []
    for identifier, raw_sequence in _read_fasta(raw_path):
        if identifier not in allowed_ids:
            if identifier in excluded_ids:
                ignored_ids.append(identifier)
                continue
            outside_hierarchy_ids.append(identifier)
            continue
        if identifier in sequences:
            raise ScopeBuildError(
                f"{raw_path}: duplicate protein ID {identifier!r}"
            )
        sequence = raw_sequence.upper().rstrip("*")
        if (
            not sequence
            or "*" in sequence
            or any(not ("A" <= residue <= "Z") for residue in sequence)
        ):
            invalid_ids.append(identifier)
            continue
        sequences[identifier] = sequence
    if not sequences:
        raise ScopeBuildError(f"{raw_path}: gffread extracted no proteins")
    lines = []
    for identifier in sorted(sequences):
        lines.append(f">{identifier}")
        sequence = sequences[identifier]
        lines.extend(
            sequence[offset:offset + 60]
            for offset in range(0, len(sequence), 60)
        )
    _write_bytes(output_path, ("\n".join(lines) + "\n").encode("ascii"))
    return {
        "proteins": len(sequences),
        "residues": sum(len(sequence) for sequence in sequences.values()),
        "missing_transcript_proteins": sorted(allowed_ids - set(sequences)),
        "ids": sorted(sequences),
        "ignored_ambiguous_protein_ids": sorted(ignored_ids),
        "ignored_outside_hierarchy_protein_ids": sorted(
            outside_hierarchy_ids
        ),
        "excluded_invalid_protein_ids": sorted(invalid_ids),
    }


def _protein_summary(statistics: Mapping[str, Any]) -> dict[str, Any]:
    missing = statistics["missing_transcript_proteins"]
    return {
        "proteins": statistics["proteins"],
        "residues": statistics["residues"],
        "missing_transcript_proteins": len(missing),
        "missing_transcript_ids_sha256": canonical_sha256(missing),
        "ignored_ambiguous_proteins": len(
            statistics["ignored_ambiguous_protein_ids"]
        ),
        "ignored_outside_hierarchy_proteins": len(
            statistics["ignored_outside_hierarchy_protein_ids"]
        ),
        "ignored_outside_hierarchy_ids_sha256": canonical_sha256(
            statistics["ignored_outside_hierarchy_protein_ids"]
        ),
        "excluded_invalid_proteins": len(
            statistics["excluded_invalid_protein_ids"]
        ),
        "excluded_invalid_protein_ids_sha256": canonical_sha256(
            statistics["excluded_invalid_protein_ids"]
        ),
    }


def _logical_argv(argv: Sequence[str | Path], work: Path) -> list[str]:
    prefix = str(work.resolve())
    result = []
    for value in argv:
        text = str(value)
        if text == prefix:
            text = "$OUTPUT"
        elif text.startswith(prefix + os.sep):
            text = "$OUTPUT/" + text[len(prefix + os.sep):]
        result.append(text)
    return result


def _prepare_in_workspace(
    work: Path,
    *,
    source_annotation: Path,
    source_genome: Path,
    target_annotation: Path,
    target_genome: Path,
    gffread_tool: Mapping[str, Any],
    runner: Runner,
) -> dict[str, Any]:
    source_hierarchy = parse_hierarchy(source_annotation)
    target_hierarchy = parse_hierarchy(target_annotation)
    hierarchies = {
        "source": source_hierarchy,
        "target": target_hierarchy,
    }
    annotations = {
        "source": (source_annotation, source_genome),
        "target": (target_annotation, target_genome),
    }
    protein_stats = {}
    commands = []
    for label in ("source", "target"):
        annotation, genome = annotations[label]
        raw_fasta = work / f"{label}.proteins.raw.fa"
        normalized_fasta = work / f"{label}.proteins.fa"
        argv = [
            gffread_tool["executable"],
            "-y",
            raw_fasta,
            "-g",
            genome,
            annotation,
        ]
        _run_checked(
            argv,
            runner=runner,
            cwd=work,
            label=f"gffread {label} protein extraction",
        )
        if not raw_fasta.is_file() or raw_fasta.stat().st_size <= 0:
            raise ScopeBuildError(
                f"gffread did not create {label} protein FASTA"
            )
        protein_stats[label] = normalize_protein_fasta(
            raw_fasta,
            normalized_fasta,
            hierarchies[label],
        )
        raw_fasta.unlink()
        commands.append({
            "label": f"extract_{label}_proteins",
            "argv": _logical_argv(argv, work),
            "shell": False,
        })
        _write_json(
            work / f"{label}.hierarchy.json",
            hierarchies[label].as_document(),
        )
    return {
        "hierarchies": hierarchies,
        "protein_stats": protein_stats,
        "commands": commands,
    }


def _parameter_document(
    *,
    minimum_identity: float,
    minimum_reciprocal_coverage: float,
    maximum_evalue: float,
    threads: int,
    minimum_gene_groups: int,
    minimum_transcript_groups: int,
) -> dict[str, Any]:
    values = {
        "minimum_identity": minimum_identity,
        "minimum_reciprocal_coverage": minimum_reciprocal_coverage,
        "maximum_evalue": maximum_evalue,
        "threads": threads,
        "minimum_gene_groups": minimum_gene_groups,
        "minimum_transcript_groups": minimum_transcript_groups,
    }
    if (
        not isinstance(threads, int)
        or isinstance(threads, bool)
        or threads < 1
    ):
        raise ValueError("threads must be a positive integer")
    for name in ("minimum_gene_groups", "minimum_transcript_groups"):
        value = values[name]
        if (
            not isinstance(value, int)
            or isinstance(value, bool)
            or value < 0
        ):
            raise ValueError(f"{name} must be a non-negative integer")
    for name in ("minimum_identity", "minimum_reciprocal_coverage"):
        value = values[name]
        if (
            not isinstance(value, (int, float))
            or isinstance(value, bool)
            or not math.isfinite(float(value))
            or not 0.0 <= float(value) <= 1.0
        ):
            raise ValueError(f"{name} must be finite in [0, 1]")
        values[name] = float(value)
    if (
        not isinstance(maximum_evalue, (int, float))
        or isinstance(maximum_evalue, bool)
        or not math.isfinite(float(maximum_evalue))
        or float(maximum_evalue) < 0.0
    ):
        raise ValueError("maximum_evalue must be finite and non-negative")
    values["maximum_evalue"] = float(maximum_evalue)
    return values


def parse_mmseqs_hits(
    path: str | Path,
    *,
    source_proteins: set[str],
    target_proteins: set[str],
) -> list[Hit]:
    path = Path(path)
    if not path.is_file():
        raise ScopeBuildError(f"MMseqs output is missing: {path}")
    hits = []
    seen = set()
    with path.open("r", encoding="utf-8", errors="strict") as handle:
        for line_number, raw in enumerate(handle, start=1):
            if not raw.strip() or raw.startswith("#"):
                continue
            columns = raw.rstrip("\r\n").split("\t")
            if len(columns) != 8:
                raise ScopeBuildError(
                    f"{path}: line {line_number} has {len(columns)} columns, "
                    "expected 8"
                )
            source_id, target_id = columns[:2]
            if source_id not in source_proteins:
                raise ScopeBuildError(
                    f"{path}: line {line_number} has unknown source protein "
                    f"{source_id!r}"
                )
            if target_id not in target_proteins:
                raise ScopeBuildError(
                    f"{path}: line {line_number} has unknown target protein "
                    f"{target_id!r}"
                )
            try:
                identity = float(columns[2])
                alignment_length = int(columns[3])
                query_coverage = float(columns[4])
                target_coverage = float(columns[5])
                evalue = float(columns[6])
                bits = float(columns[7])
            except ValueError as exc:
                raise ScopeBuildError(
                    f"{path}: line {line_number} has non-numeric metrics"
                ) from exc
            metrics = (
                identity,
                query_coverage,
                target_coverage,
                evalue,
                bits,
            )
            if (
                alignment_length <= 0
                or any(not math.isfinite(value) for value in metrics)
                or not 0.0 <= identity <= 1.0
                or not 0.0 <= query_coverage <= 1.0
                or not 0.0 <= target_coverage <= 1.0
                or evalue < 0.0
                or bits < 0.0
            ):
                raise ScopeBuildError(
                    f"{path}: line {line_number} has out-of-range metrics"
                )
            key = (source_id, target_id)
            if key in seen:
                raise ScopeBuildError(
                    f"{path}: duplicate RBH pair {source_id!r}/{target_id!r}"
                )
            seen.add(key)
            hits.append(Hit(
                source_id=source_id,
                target_id=target_id,
                identity=identity,
                alignment_length=alignment_length,
                query_coverage=query_coverage,
                target_coverage=target_coverage,
                evalue=evalue,
                bits=bits,
                line_number=line_number,
            ))
    return hits


def _hit_filter_reason(
    hit: Hit,
    parameters: Mapping[str, Any],
) -> str | None:
    if hit.identity < parameters["minimum_identity"]:
        return "identity_below_minimum"
    if hit.query_coverage < parameters["minimum_reciprocal_coverage"]:
        return "query_coverage_below_minimum"
    if hit.target_coverage < parameters["minimum_reciprocal_coverage"]:
        return "target_coverage_below_minimum"
    if hit.evalue > parameters["maximum_evalue"]:
        return "evalue_above_maximum"
    return None


def _hit_quality_key(hit: Hit) -> tuple[Any, ...]:
    return (
        -hit.bits,
        hit.evalue,
        -hit.identity,
        -min(hit.query_coverage, hit.target_coverage),
        -hit.query_coverage,
        -hit.target_coverage,
        hit.source_id,
        hit.target_id,
    )


def _maximum_cardinality_selection(
    candidates: Sequence[Any],
    *,
    left_key: Callable[[Any], str],
    right_key: Callable[[Any], str],
    quality_key: Callable[[Any], tuple[Any, ...]],
) -> list[Any]:
    """Return a deterministic maximum-cardinality bipartite selection."""

    by_left: dict[str, list[Any]] = defaultdict(list)
    for candidate in candidates:
        by_left[left_key(candidate)].append(candidate)
    for values in by_left.values():
        values.sort(key=quality_key)
    right_to_candidate: dict[str, Any] = {}

    def augment(
        left: str,
        seen_left: set[str],
        seen_right: set[str],
    ) -> bool:
        if left in seen_left:
            return False
        seen_left.add(left)
        for candidate in by_left[left]:
            right = right_key(candidate)
            if right in seen_right:
                continue
            seen_right.add(right)
            previous = right_to_candidate.get(right)
            if previous is None or augment(
                left_key(previous),
                seen_left,
                seen_right,
            ):
                right_to_candidate[right] = candidate
                return True
        return False

    for left in sorted(
        by_left,
        key=lambda value: (
            len(by_left[value]),
            quality_key(by_left[value][0]),
            value,
        ),
    ):
        augment(left, set(), set())
    return sorted(right_to_candidate.values(), key=quality_key)


def _resolve_transcript_hits(hits: Sequence[Hit]) -> list[Hit]:
    selected = _maximum_cardinality_selection(
        hits,
        left_key=lambda hit: hit.source_id,
        right_key=lambda hit: hit.target_id,
        quality_key=_hit_quality_key,
    )
    return sorted(selected, key=lambda hit: (hit.source_id, hit.target_id))


def _resolve_gene_pairs(
    transcript_hits: Sequence[Hit],
    source_hierarchy: Hierarchy,
    target_hierarchy: Hierarchy,
) -> tuple[list[tuple[str, str]], dict[tuple[str, str], list[Hit]]]:
    evidence: dict[tuple[str, str], list[Hit]] = defaultdict(list)
    for hit in transcript_hits:
        pair = (
            source_hierarchy.transcript_to_gene[hit.source_id],
            target_hierarchy.transcript_to_gene[hit.target_id],
        )
        evidence[pair].append(hit)
    def quality(pair: tuple[str, str]) -> tuple[Any, ...]:
        return (
            -len(evidence[pair]),
            -sum(hit.bits for hit in evidence[pair]),
            -max(hit.bits for hit in evidence[pair]),
            min(hit.evalue for hit in evidence[pair]),
            pair[0],
            pair[1],
        )

    selected = _maximum_cardinality_selection(
        list(evidence),
        left_key=lambda pair: pair[0],
        right_key=lambda pair: pair[1],
        quality_key=quality,
    )
    return sorted(selected), evidence


def _unscored_mapping(
    source_id: str,
    feature_type: str,
    reason: str,
) -> dict[str, Any]:
    return {
        "source_id": source_id,
        "truth_ids": [],
        "feature_type": feature_type,
        "status": "unscored",
        "reason": reason,
    }


def build_mapping_document(
    hits: Sequence[Hit],
    *,
    source_hierarchy: Hierarchy,
    target_hierarchy: Hierarchy,
    source_protein_ids: set[str],
    target_protein_ids: set[str],
    parameters: Mapping[str, Any],
    provenance: Mapping[str, Any],
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    """Filter RBHs and build deterministic one-to-one mapping entries."""

    filter_reasons = {
        (hit.source_id, hit.target_id): _hit_filter_reason(hit, parameters)
        for hit in hits
    }
    passing = [
        hit for hit in hits
        if filter_reasons[(hit.source_id, hit.target_id)] is None
    ]
    transcript_selected = _resolve_transcript_hits(passing)
    transcript_selected_keys = {
        (hit.source_id, hit.target_id) for hit in transcript_selected
    }
    gene_pairs, gene_evidence = _resolve_gene_pairs(
        transcript_selected,
        source_hierarchy,
        target_hierarchy,
    )
    gene_pair_set = set(gene_pairs)
    retained_hits = [
        hit for hit in transcript_selected
        if (
            source_hierarchy.transcript_to_gene[hit.source_id],
            target_hierarchy.transcript_to_gene[hit.target_id],
        ) in gene_pair_set
    ]
    retained_by_source = {hit.source_id: hit for hit in retained_hits}
    retained_genes = {source: target for source, target in gene_pairs}
    passing_by_source: dict[str, list[Hit]] = defaultdict(list)
    selected_by_source = {hit.source_id: hit for hit in transcript_selected}
    for hit in passing:
        passing_by_source[hit.source_id].append(hit)

    entries = []
    source_transcripts_by_gene: dict[str, list[str]] = defaultdict(list)
    for transcript_id, gene_id in source_hierarchy.transcript_to_gene.items():
        source_transcripts_by_gene[gene_id].append(transcript_id)
    for source_gene in source_hierarchy.genes:
        target_gene = retained_genes.get(source_gene)
        if target_gene is not None:
            support = gene_evidence[(source_gene, target_gene)]
            entries.append({
                "source_id": source_gene,
                "truth_ids": [target_gene],
                "feature_type": "gene",
                "status": "retained",
                "evidence": {
                    "supporting_transcript_pairs": [
                        [hit.source_id, hit.target_id]
                        for hit in sorted(
                            support,
                            key=lambda value: (
                                value.source_id, value.target_id,
                            ),
                        )
                    ],
                    "support_count": len(support),
                    "total_bits": sum(hit.bits for hit in support),
                },
            })
            continue
        transcript_ids = source_transcripts_by_gene.get(source_gene, [])
        if not any(
            transcript_id in source_protein_ids
            for transcript_id in transcript_ids
        ):
            reason = "no_extracted_protein"
        elif not any(
            passing_by_source.get(transcript_id)
            for transcript_id in transcript_ids
        ):
            reason = "no_passing_rbh"
        else:
            reason = "gene_pair_conflict"
        entries.append(_unscored_mapping(source_gene, "gene", reason))

    for source_transcript in source_hierarchy.transcripts:
        retained = retained_by_source.get(source_transcript)
        if retained is not None:
            entries.append({
                "source_id": source_transcript,
                "truth_ids": [retained.target_id],
                "feature_type": "transcript",
                "status": "retained",
                "evidence": retained.metrics(),
            })
            continue
        if source_transcript not in source_protein_ids:
            reason = "no_extracted_protein"
        elif not passing_by_source.get(source_transcript):
            reason = "no_passing_rbh"
        elif source_transcript not in selected_by_source:
            reason = "transcript_pair_conflict"
        else:
            reason = "gene_pair_conflict"
        entries.append(
            _unscored_mapping(source_transcript, "transcript", reason)
        )

    entries.sort(key=lambda entry: (
        MAPPING_FEATURE_TYPES.index(entry["feature_type"]),
        entry["source_id"],
    ))
    hit_rows = []
    retained_keys = {
        (hit.source_id, hit.target_id) for hit in retained_hits
    }
    for hit in sorted(hits, key=lambda value: (
        value.source_id, value.target_id,
    )):
        key = (hit.source_id, hit.target_id)
        hit_rows.append({
            "source_id": hit.source_id,
            "target_id": hit.target_id,
            **hit.metrics(),
            "filter_reason": filter_reasons[key] or "passed",
            "transcript_pair_selected": key in transcript_selected_keys,
            "retained_mapping": key in retained_keys,
        })

    counts = {
        "source_genes": len(source_hierarchy.genes),
        "source_transcripts": len(source_hierarchy.transcripts),
        "source_proteins": len(source_protein_ids),
        "source_genes_excluded_ambiguous": len(
            source_hierarchy.excluded_genes
        ),
        "source_transcripts_excluded_ambiguous": len(
            source_hierarchy.excluded_transcripts
        ),
        "target_genes": len(target_hierarchy.genes),
        "target_transcripts": len(target_hierarchy.transcripts),
        "target_proteins": len(target_protein_ids),
        "target_genes_excluded_ambiguous": len(
            target_hierarchy.excluded_genes
        ),
        "target_transcripts_excluded_ambiguous": len(
            target_hierarchy.excluded_transcripts
        ),
        "rbh_hits_raw": len(hits),
        "rbh_hits_passing": len(passing),
        "transcript_pairs_resolved": len(transcript_selected),
        "gene_pair_candidates": len(gene_evidence),
        "gene_groups_retained": len(retained_genes),
        "transcript_groups_retained": len(retained_hits),
        "gene_groups_unscored": sum(
            entry["feature_type"] == "gene"
            and entry["status"] == "unscored"
            for entry in entries
        ),
        "transcript_groups_unscored": sum(
            entry["feature_type"] == "transcript"
            and entry["status"] == "unscored"
            for entry in entries
        ),
    }
    document = {
        "schema_version": SCHEMA_VERSION,
        "method": METHOD,
        "metadata": {
            "scope": "full",
            "parameters": dict(parameters),
            "provenance": dict(provenance),
            "counts": counts,
        },
        "mappings": entries,
    }
    validate_mapping_against_hierarchies(
        document,
        source_hierarchy=source_hierarchy,
        target_hierarchy=target_hierarchy,
    )
    return document, hit_rows


def validate_mapping_document(document: Any) -> dict[str, Any]:
    if not isinstance(document, Mapping):
        raise ScopeBuildError("ortholog mapping document is not an object")
    if (
        document.get("schema_version") != SCHEMA_VERSION
        or document.get("method") != METHOD
        or not isinstance(document.get("metadata"), Mapping)
        or not isinstance(document.get("mappings"), list)
    ):
        raise ScopeBuildError("ortholog mapping header is unsupported")
    seen_source = set()
    seen_truth: dict[str, set[str]] = {
        feature_type: set() for feature_type in MAPPING_FEATURE_TYPES
    }
    for index, entry in enumerate(document["mappings"], start=1):
        if not isinstance(entry, Mapping):
            raise ScopeBuildError(f"mapping entry {index} is not an object")
        source_id = entry.get("source_id")
        feature_type = entry.get("feature_type")
        status = entry.get("status")
        truth_ids = entry.get("truth_ids")
        if (
            not isinstance(source_id, str)
            or not source_id
            or feature_type not in MAPPING_FEATURE_TYPES
            or status not in MAPPING_STATUSES
            or not isinstance(truth_ids, list)
            or any(not isinstance(value, str) or not value for value in truth_ids)
        ):
            raise ScopeBuildError(f"mapping entry {index} is malformed")
        key = (feature_type, source_id)
        if key in seen_source:
            raise ScopeBuildError(
                f"duplicate {feature_type} source mapping {source_id!r}"
            )
        seen_source.add(key)
        if status == "retained":
            if len(truth_ids) != 1:
                raise ScopeBuildError(
                    f"retained {feature_type} {source_id!r} must have one truth ID"
                )
            target_id = truth_ids[0]
            if target_id in seen_truth[feature_type]:
                raise ScopeBuildError(
                    f"target {feature_type} {target_id!r} is not one-to-one"
                )
            seen_truth[feature_type].add(target_id)
        elif truth_ids or not isinstance(entry.get("reason"), str):
            raise ScopeBuildError(
                f"unscored {feature_type} {source_id!r} is malformed"
            )
    return dict(document)


def _mapping_scope_counts(document: Mapping[str, Any]) -> dict[str, int]:
    return {
        "gene_groups_retained": sum(
            entry["feature_type"] == "gene"
            and entry["status"] == "retained"
            for entry in document["mappings"]
        ),
        "transcript_groups_retained": sum(
            entry["feature_type"] == "transcript"
            and entry["status"] == "retained"
            for entry in document["mappings"]
        ),
        "gene_groups_unscored": sum(
            entry["feature_type"] == "gene"
            and entry["status"] == "unscored"
            for entry in document["mappings"]
        ),
        "transcript_groups_unscored": sum(
            entry["feature_type"] == "transcript"
            and entry["status"] == "unscored"
            for entry in document["mappings"]
        ),
    }


def validate_mapping_against_hierarchies(
    document: Mapping[str, Any],
    *,
    source_hierarchy: Hierarchy,
    target_hierarchy: Hierarchy,
    require_complete_source_scope: bool = True,
) -> dict[str, Any]:
    """Validate a protein-RBH map against its exact annotation scope.

    Structural validation alone cannot detect a well-formed but fabricated or
    stale identifier.  This check binds every retained *and* unscored source
    row to the declared source hierarchy, every retained truth identifier to
    the target hierarchy, and verifies the redundant count metadata.
    """

    validated = validate_mapping_document(document)
    source_ids = {
        "gene": set(source_hierarchy.genes),
        "transcript": set(source_hierarchy.transcripts),
    }
    target_ids = {
        "gene": set(target_hierarchy.genes),
        "transcript": set(target_hierarchy.transcripts),
    }
    observed_source = {
        feature_type: set() for feature_type in MAPPING_FEATURE_TYPES
    }
    retained_genes: dict[str, str] = {}
    for entry in validated["mappings"]:
        feature_type = entry["feature_type"]
        source_id = entry["source_id"]
        if source_id not in source_ids[feature_type]:
            raise ScopeBuildError(
                f"unknown source {feature_type} ID {source_id!r} in ortholog map"
            )
        observed_source[feature_type].add(source_id)
        if entry["status"] != "retained":
            continue
        target_id = entry["truth_ids"][0]
        if target_id not in target_ids[feature_type]:
            raise ScopeBuildError(
                f"unknown target {feature_type} ID {target_id!r} in ortholog map"
            )
        if feature_type == "gene":
            retained_genes[source_id] = target_id

    if require_complete_source_scope:
        for feature_type in MAPPING_FEATURE_TYPES:
            missing = source_ids[feature_type] - observed_source[feature_type]
            if missing:
                raise ScopeBuildError(
                    f"ortholog map omits source {feature_type} ID "
                    f"{sorted(missing)[0]!r}"
                )

    for entry in validated["mappings"]:
        if (
            entry["feature_type"] != "transcript"
            or entry["status"] != "retained"
        ):
            continue
        source_gene = source_hierarchy.transcript_to_gene[entry["source_id"]]
        target_gene = target_hierarchy.transcript_to_gene[
            entry["truth_ids"][0]
        ]
        if retained_genes.get(source_gene) != target_gene:
            raise ScopeBuildError(
                f"retained transcript {entry['source_id']!r} maps outside its "
                "retained gene pair"
            )

    counts = validated["metadata"].get("counts")
    if not isinstance(counts, Mapping):
        raise ScopeBuildError("ortholog mapping count metadata is missing")
    expected_counts = {
        "source_genes": len(source_ids["gene"]),
        "source_transcripts": len(source_ids["transcript"]),
        "target_genes": len(target_ids["gene"]),
        "target_transcripts": len(target_ids["transcript"]),
        **_mapping_scope_counts(validated),
    }
    for name, expected in expected_counts.items():
        if name.endswith("_unscored") and name not in counts:
            # Early subset maps did not publish redundant unscored counts.
            continue
        observed = counts.get(name)
        if (
            not isinstance(observed, int)
            or isinstance(observed, bool)
            or observed != expected
        ):
            raise ScopeBuildError(
                f"ortholog mapping count {name} is stale: "
                f"declared={observed!r}, observed={expected}"
            )
    return validated


def _scope_annotation_records(
    document: Mapping[str, Any],
) -> Mapping[str, Any]:
    metadata = document.get("metadata")
    if not isinstance(metadata, Mapping):
        raise ScopeBuildError("ortholog mapping metadata is missing")
    scope = metadata.get("scope")
    if scope == "full":
        provenance = metadata.get("provenance")
        records = (
            provenance.get("inputs")
            if isinstance(provenance, Mapping)
            else None
        )
    elif scope == "subset":
        records = metadata.get("subset_inputs")
    else:
        raise ScopeBuildError(f"unsupported ortholog mapping scope {scope!r}")
    if not isinstance(records, Mapping):
        raise ScopeBuildError(
            f"{scope} ortholog mapping has no annotation input provenance"
        )
    return records


def _verify_annotation_fingerprint(
    record: Any,
    path: Path,
    *,
    label: str,
) -> None:
    if not isinstance(record, Mapping):
        raise ScopeBuildError(f"ortholog mapping lacks {label} fingerprint")
    observed = _input_record(path)
    if (
        record.get("size") != observed["size"]
        or record.get("sha256") != observed["sha256"]
    ):
        raise ScopeBuildError(
            f"ortholog mapping {label} fingerprint does not match {path}"
        )


def validate_mapping_against_annotations(
    document: Mapping[str, Any],
    *,
    source_annotation: str | Path,
    target_annotation: str | Path,
    require_complete_source_scope: bool = True,
    verify_input_fingerprints: bool = True,
) -> dict[str, Any]:
    """Validate map membership and, by default, its frozen input hashes."""

    source_path = Path(source_annotation).expanduser().resolve()
    target_path = Path(target_annotation).expanduser().resolve()
    validated = validate_mapping_against_hierarchies(
        document,
        source_hierarchy=parse_hierarchy(source_path),
        target_hierarchy=parse_hierarchy(target_path),
        require_complete_source_scope=require_complete_source_scope,
    )
    if verify_input_fingerprints:
        records = _scope_annotation_records(validated)
        _verify_annotation_fingerprint(
            records.get("source_annotation"),
            source_path,
            label="source annotation",
        )
        _verify_annotation_fingerprint(
            records.get("target_annotation"),
            target_path,
            label="target annotation",
        )
    return validated


def validate_mapping_from_recorded_annotations(
    document: Mapping[str, Any],
) -> dict[str, Any]:
    """Resolve frozen annotation paths and fully validate a map artifact."""

    validated = validate_mapping_document(document)
    records = _scope_annotation_records(validated)
    paths = {}
    for name in ("source_annotation", "target_annotation"):
        record = records.get(name)
        if not isinstance(record, Mapping) or not isinstance(
            record.get("path"), str,
        ):
            raise ScopeBuildError(
                f"ortholog mapping lacks recorded {name.replace('_', ' ')} path"
            )
        paths[name] = Path(record["path"]).expanduser().resolve()
    return validate_mapping_against_annotations(
        validated,
        source_annotation=paths["source_annotation"],
        target_annotation=paths["target_annotation"],
    )


def _enforce_minimum_groups(
    document: Mapping[str, Any],
    *,
    minimum_gene_groups: int,
    minimum_transcript_groups: int,
) -> None:
    retained = {
        feature_type: sum(
            entry["feature_type"] == feature_type
            and entry["status"] == "retained"
            for entry in document["mappings"]
        )
        for feature_type in MAPPING_FEATURE_TYPES
    }
    if retained["gene"] < minimum_gene_groups:
        raise ScopeBuildError(
            f"retained gene groups {retained['gene']} are below required "
            f"minimum {minimum_gene_groups}"
        )
    if retained["transcript"] < minimum_transcript_groups:
        raise ScopeBuildError(
            f"retained transcript groups {retained['transcript']} are below "
            f"required minimum {minimum_transcript_groups}"
        )


def _write_hit_table(path: Path, rows: Sequence[Mapping[str, Any]]) -> None:
    fields = (
        "source_id",
        "target_id",
        "identity",
        "alignment_length",
        "query_coverage",
        "target_coverage",
        "evalue",
        "bits",
        "filter_reason",
        "transcript_pair_selected",
        "retained_mapping",
    )
    lines = ["\t".join(fields)]
    for row in rows:
        lines.append("\t".join(str(row[field]) for field in fields))
    _write_bytes(path, ("\n".join(lines) + "\n").encode("utf-8"))


def filter_scope_document(
    full_document: Mapping[str, Any],
    *,
    source_subset_annotation: str | Path,
    target_subset_annotation: str | Path,
    minimum_gene_groups: int = 1,
    minimum_transcript_groups: int = 1,
) -> dict[str, Any]:
    """Restrict a full map to the exact source/target annotation search space."""

    validate_mapping_from_recorded_annotations(full_document)
    if minimum_gene_groups < 0 or minimum_transcript_groups < 0:
        raise ValueError("minimum retained groups must be non-negative")
    source_path = Path(source_subset_annotation).expanduser().resolve()
    target_path = Path(target_subset_annotation).expanduser().resolve()
    source = parse_hierarchy(source_path)
    target = parse_hierarchy(target_path)
    index = {
        (entry["feature_type"], entry["source_id"]): entry
        for entry in full_document["mappings"]
    }
    allowed_targets = {
        "gene": set(target.genes),
        "transcript": set(target.transcripts),
    }
    entries = []
    for feature_type, source_ids in (
        ("gene", source.genes),
        ("transcript", source.transcripts),
    ):
        for source_id in source_ids:
            original = index.get((feature_type, source_id))
            if original is None:
                entries.append(_unscored_mapping(
                    source_id,
                    feature_type,
                    "source_not_in_full_scope",
                ))
                continue
            truth_ids = original.get("truth_ids", [])
            if (
                original.get("status") == "retained"
                and len(truth_ids) == 1
                and truth_ids[0] in allowed_targets[feature_type]
            ):
                entries.append(dict(original))
            elif original.get("status") == "retained":
                entries.append(_unscored_mapping(
                    source_id,
                    feature_type,
                    "target_outside_subset",
                ))
            else:
                entries.append(dict(original))
    entries.sort(key=lambda entry: (
        MAPPING_FEATURE_TYPES.index(entry["feature_type"]),
        entry["source_id"],
    ))
    retained_counts = {
        f"{feature_type}_groups_retained": sum(
            entry["feature_type"] == feature_type
            and entry["status"] == "retained"
            for entry in entries
        )
        for feature_type in MAPPING_FEATURE_TYPES
    }
    document = {
        "schema_version": SCHEMA_VERSION,
        "method": METHOD,
        "metadata": {
            "scope": "subset",
            "parent_mapping_sha256": canonical_sha256(full_document),
            "parent_parameters": dict(
                full_document["metadata"].get("parameters", {})
            ),
            "parent_provenance": dict(
                full_document["metadata"].get("provenance", {})
            ),
            "subset_inputs": {
                "source_annotation": _input_record(source_path),
                "target_annotation": _input_record(target_path),
            },
            "counts": {
                "source_genes": len(source.genes),
                "source_transcripts": len(source.transcripts),
                "target_genes": len(target.genes),
                "target_transcripts": len(target.transcripts),
                **retained_counts,
            },
            "parameters": {
                "minimum_gene_groups": minimum_gene_groups,
                "minimum_transcript_groups": minimum_transcript_groups,
            },
        },
        "mappings": entries,
    }
    validate_mapping_against_hierarchies(
        document,
        source_hierarchy=source,
        target_hierarchy=target,
    )
    _enforce_minimum_groups(
        document,
        minimum_gene_groups=minimum_gene_groups,
        minimum_transcript_groups=minimum_transcript_groups,
    )
    return document


def write_filtered_scope(
    full_mapping: str | Path,
    source_subset_annotation: str | Path,
    target_subset_annotation: str | Path,
    output: str | Path,
    *,
    minimum_gene_groups: int = 1,
    minimum_transcript_groups: int = 1,
) -> dict[str, Any]:
    full_path = Path(full_mapping).expanduser().resolve()
    try:
        document = json.loads(full_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise ScopeBuildError(
            f"cannot read full ortholog mapping {full_path}: {exc}"
        ) from exc
    filtered = filter_scope_document(
        document,
        source_subset_annotation=source_subset_annotation,
        target_subset_annotation=target_subset_annotation,
        minimum_gene_groups=minimum_gene_groups,
        minimum_transcript_groups=minimum_transcript_groups,
    )
    _atomic_write_json(Path(output), filtered)
    return filtered


def _common_inputs(
    source_annotation: str | Path,
    source_genome: str | Path,
    target_annotation: str | Path,
    target_genome: str | Path,
) -> tuple[dict[str, Path], dict[str, Any]]:
    paths = {
        "source_annotation": Path(source_annotation).expanduser().resolve(),
        "source_genome": Path(source_genome).expanduser().resolve(),
        "target_annotation": Path(target_annotation).expanduser().resolve(),
        "target_genome": Path(target_genome).expanduser().resolve(),
    }
    return paths, {
        name: _input_record(path) for name, path in paths.items()
    }


def _dry_run_document(
    *,
    action: str,
    paths: Mapping[str, Path],
    input_records: Mapping[str, Any],
    output_dir: Path,
    gffread: str | Path,
    mmseqs: str | Path | None,
    parameters: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    gffread_path = _resolve_executable(gffread)
    commands = []
    for label in ("source", "target"):
        commands.append({
            "argv": [
                str(gffread_path),
                "-y",
                f"$OUTPUT/{label}.proteins.raw.fa",
                "-g",
                str(paths[f"{label}_genome"]),
                str(paths[f"{label}_annotation"]),
            ],
            "shell": False,
        })
    tools: dict[str, Any] = {
        "gffread": {"executable": str(gffread_path)},
    }
    if mmseqs is not None:
        mmseqs_path = _resolve_executable(mmseqs)
        tools["mmseqs"] = {"executable": str(mmseqs_path)}
        commands.append({
            "argv": [
                str(mmseqs_path),
                "easy-rbh",
                "$OUTPUT/source.proteins.fa",
                "$OUTPUT/target.proteins.fa",
                "$OUTPUT/mmseqs.easy-rbh.tsv",
                "$OUTPUT/mmseqs-tmp",
                "--format-output",
                MMSEQS_FORMAT,
                "--threads",
                str(parameters["threads"]),
                "--min-seq-id",
                str(parameters["minimum_identity"]),
                "-c",
                str(parameters["minimum_reciprocal_coverage"]),
                "--cov-mode",
                "0",
                "-e",
                str(parameters["maximum_evalue"]),
            ],
            "shell": False,
        })
    return {
        "schema_version": SCHEMA_VERSION,
        "action": action,
        "dry_run": True,
        "output_dir": str(Path(output_dir).resolve()),
        "inputs": dict(input_records),
        "tools": tools,
        "parameters": dict(parameters or {}),
        "commands": commands,
    }


def prepare_proteins(
    source_annotation: str | Path,
    source_genome: str | Path,
    target_annotation: str | Path,
    target_genome: str | Path,
    output_dir: str | Path,
    *,
    gffread: str | Path = "gffread",
    runner: Runner = subprocess.run,
) -> dict[str, Any]:
    paths, input_records = _common_inputs(
        source_annotation,
        source_genome,
        target_annotation,
        target_genome,
    )
    gffread_tool = _tool_record(
        gffread,
        version_arguments=("--version",),
        runner=runner,
    )
    result: dict[str, Any] = {}

    def build(work: Path) -> None:
        prepared = _prepare_in_workspace(
            work,
            source_annotation=paths["source_annotation"],
            source_genome=paths["source_genome"],
            target_annotation=paths["target_annotation"],
            target_genome=paths["target_genome"],
            gffread_tool=gffread_tool,
            runner=runner,
        )
        outputs = {
            name: _artifact_record(work / name)
            for name in (
                "source.proteins.fa",
                "target.proteins.fa",
                "source.hierarchy.json",
                "target.hierarchy.json",
            )
        }
        _assert_frozen_records(input_records)
        _assert_frozen_tools({"gffread": gffread_tool})
        manifest = {
            "schema_version": SCHEMA_VERSION,
            "method": PREPARE_METHOD,
            "inputs": input_records,
            "tools": {"gffread": gffread_tool},
            "commands": prepared["commands"],
            "counts": {
                label: _protein_summary(prepared["protein_stats"][label])
                for label in ("source", "target")
            },
            "outputs": outputs,
        }
        _write_json(work / "prepare.manifest.json", manifest)
        result.update(manifest)

    _publish_directory(Path(output_dir), build)
    return result


def build_ortholog_scope(
    source_annotation: str | Path,
    source_genome: str | Path,
    target_annotation: str | Path,
    target_genome: str | Path,
    output_dir: str | Path,
    *,
    gffread: str | Path = "gffread",
    mmseqs: str | Path = "mmseqs",
    minimum_identity: float = 0.30,
    minimum_reciprocal_coverage: float = 0.50,
    maximum_evalue: float = 1e-5,
    threads: int = 1,
    minimum_gene_groups: int = DEFAULT_MINIMUM_RETAINED_GROUPS,
    minimum_transcript_groups: int = DEFAULT_MINIMUM_RETAINED_GROUPS,
    source_subset_annotation: str | Path | None = None,
    target_subset_annotation: str | Path | None = None,
    runner: Runner = subprocess.run,
) -> dict[str, Any]:
    """Run gffread and MMseqs2 and atomically publish an ortholog scope."""

    if (source_subset_annotation is None) != (
        target_subset_annotation is None
    ):
        raise ValueError(
            "source and target subset annotations must be provided together"
        )
    parameters = _parameter_document(
        minimum_identity=minimum_identity,
        minimum_reciprocal_coverage=minimum_reciprocal_coverage,
        maximum_evalue=maximum_evalue,
        threads=threads,
        minimum_gene_groups=minimum_gene_groups,
        minimum_transcript_groups=minimum_transcript_groups,
    )
    paths, input_records = _common_inputs(
        source_annotation,
        source_genome,
        target_annotation,
        target_genome,
    )
    if source_subset_annotation is not None:
        input_records = {
            **input_records,
            "source_subset_annotation": _input_record(
                source_subset_annotation
            ),
            "target_subset_annotation": _input_record(
                target_subset_annotation
            ),
        }
    gffread_tool = _tool_record(
        gffread,
        version_arguments=("--version",),
        runner=runner,
    )
    mmseqs_tool = _tool_record(
        mmseqs,
        version_arguments=("version",),
        runner=runner,
    )
    result: dict[str, Any] = {}

    def build(work: Path) -> None:
        prepared = _prepare_in_workspace(
            work,
            source_annotation=paths["source_annotation"],
            source_genome=paths["source_genome"],
            target_annotation=paths["target_annotation"],
            target_genome=paths["target_genome"],
            gffread_tool=gffread_tool,
            runner=runner,
        )
        source_ids = set(prepared["protein_stats"]["source"]["ids"])
        target_ids = set(prepared["protein_stats"]["target"]["ids"])
        raw_hits = work / "mmseqs.easy-rbh.tsv"
        mmseqs_tmp = work / "mmseqs-tmp"
        mmseqs_argv = [
            mmseqs_tool["executable"],
            "easy-rbh",
            work / "source.proteins.fa",
            work / "target.proteins.fa",
            raw_hits,
            mmseqs_tmp,
            "--format-output",
            MMSEQS_FORMAT,
            "--threads",
            str(parameters["threads"]),
            "--min-seq-id",
            str(parameters["minimum_identity"]),
            "-c",
            str(parameters["minimum_reciprocal_coverage"]),
            "--cov-mode",
            "0",
            "-e",
            str(parameters["maximum_evalue"]),
        ]
        _run_checked(
            mmseqs_argv,
            runner=runner,
            cwd=work,
            label="MMseqs easy-rbh",
        )
        hits = parse_mmseqs_hits(
            raw_hits,
            source_proteins=source_ids,
            target_proteins=target_ids,
        )
        provenance = {
            "inputs": input_records,
            "tools": {
                "gffread": gffread_tool,
                "mmseqs": mmseqs_tool,
            },
            "commands": [
                *prepared["commands"],
                {
                    "label": "protein_reciprocal_best_hits",
                    "argv": _logical_argv(mmseqs_argv, work),
                    "shell": False,
                    "format_output": MMSEQS_FORMAT,
                },
            ],
            "normalized_proteins": {
                "source": _artifact_record(work / "source.proteins.fa"),
                "target": _artifact_record(work / "target.proteins.fa"),
            },
        }
        full_document, hit_rows = build_mapping_document(
            hits,
            source_hierarchy=prepared["hierarchies"]["source"],
            target_hierarchy=prepared["hierarchies"]["target"],
            source_protein_ids=source_ids,
            target_protein_ids=target_ids,
            parameters=parameters,
            provenance=provenance,
        )
        validate_mapping_against_annotations(
            full_document,
            source_annotation=paths["source_annotation"],
            target_annotation=paths["target_annotation"],
        )
        _assert_frozen_records(input_records)
        _assert_frozen_tools({
            "gffread": gffread_tool,
            "mmseqs": mmseqs_tool,
        })
        _write_hit_table(work / "rbh.hits.tsv", hit_rows)
        if source_subset_annotation is not None:
            _write_json(work / "ortholog_map.full.json", full_document)
            final_document = filter_scope_document(
                full_document,
                source_subset_annotation=source_subset_annotation,
                target_subset_annotation=target_subset_annotation,
                minimum_gene_groups=minimum_gene_groups,
                minimum_transcript_groups=minimum_transcript_groups,
            )
        else:
            final_document = full_document
            _enforce_minimum_groups(
                final_document,
                minimum_gene_groups=minimum_gene_groups,
                minimum_transcript_groups=minimum_transcript_groups,
            )
        _write_json(work / "ortholog_map.json", final_document)
        raw_hits.unlink()
        shutil.rmtree(mmseqs_tmp, ignore_errors=True)
        output_names = [
            "source.proteins.fa",
            "target.proteins.fa",
            "source.hierarchy.json",
            "target.hierarchy.json",
            "rbh.hits.tsv",
            "ortholog_map.json",
        ]
        if (work / "ortholog_map.full.json").is_file():
            output_names.append("ortholog_map.full.json")
        manifest = {
            "schema_version": SCHEMA_VERSION,
            "method": METHOD,
            "parameters": parameters,
            "inputs": input_records,
            "tools": {
                "gffread": gffread_tool,
                "mmseqs": mmseqs_tool,
            },
            "commands": provenance["commands"],
            "counts": {
                "proteins": {
                    label: _protein_summary(
                        prepared["protein_stats"][label]
                    )
                    for label in ("source", "target")
                },
                "mapping": final_document["metadata"]["counts"],
            },
            "outputs": {
                name: _artifact_record(work / name)
                for name in output_names
            },
        }
        _write_json(work / "ortholog_scope.manifest.json", manifest)
        result.update({
            "mapping": final_document,
            "manifest": manifest,
        })

    _publish_directory(Path(output_dir), build)
    return result


def validate_scope_bundle(bundle: str | Path) -> dict[str, Any]:
    """Validate every immutable artifact in one published full-scope bundle."""

    requested_path = Path(bundle).expanduser()
    if requested_path.is_symlink():
        raise ScopeBuildError(
            f"ortholog scope bundle must not be a symlink: {requested_path}"
        )
    bundle_path = requested_path.resolve()
    if bundle_path.name == "ortholog_map.json":
        bundle_path = bundle_path.parent
    if not bundle_path.is_dir():
        raise ScopeBuildError(f"ortholog scope bundle is not a directory: {bundle_path}")
    mapping_path = bundle_path / "ortholog_map.json"
    scope_manifest_path = bundle_path / "ortholog_scope.manifest.json"
    try:
        mapping_document = json.loads(
            mapping_path.read_text(encoding="utf-8")
        )
        scope_manifest = json.loads(
            scope_manifest_path.read_text(encoding="utf-8")
        )
    except (OSError, json.JSONDecodeError) as exc:
        raise ScopeBuildError(
            f"cannot read ortholog scope bundle {bundle_path}: {exc}"
        ) from exc
    validated = validate_mapping_from_recorded_annotations(mapping_document)
    if validated["metadata"].get("scope") != "full":
        raise ScopeBuildError("registry mapping must have full scope")
    manifest_outputs = (
        scope_manifest.get("outputs")
        if isinstance(scope_manifest, Mapping)
        else None
    )
    manifest_counts = (
        scope_manifest.get("counts")
        if isinstance(scope_manifest, Mapping)
        else None
    )
    mapping_provenance = validated["metadata"].get("provenance")
    required_outputs = {
        "source.proteins.fa",
        "target.proteins.fa",
        "source.hierarchy.json",
        "target.hierarchy.json",
        "rbh.hits.tsv",
        "ortholog_map.json",
    }
    if (
        not isinstance(scope_manifest, Mapping)
        or scope_manifest.get("schema_version") != SCHEMA_VERSION
        or scope_manifest.get("method") != METHOD
        or scope_manifest.get("parameters")
        != validated["metadata"].get("parameters")
        or not isinstance(mapping_provenance, Mapping)
        or scope_manifest.get("inputs") != mapping_provenance.get("inputs")
        or scope_manifest.get("tools") != mapping_provenance.get("tools")
        or scope_manifest.get("commands") != mapping_provenance.get("commands")
        or not isinstance(manifest_counts, Mapping)
        or manifest_counts.get("mapping") != validated["metadata"].get("counts")
        or not isinstance(manifest_outputs, Mapping)
        or not required_outputs.issubset(manifest_outputs)
    ):
        raise ScopeBuildError("ortholog scope manifest is stale or tampered")
    for name, record in sorted(manifest_outputs.items()):
        artifact = bundle_path / name
        if (
            not isinstance(name, str)
            or not name
            or Path(name).name != name
            or not isinstance(record, Mapping)
            or record.get("name") != name
            or not artifact.is_file()
            or artifact.is_symlink()
            or artifact.stat().st_size <= 0
            or record.get("size") != artifact.stat().st_size
            or record.get("sha256") != sha256_file(artifact)
        ):
            raise ScopeBuildError(
                "ortholog scope manifest is stale or tampered: "
                f"artifact {name!r}"
            )
    if scope_manifest_path.is_symlink():
        raise ScopeBuildError("ortholog scope manifest must not be a symlink")
    allowed_paths = {
        scope_manifest_path,
        *(bundle_path / name for name in manifest_outputs),
    }
    unexpected = sorted(
        str(path.relative_to(bundle_path))
        for path in bundle_path.rglob("*")
        if path not in allowed_paths
    )
    if unexpected:
        raise ScopeBuildError(
            f"ortholog scope bundle has unexpected artifact {unexpected[0]!r}"
        )
    return {
        "bundle_path": bundle_path,
        "mapping_path": mapping_path,
        "mapping_sha256": sha256_file(mapping_path),
        "mapping": validated,
        "manifest": scope_manifest,
    }


def _registry_manifest_document(manifest: str | Path) -> dict[str, Any]:
    manifest_path = Path(manifest).expanduser().resolve()
    try:
        document = json.loads(manifest_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise ScopeBuildError(
            f"cannot read canonical dataset manifest {manifest_path}: {exc}"
        ) from exc
    if not isinstance(document, Mapping):
        raise ScopeBuildError("canonical dataset manifest is not an object")
    return dict(document)


def _biological_scenario_ids(manifest_document: Mapping[str, Any]) -> tuple[str, ...]:
    scenarios = (
        manifest_document.get("scenarios")
        if isinstance(manifest_document, Mapping)
        else None
    )
    if not isinstance(scenarios, list):
        raise ScopeBuildError("canonical dataset manifest has no scenarios")
    biological_ids = []
    for scenario in scenarios:
        if not isinstance(scenario, Mapping) or scenario.get("kind") != "biological":
            continue
        scenario_id = scenario.get("id")
        if not isinstance(scenario_id, str) or not scenario_id:
            raise ScopeBuildError("biological scenario has no valid ID")
        biological_ids.append(scenario_id)
    if not biological_ids or len(set(biological_ids)) != len(biological_ids):
        raise ScopeBuildError(
            "canonical dataset manifest biological scenario IDs are empty or duplicate"
        )
    return tuple(sorted(biological_ids))


def finalize_mapping_registry(
    manifest: str | Path,
    mapping_root: str | Path,
    output: str | Path,
) -> dict[str, Any]:
    """Publish a deterministic, manifest-bound registry of validated maps."""

    manifest_document = _registry_manifest_document(manifest)
    mapping_root = Path(mapping_root).expanduser().resolve()
    output_path = Path(output).expanduser().resolve()
    expected = set(_biological_scenario_ids(manifest_document))
    discovered = {
        path.parent.name
        for path in mapping_root.glob("*/ortholog_map.json")
        if path.is_file()
    }
    missing = expected - discovered
    unknown = discovered - expected
    if missing:
        raise ScopeBuildError(
            f"ortholog mapping root is missing scenario {sorted(missing)[0]!r}"
        )
    if unknown:
        raise ScopeBuildError(
            f"ortholog mapping root has unexpected scenario {sorted(unknown)[0]!r}"
        )

    mappings = {}
    for scenario_id in sorted(expected):
        bundle = validate_scope_bundle(mapping_root / scenario_id)
        mapping_path = bundle["mapping_path"]
        relative_path = Path(os.path.relpath(
            mapping_path,
            start=output_path.parent,
        )).as_posix()
        mappings[scenario_id] = {
            "id_policy": "ortholog-map",
            "path": relative_path,
            "sha256": bundle["mapping_sha256"],
        }
    registry = {
        "schema_version": REGISTRY_SCHEMA_VERSION,
        "manifest_sha256": canonical_sha256(manifest_document),
        "mappings": mappings,
    }
    _atomic_write_json(output_path, registry)
    return registry


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="action", required=True)

    def add_inputs(command: argparse.ArgumentParser) -> None:
        command.add_argument("source_annotation", type=Path)
        command.add_argument("source_genome", type=Path)
        command.add_argument("target_annotation", type=Path)
        command.add_argument("target_genome", type=Path)
        command.add_argument("output_dir", type=Path)
        command.add_argument("--gffread", default="gffread")
        command.add_argument("--dry-run", action="store_true")

    prepare = subparsers.add_parser(
        "prepare",
        help="extract and normalize transcript proteins only",
    )
    add_inputs(prepare)

    build = subparsers.add_parser(
        "build",
        help="extract proteins, run MMseqs easy-rbh, and publish a map",
    )
    add_inputs(build)
    build.add_argument("--mmseqs", default="mmseqs")
    build.add_argument("--minimum-identity", type=float, default=0.30)
    build.add_argument(
        "--minimum-reciprocal-coverage",
        type=float,
        default=0.50,
    )
    build.add_argument("--maximum-evalue", type=float, default=1e-5)
    build.add_argument("--threads", type=int, default=1)
    build.add_argument(
        "--minimum-gene-groups",
        type=int,
        default=DEFAULT_MINIMUM_RETAINED_GROUPS,
    )
    build.add_argument(
        "--minimum-transcript-groups",
        type=int,
        default=DEFAULT_MINIMUM_RETAINED_GROUPS,
    )
    build.add_argument("--source-subset-annotation", type=Path)
    build.add_argument("--target-subset-annotation", type=Path)

    filter_parser = subparsers.add_parser(
        "filter",
        help="restrict a full map to source and target subset annotations",
    )
    filter_parser.add_argument("full_mapping", type=Path)
    filter_parser.add_argument("source_subset_annotation", type=Path)
    filter_parser.add_argument("target_subset_annotation", type=Path)
    filter_parser.add_argument("output", type=Path)
    filter_parser.add_argument("--minimum-gene-groups", type=int, default=1)
    filter_parser.add_argument(
        "--minimum-transcript-groups",
        type=int,
        default=1,
    )

    finalize = subparsers.add_parser(
        "finalize-registry",
        help="validate completed full maps and publish their registry",
    )
    finalize.add_argument("manifest", type=Path)
    finalize.add_argument("mapping_root", type=Path)
    finalize.add_argument("output", type=Path)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    arguments = parser.parse_args(argv)
    try:
        if arguments.action == "finalize-registry":
            finalize_mapping_registry(
                arguments.manifest,
                arguments.mapping_root,
                arguments.output,
            )
            return 0
        if arguments.action == "filter":
            write_filtered_scope(
                arguments.full_mapping,
                arguments.source_subset_annotation,
                arguments.target_subset_annotation,
                arguments.output,
                minimum_gene_groups=arguments.minimum_gene_groups,
                minimum_transcript_groups=(
                    arguments.minimum_transcript_groups
                ),
            )
            return 0
        paths, inputs = _common_inputs(
            arguments.source_annotation,
            arguments.source_genome,
            arguments.target_annotation,
            arguments.target_genome,
        )
        parameters = None
        mmseqs = None
        if arguments.action == "build":
            parameters = _parameter_document(
                minimum_identity=arguments.minimum_identity,
                minimum_reciprocal_coverage=(
                    arguments.minimum_reciprocal_coverage
                ),
                maximum_evalue=arguments.maximum_evalue,
                threads=arguments.threads,
                minimum_gene_groups=arguments.minimum_gene_groups,
                minimum_transcript_groups=(
                    arguments.minimum_transcript_groups
                ),
            )
            mmseqs = arguments.mmseqs
            if (
                (arguments.source_subset_annotation is None)
                != (arguments.target_subset_annotation is None)
            ):
                raise ValueError(
                    "source and target subset annotations must be provided "
                    "together"
                )
            if arguments.source_subset_annotation is not None:
                inputs = {
                    **inputs,
                    "source_subset_annotation": _input_record(
                        arguments.source_subset_annotation
                    ),
                    "target_subset_annotation": _input_record(
                        arguments.target_subset_annotation
                    ),
                }
        if arguments.dry_run:
            print(json.dumps(
                _dry_run_document(
                    action=arguments.action,
                    paths=paths,
                    input_records=inputs,
                    output_dir=arguments.output_dir,
                    gffread=arguments.gffread,
                    mmseqs=mmseqs,
                    parameters=parameters,
                ),
                indent=2,
                sort_keys=True,
            ))
            return 0
        if arguments.action == "prepare":
            prepare_proteins(
                arguments.source_annotation,
                arguments.source_genome,
                arguments.target_annotation,
                arguments.target_genome,
                arguments.output_dir,
                gffread=arguments.gffread,
            )
        else:
            build_ortholog_scope(
                arguments.source_annotation,
                arguments.source_genome,
                arguments.target_annotation,
                arguments.target_genome,
                arguments.output_dir,
                gffread=arguments.gffread,
                mmseqs=arguments.mmseqs,
                minimum_identity=arguments.minimum_identity,
                minimum_reciprocal_coverage=(
                    arguments.minimum_reciprocal_coverage
                ),
                maximum_evalue=arguments.maximum_evalue,
                threads=arguments.threads,
                minimum_gene_groups=arguments.minimum_gene_groups,
                minimum_transcript_groups=(
                    arguments.minimum_transcript_groups
                ),
                source_subset_annotation=(
                    arguments.source_subset_annotation
                ),
                target_subset_annotation=(
                    arguments.target_subset_annotation
                ),
            )
        return 0
    except (FileExistsError, OSError, ScopeBuildError, ValueError) as exc:
        print(f"ortholog-scope: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
