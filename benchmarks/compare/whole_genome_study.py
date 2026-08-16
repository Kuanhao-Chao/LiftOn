#!/usr/bin/env python3
"""Acquire and prepare the v1.0.11 biology-first whole-genome study.

The module is intentionally separate from the sealed release campaign. It
pins exact target assemblies and released annotations, stores downloaded
bytes content-addressably, validates the honey-bee chromosome-name transform,
builds uniform protein-RBH scopes, and exports a runtime-only dataset registry.
"""
from __future__ import annotations

import argparse
import concurrent.futures
import datetime as dt
import gzip
import hashlib
import json
import os
import shutil
import stat
import subprocess
import tempfile
import urllib.parse
import urllib.request
import zipfile
from pathlib import Path, PurePosixPath
from typing import Any, BinaryIO, Mapping, Sequence

from lifton.gff3_validator import GENE_TYPES, TRANSCRIPT_TYPES

from . import ortholog_scope


HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parent.parent
DEFAULT_STUDY = (
    REPO_ROOT / "benchmarks" / "manifests" /
    "lifton_v111_biology_study.json"
)
DEFAULT_BENCHMARK_REGISTRY = HERE / "benchmarks.json"
LOCK_NAME = "annotation-lock.json"
PREFLIGHT_NAME = "study-preflight.json"
DATASET_OVERLAY = "runtime/datasets.json"
MODEL_VIEWS_DIR = "runtime/model-views"
SCHEMA_VERSION = 1
EXPECTED_PAIR_IDS = (
    "bee",
    "drosophila",
    "t2_human_to_gorilla",
    "t2_mouse_to_caroli",
    "t3_dog_to_cat",
    "t3_human_to_macaque",
    "t3_human_to_marmoset",
)
FASTA_CHARACTERS = frozenset("ACGTURYSWKMBDHVN.-")
MODEL_FEATURE_TYPES = GENE_TYPES | TRANSCRIPT_TYPES | {"exon", "CDS"}


class StudyError(RuntimeError):
    """Raised when study acquisition or preparation cannot fail closed."""


def utc_now() -> str:
    return dt.datetime.now(dt.timezone.utc).isoformat().replace("+00:00", "Z")


def canonical_sha256(value: Any) -> str:
    payload = json.dumps(
        value,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=True,
        allow_nan=False,
    ).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def sha256_file(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def md5_file(path: str | Path) -> str:
    digest = hashlib.md5(usedforsecurity=False)
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def file_record(path: str | Path) -> dict[str, Any]:
    resolved = Path(path).resolve()
    if not resolved.is_file() or resolved.stat().st_size <= 0:
        raise StudyError(f"required file is missing or empty: {resolved}")
    return {
        "path": str(resolved),
        "size": resolved.stat().st_size,
        "sha256": sha256_file(resolved),
    }


def _atomic_json(path: Path, document: Any) -> None:
    path = path.resolve()
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists() or path.is_symlink():
        raise FileExistsError(f"refusing to replace immutable study file: {path}")
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    try:
        with temporary.open("x", encoding="utf-8") as handle:
            json.dump(document, handle, indent=2, sort_keys=True)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    except Exception:
        if temporary.exists():
            temporary.unlink()
        raise


def load_study(path: str | Path = DEFAULT_STUDY) -> dict[str, Any]:
    study_path = Path(path).resolve()
    try:
        document = json.loads(study_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise StudyError(f"cannot read study manifest {study_path}: {exc}") from exc
    if not isinstance(document, Mapping) or document.get("schema_version") != 1:
        raise StudyError("study manifest schema_version must be 1")
    pairs = document.get("pairs")
    if (
        not isinstance(pairs, list)
        or tuple(pair.get("id") for pair in pairs if isinstance(pair, Mapping))
        != EXPECTED_PAIR_IDS
    ):
        raise StudyError("study pair IDs/order do not match the approved cohort")
    if document.get("candidate", {}).get("sha") != (
        "c623f0bddc5b5051a3670a3b0064c52cc61bb719"
    ):
        raise StudyError("study candidate SHA is not exact v1.0.11")
    if document.get("reference", {}).get("sha") != (
        "e503643d8346c600fedabcd3a4dff5c0873a4a37"
    ):
        raise StudyError("study reference SHA is not exact v1.0.8")
    for pair in pairs:
        annotation = pair.get("target_annotation")
        if not isinstance(annotation, Mapping):
            raise StudyError(f"{pair.get('id')}: target annotation is missing")
        required = {
            "provider", "assembly_accession", "assembly_name", "release",
            "release_date", "reported_gene_count", "transport",
        }
        if not required.issubset(annotation):
            raise StudyError(
                f"{pair.get('id')}: target annotation metadata is incomplete"
            )
    return dict(document)


def _safe_cache_root(path: str | Path) -> Path:
    root = Path(path).expanduser().resolve()
    if root in {Path("/"), Path.home().resolve(), REPO_ROOT.resolve()}:
        raise StudyError(f"unsafe study cache root: {root}")
    root.mkdir(parents=True, exist_ok=True)
    return root


def validate_zip_archive(path: str | Path) -> list[zipfile.ZipInfo]:
    archive = Path(path)
    try:
        handle = zipfile.ZipFile(archive)
    except (OSError, zipfile.BadZipFile) as exc:
        raise StudyError(f"unsafe or unreadable ZIP archive {archive}: {exc}") from exc
    with handle:
        members = handle.infolist()
        if not members:
            raise StudyError(f"ZIP archive is empty: {archive}")
        for member in members:
            pure = PurePosixPath(member.filename)
            mode = member.external_attr >> 16
            if (
                pure.is_absolute()
                or ".." in pure.parts
                or "" in pure.parts
                or stat.S_ISLNK(mode)
            ):
                raise StudyError(
                    f"unsafe ZIP member in {archive}: {member.filename!r}"
                )
        return members


def _object_path(cache_root: Path, digest: str, name: str) -> Path:
    return cache_root / "objects" / "sha256" / digest[:2] / digest / name


def _publish_file(cache_root: Path, source: Path, name: str) -> dict[str, Any]:
    digest = sha256_file(source)
    destination = _object_path(cache_root, digest, name)
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists():
        if sha256_file(destination) != digest:
            raise StudyError(f"content-addressed object is corrupted: {destination}")
        source.unlink(missing_ok=True)
    else:
        os.replace(source, destination)
    return {
        "path": str(destination.relative_to(cache_root)),
        "size": destination.stat().st_size,
        "sha256": digest,
    }


def _publish_stream(
    cache_root: Path,
    source: BinaryIO,
    *,
    name: str,
) -> dict[str, Any]:
    staging_root = cache_root / "staging"
    staging_root.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        prefix="extract-", suffix=".part", dir=staging_root,
    )
    temporary = Path(temporary_name)
    digest = hashlib.sha256()
    try:
        with os.fdopen(descriptor, "wb") as destination:
            while True:
                block = source.read(1024 * 1024)
                if not block:
                    break
                destination.write(block)
                digest.update(block)
            destination.flush()
            os.fsync(destination.fileno())
        return _publish_file(cache_root, temporary, name)
    except Exception:
        temporary.unlink(missing_ok=True)
        raise


def _resolved_object(cache_root: Path, record: Mapping[str, Any]) -> Path:
    relative = Path(str(record.get("path", "")))
    if relative.is_absolute() or ".." in relative.parts:
        raise StudyError("annotation lock contains an unsafe object path")
    path = (cache_root / relative).resolve()
    if not path.is_relative_to(cache_root):
        raise StudyError("annotation lock object escapes the cache root")
    observed = file_record(path)
    if (
        observed["size"] != record.get("size")
        or observed["sha256"] != record.get("sha256")
    ):
        raise StudyError(f"annotation lock object changed: {path}")
    return path


def _download_url(url: str, destination: Path) -> None:
    request = urllib.request.Request(url, headers={"User-Agent": "LiftOn-study/1"})
    with urllib.request.urlopen(request, timeout=120) as response, destination.open(
        "xb"
    ) as handle:
        shutil.copyfileobj(response, handle, length=1024 * 1024)
        handle.flush()
        os.fsync(handle.fileno())


def _download_ncbi(accession: str, destination: Path, datasets: Path) -> None:
    command = [
        str(datasets),
        "download", "genome", "accession", accession,
        "--include", "gff3,protein",
        "--filename", str(destination),
        "--no-progressbar",
    ]
    result = subprocess.run(command, text=True, capture_output=True, check=False)
    if result.returncode != 0:
        raise StudyError(
            f"NCBI Datasets failed for {accession}: "
            f"{result.stderr[-2000:].strip()}"
        )


def _select_member(
    members: Sequence[zipfile.ZipInfo],
    *,
    suffixes: Sequence[str],
    label: str,
) -> zipfile.ZipInfo:
    selected = [
        member for member in members
        if not member.is_dir()
        and any(member.filename.endswith(suffix) for suffix in suffixes)
    ]
    if len(selected) != 1:
        raise StudyError(
            f"{label}: expected one ZIP member matching {suffixes}, "
            f"observed {[member.filename for member in selected]}"
        )
    return selected[0]


def _extract_ncbi_package(
    cache_root: Path,
    archive: Path,
    accession: str,
) -> dict[str, Any]:
    members = validate_zip_archive(archive)
    gff = _select_member(
        members,
        suffixes=(
            "/genomic.gff", "/genomic.gff.gz",
            "_genomic.gff", "_genomic.gff.gz",
        ),
        label=f"{accession} annotation",
    )
    protein = _select_member(
        members,
        suffixes=(
            "/protein.faa", "/protein.faa.gz",
            "_protein.faa", "_protein.faa.gz",
        ),
        label=f"{accession} proteins",
    )
    report = _select_member(
        members,
        suffixes=("assembly_data_report.jsonl",),
        label=f"{accession} assembly report",
    )
    with zipfile.ZipFile(archive) as handle:
        records = {}
        for role, member in (
            ("truth_gff", gff),
            ("target_proteins", protein),
            ("assembly_report", report),
        ):
            name = Path(member.filename).name
            with handle.open(member) as source:
                records[role] = _publish_stream(
                    cache_root,
                    source,
                    name=f"{accession}-{name}",
                )
            records[role]["archive_member"] = member.filename
    return records


def _benchmark_entries(path: str | Path) -> dict[str, Mapping[str, Any]]:
    document = json.loads(Path(path).read_text(encoding="utf-8"))
    entries = {
        row.get("id"): row
        for row in document.get("benchmarks", [])
        if isinstance(row, Mapping) and isinstance(row.get("id"), str)
    }
    missing = set(EXPECTED_PAIR_IDS) - set(entries)
    if missing:
        raise StudyError(f"benchmark registry is missing study IDs: {sorted(missing)}")
    return entries


def _open_text(path: Path):
    return gzip.open(path, "rt", encoding="utf-8", errors="strict") \
        if path.name.endswith(".gz") else path.open(
            "r", encoding="utf-8", errors="strict",
        )


def fasta_lengths(path: str | Path) -> dict[str, int]:
    fasta = Path(path)
    lengths: dict[str, int] = {}
    identifier = None
    length = 0
    with _open_text(fasta) as handle:
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
                    raise StudyError(
                        f"{fasta}: invalid FASTA ID at line {line_number}"
                    )
                if identifier in lengths:
                    raise StudyError(
                        f"{fasta}: duplicate FASTA ID {identifier!r}"
                    )
                length = 0
            elif identifier is None:
                raise StudyError(f"{fasta}: sequence precedes first FASTA header")
            else:
                length += len(line)
    if identifier is not None:
        if identifier in lengths:
            raise StudyError(f"{fasta}: duplicate FASTA ID {identifier!r}")
        lengths[identifier] = length
    if not lengths or any(value <= 0 for value in lengths.values()):
        raise StudyError(f"{fasta}: FASTA has no non-empty records")
    return lengths


def normalize_bee_annotation(
    raw_gff: str | Path,
    target_fasta: str | Path,
    output: str | Path,
    seqid_map: Mapping[str, str],
    expected_lengths: Mapping[str, int],
) -> dict[str, Any]:
    """Rewrite only column 1 after exact chromosome-length validation."""

    raw_path = Path(raw_gff).resolve()
    target_path = Path(target_fasta).resolve()
    output_path = Path(output).resolve()
    target_lengths = fasta_lengths(target_path)
    if len(set(seqid_map.values())) != len(seqid_map):
        raise StudyError("bee sequence-ID mapping is not one-to-one")
    for source_id, target_id in seqid_map.items():
        expected = expected_lengths.get(source_id)
        if expected is None or target_lengths.get(target_id) != expected:
            raise StudyError(
                f"bee chromosome length mismatch: {source_id} -> {target_id}"
            )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    directives: dict[str, int] = {}
    embedded_lengths: dict[str, int] = {}
    embedded_id: str | None = None
    feature_counts: dict[str, int] = {}
    rows = 0
    with _open_text(raw_path) as source, output_path.open(
        "x", encoding="utf-8", newline="\n",
    ) as destination:
        destination.write("##gff-version 3\n")
        for source_id, target_id in seqid_map.items():
            destination.write(
                f"##sequence-region {target_id} 1 {expected_lengths[source_id]}\n"
            )
        for line_number, raw in enumerate(source, start=1):
            if raw.startswith("##FASTA"):
                continue
            if raw.startswith(">"):
                embedded_id = raw[1:].strip().split()[0]
                if (
                    embedded_id not in seqid_map
                    or embedded_id in embedded_lengths
                ):
                    raise StudyError(
                        f"{raw_path}: invalid embedded FASTA ID "
                        f"{embedded_id!r} at line {line_number}"
                    )
                embedded_lengths[embedded_id] = 0
                continue
            if embedded_id is not None and "\t" not in raw:
                sequence = raw.strip().upper()
                if sequence and set(sequence) <= FASTA_CHARACTERS:
                    embedded_lengths[embedded_id] += len(sequence)
                    continue
                if sequence:
                    raise StudyError(
                        f"{raw_path}: invalid embedded FASTA sequence at "
                        f"line {line_number}"
                    )
            if not raw.strip():
                continue
            if raw.startswith("##sequence-region"):
                fields = raw.split()
                if len(fields) == 4 and fields[1] in seqid_map:
                    directives[fields[1]] = int(fields[3])
                continue
            if raw.startswith("#"):
                continue
            columns = raw.rstrip("\r\n").split("\t")
            if len(columns) != 9:
                raise StudyError(
                    f"{raw_path}: line {line_number} has {len(columns)} columns"
                )
            source_id = columns[0]
            embedded_id = None
            if source_id not in seqid_map:
                raise StudyError(
                    f"{raw_path}: line {line_number} has unmapped seqid {source_id!r}"
                )
            if columns[8] not in {"", "."} and "=" not in columns[8]:
                raise StudyError(
                    f"{raw_path}: line {line_number} is not GFF3"
                )
            columns[0] = seqid_map[source_id]
            destination.write("\t".join(columns) + "\n")
            rows += 1
            feature_counts[columns[2]] = feature_counts.get(columns[2], 0) + 1
        destination.flush()
        os.fsync(destination.fileno())
    if rows == 0:
        raise StudyError("bee annotation normalization emitted no features")
    for source_id, observed in directives.items():
        if observed != expected_lengths[source_id]:
            raise StudyError(f"bee source directive length mismatch: {source_id}")
    if embedded_lengths and embedded_lengths != dict(expected_lengths):
        raise StudyError("bee embedded FASTA lengths do not match the target assembly")
    return {
        "method": "bee-gff3-seqid-column1-normalization-v1",
        "input": file_record(raw_path),
        "target_fasta": file_record(target_path),
        "seqid_map": dict(seqid_map),
        "expected_lengths": dict(expected_lengths),
        "source_sequence_directives": directives,
        "embedded_fasta_lengths": dict(sorted(embedded_lengths.items())),
        "feature_rows": rows,
        "feature_counts": dict(sorted(feature_counts.items())),
        "output": file_record(output_path),
    }


def build_model_view(
    annotation: str | Path,
    output: str | Path,
) -> dict[str, Any]:
    """Create a lossless gene-model-only view of a released annotation."""

    source_path = Path(annotation).resolve()
    output_path = Path(output).resolve()
    hierarchy = ortholog_scope.parse_hierarchy(source_path)
    excluded_genes = set(hierarchy.excluded_genes)
    excluded_transcripts = set(hierarchy.excluded_transcripts)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    feature_counts: dict[str, int] = {}
    excluded_feature_counts: dict[str, int] = {}
    sequence_regions = 0
    with _open_text(source_path) as source, output_path.open(
        "x", encoding="utf-8", newline="\n",
    ) as destination:
        destination.write("##gff-version 3\n")
        for line_number, raw in enumerate(source, start=1):
            if raw.startswith("##FASTA") or raw.startswith(">"):
                break
            if raw.startswith("##sequence-region"):
                destination.write(raw.rstrip("\r\n") + "\n")
                sequence_regions += 1
                continue
            if not raw.strip() or raw.startswith("#"):
                continue
            columns = raw.rstrip("\r\n").split("\t")
            if len(columns) != 9:
                raise StudyError(
                    f"{source_path}: line {line_number} has "
                    f"{len(columns)} columns"
                )
            feature_type = columns[2]
            if feature_type not in MODEL_FEATURE_TYPES:
                continue
            if columns[8] not in {"", "."} and "=" not in columns[8]:
                raise StudyError(
                    f"{source_path}: line {line_number} is not GFF3"
                )
            attributes = columns[8].split(";")
            values = {}
            for index, item in enumerate(attributes):
                if "=" not in item:
                    continue
                name, raw_value = item.split("=", 1)
                values[name] = (
                    index,
                    tuple(
                        urllib.parse.unquote(value, errors="strict")
                        for value in raw_value.split(",")
                        if value
                    ),
                    tuple(value for value in raw_value.split(",") if value),
                )
            excluded = False
            if feature_type in GENE_TYPES:
                identifiers = values.get("ID", (-1, (), ()))[1]
                excluded = len(identifiers) == 1 and (
                    identifiers[0] in excluded_genes
                )
            elif feature_type in TRANSCRIPT_TYPES:
                identifiers = values.get("ID", (-1, (), ()))[1]
                parents = values.get("Parent", (-1, (), ()))[1]
                excluded = (
                    len(identifiers) == 1
                    and identifiers[0] in excluded_transcripts
                ) or any(parent in excluded_genes for parent in parents)
            else:
                parent_record = values.get("Parent")
                if parent_record is not None:
                    retained = [
                        raw_value
                        for decoded, raw_value in zip(
                            parent_record[1], parent_record[2]
                        )
                        if decoded not in excluded_transcripts
                    ]
                    if not retained:
                        excluded = True
                    elif len(retained) != len(parent_record[1]):
                        attributes[parent_record[0]] = (
                            "Parent=" + ",".join(retained)
                        )
                        columns[8] = ";".join(attributes)
            if excluded:
                excluded_feature_counts[feature_type] = (
                    excluded_feature_counts.get(feature_type, 0) + 1
                )
                continue
            destination.write("\t".join(columns) + "\n")
            feature_counts[feature_type] = (
                feature_counts.get(feature_type, 0) + 1
            )
        destination.flush()
        os.fsync(destination.fileno())
    if not any(
        feature_counts.get(feature_type, 0)
        for feature_type in GENE_TYPES
    ) or not any(
        feature_counts.get(feature_type, 0)
        for feature_type in TRANSCRIPT_TYPES
    ):
        raise StudyError(f"{source_path}: model view has no gene hierarchy")
    return {
        "method": "gene-model-only-gff3-view-v1",
        "input": file_record(source_path),
        "output": file_record(output_path),
        "feature_counts": dict(sorted(feature_counts.items())),
        "excluded_feature_counts": dict(
            sorted(excluded_feature_counts.items())
        ),
        "excluded_model_ids": {
            "genes": len(excluded_genes),
            "gene_ids_sha256": canonical_sha256(sorted(excluded_genes)),
            "transcripts": len(excluded_transcripts),
            "transcript_ids_sha256": canonical_sha256(
                sorted(excluded_transcripts)
            ),
        },
        "sequence_regions": sequence_regions,
    }


def _annotation_info(report_path: Path) -> Mapping[str, Any]:
    first = report_path.read_text(encoding="utf-8").splitlines()[0]
    report = json.loads(first)
    info = report.get("annotationInfo")
    if not isinstance(info, Mapping):
        raise StudyError(f"assembly report has no annotationInfo: {report_path}")
    return info


def _verify_annotation_metadata(
    pair: Mapping[str, Any],
    report_path: Path,
) -> dict[str, Any]:
    expected = pair["target_annotation"]
    info = _annotation_info(report_path)
    observed_release = str(info.get("name", ""))
    pinned_release = str(expected["release"])
    release_matches = (
        observed_release == pinned_release
        or (
            "Annotation Release " in pinned_release
            and pinned_release.split("Annotation Release ", 1)[1]
            in observed_release
        )
    )
    gene_count = (
        info.get("stats", {}).get("geneCounts", {}).get("total")
        if isinstance(info.get("stats"), Mapping)
        else None
    )
    if (
        not release_matches
        or info.get("releaseDate") != expected["release_date"]
        or int(gene_count or -1) != int(expected["reported_gene_count"])
    ):
        raise StudyError(
            f"{pair['id']}: NCBI annotation metadata disagrees with the pin"
        )
    return {
        "name": observed_release,
        "provider": info.get("provider"),
        "release_date": info.get("releaseDate"),
        "gene_counts": info.get("stats", {}).get("geneCounts"),
        "method": info.get("method"),
        "pipeline": info.get("pipeline"),
        "software_version": info.get("softwareVersion"),
        "status": info.get("status"),
        "report_url": info.get("reportUrl"),
    }


def acquire_annotations(
    study_path: str | Path,
    cache_root: str | Path,
    *,
    benchmark_registry: str | Path = DEFAULT_BENCHMARK_REGISTRY,
    datasets: str | Path = "datasets",
) -> dict[str, Any]:
    study = load_study(study_path)
    root = _safe_cache_root(cache_root)
    lock_path = root / LOCK_NAME
    if lock_path.exists():
        return verify_annotation_lock(study_path, lock_path, root)
    entries = _benchmark_entries(benchmark_registry)
    staging = root / "staging"
    staging.mkdir(parents=True, exist_ok=True)
    targets = {}
    for pair in study["pairs"]:
        pair_id = pair["id"]
        annotation = pair["target_annotation"]
        archive = staging / f"{pair_id}-{annotation['assembly_accession']}.zip"
        if annotation["transport"] == "ncbi_datasets":
            if archive.exists():
                try:
                    validate_zip_archive(archive)
                except StudyError:
                    archive.unlink()
            if not archive.exists():
                _download_ncbi(
                    annotation["assembly_accession"], archive, Path(datasets),
                )
            extracted = _extract_ncbi_package(
                root, archive, annotation["assembly_accession"],
            )
            report_path = _resolved_object(root, extracted["assembly_report"])
            metadata = _verify_annotation_metadata(pair, report_path)
            archive_record = _publish_file(
                root, archive, f"{annotation['assembly_accession']}.zip",
            )
            targets[pair_id] = {
                "identity": dict(annotation),
                "archive": archive_record,
                "artifacts": extracted,
                "observed_annotation": metadata,
            }
            continue
        if annotation["transport"] != "figshare":
            raise StudyError(f"{pair_id}: unsupported annotation transport")
        if archive.exists() and (
            archive.stat().st_size != annotation["expected_size"]
            or md5_file(archive) != annotation["expected_md5"]
        ):
            archive.unlink()
        if not archive.exists():
            _download_url(annotation["download_url"], archive)
        if (
            archive.stat().st_size != annotation["expected_size"]
            or md5_file(archive) != annotation["expected_md5"]
        ):
            raise StudyError("bee Figshare archive size or MD5 does not match")
        members = validate_zip_archive(archive)
        member = _select_member(
            members,
            suffixes=(".gff", ".gff3"),
            label="bee published annotation",
        )
        with zipfile.ZipFile(archive) as handle, handle.open(member) as source:
            raw_record = _publish_stream(
                root, source, name="bee-published-raw.gff",
            )
        raw_path = _resolved_object(root, raw_record)
        normalized_temp = staging / "bee-normalized.gff3"
        normalized_temp.unlink(missing_ok=True)
        normalization = normalize_bee_annotation(
            raw_path,
            entries[pair_id]["tgt_genome"],
            normalized_temp,
            annotation["chromosome_seqid_map"],
            annotation["chromosome_lengths"],
        )
        normalized_record = _publish_file(
            root, normalized_temp, "bee-published-normalized.gff3",
        )
        normalization["output"] = normalized_record
        archive_record = _publish_file(root, archive, annotation["filename"])
        targets[pair_id] = {
            "identity": dict(annotation),
            "archive": archive_record,
            "artifacts": {
                "raw_truth_gff": raw_record,
                "truth_gff": normalized_record,
            },
            "normalization": normalization,
        }
    lock = {
        "schema_version": SCHEMA_VERSION,
        "kind": "lifton-v1.0.11-biology-study-annotation-lock",
        "created_at": utc_now(),
        "study": file_record(study_path),
        "study_sha256": canonical_sha256(study),
        "benchmark_registry": file_record(benchmark_registry),
        "targets": targets,
    }
    lock["fingerprint"] = canonical_sha256(lock)
    _atomic_json(lock_path, lock)
    return verify_annotation_lock(study_path, lock_path, root)


def verify_annotation_lock(
    study_path: str | Path,
    lock_path: str | Path,
    cache_root: str | Path,
) -> dict[str, Any]:
    study = load_study(study_path)
    root = _safe_cache_root(cache_root)
    try:
        lock = json.loads(Path(lock_path).read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise StudyError(f"cannot read annotation lock: {exc}") from exc
    fingerprint = lock.get("fingerprint") if isinstance(lock, Mapping) else None
    material = dict(lock) if isinstance(lock, Mapping) else {}
    material.pop("fingerprint", None)
    if (
        lock.get("schema_version") != SCHEMA_VERSION
        or lock.get("kind") != "lifton-v1.0.11-biology-study-annotation-lock"
        or lock.get("study_sha256") != canonical_sha256(study)
        or fingerprint != canonical_sha256(material)
        or set(lock.get("targets", {})) != set(EXPECTED_PAIR_IDS)
    ):
        raise StudyError("annotation lock schema, study binding, or fingerprint is invalid")
    if lock.get("study") != file_record(study_path):
        raise StudyError("annotation lock study file record changed")
    registry_record = lock.get("benchmark_registry")
    if not isinstance(registry_record, Mapping):
        raise StudyError("annotation lock benchmark registry record is missing")
    if registry_record != file_record(registry_record.get("path", "")):
        raise StudyError("annotation lock benchmark registry changed")
    verified_files = 0
    for pair_id, target in lock["targets"].items():
        if target.get("identity") != next(
            pair["target_annotation"] for pair in study["pairs"]
            if pair["id"] == pair_id
        ):
            raise StudyError(f"annotation identity changed in lock: {pair_id}")
        _resolved_object(root, target["archive"])
        verified_files += 1
        for record in target.get("artifacts", {}).values():
            _resolved_object(root, record)
            verified_files += 1
    return {
        "verified": True,
        "study_id": study["study_id"],
        "targets": len(lock["targets"]),
        "files": verified_files,
        "lock": file_record(lock_path),
        "document": lock,
    }


def _truth_path(cache_root: Path, target: Mapping[str, Any]) -> Path:
    return _resolved_object(cache_root, target["artifacts"]["truth_gff"])


def prepare_study(
    study_path: str | Path,
    cache_root: str | Path,
    *,
    benchmark_registry: str | Path = DEFAULT_BENCHMARK_REGISTRY,
    gffread: str | Path = "gffread",
    mmseqs: str | Path = "mmseqs",
    threads: int = 8,
    max_active: int = 2,
) -> dict[str, Any]:
    study = load_study(study_path)
    root = _safe_cache_root(cache_root)
    verification = verify_annotation_lock(
        study_path, root / LOCK_NAME, root,
    )
    lock = verification["document"]
    if lock["benchmark_registry"] != file_record(benchmark_registry):
        raise StudyError(
            "configured benchmark registry differs from the annotation lock"
        )
    entries = _benchmark_entries(benchmark_registry)
    views: dict[str, dict[str, Any]] = {}
    staging = root / "staging"
    staging.mkdir(parents=True, exist_ok=True)
    view_records_root = root / MODEL_VIEWS_DIR
    view_records_root.mkdir(parents=True, exist_ok=True)

    def derive_views(pair: Mapping[str, Any]) -> tuple[str, dict[str, Any]]:
        pair_id = pair["id"]
        benchmark = entries[pair_id]
        paths = {
            "source": Path(benchmark["ref_gff"]),
            "target": _truth_path(root, lock["targets"][pair_id]),
        }
        result = {}
        for role, source_path in paths.items():
            record_path = view_records_root / pair_id / f"{role}.json"
            if record_path.exists():
                audit = json.loads(record_path.read_text(encoding="utf-8"))
                if (
                    audit.get("method")
                    != "gene-model-only-gff3-view-v1"
                    or audit.get("input") != file_record(source_path)
                ):
                    raise StudyError(
                        f"{pair_id}: cached {role} model view changed"
                    )
                _resolved_object(root, audit["output"])
                result[role] = audit
                continue
            print(
                f"[{utc_now()}] building {role} model view: {pair_id}",
                file=os.sys.stderr,
                flush=True,
            )
            temporary = staging / f"{pair_id}-{role}-model-view.gff3"
            temporary.unlink(missing_ok=True)
            audit = build_model_view(source_path, temporary)
            artifact = _publish_file(
                root, temporary, f"{pair_id}-{role}-model-view.gff3",
            )
            audit["output"] = artifact
            _atomic_json(record_path, audit)
            result[role] = audit
        return pair_id, result

    with concurrent.futures.ThreadPoolExecutor(max_workers=max_active) as executor:
        futures = {
            executor.submit(derive_views, pair): pair["id"]
            for pair in study["pairs"]
        }
        for future in concurrent.futures.as_completed(futures):
            try:
                pair_id, record = future.result()
            except Exception as exc:
                failed = futures[future]
                for pending in futures:
                    pending.cancel()
                print(
                    f"[{utc_now()}] model-view preparation failed: "
                    f"{failed}: {exc}",
                    file=os.sys.stderr,
                    flush=True,
                )
                raise StudyError(
                    f"{failed}: model-view preparation failed: {exc}"
                ) from exc
            views[pair_id] = record

    scopes_root = root / "runtime" / "ortholog-scopes"
    scopes_root.mkdir(parents=True, exist_ok=True)

    def build(pair: Mapping[str, Any]) -> tuple[str, dict[str, Any]]:
        pair_id = pair["id"]
        destination = scopes_root / pair_id
        print(
            f"[{utc_now()}] preparing ortholog scope: {pair_id}",
            file=os.sys.stderr,
            flush=True,
        )
        if destination.exists():
            bundle = ortholog_scope.validate_scope_bundle(destination)
        else:
            benchmark = entries[pair_id]
            ortholog_scope.build_ortholog_scope(
                _resolved_object(root, views[pair_id]["source"]["output"]),
                benchmark["ref_genome"],
                _resolved_object(root, views[pair_id]["target"]["output"]),
                benchmark["tgt_genome"],
                destination,
                gffread=gffread,
                mmseqs=mmseqs,
                threads=threads,
                minimum_gene_groups=100,
                minimum_transcript_groups=100,
            )
            bundle = ortholog_scope.validate_scope_bundle(destination)
        print(
            f"[{utc_now()}] ortholog scope ready: {pair_id}",
            file=os.sys.stderr,
            flush=True,
        )
        return pair_id, {
            "mapping": file_record(bundle["mapping_path"]),
            "manifest": file_record(
                destination / "ortholog_scope.manifest.json"
            ),
            "counts": bundle["mapping"]["metadata"]["counts"],
        }

    scopes = {}
    with concurrent.futures.ThreadPoolExecutor(max_workers=max_active) as executor:
        futures = {executor.submit(build, pair): pair["id"] for pair in study["pairs"]}
        for future in concurrent.futures.as_completed(futures):
            try:
                pair_id, record = future.result()
            except Exception as exc:
                failed = futures[future]
                for pending in futures:
                    pending.cancel()
                print(
                    f"[{utc_now()}] ortholog-scope preparation failed: "
                    f"{failed}: {exc}",
                    file=os.sys.stderr,
                    flush=True,
                )
                raise StudyError(
                    f"{failed}: ortholog-scope preparation failed: {exc}"
                ) from exc
            scopes[pair_id] = record

    datasets = []
    for pair in study["pairs"]:
        pair_id = pair["id"]
        benchmark = entries[pair_id]
        datasets.append({
            "id": pair_id,
            "species": pair["public_label"],
            "reference_fa": Path(benchmark["ref_genome"]).resolve().as_uri(),
            "target_fa": Path(benchmark["tgt_genome"]).resolve().as_uri(),
            "reference_gff": Path(benchmark["ref_gff"]).resolve().as_uri(),
            "approx_size_gb": round(
                (
                    Path(benchmark["ref_genome"]).stat().st_size
                    + Path(benchmark["tgt_genome"]).stat().st_size
                ) / 1_000_000_000,
                3,
            ),
            "cross_species": bool(benchmark.get("cross_species")),
            "annotation_database": benchmark.get("annotation_database", "RefSeq"),
            "truth_gff": _resolved_object(
                root, views[pair_id]["target"]["output"],
            ).as_uri(),
            "truth_source_gff": _resolved_object(
                root, views[pair_id]["source"]["output"],
            ).as_uri(),
            "ortholog_map": Path(scopes[pair_id]["mapping"]["path"]).as_uri(),
            "truth_id_policy": "ortholog-map",
        })
    overlay = {
        "_study": {
            "id": study["study_id"],
            "manifest": file_record(study_path),
            "annotation_lock": file_record(root / LOCK_NAME),
        },
        "datasets": datasets,
        "lifton_flags": [
            "--stream", "--inmemory-liftoff", "--locus-pipeline",
            "-t", "8", "--native", "-copies",
        ],
        "evaluation_flags": [],
    }
    overlay_path = root / DATASET_OVERLAY
    if overlay_path.exists():
        observed = json.loads(overlay_path.read_text(encoding="utf-8"))
        if observed != overlay:
            raise StudyError("runtime dataset registry exists with different content")
    else:
        _atomic_json(overlay_path, overlay)
    preflight = {
        "schema_version": SCHEMA_VERSION,
        "kind": "lifton-v1.0.11-biology-study-preflight",
        "created_at": utc_now(),
        "campaign_ready": True,
        "study": file_record(study_path),
        "annotation_lock": file_record(root / LOCK_NAME),
        "benchmark_registry": file_record(benchmark_registry),
        "dataset_registry": file_record(overlay_path),
        "ortholog_scopes": dict(sorted(scopes.items())),
        "model_views": dict(sorted(views.items())),
        "blockers": [],
    }
    preflight["fingerprint"] = canonical_sha256(preflight)
    preflight_path = root / PREFLIGHT_NAME
    if preflight_path.exists():
        existing = json.loads(preflight_path.read_text(encoding="utf-8"))
        comparable = dict(existing)
        comparable["created_at"] = preflight["created_at"]
        comparable["fingerprint"] = preflight["fingerprint"]
        if comparable != preflight:
            raise StudyError("study preflight exists with different evidence")
        preflight = existing
    else:
        _atomic_json(preflight_path, preflight)
    return preflight


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--study", type=Path, default=DEFAULT_STUDY)
    parser.add_argument("--cache-root", type=Path, required=True)
    parser.add_argument(
        "--benchmark-registry", type=Path, default=DEFAULT_BENCHMARK_REGISTRY,
    )
    subparsers = parser.add_subparsers(dest="action", required=True)
    acquire = subparsers.add_parser("acquire")
    acquire.add_argument("--datasets", default="datasets")
    subparsers.add_parser("verify")
    prepare = subparsers.add_parser("prepare")
    prepare.add_argument("--gffread", default="gffread")
    prepare.add_argument("--mmseqs", default="mmseqs")
    prepare.add_argument("--threads", type=int, default=8)
    prepare.add_argument("--max-active", type=int, default=2)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    arguments = build_parser().parse_args(argv)
    try:
        if arguments.action == "acquire":
            result = acquire_annotations(
                arguments.study,
                arguments.cache_root,
                benchmark_registry=arguments.benchmark_registry,
                datasets=arguments.datasets,
            )
        elif arguments.action == "verify":
            result = verify_annotation_lock(
                arguments.study,
                arguments.cache_root / LOCK_NAME,
                arguments.cache_root,
            )
        else:
            result = prepare_study(
                arguments.study,
                arguments.cache_root,
                benchmark_registry=arguments.benchmark_registry,
                gffread=arguments.gffread,
                mmseqs=arguments.mmseqs,
                threads=arguments.threads,
                max_active=arguments.max_active,
            )
        print(json.dumps(result, indent=2, sort_keys=True, default=str))
        return 0
    except (FileExistsError, OSError, StudyError, ValueError) as exc:
        print(f"whole-genome-study: {exc}", file=os.sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
