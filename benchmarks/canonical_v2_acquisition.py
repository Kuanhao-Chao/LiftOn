#!/usr/bin/env python3
"""Acquire and prepare the canonical-v2 benchmark inputs.

The manifest fixes source identities. This worker turns the first observed
bytes into an explicit, content-addressed lock only after the operator
acknowledges sources whose upstream provider does not publish a checksum.
Downloads are executed as argument vectors, never through a shell.

Materialization is a separate fail-closed step. It decompresses runtime
assets, creates verified FASTA indexes, generates synthetic truth, validates
curated ortholog maps, and exports merged registry overlays. Missing mappings
or incompatible sequence identifiers are reported as blockers; no runnable
truth-required registry row is emitted in that state.
"""
from __future__ import annotations

import argparse
import gzip
import json
import os
import shutil
import stat
import subprocess
import sys
import tempfile
import zipfile
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path, PurePosixPath
from typing import Any, Callable, Mapping, Sequence

from benchmarks import manifest_tools


HERE = Path(__file__).resolve().parent
DEFAULT_CACHE_ROOT = HERE / "data" / "canonical-v2-cache"
DEFAULT_BENCHMARK_REGISTRY = HERE / "compare" / "benchmarks.json"
DEFAULT_DATASET_REGISTRY = HERE / "datasets.json"
DEFAULT_WORK_ROOT = HERE / "compare" / "work"
LOCK_NAME = "acquisition.lock.json"
PREFLIGHT_NAME = "preflight.json"
ORTHOLOG_SCHEMA_VERSION = 1
MINIMUM_FULL_TRUTH_GROUPS = 100
REMOTE_PENDING_PIN = "identity_pinned_bytes_pending"
FASTA_ROLES = {
    "genome",
    "reference_genome",
    "target_genome",
}
BIOLOGICAL_IDS = tuple(manifest_tools.EXPECTED_SCENARIO_IDS[:8])


class AcquisitionError(RuntimeError):
    """Raised when acquisition or materialization cannot safely continue."""


CommandRunner = Callable[[Sequence[str]], int]
SyntheticBuilder = Callable[
    [Mapping[str, Any], Path, Path, Path],
    Mapping[str, Path],
]


def _atomic_write_json(path: Path, document: Mapping[str, Any]) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{path.name}.",
        suffix=".tmp",
        dir=path.parent,
    )
    temporary = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as handle:
            json.dump(document, handle, indent=2, sort_keys=True)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    finally:
        if temporary.exists():
            temporary.unlink()


def _atomic_write_bytes(path: Path, payload: bytes) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{path.name}.",
        suffix=".tmp",
        dir=path.parent,
    )
    temporary = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "wb") as handle:
            handle.write(payload)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    finally:
        if temporary.exists():
            temporary.unlink()


def _default_runner(argv: Sequence[str]) -> int:
    return subprocess.run(list(argv), check=False).returncode


def _require_nonempty(path: Path, label: str) -> Path:
    path = Path(path).resolve()
    if not path.is_file() or path.stat().st_size == 0:
        raise AcquisitionError(f"{label} is missing or empty: {path}")
    return path


def _file_record(path: Path) -> dict[str, Any]:
    path = _require_nonempty(path, "fingerprinted file")
    return {
        "path": str(path),
        "size": path.stat().st_size,
        "sha256": manifest_tools.sha256_file(path),
    }


def _safe_cache_root(cache_root: Path) -> Path:
    """Keep CLI mutations below the repository's already ignored data tree."""

    cache_root = Path(cache_root).expanduser().resolve()
    data_root = (HERE / "data").resolve()
    if cache_root == data_root or not cache_root.is_relative_to(data_root):
        raise AcquisitionError(
            f"cache root must be a child of ignored {data_root}: {cache_root}"
        )
    return cache_root


def _safe_lock_path(lock_path: Path | None, cache_root: Path) -> Path:
    """Resolve a lock below its cache root, rejecting external mutations."""

    cache_root = Path(cache_root).expanduser().resolve()
    resolved = Path(lock_path or cache_root / LOCK_NAME).expanduser().resolve()
    if resolved == cache_root or not resolved.is_relative_to(cache_root):
        raise AcquisitionError(
            f"lock path must be a child of cache root {cache_root}: {resolved}"
        )
    return resolved


def _safe_zip_name(name: str) -> PurePosixPath:
    if not name or "\0" in name or "\\" in name:
        raise AcquisitionError(f"unsafe archive member name: {name!r}")
    path = PurePosixPath(name)
    if (
        path.is_absolute()
        or ":" in path.parts[0]
        or any(part in {"", ".", ".."} for part in path.parts)
    ):
        raise AcquisitionError(f"unsafe archive member name: {name!r}")
    return path


def validate_zip_archive(archive: Path) -> dict[str, zipfile.ZipInfo]:
    """Reject traversal, links, encryption, and duplicate member names."""

    archive = _require_nonempty(archive, "NCBI archive")
    try:
        with zipfile.ZipFile(archive) as handle:
            members: dict[str, zipfile.ZipInfo] = {}
            for info in handle.infolist():
                _safe_zip_name(info.filename)
                if info.filename in members:
                    raise AcquisitionError(
                        f"duplicate archive member: {info.filename}"
                    )
                mode = info.external_attr >> 16
                if stat.S_ISLNK(mode):
                    raise AcquisitionError(
                        f"archive links are not allowed: {info.filename}"
                    )
                if info.flag_bits & 0x1:
                    raise AcquisitionError(
                        f"encrypted archive member is not allowed: {info.filename}"
                    )
                members[info.filename] = info
            return members
    except (OSError, zipfile.BadZipFile) as exc:
        raise AcquisitionError(f"invalid ZIP archive {archive}: {exc}") from exc


def _publish_content(
    source: Path,
    *,
    filename: str,
    cache_root: Path,
    expected_sha256: str | None,
) -> dict[str, Any]:
    source = _require_nonempty(source, "acquired source file")
    if Path(filename).name != filename or filename in {"", ".", ".."}:
        raise AcquisitionError(f"unsafe content filename: {filename!r}")
    digest = manifest_tools.sha256_file(source)
    if expected_sha256 is not None and digest != expected_sha256:
        raise AcquisitionError(
            f"declared SHA-256 mismatch for {filename}: "
            f"expected {expected_sha256}, observed {digest}"
        )
    relative = Path("sha256") / digest[:2] / digest / filename
    destination = cache_root / relative
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists():
        observed = manifest_tools.sha256_file(destination)
        if observed != digest:
            raise AcquisitionError(
                f"content-addressed destination changed: {destination}"
            )
    else:
        descriptor, temporary_name = tempfile.mkstemp(
            prefix=f".{filename}.",
            suffix=".tmp",
            dir=destination.parent,
        )
        temporary = Path(temporary_name)
        try:
            with os.fdopen(descriptor, "wb") as output, source.open("rb") as input_:
                shutil.copyfileobj(input_, output, length=1024 * 1024)
                output.flush()
                os.fsync(output.fileno())
            os.replace(temporary, destination)
        finally:
            if temporary.exists():
                temporary.unlink()
    return {
        "path": str(relative),
        "bytes": destination.stat().st_size,
        "sha256": digest,
    }


def _run_command(argv: Sequence[str], runner: CommandRunner) -> None:
    if not argv or any(not isinstance(item, str) or "\0" in item for item in argv):
        raise AcquisitionError("refusing a malformed acquisition argument vector")
    try:
        return_code = runner(tuple(argv))
    except OSError as exc:
        raise AcquisitionError(
            f"cannot execute acquisition command {argv[0]!r}: {exc}"
        ) from exc
    if return_code != 0:
        raise AcquisitionError(
            f"acquisition command failed with exit {return_code}: "
            + " ".join(json.dumps(item) for item in argv)
        )


def _download_https_file(
    source: Mapping[str, Any],
    file_record: Mapping[str, Any],
    request_dir: Path,
    runner: CommandRunner,
) -> Path:
    locator = file_record["locator"]
    destination = request_dir / locator["filename"]
    if destination.is_file() and destination.stat().st_size:
        return destination
    partial = destination.with_name(destination.name + ".part")
    argv = [
        "curl",
        "--fail",
        "--location",
        "--retry",
        "3",
        "--continue-at",
        "-",
        "--output",
        str(partial),
        locator["url"],
    ]
    _run_command(argv, runner)
    _require_nonempty(partial, f"{source['id']}.{file_record['role']} partial")
    os.replace(partial, destination)
    return destination


def _download_ncbi_archive(
    source: Mapping[str, Any],
    request_dir: Path,
    runner: CommandRunner,
) -> Path:
    acquisition = source["acquisition"]
    destination = request_dir / acquisition["filename"]
    if destination.is_file() and destination.stat().st_size:
        return destination
    partial = destination.with_name(destination.name + ".part")
    if partial.exists():
        partial.unlink()
    argv = [
        "datasets",
        "download",
        "genome",
        "accession",
        source["identity"]["accession"],
        "--include",
        ",".join(acquisition["include"]),
        "--filename",
        str(partial),
    ]
    _run_command(argv, runner)
    _require_nonempty(partial, f"{source['id']} archive partial")
    validate_zip_archive(partial)
    os.replace(partial, destination)
    return destination


def _extract_member(
    archive: Path,
    member: str,
    destination: Path,
) -> Path:
    members = validate_zip_archive(archive)
    if member not in members or members[member].is_dir():
        raise AcquisitionError(f"declared archive member is missing: {member}")
    destination.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{destination.name}.",
        suffix=".tmp",
        dir=destination.parent,
    )
    temporary = Path(temporary_name)
    try:
        with (
            os.fdopen(descriptor, "wb") as output,
            zipfile.ZipFile(archive) as handle,
            handle.open(member) as input_,
        ):
            shutil.copyfileobj(input_, output, length=1024 * 1024)
            output.flush()
            os.fsync(output.fileno())
        _require_nonempty(temporary, f"extracted archive member {member}")
        os.replace(temporary, destination)
    except (OSError, zipfile.BadZipFile) as exc:
        raise AcquisitionError(
            f"cannot extract archive member {member}: {exc}"
        ) from exc
    finally:
        if temporary.exists():
            temporary.unlink()
    return destination


def _acquire_source(
    source: Mapping[str, Any],
    cache_root: Path,
    runner: CommandRunner,
) -> dict[str, Any]:
    request_sha = manifest_tools.source_request_sha256(source)
    request_dir = cache_root / "requests" / request_sha
    request_dir.mkdir(parents=True, exist_ok=True)
    locked_files = {}
    if source["transport"] == "https_files":
        for file_record in source["files"]:
            acquired = _download_https_file(
                source, file_record, request_dir, runner,
            )
            locked_files[file_record["role"]] = _publish_content(
                acquired,
                filename=file_record["locator"]["filename"],
                cache_root=cache_root,
                expected_sha256=file_record["expected_sha256"],
            )
    elif source["transport"] == "ncbi_datasets":
        archive = _download_ncbi_archive(source, request_dir, runner)
        extracted_dir = request_dir / "extracted"
        for file_record in source["files"]:
            member = file_record["locator"]["archive_member"]
            filename = PurePosixPath(member).name
            extracted = _extract_member(
                archive,
                member,
                extracted_dir / f"{file_record['role']}-{filename}",
            )
            locked_files[file_record["role"]] = _publish_content(
                extracted,
                filename=filename,
                cache_root=cache_root,
                expected_sha256=file_record["expected_sha256"],
            )
    else:
        raise AcquisitionError(
            f"source {source['id']} is not remotely acquirable"
        )
    return {
        "request_sha256": request_sha,
        "files": locked_files,
    }


def acquire(
    manifest: Mapping[str, Any],
    *,
    cache_root: Path,
    lock_path: Path | None = None,
    accept_identity_pinned_bytes: bool = False,
    runner: CommandRunner = _default_runner,
) -> tuple[dict[str, Any], dict[str, Any]]:
    """Acquire every remote source and atomically publish a verified lock."""

    manifest = manifest_tools.validate_manifest(manifest)
    cache_root = Path(cache_root).resolve()
    lock_path = _safe_lock_path(lock_path, cache_root)
    if lock_path.exists():
        lock = manifest_tools.load_json(lock_path)
        verification = manifest_tools.verify_acquisition_lock(
            manifest, lock, cache_root,
        )
        return lock, {**verification, "resumed_from_verified_lock": True}
    pending = sorted(
        f"{source['id']}:{item['role']}"
        for source in manifest["sources"]
        if source["transport"] != "existing_registry"
        for item in source["files"]
        if item["pin_state"] == REMOTE_PENDING_PIN
    )
    if pending and not accept_identity_pinned_bytes:
        raise AcquisitionError(
            "upstream checksums are unavailable for identity-pinned bytes; "
            "review the dry-run plan and pass "
            "--accept-identity-pinned-bytes to lock the first observed bytes. "
            f"Pending roles: {pending}"
        )
    cache_root.mkdir(parents=True, exist_ok=True)
    sources = {}
    for source in manifest["sources"]:
        if source["transport"] == "existing_registry":
            continue
        sources[source["id"]] = _acquire_source(source, cache_root, runner)
    lock = {
        "schema_version": manifest_tools.LOCK_SCHEMA_VERSION,
        "manifest_sha256": manifest_tools.canonical_sha256(manifest),
        "sources": sources,
    }
    verification = manifest_tools.verify_acquisition_lock(
        manifest, lock, cache_root,
    )
    _atomic_write_json(lock_path, lock)
    return lock, {**verification, "resumed_from_verified_lock": False}


def acquisition_dry_run(
    manifest: Mapping[str, Any],
    cache_root: Path,
    *,
    lock_path: Path | None = None,
) -> dict[str, Any]:
    manifest = manifest_tools.validate_manifest(manifest)
    cache_root = Path(cache_root).resolve()
    resolved_lock = _safe_lock_path(lock_path, cache_root)
    plan = manifest_tools.build_acquisition_plan(manifest, cache_root)
    by_id = {source["id"]: source for source in manifest["sources"]}
    for step in plan["steps"]:
        source = by_id[step["source_id"]]
        request_dir = Path(step["staging_directory"])
        if source["transport"] == "ncbi_datasets":
            acquisition = source["acquisition"]
            step["commands"] = [[
                "datasets",
                "download",
                "genome",
                "accession",
                source["identity"]["accession"],
                "--include",
                ",".join(acquisition["include"]),
                "--filename",
                str(request_dir / f"{acquisition['filename']}.part"),
            ]]
        else:
            step["commands"] = [
                [
                    "curl",
                    "--fail",
                    "--location",
                    "--retry",
                    "3",
                    "--continue-at",
                    "-",
                    "--output",
                    str(request_dir / f"{item['locator']['filename']}.part"),
                    item["locator"]["url"],
                ]
                for item in source["files"]
            ]
    plan.pop("plan_sha256")
    plan["plan_sha256"] = manifest_tools.canonical_sha256(plan)
    pending = sorted(
        f"{source['id']}:{item['role']}"
        for source in manifest["sources"]
        if source["transport"] != "existing_registry"
        for item in source["files"]
        if item["pin_state"] == REMOTE_PENDING_PIN
    )
    return {
        **plan,
        "execution": {
            "shell": False,
            "atomic_lock": str(resolved_lock),
            "resume_policy": (
                "verified lock and completed request files are reused; "
                "HTTPS .part files resume; incomplete NCBI archives restart"
            ),
            "requires_identity_pin_acknowledgement": bool(pending),
            "identity_pinned_roles": pending,
        },
    }


def _locked_source_paths(
    manifest: Mapping[str, Any],
    lock: Mapping[str, Any],
    cache_root: Path,
) -> dict[str, dict[str, Path]]:
    manifest_tools.verify_acquisition_lock(manifest, lock, cache_root)
    result = {}
    for source in manifest["sources"]:
        if source["transport"] == "existing_registry":
            continue
        result[source["id"]] = {
            role: (cache_root / record["path"]).resolve()
            for role, record in lock["sources"][source["id"]]["files"].items()
        }
    return result


def _load_benchmark_registry(path: Path) -> dict[str, Any]:
    try:
        document = json.loads(Path(path).read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise AcquisitionError(f"cannot load benchmark registry {path}: {exc}") from exc
    entries = document.get("benchmarks") if isinstance(document, dict) else None
    if not isinstance(entries, list):
        raise AcquisitionError("benchmark registry has no benchmarks list")
    ids = [
        row.get("id") for row in entries
        if isinstance(row, dict) and row.get("id")
    ]
    if len(ids) != len(set(ids)):
        raise AcquisitionError("benchmark registry IDs are duplicated")
    return document


def _existing_source_paths(
    manifest: Mapping[str, Any],
    benchmark_registry: Mapping[str, Any],
) -> dict[str, dict[str, Path]]:
    by_id = {
        row["id"]: row
        for row in benchmark_registry["benchmarks"]
        if isinstance(row, dict) and row.get("id")
    }
    result = {}
    for source in manifest["sources"]:
        if source["transport"] != "existing_registry":
            continue
        benchmark_id = source["identity"]["benchmark_id"]
        if benchmark_id not in by_id:
            raise AcquisitionError(
                f"existing benchmark is missing from registry: {benchmark_id}"
            )
        benchmark = by_id[benchmark_id]
        roles = {}
        for file_record in source["files"]:
            key = file_record["locator"]["registry_key"]
            value = benchmark.get(key)
            if not value:
                raise AcquisitionError(
                    f"{benchmark_id}.{key} is not a usable registry path"
                )
            roles[file_record["role"]] = _require_nonempty(
                Path(value), f"{benchmark_id}.{key}",
            )
        result[source["id"]] = roles
    return result


def _runtime_name(source: Path) -> tuple[str, bool]:
    if source.name.lower().endswith(".gz"):
        name = source.name[:-3]
        if not name:
            raise AcquisitionError(f"compressed input has no runtime name: {source}")
        return name, True
    return source.name, False


def _atomic_symlink(source: Path, destination: Path) -> None:
    source = Path(source).resolve()
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.is_symlink():
        if destination.resolve() != source:
            raise AcquisitionError(
                f"refusing to replace materialized symlink: {destination}"
            )
        return
    if destination.exists():
        if (
            destination.is_file()
            and destination.stat().st_size == source.stat().st_size
            and manifest_tools.sha256_file(destination)
            == manifest_tools.sha256_file(source)
        ):
            return
        raise AcquisitionError(
            f"refusing to replace materialized path: {destination}"
        )
    temporary = destination.with_name(
        f".{destination.name}.{os.getpid()}.tmp"
    )
    if temporary.exists() or temporary.is_symlink():
        temporary.unlink()
    temporary.symlink_to(source)
    os.replace(temporary, destination)


def _decompress_gzip(source: Path, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{destination.name}.",
        suffix=".tmp",
        dir=destination.parent,
    )
    temporary = Path(temporary_name)
    try:
        with (
            os.fdopen(descriptor, "wb") as output,
            gzip.open(source, "rb") as input_,
        ):
            shutil.copyfileobj(input_, output, length=1024 * 1024)
            output.flush()
            os.fsync(output.fileno())
        _require_nonempty(temporary, f"decompressed runtime input {source}")
        if destination.exists():
            if (
                destination.stat().st_size == temporary.stat().st_size
                and manifest_tools.sha256_file(destination)
                == manifest_tools.sha256_file(temporary)
            ):
                return
            raise AcquisitionError(
                f"decompressed runtime input changed: {destination}"
            )
        os.replace(temporary, destination)
    except (OSError, gzip.BadGzipFile, EOFError) as exc:
        raise AcquisitionError(f"cannot decompress {source}: {exc}") from exc
    finally:
        if temporary.exists():
            temporary.unlink()


def _fasta_index_payload(path: Path) -> bytes:
    records = []
    seen = set()
    with Path(path).open("rb") as handle:
        current: dict[str, Any] | None = None
        previous_bases = 0

        def finish() -> None:
            if current is None:
                return
            if current["length"] == 0:
                raise AcquisitionError(
                    f"FASTA record {current['name']!r} has no sequence: {path}"
                )
            records.append(
                (
                    current["name"],
                    current["length"],
                    current["offset"],
                    current["line_bases"],
                    current["line_width"],
                )
            )

        while True:
            offset = handle.tell()
            raw = handle.readline()
            if not raw:
                break
            if raw.startswith(b">"):
                finish()
                name = raw[1:].strip().split(None, 1)[0]
                try:
                    decoded = name.decode("ascii")
                except UnicodeDecodeError as exc:
                    raise AcquisitionError(
                        f"non-ASCII FASTA identifier in {path}"
                    ) from exc
                if not decoded or decoded in seen:
                    raise AcquisitionError(
                        f"empty or duplicate FASTA identifier {decoded!r}: {path}"
                    )
                seen.add(decoded)
                current = {
                    "name": decoded,
                    "length": 0,
                    "offset": handle.tell(),
                    "line_bases": 0,
                    "line_width": 0,
                }
                previous_bases = 0
                continue
            if current is None:
                if raw.strip():
                    raise AcquisitionError(
                        f"FASTA sequence precedes the first header: {path}"
                    )
                continue
            sequence = raw.rstrip(b"\r\n")
            if not sequence:
                raise AcquisitionError(f"blank FASTA sequence line: {path}")
            if any(byte in b" \t" for byte in sequence):
                raise AcquisitionError(f"whitespace inside FASTA sequence: {path}")
            bases = len(sequence)
            width = len(raw)
            if current["line_bases"] == 0:
                current["line_bases"] = bases
                current["line_width"] = width
            elif previous_bases < current["line_bases"]:
                raise AcquisitionError(
                    f"short FASTA line is not final in its record: {path}"
                )
            elif (
                bases > current["line_bases"]
                or (
                    bases == current["line_bases"]
                    and width != current["line_width"]
                )
            ):
                raise AcquisitionError(
                    f"inconsistent FASTA line width at byte {offset}: {path}"
                )
            current["length"] += bases
            previous_bases = bases
        finish()
    if not records:
        raise AcquisitionError(f"FASTA contains no records: {path}")
    return "".join(
        f"{name}\t{length}\t{offset}\t{line_bases}\t{line_width}\n"
        for name, length, offset, line_bases, line_width in records
    ).encode("ascii")


def ensure_fasta_index(path: Path) -> Path:
    payload = _fasta_index_payload(path)
    destination = Path(str(path) + ".fai")
    if destination.exists():
        if destination.read_bytes() != payload:
            raise AcquisitionError(f"stale or invalid FASTA index: {destination}")
    else:
        _atomic_write_bytes(destination, payload)
    return destination


def _materialize_asset(
    source: Path,
    destination_dir: Path,
    *,
    fasta: bool,
) -> tuple[Path, dict[str, Any]]:
    source_record = _file_record(source)
    runtime_name, compressed = _runtime_name(source)
    destination = destination_dir / runtime_name
    if compressed:
        _decompress_gzip(source, destination)
        transform = "gzip-decompress"
    else:
        _atomic_symlink(source, destination)
        transform = "verified-symlink"
    runtime_record = _file_record(destination)
    record = {
        "source": source_record,
        "runtime": runtime_record,
        "transform": transform,
    }
    if fasta:
        index = ensure_fasta_index(destination)
        record["fasta_index"] = _file_record(index)
    return destination, record


def materialize_source_assets(
    manifest: Mapping[str, Any],
    source_paths: Mapping[str, Mapping[str, Path]],
    runtime_root: Path,
) -> tuple[dict[str, dict[str, Path]], dict[str, Any]]:
    runtime_paths = {}
    provenance = {}
    for source in manifest["sources"]:
        source_id = source["id"]
        if source_id not in source_paths:
            raise AcquisitionError(f"source bytes are unresolved: {source_id}")
        runtime_paths[source_id] = {}
        provenance[source_id] = {}
        for file_record in source["files"]:
            role = file_record["role"]
            path, record = _materialize_asset(
                source_paths[source_id][role],
                runtime_root / "sources" / source_id / role,
                fasta=role in FASTA_ROLES or role == "genome",
            )
            runtime_paths[source_id][role] = path
            provenance[source_id][role] = record
    return runtime_paths, provenance


def _resolve_scenario_inputs(
    scenario: Mapping[str, Any],
    runtime_paths: Mapping[str, Mapping[str, Path]],
) -> dict[str, Path]:
    resolved = {}
    for role, binding in scenario["inputs"].items():
        source_id, source_role = binding.split(":", 1)
        try:
            resolved[role] = runtime_paths[source_id][source_role]
        except KeyError as exc:
            raise AcquisitionError(
                f"{scenario['id']}.{role} is unresolved: {binding}"
            ) from exc
    return resolved


def _ortholog_mappings(
    path: Path | None,
    manifest: Mapping[str, Any],
    scenario_inputs: Mapping[str, Mapping[str, Path]],
) -> tuple[dict[str, dict[str, Any]], list[dict[str, str]]]:
    biological = {
        scenario["id"] for scenario in manifest["scenarios"]
        if scenario["kind"] == "biological"
    }
    if path is None:
        return {}, [
            {
                "code": "missing_ortholog_map",
                "scenario_id": scenario_id,
                "message": (
                    "independent target truth requires a curated source-to-target "
                    "gene/transcript ortholog map"
                ),
            }
            for scenario_id in sorted(biological)
        ]
    path = Path(path).resolve()
    try:
        registry_document = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise AcquisitionError(f"cannot load ortholog registry {path}: {exc}") from exc
    if (
        not isinstance(registry_document, dict)
        or set(registry_document) != {
            "schema_version", "manifest_sha256", "mappings",
        }
        or registry_document["schema_version"] != ORTHOLOG_SCHEMA_VERSION
        or registry_document["manifest_sha256"]
        != manifest_tools.canonical_sha256(manifest)
        or not isinstance(registry_document["mappings"], dict)
    ):
        raise AcquisitionError("ortholog registry schema or manifest SHA is invalid")
    unknown = set(registry_document["mappings"]) - biological
    if unknown:
        raise AcquisitionError(
            f"ortholog registry has unexpected scenarios: {sorted(unknown)}"
        )
    mappings = {}
    blockers = []
    for scenario_id in sorted(biological):
        raw = registry_document["mappings"].get(scenario_id)
        if not isinstance(raw, dict) or set(raw) != {
            "path", "sha256", "id_policy",
        }:
            blockers.append({
                "code": "missing_ortholog_map",
                "scenario_id": scenario_id,
                "message": "ortholog registry has no complete mapping record",
            })
            continue
        if raw["id_policy"] != "ortholog-map":
            raise AcquisitionError(
                f"{scenario_id}: biological truth requires ortholog-map policy"
            )
        mapping_path = Path(str(raw["path"])).expanduser()
        if not mapping_path.is_absolute():
            mapping_path = path.parent / mapping_path
        mapping_path = _require_nonempty(
            mapping_path, f"{scenario_id} ortholog map",
        )
        digest = manifest_tools.sha256_file(mapping_path)
        if digest != raw["sha256"]:
            raise AcquisitionError(
                f"{scenario_id}: ortholog map SHA-256 changed"
            )
        from benchmarks.compare import ortholog_scope, target_truth

        try:
            mapping_document = (
                ortholog_scope.validate_mapping_from_recorded_annotations(
                    json.loads(mapping_path.read_text(encoding="utf-8"))
                )
            )
        except (
            OSError,
            json.JSONDecodeError,
            ortholog_scope.ScopeBuildError,
        ) as exc:
            raise AcquisitionError(
                f"{scenario_id}: unsupported ortholog scope: {exc}"
            ) from exc
        metadata = mapping_document.get("metadata")
        provenance = (
            metadata.get("provenance")
            if isinstance(metadata, Mapping)
            else None
        )
        provenance_inputs = (
            provenance.get("inputs")
            if isinstance(provenance, Mapping)
            else None
        )
        expected_inputs = {
            "source_annotation": "reference_annotation",
            "source_genome": "reference_genome",
            "target_annotation": "target_truth",
            "target_genome": "target_genome",
        }
        if (
            not isinstance(metadata, Mapping)
            or metadata.get("scope") != "full"
            or not isinstance(provenance_inputs, Mapping)
        ):
            raise AcquisitionError(
                f"{scenario_id}: ortholog scope lacks full-scope provenance"
            )
        for provenance_role, scenario_role in expected_inputs.items():
            record = provenance_inputs.get(provenance_role)
            expected_path = scenario_inputs[scenario_id][scenario_role]
            if (
                not isinstance(record, Mapping)
                or record.get("sha256")
                != manifest_tools.sha256_file(expected_path)
                or record.get("size") != expected_path.stat().st_size
            ):
                raise AcquisitionError(
                    f"{scenario_id}: ortholog scope input fingerprint "
                    f"does not match {scenario_role}"
                )
        commands = provenance.get("commands")
        tools = provenance.get("tools")
        if (
            not isinstance(commands, list)
            or not commands
            or any(
                not isinstance(command, Mapping)
                or command.get("shell") is not False
                or not isinstance(command.get("argv"), list)
                or not command["argv"]
                for command in commands
            )
            or not isinstance(tools, Mapping)
            or set(tools) != {"gffread", "mmseqs"}
        ):
            raise AcquisitionError(
                f"{scenario_id}: ortholog scope tool/command provenance "
                "is incomplete"
            )
        counts = metadata.get("counts")
        if not isinstance(counts, Mapping):
            raise AcquisitionError(
                f"{scenario_id}: ortholog scope counts are missing"
            )
        for key in (
            "gene_groups_retained",
            "transcript_groups_retained",
        ):
            value = counts.get(key)
            feature_type = key.removesuffix("_groups_retained")
            observed = sum(
                entry.get("feature_type") == feature_type
                and entry.get("status") == "retained"
                for entry in mapping_document["mappings"]
            )
            if (
                not isinstance(value, int)
                or isinstance(value, bool)
                or value != observed
                or value < MINIMUM_FULL_TRUTH_GROUPS
            ):
                raise AcquisitionError(
                    f"{scenario_id}: {key} must equal the retained entries "
                    f"and be at least {MINIMUM_FULL_TRUTH_GROUPS}; "
                    f"declared={value!r}, observed={observed}"
                )
        entries, _mapping_document = target_truth.load_ortholog_map(
            mapping_path
        )
        mappings[scenario_id] = {
            "path": mapping_path,
            "sha256": digest,
            "id_policy": "ortholog-map",
            "entries": len(entries),
            "method": mapping_document["method"],
            "scope": metadata["scope"],
            "counts": {
                key: counts[key]
                for key in (
                    "gene_groups_retained",
                    "transcript_groups_retained",
                )
            },
        }
    return mappings, blockers


def _default_synthetic_builder(
    scenario: Mapping[str, Any],
    source_fasta: Path,
    source_gff: Path,
    output_dir: Path,
) -> Mapping[str, Path]:
    from benchmarks.compare import synthetic_truth

    scenario_id = scenario["id"]
    kwargs = manifest_tools.synthetic_builder_kwargs(scenario)
    if scenario_id == "v2_synth_chr22_fragmented":
        return synthetic_truth.build_fragmented_case(
            source_fasta, source_gff, output_dir, **kwargs,
        )
    if scenario_id == "v2_synth_chr22_sv":
        return synthetic_truth.build_sv_case(
            source_fasta, source_gff, output_dir, **kwargs,
        )
    raise AcquisitionError(f"unknown synthetic scenario: {scenario_id}")


def _verified_synthetic_outputs(
    output_dir: Path,
    scenario: Mapping[str, Any],
    source_fasta: Path,
    source_gff: Path,
) -> tuple[dict[str, Path], dict[str, Any]] | None:
    manifest_path = output_dir / "transform.manifest.json"
    if not manifest_path.is_file():
        return None
    try:
        document = json.loads(manifest_path.read_text(encoding="utf-8"))
        expected_transform = manifest_tools.expected_synthetic_transform(scenario)
        if (
            not isinstance(document, dict)
            or document.get("schema_version") != 1
            or document.get("kind") != "lifton-synthetic-target-truth"
            or document.get("transform") != expected_transform
            or document.get("source", {}).get("seqid")
            != scenario["design"]["chromosome"]
            or document["source"]["fasta"]["sha256"]
            != manifest_tools.sha256_file(source_fasta)
            or document["source"]["gff"]["sha256"]
            != manifest_tools.sha256_file(source_gff)
        ):
            return None
        declared_outputs = document.get("outputs")
        expected_names = {
            "target_fasta": "target.fa",
            "truth_gff": "target.truth.gff3",
            "ortholog_map": "ortholog_map.json",
        }
        if (
            not isinstance(declared_outputs, dict)
            or set(declared_outputs) != set(expected_names)
        ):
            return None
        outputs = {
            role: output_dir / expected_name
            for role, expected_name in expected_names.items()
        }
        for role, expected_name in expected_names.items():
            record = declared_outputs[role]
            if (
                not isinstance(record, dict)
                or record.get("name") != expected_name
                or not isinstance(record.get("sha256"), str)
                or len(record["sha256"]) != 64
            ):
                return None
        for role in expected_names:
            if (
                not outputs[role].is_file()
                or outputs[role].stat().st_size <= 0
                or manifest_tools.sha256_file(outputs[role])
                != declared_outputs[role]["sha256"]
            ):
                return None
        provenance = {
            "design_sha256": manifest_tools.canonical_sha256(
                scenario["design"]
            ),
            "source": {
                "fasta": _file_record(source_fasta),
                "gff": _file_record(source_gff),
            },
            "transform": expected_transform,
            "transform_manifest": _file_record(manifest_path),
            "outputs": {
                role: _file_record(path)
                for role, path in sorted(outputs.items())
            },
        }
        return outputs, provenance
    except (AttributeError, KeyError, OSError, ValueError, TypeError):
        return None


def generate_synthetic_scenarios(
    manifest: Mapping[str, Any],
    runtime_paths: Mapping[str, Mapping[str, Path]],
    runtime_root: Path,
    *,
    builder: SyntheticBuilder = _default_synthetic_builder,
) -> tuple[dict[str, dict[str, Path]], dict[str, Any]]:
    outputs = {}
    provenance = {}
    for scenario in manifest["scenarios"]:
        if scenario["kind"] != "synthetic":
            continue
        inputs = _resolve_scenario_inputs(scenario, runtime_paths)
        destination = runtime_root / "scenarios" / scenario["id"]
        verified_bundle = _verified_synthetic_outputs(
            destination,
            scenario,
            inputs["reference_genome"],
            inputs["reference_annotation"],
        )
        if verified_bundle is None:
            if destination.exists() and any(destination.iterdir()):
                raise AcquisitionError(
                    f"stale synthetic output requires manual quarantine: {destination}"
                )
            destination.parent.mkdir(parents=True, exist_ok=True)
            temporary = Path(tempfile.mkdtemp(
                prefix=f".{scenario['id']}.",
                dir=destination.parent,
            ))
            try:
                builder(
                    scenario,
                    inputs["reference_genome"],
                    inputs["reference_annotation"],
                    temporary,
                )
                if destination.exists():
                    destination.rmdir()
                os.replace(temporary, destination)
            finally:
                if temporary.exists():
                    shutil.rmtree(temporary)
            verified_bundle = _verified_synthetic_outputs(
                destination,
                scenario,
                inputs["reference_genome"],
                inputs["reference_annotation"],
            )
        if verified_bundle is None:
            raise AcquisitionError(
                f"synthetic scenario did not publish verified outputs: "
                f"{scenario['id']}"
            )
        verified, verified_provenance = verified_bundle
        ensure_fasta_index(verified["target_fasta"])
        outputs[scenario["id"]] = {
            **inputs,
            **verified,
        }
        provenance[scenario["id"]] = verified_provenance
    return outputs, provenance


def _fasta_seqids(path: Path) -> set[str]:
    index = ensure_fasta_index(path)
    return {
        line.split("\t", 1)[0]
        for line in index.read_text(encoding="ascii").splitlines()
        if line
    }


def _annotation_seqids(path: Path) -> set[str]:
    result = set()
    with Path(path).open("rt", encoding="utf-8", errors="strict") as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.split("\t", 1)
            if len(fields) == 2:
                result.add(fields[0])
    if not result:
        raise AcquisitionError(f"annotation contains no feature records: {path}")
    return result


def _seqid_blockers(
    scenario: Mapping[str, Any],
    inputs: Mapping[str, Path],
) -> list[dict[str, str]]:
    blockers = []
    pairs = (
        ("reference_annotation", "reference_genome"),
        ("target_truth", "target_genome"),
    )
    for annotation_role, genome_role in pairs:
        if annotation_role not in inputs or genome_role not in inputs:
            continue
        overlap = (
            _annotation_seqids(inputs[annotation_role])
            & _fasta_seqids(inputs[genome_role])
        )
        if not overlap:
            blockers.append({
                "code": "sequence_id_mismatch",
                "scenario_id": scenario["id"],
                "message": (
                    f"{annotation_role} and {genome_role} have no exact "
                    "sequence identifier in common"
                ),
            })
    if (
        "chr22" in str(scenario["design"].get("reference_scope", ""))
        and _reference_subset_seqid(inputs) is None
    ):
        blockers.append({
            "code": "reference_subset_seqid_unresolved",
            "scenario_id": scenario["id"],
            "message": (
                "the chr22 subset scope cannot be resolved to one exact "
                "reference annotation/FASTA sequence identifier"
            ),
        })
    return blockers


def _reference_subset_seqid(
    inputs: Mapping[str, Path],
) -> str | None:
    annotation = inputs["reference_annotation"]
    common = (
        _annotation_seqids(annotation)
        & _fasta_seqids(inputs["reference_genome"])
    )
    for candidate in ("chr22", "22"):
        if candidate in common:
            return candidate
    with Path(annotation).open(
        "rt", encoding="utf-8", errors="strict",
    ) as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            columns = line.rstrip("\r\n").split("\t")
            if len(columns) != 9 or columns[0] not in common:
                continue
            attributes = {}
            for item in columns[8].split(";"):
                key, separator, value = item.strip().partition("=")
                if separator and key not in attributes:
                    attributes[key] = value
            if (
                attributes.get("chromosome") == "22"
                or attributes.get("Name") == "22"
            ):
                return columns[0]
    return None


def _annotation_database(scenario: Mapping[str, Any]) -> str:
    if scenario["id"] == "v2_dialect_ensembl116_gtf":
        return "ENSEMBL"
    if scenario["id"].startswith("v2_dialect_"):
        return "OTHER"
    return "RefSeq"


def _cross_species(scenario: Mapping[str, Any]) -> bool:
    return scenario["id"] not in {
        "v2_truth_human_grch38_chm13",
        "v2_dialect_ensembl116_gtf",
        "v2_truth_soybean_w82_lee",
    }


def _benchmark_row(
    scenario: Mapping[str, Any],
    inputs: Mapping[str, Path],
    mapping: Mapping[str, Any],
) -> dict[str, Any]:
    design = scenario["design"]
    scope = str(design.get("reference_scope", ""))
    ref_chrom = (
        _reference_subset_seqid(inputs)
        if "chr22" in scope
        else "AUTO_LARGEST_CODING"
    )
    if ref_chrom is None:
        raise AcquisitionError(
            f"{scenario['id']}: reference chr22 sequence ID is unresolved"
        )
    target_chrom = (
        "WHOLE"
        if scenario["id"] in {
            "v2_deep_zebrafish_xenopus",
            "v2_deep_tomato_rice",
        }
        else "AUTO_SYNTENIC"
    )
    return {
        "id": scenario["id"],
        "species": (
            f"{design.get('reference_label', scenario['id'])} -> "
            f"{design.get('target_label', scenario['id'])}"
        ),
        "cross_species": _cross_species(scenario),
        "dimension": "canonical_v2_truth",
        "ref_genome": str(inputs["reference_genome"].resolve()),
        "ref_gff": str(inputs["reference_annotation"].resolve()),
        "ref_proteins": (
            str(inputs["reference_proteins"].resolve())
            if inputs.get("reference_proteins")
            else None
        ),
        "tgt_genome": str(inputs["target_genome"].resolve()),
        "ref_chrom": ref_chrom,
        "tgt_chrom": target_chrom,
        "full_input_mode": "raw",
        "miniprot_target_space": "transcript",
        "annotation_database": _annotation_database(scenario),
        "input_format": design.get("input_format", "GFF3"),
        "target_truth": {
            "gff": str(inputs["target_truth"].resolve()),
            "ortholog_map": str(mapping["path"].resolve()),
            "id_policy": mapping["id_policy"],
        },
    }


def _synthetic_benchmark_row(
    scenario: Mapping[str, Any],
    outputs: Mapping[str, Path],
) -> dict[str, Any]:
    truth = {
        "gff": str(outputs["truth_gff"].resolve()),
        "ortholog_map": str(outputs["ortholog_map"].resolve()),
        "id_policy": "ortholog-map",
    }
    return {
        "id": scenario["id"],
        "species": f"Synthetic human chr22: {scenario['id']}",
        "cross_species": False,
        "dimension": "canonical_v2_synthetic_truth",
        "ref_genome": str(outputs["reference_genome"].resolve()),
        "ref_gff": str(outputs["reference_annotation"].resolve()),
        "ref_proteins": None,
        "tgt_genome": str(outputs["target_fasta"].resolve()),
        "ref_chrom": "chr22",
        "tgt_chrom": "WHOLE",
        "full_input_mode": "raw",
        "miniprot_target_space": "transcript",
        "annotation_database": "RefSeq",
        "input_format": "GFF3",
        "target_truth": dict(truth),
        "target_truth_by_panel": {
            "subset": dict(truth),
            "full": dict(truth),
        },
    }


def _materialize_subset_truth(
    row: Mapping[str, Any],
    mapping: Mapping[str, Any],
    *,
    work_root: Path,
    runtime_root: Path,
) -> tuple[dict[str, Any], dict[str, Any]]:
    """Create truth scoped to the exact prepared source/target subset."""

    from benchmarks.compare import ortholog_scope, target_truth

    scenario_id = str(row["id"])
    subset_manifest_path = (
        Path(work_root) / scenario_id / "subset" / "subset.manifest.json"
    ).resolve()
    try:
        subset_manifest = json.loads(
            subset_manifest_path.read_text(encoding="utf-8")
        )
    except (OSError, json.JSONDecodeError) as exc:
        raise AcquisitionError(
            f"{scenario_id}: cannot read verified subset manifest: {exc}"
        ) from exc
    if (
        subset_manifest.get("schema_version") != 2
        or subset_manifest.get("id") != scenario_id
        or not isinstance(subset_manifest.get("outputs"), Mapping)
        or not isinstance(subset_manifest.get("paths"), Mapping)
    ):
        raise AcquisitionError(
            f"{scenario_id}: subset manifest is not schema-v2 provenance"
        )
    for label, record in subset_manifest["outputs"].items():
        if not isinstance(record, Mapping):
            raise AcquisitionError(
                f"{scenario_id}: malformed subset output record {label!r}"
            )
        try:
            live = _file_record(Path(str(record["path"])))
        except (KeyError, OSError) as exc:
            raise AcquisitionError(
                f"{scenario_id}: subset output is unavailable: {label}"
            ) from exc
        if dict(record) != live:
            raise AcquisitionError(
                f"{scenario_id}: subset output changed after preparation: {label}"
            )
    paths = subset_manifest["paths"]
    source_subset = Path(str(paths["ref_gff"])).resolve()
    target_subset_fasta = Path(str(paths["tgt_fa"])).resolve()
    output_root = runtime_root / "truth_scopes" / scenario_id / "subset"
    output_root.mkdir(parents=True, exist_ok=True)
    truth_path = output_root / "target.truth.gff3"
    map_path = output_root / "ortholog_map.json"
    filter_record = target_truth.filter_truth_to_target_fasta(
        row["target_truth"]["gff"],
        target_subset_fasta,
        truth_path,
        force=True,
    )
    try:
        filtered_mapping = ortholog_scope.write_filtered_scope(
            mapping["path"],
            source_subset,
            truth_path,
            map_path,
            minimum_gene_groups=10,
            minimum_transcript_groups=10,
        )
    except (OSError, ValueError, ortholog_scope.ScopeBuildError) as exc:
        raise AcquisitionError(
            f"{scenario_id}: cannot materialize subset ortholog truth: {exc}"
        ) from exc
    truth = {
        "gff": str(truth_path.resolve()),
        "ortholog_map": str(map_path.resolve()),
        "id_policy": "ortholog-map",
    }
    provenance = {
        "schema_version": 2,
        "method": "canonical-v2-panel-truth-v2",
        "scenario_id": scenario_id,
        "inputs": {
            "subset_manifest": _file_record(subset_manifest_path),
            "source_subset_annotation": _file_record(source_subset),
            "target_subset_fasta": _file_record(target_subset_fasta),
            "full_truth_gff": _file_record(Path(row["target_truth"]["gff"])),
            "full_ortholog_map": _file_record(Path(mapping["path"])),
        },
        "filter": filter_record,
        "mapping": {
            "sha256": manifest_tools.canonical_sha256(filtered_mapping),
            "counts": filtered_mapping["metadata"]["counts"],
        },
        "outputs": {
            "truth_gff": _file_record(truth_path),
            "ortholog_map": _file_record(map_path),
        },
    }
    provenance_path = output_root / "panel_truth.manifest.json"
    _atomic_write_json(provenance_path, provenance)
    provenance["manifest"] = _file_record(provenance_path)
    return truth, provenance


def _merge_rows(
    base: Sequence[Mapping[str, Any]],
    additions: Sequence[Mapping[str, Any]],
) -> list[dict[str, Any]]:
    addition_ids = {row["id"] for row in additions}
    if len(addition_ids) != len(additions):
        raise AcquisitionError("generated registry IDs are duplicated")
    overlap = addition_ids & {
        row.get("id") for row in base if isinstance(row, Mapping)
    }
    if overlap:
        raise AcquisitionError(
            f"generated registry IDs already exist: {sorted(overlap)}"
        )
    return [dict(row) for row in base] + [dict(row) for row in additions]


def _e2e_dataset_row(
    row: Mapping[str, Any],
    data_root: Path,
) -> dict[str, Any]:
    scenario_id = row["id"]
    destination = data_root / scenario_id
    destination.mkdir(parents=True, exist_ok=True)
    values = {
        "reference.fa": Path(row["ref_genome"]),
        "reference.gff": Path(row["ref_gff"]),
        "target.fa": Path(row["tgt_genome"]),
        "target.truth.gff3": Path(row["target_truth"]["gff"]),
        "ortholog_map.json": Path(row["target_truth"]["ortholog_map"]),
    }
    for filename, source in values.items():
        _atomic_symlink(source, destination / filename)
    return {
        "id": scenario_id,
        "species": row["species"],
        "reference_fa": f"file://{destination / 'reference.fa'}",
        "target_fa": f"file://{destination / 'target.fa'}",
        "reference_gff": f"file://{destination / 'reference.gff'}",
        "target_gff": None,
        "approx_size_gb": 0.0,
        "cross_species": bool(row["cross_species"]),
        "annotation_database": row["annotation_database"],
        "truth_gff": f"file://{destination / 'target.truth.gff3'}",
        "ortholog_map": f"file://{destination / 'ortholog_map.json'}",
        "truth_id_policy": row["target_truth"]["id_policy"],
    }


def _subset_preflight(
    benchmark_ids: Sequence[str],
    *,
    benchmark_registry: Path,
    work_root: Path,
    threads: int,
) -> tuple[list[dict[str, str]], list[list[str]]]:
    from benchmarks.compare import subset_builder, tool_runners

    try:
        registry_document = json.loads(
            Path(benchmark_registry).read_text(encoding="utf-8")
        )
        tools = registry_document["tools"]
        benchmark_rows = {
            row["id"]: row for row in registry_document["benchmarks"]
            if isinstance(row, Mapping) and isinstance(row.get("id"), str)
        }
    except (KeyError, OSError, TypeError, json.JSONDecodeError) as exc:
        raise AcquisitionError(
            f"cannot verify subset preparation registry: {exc}"
        ) from exc
    blockers = []
    actions = []
    for benchmark_id in sorted(benchmark_ids):
        work = Path(work_root) / benchmark_id
        bench = benchmark_rows.get(benchmark_id)
        verified = {
            "subset_manifest": (
                isinstance(bench, Mapping)
                and subset_builder.verify_cached_subset(
                    bench, work, tools, threads=threads,
                )
            ),
            "liftoff_cache": tool_runners.verify_cached_tool_run(
                work, "liftoff",
            ),
            "miniprot_cache": tool_runners.verify_cached_tool_run(
                work, "miniprot",
            ),
        }
        missing = sorted(name for name, valid in verified.items() if not valid)
        if not missing:
            continue
        blockers.append({
            "code": "subset_inputs_not_prepared",
            "scenario_id": benchmark_id,
            "message": f"missing cached subset artifacts: {missing}",
        })
        actions.append([
            sys.executable,
            "-B",
            "-m",
            "benchmarks.compare.build_inputs",
            benchmark_id,
            "--registry",
            str(Path(benchmark_registry).resolve()),
            "--work-root",
            str(Path(work_root).resolve()),
            "--threads",
            str(threads),
            "--no-force",
        ])
    return blockers, actions


def _ortholog_build_actions(
    mapping_blockers: Sequence[Mapping[str, str]],
    scenario_inputs: Mapping[str, Mapping[str, Path]],
    runtime_root: Path,
    *,
    threads: int,
    excluded_scenarios: set[str] | None = None,
) -> list[list[str]]:
    actions = []
    excluded_scenarios = set(excluded_scenarios or ())
    for blocker in sorted(
        mapping_blockers,
        key=lambda item: item["scenario_id"],
    ):
        if blocker.get("code") != "missing_ortholog_map":
            continue
        scenario_id = blocker["scenario_id"]
        if scenario_id in excluded_scenarios:
            continue
        inputs = scenario_inputs[scenario_id]
        actions.append([
            sys.executable,
            "-B",
            "-m",
            "benchmarks.compare.ortholog_scope",
            "build",
            str(inputs["reference_annotation"]),
            str(inputs["reference_genome"]),
            str(inputs["target_truth"]),
            str(inputs["target_genome"]),
            str(runtime_root / "ortholog_scopes" / scenario_id),
            "--threads",
            str(threads),
            "--minimum-gene-groups",
            str(MINIMUM_FULL_TRUTH_GROUPS),
            "--minimum-transcript-groups",
            str(MINIMUM_FULL_TRUTH_GROUPS),
        ])
    return actions


def _ortholog_handoff(
    manifest: Mapping[str, Any],
    runtime_root: Path,
) -> dict[str, Any]:
    return {
        "schema_version": ORTHOLOG_SCHEMA_VERSION,
        "manifest_sha256": manifest_tools.canonical_sha256(manifest),
        "mappings": {
            scenario_id: {
                "path": str(
                    runtime_root
                    / "ortholog_scopes"
                    / scenario_id
                    / "ortholog_map.json"
                ),
                "sha256": "required-after-build",
                "id_policy": "ortholog-map",
            }
            for scenario_id in BIOLOGICAL_IDS
        },
        "instruction": (
            "Replace every required-after-build marker with the observed "
            "lowercase SHA-256. The prepare command rejects placeholders, "
            "missing scenarios, stale hashes, weak scope, or provenance drift."
        ),
    }


def materialize(
    manifest: Mapping[str, Any],
    lock: Mapping[str, Any],
    *,
    cache_root: Path,
    benchmark_registry: Path,
    dataset_registry: Path,
    ortholog_registry: Path | None = None,
    work_root: Path = DEFAULT_WORK_ROOT,
    data_root: Path = HERE / "data",
    threads: int = 8,
    synthetic_builder: SyntheticBuilder = _default_synthetic_builder,
) -> dict[str, Any]:
    """Prepare runtime assets and publish overlays only when truth is resolved."""

    manifest = manifest_tools.validate_manifest(manifest)
    cache_root = Path(cache_root).resolve()
    # Invalidate any previously ready preflight before touching derived
    # artifacts.  A killed or failed rematerialization must never leave an
    # older campaign-ready document authorizing bytes that may have been
    # replaced in place.
    _atomic_write_json(cache_root / PREFLIGHT_NAME, {
        "schema_version": 2,
        "manifest_sha256": manifest_tools.canonical_sha256(manifest),
        "acquisition_lock_sha256": manifest_tools.canonical_sha256(lock),
        "campaign_ready": False,
        "registries_exported": False,
        "blockers": [{
            "code": "materialization_in_progress",
            "scenario_id": "canonical-v2",
            "message": (
                "materialization has started but final derived artifacts and "
                "registry overlays have not been reverified"
            ),
        }],
        "remaining_actions": [],
    })
    runtime_root = cache_root / "runtime"
    base_benchmark = _load_benchmark_registry(benchmark_registry)
    try:
        base_dataset = json.loads(
            Path(dataset_registry).read_text(encoding="utf-8")
        )
    except (OSError, json.JSONDecodeError) as exc:
        raise AcquisitionError(f"cannot load dataset registry: {exc}") from exc
    if not isinstance(base_dataset.get("datasets"), list):
        raise AcquisitionError("dataset registry has no datasets list")

    source_paths = _locked_source_paths(
        manifest, lock, cache_root,
    )
    source_paths.update(_existing_source_paths(manifest, base_benchmark))
    runtime_paths, source_provenance = materialize_source_assets(
        manifest, source_paths, runtime_root,
    )
    synthetic, synthetic_provenance = generate_synthetic_scenarios(
        manifest,
        runtime_paths,
        runtime_root,
        builder=synthetic_builder,
    )
    resolved_scenarios = {}
    biological_inputs = {}
    blockers = []
    for scenario in manifest["scenarios"]:
        if scenario["kind"] != "biological":
            continue
        inputs = _resolve_scenario_inputs(scenario, runtime_paths)
        biological_inputs[scenario["id"]] = inputs
        resolved_scenarios[scenario["id"]] = {
            role: str(path.resolve()) for role, path in sorted(inputs.items())
        }
        blockers.extend(_seqid_blockers(scenario, inputs))
    mappings, mapping_blockers = _ortholog_mappings(
        ortholog_registry,
        manifest,
        biological_inputs,
    )
    blockers.extend(mapping_blockers)
    biological_rows = []
    for scenario in manifest["scenarios"]:
        if scenario["kind"] != "biological":
            continue
        inputs = biological_inputs[scenario["id"]]
        if scenario["id"] in mappings:
            biological_rows.append(
                _benchmark_row(scenario, inputs, mappings[scenario["id"]])
            )
    if blockers:
        actions = _ortholog_build_actions(
            mapping_blockers,
            biological_inputs,
            runtime_root,
            threads=threads,
            excluded_scenarios={
                blocker["scenario_id"]
                for blocker in blockers
                if blocker["code"] != "missing_ortholog_map"
            },
        )
        report = {
            "schema_version": 2,
            "manifest_sha256": manifest_tools.canonical_sha256(manifest),
            "acquisition_lock_sha256": manifest_tools.canonical_sha256(lock),
            "campaign_ready": False,
            "registries_exported": False,
            "blockers": sorted(
                blockers,
                key=lambda item: (
                    item["scenario_id"], item["code"], item["message"],
                ),
            ),
            "remaining_actions": actions,
            "ortholog_registry_handoff": _ortholog_handoff(
                manifest, runtime_root,
            ),
            "resolved_scenarios": resolved_scenarios,
            "source_provenance": source_provenance,
            "synthetic_provenance": synthetic_provenance,
        }
        _atomic_write_json(cache_root / PREFLIGHT_NAME, report)
        return report

    synthetic_rows = [
        _synthetic_benchmark_row(scenario, synthetic[scenario["id"]])
        for scenario in manifest["scenarios"]
        if scenario["kind"] == "synthetic"
    ]
    additions = biological_rows + synthetic_rows
    overlay_root = runtime_root / "registries"
    benchmark_preparation = overlay_root / "benchmarks.preparation.json"
    dataset_overlay = overlay_root / "datasets.json"
    benchmark_document = {
        **base_benchmark,
        "_canonical_v2": {
            "manifest_sha256": manifest_tools.canonical_sha256(manifest),
            "acquisition_lock_sha256": manifest_tools.canonical_sha256(lock),
        },
        "benchmarks": _merge_rows(
            base_benchmark["benchmarks"], additions,
        ),
    }
    e2e_ids = {
        scenario["id"]
        for scenario in manifest["scenarios"]
        if "cross_e2e" in scenario["panels"]
    }
    e2e_rows = [
        _e2e_dataset_row(row, Path(data_root))
        for row in biological_rows
        if row["id"] in e2e_ids
    ]
    dataset_document = {
        **base_dataset,
        "_canonical_v2": {
            "manifest_sha256": manifest_tools.canonical_sha256(manifest),
            "acquisition_lock_sha256": manifest_tools.canonical_sha256(lock),
        },
        "datasets": _merge_rows(base_dataset["datasets"], e2e_rows),
    }
    _atomic_write_json(benchmark_preparation, benchmark_document)
    _atomic_write_json(dataset_overlay, dataset_document)

    subset_ids = [
        scenario["id"]
        for scenario in manifest["scenarios"]
        if (
            scenario["kind"] in {"biological", "synthetic"}
            and (
                "subset" in scenario["panels"]
                or "synthetic" in scenario["panels"]
            )
        )
    ]
    prep_blockers, actions = _subset_preflight(
        subset_ids,
        benchmark_registry=benchmark_preparation,
        work_root=work_root,
        threads=threads,
    )
    prep_blockers.sort(
        key=lambda item: (
            item["scenario_id"], item["code"], item["message"],
        )
    )
    common_report = {
        "schema_version": 2,
        "manifest_sha256": manifest_tools.canonical_sha256(manifest),
        "acquisition_lock_sha256": manifest_tools.canonical_sha256(lock),
        "resolved_scenarios": resolved_scenarios,
        "source_provenance": source_provenance,
        "synthetic_provenance": synthetic_provenance,
        "ortholog_mappings": {
            scenario_id: {
                key: (
                    str(value.resolve()) if isinstance(value, Path) else value
                )
                for key, value in record.items()
            }
            for scenario_id, record in sorted(mappings.items())
        },
    }
    if prep_blockers or actions:
        report = {
            **common_report,
            "campaign_ready": False,
            "registries_exported": False,
            "preparation_registries_exported": True,
            "registries": {
                "benchmark_preparation": _file_record(benchmark_preparation),
                "dataset_preparation": _file_record(dataset_overlay),
            },
            "blockers": prep_blockers,
            "remaining_actions": actions,
        }
        _atomic_write_json(cache_root / PREFLIGHT_NAME, report)
        return report

    panel_truth_provenance = {}
    final_biological_rows = []
    for row in biological_rows:
        scenario_id = row["id"]
        subset_truth, provenance = _materialize_subset_truth(
            row,
            mappings[scenario_id],
            work_root=work_root,
            runtime_root=runtime_root,
        )
        final_row = dict(row)
        final_row["target_truth_by_panel"] = {
            "subset": subset_truth,
            "full": dict(row["target_truth"]),
        }
        final_biological_rows.append(final_row)
        panel_truth_provenance[scenario_id] = provenance

    final_additions = final_biological_rows + synthetic_rows
    benchmark_overlay = overlay_root / "benchmarks.json"
    final_benchmark_document = {
        **base_benchmark,
        "_canonical_v2": {
            "manifest_sha256": manifest_tools.canonical_sha256(manifest),
            "acquisition_lock_sha256": manifest_tools.canonical_sha256(lock),
        },
        "benchmarks": _merge_rows(
            base_benchmark["benchmarks"], final_additions,
        ),
    }
    final_e2e_rows = [
        _e2e_dataset_row(row, Path(data_root))
        for row in final_biological_rows
        if row["id"] in e2e_ids
    ]
    final_dataset_document = {
        **base_dataset,
        "_canonical_v2": {
            "manifest_sha256": manifest_tools.canonical_sha256(manifest),
            "acquisition_lock_sha256": manifest_tools.canonical_sha256(lock),
        },
        "datasets": _merge_rows(base_dataset["datasets"], final_e2e_rows),
    }
    _atomic_write_json(benchmark_overlay, final_benchmark_document)
    _atomic_write_json(dataset_overlay, final_dataset_document)
    registry_records = {
        "benchmark": _file_record(benchmark_overlay),
        "dataset": _file_record(dataset_overlay),
    }
    if ortholog_registry is not None:
        registry_records["ortholog"] = _file_record(
            Path(ortholog_registry).resolve()
        )
    report = {
        **common_report,
        "campaign_ready": True,
        "registries_exported": True,
        "preparation_registries_exported": True,
        "registries": registry_records,
        "blockers": [],
        "remaining_actions": [],
        "panel_truth_provenance": panel_truth_provenance,
    }
    _atomic_write_json(cache_root / PREFLIGHT_NAME, report)
    return report


def _load_manifest(path: Path) -> dict[str, Any]:
    return manifest_tools.validate_manifest(manifest_tools.load_json(path))


def _bytecode_safe_command(argv: Sequence[str]) -> list[str]:
    command = [str(value) for value in argv]
    if (
        command
        and Path(command[0]).resolve() == Path(sys.executable).resolve()
        and "-B" not in command[1:3]
    ):
        command.insert(1, "-B")
    return command


def _bootstrap_ortholog_state(
    manifest: Mapping[str, Any],
    *,
    manifest_path: Path,
    ortholog_root: Path,
    final_registry: Path,
) -> dict[str, Any]:
    """Validate completed bundles and select a full or partial registry."""

    from benchmarks.compare import ortholog_scope

    expected = tuple(sorted(
        scenario["id"] for scenario in manifest["scenarios"]
        if scenario["kind"] == "biological"
    ))
    if not expected:
        raise AcquisitionError("bootstrap manifest has no biological scenarios")
    valid = {}
    missing = []
    invalid = []
    ortholog_root = Path(ortholog_root).resolve()
    final_registry = Path(final_registry).resolve()
    partial_registry = ortholog_root / "ortholog_registry.partial.json"
    if ortholog_root.is_dir():
        allowed_names = {
            *expected,
            final_registry.name,
            partial_registry.name,
        }
        for child in sorted(ortholog_root.iterdir()):
            if child.name not in allowed_names and not any(
                child.name.startswith(f".{scenario_id}.tmp-")
                for scenario_id in expected
            ):
                invalid.append((child, "unexpected bootstrap artifact"))
    for scenario_id in expected:
        bundle_path = ortholog_root / scenario_id
        interrupted = sorted(
            ortholog_root.glob(f".{scenario_id}.tmp-*")
        )
        if interrupted:
            invalid.append((interrupted[0], "interrupted temporary bundle"))
            continue
        if not bundle_path.exists() and not bundle_path.is_symlink():
            missing.append(scenario_id)
            continue
        if bundle_path.is_symlink() or not bundle_path.is_dir():
            invalid.append((bundle_path, "bundle path is not a real directory"))
            continue
        try:
            bundle = ortholog_scope.validate_scope_bundle(bundle_path)
        except (OSError, ValueError, ortholog_scope.ScopeBuildError) as exc:
            invalid.append((bundle_path, str(exc)))
            continue
        mapping_path = Path(bundle["mapping_path"])
        valid[scenario_id] = {
            "id_policy": "ortholog-map",
            "path": Path(os.path.relpath(
                mapping_path,
                start=final_registry.parent,
            )).as_posix(),
            "sha256": bundle["mapping_sha256"],
        }
    if invalid:
        path, reason = invalid[0]
        raise AcquisitionError(
            "invalid partial ortholog scope; quarantine the path outside "
            f"{ortholog_root} and rerun bootstrap: {path}: {reason}"
        )

    expected_registry = {
        "schema_version": ORTHOLOG_SCHEMA_VERSION,
        "manifest_sha256": manifest_tools.canonical_sha256(manifest),
        "mappings": dict(sorted(valid.items())),
    }
    if final_registry.exists() or final_registry.is_symlink():
        if final_registry.is_symlink() or not final_registry.is_file():
            raise AcquisitionError(
                "final ortholog registry is not a regular file; quarantine it "
                f"outside {ortholog_root} and rerun bootstrap"
            )
        if missing:
            raise AcquisitionError(
                "final ortholog registry exists but scenario bundles are "
                "incomplete; quarantine the registry and affected scope "
                f"outside {ortholog_root}, then rerun bootstrap"
            )
        try:
            observed_registry = json.loads(
                final_registry.read_text(encoding="utf-8")
            )
        except (OSError, json.JSONDecodeError) as exc:
            raise AcquisitionError(
                "final ortholog registry is unreadable; quarantine it outside "
                f"{ortholog_root} and rerun bootstrap: {exc}"
            ) from exc
        if observed_registry != expected_registry:
            raise AcquisitionError(
                "final ortholog registry is stale or tampered; quarantine it "
                f"outside {ortholog_root} and rerun bootstrap"
            )
        registry_path: Path | None = final_registry
    elif not missing:
        ortholog_root.mkdir(parents=True, exist_ok=True)
        try:
            ortholog_scope.finalize_mapping_registry(
                manifest_path,
                ortholog_root,
                final_registry,
            )
        except (OSError, ValueError, ortholog_scope.ScopeBuildError) as exc:
            raise AcquisitionError(
                f"cannot finalize validated ortholog bundles: {exc}"
            ) from exc
        observed_registry = manifest_tools.load_json(final_registry)
        if observed_registry != expected_registry:
            raise AcquisitionError(
                "finalized ortholog registry does not match validated bundles"
            )
        registry_path = final_registry
    elif valid:
        ortholog_root.mkdir(parents=True, exist_ok=True)
        registry_path = partial_registry
        _atomic_write_json(registry_path, expected_registry)
    else:
        registry_path = None
    return {
        "valid_scenarios": sorted(valid),
        "missing_scenarios": sorted(missing),
        "registry_path": registry_path,
    }


def _ortholog_action_scenario(action: Sequence[str]) -> str | None:
    command = [str(value) for value in action]
    try:
        module_index = command.index("benchmarks.compare.ortholog_scope")
    except ValueError:
        return None
    if command[module_index + 1:module_index + 2] != ["build"]:
        return None
    output_index = module_index + 6
    if output_index >= len(command):
        raise AcquisitionError(f"malformed ortholog build action: {command}")
    return Path(command[output_index]).name


def _resume_safe_actions(
    actions: Sequence[Sequence[str]],
    ortholog_state: Mapping[str, Any],
) -> list[list[str]]:
    """Drop completed ortholog actions and reject undeclared scenarios."""

    valid = set(ortholog_state["valid_scenarios"])
    missing = set(ortholog_state["missing_scenarios"])
    selected = []
    for raw in actions:
        if (
            isinstance(raw, (str, bytes))
            or not isinstance(raw, Sequence)
            or not raw
        ):
            raise AcquisitionError(f"malformed bootstrap action: {raw!r}")
        action = [str(value) for value in raw]
        scenario_id = _ortholog_action_scenario(action)
        if scenario_id is None:
            selected.append(action)
        elif scenario_id in missing:
            selected.append(action)
        elif scenario_id not in valid:
            raise AcquisitionError(
                f"ortholog action targets undeclared scenario {scenario_id!r}"
            )
    return selected


def _run_bootstrap_actions(
    actions: Sequence[Sequence[str]],
    *,
    cache_root: Path,
    max_active: int,
) -> list[dict[str, Any]]:
    if max_active < 1 or max_active > 4:
        raise AcquisitionError("bootstrap max-active must be between 1 and 4")
    action_root = Path(cache_root) / "bootstrap" / "actions"
    action_root.mkdir(parents=True, exist_ok=True)
    environment = os.environ.copy()
    environment.update({
        "PYTHONDONTWRITEBYTECODE": "1",
        "PYTHONNOUSERSITE": "1",
        "PYTHONHASHSEED": "0",
    })

    def run(index: int, raw: Sequence[str]) -> dict[str, Any]:
        command = _bytecode_safe_command(raw)
        fingerprint = manifest_tools.canonical_sha256(command)[:16]
        prefix = action_root / f"{index:03d}-{fingerprint}"
        stdout_path = prefix.with_suffix(".stdout.log")
        stderr_path = prefix.with_suffix(".stderr.log")
        with stdout_path.open("w", encoding="utf-8") as stdout, \
                stderr_path.open("w", encoding="utf-8") as stderr:
            completed = subprocess.run(
                command,
                cwd=HERE.parent,
                env=environment,
                stdout=stdout,
                stderr=stderr,
                check=False,
            )
        result = {
            "schema_version": 1,
            "index": index,
            "command": command,
            "fingerprint": fingerprint,
            "returncode": completed.returncode,
            "stdout": str(stdout_path.resolve()),
            "stderr": str(stderr_path.resolve()),
        }
        _atomic_write_json(prefix.with_suffix(".result.json"), result)
        return result

    results = []
    with ThreadPoolExecutor(max_workers=max_active) as executor:
        futures = {
            executor.submit(run, index, action): index
            for index, action in enumerate(actions, start=1)
        }
        for future in as_completed(futures):
            results.append(future.result())
    results.sort(key=lambda item: item["index"])
    failures = [item for item in results if item["returncode"] != 0]
    if failures:
        first = failures[0]
        raise AcquisitionError(
            "bootstrap action failed; inspect "
            f"{first['stderr']}: {first['command']}"
        )
    return results


def bootstrap(
    manifest: Mapping[str, Any],
    *,
    manifest_path: Path,
    cache_root: Path,
    lock_path: Path | None,
    benchmark_registry: Path,
    dataset_registry: Path,
    threads: int,
    max_active: int,
    accept_identity_pinned_bytes: bool,
) -> dict[str, Any]:
    """Acquire and resumably drive every preparation action to readiness."""

    if threads < 1 or threads > 32:
        raise AcquisitionError("bootstrap threads must be between 1 and 32")
    cache_root = _safe_cache_root(cache_root)
    lock_path = _safe_lock_path(lock_path, cache_root)
    lock, verification = acquire(
        manifest,
        cache_root=cache_root,
        lock_path=lock_path,
        accept_identity_pinned_bytes=accept_identity_pinned_bytes,
    )
    runtime_root = cache_root / "runtime"
    work_root = cache_root / "work"
    ortholog_root = runtime_root / "ortholog_scopes"
    ortholog_registry = ortholog_root / "ortholog_registry.json"
    action_rounds = []
    for round_number in range(1, 6):
        ortholog_state = _bootstrap_ortholog_state(
            manifest,
            manifest_path=manifest_path,
            ortholog_root=ortholog_root,
            final_registry=ortholog_registry,
        )
        report = materialize(
            manifest,
            lock,
            cache_root=cache_root,
            benchmark_registry=benchmark_registry,
            dataset_registry=dataset_registry,
            ortholog_registry=ortholog_state["registry_path"],
            work_root=work_root,
            data_root=runtime_root / "datasets",
            threads=threads,
        )
        if report["campaign_ready"]:
            return {
                "schema_version": 1,
                "campaign_ready": True,
                "manifest": _file_record(manifest_path),
                "lock": _file_record(lock_path),
                "lock_verification": verification,
                "preflight": _file_record(cache_root / PREFLIGHT_NAME),
                "action_rounds": action_rounds,
                "ortholog_bootstrap": {
                    "valid_scenarios": ortholog_state["valid_scenarios"],
                    "missing_scenarios": ortholog_state["missing_scenarios"],
                    "registry": _file_record(ortholog_registry),
                },
                "registries": report["registries"],
            }
        raw_actions = report.get("remaining_actions")
        actions = (
            _resume_safe_actions(raw_actions, ortholog_state)
            if isinstance(raw_actions, list)
            else None
        )
        if not isinstance(actions, list) or not actions:
            raise AcquisitionError(
                "bootstrap cannot resolve acquisition blockers: "
                + json.dumps(report.get("blockers", []), sort_keys=True)
            )
        results = _run_bootstrap_actions(
            actions, cache_root=cache_root, max_active=max_active,
        )
        action_rounds.append({
            "round": round_number,
            "actions": results,
        })
    raise AcquisitionError("bootstrap exceeded five preparation rounds")


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--manifest", type=Path, default=manifest_tools.DEFAULT_MANIFEST,
    )
    parser.add_argument("--cache-root", type=Path, default=DEFAULT_CACHE_ROOT)
    subparsers = parser.add_subparsers(dest="command", required=True)

    acquire_parser = subparsers.add_parser(
        "acquire", help="download and lock all remote source bytes",
    )
    acquire_parser.add_argument("--lock", type=Path)
    acquire_parser.add_argument("--dry-run", action="store_true")
    acquire_parser.add_argument(
        "--accept-identity-pinned-bytes",
        action="store_true",
        help=(
            "explicitly lock first-observed bytes for identities whose provider "
            "does not declare a checksum"
        ),
    )

    prepare = subparsers.add_parser(
        "prepare", help="materialize assets and export registry overlays",
    )
    prepare.add_argument("--lock", type=Path)
    prepare.add_argument(
        "--benchmark-registry",
        type=Path,
        default=DEFAULT_BENCHMARK_REGISTRY,
    )
    prepare.add_argument(
        "--dataset-registry",
        type=Path,
        default=DEFAULT_DATASET_REGISTRY,
    )
    prepare.add_argument("--ortholog-registry", type=Path)
    prepare.add_argument("--work-root", type=Path)
    prepare.add_argument("--threads", type=int, default=8)

    bootstrap_parser = subparsers.add_parser(
        "bootstrap",
        help="resumably acquire, prepare, validate, and export final overlays",
    )
    bootstrap_parser.add_argument("--lock", type=Path)
    bootstrap_parser.add_argument(
        "--benchmark-registry",
        type=Path,
        default=DEFAULT_BENCHMARK_REGISTRY,
    )
    bootstrap_parser.add_argument(
        "--dataset-registry",
        type=Path,
        default=DEFAULT_DATASET_REGISTRY,
    )
    bootstrap_parser.add_argument("--threads", type=int, default=8)
    bootstrap_parser.add_argument("--max-active", type=int, default=2)
    bootstrap_parser.add_argument(
        "--accept-identity-pinned-bytes",
        action="store_true",
        help="acknowledge and lock source bytes lacking provider checksums",
    )
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    arguments = _parser().parse_args(argv)
    try:
        cache_root = _safe_cache_root(arguments.cache_root)
        manifest = _load_manifest(arguments.manifest)
        if arguments.command == "acquire":
            if arguments.dry_run:
                result = acquisition_dry_run(
                    manifest,
                    cache_root,
                    lock_path=arguments.lock,
                )
            else:
                lock, verification = acquire(
                    manifest,
                    cache_root=cache_root,
                    lock_path=arguments.lock,
                    accept_identity_pinned_bytes=(
                        arguments.accept_identity_pinned_bytes
                    ),
                )
                result = {
                    "lock": str(
                        Path(arguments.lock or cache_root / LOCK_NAME).resolve()
                    ),
                    "lock_sha256": manifest_tools.canonical_sha256(lock),
                    "verification": verification,
                }
        elif arguments.command == "prepare":
            lock_path = _safe_lock_path(arguments.lock, cache_root)
            work_root = Path(
                arguments.work_root or cache_root / "work"
            ).expanduser().resolve()
            if work_root == cache_root or not work_root.is_relative_to(cache_root):
                raise AcquisitionError(
                    "prepare work root must be a child of cache root "
                    f"{cache_root}: {work_root}"
                )
            report = materialize(
                manifest,
                manifest_tools.load_json(lock_path),
                cache_root=cache_root,
                benchmark_registry=arguments.benchmark_registry,
                dataset_registry=arguments.dataset_registry,
                ortholog_registry=arguments.ortholog_registry,
                work_root=work_root,
                data_root=cache_root / "runtime" / "datasets",
                threads=arguments.threads,
            )
            result = report
        else:
            result = bootstrap(
                manifest,
                manifest_path=Path(arguments.manifest).resolve(),
                cache_root=cache_root,
                lock_path=arguments.lock,
                benchmark_registry=arguments.benchmark_registry,
                dataset_registry=arguments.dataset_registry,
                threads=arguments.threads,
                max_active=arguments.max_active,
                accept_identity_pinned_bytes=(
                    arguments.accept_identity_pinned_bytes
                ),
            )
        print(json.dumps(result, indent=2, sort_keys=True))
        if arguments.command in {"prepare", "bootstrap"} and not result[
            "campaign_ready"
        ]:
            return 3
        return 0
    except (
        AcquisitionError,
        manifest_tools.ManifestError,
        OSError,
        ValueError,
    ) as exc:
        print(f"canonical-v2 acquisition error: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
