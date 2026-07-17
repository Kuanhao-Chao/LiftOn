"""Pinned, paired LiftOn release evaluation helpers.

This module compares two exact LiftOn source trees while deliberately sharing
the Python environment and external binaries.  It is designed to be called by
``build_controller`` workers: one cell runs one AB/BA pair, writes all outputs
below its cell directory, validates both GFF3 files with the current validator,
and re-scores both annotations with the current neutral evaluator.

The canonical benchmark baseline is read-only.  Paired results are independent
JSON documents that can be aggregated by :mod:`release_report`.
"""
from __future__ import annotations

import argparse
import contextlib
import dataclasses
import datetime as dt
import hashlib
import importlib.machinery
import json
import math
import os
import re
import shutil
import statistics
import subprocess
import sys
import urllib.parse
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping, Sequence

from lifton.gff3_validator import validate_gff3_file

from benchmarks import run_benchmarks

from . import devel_refresh
from . import evaluator
from . import fourway_compare
from .profiling import run_profiled


HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parent.parent
DEFAULT_BENCHMARK_REGISTRY = HERE / "benchmarks.json"
DEFAULT_REFERENCE_SHA = "3739dfc8f73396fccab7d7be0f95e008179cea5d"
SCHEMA_VERSION = 3
DEFAULT_BOOTSTRAP_SEED = 20260717
DEFAULT_BOOTSTRAP_REPLICATES = 10_000
E2E_MODE_FEATURES = {
    "safe": frozenset(),
    "stream": frozenset({"--stream"}),
    "inmemory": frozenset({"--inmemory-liftoff"}),
    "native": frozenset({"--native"}),
    "stream-inmemory": frozenset({"--stream", "--inmemory-liftoff"}),
    "stream-native": frozenset({"--stream", "--native"}),
    "inmemory-native": frozenset({"--inmemory-liftoff", "--native"}),
    "fast": frozenset({"--stream", "--inmemory-liftoff", "--native"}),
}
E2E_MODES = tuple(E2E_MODE_FEATURES)
SEMANTIC_HASH_ALGORITHM = "sha256-multiset-sum2-v1"
_SEMANTIC_MODULUS = 1 << 256
STABLE_ID_FEATURE_TYPES = ("CDS", "exon")
STABLE_ID_METHOD = "declared-gff3-id-same-type-copy-aware-v2"
_COPY_SUFFIX_RE = re.compile(r"_(\d+)$")


@dataclass(frozen=True)
class SourceSpec:
    """One exact LiftOn implementation executed in the shared environment."""

    label: str
    root: Path
    sha: str
    lifton_executable: Path

    def environment(self) -> dict[str, str]:
        current = os.environ.get("PYTHONPATH", "")
        pythonpath = str(self.root) if not current else f"{self.root}{os.pathsep}{current}"
        current_path = os.environ.get("PATH", "")
        executable_dir = str(self.lifton_executable.parent)
        path = (
            executable_dir
            if not current_path
            else f"{executable_dir}{os.pathsep}{current_path}"
        )
        return {
            "PATH": path,
            "PYTHONPATH": pythonpath,
            "PYTHONDONTWRITEBYTECODE": "1",
            "PYTHONNOUSERSITE": "1",
            "PYTHONHASHSEED": "0",
        }


@dataclass(frozen=True)
class PanelInputs:
    """Resolved, immutable inputs for one benchmark cell."""

    benchmark: str
    panel: str
    species: str
    cross_species: bool
    annotation_database: str
    ref_fa: Path
    ref_gff: Path
    tgt_fa: Path
    liftoff_gff: Path | None = None
    miniprot_gff: Path | None = None
    transcripts_fa: Path | None = None
    proteins_fa: Path | None = None
    truth_gff: Path | None = None

    def required_paths(self) -> dict[str, Path]:
        paths = {
            "ref_fa": self.ref_fa,
            "ref_gff": self.ref_gff,
            "tgt_fa": self.tgt_fa,
        }
        optional = {
            "liftoff_gff": self.liftoff_gff,
            "miniprot_gff": self.miniprot_gff,
            "transcripts_fa": self.transcripts_fa,
            "proteins_fa": self.proteins_fa,
        }
        paths.update({name: path for name, path in optional.items() if path is not None})
        return paths


def _git(root: Path, *arguments: str) -> str:
    result = subprocess.run(
        ["git", "-C", str(root), *arguments],
        text=True,
        capture_output=True,
        check=False,
    )
    if result.returncode != 0:
        raise RuntimeError(
            f"git {' '.join(arguments)} failed in {root}: {result.stderr.strip()}"
        )
    return result.stdout.strip()


def _git_other_paths(root: Path, *, ignored: bool) -> list[str]:
    arguments = ["ls-files", "-z", "--others", "--exclude-standard"]
    if ignored:
        arguments.insert(3, "--ignored")
    pathspecs = sorted({
        *(f"*{suffix}" for suffix in importlib.machinery.all_suffixes()),
        "*.pth",
        "*.egg-link",
        "sitecustomize",
        "usercustomize",
    })
    arguments.extend(["--", *pathspecs])
    return sorted(path for path in _git(root, *arguments).split("\0") if path)


def _can_affect_python_import(path: str) -> bool:
    suffixes = tuple(importlib.machinery.all_suffixes())
    name = Path(path).name
    return (
        name.endswith(suffixes)
        or name.endswith((".pth", ".egg-link"))
        or name in {"sitecustomize", "usercustomize"}
    )


def _unversioned_import_paths(root: Path) -> list[dict[str, str]]:
    paths: list[dict[str, str]] = []
    for category, ignored in (("untracked", False), ("ignored", True)):
        paths.extend(
            {"path": path, "category": category}
            for path in _git_other_paths(root, ignored=ignored)
            if _can_affect_python_import(path)
        )
    return sorted(paths, key=lambda item: (item["path"], item["category"]))


def verify_source(spec: SourceSpec) -> dict[str, Any]:
    """Fail closed if a source root or its import path is not the pinned SHA."""

    root = spec.root.resolve()
    if not root.is_dir():
        raise RuntimeError(f"{spec.label} source root does not exist: {root}")
    actual_sha = _git(root, "rev-parse", "HEAD")
    if actual_sha != spec.sha:
        raise RuntimeError(
            f"{spec.label} source SHA is {actual_sha}, expected {spec.sha}"
        )
    dirty = _git(root, "status", "--porcelain", "--untracked-files=no")
    if dirty:
        raise RuntimeError(f"{spec.label} source tree has tracked changes")
    unversioned_imports = _unversioned_import_paths(root)
    if unversioned_imports:
        examples = ", ".join(
            f"{item['category']}:{item['path']}"
            for item in unversioned_imports[:20]
        )
        if len(unversioned_imports) > 20:
            examples += f", ... ({len(unversioned_imports)} total)"
        raise RuntimeError(
            f"{spec.label} source tree has unversioned import-affecting "
            f"files: {examples}"
        )
    command = [
        str(spec.lifton_executable.parent / "python"),
        "-c",
        "import lifton; print(lifton.__file__)",
    ]
    probe = subprocess.run(
        command,
        cwd="/tmp",
        env={**os.environ, **spec.environment()},
        text=True,
        capture_output=True,
        check=False,
    )
    imported = Path(probe.stdout.strip()).resolve() if probe.returncode == 0 else None
    if probe.returncode != 0 or imported is None or not imported.is_relative_to(root):
        raise RuntimeError(
            f"{spec.label} import probe resolved to {imported}; stderr={probe.stderr.strip()}"
        )
    return {
        "label": spec.label,
        "root": str(root),
        "sha": actual_sha,
        "lifton_executable": str(spec.lifton_executable),
        "imported_package": str(imported),
    }


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _canonical_attributes(value: str) -> str:
    fields = [field for field in value.split(";") if field]
    normalized = []
    for field in fields:
        key, separator, raw = field.partition("=")
        if not separator:
            normalized.append((key, ""))
            continue
        values = ",".join(sorted(part.strip() for part in raw.split(",")))
        normalized.append((key.strip(), values))
    return ";".join(f"{key}={raw}" if raw else key for key, raw in sorted(normalized))


class _SemanticMultisetHash:
    """Constant-row-memory commutative digest with duplicate multiplicity."""

    __slots__ = ("count", "_sum", "_sum_squares")

    def __init__(self) -> None:
        self.count = 0
        self._sum = 0
        self._sum_squares = 0

    def add(self, row: str) -> None:
        digest = hashlib.sha256()
        digest.update(b"lifton-gff3-semantic-row-v1\0")
        digest.update(row.encode("utf-8", "surrogateescape"))
        value = int.from_bytes(digest.digest(), "big")
        self.count += 1
        self._sum = (self._sum + value) % _SEMANTIC_MODULUS
        self._sum_squares = (
            self._sum_squares + (value * value)
        ) % _SEMANTIC_MODULUS

    def hexdigest(self) -> str:
        digest = hashlib.sha256()
        digest.update(SEMANTIC_HASH_ALGORITHM.encode("ascii"))
        digest.update(b"\0")
        digest.update(str(self.count).encode("ascii"))
        digest.update(b"\0")
        digest.update(self._sum.to_bytes(32, "big"))
        digest.update(self._sum_squares.to_bytes(32, "big"))
        return digest.hexdigest()


def gff3_fingerprints(path: Path) -> dict[str, Any]:
    """Return byte and order-independent semantic hashes plus feature counts."""

    byte_sha = sha256_file(path)
    semantic = _SemanticMultisetHash()
    feature_counts: Counter[str] = Counter()
    with Path(path).open(encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            line = raw.rstrip("\r\n")
            if not line or line.startswith("#"):
                continue
            columns = line.split("\t")
            if len(columns) != 9:
                semantic.add(line)
                continue
            columns[8] = _canonical_attributes(columns[8])
            feature_counts[columns[2]] += 1
            semantic.add("\t".join(columns))
    return {
        "byte_sha256": byte_sha,
        "semantic_sha256": semantic.hexdigest(),
        "semantic_algorithm": SEMANTIC_HASH_ALGORITHM,
        "feature_records": semantic.count,
        "feature_counts": dict(sorted(feature_counts.items())),
    }


def _declared_id_index(
    path: Path,
    feature_types: Sequence[str] = STABLE_ID_FEATURE_TYPES,
) -> dict[str, dict[str, Any]]:
    """Index explicit column-9 IDs without inventing database row IDs.

    Discontinuous CDS segments legitimately repeat one ID, so the preservation
    denominator is the set of declared logical IDs rather than the row count.
    """

    selected = tuple(feature_types)
    index = {
        feature_type: {
            "n_records": 0,
            "n_records_with_id": 0,
            "ids": set(),
            "parents_by_id": {},
        }
        for feature_type in selected
    }
    with Path(path).open(encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            if not raw or raw.startswith("#"):
                continue
            columns = raw.rstrip("\r\n").split("\t")
            if len(columns) != 9 or columns[2] not in index:
                continue
            row = index[columns[2]]
            row["n_records"] += 1
            declared_id = None
            declared_parents: set[str] = set()
            for field in columns[8].split(";"):
                key, separator, value = field.strip().partition("=")
                if not separator:
                    continue
                key = key.strip()
                if key == "ID" and declared_id is None:
                    # ID is scalar in GFF3. Taking the first non-empty value
                    # keeps malformed comma lists from inflating the denominator.
                    declared_id = next(
                        (
                            urllib.parse.unquote(item.strip())
                            for item in value.split(",")
                            if item.strip()
                        ),
                        None,
                    )
                elif key == "Parent":
                    declared_parents.update(
                        urllib.parse.unquote(item.strip())
                        for item in value.split(",")
                        if item.strip()
                    )
            if declared_id:
                row["n_records_with_id"] += 1
                row["ids"].add(declared_id)
                row["parents_by_id"].setdefault(declared_id, set()).update(
                    declared_parents
                )
    return index


def stable_id_preservation(
    reference_gff: Path,
    output_gff: Path,
    *,
    reference_index: Mapping[str, Mapping[str, Any]] | None = None,
) -> dict[str, Any]:
    """Measure declared CDS/exon ID continuity, not biological completeness."""

    reference = (
        dict(reference_index)
        if reference_index is not None
        else _declared_id_index(reference_gff)
    )
    output = _declared_id_index(output_gff)
    by_type: dict[str, dict[str, Any]] = {}
    for feature_type in STABLE_ID_FEATURE_TYPES:
        reference_row = reference[feature_type]
        output_row = output[feature_type]
        reference_ids = set(reference_row["ids"])
        output_ids = set(output_row["ids"])
        reference_parents = reference_row["parents_by_id"]
        output_parents = output_row["parents_by_id"]
        preserved: set[str] = set()
        for output_id in output_ids:
            if output_id in reference_ids:
                preserved.add(output_id)
                continue
            match = _COPY_SUFFIX_RE.search(output_id)
            if match:
                base = output_id[:match.start()]
                suffix = match.group(0)
                copy_parent_matches = any(
                    output_parent.endswith(suffix)
                    and output_parent[:-len(suffix)]
                    in reference_parents.get(base, set())
                    for output_parent in output_parents.get(output_id, set())
                )
                if base in reference_ids and copy_parent_matches:
                    preserved.add(base)
        n_reference_ids = len(reference_ids)
        applicable = n_reference_ids > 0
        if applicable:
            reason = None
        elif reference_row["n_records"] == 0:
            reason = "reference_feature_type_absent"
        else:
            reason = "no_declared_reference_ids"
        by_type[feature_type] = {
            "applicable": applicable,
            "reason": reason,
            "n_reference_records": reference_row["n_records"],
            "n_reference_records_with_id": (
                reference_row["n_records_with_id"]
            ),
            "n_reference_ids": n_reference_ids,
            "n_preserved_ids": len(preserved),
            "n_output_records": output_row["n_records"],
            "n_output_records_with_id": output_row["n_records_with_id"],
            "n_output_ids": len(output_ids),
            "preservation_rate": (
                len(preserved) / n_reference_ids
                if applicable else None
            ),
        }
    return {
        "method": STABLE_ID_METHOD,
        "by_type": by_type,
    }


def _require_nonempty(paths: Mapping[str, Path]) -> None:
    missing = [
        f"{name}={path}" for name, path in paths.items()
        if not path.is_file() or path.stat().st_size == 0
    ]
    if missing:
        raise RuntimeError("missing or empty paired input: " + ", ".join(missing))


def _benchmark_entry(path: Path, benchmark: str) -> dict[str, Any]:
    raw = json.loads(Path(path).read_text())
    entries = raw.get("benchmarks", []) if isinstance(raw, Mapping) else []
    matches = [
        item for item in entries
        if isinstance(item, Mapping) and item.get("id") == benchmark
    ]
    if len(matches) != 1:
        raise ValueError(
            f"expected one benchmark {benchmark!r} in {Path(path).resolve()}"
        )
    return dict(matches[0])


def resolve_panel_inputs(
    panel: str,
    benchmark: str,
    *,
    benchmark_registry: Path = DEFAULT_BENCHMARK_REGISTRY,
    dataset_registry: Path = run_benchmarks.DEFAULT_REGISTRY,
) -> PanelInputs:
    """Resolve canonical cached inputs without creating or refreshing them."""

    if panel == "e2e":
        registry = run_benchmarks.load_registry(
            Path(dataset_registry),
        )
        matches = [item for item in registry.datasets if item.id == benchmark]
        if len(matches) != 1:
            raise ValueError(f"unknown end-to-end dataset: {benchmark}")
        dataset = matches[0]
        data_dir = run_benchmarks.DEFAULT_DATA_DIR / benchmark
        result = PanelInputs(
            benchmark=benchmark,
            panel=panel,
            species=dataset.species,
            cross_species=False,
            annotation_database="RefSeq",
            ref_fa=data_dir / run_benchmarks._filename_for(dataset.reference_fa),
            ref_gff=data_dir / run_benchmarks._filename_for(dataset.reference_gff),
            tgt_fa=data_dir / run_benchmarks._filename_for(dataset.target_fa),
        )
        _require_nonempty(result.required_paths())
        return result

    benchmark_registry = Path(benchmark_registry).resolve()
    bench = _benchmark_entry(benchmark_registry, benchmark)
    common = {
        "benchmark": benchmark,
        "panel": panel,
        "species": str(bench["species"]),
        "cross_species": bool(bench["cross_species"]),
        "annotation_database": str(bench.get("annotation_database", "RefSeq")),
    }
    if panel == "subset":
        work = fourway_compare.WORK / benchmark
        manifest_path = work / "subset" / "subset.manifest.json"
        liftoff_gff = work / "tools" / "liftoff" / "liftoff.gff3"
        miniprot_gff = work / "tools" / "miniprot" / "miniprot.gff3"
        _require_nonempty({
            "subset_manifest": manifest_path,
            "liftoff_gff": liftoff_gff,
            "miniprot_gff": miniprot_gff,
        })
        manifest = json.loads(manifest_path.read_text())
        paths = manifest["paths"]
        result = PanelInputs(
            **common,
            ref_fa=Path(paths["ref_fa"]),
            ref_gff=Path(paths["ref_gff"]),
            tgt_fa=Path(paths["tgt_fa"]),
            liftoff_gff=Path(liftoff_gff),
            miniprot_gff=Path(miniprot_gff),
        )
    elif panel == "full":
        if benchmark_registry == DEFAULT_BENCHMARK_REGISTRY.resolve():
            paths = fourway_compare._full_paths(benchmark)
        else:
            paths = {
                "ref_fa": Path(bench["ref_genome"]),
                "ref_gff": Path(bench["ref_gff"]),
                "tgt_fa": Path(bench["tgt_genome"]),
            }
        liftoff_gff, miniprot_gff, transcripts_fa, proteins_fa = (
            devel_refresh._cached_inputs(benchmark)
        )
        result = PanelInputs(
            **common,
            ref_fa=Path(paths["ref_fa"]),
            ref_gff=Path(paths["ref_gff"]),
            tgt_fa=Path(paths["tgt_fa"]),
            liftoff_gff=Path(liftoff_gff),
            miniprot_gff=Path(miniprot_gff),
            transcripts_fa=Path(transcripts_fa),
            proteins_fa=Path(proteins_fa),
        )
    else:
        raise ValueError(
            f"paired panel must be subset, full, or e2e, not {panel!r}"
        )
    _require_nonempty(result.required_paths())
    return result


def input_fingerprints(inputs: PanelInputs) -> dict[str, Any]:
    return {
        name: {
            "path": str(path.resolve()),
            "size": path.stat().st_size,
            "sha256": sha256_file(path),
        }
        for name, path in sorted(inputs.required_paths().items())
    }


def _isolated_link(source: Path, destination: Path) -> Path:
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists() or destination.is_symlink():
        raise FileExistsError(f"refusing to replace paired input link: {destination}")
    destination.symlink_to(source.resolve())
    return destination


def _isolated_fasta(source: Path, destination: Path) -> Path:
    """Link a read-only FASTA and copy mutable indexes into the cell."""

    linked = _isolated_link(source, destination)
    sidecars = (
        (Path(str(source) + ".fai"), Path(str(destination) + ".fai")),
        (Path(str(source) + ".gzi"), Path(str(destination) + ".gzi")),
        (source.with_suffix(".dict"), destination.with_suffix(".dict")),
    )
    for canonical, isolated in sidecars:
        if not canonical.is_file():
            continue
        isolated.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(canonical, isolated)
    return linked


def _isolated_inputs(inputs: PanelInputs, root: Path) -> PanelInputs:
    """Give one version private database/index sidecar names via symlinks."""

    linked = {
        name: _isolated_link(path, root / "inputs" / f"{name}{path.suffix}")
        for name, path in inputs.required_paths().items()
    }
    return PanelInputs(
        benchmark=inputs.benchmark,
        panel=inputs.panel,
        species=inputs.species,
        cross_species=inputs.cross_species,
        annotation_database=inputs.annotation_database,
        ref_fa=linked["ref_fa"],
        ref_gff=linked["ref_gff"],
        tgt_fa=linked["tgt_fa"],
        liftoff_gff=linked.get("liftoff_gff"),
        miniprot_gff=linked.get("miniprot_gff"),
        transcripts_fa=linked.get("transcripts_fa"),
        proteins_fa=linked.get("proteins_fa"),
        truth_gff=linked.get("truth_gff"),
    )


def build_lifton_argv(
    source: SourceSpec,
    inputs: PanelInputs,
    output: Path,
    *,
    threads: int,
) -> list[str]:
    """Build the frozen subset/full protocol for either exact source tree."""

    argv = [
        str(source.lifton_executable),
        "-t",
        str(1 if inputs.panel == "subset" else threads),
        "-copies",
        "-ad",
        inputs.annotation_database,
        "-g",
        str(inputs.ref_gff),
    ]
    if inputs.liftoff_gff is not None:
        argv.extend(["-L", str(inputs.liftoff_gff)])
    if inputs.miniprot_gff is not None:
        argv.extend(["-M", str(inputs.miniprot_gff)])
    if inputs.transcripts_fa is not None:
        argv.extend(["-T", str(inputs.transcripts_fa)])
    if inputs.proteins_fa is not None:
        argv.extend(["-P", str(inputs.proteins_fa)])
    if inputs.panel == "full":
        argv.append("--locus-pipeline")
    argv.extend([
        "--no-miniprot-rescue",
        "-o",
        str(output),
        str(inputs.tgt_fa),
        str(inputs.ref_fa),
    ])
    return argv


def _validation_document(path: Path) -> dict[str, Any]:
    result = validate_gff3_file(str(path))
    return {
        "is_valid": result.is_valid,
        "stats": dict(sorted(result.stats.items())),
        "n_errors": len(result.errors),
        "n_warnings": len(result.warnings),
        "issues": [
            {
                "severity": str(issue.severity),
                "lineno": issue.lineno,
                "feature_id": issue.feature_id,
                "check": issue.check,
                "message": issue.message,
            }
            for issue in result.issues
        ],
    }


def _native_manifest_record(path: Path) -> dict[str, Any]:
    record: dict[str, Any] = {
        "path": str(path.resolve()),
        "present": path.is_file(),
    }
    if path.is_file():
        record.update({
            "size": path.stat().st_size,
            "sha256": sha256_file(path),
        })
    return record


def _write_release_arm_manifest(
    version_dir: Path,
    *,
    source: Mapping[str, Any],
    protocol: Mapping[str, Any],
    profile: Mapping[str, Any],
    output: Path,
    fingerprints: Mapping[str, Any],
    validation: Mapping[str, Any],
    native_manifests: Mapping[str, Path],
) -> tuple[Path, dict[str, Any]]:
    """Publish source-neutral evidence after one arm passes every local check."""

    manifest_path = version_dir / "release_run_manifest.json"
    native = {
        name: _native_manifest_record(path)
        for name, path in sorted(native_manifests.items())
    }
    payload = {
        "schema_version": SCHEMA_VERSION,
        "kind": "paired_release_arm",
        "run": {
            "status": "success",
            "finished_at": (
                dt.datetime.now(dt.timezone.utc)
                .isoformat()
                .replace("+00:00", "Z")
            ),
        },
        "source": dict(source),
        "protocol": dict(protocol),
        "profile": dict(profile),
        "artifacts": {
            "output_gff": {
                "path": str(output.resolve()),
                "size": output.stat().st_size,
                "byte_sha256": fingerprints["byte_sha256"],
                "semantic_sha256": fingerprints["semantic_sha256"],
            },
            "native_manifests": native,
        },
        "validation": dict(validation),
    }
    temporary = manifest_path.with_suffix(".json.tmp")
    with temporary.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, default=str)
        handle.write("\n")
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(temporary, manifest_path)
    return manifest_path, payload


def _run_one(
    source: SourceSpec,
    inputs: PanelInputs,
    version_dir: Path,
    *,
    threads: int,
) -> tuple[Path, dict[str, Any]]:
    version_dir.mkdir(parents=True, exist_ok=False)
    isolated = _isolated_inputs(inputs, version_dir)
    output = version_dir / f"{source.label}.gff3"
    argv = build_lifton_argv(source, isolated, output, threads=threads)
    profile = run_profiled(
        argv,
        label=source.label,
        log_dir=version_dir / "logs",
        env=source.environment(),
        cwd=version_dir,
        log=lambda message: print(message, flush=True),
    )
    if profile.exit_code != 0:
        raise RuntimeError(
            f"{source.label} exited {profile.exit_code}; see {profile.stderr_path}"
        )
    if not output.is_file() or output.stat().st_size == 0:
        raise RuntimeError(f"{source.label} did not publish a non-empty GFF3")
    validation = _validation_document(output)
    if not validation["is_valid"]:
        raise RuntimeError(
            f"{source.label} produced invalid GFF3 with {validation['n_errors']} errors"
        )
    document = {
        "source": dataclasses.asdict(source),
        "argv": argv,
        "profile": dataclasses.asdict(profile),
        "output_gff": str(output),
        "fingerprints": gff3_fingerprints(output),
        "validation": validation,
    }
    document["source"]["root"] = str(source.root)
    document["source"]["lifton_executable"] = str(source.lifton_executable)
    release_manifest, release_manifest_document = _write_release_arm_manifest(
        version_dir,
        source=document["source"],
        protocol={"kind": inputs.panel, "argv": argv},
        profile=document["profile"],
        output=output,
        fingerprints=document["fingerprints"],
        validation=validation,
        native_manifests={
            "lift": version_dir / "lifton_output" / "run_manifest.json",
        },
    )
    document["release_manifest"] = str(release_manifest)
    document["native_manifests"] = release_manifest_document["artifacts"][
        "native_manifests"
    ]
    return output, document


@contextlib.contextmanager
def _temporary_environment(values: Mapping[str, str]):
    previous = {name: os.environ.get(name) for name in values}
    os.environ.update({str(name): str(value) for name, value in values.items()})
    try:
        yield
    finally:
        for name, value in previous.items():
            if value is None:
                os.environ.pop(name, None)
            else:
                os.environ[name] = value


def _e2e_dataset(
    benchmark: str,
    dataset_registry: Path,
) -> tuple[run_benchmarks.Registry, run_benchmarks.Dataset]:
    registry = run_benchmarks.load_registry(Path(dataset_registry))
    matches = [item for item in registry.datasets if item.id == benchmark]
    if len(matches) != 1:
        raise ValueError(f"unknown end-to-end dataset: {benchmark}")
    return registry, matches[0]


def e2e_flags(
    mode: str,
    *,
    threads: int,
    dataset_registry: Path = run_benchmarks.DEFAULT_REGISTRY,
) -> list[str]:
    """Return one explicit fast-path combination on the production protocol."""

    registry = run_benchmarks.load_registry(Path(dataset_registry))
    try:
        enabled = E2E_MODE_FEATURES[mode]
    except KeyError as exc:
        raise ValueError(
            f"unknown end-to-end mode {mode!r}; choose from {', '.join(E2E_MODES)}"
        ) from exc
    accelerated = E2E_MODE_FEATURES["fast"]
    flags = [
        flag for flag in registry.lifton_flags
        if flag not in accelerated or flag in enabled
    ]
    try:
        index = flags.index("-t") + 1
        flags[index] = str(threads)
    except (ValueError, IndexError):
        flags.extend(["-t", str(threads)])
    return flags


def validate_e2e_biology(row: Mapping[str, Any]) -> dict[str, Any]:
    """Fail closed when the E2E harness did not produce meaningful biology."""

    biological = row.get("biological_summary")
    if not isinstance(biological, Mapping):
        raise RuntimeError("end-to-end result has no biological_summary")
    errors = run_benchmarks.validate_biological_result(
        dict(row),
        require_evaluation=True,
        require_identity=True,
    )
    if errors:
        raise RuntimeError(
            "end-to-end biological result is invalid: " + "; ".join(errors)
        )
    return dict(biological)


def _run_e2e_one(
    source: SourceSpec,
    inputs: PanelInputs,
    version_dir: Path,
    *,
    threads: int,
    mode: str,
    dataset_registry: Path,
) -> tuple[Path, dict[str, Any]]:
    version_dir.mkdir(parents=True, exist_ok=False)
    registry, dataset = _e2e_dataset(inputs.benchmark, dataset_registry)
    dataset = dataclasses.replace(dataset, target_gff=None)
    dataset_dir = version_dir / "data" / dataset.id
    source_paths = {
        dataset.reference_fa: inputs.ref_fa,
        dataset.target_fa: inputs.tgt_fa,
        dataset.reference_gff: inputs.ref_gff,
    }
    for url, source_path in source_paths.items():
        _isolated_link(
            source_path,
            dataset_dir / run_benchmarks._filename_for(url),
        )
    flags = e2e_flags(
        mode, threads=threads, dataset_registry=dataset_registry,
    )
    with _temporary_environment(source.environment()):
        row = run_benchmarks.run_dataset(
            dataset,
            data_dir=version_dir / "data",
            results_dir=version_dir / "artifacts",
            lifton_flags=flags,
            evaluation_flags=list(registry.evaluation_flags),
            do_download=False,
            do_lift=True,
            do_evaluate=True,
            force=True,
            log=lambda message: print(message, flush=True),
        )
    if row.get("error"):
        raise RuntimeError(f"{source.label} end-to-end error: {row['error']}")
    profile = row.get("lift_profile") or {}
    eval_profile = row.get("eval_profile") or {}
    if profile.get("exit_code") != 0 or eval_profile.get("exit_code") != 0:
        raise RuntimeError(f"{source.label} end-to-end profile is unsuccessful")
    biological_summary = validate_e2e_biology(row)
    output = Path(row["out_gff"])
    if not output.is_file() or output.stat().st_size == 0:
        raise RuntimeError(f"{source.label} did not publish a non-empty GFF3")
    validation = _validation_document(output)
    if not validation["is_valid"]:
        raise RuntimeError(
            f"{source.label} produced invalid GFF3 with {validation['n_errors']} errors"
        )
    source_document = dataclasses.asdict(source)
    source_document["root"] = str(source.root)
    source_document["lifton_executable"] = str(source.lifton_executable)
    document = {
        "source": source_document,
        "e2e_mode": mode,
        "lifton_flags": flags,
        "evaluation_flags": list(registry.evaluation_flags),
        "profile": profile,
        "evaluation_profile": eval_profile,
        "biological_summary": biological_summary,
        "score_summary": row.get("score_summary"),
        "evaluation_summary": row.get("evaluation_summary"),
        "output_gff": str(output),
        "fingerprints": gff3_fingerprints(output),
        "validation": validation,
    }
    release_manifest, release_manifest_document = _write_release_arm_manifest(
        version_dir,
        source=source_document,
        protocol={
            "kind": "e2e",
            "mode": mode,
            "lifton_flags": flags,
            "evaluation_flags": list(registry.evaluation_flags),
        },
        profile=profile,
        output=output,
        fingerprints=document["fingerprints"],
        validation=validation,
        native_manifests={
            "lift": output.parent / "lifton_output" / "run_manifest.json",
            "evaluation": Path(
                row.get("evaluation_manifest")
                or (
                    output.parent
                    / "evaluation"
                    / "lifton_output"
                    / "run_manifest.json"
                )
            ),
        },
    )
    document["release_manifest"] = str(release_manifest)
    document["native_manifests"] = release_manifest_document["artifacts"][
        "native_manifests"
    ]
    return output, document


def _score_pair(
    inputs: PanelInputs,
    outputs: Mapping[str, Path],
    documents: dict[str, dict[str, Any]],
    cell_dir: Path,
    *,
    threads: int,
) -> None:
    score_root = cell_dir / "score-input"
    score_input = score_root / f"reference-annotation{inputs.ref_gff.suffix}"
    isolated_ref_gff = _isolated_link(inputs.ref_gff, score_input)
    isolated_ref_fa = _isolated_fasta(
        inputs.ref_fa,
        score_root / f"reference-genome{''.join(inputs.ref_fa.suffixes)}",
    )
    isolated_tgt_fa = _isolated_fasta(
        inputs.tgt_fa,
        score_root / f"target-genome{''.join(inputs.tgt_fa.suffixes)}",
    )
    reference, reference_index = evaluator.build_reference(
        str(isolated_ref_gff),
        str(isolated_ref_fa),
        log=lambda message: print(message, flush=True),
    )
    reference_id_index = _declared_id_index(isolated_ref_gff)
    manifest = {
        "id": inputs.benchmark,
        "species": inputs.species,
        "cross_species": inputs.cross_species,
        "miniprot_target_space": "transcript",
        "protein_acc_to_mrna": {},
    }
    for label, output in outputs.items():
        summary = dict(evaluator.evaluate_tool(
            label,
            str(output),
            str(isolated_tgt_fa),
            reference,
            manifest,
            cell_dir / "evaluation",
            documents[label]["profile"],
            log=lambda message: print(message, flush=True),
            ref_index=reference_index,
            threads=threads,
        ))
        expected_tsv = (
            cell_dir / "evaluation" / f"{label}.transcripts.tsv"
        ).resolve()
        transcript_value = summary.get("transcripts_tsv")
        if not isinstance(transcript_value, str):
            raise RuntimeError(
                f"{label} neutral evaluator did not report transcripts_tsv"
            )
        transcript_path = Path(transcript_value).resolve()
        if transcript_path != expected_tsv:
            raise RuntimeError(
                f"{label} neutral evaluator transcripts_tsv is outside the "
                f"canonical cell path: {transcript_path}"
            )
        try:
            transcript_stat = transcript_path.stat()
        except OSError as exc:
            raise RuntimeError(
                f"{label} neutral evaluator TSV is unavailable: {exc}"
            ) from exc
        if not transcript_path.is_file() or transcript_stat.st_size <= 0:
            raise RuntimeError(
                f"{label} neutral evaluator TSV is missing or empty: "
                f"{transcript_path}"
            )
        summary["transcripts_tsv"] = str(transcript_path)
        summary["stable_id_preservation"] = stable_id_preservation(
            isolated_ref_gff,
            output,
            reference_index=reference_id_index,
        )
        documents[label]["summary"] = summary
        documents[label]["evaluation_artifacts"] = {
            "transcripts_tsv": {
                "path": str(transcript_path),
                "size": transcript_stat.st_size,
                "sha256": sha256_file(transcript_path),
            },
        }


def run_paired_cell(
    *,
    panel: str,
    benchmark: str,
    repetition: int,
    candidate: SourceSpec,
    reference: SourceSpec,
    cell_dir: Path,
    threads: int = 8,
    order: Sequence[str] | None = None,
    candidate_e2e_mode: str = "fast",
    reference_e2e_mode: str = "fast",
    benchmark_registry: Path = DEFAULT_BENCHMARK_REGISTRY,
    dataset_registry: Path = run_benchmarks.DEFAULT_REGISTRY,
) -> dict[str, Any]:
    """Run one isolated AB/BA repetition and write ``pair_result.json``."""

    if not 1 <= repetition <= 5:
        raise ValueError("repetition must be between 1 and 5")
    if threads < 1:
        raise ValueError("threads must be positive")
    if cell_dir.exists() and any(cell_dir.iterdir()):
        raise FileExistsError(f"paired cell directory is not empty: {cell_dir}")
    cell_dir.mkdir(parents=True, exist_ok=True)
    sources = {"candidate": candidate, "reference": reference}
    selected_order = list(order or (
        ("reference", "candidate") if repetition % 2 else ("candidate", "reference")
    ))
    if sorted(selected_order) != ["candidate", "reference"]:
        raise ValueError("paired order must contain candidate and reference exactly once")
    provenance = {label: verify_source(source) for label, source in sources.items()}
    benchmark_registry = Path(benchmark_registry).resolve()
    dataset_registry = Path(dataset_registry).resolve()
    inputs = resolve_panel_inputs(
        panel,
        benchmark,
        benchmark_registry=benchmark_registry,
        dataset_registry=dataset_registry,
    )
    inputs_document = input_fingerprints(inputs)
    outputs: dict[str, Path] = {}
    documents: dict[str, dict[str, Any]] = {}
    for label in selected_order:
        if panel == "e2e":
            mode = (
                candidate_e2e_mode if label == "candidate"
                else reference_e2e_mode
            )
            output, document = _run_e2e_one(
                sources[label],
                inputs,
                cell_dir / label,
                threads=threads,
                mode=mode,
                dataset_registry=dataset_registry,
            )
        else:
            output, document = _run_one(
                sources[label],
                inputs,
                cell_dir / label,
                threads=threads,
            )
        outputs[label] = output
        documents[label] = document
    _score_pair(inputs, outputs, documents, cell_dir, threads=threads)
    reference_profile = documents["reference"]["profile"]
    candidate_profile = documents["candidate"]["profile"]
    wall_ratio = (
        candidate_profile["wall_clock_seconds"]
        / reference_profile["wall_clock_seconds"]
    )
    rss_ratio = (
        candidate_profile["peak_rss_mb"]
        / reference_profile["peak_rss_mb"]
    )
    payload = {
        "schema_version": SCHEMA_VERSION,
        "panel": panel,
        "benchmark": benchmark,
        "repetition": repetition,
        "order": selected_order,
        "threads": threads,
        "modes": {
            "candidate": candidate_e2e_mode if panel == "e2e" else panel,
            "reference": reference_e2e_mode if panel == "e2e" else panel,
        },
        "registries": {
            "benchmark": str(benchmark_registry),
            "dataset": str(dataset_registry),
        },
        "provenance": provenance,
        "inputs": inputs_document,
        "versions": documents,
        "ratios": {
            "wall": wall_ratio,
            "peak_rss": rss_ratio,
        },
    }
    output_path = cell_dir / "pair_result.json"
    temporary = output_path.with_suffix(".json.tmp")
    temporary.write_text(json.dumps(payload, indent=2, default=str) + "\n")
    os.replace(temporary, output_path)
    return payload


def _percentile(values: Sequence[float], quantile: float) -> float | None:
    ordered = sorted(values)
    if not ordered:
        return None
    if len(ordered) == 1:
        return ordered[0]
    position = (len(ordered) - 1) * quantile
    lower = int(position)
    upper = min(lower + 1, len(ordered) - 1)
    fraction = position - lower
    return ordered[lower] * (1.0 - fraction) + ordered[upper] * fraction


def bootstrap_geomean_ratio(
    ratios: Sequence[float],
    *,
    replicates: int = DEFAULT_BOOTSTRAP_REPLICATES,
    seed: int = DEFAULT_BOOTSTRAP_SEED,
) -> dict[str, Any] | None:
    """Deterministic cell-level bootstrap CI for a geometric mean ratio."""

    values = [float(value) for value in ratios if value and value > 0]
    if not values:
        return None
    logs = [math.log(value) for value in values]
    if len(logs) == 1:
        estimate = values[0]
        return {
            "estimate": estimate,
            "low": estimate,
            "high": estimate,
            "replicates": 1,
            "seed": seed,
        }
    import numpy as np

    rng = np.random.default_rng(seed)
    array = np.asarray(logs, dtype=float)
    means = []
    remaining = max(1, int(replicates))
    while remaining:
        batch = min(256, remaining)
        indices = rng.integers(0, len(array), size=(batch, len(array)))
        means.extend(np.exp(array[indices].mean(axis=1)).tolist())
        remaining -= batch
    return {
        "estimate": math.exp(statistics.fmean(logs)),
        "low": _percentile(means, 0.025),
        "high": _percentile(means, 0.975),
        "replicates": len(means),
        "seed": seed,
    }


def _source_from_args(label: str, root: str, sha: str, executable: str) -> SourceSpec:
    return SourceSpec(
        label=label,
        root=Path(root).resolve(),
        sha=sha,
        lifton_executable=Path(executable).resolve(),
    )


def _positive_int(value: str) -> int:
    try:
        parsed = int(value)
    except ValueError as exc:
        raise argparse.ArgumentTypeError("must be a positive integer") from exc
    if parsed < 1:
        raise argparse.ArgumentTypeError("must be a positive integer")
    return parsed


def _repetition_count(value: str) -> int:
    parsed = _positive_int(value)
    if parsed > 5:
        raise argparse.ArgumentTypeError("must be between 1 and 5")
    return parsed


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    run = subparsers.add_parser("run-pair", help="run one paired benchmark cell")
    run.add_argument("--panel", choices=("subset", "full", "e2e"), required=True)
    run.add_argument("--benchmark", required=True)
    run.add_argument("--repetition", type=_repetition_count, required=True)
    run.add_argument("--candidate-root", required=True)
    run.add_argument("--candidate-sha", required=True)
    run.add_argument("--reference-root", required=True)
    run.add_argument("--reference-sha", default=DEFAULT_REFERENCE_SHA)
    run.add_argument(
        "--lifton-executable",
        default=str(Path(sys.executable).with_name("lifton")),
    )
    run.add_argument("--cell-dir", required=True)
    run.add_argument("--threads", type=_positive_int, default=8)
    run.add_argument(
        "--benchmark-registry",
        default=str(DEFAULT_BENCHMARK_REGISTRY),
    )
    run.add_argument(
        "--dataset-registry",
        default=str(run_benchmarks.DEFAULT_REGISTRY),
    )
    run.add_argument(
        "--candidate-e2e-mode",
        choices=E2E_MODES,
        default="fast",
    )
    run.add_argument(
        "--reference-e2e-mode",
        choices=E2E_MODES,
        default="fast",
    )
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if args.command == "run-pair":
        executable = str(Path(args.lifton_executable).resolve())
        candidate = _source_from_args(
            "candidate", args.candidate_root, args.candidate_sha, executable,
        )
        reference = _source_from_args(
            "reference", args.reference_root, args.reference_sha, executable,
        )
        run_paired_cell(
            panel=args.panel,
            benchmark=args.benchmark,
            repetition=args.repetition,
            candidate=candidate,
            reference=reference,
            cell_dir=Path(args.cell_dir).resolve(),
            threads=args.threads,
            candidate_e2e_mode=args.candidate_e2e_mode,
            reference_e2e_mode=args.reference_e2e_mode,
            benchmark_registry=Path(args.benchmark_registry).resolve(),
            dataset_registry=Path(args.dataset_registry).resolve(),
        )
        return 0
    raise AssertionError(args.command)


if __name__ == "__main__":
    raise SystemExit(main())
