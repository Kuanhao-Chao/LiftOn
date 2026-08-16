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
import functools
import hashlib
import importlib.machinery
import json
import math
import os
import re
import resource
import shutil
import statistics
import subprocess
import sys
import time
import urllib.parse
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping, Sequence

from lifton.gff3_validator import (
    validate_gff3_file,
    validate_gff3_structure,
    validate_gff3_target_bounds,
)

from benchmarks import run_benchmarks

from . import devel_refresh
from . import evaluator
from . import fourway_compare
from . import target_truth
from .profiling import run_profiled


HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parent.parent
DEFAULT_BENCHMARK_REGISTRY = HERE / "benchmarks.json"
DEFAULT_REFERENCE_SHA = "3739dfc8f73396fccab7d7be0f95e008179cea5d"
SCHEMA_VERSION = 4
NEUTRAL_EVALUATION_FORMAT = "neutral_evaluator_v1"
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
OPTIONAL_EXECUTION_FLAGS = frozenset({
    "--stream",
    "--inmemory-liftoff",
    "--locus-pipeline",
    "--native",
})
PROBED_CLI_OPTIONS = tuple(sorted({
    *OPTIONAL_EXECUTION_FLAGS,
    "--no-miniprot-rescue",
}))
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


@functools.lru_cache(maxsize=None)
def _source_cli_options(
    root: str,
    sha: str,
    executable: str,
) -> frozenset[str]:
    source = SourceSpec(
        label="probe",
        root=Path(root),
        sha=sha,
        lifton_executable=Path(executable),
    )
    try:
        result = subprocess.run(
            [str(source.lifton_executable), "--help"],
            cwd=source.root,
            env=source.environment(),
            text=True,
            capture_output=True,
            check=False,
            timeout=30,
        )
    except (OSError, subprocess.TimeoutExpired) as exc:
        raise RuntimeError(
            f"could not probe LiftOn CLI options for {source.root}: {exc}"
        ) from exc
    if result.returncode != 0:
        raise RuntimeError(
            f"could not probe LiftOn CLI options for {source.root}: "
            f"--help exited {result.returncode}"
        )
    help_text = f"{result.stdout}\n{result.stderr}"
    return frozenset(
        option for option in PROBED_CLI_OPTIONS
        if option in help_text
    )


def source_cli_options(source: SourceSpec) -> frozenset[str]:
    return _source_cli_options(
        str(source.root.resolve()),
        source.sha,
        str(source.lifton_executable.resolve()),
    )


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
    truth_source_gff: Path | None = None
    ortholog_map: Path | None = None
    truth_id_policy: str = "ortholog-map"

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

    def evaluation_paths(self) -> dict[str, Path]:
        """Inputs consumed only after both LiftOn arms have completed."""

        optional = {
            "truth_gff": self.truth_gff,
            "truth_source_gff": self.truth_source_gff,
            "ortholog_map": self.ortholog_map,
        }
        return {
            name: path for name, path in optional.items() if path is not None
        }


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
        "cli_options": sorted(source_cli_options(spec)),
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


def _evaluation_input_paths(
    bench: Mapping[str, Any],
    registry_path: Path,
    panel: str | None = None,
) -> dict[str, Any]:
    by_panel = bench.get("target_truth_by_panel")
    if by_panel is not None:
        if not isinstance(by_panel, Mapping):
            raise ValueError("target_truth_by_panel must be an object")
        unknown = set(by_panel) - {"subset", "full", "e2e"}
        if unknown:
            raise ValueError(
                f"target_truth_by_panel has unsupported panel "
                f"{sorted(unknown)[0]!r}"
            )
        if panel is None:
            raise ValueError(
                "panel is required when target_truth_by_panel is declared"
            )
        if panel not in by_panel:
            raise ValueError(
                f"target_truth_by_panel has no exact {panel!r} selection"
            )
        truth = by_panel[panel]
        if not isinstance(truth, Mapping):
            raise ValueError(
                f"target_truth_by_panel.{panel} must be an object"
            )
        truth_document = truth
        legacy_values = {}
    else:
        truth = bench.get("target_truth")
        truth_document = truth if isinstance(truth, Mapping) else {}
        legacy_values = bench
    values = {
        "truth_gff": (
            truth_document.get("gff")
            or truth_document.get("truth_gff")
            or legacy_values.get("truth_gff")
        ),
        "ortholog_map": (
            truth_document.get("ortholog_map")
            or legacy_values.get("ortholog_map")
        ),
    }
    resolved = {}
    for name, value in values.items():
        if not value:
            continue
        path = Path(str(value)).expanduser()
        if not path.is_absolute():
            path = registry_path.parent / path
        resolved[name] = path.resolve()
    if "ortholog_map" in resolved and "truth_gff" not in resolved:
        raise ValueError("ortholog_map requires a target truth GFF3")
    if "truth_gff" in resolved:
        policy = str(
            truth_document.get("id_policy", "ortholog-map")
        ).strip().lower()
        if policy not in {"ortholog-map", "exact-id"}:
            raise ValueError(
                "target_truth.id_policy must be ortholog-map or exact-id"
            )
        if policy == "ortholog-map" and "ortholog_map" not in resolved:
            raise ValueError(
                "independent target truth requires a non-empty ortholog_map"
            )
        if policy == "exact-id" and "ortholog_map" in resolved:
            raise ValueError(
                "exact-id target truth cannot also declare ortholog_map"
            )
        resolved["truth_id_policy"] = policy
    return resolved


def _e2e_dataset_asset_path(value: str, benchmark: str) -> Path:
    """Resolve a local canonical overlay or the legacy download-cache path."""

    parsed = urllib.parse.urlsplit(str(value))
    if parsed.scheme.lower() == "file":
        if parsed.netloc not in {"", "localhost"}:
            raise ValueError(
                f"{benchmark}: file dataset assets must use a local authority"
            )
        if parsed.query or parsed.fragment:
            raise ValueError(
                f"{benchmark}: file dataset assets cannot contain a query or "
                "fragment"
            )
        path = Path(urllib.parse.unquote(parsed.path)).expanduser()
        if not path.is_absolute():
            raise ValueError(
                f"{benchmark}: file dataset assets must use an absolute path"
            )
        return path.resolve()
    return (
        run_benchmarks.DEFAULT_DATA_DIR
        / benchmark
        / run_benchmarks._filename_for(str(value))
    )


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
        truth_gff = (
            _e2e_dataset_asset_path(dataset.truth_gff, benchmark)
            if dataset.truth_gff
            else None
        )
        ortholog_map = (
            _e2e_dataset_asset_path(dataset.ortholog_map, benchmark)
            if dataset.ortholog_map
            else None
        )
        truth_source_gff = (
            _e2e_dataset_asset_path(dataset.truth_source_gff, benchmark)
            if dataset.truth_source_gff
            else None
        )
        truth_id_policy = str(dataset.truth_id_policy).strip().lower()
        if truth_id_policy not in {"ortholog-map", "exact-id"}:
            raise ValueError(
                f"{benchmark}: truth_id_policy must be ortholog-map or exact-id"
            )
        if ortholog_map is not None and truth_gff is None:
            raise ValueError(
                f"{benchmark}: ortholog_map requires an independent truth_gff"
            )
        if truth_gff is not None:
            if truth_id_policy == "ortholog-map" and ortholog_map is None:
                raise ValueError(
                    f"{benchmark}: independent E2E target truth requires a "
                    "non-empty ortholog_map"
                )
            if truth_id_policy == "exact-id" and ortholog_map is not None:
                raise ValueError(
                    f"{benchmark}: exact-id E2E truth cannot also declare "
                    "ortholog_map"
                )
        result = PanelInputs(
            benchmark=benchmark,
            panel=panel,
            species=dataset.species,
            cross_species=bool(dataset.cross_species),
            annotation_database=str(dataset.annotation_database),
            ref_fa=_e2e_dataset_asset_path(dataset.reference_fa, benchmark),
            ref_gff=_e2e_dataset_asset_path(dataset.reference_gff, benchmark),
            tgt_fa=_e2e_dataset_asset_path(dataset.target_fa, benchmark),
            truth_gff=truth_gff,
            truth_source_gff=truth_source_gff,
            ortholog_map=ortholog_map,
            truth_id_policy=truth_id_policy,
        )
        _require_nonempty(result.required_paths())
        _require_nonempty(result.evaluation_paths())
        return result

    benchmark_registry = Path(benchmark_registry).resolve()
    bench = _benchmark_entry(benchmark_registry, benchmark)
    common = {
        "benchmark": benchmark,
        "panel": panel,
        "species": str(bench["species"]),
        "cross_species": bool(bench["cross_species"]),
        "annotation_database": str(bench.get("annotation_database", "RefSeq")),
        **_evaluation_input_paths(bench, benchmark_registry, panel),
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
        full_input_mode = str(
            bench.get("full_input_mode", "cached")
        ).strip().lower()
        if full_input_mode == "cached":
            liftoff_gff, miniprot_gff, transcripts_fa, proteins_fa = (
                devel_refresh._cached_inputs(benchmark)
            )
        elif full_input_mode == "raw":
            liftoff_gff = miniprot_gff = None
            transcripts_fa = proteins_fa = None
        else:
            raise ValueError(
                f"{benchmark}: full_input_mode must be cached or raw"
            )
        result = PanelInputs(
            **common,
            ref_fa=Path(paths["ref_fa"]),
            ref_gff=Path(paths["ref_gff"]),
            tgt_fa=Path(paths["tgt_fa"]),
            liftoff_gff=(
                Path(liftoff_gff) if liftoff_gff is not None else None
            ),
            miniprot_gff=(
                Path(miniprot_gff) if miniprot_gff is not None else None
            ),
            transcripts_fa=(
                Path(transcripts_fa) if transcripts_fa is not None else None
            ),
            proteins_fa=(
                Path(proteins_fa) if proteins_fa is not None else None
            ),
        )
    else:
        raise ValueError(
            f"paired panel must be subset, full, or e2e, not {panel!r}"
        )
    _require_nonempty(result.required_paths())
    _require_nonempty(result.evaluation_paths())
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


def evaluation_input_fingerprints(inputs: PanelInputs) -> dict[str, Any]:
    """Fingerprint truth inputs separately from LiftOn execution inputs."""

    return {
        name: {
            "path": str(path.resolve()),
            "size": path.stat().st_size,
            "sha256": sha256_file(path),
        }
        for name, path in sorted(inputs.evaluation_paths().items())
    }


def _isolated_link(source: Path, destination: Path, *, reuse: bool = False) -> Path:
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists() or destination.is_symlink():
        if (
            reuse
            and destination.is_symlink()
            and destination.resolve() == source.resolve()
        ):
            # A re-score runs inside a cell that is already populated. Reusing
            # a link that resolves to the very same input keeps the cell's
            # cached annotation databases valid; anything else still fails.
            return destination
        raise FileExistsError(f"refusing to replace paired input link: {destination}")
    destination.symlink_to(source.resolve())
    return destination


def _isolated_fasta(source: Path, destination: Path, *, reuse: bool = False) -> Path:
    """Link a read-only FASTA and copy mutable indexes into the cell."""

    linked = _isolated_link(source, destination, reuse=reuse)
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
    if (
        inputs.panel == "full"
        and "--locus-pipeline" in source_cli_options(source)
    ):
        argv.append("--locus-pipeline")
    if "--no-miniprot-rescue" in source_cli_options(source):
        argv.append("--no-miniprot-rescue")
    argv.extend([
        "-o",
        str(output),
        str(inputs.tgt_fa),
        str(inputs.ref_fa),
    ])
    return argv


def _issue_document(issue: Any, validator: str) -> dict[str, Any]:
    return {
        "validator": validator,
        "severity": str(issue.severity),
        "lineno": issue.lineno,
        "feature_id": issue.feature_id,
        "check": issue.check,
        "message": issue.message,
    }


def _validation_pass_document(result: Any) -> dict[str, Any]:
    return {
        "is_valid": result.is_valid,
        "n_errors": result.severity_totals.get(
            "ERROR", len(result.errors),
        ),
        "n_warnings": result.severity_totals.get(
            "WARNING", len(result.warnings),
        ),
        "n_issues_reported": len(result.issues),
        "issue_totals": dict(sorted(result.issue_totals.items())),
    }


def _validation_document(
    path: Path,
    target_fasta: Path | None = None,
) -> dict[str, Any]:
    full = validate_gff3_file(str(path))
    structural = validate_gff3_structure(str(path))
    passes = {
        "full_semantic": _validation_pass_document(full),
        "streaming_structure": _validation_pass_document(structural),
    }
    results = (("full_semantic", full), ("streaming_structure", structural))
    if target_fasta is not None:
        target_bounds = validate_gff3_target_bounds(
            str(path), str(target_fasta), strict_sequence_regions=True,
        )
        passes["target_fasta_bounds"] = {
            **_validation_pass_document(target_bounds),
            "target_fasta": str(target_fasta.resolve()),
        }
        results = (*results, ("target_fasta_bounds", target_bounds))
    issues = [
        _issue_document(issue, validator)
        for validator, result in results
        for issue in result.issues
    ]
    return {
        "is_valid": all(item["is_valid"] for item in passes.values()),
        "stats": dict(sorted(full.stats.items())),
        "n_errors": sum(item["n_errors"] for item in passes.values()),
        "n_warnings": sum(item["n_warnings"] for item in passes.values()),
        "passes": passes,
        "issues": issues,
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
    validation = _validation_document(output, isolated.tgt_fa)
    if not validation["is_valid"]:
        print(
            f"[bench] {source.label} validation reported "
            f"{validation['n_errors']} errors",
            flush=True,
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
    source: SourceSpec | None = None,
    dataset_registry: Path = run_benchmarks.DEFAULT_REGISTRY,
) -> list[str]:
    """Return the requested fast-path combination supported by one release."""

    registry = run_benchmarks.load_registry(Path(dataset_registry))
    try:
        enabled = E2E_MODE_FEATURES[mode]
    except KeyError as exc:
        raise ValueError(
            f"unknown end-to-end mode {mode!r}; choose from {', '.join(E2E_MODES)}"
        ) from exc
    requested = frozenset({"--locus-pipeline", *enabled})
    supported = source_cli_options(source) if source is not None else requested
    flags = [
        flag for flag in registry.lifton_flags
        if (
            flag not in OPTIONAL_EXECUTION_FLAGS
            or (flag in requested and flag in supported)
        )
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
    dataset = dataclasses.replace(
        dataset,
        target_gff=None,
        truth_gff=None,
        ortholog_map=None,
    )
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
        mode,
        threads=threads,
        source=source,
        dataset_registry=dataset_registry,
    )
    source_environment = source.environment()
    source_environment[run_benchmarks.PROCESS_GROUP_PROFILE_ENV] = "1"
    with _temporary_environment(source_environment):
        row = run_benchmarks.run_dataset(
            dataset,
            data_dir=version_dir / "data",
            results_dir=version_dir / "artifacts",
            lifton_flags=flags,
            evaluation_flags=list(registry.evaluation_flags),
            do_download=False,
            do_lift=True,
            do_evaluate=False,
            force=True,
            log=lambda message: print(message, flush=True),
        )
    if row.get("error"):
        raise RuntimeError(f"{source.label} end-to-end error: {row['error']}")
    profile = row.get("lift_profile") or {}
    if profile.get("exit_code") != 0:
        raise RuntimeError(f"{source.label} end-to-end lift is unsuccessful")
    output = Path(row["out_gff"])
    if not output.is_file() or output.stat().st_size == 0:
        raise RuntimeError(f"{source.label} did not publish a non-empty GFF3")
    validation = _validation_document(output, inputs.tgt_fa)
    if not validation["is_valid"]:
        print(
            f"[bench] {source.label} validation reported "
            f"{validation['n_errors']} errors",
            flush=True,
        )
    source_document = dataclasses.asdict(source)
    source_document["root"] = str(source.root)
    source_document["lifton_executable"] = str(source.lifton_executable)
    document = {
        "source": source_document,
        "e2e_mode": mode,
        "lifton_flags": flags,
        "evaluation_flags": [],
        "profile": profile,
        "evaluation_profile": None,
        "evaluation_method": NEUTRAL_EVALUATION_FORMAT,
        "biological_summary": None,
        "score_summary": row.get("score_summary"),
        "evaluation_summary": None,
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
            "evaluation_flags": [],
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


def _neutral_evaluation_profile(started: float) -> dict[str, Any]:
    usage = resource.getrusage(resource.RUSAGE_SELF)
    return {
        "exit_code": 0,
        "wall_clock_seconds": max(time.perf_counter() - started, 1e-9),
        "peak_rss_mb": max(float(usage.ru_maxrss) / 1024.0, 1e-9),
        "method": NEUTRAL_EVALUATION_FORMAT,
    }


def _neutral_e2e_summary(summary: Mapping[str, Any]) -> dict[str, Any]:
    reference_total = int(summary["n_reference_total"])
    reference_coding = int(summary["n_reference_coding"])
    protein_identity = summary.get("protein_identity") or {}
    return {
        "format": NEUTRAL_EVALUATION_FORMAT,
        "records": reference_total,
        "coding_records": reference_coding,
        "noncoding_records": reference_total - reference_coding,
        "malformed_records": 0,
        "avg_identity": protein_identity.get("mean"),
        "score_file": summary.get("transcripts_tsv"),
    }


def _neutral_e2e_biology(
    summary: Mapping[str, Any], score_summary: Mapping[str, Any] | None,
) -> dict[str, Any]:
    reference_total = int(summary["n_reference_total"])
    recovered = int(summary["n_recovered_any"])
    emitted = (
        int(score_summary["records"])
        if isinstance(score_summary, Mapping)
        and isinstance(score_summary.get("records"), int)
        and score_summary["records"] > 0
        else int(summary["n_tool_features"])
    )
    return {
        "schema_version": 1,
        "reference_features": reference_total,
        "mapped_features": recovered,
        "lost_features": reference_total - recovered,
        "extra_copy_features": int(summary["n_extra_copies"]),
        "feature_completeness": recovered / reference_total,
        "emitted_transcript_records": emitted,
        "mapped_transcripts_reported": recovered,
        "evaluated_transcript_records": reference_total,
        "evaluated_coding_records": int(summary["n_reference_coding"]),
        "mean_protein_identity": (
            (summary.get("protein_identity") or {}).get("mean")
        ),
    }


def _score_pair(
    inputs: PanelInputs,
    outputs: Mapping[str, Path],
    documents: dict[str, dict[str, Any]],
    cell_dir: Path,
    *,
    threads: int,
    sequence_scoring: bool = False,
    reuse_isolated_inputs: bool = False,
) -> None:
    score_root = cell_dir / "score-input"
    score_input = score_root / f"reference-annotation{inputs.ref_gff.suffix}"
    isolated_ref_gff = _isolated_link(
        inputs.ref_gff, score_input, reuse=reuse_isolated_inputs,
    )
    isolated_ref_fa = _isolated_fasta(
        inputs.ref_fa,
        score_root / f"reference-genome{''.join(inputs.ref_fa.suffixes)}",
        reuse=reuse_isolated_inputs,
    )
    isolated_tgt_fa = _isolated_fasta(
        inputs.tgt_fa,
        score_root / f"target-genome{''.join(inputs.tgt_fa.suffixes)}",
        reuse=reuse_isolated_inputs,
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
        neutral_started = time.perf_counter()
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
        if inputs.panel == "e2e":
            documents[label]["evaluation_profile"] = (
                _neutral_evaluation_profile(neutral_started)
            )
            documents[label]["evaluation_method"] = NEUTRAL_EVALUATION_FORMAT
            documents[label]["evaluation_summary"] = _neutral_e2e_summary(
                summary
            )
            documents[label]["biological_summary"] = _neutral_e2e_biology(
                summary, documents[label].get("score_summary"),
            )
        documents[label]["summary"] = summary
        documents[label]["evaluation_artifacts"] = {
            "transcripts_tsv": {
                "path": str(transcript_path),
                "size": transcript_stat.st_size,
                "sha256": sha256_file(transcript_path),
            },
        }

    if inputs.truth_gff is not None:
        truth_path = _isolated_link(
            inputs.truth_gff,
            score_root / f"target-truth{''.join(inputs.truth_gff.suffixes)}",
            reuse=reuse_isolated_inputs,
        )
        mapping_path = (
            _isolated_link(
                inputs.ortholog_map,
                score_root / (
                    f"ortholog-map{''.join(inputs.ortholog_map.suffixes)}"
                ),
                reuse=reuse_isolated_inputs,
            )
            if inputs.ortholog_map is not None
            else None
        )
        truth_source_path = (
            _isolated_link(
                inputs.truth_source_gff,
                score_root / (
                    "truth-source"
                    f"{''.join(inputs.truth_source_gff.suffixes)}"
                ),
                reuse=reuse_isolated_inputs,
            )
            if inputs.truth_source_gff is not None
            else isolated_ref_gff
        )
        truth_validation = validate_gff3_target_bounds(
            str(truth_path), str(isolated_tgt_fa),
            strict_sequence_regions=True,
        )
        if not truth_validation.is_valid:
            checks = ", ".join(
                sorted({issue.check for issue in truth_validation.errors})
            )
            raise RuntimeError(
                f"target truth annotation is outside the target FASTA: {checks}"
            )
        truth_evidence = {
            "id_policy": inputs.truth_id_policy,
            "mapping_required": inputs.truth_id_policy == "ortholog-map",
            "gff": {
                "path": str(truth_path.resolve()),
                "size": truth_path.stat().st_size,
                "sha256": sha256_file(truth_path),
            },
            "ortholog_map": (
                {
                    "path": str(mapping_path.resolve()),
                    "size": mapping_path.stat().st_size,
                    "sha256": sha256_file(mapping_path),
                }
                if mapping_path is not None
                else None
            ),
            "source_model_view": {
                "path": str(truth_source_path.resolve()),
                "size": truth_source_path.stat().st_size,
                "sha256": sha256_file(truth_source_path),
            },
            "target_bounds_validation": {
                "is_valid": True,
                "n_errors": truth_validation.severity_totals.get("ERROR", 0),
                "n_warnings": truth_validation.severity_totals.get("WARNING", 0),
                "n_issues_reported": len(truth_validation.issues),
                "issue_totals": dict(sorted(truth_validation.issue_totals.items())),
                "stats": dict(sorted(truth_validation.stats.items())),
                "issues": [
                    _issue_document(issue, "target_truth_bounds")
                    for issue in truth_validation.issues
                ],
            },
        }
        for label, output in outputs.items():
            metrics = target_truth.score_target_truth(
                output,
                truth_path,
                ortholog_map=mapping_path,
                source_gff=truth_source_path,
                target_fasta=(isolated_tgt_fa if sequence_scoring else None),
                id_policy=inputs.truth_id_policy,
            )
            metrics_path = target_truth.write_target_truth_metrics(
                cell_dir / "evaluation" / f"{label}.target_truth.json",
                metrics,
            )
            documents[label]["summary"]["target_truth"] = metrics
            documents[label]["target_truth_evidence"] = truth_evidence
            documents[label]["evaluation_artifacts"]["target_truth"] = {
                "path": str(metrics_path.resolve()),
                "size": metrics_path.stat().st_size,
                "sha256": sha256_file(metrics_path),
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

    if not 1 <= repetition <= 10:
        raise ValueError("repetition must be between 1 and 10")
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
    evaluation_inputs_document = evaluation_input_fingerprints(inputs)
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
    return _score_and_publish_pair(
        inputs,
        outputs,
        documents,
        cell_dir,
        panel=panel,
        benchmark=benchmark,
        repetition=repetition,
        selected_order=selected_order,
        threads=threads,
        candidate_e2e_mode=candidate_e2e_mode,
        reference_e2e_mode=reference_e2e_mode,
        benchmark_registry=benchmark_registry,
        dataset_registry=dataset_registry,
        provenance=provenance,
        inputs_document=inputs_document,
        evaluation_inputs_document=evaluation_inputs_document,
    )


def _score_and_publish_pair(
    inputs: PanelInputs,
    outputs: Mapping[str, Path],
    documents: dict[str, dict[str, Any]],
    cell_dir: Path,
    *,
    panel: str,
    benchmark: str,
    repetition: int,
    selected_order: list[str],
    threads: int,
    candidate_e2e_mode: str,
    reference_e2e_mode: str,
    benchmark_registry: Path,
    dataset_registry: Path,
    provenance: dict[str, Any],
    inputs_document: dict[str, Any],
    evaluation_inputs_document: dict[str, Any],
    reuse_isolated_inputs: bool = False,
) -> dict[str, Any]:
    """Score both arms and publish ``pair_result.json``.

    Shared by ``run-pair`` (which has just produced the arm outputs) and
    ``score-pair`` (which reuses arm outputs that already exist), so both
    entry points build byte-identical payloads from the same code.
    """

    _score_pair(
        inputs,
        outputs,
        documents,
        cell_dir,
        threads=threads,
        sequence_scoring=repetition == 1,
        reuse_isolated_inputs=reuse_isolated_inputs,
    )
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
        "evaluation_inputs": evaluation_inputs_document,
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


def rebuild_arm_document(cell_dir: Path, label: str) -> tuple[Path, dict[str, Any]]:
    """Rebuild one arm's lift-side document from artifacts already on disk.

    Every field is taken from the arm's own ``release_run_manifest.json`` or
    re-derived deterministically from the published GFF3 and score file, so
    the result is what ``_run_e2e_one`` returned for that execution. The
    published output is re-hashed and compared with the manifest, which fails
    closed if the artifact has drifted since it was written.
    """

    arm_dir = Path(cell_dir) / label
    manifest_path = arm_dir / "release_run_manifest.json"
    if not manifest_path.is_file():
        raise RuntimeError(f"{label} arm has no release_run_manifest.json")
    manifest = json.loads(manifest_path.read_text())
    if manifest.get("kind") != "paired_release_arm":
        raise RuntimeError(f"{label} arm manifest has an unexpected kind")
    if (manifest.get("run") or {}).get("status") != "success":
        raise RuntimeError(f"{label} arm did not complete successfully")

    artifacts = manifest["artifacts"]["output_gff"]
    output = Path(artifacts["path"])
    if not output.is_file() or output.stat().st_size == 0:
        raise RuntimeError(f"{label} arm output GFF3 is missing or empty")
    fingerprints = gff3_fingerprints(output)
    for field in ("byte_sha256", "semantic_sha256"):
        if fingerprints[field] != artifacts[field]:
            raise RuntimeError(
                f"{label} arm output GFF3 {field} no longer matches its "
                f"manifest; the artifact has drifted since it was published"
            )

    score_file = output.parent / "lifton_output" / "score.txt"
    if not score_file.is_file():
        raise RuntimeError(f"{label} arm has no lifton_output/score.txt")
    score_summary = dict(
        run_benchmarks.parse_score_txt(score_file).__dict__
    )

    protocol = manifest["protocol"]
    # Key order matches _run_e2e_one exactly; pair_result.json is emitted with
    # json.dumps, which preserves insertion order, so this is load-bearing.
    document: dict[str, Any] = {
        "source": manifest["source"],
        "e2e_mode": protocol["mode"],
        "lifton_flags": protocol["lifton_flags"],
        "evaluation_flags": protocol["evaluation_flags"],
        "profile": manifest["profile"],
        "evaluation_profile": None,
        "evaluation_method": NEUTRAL_EVALUATION_FORMAT,
        "biological_summary": None,
        "score_summary": score_summary,
        "evaluation_summary": None,
        "output_gff": str(output),
        "fingerprints": fingerprints,
        "validation": manifest["validation"],
    }
    document["release_manifest"] = str(manifest_path)
    document["native_manifests"] = manifest["artifacts"]["native_manifests"]
    return output, document


def score_existing_pair(
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
    """Re-score a paired cell whose two arms have already been executed.

    This repeats only the evaluation half of ``run_paired_cell``: both arms'
    GFF3 outputs are reused as published, and every metric is recomputed by
    the same ``_score_pair`` the original run used. Timings and resource
    measurements therefore remain those of the original execution.
    """

    cell_dir = Path(cell_dir)
    if not cell_dir.is_dir():
        raise FileNotFoundError(f"paired cell directory does not exist: {cell_dir}")
    if not 1 <= repetition <= 10:
        raise ValueError("repetition must be between 1 and 10")
    if threads < 1:
        raise ValueError("threads must be positive")
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
    evaluation_inputs_document = evaluation_input_fingerprints(inputs)
    outputs: dict[str, Path] = {}
    documents: dict[str, dict[str, Any]] = {}
    for label in selected_order:
        output, document = rebuild_arm_document(cell_dir, label)
        recorded_sha = document["source"].get("sha")
        expected_sha = sources[label].sha
        if recorded_sha != expected_sha:
            raise RuntimeError(
                f"{label} arm was produced by {recorded_sha}, not {expected_sha}"
            )
        outputs[label] = output
        documents[label] = document
    return _score_and_publish_pair(
        inputs,
        outputs,
        documents,
        cell_dir,
        panel=panel,
        benchmark=benchmark,
        repetition=repetition,
        selected_order=selected_order,
        threads=threads,
        candidate_e2e_mode=candidate_e2e_mode,
        reference_e2e_mode=reference_e2e_mode,
        benchmark_registry=benchmark_registry,
        dataset_registry=dataset_registry,
        provenance=provenance,
        inputs_document=inputs_document,
        evaluation_inputs_document=evaluation_inputs_document,
        reuse_isolated_inputs=True,
    )


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
    if parsed > 10:
        raise argparse.ArgumentTypeError("must be between 1 and 10")
    return parsed


def _add_paired_cell_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--panel", choices=("subset", "full", "e2e"), required=True)
    parser.add_argument("--benchmark", required=True)
    parser.add_argument("--repetition", type=_repetition_count, required=True)
    parser.add_argument("--candidate-root", required=True)
    parser.add_argument("--candidate-sha", required=True)
    parser.add_argument("--reference-root", required=True)
    parser.add_argument("--reference-sha", default=DEFAULT_REFERENCE_SHA)
    parser.add_argument(
        "--lifton-executable",
        default=str(Path(sys.executable).with_name("lifton")),
    )
    parser.add_argument("--cell-dir", required=True)
    parser.add_argument("--threads", type=_positive_int, default=8)
    parser.add_argument(
        "--benchmark-registry",
        default=str(DEFAULT_BENCHMARK_REGISTRY),
    )
    parser.add_argument(
        "--dataset-registry",
        default=str(run_benchmarks.DEFAULT_REGISTRY),
    )
    parser.add_argument(
        "--candidate-e2e-mode",
        choices=E2E_MODES,
        default="fast",
    )
    parser.add_argument(
        "--reference-e2e-mode",
        choices=E2E_MODES,
        default="fast",
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    _add_paired_cell_arguments(
        subparsers.add_parser("run-pair", help="run one paired benchmark cell")
    )
    _add_paired_cell_arguments(
        subparsers.add_parser(
            "score-pair",
            help=(
                "re-score a paired cell whose two arms have already been "
                "executed, reusing their published outputs"
            ),
        )
    )
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if args.command not in {"run-pair", "score-pair"}:
        raise AssertionError(args.command)
    executable = str(Path(args.lifton_executable).resolve())
    candidate = _source_from_args(
        "candidate", args.candidate_root, args.candidate_sha, executable,
    )
    reference = _source_from_args(
        "reference", args.reference_root, args.reference_sha, executable,
    )
    entry = run_paired_cell if args.command == "run-pair" else score_existing_pair
    entry(
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


if __name__ == "__main__":
    raise SystemExit(main())
