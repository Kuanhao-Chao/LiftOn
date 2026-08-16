#!/usr/bin/env python3
"""Guarded, resumable tmux controller for LiftOn benchmark builds.

The controller deliberately keeps all scheduler state below
``benchmarks/compare/_runs/<run-id>`` and never writes the canonical benchmark
baseline.  A cell is successful only after its command, result schema, LiftOn
run manifest, and final GFF3 have all passed independent checks.

Typical use::

    python -m benchmarks.compare.build_controller start --stage subset-canary
    python -m benchmarks.compare.build_controller status <run-id>
    python -m benchmarks.compare.build_controller retry <run-id>
    python -m benchmarks.compare.build_controller reconcile <run-id>

``start --dry-run`` builds the immutable plan without contacting tmux.  This is
also useful for reviewing exact commands and provenance before a long run.
"""
from __future__ import annotations

import argparse
import contextlib
import datetime as dt
import hashlib
import importlib.metadata as importlib_metadata
import json
import math
import os
import re
import signal
import shlex
import shutil
import subprocess
import sys
import threading
import time
from collections import Counter
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Callable, Mapping, Sequence

try:
    import fcntl
except ImportError:  # pragma: no cover - tmux controller targets POSIX hosts
    fcntl = None


HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parent.parent
# Detached tmux workers execute this file by absolute path. In that mode
# Python puts ``benchmarks/compare`` (not the repository root) on ``sys.path``,
# so late imports such as ``benchmarks.compare.fourway_compare`` would fail.
# Make script-mode execution equivalent to the documented ``python -m`` form.
if __package__ in (None, "") and str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
DEFAULT_RUNS_ROOT = HERE / "_runs"
DEFAULT_REGISTRY = HERE / "benchmarks.json"
DEFAULT_DATASET_REGISTRY = HERE.parent / "datasets.json"
DEFAULT_BASELINE = HERE / "fourway_results.json"
DEFAULT_PROFILE_REGISTRY = HERE / "campaign_profiles.json"
CONTROLLER_SCHEMA_VERSION = 3
DEFAULT_PAIRED_REPETITIONS = 4
STATUS_EXIT_SUCCESS = 0
STATUS_EXIT_FAILED = 1
STATUS_EXIT_INVALID = 2
STATUS_EXIT_INCOMPLETE = 3
REQUIRED_RUNTIME_DISTRIBUTIONS = (
    "duckdb",
    "gffutils",
    "mappy",
    "numpy",
    "parasail",
    "pyarrow",
    "pyfaidx",
    "pysam",
)
OPTIONAL_RUNTIME_DISTRIBUTIONS = ("lifton",)
PROVENANCE_TOOLING_FILES = {
    "tooling_build_controller": Path(__file__).resolve(),
    "tooling_evaluator": HERE / "evaluator.py",
    "tooling_release_evaluation": HERE / "release_evaluation.py",
    "tooling_release_report": HERE / "release_report.py",
    "tooling_run_benchmarks": HERE.parent / "run_benchmarks.py",
    "tooling_gff3_validator": REPO_ROOT / "lifton" / "gff3_validator.py",
}
WATCHDOG_LOW_CPU_FRACTION = 0.05
WATCHDOG_CONTROL_FILES = {
    ".failed.json", ".success", ".terminal.lock", "cell.json", "exit.json",
    "performance_retry.json", "status.json",
}

SUBSET_CANARIES = ("human_mane", "human_to_zebrafish")
FULL_CANARIES = ("arabidopsis", "human_to_zebrafish")
E2E_DATASETS = ("bee", "human", "arabidopsis", "rice", "mouse")
E2E_CANARIES = ("bee", "human")
PAIRED_STAGE_PREFIX = "paired-"
PAIRED_STAGES = (
    "paired-subset-canary", "paired-subset",
    "paired-full-canary", "paired-full",
    "paired-e2e-canary", "paired-e2e",
)
PAIRED_E2E_MODES = (
    "safe", "stream", "inmemory", "native",
    "stream-inmemory", "stream-native", "inmemory-native", "fast",
)

ACTIVE_STATES = {"launching", "running"}
PENDING_STATES = {"pending", "retry_pending"}
# How long a launched worker may remain unobservable before liveness probes
# stop giving it the benefit of the doubt. Shared by admission control and
# orphan detection so the two can never disagree about whether a cell counts.
LAUNCH_GRACE_SECONDS = 60.0


@dataclass(frozen=True)
class Policy:
    """Resource policy applied before every cell launch."""

    threads_per_cell: int = 8
    scheduler_threads_per_cell: int | None = None
    max_active: int = 4
    max_full: int = 2
    max_worker_threads: int = 32
    load1_limit: float = 32.0
    min_available_gib: float = 256.0
    stagger_seconds: float = 15.0
    poll_seconds: float = 30.0
    subset_timeout_seconds: float = 3.0 * 60.0 * 60.0
    full_timeout_seconds: float = 8.0 * 60.0 * 60.0
    e2e_timeout_seconds: float = 24.0 * 60.0 * 60.0
    stall_timeout_seconds: float = 2.0 * 60.0 * 60.0
    watchdog_poll_seconds: float = 30.0
    terminate_grace_seconds: float = 30.0

    def validate(self) -> None:
        numeric_positive = {
            "threads_per_cell": self.threads_per_cell,
            "max_active": self.max_active,
            "max_full": self.max_full,
            "max_worker_threads": self.max_worker_threads,
            "load1_limit": self.load1_limit,
            "min_available_gib": self.min_available_gib,
            "subset_timeout_seconds": self.subset_timeout_seconds,
            "full_timeout_seconds": self.full_timeout_seconds,
            "e2e_timeout_seconds": self.e2e_timeout_seconds,
            "stall_timeout_seconds": self.stall_timeout_seconds,
            "watchdog_poll_seconds": self.watchdog_poll_seconds,
            "terminate_grace_seconds": self.terminate_grace_seconds,
        }
        if self.scheduler_threads_per_cell is not None:
            numeric_positive["scheduler_threads_per_cell"] = (
                self.scheduler_threads_per_cell
            )
        bad = [name for name, value in numeric_positive.items() if value <= 0]
        if bad:
            raise ValueError(f"policy values must be positive: {', '.join(bad)}")
        if self.max_full > self.max_active:
            raise ValueError("max_full cannot exceed max_active")
        if self.threads_per_cell > self.max_worker_threads:
            raise ValueError("threads_per_cell cannot exceed max_worker_threads")
        if (
            self.scheduler_threads_per_cell is not None
            and self.scheduler_threads_per_cell > self.max_worker_threads
        ):
            raise ValueError(
                "scheduler_threads_per_cell cannot exceed max_worker_threads"
            )


def utc_now() -> str:
    return dt.datetime.now(dt.timezone.utc).isoformat().replace("+00:00", "Z")


def canonical_hash(value: Any) -> str:
    payload = json.dumps(value, sort_keys=True, separators=(",", ":"),
                         ensure_ascii=True).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def atomic_write_json(path: Path, payload: Any) -> None:
    """Durably replace *path* with a JSON document from the same directory."""

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(
        f".{path.name}.{os.getpid()}.{time.time_ns()}.tmp"
    )
    descriptor = -1
    try:
        with temporary.open("w", encoding="utf-8") as handle:
            json.dump(payload, handle, indent=2, sort_keys=True, default=str)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
        try:
            descriptor = os.open(path.parent, os.O_RDONLY)
            os.fsync(descriptor)
        except OSError:
            pass
    finally:
        if descriptor >= 0:
            os.close(descriptor)
        try:
            temporary.unlink()
        except FileNotFoundError:
            pass


def read_json(path: Path) -> Any:
    with Path(path).open(encoding="utf-8") as handle:
        return json.load(handle)


def _run_capture(argv: Sequence[str], *, cwd: Path = REPO_ROOT,
                 timeout: float = 15.0) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        list(argv), cwd=str(cwd), text=True, capture_output=True,
        check=False, timeout=timeout,
    )


def _git_text(args: Sequence[str], repo_root: Path) -> str:
    result = _run_capture(["git", *args], cwd=repo_root, timeout=60.0)
    if result.returncode != 0:
        message = (result.stderr or result.stdout).strip()
        raise RuntimeError(f"git {' '.join(args)} failed: {message}")
    return result.stdout


def _git_diff_hash(repo_root: Path) -> str:
    """Hash the complete tracked working-tree diff without retaining it."""

    process = subprocess.Popen(
        ["git", "diff", "--binary", "HEAD", "--"], cwd=str(repo_root),
        stdout=subprocess.PIPE, stderr=subprocess.PIPE,
    )
    digest = hashlib.sha256()
    assert process.stdout is not None
    for block in iter(lambda: process.stdout.read(1024 * 1024), b""):
        digest.update(block)
    stderr = process.stderr.read().decode("utf-8", "replace") if process.stderr else ""
    returncode = process.wait()
    if returncode != 0:
        raise RuntimeError(f"git diff failed: {stderr.strip()}")
    return digest.hexdigest()


def collect_git_state(repo_root: Path = REPO_ROOT) -> dict[str, Any]:
    head = _git_text(["rev-parse", "HEAD"], repo_root).strip()
    name_status = _git_text(["diff", "--name-status", "HEAD", "--"], repo_root)
    return {
        "head": head,
        "tracked_diff_sha256": _git_diff_hash(repo_root),
        "dirty_tracked_paths": [line for line in name_status.splitlines() if line],
    }


def _resolve_executable(value: str) -> str | None:
    path = Path(value).expanduser()
    if path.is_absolute() or "/" in value:
        return str(path.resolve()) if path.exists() else None
    found = shutil.which(value)
    return str(Path(found).resolve()) if found else None


def probe_tool(name: str, candidate: str) -> dict[str, Any]:
    resolved = _resolve_executable(candidate)
    record: dict[str, Any] = {
        "requested": candidate,
        "path": resolved,
        "available": bool(resolved),
        "version": None,
    }
    if not resolved:
        return record
    try:
        stat = Path(resolved).stat()
        record.update({"size": stat.st_size, "mtime_ns": stat.st_mtime_ns})
    except OSError as exc:
        record["stat_error"] = str(exc)
    try:
        result = _run_capture([resolved, "--version"], timeout=10.0)
        stdout = (result.stdout or "").strip()
        stderr = (result.stderr or "").strip()
        # Successful version commands commonly emit the version on stdout and
        # environment warnings on stderr (for example a random Matplotlib cache
        # path). Prefer stdout so those nondeterministic warnings cannot change
        # an otherwise immutable run fingerprint. Retain stderr-only tools such
        # as tmux by falling back when stdout is empty.
        text = stdout or stderr
        record["version"] = text[:1000]
        record["version_exit_code"] = result.returncode
    except (OSError, subprocess.TimeoutExpired) as exc:
        record["version_error"] = str(exc)
    return record


def collect_runtime_dependencies() -> dict[str, Any]:
    """Fingerprint the shared Python runtime used by every benchmark arm."""

    distributions: dict[str, Any] = {}
    required = set(REQUIRED_RUNTIME_DISTRIBUTIONS)
    for requested in sorted(required | set(OPTIONAL_RUNTIME_DISTRIBUTIONS)):
        try:
            distribution = importlib_metadata.distribution(requested)
        except importlib_metadata.PackageNotFoundError as exc:
            if requested in required:
                raise RuntimeError(
                    "required runtime distribution cannot be fingerprinted: "
                    f"{requested}"
                ) from exc
            continue
        version = str(distribution.version or "").strip()
        if not version:
            raise RuntimeError(
                "runtime distribution has no fingerprintable version: "
                f"{requested}"
            )
        metadata_hashes: dict[str, Any] = {}
        for filename in ("METADATA", "WHEEL", "INSTALLER", "direct_url.json"):
            try:
                text = distribution.read_text(filename)
            except (OSError, UnicodeError) as exc:
                raise RuntimeError(
                    "runtime distribution metadata cannot be fingerprinted: "
                    f"{requested}/{filename}: {exc}"
                ) from exc
            if text is None:
                continue
            encoded = text.encode("utf-8")
            metadata_hashes[filename] = {
                "size": len(encoded),
                "sha256": hashlib.sha256(encoded).hexdigest(),
            }
        distributions[requested] = {
            "name": str(distribution.metadata.get("Name") or requested),
            "version": version,
            "metadata": metadata_hashes,
        }
    return {
        "python": {
            "implementation": sys.implementation.name,
            "version": ".".join(str(part) for part in sys.version_info[:3]),
            "cache_tag": sys.implementation.cache_tag,
            "abi_flags": getattr(sys, "abiflags", ""),
            "executable": str(Path(sys.executable).resolve()),
        },
        "distributions": distributions,
    }


def _registry_tools(registry: Path) -> dict[str, str]:
    try:
        raw = read_json(registry)
    except (OSError, ValueError, TypeError):
        raw = {}
    configured = raw.get("tools", {}) if isinstance(raw, dict) else {}
    tools = {
        name: str(value) for name, value in configured.items()
        if isinstance(value, str) and (name.endswith("_bin") or name.endswith("_python"))
    }
    tools.update({
        "controller_python": sys.executable,
        "tmux": "tmux",
        "git": "git",
        "make": "make",
        "gff3_validate": "gff3-validate",
    })
    return tools


def collect_provenance(
    *,
    repo_root: Path = REPO_ROOT,
    registry: Path = DEFAULT_REGISTRY,
    dataset_registry: Path = DEFAULT_DATASET_REGISTRY,
    baseline: Path = DEFAULT_BASELINE,
    profile_registry: Path | None = None,
) -> dict[str, Any]:
    files: dict[str, Any] = {}
    required_files = [
        ("benchmark_registry", registry),
        ("dataset_registry", dataset_registry),
        ("baseline", baseline),
        *sorted(PROVENANCE_TOOLING_FILES.items()),
    ]
    if profile_registry is not None:
        required_files.append(
            ("campaign_profile_registry", Path(profile_registry))
        )
    for label, path in required_files:
        path = Path(path).resolve()
        if not path.is_file():
            raise FileNotFoundError(f"required {label} does not exist: {path}")
        files[label] = {
            "path": str(path),
            "size": path.stat().st_size,
            "sha256": sha256_file(path),
        }
    tools = {
        name: probe_tool(name, command)
        for name, command in sorted(_registry_tools(Path(registry)).items())
    }
    document = {
        "git": collect_git_state(repo_root),
        "files": files,
        "tools": tools,
        "runtime": collect_runtime_dependencies(),
    }
    document["fingerprint"] = canonical_hash(document)
    return document


def _baseline_keys(path: Path, prefix: str) -> list[str]:
    raw = read_json(path)
    if not isinstance(raw, dict):
        raise ValueError(f"baseline must be a JSON object: {path}")
    return [key.split(":", 1)[1] for key in raw if key.startswith(prefix + ":")]


def _dataset_ids(path: Path) -> list[str]:
    raw = read_json(path)
    entries = raw.get("datasets", []) if isinstance(raw, dict) else []
    return [
        str(item["id"]) for item in entries
        if isinstance(item, dict) and item.get("id")
    ]


def _benchmark_ids(path: Path) -> list[str]:
    raw = read_json(path)
    entries = raw.get("benchmarks", []) if isinstance(raw, dict) else []
    return [
        str(item["id"]) for item in entries
        if isinstance(item, dict) and item.get("id")
    ]


def _base_stage(stage: str) -> str:
    return stage[len(PAIRED_STAGE_PREFIX):] if stage.startswith(PAIRED_STAGE_PREFIX) else stage


def _paired_panel(stage: str) -> str | None:
    if not stage.startswith(PAIRED_STAGE_PREFIX):
        return None
    base = _base_stage(stage)
    return base.split("-", 1)[0]


def select_ids(stage: str, *, baseline: Path, dataset_registry: Path,
               requested: Sequence[str] | None = None,
               benchmark_registry: Path | None = None) -> list[str]:
    stage = _base_stage(stage)
    if requested:
        selected = list(dict.fromkeys(requested))
    elif stage == "subset-canary":
        selected = list(SUBSET_CANARIES)
    elif stage == "subset":
        selected = _baseline_keys(baseline, "subset")
    elif stage == "full-canary":
        selected = list(FULL_CANARIES)
    elif stage == "full":
        selected = _baseline_keys(baseline, "full")
    elif stage == "e2e-canary":
        selected = list(E2E_CANARIES)
    elif stage == "e2e":
        selected = list(E2E_DATASETS)
    elif stage == "gates":
        selected = ["python-suite", "fast-gate", "benchmark-gate"]
    else:
        raise ValueError(f"unsupported stage: {stage}")

    if stage.startswith("subset"):
        available = set(
            _benchmark_ids(benchmark_registry)
            if requested and benchmark_registry is not None
            else _baseline_keys(baseline, "subset")
        )
    elif stage.startswith("full"):
        available = set(
            _benchmark_ids(benchmark_registry)
            if requested and benchmark_registry is not None
            else _baseline_keys(baseline, "full")
        )
    elif stage.startswith("e2e"):
        available = set(_dataset_ids(dataset_registry))
    else:
        available = set(selected)
    unknown = [item for item in selected if item not in available]
    if unknown:
        raise ValueError(f"ids are not supported by stage {stage}: {unknown}")
    return selected


def safe_name(value: str, *, limit: int = 70) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9_.-]+", "-", value).strip("-.") or "cell"
    if len(cleaned) <= limit:
        return cleaned
    suffix = hashlib.sha256(cleaned.encode()).hexdigest()[:10]
    return f"{cleaned[:limit - 11]}-{suffix}"


def paired_configuration(
    *,
    stage: str,
    candidate_root: Path,
    candidate_sha: str,
    reference_root: Path,
    reference_sha: str,
    repetitions: int,
    lifton_executable: Path,
    candidate_e2e_mode: str,
    reference_e2e_mode: str,
    benchmark_registry: Path = DEFAULT_REGISTRY,
    dataset_registry: Path = DEFAULT_DATASET_REGISTRY,
) -> dict[str, Any]:
    """Validate and normalize the immutable paired-source configuration."""

    panel = _paired_panel(stage)
    if panel not in {"subset", "full", "e2e"}:
        raise ValueError(f"stage {stage!r} is not a paired benchmark stage")
    if not 1 <= repetitions <= 10:
        raise ValueError("paired repetitions must be between 1 and 10")
    if not _base_stage(stage).endswith("-canary") and repetitions % 2:
        raise ValueError(
            "non-canary paired stages require an even repetition count for "
            "balanced reference/candidate order"
        )
    for label, sha in (("candidate", candidate_sha), ("reference", reference_sha)):
        if not re.fullmatch(r"[0-9a-fA-F]{40}", str(sha)):
            raise ValueError(f"{label} SHA must be an exact 40-character Git SHA")
    for label, mode in (
        ("candidate", candidate_e2e_mode),
        ("reference", reference_e2e_mode),
    ):
        if mode not in PAIRED_E2E_MODES:
            raise ValueError(
                f"{label} E2E mode must be one of {', '.join(PAIRED_E2E_MODES)}"
            )
    return {
        "panel": panel,
        "repetitions": int(repetitions),
        "exclusive": False,
        "validation_policy": "record_invalid",
        "lifton_executable": str(Path(lifton_executable).resolve()),
        "registries": {
            "benchmark": str(Path(benchmark_registry).resolve()),
            "dataset": str(Path(dataset_registry).resolve()),
        },
        "candidate": {
            "root": str(Path(candidate_root).resolve()),
            "sha": str(candidate_sha).lower(),
            "e2e_mode": candidate_e2e_mode,
        },
        "reference": {
            "root": str(Path(reference_root).resolve()),
            "sha": str(reference_sha).lower(),
            "e2e_mode": reference_e2e_mode,
        },
    }


def _paired_source_specs(configuration: Mapping[str, Any]) -> dict[str, Any]:
    from benchmarks.compare import release_evaluation

    executable = Path(configuration["lifton_executable"])
    return {
        label: release_evaluation.SourceSpec(
            label=label,
            root=Path(configuration[label]["root"]),
            sha=str(configuration[label]["sha"]),
            lifton_executable=executable,
        )
        for label in ("candidate", "reference")
    }


def _prepare_paired_provenance(
    configuration: Mapping[str, Any],
    benchmark_ids: Sequence[str],
) -> dict[str, Any]:
    """Verify both sources and hash every canonical input before planning."""

    from benchmarks.compare import release_evaluation

    sources = {
        label: release_evaluation.verify_source(spec)
        for label, spec in _paired_source_specs(configuration).items()
    }
    inputs: dict[str, Any] = {}
    evaluation_inputs: dict[str, Any] = {}
    for benchmark in benchmark_ids:
        resolved = release_evaluation.resolve_panel_inputs(
            str(configuration["panel"]), benchmark,
            benchmark_registry=Path(configuration["registries"]["benchmark"]),
            dataset_registry=Path(configuration["registries"]["dataset"]),
        )
        record_groups = {
            "inputs": release_evaluation.input_fingerprints(resolved),
            "evaluation_inputs": (
                release_evaluation.evaluation_input_fingerprints(resolved)
            ),
        }
        for records in record_groups.values():
            for record in records.values():
                stable = _stable_paired_input_record(
                    Path(record["path"]),
                    sha256_file=release_evaluation.sha256_file,
                )
                record.clear()
                record.update(stable)
        inputs[benchmark] = record_groups["inputs"]
        evaluation_inputs[benchmark] = record_groups["evaluation_inputs"]
    return {
        "configuration": dict(configuration),
        "sources": sources,
        "inputs": inputs,
        "evaluation_inputs": evaluation_inputs,
    }


def _current_paired_provenance(plan: Mapping[str, Any]) -> dict[str, Any]:
    """Recheck sources and content-address every frozen paired input."""

    from benchmarks.compare import release_evaluation

    configuration = plan["paired"]
    sources = {
        label: release_evaluation.verify_source(spec)
        for label, spec in _paired_source_specs(configuration).items()
    }
    frozen = plan["provenance"]["paired"]
    current_groups: dict[str, Any] = {}
    for group_name in ("inputs", "evaluation_inputs"):
        expected_group = frozen[group_name]
        current_group: dict[str, Any] = {}
        for benchmark, expected_records in expected_group.items():
            current_records: dict[str, Any] = {}
            for name, expected in expected_records.items():
                path = Path(expected["path"])
                try:
                    current_records[name] = _stable_paired_input_record(
                        path,
                        sha256_file=release_evaluation.sha256_file,
                    )
                except (OSError, RuntimeError, ValueError) as exc:
                    current_records[name] = {
                        "path": str(path), "error": str(exc),
                    }
            current_group[benchmark] = current_records
        current_groups[group_name] = current_group
    return {
        "configuration": dict(configuration),
        "sources": sources,
        **current_groups,
    }


def _stable_paired_input_record(
    path: Path,
    *,
    sha256_file: Callable[[Path], str],
) -> dict[str, Any]:
    """Hash one paired input while rejecting a concurrent replacement."""

    resolved = Path(path).resolve()
    before = resolved.stat()
    digest = sha256_file(resolved)
    after = resolved.stat()
    if any(
        getattr(before, field_name) != getattr(after, field_name)
        for field_name in (
            "st_size", "st_mtime_ns", "st_ctime_ns", "st_dev", "st_ino",
        )
    ):
        raise RuntimeError(f"paired input changed while it was hashed: {resolved}")
    return {
        "path": str(resolved),
        "size": after.st_size,
        "sha256": digest,
    }


def _add_paired_provenance(
    provenance: Mapping[str, Any],
    paired: Mapping[str, Any],
) -> dict[str, Any]:
    document = {
        key: value for key, value in provenance.items() if key != "fingerprint"
    }
    document["paired"] = dict(paired)
    document["fingerprint"] = canonical_hash(document)
    return document


def _validator_command(gff_path: Path) -> list[str]:
    return [
        sys.executable, "-m", "lifton.gff3_validator", str(gff_path), "--json",
    ]


def _gate_cell(name: str, cell_dir: Path) -> dict[str, Any]:
    if name == "python-suite":
        command = [sys.executable, "-m", "pytest", "tests/", "-q"]
    elif name == "fast-gate":
        command = ["make", "test-fast", f"LIFTON_PY={sys.executable}"]
    elif name == "benchmark-gate":
        command = ["make", "benchmark-gate", f"LIFTON_PY={sys.executable}"]
    else:
        raise ValueError(f"unknown gate: {name}")
    return {
        "id": f"gate__{safe_name(name)}",
        "kind": "gate",
        "benchmark": name,
        "mode": "gate",
        "threads": 1,
        "full_job": False,
        "command": command,
        "environment": {},
        "artifacts": {},
        "cell_dir": str(cell_dir),
    }


def _subset_cell(benchmark: str, cell_dir: Path, threads: int) -> dict[str, Any]:
    output_root = cell_dir / "work" / benchmark / "_fourway" / "lifton_devel"
    result_json = cell_dir / "result.json"
    gff = output_root / "devel.gff3"
    return {
        "id": f"subset__{safe_name(benchmark)}",
        "kind": "subset",
        "benchmark": benchmark,
        "mode": "subset",
        "threads": threads,
        "full_job": False,
        "command": [
            sys.executable, str(Path(__file__).resolve()), "_run-subset",
            "--benchmark", benchmark, "--cell-dir", str(cell_dir),
            "--threads", str(threads),
        ],
        "environment": {},
        "artifacts": {
            "result_json": str(result_json),
            "result_key": f"subset:{benchmark}",
            "gff": str(gff),
            "manifest": str(output_root / "lifton_output" / "run_manifest.json"),
            "gff_validator": _validator_command(gff),
        },
        "cell_dir": str(cell_dir),
    }


def _full_cell(benchmark: str, cell_dir: Path, threads: int) -> dict[str, Any]:
    output_root = cell_dir / "work" / benchmark / "_devel_refresh"
    result_json = cell_dir / f"{benchmark}.json"
    gff = output_root / "devel_refresh.gff3"
    return {
        "id": f"full__{safe_name(benchmark)}",
        "kind": "full_refresh",
        "benchmark": benchmark,
        "mode": "full",
        "threads": threads,
        "full_job": True,
        "command": [
            sys.executable, str(Path(__file__).resolve()), "_run-refresh",
            "--benchmark", benchmark, "--result-dir", str(cell_dir),
            "--threads", str(threads),
        ],
        "environment": {},
        "artifacts": {
            "result_json": str(result_json),
            "result_key": f"full:{benchmark}",
            "gff": str(gff),
            "manifest": str(output_root / "lifton_output" / "run_manifest.json"),
            "gff_validator": _validator_command(gff),
        },
        "cell_dir": str(cell_dir),
    }


def _e2e_cell(dataset: str, cell_dir: Path, threads: int,
              dataset_registry: Path) -> dict[str, Any]:
    results_dir = cell_dir / "artifacts"
    output_root = results_dir / dataset
    result_json = cell_dir / "result.json"
    gff = output_root / "lifton.gff3"
    return {
        "id": f"e2e__{safe_name(dataset)}",
        "kind": "end_to_end",
        "benchmark": dataset,
        "mode": "end_to_end",
        "threads": threads,
        "full_job": True,
        "command": [
            sys.executable, str(Path(__file__).resolve()), "_run-e2e",
            "--dataset", dataset, "--cell-dir", str(cell_dir),
            "--registry", str(dataset_registry), "--threads", str(threads),
        ],
        "environment": {},
        "artifacts": {
            "result_json": str(result_json),
            "gff": str(gff),
            "manifest": str(output_root / "lifton_output" / "run_manifest.json"),
            "gff_validator": _validator_command(gff),
        },
        "cell_dir": str(cell_dir),
    }


def _paired_cell(
    benchmark: str,
    cell_dir: Path,
    threads: int,
    *,
    configuration: Mapping[str, Any],
    repetition: int,
    input_fingerprints: Mapping[str, Any],
    evaluation_input_fingerprints: Mapping[str, Any],
) -> dict[str, Any]:
    panel = str(configuration["panel"])
    pair_root = cell_dir / "pair"
    if panel == "e2e":
        candidate_gff = pair_root / "candidate" / "artifacts" / benchmark / "lifton.gff3"
        reference_gff = pair_root / "reference" / "artifacts" / benchmark / "lifton.gff3"
    else:
        candidate_gff = pair_root / "candidate" / "candidate.gff3"
        reference_gff = pair_root / "reference" / "reference.gff3"
    candidate_manifest = pair_root / "candidate" / "release_run_manifest.json"
    reference_manifest = pair_root / "reference" / "release_run_manifest.json"
    expected_order = (
        ["reference", "candidate"]
        if repetition % 2
        else ["candidate", "reference"]
    )
    executable = str(configuration["lifton_executable"])
    command = [
        sys.executable, "-m", "benchmarks.compare.release_evaluation", "run-pair",
        "--panel", panel,
        "--benchmark", benchmark,
        "--repetition", str(repetition),
        "--candidate-root", str(configuration["candidate"]["root"]),
        "--candidate-sha", str(configuration["candidate"]["sha"]),
        "--reference-root", str(configuration["reference"]["root"]),
        "--reference-sha", str(configuration["reference"]["sha"]),
        "--lifton-executable", executable,
        "--cell-dir", str(pair_root),
        "--threads", str(threads),
        "--benchmark-registry", str(configuration["registries"]["benchmark"]),
        "--dataset-registry", str(configuration["registries"]["dataset"]),
        "--candidate-e2e-mode", str(configuration["candidate"]["e2e_mode"]),
        "--reference-e2e-mode", str(configuration["reference"]["e2e_mode"]),
    ]
    return {
        "id": (
            f"paired_{panel}__{safe_name(benchmark)}"
            f"__repetition_{repetition:02d}"
        ),
        "kind": "paired_release",
        "benchmark": benchmark,
        "mode": f"paired_{panel}",
        "panel": panel,
        "repetition": repetition,
        "expected_order": expected_order,
        "threads": threads,
        "full_job": panel in {"full", "e2e"},
        "exclusive": bool(configuration.get("exclusive", False)),
        "command": command,
        "environment": {},
        "paired": dict(configuration),
        "input_fingerprints": dict(input_fingerprints),
        "evaluation_input_fingerprints": dict(
            evaluation_input_fingerprints
        ),
        "artifacts": {
            "result_json": str(pair_root / "pair_result.json"),
            "result_key": f"paired:{panel}:{benchmark}:repetition-{repetition:02d}",
            "candidate_gff": str(candidate_gff),
            "reference_gff": str(reference_gff),
            "candidate_manifest": str(candidate_manifest),
            "reference_manifest": str(reference_manifest),
            "candidate_gff_validator": _validator_command(candidate_gff),
            "reference_gff_validator": _validator_command(reference_gff),
        },
        "cell_dir": str(cell_dir),
    }


def build_cells(stage: str, ids: Sequence[str], *, run_dir: Path,
                policy: Policy, dataset_registry: Path,
                paired: Mapping[str, Any] | None = None,
                paired_inputs: Mapping[str, Any] | None = None,
                paired_evaluation_inputs: Mapping[str, Any] | None = None,
                ) -> list[dict[str, Any]]:
    cells: list[dict[str, Any]] = []
    panel = _paired_panel(stage)
    if panel is not None:
        if (
            paired is None
            or paired_inputs is None
            or paired_evaluation_inputs is None
        ):
            raise ValueError("paired stages require immutable source and input configuration")
        for item in ids:
            for repetition in range(1, int(paired["repetitions"]) + 1):
                cell_id = (
                    f"paired_{panel}__{safe_name(item)}"
                    f"__repetition_{repetition:02d}"
                )
                cells.append(_paired_cell(
                    item,
                    run_dir / "cells" / cell_id,
                    policy.threads_per_cell,
                    configuration=paired,
                    repetition=repetition,
                    input_fingerprints=paired_inputs[item],
                    evaluation_input_fingerprints=(
                        paired_evaluation_inputs[item]
                    ),
                ))
        for cell in cells:
            cell["scheduler_thread_cost"] = (
                policy.scheduler_threads_per_cell or cell["threads"]
            )
        return cells
    for item in ids:
        prefix = "gate" if stage == "gates" else stage.split("-", 1)[0]
        cell_id = f"{prefix}__{safe_name(item)}"
        cell_dir = run_dir / "cells" / cell_id
        if stage == "gates":
            cell = _gate_cell(item, cell_dir)
        elif stage.startswith("subset"):
            cell = _subset_cell(item, cell_dir, policy.threads_per_cell)
        elif stage.startswith("full"):
            cell = _full_cell(item, cell_dir, policy.threads_per_cell)
        elif stage.startswith("e2e"):
            cell = _e2e_cell(item, cell_dir, policy.threads_per_cell,
                             dataset_registry)
        else:
            raise ValueError(stage)
        cells.append(cell)
    for cell in cells:
        cell["scheduler_thread_cost"] = (
            policy.scheduler_threads_per_cell or cell["threads"]
        )
    return cells


def _cell_fingerprint(cell: Mapping[str, Any], provenance_fingerprint: str) -> str:
    material = {
        "provenance": provenance_fingerprint,
        "cell": {
            key: value for key, value in cell.items()
            if key != "fingerprint"
        },
    }
    return canonical_hash(material)


def _plan_fingerprint(plan: Mapping[str, Any]) -> str:
    return canonical_hash({
        key: value for key, value in plan.items()
        if key != "fingerprint"
    })


def validate_plan_integrity(plan: Mapping[str, Any]) -> None:
    """Recompute every immutable digest and paired redundancy contract."""

    provenance = plan.get("provenance")
    if not isinstance(provenance, Mapping):
        raise ValueError("plan provenance is missing or malformed")
    provenance_fingerprint = provenance.get("fingerprint")
    if not isinstance(provenance_fingerprint, str):
        raise ValueError("plan provenance fingerprint is missing")
    expected_provenance = canonical_hash({
        key: value for key, value in provenance.items()
        if key != "fingerprint"
    })
    if provenance_fingerprint != expected_provenance:
        raise ValueError("plan provenance fingerprint does not match its content")

    cells = plan.get("cells")
    if not isinstance(cells, list) or not cells:
        raise ValueError("plan cells are missing or empty")
    for cell in cells:
        if not isinstance(cell, Mapping):
            raise ValueError("plan contains a malformed cell")
        expected_cell = _cell_fingerprint(cell, provenance_fingerprint)
        if cell.get("fingerprint") != expected_cell:
            raise ValueError(
                f"cell fingerprint does not match immutable content: "
                f"{cell.get('id', '<unknown>')}"
            )

    paired = plan.get("paired")
    frozen_paired = provenance.get("paired")
    if paired is None:
        if frozen_paired is not None:
            raise ValueError("unpaired plan unexpectedly contains paired provenance")
    else:
        if not isinstance(frozen_paired, Mapping):
            raise ValueError("paired plan lacks frozen paired provenance")
        if frozen_paired.get("configuration") != paired:
            raise ValueError(
                "plan paired configuration disagrees with frozen provenance"
            )
        frozen_inputs = frozen_paired.get("inputs")
        frozen_evaluation_inputs = frozen_paired.get("evaluation_inputs")
        if not isinstance(frozen_inputs, Mapping):
            raise ValueError("paired plan lacks frozen input fingerprints")
        if not isinstance(frozen_evaluation_inputs, Mapping):
            raise ValueError(
                "paired plan lacks frozen evaluation-input fingerprints"
            )
        for cell in cells:
            if cell.get("kind") != "paired_release":
                raise ValueError("paired plan contains a non-paired cell")
            if cell.get("paired") != paired:
                raise ValueError(
                    f"cell paired configuration disagrees with plan: {cell['id']}"
                )
            if cell.get("input_fingerprints") != frozen_inputs.get(
                cell.get("benchmark")
            ):
                raise ValueError(
                    f"cell input fingerprints disagree with frozen provenance: "
                    f"{cell['id']}"
                )
            if (
                cell.get("evaluation_input_fingerprints")
                != frozen_evaluation_inputs.get(cell.get("benchmark"))
            ):
                raise ValueError(
                    "cell evaluation-input fingerprints disagree with frozen "
                    f"provenance: {cell['id']}"
                )
            repetition = cell.get("repetition")
            if (
                not isinstance(repetition, int)
                or isinstance(repetition, bool)
                or not 1 <= repetition <= int(paired["repetitions"])
            ):
                raise ValueError(f"cell repetition is invalid: {cell['id']}")
            expected_order = (
                ["reference", "candidate"]
                if repetition % 2 else ["candidate", "reference"]
            )
            if cell.get("expected_order") != expected_order:
                raise ValueError(
                    f"cell expected order is inconsistent: {cell['id']}"
                )
            if cell.get("panel") != paired.get("panel"):
                raise ValueError(f"cell panel is inconsistent: {cell['id']}")

    campaign_case = plan.get("campaign_case")
    frozen_campaign_case = provenance.get("campaign_case")
    if campaign_case is None:
        if frozen_campaign_case is not None:
            raise ValueError(
                "unprofiled plan unexpectedly contains campaign-case provenance"
            )
    else:
        if not isinstance(campaign_case, Mapping):
            raise ValueError("plan campaign_case is malformed")
        if frozen_campaign_case != campaign_case:
            raise ValueError(
                "plan campaign_case disagrees with frozen provenance"
            )
        case = campaign_case.get("case")
        if not isinstance(case, Mapping):
            raise ValueError("plan campaign_case definition is malformed")
        if plan.get("stage") != case.get("stage"):
            raise ValueError("plan stage disagrees with campaign_case")
        if plan.get("ids") != case.get("ids"):
            raise ValueError("plan ids disagree with campaign_case")
        if plan.get("policy", {}).get("threads_per_cell") != case.get("threads"):
            raise ValueError("plan thread policy disagrees with campaign_case")
        for cell in cells:
            if cell.get("campaign_case") != campaign_case:
                raise ValueError(
                    f"cell campaign_case disagrees with plan: {cell['id']}"
                )

    expected_plan = _plan_fingerprint(plan)
    if plan.get("fingerprint") != expected_plan:
        raise ValueError("plan fingerprint does not match immutable content")


def create_plan(
    *,
    run_id: str | None,
    stage: str,
    requested_ids: Sequence[str] | None,
    runs_root: Path,
    repo_root: Path,
    registry: Path,
    dataset_registry: Path,
    baseline: Path,
    policy: Policy,
    paired: Mapping[str, Any] | None = None,
    campaign_case: Mapping[str, Any] | None = None,
    profile_registry: Path | None = None,
) -> tuple[Path, dict[str, Any]]:
    policy.validate()
    panel = _paired_panel(stage)
    if panel is None and paired is not None:
        raise ValueError("paired source configuration is only valid for paired stages")
    if panel is not None:
        if paired is None:
            raise ValueError(
                "paired stages require candidate/reference roots and exact SHAs"
            )
        try:
            paired = paired_configuration(
                stage=stage,
                candidate_root=Path(paired["candidate"]["root"]),
                candidate_sha=str(paired["candidate"]["sha"]),
                reference_root=Path(paired["reference"]["root"]),
                reference_sha=str(paired["reference"]["sha"]),
                repetitions=int(paired["repetitions"]),
                lifton_executable=Path(paired["lifton_executable"]),
                candidate_e2e_mode=str(paired["candidate"]["e2e_mode"]),
                reference_e2e_mode=str(paired["reference"]["e2e_mode"]),
                benchmark_registry=Path(
                    paired.get("registries", {}).get(
                        "benchmark", DEFAULT_REGISTRY,
                    )
                ),
                dataset_registry=Path(
                    paired.get("registries", {}).get(
                        "dataset", DEFAULT_DATASET_REGISTRY,
                    )
                ),
            )
        except (AttributeError, KeyError, TypeError) as exc:
            raise ValueError(
                f"paired source configuration is incomplete: {exc}"
            ) from exc
        expected_registries = {
            "benchmark": str(Path(registry).resolve()),
            "dataset": str(Path(dataset_registry).resolve()),
        }
        if paired["registries"] != expected_registries:
            raise ValueError(
                "paired registry configuration does not match the plan inputs"
            )
    ids = select_ids(
        stage, baseline=baseline, dataset_registry=dataset_registry,
        requested=requested_ids,
        benchmark_registry=registry if campaign_case is not None else None,
    )
    provenance = collect_provenance(
        repo_root=repo_root, registry=registry,
        dataset_registry=dataset_registry, baseline=baseline,
        profile_registry=profile_registry,
    )
    if campaign_case is not None:
        campaign_case = json.loads(json.dumps(campaign_case))
        case = campaign_case.get("case")
        if not isinstance(case, Mapping):
            raise ValueError("campaign_case must contain a normalized case")
        if (
            case.get("stage") != stage
            or case.get("ids") != ids
            or case.get("threads") != policy.threads_per_cell
        ):
            raise ValueError(
                "campaign_case stage, ids, or threads disagree with the plan"
            )
        provenance["campaign_case"] = campaign_case
        provenance["fingerprint"] = canonical_hash({
            key: value for key, value in provenance.items()
            if key != "fingerprint"
        })
    paired_provenance = None
    if paired is not None:
        paired_provenance = _prepare_paired_provenance(paired, ids)
        provenance = _add_paired_provenance(provenance, paired_provenance)
    if run_id is None:
        stamp = dt.datetime.now(dt.timezone.utc).strftime("%Y%m%dT%H%M%SZ")
        run_id = f"{stamp}-{stage}-{provenance['fingerprint'][:8]}"
    run_id = safe_name(run_id, limit=100)
    run_dir = Path(runs_root).resolve() / run_id
    cells = build_cells(
        stage, ids, run_dir=run_dir, policy=policy,
        dataset_registry=dataset_registry, paired=paired,
        paired_inputs=(
            paired_provenance["inputs"]
            if paired_provenance is not None else None
        ),
        paired_evaluation_inputs=(
            paired_provenance["evaluation_inputs"]
            if paired_provenance is not None else None
        ),
    )
    for cell in cells:
        if campaign_case is not None:
            cell["campaign_case"] = campaign_case
        cell["fingerprint"] = _cell_fingerprint(cell, provenance["fingerprint"])
    plan: dict[str, Any] = {
        "schema_version": CONTROLLER_SCHEMA_VERSION,
        "run_id": run_id,
        "created_at": utc_now(),
        "repo_root": str(Path(repo_root).resolve()),
        "run_dir": str(run_dir),
        "stage": stage,
        "ids": ids,
        "policy": asdict(policy),
        "inputs": {
            "registry": str(Path(registry).resolve()),
            "dataset_registry": str(Path(dataset_registry).resolve()),
            "baseline": str(Path(baseline).resolve()),
        },
        "provenance": provenance,
        "cells": cells,
    }
    if paired is not None:
        plan["paired"] = paired
    if campaign_case is not None:
        plan["campaign_case"] = campaign_case
        plan["inputs"]["profile_registry"] = str(
            Path(profile_registry or DEFAULT_PROFILE_REGISTRY).resolve()
        )
    plan["fingerprint"] = _plan_fingerprint(plan)
    return run_dir, plan


def initialize_run(run_dir: Path, plan: Mapping[str, Any]) -> None:
    validate_plan_integrity(plan)
    validate_plan_layout(plan)
    if run_dir.exists() and any(run_dir.iterdir()):
        raise FileExistsError(f"run directory already exists and is not empty: {run_dir}")
    run_dir.mkdir(parents=True, exist_ok=True)
    atomic_write_json(run_dir / "plan.json", plan)
    for cell in plan["cells"]:
        cell_dir = Path(cell["cell_dir"])
        cell_dir.mkdir(parents=True, exist_ok=True)
        atomic_write_json(cell_dir / "cell.json", cell)
        atomic_write_json(cell_dir / "status.json", {
            "schema_version": CONTROLLER_SCHEMA_VERSION,
            "cell_id": cell["id"],
            "fingerprint": cell["fingerprint"],
            "state": "pending",
            "attempts": 0,
            "updated_at": utc_now(),
        })
    atomic_write_json(run_dir / "controller.json", {
        "schema_version": CONTROLLER_SCHEMA_VERSION,
        "run_id": plan["run_id"],
        "fingerprint": plan["fingerprint"],
        "state": "planned",
        "created_at": utc_now(),
    })


def load_plan(run_dir: Path) -> dict[str, Any]:
    plan = read_json(Path(run_dir) / "plan.json")
    if plan.get("schema_version") != CONTROLLER_SCHEMA_VERSION:
        raise ValueError("unsupported or missing controller plan schema")
    validate_plan_integrity(plan)
    validate_plan_layout(plan)
    return plan


def validate_plan_layout(plan: Mapping[str, Any]) -> None:
    """Fail closed if a generated output escapes its immutable cell directory."""

    run_dir = Path(plan["run_dir"]).resolve()
    cells_root = run_dir / "cells"
    for cell in plan["cells"]:
        cell_dir = Path(cell["cell_dir"]).resolve()
        if not cell_dir.is_relative_to(cells_root):
            raise ValueError(f"cell directory escapes the run root: {cell_dir}")
        if cell["kind"] == "gate":
            continue
        if cell["kind"] == "paired_release":
            output_names = (
                "result_json",
                "candidate_gff", "reference_gff",
                "candidate_manifest", "reference_manifest",
            )
        else:
            output_names = ("result_json", "gff", "manifest")
        for name in output_names:
            output = Path(cell["artifacts"][name]).resolve()
            if not output.is_relative_to(cell_dir):
                raise ValueError(
                    f"cell {cell['id']} artifact {name} escapes isolation: {output}"
                )


def collect_current_provenance(plan: Mapping[str, Any]) -> dict[str, Any]:
    inputs = plan["inputs"]
    provenance = collect_provenance(
        repo_root=Path(plan["repo_root"]),
        registry=Path(inputs["registry"]),
        dataset_registry=Path(inputs["dataset_registry"]),
        baseline=Path(inputs["baseline"]),
        profile_registry=(
            Path(inputs["profile_registry"])
            if inputs.get("profile_registry") else None
        ),
    )
    if plan.get("campaign_case") is not None:
        provenance["campaign_case"] = plan["campaign_case"]
        provenance["fingerprint"] = canonical_hash({
            key: value for key, value in provenance.items()
            if key != "fingerprint"
        })
    if plan.get("paired") is not None:
        provenance = _add_paired_provenance(
            provenance, _current_paired_provenance(plan),
        )
    return provenance


def assert_matching_provenance(plan: Mapping[str, Any]) -> None:
    current = collect_current_provenance(plan)
    expected = plan["provenance"]["fingerprint"]
    if current["fingerprint"] != expected:
        raise RuntimeError(
            "run provenance no longer matches its immutable plan; start a new "
            f"run (expected {expected[:12]}, current {current['fingerprint'][:12]})"
        )


def _cell_for(plan: Mapping[str, Any], cell_id: str) -> dict[str, Any]:
    for cell in plan["cells"]:
        if cell["id"] == cell_id:
            return cell
    raise KeyError(f"unknown cell: {cell_id}")


def _read_status(cell: Mapping[str, Any]) -> dict[str, Any]:
    path = Path(cell["cell_dir"]) / "status.json"
    try:
        status = read_json(path)
    except (OSError, ValueError, TypeError):
        status = {
            "cell_id": cell["id"], "fingerprint": cell["fingerprint"],
            "state": "pending", "attempts": 0,
        }
    return status


def _write_status(
    cell: Mapping[str, Any],
    state: str,
    *,
    clear_fields: Sequence[str] = (),
    **fields: Any,
) -> None:
    old = _read_status(cell)
    payload = {
        **old,
        "schema_version": CONTROLLER_SCHEMA_VERSION,
        "cell_id": cell["id"],
        "fingerprint": cell["fingerprint"],
        "state": state,
        "updated_at": utc_now(),
        **fields,
    }
    for field_name in clear_fields:
        payload.pop(field_name, None)
    atomic_write_json(Path(cell["cell_dir"]) / "status.json", payload)


@contextlib.contextmanager
def _terminal_transition(cell: Mapping[str, Any]):
    lock_path = Path(cell["cell_dir"]) / ".terminal.lock"
    lock_path.parent.mkdir(parents=True, exist_ok=True)
    with lock_path.open("a+", encoding="utf-8") as handle:
        if fcntl is not None:
            fcntl.flock(handle.fileno(), fcntl.LOCK_EX)
        try:
            yield
        finally:
            if fcntl is not None:
                fcntl.flock(handle.fileno(), fcntl.LOCK_UN)


def _read_proc_stat(pid: int, proc_root: Path = Path("/proc")) -> dict[str, int] | None:
    """Return the process identity/accounting fields needed by the watchdog."""

    try:
        text = (Path(proc_root) / str(pid) / "stat").read_text(encoding="utf-8")
        # ``comm`` is parenthesized and may itself contain spaces or ``)``.
        fields = text[text.rfind(")") + 2:].split()
        return {
            "pgrp": int(fields[2]),
            "session": int(fields[3]),
            "utime": int(fields[11]),
            "stime": int(fields[12]),
            "start_ticks": int(fields[19]),
        }
    except (OSError, ValueError, IndexError):
        return None


def _pid_identity_matches(pid: int, pgid: int, start_ticks: int | None,
                          proc_root: Path = Path("/proc")) -> bool:
    """Protect orphan cleanup from signalling a recycled PID/process group."""

    if start_ticks is None:
        return False
    record = _read_proc_stat(pid, proc_root)
    return bool(
        record
        and record["pgrp"] == pgid
        and record["session"] == pgid
        and record["start_ticks"] == start_ticks
    )


def _worker_identity_matches(status: Mapping[str, Any]) -> bool:
    """Return whether the recorded cell worker is still the same process."""

    try:
        pid = int(status["worker_pid"])
        start_ticks = int(status["worker_start_ticks"])
    except (KeyError, TypeError, ValueError):
        return False
    record = _read_proc_stat(pid)
    return bool(record and record["start_ticks"] == start_ticks)


def _process_group_cpu_seconds(pgid: int, proc_root: Path = Path("/proc")) -> float | None:
    """Best-effort live CPU usage for every member of a POSIX process group."""

    try:
        clock_ticks = float(os.sysconf("SC_CLK_TCK"))
        total_ticks = 0
        found = False
        for entry in Path(proc_root).iterdir():
            if not entry.name.isdigit():
                continue
            record = _read_proc_stat(int(entry.name), proc_root)
            if record and record["pgrp"] == pgid:
                total_ticks += record["utime"] + record["stime"]
                found = True
        return total_ticks / clock_ticks if found else None
    except (OSError, ValueError):
        return None


def _progress_snapshot(cell_dir: Path) -> tuple[int, int, int]:
    """Cheaply summarize output/log progress without reading large artifacts."""

    newest_ns = 0
    total_size = 0
    files = 0
    try:
        for root, _dirs, names in os.walk(cell_dir):
            for name in names:
                if name in WATCHDOG_CONTROL_FILES or name.endswith(".watchdog.json"):
                    continue
                try:
                    stat = (Path(root) / name).stat()
                except OSError:
                    continue
                newest_ns = max(newest_ns, stat.st_mtime_ns)
                total_size += stat.st_size
                files += 1
    except OSError:
        pass
    return newest_ns, total_size, files


def _terminate_process_group(
    *,
    pid: int,
    pgid: int,
    start_ticks: int | None,
    grace_seconds: float,
    process: subprocess.Popen[Any] | None = None,
) -> dict[str, Any]:
    """Terminate a verified process group, escalating TERM to KILL if needed."""

    cleanup = {
        "identity_verified": False,
        "term_sent": False,
        "kill_sent": False,
        "error": None,
    }
    owned_handle = process is not None
    verified = _pid_identity_matches(pid, pgid, start_ticks)
    # A live Popen handle is itself a safe identity token even if /proc is
    # temporarily unavailable. Orphan cleanup, which has no handle, fails
    # closed unless the persisted PID/start token still matches.
    if not verified and not owned_handle:
        cleanup["error"] = "process identity no longer matches; signal refused"
        return cleanup
    cleanup["identity_verified"] = True

    try:
        os.killpg(pgid, signal.SIGTERM)
        cleanup["term_sent"] = True
    except ProcessLookupError:
        return cleanup
    except OSError as exc:
        cleanup["error"] = f"SIGTERM failed: {exc}"
        return cleanup

    deadline = time.monotonic() + grace_seconds
    while time.monotonic() < deadline:
        if process is not None:
            process.poll()
        try:
            os.killpg(pgid, 0)
        except ProcessLookupError:
            return cleanup
        except PermissionError:
            break
        time.sleep(min(0.1, max(0.0, deadline - time.monotonic())))

    try:
        os.killpg(pgid, signal.SIGKILL)
        cleanup["kill_sent"] = True
    except ProcessLookupError:
        pass
    except OSError as exc:
        cleanup["error"] = f"SIGKILL failed: {exc}"
    if process is not None:
        try:
            process.wait(timeout=max(1.0, grace_seconds))
        except subprocess.TimeoutExpired:
            cleanup["error"] = cleanup["error"] or "process remained after SIGKILL"
    return cleanup


def _hard_timeout_seconds(cell: Mapping[str, Any], policy: Policy) -> float:
    if cell["kind"] == "paired_release":
        panel = cell.get("panel")
        if panel == "e2e":
            base = policy.e2e_timeout_seconds
        elif panel == "full":
            base = policy.full_timeout_seconds
        elif panel == "subset":
            base = policy.subset_timeout_seconds
        else:
            raise ValueError(f"unknown paired panel: {panel!r}")
        return 2.0 * base
    if cell["kind"] == "end_to_end":
        return policy.e2e_timeout_seconds
    if cell["kind"] == "full_refresh":
        return policy.full_timeout_seconds
    return policy.subset_timeout_seconds


class _SupervisorSignal(RuntimeError):
    def __init__(self, signum: int):
        self.signum = signum
        super().__init__(f"worker received signal {signum}")


def _run_supervised(
    *,
    plan: Mapping[str, Any],
    cell: Mapping[str, Any],
    attempt: int,
    environment: Mapping[str, str],
    stdout_path: Path,
    stderr_path: Path,
) -> dict[str, Any]:
    """Run one cell in its own session with hard/stall watchdog enforcement."""

    policy = Policy(**plan["policy"])
    hard_timeout = _hard_timeout_seconds(cell, policy)
    watchdog_path = Path(cell["cell_dir"]) / f"attempt-{attempt:02d}.watchdog.json"
    started = time.monotonic()
    launch_error: str | None = None
    process: subprocess.Popen[Any] | None = None
    pid: int | None = None
    pgid: int | None = None
    start_ticks: int | None = None
    reason: str | None = None
    cleanup: dict[str, Any] = {}
    previous_handlers: dict[int, Any] = {}

    def _signal_handler(signum: int, _frame: Any) -> None:
        raise _SupervisorSignal(signum)

    if threading.current_thread() is threading.main_thread():
        for signum in (signal.SIGHUP, signal.SIGINT, signal.SIGTERM):
            previous_handlers[signum] = signal.getsignal(signum)
            signal.signal(signum, _signal_handler)

    try:
        with stdout_path.open("wb") as stdout, stderr_path.open("wb") as stderr:
            try:
                process = subprocess.Popen(
                    cell["command"], cwd=plan["repo_root"], env=dict(environment),
                    stdout=stdout, stderr=stderr, start_new_session=True,
                )
            except OSError as exc:
                launch_error = str(exc)
                returncode = 127
            else:
                pid = process.pid
                pgid = process.pid
                record = None
                for _ in range(5):
                    record = _read_proc_stat(pid)
                    if record is not None or process.poll() is not None:
                        break
                    time.sleep(0.01)
                start_ticks = record["start_ticks"] if record else None
                _write_status(
                    cell, "running", attempts=attempt,
                    command_pid=pid, command_pgid=pgid,
                    command_start_ticks=start_ticks,
                    hard_timeout_seconds=hard_timeout,
                    stall_timeout_seconds=policy.stall_timeout_seconds,
                )
                last_snapshot = _progress_snapshot(Path(cell["cell_dir"]))
                last_activity = time.monotonic()
                last_sample = last_activity
                last_cpu = _process_group_cpu_seconds(pgid)
                cpu_observable = last_cpu is not None

                while True:
                    now = time.monotonic()
                    remaining_hard = hard_timeout - (now - started)
                    remaining_stall = policy.stall_timeout_seconds - (now - last_activity)
                    wait_limits = [policy.watchdog_poll_seconds, remaining_hard]
                    if cpu_observable:
                        wait_limits.append(remaining_stall)
                    wait_seconds = max(0.01, min(wait_limits))
                    try:
                        returncode = process.wait(timeout=wait_seconds)
                        break
                    except subprocess.TimeoutExpired:
                        pass

                    now = time.monotonic()
                    snapshot = _progress_snapshot(Path(cell["cell_dir"]))
                    current_cpu = _process_group_cpu_seconds(pgid)
                    cpu_fraction = None
                    if current_cpu is None:
                        cpu_observable = False
                    elif last_cpu is not None:
                        interval = max(now - last_sample, 1e-9)
                        cpu_fraction = max(0.0, current_cpu - last_cpu) / interval
                    if snapshot != last_snapshot or (
                        cpu_fraction is not None
                        and cpu_fraction >= WATCHDOG_LOW_CPU_FRACTION
                    ):
                        last_activity = now
                    last_snapshot = snapshot
                    last_sample = now
                    last_cpu = current_cpu
                    atomic_write_json(watchdog_path, {
                        "schema_version": CONTROLLER_SCHEMA_VERSION,
                        "cell_id": cell["id"],
                        "attempt": attempt,
                        "pid": pid,
                        "pgid": pgid,
                        "process_start_ticks": start_ticks,
                        "elapsed_seconds": now - started,
                        "seconds_since_activity": now - last_activity,
                        "group_cpu_seconds": current_cpu,
                        "cpu_fraction": cpu_fraction,
                        "cpu_observable": cpu_observable,
                        "progress": {
                            "newest_mtime_ns": snapshot[0],
                            "total_size": snapshot[1],
                            "file_count": snapshot[2],
                        },
                        "hard_timeout_seconds": hard_timeout,
                        "stall_timeout_seconds": policy.stall_timeout_seconds,
                        "updated_at": utc_now(),
                    })
                    if now - started >= hard_timeout:
                        reason = "hard_timeout"
                        cleanup = _terminate_process_group(
                            pid=pid, pgid=pgid, start_ticks=start_ticks,
                            grace_seconds=policy.terminate_grace_seconds,
                            process=process,
                        )
                        returncode = process.wait()
                        break
                    if (
                        cpu_observable
                        and now - last_activity >= policy.stall_timeout_seconds
                    ):
                        reason = "stall_timeout"
                        cleanup = _terminate_process_group(
                            pid=pid, pgid=pgid, start_ticks=start_ticks,
                            grace_seconds=policy.terminate_grace_seconds,
                            process=process,
                        )
                        returncode = process.wait()
                        break
                if reason is None and pid is not None and pgid is not None:
                    # A well-behaved cell command waits for its descendants.
                    # Clean up any process that escaped that contract before
                    # artifact validation can publish success.
                    try:
                        os.killpg(pgid, 0)
                    except ProcessLookupError:
                        pass
                    else:
                        cleanup = _terminate_process_group(
                            pid=pid, pgid=pgid, start_ticks=start_ticks,
                            grace_seconds=policy.terminate_grace_seconds,
                            process=process,
                        )
    except _SupervisorSignal as exc:
        reason = f"worker_signal_{signal.Signals(exc.signum).name.lower()}"
        if process is not None and pid is not None and pgid is not None:
            cleanup = _terminate_process_group(
                pid=pid, pgid=pgid, start_ticks=start_ticks,
                grace_seconds=policy.terminate_grace_seconds,
                process=process,
            )
            returncode = process.wait()
        else:
            returncode = 128 + exc.signum
    except BaseException:
        if process is not None and pid is not None and pgid is not None:
            _terminate_process_group(
                pid=pid, pgid=pgid, start_ticks=start_ticks,
                grace_seconds=policy.terminate_grace_seconds,
                process=process,
            )
        raise
    finally:
        for signum, handler in previous_handlers.items():
            signal.signal(signum, handler)

    elapsed = time.monotonic() - started
    watchdog = {
        "reason": reason,
        "hard_timeout_seconds": hard_timeout,
        "stall_timeout_seconds": policy.stall_timeout_seconds,
        "elapsed_seconds": elapsed,
        "pid": pid,
        "pgid": pgid,
        "process_start_ticks": start_ticks,
        "cleanup": cleanup,
    }
    atomic_write_json(watchdog_path, {
        "schema_version": CONTROLLER_SCHEMA_VERSION,
        "cell_id": cell["id"],
        "attempt": attempt,
        **watchdog,
        "finished_at": utc_now(),
    })
    return {
        "returncode": returncode,
        "elapsed_seconds": elapsed,
        "launch_error": launch_error,
        "watchdog": watchdog,
        "watchdog_path": str(watchdog_path),
    }


def validate_gff3_structure(path: Path) -> list[str]:
    """Fast streaming sanity check before the full project validator."""

    errors: list[str] = []
    feature_count = 0
    saw_version = False
    try:
        with Path(path).open(encoding="utf-8", errors="replace") as handle:
            for line_number, raw in enumerate(handle, 1):
                line = raw.rstrip("\r\n")
                if not line:
                    continue
                if line.startswith("#"):
                    if feature_count == 0 and line == "##gff-version 3":
                        saw_version = True
                    continue
                if feature_count == 0 and not saw_version:
                    errors.append("first feature appears before '##gff-version 3'")
                feature_count += 1
                columns = line.split("\t")
                if len(columns) != 9:
                    errors.append(f"line {line_number}: expected 9 columns, got {len(columns)}")
                    if len(errors) >= 20:
                        break
                    continue
                try:
                    start, end = int(columns[3]), int(columns[4])
                    if start < 1 or end < start:
                        raise ValueError
                except ValueError:
                    errors.append(f"line {line_number}: invalid coordinates")
                if columns[6] not in {"+", "-", ".", "?"}:
                    errors.append(f"line {line_number}: invalid strand {columns[6]!r}")
                if columns[2] == "CDS" and columns[7] not in {"0", "1", "2"}:
                    errors.append(f"line {line_number}: invalid CDS phase {columns[7]!r}")
                if len(errors) >= 20:
                    break
    except OSError as exc:
        return [f"cannot read GFF3: {exc}"]
    if not saw_version:
        errors.append("missing '##gff-version 3' directive")
    if feature_count == 0:
        errors.append("GFF3 contains no feature records")
    return errors


def _artifact_is_fresh(path: Path, started_ns: int) -> bool:
    try:
        return path.is_file() and path.stat().st_size > 0 and path.stat().st_mtime_ns >= started_ns
    except OSError:
        return False


def _success_artifact_record(
    path: Path,
    *,
    sha256: str | None = None,
) -> dict[str, Any]:
    stat = Path(path).stat()
    return {
        "path": str(Path(path).resolve()),
        "size": stat.st_size,
        "mtime_ns": stat.st_mtime_ns,
        "sha256": sha256 or sha256_file(Path(path)),
    }


def _number(value: Any) -> bool:
    return isinstance(value, (int, float)) and not isinstance(value, bool)


def _finite_number(
    value: Any,
    *,
    positive: bool = False,
    unit_interval: bool = False,
) -> bool:
    if not _number(value) or not math.isfinite(float(value)):
        return False
    if positive and float(value) <= 0:
        return False
    if unit_interval and not 0.0 <= float(value) <= 1.0:
        return False
    return True


def _same_resolved_path(value: Any, expected: Any) -> bool:
    if not isinstance(value, (str, os.PathLike)):
        return False
    try:
        return Path(value).resolve() == Path(expected).resolve()
    except (OSError, TypeError, ValueError):
        return False


def _sha256_text(value: Any) -> bool:
    return isinstance(value, str) and bool(re.fullmatch(r"[0-9a-f]{64}", value))


def _validate_successful_profile(
    profile: Any,
    *,
    label: str,
) -> list[str]:
    if not isinstance(profile, Mapping):
        return [f"{label} is missing or not an object"]
    errors = []
    if profile.get("exit_code") != 0:
        errors.append(f"{label}.exit_code is not 0")
    for field_name in ("wall_clock_seconds", "peak_rss_mb"):
        if not _finite_number(profile.get(field_name), positive=True):
            errors.append(f"{label}.{field_name} is not finite and positive")
    return errors


def _validate_e2e_payload(payload: Mapping[str, Any], *, label: str) -> list[str]:
    from benchmarks.compare import release_evaluation

    try:
        release_evaluation.validate_e2e_biology(payload)
    except (RuntimeError, TypeError, ValueError) as exc:
        return [f"{label}: {exc}"]
    return []


def _validate_paired_source(
    document: Any,
    expected: Mapping[str, Any],
    *,
    label: str,
    require_import: bool,
) -> list[str]:
    if not isinstance(document, Mapping):
        return [f"{label} source provenance is missing or not an object"]
    errors = []
    if document.get("label") != label:
        errors.append(f"{label} source label is inconsistent")
    if not _same_resolved_path(document.get("root"), expected["root"]):
        errors.append(f"{label} source root does not match the plan")
    if document.get("sha") != expected["sha"]:
        errors.append(f"{label} source SHA does not match the plan")
    if not _same_resolved_path(
        document.get("lifton_executable"), expected["lifton_executable"],
    ):
        errors.append(f"{label} LiftOn executable does not match the plan")
    if require_import:
        imported = document.get("imported_package")
        if not isinstance(imported, str):
            errors.append(f"{label} imported package provenance is missing")
        else:
            try:
                within_root = Path(imported).resolve().is_relative_to(
                    Path(expected["root"]).resolve()
                )
            except (OSError, TypeError, ValueError):
                within_root = False
            if not within_root:
                errors.append(f"{label} imported package is outside its source root")
    return errors


def _validate_gff_fingerprint(document: Any, *, label: str) -> list[str]:
    from benchmarks.compare import release_evaluation

    if not isinstance(document, Mapping):
        return [f"{label} GFF3 fingerprints are missing or not an object"]
    errors = []
    for field_name in ("byte_sha256", "semantic_sha256"):
        if not _sha256_text(document.get(field_name)):
            errors.append(f"{label} fingerprint {field_name} is not a SHA-256")
    if document.get("semantic_algorithm") != (
        release_evaluation.SEMANTIC_HASH_ALGORITHM
    ):
        errors.append(f"{label} fingerprint semantic_algorithm is unsupported")
    feature_records = document.get("feature_records")
    if (
        not isinstance(feature_records, int)
        or isinstance(feature_records, bool)
        or feature_records <= 0
    ):
        errors.append(f"{label} fingerprint feature_records is not positive")
    feature_counts = document.get("feature_counts")
    if not isinstance(feature_counts, Mapping) or not feature_counts:
        errors.append(f"{label} fingerprint feature_counts is empty or malformed")
    else:
        malformed = [
            name for name, count in feature_counts.items()
            if (
                not isinstance(name, str)
                or not name
                or not isinstance(count, int)
                or isinstance(count, bool)
                or count < 0
            )
        ]
        if malformed:
            errors.append(f"{label} fingerprint feature_counts has invalid entries")
        elif isinstance(feature_records, int) and sum(feature_counts.values()) != feature_records:
            errors.append(f"{label} feature counts do not sum to feature_records")
    return errors


def _validate_evaluation_artifact(
    cell: Mapping[str, Any],
    version: Mapping[str, Any],
    *,
    label: str,
) -> list[str]:
    artifacts = version.get("evaluation_artifacts")
    if not isinstance(artifacts, Mapping):
        return [f"{label} evaluation_artifacts are missing or malformed"]
    record = artifacts.get("transcripts_tsv")
    if not isinstance(record, Mapping):
        return [f"{label} transcripts_tsv evidence is missing or malformed"]
    errors = []
    expected = (
        Path(cell["artifacts"]["result_json"]).parent
        / "evaluation"
        / f"{label}.transcripts.tsv"
    ).resolve()
    path_value = record.get("path")
    try:
        is_absolute = (
            isinstance(path_value, (str, os.PathLike))
            and Path(path_value).is_absolute()
        )
    except (OSError, TypeError, ValueError):
        is_absolute = False
    if not is_absolute:
        errors.append(f"{label} transcripts_tsv path is not absolute")
    if not _same_resolved_path(path_value, expected):
        errors.append(
            f"{label} transcripts_tsv path does not match the canonical cell path"
        )
    size = record.get("size")
    if (
        not isinstance(size, int)
        or isinstance(size, bool)
        or size <= 0
    ):
        errors.append(f"{label} transcripts_tsv size is not positive")
    if not _sha256_text(record.get("sha256")):
        errors.append(f"{label} transcripts_tsv sha256 is malformed")
    summary = version.get("summary")
    summary_path = (
        summary.get("transcripts_tsv")
        if isinstance(summary, Mapping) else None
    )
    if not _same_resolved_path(summary_path, path_value):
        errors.append(
            f"{label} evaluator summary transcripts_tsv path disagrees with "
            "evaluation_artifacts"
        )
    truth_summary = (
        summary.get("target_truth")
        if isinstance(summary, Mapping) else None
    )
    truth_record = artifacts.get("target_truth")
    if truth_summary is not None or truth_record is not None:
        expected_truth = (
            Path(cell["artifacts"]["result_json"]).parent
            / "evaluation"
            / f"{label}.target_truth.json"
        ).resolve()
        if not isinstance(truth_summary, Mapping):
            errors.append(f"{label} target-truth summary is malformed")
        if not isinstance(truth_record, Mapping):
            errors.append(f"{label} target-truth artifact evidence is missing")
        else:
            if not _same_resolved_path(
                truth_record.get("path"), expected_truth,
            ):
                errors.append(
                    f"{label} target-truth path does not match the canonical "
                    "cell path"
                )
            if (
                not isinstance(truth_record.get("size"), int)
                or isinstance(truth_record.get("size"), bool)
                or truth_record["size"] <= 0
            ):
                errors.append(f"{label} target-truth size is not positive")
            if not _sha256_text(truth_record.get("sha256")):
                errors.append(f"{label} target-truth sha256 is malformed")
    return errors


def _validate_paired_summary(
    summary: Any,
    *,
    benchmark: str,
    label: str,
    profile: Mapping[str, Any],
) -> list[str]:
    if not isinstance(summary, Mapping) or not summary:
        return [f"{label} evaluator summary is missing or empty"]
    errors = []
    if summary.get("benchmark") != benchmark:
        errors.append(f"{label} evaluator benchmark does not match the cell")
    if summary.get("tool") != label:
        errors.append(f"{label} evaluator tool label is inconsistent")
    total = summary.get("n_reference_total")
    if not isinstance(total, int) or isinstance(total, bool) or total <= 0:
        errors.append(f"{label} evaluator has no reference records")
    for field_name in ("completeness_all", "completeness_coding"):
        value = summary.get(field_name)
        if value is not None and not _finite_number(value, unit_interval=True):
            errors.append(f"{label} evaluator {field_name} is not finite in [0, 1]")
    identity = summary.get("protein_identity")
    if not isinstance(identity, Mapping):
        errors.append(f"{label} evaluator protein_identity is malformed")
    else:
        count = identity.get("n")
        mean = identity.get("mean")
        if not isinstance(count, int) or isinstance(count, bool) or count < 0:
            errors.append(f"{label} evaluator protein_identity.n is invalid")
        elif count > 0 and not _finite_number(mean, unit_interval=True):
            errors.append(f"{label} evaluator protein identity is not finite in [0, 1]")
    summary_profile = summary.get("profile")
    if not isinstance(summary_profile, Mapping):
        errors.append(f"{label} evaluator profile is missing")
    else:
        for field_name in ("wall_clock_seconds", "peak_rss_mb"):
            if summary_profile.get(field_name) != profile.get(field_name):
                errors.append(
                    f"{label} evaluator profile {field_name} disagrees with the run profile"
                )
    return errors


def _validate_paired_result(
    cell: Mapping[str, Any],
    raw: Any,
    *,
    allow_invalid_validation: bool = False,
) -> list[str]:
    from benchmarks.compare import release_evaluation

    if not isinstance(raw, Mapping):
        return ["paired result JSON is not an object"]
    errors: list[str] = []
    configuration = cell["paired"]
    panel = cell["panel"]
    if raw.get("schema_version") != release_evaluation.SCHEMA_VERSION:
        errors.append(
            "paired result schema_version is not "
            f"{release_evaluation.SCHEMA_VERSION}"
        )
    if raw.get("panel") != panel:
        errors.append("paired result panel does not match the cell")
    if raw.get("benchmark") != cell["benchmark"]:
        errors.append("paired result benchmark does not match the cell")
    if raw.get("repetition") != cell["repetition"]:
        errors.append("paired result repetition does not match the cell")
    if raw.get("order") != cell["expected_order"]:
        errors.append("paired result AB/BA order does not match the plan")
    if raw.get("threads") != cell["threads"]:
        errors.append("paired result thread count does not match the plan")

    expected_modes = {
        label: (
            configuration[label]["e2e_mode"] if panel == "e2e" else panel
        )
        for label in ("candidate", "reference")
    }
    if raw.get("modes") != expected_modes:
        errors.append("paired result modes do not match the plan")
    result_registries = raw.get("registries")
    expected_registries = configuration["registries"]
    if not isinstance(result_registries, Mapping) or any(
        not _same_resolved_path(
            result_registries.get(name), expected_registries[name],
        )
        for name in ("benchmark", "dataset")
    ):
        errors.append("paired result registries do not match the plan")

    source_expectations = {
        label: {
            **configuration[label],
            "lifton_executable": configuration["lifton_executable"],
        }
        for label in ("candidate", "reference")
    }
    provenance = raw.get("provenance")
    if not isinstance(provenance, Mapping):
        errors.append("paired result provenance is missing or not an object")
    else:
        for label in ("candidate", "reference"):
            errors.extend(_validate_paired_source(
                provenance.get(label), source_expectations[label],
                label=label, require_import=True,
            ))

    expected_inputs = cell["input_fingerprints"]
    inputs = raw.get("inputs")
    if not isinstance(inputs, Mapping):
        errors.append("paired result inputs are missing or not an object")
    elif set(inputs) != set(expected_inputs):
        errors.append("paired result input set does not match the plan")
    else:
        for name, expected in expected_inputs.items():
            record = inputs.get(name)
            if not isinstance(record, Mapping):
                errors.append(f"paired input {name!r} is not an object")
                continue
            if not _same_resolved_path(record.get("path"), expected["path"]):
                errors.append(f"paired input {name!r} path does not match the plan")
            if record.get("size") != expected["size"]:
                errors.append(f"paired input {name!r} size does not match the plan")
            if (
                not _sha256_text(record.get("sha256"))
                or record.get("sha256") != expected["sha256"]
            ):
                errors.append(f"paired input {name!r} hash does not match the plan")

    expected_evaluation_inputs = cell["evaluation_input_fingerprints"]
    evaluation_inputs = raw.get("evaluation_inputs")
    if not isinstance(evaluation_inputs, Mapping):
        errors.append(
            "paired result evaluation_inputs are missing or not an object"
        )
    elif set(evaluation_inputs) != set(expected_evaluation_inputs):
        errors.append(
            "paired result evaluation-input set does not match the plan"
        )
    else:
        for name, expected in expected_evaluation_inputs.items():
            record = evaluation_inputs.get(name)
            if not isinstance(record, Mapping):
                errors.append(
                    f"paired evaluation input {name!r} is not an object"
                )
                continue
            if not _same_resolved_path(record.get("path"), expected["path"]):
                errors.append(
                    f"paired evaluation input {name!r} path does not match "
                    "the plan"
                )
            if record.get("size") != expected["size"]:
                errors.append(
                    f"paired evaluation input {name!r} size does not match "
                    "the plan"
                )
            if (
                not _sha256_text(record.get("sha256"))
                or record.get("sha256") != expected["sha256"]
            ):
                errors.append(
                    f"paired evaluation input {name!r} hash does not match "
                    "the plan"
                )

    versions = raw.get("versions")
    profiles: dict[str, Mapping[str, Any]] = {}
    if not isinstance(versions, Mapping):
        errors.append("paired result versions are missing or not an object")
        versions = {}
    for label in ("candidate", "reference"):
        version = versions.get(label)
        if not isinstance(version, Mapping):
            errors.append(f"paired result lacks {label} version details")
            continue
        errors.extend(_validate_paired_source(
            version.get("source"), source_expectations[label],
            label=label, require_import=False,
        ))
        profile = version.get("profile")
        errors.extend(_validate_successful_profile(
            profile, label=f"{label} profile",
        ))
        if isinstance(profile, Mapping):
            profiles[label] = profile
        if panel == "e2e":
            if version.get("e2e_mode") != configuration[label]["e2e_mode"]:
                errors.append(f"{label} E2E mode does not match the plan")
            errors.extend(_validate_successful_profile(
                version.get("evaluation_profile"),
                label=f"{label} evaluation profile",
            ))
            errors.extend(_validate_e2e_payload(version, label=f"{label} E2E biology"))
            for field_name in ("lifton_flags", "evaluation_flags"):
                if not isinstance(version.get(field_name), list):
                    errors.append(f"{label} E2E {field_name} are missing")
        elif not isinstance(version.get("argv"), list) or not version["argv"]:
            errors.append(f"{label} paired argv is missing")
        expected_gff = cell["artifacts"][f"{label}_gff"]
        if not _same_resolved_path(version.get("output_gff"), expected_gff):
            errors.append(f"{label} output GFF3 path does not match the cell")
        expected_manifest = cell["artifacts"][f"{label}_manifest"]
        if not _same_resolved_path(
            version.get("release_manifest"), expected_manifest,
        ):
            errors.append(f"{label} release manifest path does not match the cell")
        if not isinstance(version.get("native_manifests"), Mapping):
            errors.append(f"{label} native-manifest evidence is missing")
        errors.extend(_validate_gff_fingerprint(
            version.get("fingerprints"), label=label,
        ))
        validation = version.get("validation")
        if not isinstance(validation, Mapping):
            errors.append(f"{label} GFF3 validation is missing")
        elif not allow_invalid_validation:
            if validation.get("is_valid") is not True:
                errors.append(f"{label} GFF3 validation is not valid")
            if validation.get("n_errors") != 0:
                errors.append(f"{label} GFF3 validation reports errors")
        if isinstance(profile, Mapping):
            errors.extend(_validate_paired_summary(
                version.get("summary"),
                benchmark=cell["benchmark"], label=label, profile=profile,
            ))
        errors.extend(_validate_evaluation_artifact(
            cell, version, label=label,
        ))

    # A same-commit pair is a mode/scheduling qualification, not a version
    # comparison. LiftOn's I/O and parallel flags promise byte-identical
    # annotations, so accepting two individually valid but divergent outputs
    # would hide exactly the class of real-genome regression this canary is
    # meant to detect.
    if (
        configuration["candidate"]["sha"]
        == configuration["reference"]["sha"]
        and all(isinstance(versions.get(label), Mapping)
                for label in ("candidate", "reference"))
    ):
        candidate_version = versions["candidate"]
        reference_version = versions["reference"]
        candidate_fingerprints = candidate_version.get("fingerprints")
        reference_fingerprints = reference_version.get("fingerprints")
        if (isinstance(candidate_fingerprints, Mapping)
                and isinstance(reference_fingerprints, Mapping)):
            for field_name in ("byte_sha256", "semantic_sha256"):
                if candidate_fingerprints.get(field_name) != (
                    reference_fingerprints.get(field_name)
                ):
                    errors.append(
                        "same-SHA paired outputs disagree on " + field_name
                    )
        candidate_evidence = candidate_version.get("evaluation_artifacts")
        reference_evidence = reference_version.get("evaluation_artifacts")
        candidate_tsv = (
            candidate_evidence.get("transcripts_tsv")
            if isinstance(candidate_evidence, Mapping) else None
        )
        reference_tsv = (
            reference_evidence.get("transcripts_tsv")
            if isinstance(reference_evidence, Mapping) else None
        )
        if (isinstance(candidate_tsv, Mapping)
                and isinstance(reference_tsv, Mapping)
                and candidate_tsv.get("sha256") != reference_tsv.get("sha256")):
            errors.append(
                "same-SHA paired neutral-evaluator TSVs disagree"
            )

    ratios = raw.get("ratios")
    if not isinstance(ratios, Mapping):
        errors.append("paired result ratios are missing or not an object")
    else:
        ratio_fields = {
            "wall": "wall_clock_seconds",
            "peak_rss": "peak_rss_mb",
        }
        for ratio_name, profile_field in ratio_fields.items():
            value = ratios.get(ratio_name)
            if not _finite_number(value, positive=True):
                errors.append(f"paired ratio {ratio_name} is not finite and positive")
                continue
            if set(profiles) == {"candidate", "reference"}:
                expected = (
                    float(profiles["candidate"][profile_field])
                    / float(profiles["reference"][profile_field])
                )
                if not math.isclose(
                    float(value), expected, rel_tol=1e-12, abs_tol=1e-12,
                ):
                    errors.append(
                        f"paired ratio {ratio_name} disagrees with the profiles"
                    )
    return errors


def validate_result_schema(
    cell: Mapping[str, Any],
    path: Path,
    *,
    allow_invalid_validation: bool = False,
) -> list[str]:
    errors: list[str] = []
    try:
        raw = read_json(path)
    except (OSError, ValueError, TypeError) as exc:
        return [f"result JSON is unreadable: {exc}"]
    kind = cell["kind"]
    benchmark = cell["benchmark"]
    if kind == "subset":
        key = cell["artifacts"]["result_key"]
        record = raw.get(key) if isinstance(raw, dict) else None
        if not isinstance(record, dict):
            return [f"result JSON lacks object {key!r}"]
        if record.get("benchmark") != benchmark or record.get("mode") != "subset":
            errors.append("result benchmark/mode does not match the cell")
        tools = record.get("tools", [])
        if "lifton_devel" not in tools:
            errors.append("lifton_devel result is absent")
        validity = record.get("validity", {}).get("lifton_devel", {})
        if validity.get("valid") is not True:
            errors.append("result reports an invalid lifton_devel GFF3")
        for key_name in ("wall_s", "peak_rss_mb", "completeness_coding", "mean_pi"):
            if not isinstance(record.get(key_name), dict):
                errors.append(f"result field {key_name!r} is missing or not an object")
        for key_name in ("wall_s", "peak_rss_mb", "completeness_coding", "mean_pi"):
            value = (record.get(key_name) or {}).get("lifton_devel")
            if not _number(value):
                errors.append(f"result field {key_name}.lifton_devel is not numeric")
        for key_name in ("wall_s", "peak_rss_mb"):
            value = (record.get(key_name) or {}).get("lifton_devel")
            if _number(value) and value <= 0:
                errors.append(f"result field {key_name}.lifton_devel is not positive")
    elif kind == "full_refresh":
        if not isinstance(raw, dict) or raw.get("bid") != benchmark:
            return ["refresh result benchmark does not match the cell"]
        summary = raw.get("summary")
        if not isinstance(summary, dict):
            errors.append("refresh result lacks summary object")
        else:
            for key_name in ("completeness_coding", "protein_identity", "profile"):
                if key_name not in summary:
                    errors.append(f"refresh summary lacks {key_name!r}")
            profile = summary.get("profile") or {}
            protein = summary.get("protein_identity") or {}
            for label, value in (
                ("completeness_coding", summary.get("completeness_coding")),
                ("protein_identity.mean", protein.get("mean")),
                ("profile.wall_clock_seconds", profile.get("wall_clock_seconds")),
                ("profile.peak_rss_mb", profile.get("peak_rss_mb")),
            ):
                if not _number(value):
                    errors.append(f"refresh summary field {label} is not numeric")
        if raw.get("validity", {}).get("valid") is not True:
            errors.append("refresh result reports an invalid GFF3")
    elif kind == "end_to_end":
        rows = raw.get("rows") if isinstance(raw, dict) else None
        matches = [
            row for row in (rows or [])
            if isinstance(row, Mapping) and row.get("dataset") == benchmark
        ]
        if len(matches) != 1:
            return [f"expected exactly one result row for {benchmark!r}"]
        row = matches[0]
        if row.get("error"):
            errors.append(f"end-to-end harness error: {row['error']}")
        profile = row.get("lift_profile") or {}
        if profile.get("exit_code") != 0:
            errors.append("end-to-end LiftOn profile does not have exit_code 0")
        for key_name in ("wall_clock_seconds", "peak_rss_mb"):
            value = profile.get(key_name)
            if not _number(value) or value <= 0:
                errors.append(f"end-to-end profile field {key_name} is not positive numeric")
        eval_profile = row.get("eval_profile")
        if not isinstance(eval_profile, dict) or eval_profile.get("exit_code") != 0:
            errors.append("end-to-end evaluation profile does not have exit_code 0")
        elif any(
            not _finite_number(eval_profile.get(field_name), positive=True)
            for field_name in ("wall_clock_seconds", "peak_rss_mb")
        ):
            errors.append("end-to-end evaluation profile metrics are not finite and positive")
        errors.extend(_validate_e2e_payload(row, label="end-to-end biology"))
    elif kind == "paired_release":
        errors.extend(_validate_paired_result(
            cell, raw, allow_invalid_validation=allow_invalid_validation,
        ))
    elif kind != "gate":
        errors.append(f"unknown result kind {kind!r}")
    return errors


def validate_manifest(path: Path) -> list[str]:
    try:
        manifest = read_json(path)
    except (OSError, ValueError, TypeError) as exc:
        return [f"run manifest is unreadable: {exc}"]
    errors = []
    if not isinstance(manifest, dict) or not manifest.get("schema_version"):
        errors.append("run manifest lacks schema_version")
    run = manifest.get("run", {}) if isinstance(manifest, dict) else {}
    if run.get("status") != "success":
        errors.append(f"run manifest status is {run.get('status')!r}, not 'success'")
    if not run.get("finished_at"):
        errors.append("run manifest lacks finished_at")
    return errors


def validate_release_arm_manifest(
    cell: Mapping[str, Any],
    label: str,
    path: Path,
    *,
    version: Mapping[str, Any],
    observed_fingerprints: Mapping[str, Any],
    allow_invalid_validation: bool = False,
) -> list[str]:
    """Cross-check one neutral arm manifest against plan and live artifacts."""

    from benchmarks.compare import release_evaluation

    try:
        manifest = read_json(path)
    except (OSError, TypeError, ValueError) as exc:
        return [f"{label} release manifest is unreadable: {exc}"]
    if not isinstance(manifest, Mapping):
        return [f"{label} release manifest is not an object"]
    errors = []
    if manifest.get("schema_version") != release_evaluation.SCHEMA_VERSION:
        errors.append(
            f"{label} release manifest schema_version is not "
            f"{release_evaluation.SCHEMA_VERSION}"
        )
    if manifest.get("kind") != "paired_release_arm":
        errors.append(f"{label} release manifest kind is not paired_release_arm")
    configuration = cell["paired"]
    expected_source = {
        **configuration[label],
        "lifton_executable": configuration["lifton_executable"],
    }
    errors.extend(_validate_paired_source(
        manifest.get("source"), expected_source,
        label=label, require_import=False,
    ))

    protocol = manifest.get("protocol")
    if not isinstance(protocol, Mapping):
        errors.append(f"{label} release manifest protocol is missing")
    else:
        if protocol.get("kind") != cell["panel"]:
            errors.append(f"{label} release manifest protocol kind is inconsistent")
        if cell["panel"] == "e2e":
            if protocol.get("mode") != configuration[label]["e2e_mode"]:
                errors.append(f"{label} release manifest E2E mode is inconsistent")
            for field_name in ("lifton_flags", "evaluation_flags"):
                if protocol.get(field_name) != version.get(field_name):
                    errors.append(
                        f"{label} release manifest {field_name} "
                        "disagrees with pair result"
                    )
        elif protocol.get("argv") != version.get("argv"):
            errors.append(
                f"{label} release manifest argv disagrees with pair result"
            )

    profile = manifest.get("profile")
    errors.extend(_validate_successful_profile(
        profile, label=f"{label} release manifest profile",
    ))
    if profile != version.get("profile"):
        errors.append(f"{label} release manifest profile disagrees with pair result")

    manifest_artifacts = manifest.get("artifacts")
    if not isinstance(manifest_artifacts, Mapping):
        errors.append(f"{label} release manifest artifacts are missing")
        manifest_artifacts = {}
    output = manifest_artifacts.get("output_gff")
    expected_gff = Path(cell["artifacts"][f"{label}_gff"])
    if not isinstance(output, Mapping):
        errors.append(f"{label} release manifest output GFF3 is missing")
    else:
        if not _same_resolved_path(output.get("path"), expected_gff):
            errors.append(f"{label} release manifest GFF3 path is inconsistent")
        if output.get("size") != expected_gff.stat().st_size:
            errors.append(f"{label} release manifest GFF3 size is inconsistent")
        for manifest_key, fingerprint_key in (
            ("byte_sha256", "byte_sha256"),
            ("semantic_sha256", "semantic_sha256"),
        ):
            if output.get(manifest_key) != observed_fingerprints.get(fingerprint_key):
                errors.append(
                    f"{label} release manifest {manifest_key} is inconsistent"
                )

    validation = manifest.get("validation")
    if validation != version.get("validation"):
        errors.append(f"{label} release manifest validation disagrees with pair result")
    if (
        not allow_invalid_validation
        and (
            not isinstance(validation, Mapping)
            or validation.get("is_valid") is not True
            or validation.get("n_errors") != 0
        )
    ):
        errors.append(f"{label} release manifest validation is not successful")

    native = manifest_artifacts.get("native_manifests")
    if native != version.get("native_manifests"):
        errors.append(
            f"{label} release manifest native evidence disagrees with pair result"
        )
    if not isinstance(native, Mapping):
        errors.append(f"{label} release manifest native evidence is missing")
    else:
        for name, record in native.items():
            if not isinstance(name, str) or not isinstance(record, Mapping):
                errors.append(f"{label} native-manifest evidence is malformed")
                continue
            if not isinstance(record.get("path"), str):
                errors.append(f"{label} native manifest {name!r} lacks a path")
            if not isinstance(record.get("present"), bool):
                errors.append(
                    f"{label} native manifest {name!r} lacks a presence flag"
                )
            elif record["present"] and (
                not isinstance(record.get("size"), int)
                or record["size"] <= 0
                or not _sha256_text(record.get("sha256"))
            ):
                errors.append(
                    f"{label} native manifest {name!r} fingerprint is malformed"
                )
    return errors


def _run_gff_validator(cell: Mapping[str, Any], gff: Path) -> tuple[list[str], dict[str, Any]]:
    template = cell["artifacts"].get("gff_validator", [])
    command = [str(part).replace("{gff}", str(gff)) for part in template]
    if not command:
        return ["no full GFF3 validator command is configured"], {}
    try:
        result = subprocess.run(
            command, cwd=str(REPO_ROOT), text=True, capture_output=True,
            check=False,
        )
    except OSError as exc:
        return [f"could not launch full GFF3 validator: {exc}"], {}
    parsed: dict[str, Any] = {}
    try:
        parsed = json.loads(result.stdout)
    except (ValueError, TypeError):
        pass
    report = {
        "command": command,
        "exit_code": result.returncode,
        "result": parsed,
        "stdout_tail": (result.stdout or "")[-4000:],
        "stderr_tail": (result.stderr or "")[-4000:],
        "checked_at": utc_now(),
    }
    errors = []
    if result.returncode != 0:
        errors.append(f"full GFF3 validator exited {result.returncode}")
    if parsed.get("is_valid") is not True:
        errors.append("full GFF3 validator did not report is_valid=true")
    return errors, report


def _validate_paired_artifacts(
    cell: Mapping[str, Any],
    started_ns: int,
    *,
    persist_validator_reports: bool = True,
) -> tuple[list[str], dict[str, Any]]:
    from benchmarks.compare import release_evaluation

    artifacts = cell["artifacts"]
    required = {
        name: Path(artifacts[name])
        for name in (
            "result_json",
            "candidate_gff", "reference_gff",
            "candidate_manifest", "reference_manifest",
        )
    }
    errors = [
        f"{name} is missing, empty, or stale: {path}"
        for name, path in required.items()
        if not _artifact_is_fresh(path, started_ns)
    ]
    if errors:
        return errors, {
            "artifacts": {name: str(path) for name, path in required.items()},
        }

    allow_invalid_validation = (
        cell.get("paired", {}).get("validation_policy") == "record_invalid"
    )
    errors.extend(validate_result_schema(
        cell,
        required["result_json"],
        allow_invalid_validation=allow_invalid_validation,
    ))
    validation_reports: dict[str, Any] = {}
    fingerprints: dict[str, Any] = {}
    try:
        raw = read_json(required["result_json"])
    except (OSError, TypeError, ValueError):
        raw = {}
    versions = raw.get("versions", {}) if isinstance(raw, Mapping) else {}
    evaluation_artifacts: dict[str, Any] = {}
    evidence_paths = dict(required)
    result_mtime_ns = required["result_json"].stat().st_mtime_ns
    for label in ("candidate", "reference"):
        manifest_errors = validate_manifest(required[f"{label}_manifest"])
        errors.extend(f"{label} manifest: {error}" for error in manifest_errors)
        structure_errors = validate_gff3_structure(required[f"{label}_gff"])
        if not allow_invalid_validation:
            errors.extend(f"{label} GFF3: {error}" for error in structure_errors)

        validator_cell = {
            **cell,
            "artifacts": {
                "gff_validator": artifacts[f"{label}_gff_validator"],
            },
        }
        validator_errors, validator_report = _run_gff_validator(
            validator_cell, required[f"{label}_gff"],
        )
        validation_reports[label] = validator_report
        if not allow_invalid_validation:
            errors.extend(f"{label} GFF3: {error}" for error in validator_errors)
        validation_path = (
            Path(cell["cell_dir"]) / f"{label}_gff_validation.json"
        )
        if persist_validator_reports:
            atomic_write_json(validation_path, validator_report)
        if validation_path.is_file():
            evidence_paths[f"{label}_gff_validation"] = validation_path
        else:
            errors.append(
                f"{label} GFF3 validation report is missing: {validation_path}"
            )

        try:
            observed = release_evaluation.gff3_fingerprints(
                required[f"{label}_gff"],
            )
        except (OSError, TypeError, ValueError) as exc:
            errors.append(f"{label} GFF3 fingerprinting failed: {exc}")
            continue
        fingerprints[label] = observed
        version = versions.get(label) if isinstance(versions, Mapping) else None
        recorded = (
            version.get("fingerprints")
            if isinstance(version, Mapping) else None
        )
        if recorded != observed:
            errors.append(
                f"{label} GFF3 fingerprints do not match pair_result.json"
            )
        if isinstance(version, Mapping):
            errors.extend(validate_release_arm_manifest(
                cell,
                label,
                required[f"{label}_manifest"],
                version=version,
                observed_fingerprints=observed,
                allow_invalid_validation=allow_invalid_validation,
            ))
            artifact_group = version.get("evaluation_artifacts")
            record = (
                artifact_group.get("transcripts_tsv")
                if isinstance(artifact_group, Mapping) else None
            )
            expected_tsv = (
                required["result_json"].parent
                / "evaluation"
                / f"{label}.transcripts.tsv"
            ).resolve()
            if not isinstance(record, Mapping):
                errors.append(
                    f"{label} evaluator TSV evidence is missing from pair result"
                )
            elif not _same_resolved_path(record.get("path"), expected_tsv):
                errors.append(
                    f"{label} evaluator TSV path is inconsistent with pair result"
                )
            elif not _artifact_is_fresh(expected_tsv, started_ns):
                errors.append(
                    f"{label} evaluator TSV is missing, empty, or stale: "
                    f"{expected_tsv}"
                )
            else:
                try:
                    stat = expected_tsv.stat()
                    observed_sha = sha256_file(expected_tsv)
                except OSError as exc:
                    errors.append(
                        f"{label} evaluator TSV cannot be fingerprinted: {exc}"
                    )
                else:
                    observed_record = {
                        "path": str(expected_tsv),
                        "size": stat.st_size,
                        "sha256": observed_sha,
                        "mtime_ns": stat.st_mtime_ns,
                    }
                    evaluation_artifacts[label] = {
                        "transcripts_tsv": observed_record,
                    }
                    if record.get("size") != stat.st_size:
                        errors.append(
                            f"{label} evaluator TSV size disagrees with pair result"
                        )
                    if record.get("sha256") != observed_sha:
                        errors.append(
                            f"{label} evaluator TSV hash disagrees with pair result"
                        )
                    if stat.st_mtime_ns > result_mtime_ns:
                        errors.append(
                            f"{label} evaluator TSV changed after pair_result.json "
                            "was published"
                        )
            truth_summary = (
                version.get("summary", {}).get("target_truth")
                if isinstance(version.get("summary"), Mapping) else None
            )
            truth_record = (
                artifact_group.get("target_truth")
                if isinstance(artifact_group, Mapping) else None
            )
            if truth_summary is not None or truth_record is not None:
                expected_truth = (
                    required["result_json"].parent
                    / "evaluation"
                    / f"{label}.target_truth.json"
                ).resolve()
                if not isinstance(truth_summary, Mapping):
                    errors.append(f"{label} target-truth summary is malformed")
                if not isinstance(truth_record, Mapping):
                    errors.append(
                        f"{label} target-truth artifact evidence is missing"
                    )
                elif not _same_resolved_path(
                    truth_record.get("path"), expected_truth,
                ):
                    errors.append(
                        f"{label} target-truth artifact path is inconsistent"
                    )
                elif not _artifact_is_fresh(expected_truth, started_ns):
                    errors.append(
                        f"{label} target-truth artifact is missing, empty, "
                        f"or stale: {expected_truth}"
                    )
                else:
                    try:
                        truth_stat = expected_truth.stat()
                        truth_sha = sha256_file(expected_truth)
                        truth_document = read_json(expected_truth)
                    except (OSError, TypeError, ValueError) as exc:
                        errors.append(
                            f"{label} target-truth artifact cannot be "
                            f"validated: {exc}"
                        )
                    else:
                        evaluation_artifacts.setdefault(label, {})[
                            "target_truth"
                        ] = {
                            "path": str(expected_truth),
                            "size": truth_stat.st_size,
                            "sha256": truth_sha,
                            "mtime_ns": truth_stat.st_mtime_ns,
                        }
                        if (
                            truth_record.get("size") != truth_stat.st_size
                            or truth_record.get("sha256") != truth_sha
                        ):
                            errors.append(
                                f"{label} target-truth fingerprint disagrees "
                                "with live evidence"
                            )
                        if truth_document != truth_summary:
                            errors.append(
                                f"{label} target-truth summary disagrees with "
                                "its artifact"
                            )
                        if truth_stat.st_mtime_ns > result_mtime_ns:
                            errors.append(
                                f"{label} target-truth artifact changed after "
                                "pair_result.json was published"
                            )
            native = version.get("native_manifests")
            if isinstance(native, Mapping):
                arm_root = required["result_json"].parent / label
                for name, native_record in native.items():
                    if (
                        not isinstance(name, str)
                        or not isinstance(native_record, Mapping)
                        or native_record.get("present") is not True
                    ):
                        continue
                    native_path = Path(str(native_record.get("path", ""))).resolve()
                    evidence_name = (
                        f"{label}_native_manifest_{safe_name(name, limit=40)}"
                    )
                    if not native_path.is_relative_to(arm_root.resolve()):
                        errors.append(
                            f"{label} native manifest {name!r} escapes its arm"
                        )
                        continue
                    if not _artifact_is_fresh(native_path, started_ns):
                        errors.append(
                            f"{label} native manifest {name!r} is missing, "
                            f"empty, or stale: {native_path}"
                        )
                        continue
                    observed_native = _success_artifact_record(native_path)
                    if (
                        native_record.get("size") != observed_native["size"]
                        or native_record.get("sha256")
                        != observed_native["sha256"]
                    ):
                        errors.append(
                            f"{label} native manifest {name!r} fingerprint "
                            "disagrees with live evidence"
                        )
                    evidence_paths[evidence_name] = native_path

    stats = {}
    for name, path in evidence_paths.items():
        known_sha = None
        if name == "candidate_gff" and "candidate" in fingerprints:
            known_sha = fingerprints["candidate"]["byte_sha256"]
        elif name == "reference_gff" and "reference" in fingerprints:
            known_sha = fingerprints["reference"]["byte_sha256"]
        stats[name] = _success_artifact_record(path, sha256=known_sha)
    return errors, {
        "artifacts": stats,
        "evaluation_artifacts": evaluation_artifacts,
        "gff_validation": validation_reports,
        "gff_fingerprints": fingerprints,
    }


def validate_artifacts(
    cell: Mapping[str, Any],
    started_ns: int,
    *,
    persist_validator_reports: bool = True,
) -> tuple[list[str], dict[str, Any]]:
    """Validate a cell, optionally keeping validator reruns in memory only.

    Non-persisting mode is for auditing already sealed artifacts. Its returned
    validator details describe the live rerun, while artifact records continue
    to describe the existing sealed reports.
    """

    if cell["kind"] == "gate":
        return [], {"gate": "exit-code only"}
    if cell["kind"] == "paired_release":
        return _validate_paired_artifacts(
            cell,
            started_ns,
            persist_validator_reports=persist_validator_reports,
        )
    artifacts = cell["artifacts"]
    required = {
        name: Path(artifacts[name]) for name in ("result_json", "gff", "manifest")
    }
    errors = []
    for name, path in required.items():
        if not _artifact_is_fresh(path, started_ns):
            errors.append(f"{name} is missing, empty, or stale: {path}")
    if errors:
        return errors, {"artifacts": {name: str(path) for name, path in required.items()}}

    errors.extend(validate_result_schema(cell, required["result_json"]))
    errors.extend(validate_manifest(required["manifest"]))
    errors.extend(validate_gff3_structure(required["gff"]))
    validator_errors, validator_report = _run_gff_validator(cell, required["gff"])
    errors.extend(validator_errors)
    validation_path = Path(cell["cell_dir"]) / "gff_validation.json"
    if persist_validator_reports:
        atomic_write_json(validation_path, validator_report)
    evidence_paths = {
        **required,
    }
    if validation_path.is_file():
        evidence_paths["gff_validation"] = validation_path
    else:
        errors.append(f"GFF3 validation report is missing: {validation_path}")
    stats = {
        name: _success_artifact_record(path)
        for name, path in evidence_paths.items()
    }
    return errors, {"artifacts": stats, "gff_validation": validator_report}


def performance_regressions(cell: Mapping[str, Any], baseline_path: Path) -> list[dict[str, Any]]:
    """Return performance failures against the immutable canonical result.

    Subset wall time is a paired measurement: the same cell runs the frozen
    stable executable immediately before devel. Normalize devel's raw change by
    that control's change so host and input drift cannot flip a near-threshold
    result. Memory and cells without a paired control remain absolute checks.
    """

    if cell["kind"] not in {"subset", "full_refresh"}:
        return []
    baseline = read_json(baseline_path)
    base = baseline.get(cell["artifacts"]["result_key"], {})
    result = read_json(Path(cell["artifacts"]["result_json"]))
    if cell["kind"] == "subset":
        current_record = result.get(cell["artifacts"]["result_key"], {})
        baseline_wall = base.get("wall_s", {})
        current_wall_values = current_record.get("wall_s", {})
        current_wall = current_wall_values.get("lifton_devel")
        current_rss = current_record.get("peak_rss_mb", {}).get("lifton_devel")
    else:
        baseline_wall = {}
        current_wall_values = {}
        profile = result.get("summary", {}).get("profile", {})
        current_wall = profile.get("wall_clock_seconds")
        current_rss = profile.get("peak_rss_mb")
    values = {
        "wall_s": (base.get("wall_s", {}).get("lifton_devel"), current_wall),
        "peak_rss_mb": (base.get("peak_rss_mb", {}).get("lifton_devel"), current_rss),
    }
    regressions = []
    for metric, (old, new) in values.items():
        if not (_number(old) and _number(new)) or old <= 0:
            continue
        raw_ratio = float(new) / float(old)
        ratio = raw_ratio
        comparison = "absolute"
        control = None
        if metric == "wall_s" and cell["kind"] == "subset":
            control_old = baseline_wall.get("lifton_stable")
            control_new = current_wall_values.get("lifton_stable")
            if (_number(control_old) and _number(control_new)
                    and control_old > 0 and control_new > 0):
                control_ratio = float(control_new) / float(control_old)
                ratio = raw_ratio / control_ratio
                comparison = "paired_stable_normalized"
                control = {
                    "tool": "lifton_stable",
                    "baseline": control_old,
                    "current": control_new,
                    "ratio": control_ratio,
                }
        if ratio > 1.25:
            record = {
                "metric": metric, "baseline": old, "current": new,
                "raw_ratio": raw_ratio, "ratio": ratio,
                "comparison": comparison, "threshold": 1.25,
            }
            if control is not None:
                record["control"] = control
            regressions.append(record)
    return regressions


def _unlink_markers(cell_dir: Path) -> None:
    for name in (".success", ".failed.json"):
        try:
            (cell_dir / name).unlink()
        except FileNotFoundError:
            pass


def _relocate_archived_value(value: Any, old_root: Path, new_root: Path) -> Any:
    """Rewrite only absolute paths that moved with one archived pair tree."""

    if isinstance(value, Mapping):
        return {
            key: _relocate_archived_value(item, old_root, new_root)
            for key, item in value.items()
        }
    if isinstance(value, list):
        return [
            _relocate_archived_value(item, old_root, new_root)
            for item in value
        ]
    if isinstance(value, str):
        path = Path(value)
        if path.is_absolute() and path.is_relative_to(old_root):
            return str(new_root / path.relative_to(old_root))
    return value


def _relocate_archived_pair_evidence(
    archive: Path,
    *,
    original_root: Path,
) -> None:
    """Make moved pair/arm evidence self-contained without touching GFF3 bytes."""

    marker = archive / "relocation.json"
    if marker.is_file():
        return
    evidence_paths = [
        archive / "pair_result.json",
        archive / "candidate" / "release_run_manifest.json",
        archive / "reference" / "release_run_manifest.json",
    ]
    records = []
    for path in evidence_paths:
        if not path.is_file():
            continue
        before = sha256_file(path)
        document = read_json(path)
        relocated = _relocate_archived_value(
            document, original_root, archive,
        )
        atomic_write_json(path, relocated)
        records.append({
            "path": str(path),
            "original_sha256": before,
            "relocated_sha256": sha256_file(path),
        })
    atomic_write_json(marker, {
        "schema_version": CONTROLLER_SCHEMA_VERSION,
        "original_root": str(original_root),
        "archive_root": str(archive),
        "relocated_at": utc_now(),
        "documents": records,
        "gff3_bytes_modified": False,
    })


def _prepare_attempt_workspace(cell: Mapping[str, Any], attempt: int) -> None:
    """Archive a failed paired workspace so its immutable command can resume."""

    if cell["kind"] != "paired_release" or attempt <= 1:
        return
    cell_dir = Path(cell["cell_dir"])
    pair_dir = Path(cell["artifacts"]["result_json"]).parent
    archive = cell_dir / f"attempt-{attempt - 1:02d}.pair"
    if archive.exists():
        if pair_dir.exists():
            raise RuntimeError(
                f"paired retry archive already exists while workspace remains: {archive}"
            )
        _relocate_archived_pair_evidence(
            archive, original_root=pair_dir,
        )
        return
    if pair_dir.exists():
        pair_dir.rename(archive)
        _relocate_archived_pair_evidence(
            archive, original_root=pair_dir,
        )


def execute_cell(run_dir: Path, cell_id: str) -> int:
    plan = load_plan(run_dir)
    cell = _cell_for(plan, cell_id)
    cell_dir = Path(cell["cell_dir"])
    success_path = cell_dir / ".success"
    if success_path.exists():
        success = read_json(success_path)
        if success.get("fingerprint") == cell["fingerprint"]:
            return 0

    try:
        assert_matching_provenance(plan)
    except Exception as exc:
        _mark_failed(cell, [str(exc)], returncode=None, attempt=_read_status(cell).get("attempts", 0))
        return 2

    old_status = _read_status(cell)
    attempt = int(old_status.get("attempts", 0)) + 1
    try:
        _prepare_attempt_workspace(cell, attempt)
    except (OSError, RuntimeError) as exc:
        _mark_failed(
            cell,
            [f"could not archive the prior paired attempt: {exc}"],
            returncode=None,
            attempt=attempt,
        )
        return 2
    started_ns = time.time_ns()
    started_at = utc_now()
    _unlink_markers(cell_dir)
    worker_record = _read_proc_stat(os.getpid())
    _write_status(
        cell, "running", attempts=attempt, started_at=started_at,
        started_ns=started_ns, session=cell_session_name(plan, cell),
        worker_pid=os.getpid(),
        worker_start_ticks=(worker_record or {}).get("start_ticks"),
        isolated_retry=bool(old_status.get("isolated_retry")),
    )
    environment = os.environ.copy()
    environment.update({
        "PYTHONDONTWRITEBYTECODE": "1",
        "PYTHONNOUSERSITE": "1",
        "PYTHONHASHSEED": "0",
        "LIFTON_BENCH_THREADS": str(cell["threads"]),
    })
    environment.update({str(key): str(value) for key, value in cell.get("environment", {}).items()})
    stdout_path = cell_dir / f"attempt-{attempt:02d}.stdout.log"
    stderr_path = cell_dir / f"attempt-{attempt:02d}.stderr.log"
    outcome = _run_supervised(
        plan=plan, cell=cell, attempt=attempt, environment=environment,
        stdout_path=stdout_path, stderr_path=stderr_path,
    )
    returncode = int(outcome["returncode"])
    elapsed = float(outcome["elapsed_seconds"])
    launch_error = outcome.get("launch_error")
    exit_document = {
        "schema_version": CONTROLLER_SCHEMA_VERSION,
        "cell_id": cell["id"],
        "fingerprint": cell["fingerprint"],
        "attempt": attempt,
        "command": cell["command"],
        "started_at": started_at,
        "started_ns": started_ns,
        "finished_at": utc_now(),
        "elapsed_seconds": elapsed,
        "returncode": returncode,
        "launch_error": launch_error,
        "stdout": str(stdout_path),
        "stderr": str(stderr_path),
        "watchdog": outcome["watchdog"],
        "watchdog_path": outcome["watchdog_path"],
    }
    atomic_write_json(cell_dir / f"attempt-{attempt:02d}.exit.json", exit_document)
    atomic_write_json(cell_dir / "exit.json", exit_document)

    errors = []
    validation: dict[str, Any] = {}
    watchdog_reason = outcome["watchdog"].get("reason")
    if watchdog_reason:
        errors.append(
            f"watchdog stopped the command ({watchdog_reason}) after {elapsed:.1f}s"
        )
    if launch_error:
        errors.append(f"could not launch command: {launch_error}")
    if returncode != 0:
        errors.append(f"command exited with status {returncode}")
    else:
        try:
            assert_matching_provenance(plan)
        except Exception as exc:
            errors.append(f"provenance changed while the cell was running: {exc}")
        if not errors:
            errors, validation = validate_artifacts(cell, started_ns)
    regressions: list[dict[str, Any]] = []
    if not errors:
        regressions = performance_regressions(cell, Path(plan["inputs"]["baseline"]))
        performance_retry_path = cell_dir / "performance_retry.json"
        if regressions and not performance_retry_path.exists():
            retry = {
                "schema_version": CONTROLLER_SCHEMA_VERSION,
                "cell_id": cell["id"],
                "fingerprint": cell["fingerprint"],
                "attempt": attempt,
                "reason": "performance regression above 25%; one isolated rerun required",
                "regressions": regressions,
                "created_at": utc_now(),
            }
            atomic_write_json(performance_retry_path, retry)
            _write_status(
                cell, "retry_pending", attempts=attempt, isolated_retry=True,
                validation=validation, performance=regressions,
            )
            return 75
        if regressions:
            errors.append("performance regression remained above 25% on isolated rerun")

    if errors:
        _mark_failed(
            cell, errors, returncode=returncode, attempt=attempt,
            validation=validation, performance=regressions,
            watchdog=outcome["watchdog"],
        )
        return 1

    success_document = {
        "schema_version": CONTROLLER_SCHEMA_VERSION,
        "cell_id": cell["id"],
        "fingerprint": cell["fingerprint"],
        "attempt": attempt,
        "completed_at": utc_now(),
        "exit": exit_document,
        "validation": validation,
        "performance": regressions,
    }
    _mark_success(
        cell,
        success_document,
        validation=validation,
        performance=regressions,
    )
    return 0


def _mark_success(
    cell: Mapping[str, Any],
    document: Mapping[str, Any],
    *,
    validation: Mapping[str, Any],
    performance: Sequence[Mapping[str, Any]],
) -> None:
    cell_dir = Path(cell["cell_dir"])
    with _terminal_transition(cell):
        atomic_write_json(cell_dir / ".success", document)
        try:
            (cell_dir / ".failed.json").unlink()
        except FileNotFoundError:
            pass
        _write_status(
            cell,
            "success",
            clear_fields=("errors", "failed_at", "returncode", "watchdog"),
            attempts=int(document["attempt"]),
            completed_at=document["completed_at"],
            isolated_retry=False,
            validation=validation,
            performance=performance,
        )


def _success_evidence_errors(
    cell: Mapping[str, Any],
    success: Mapping[str, Any],
) -> list[str]:
    validation = success.get("validation")
    if not isinstance(validation, Mapping):
        return ["success marker validation evidence is missing"]
    evidence_records: dict[str, Any] = {}
    artifacts = validation.get("artifacts")
    if isinstance(artifacts, Mapping):
        evidence_records.update(artifacts)
    evaluation_artifacts = validation.get("evaluation_artifacts")
    if isinstance(evaluation_artifacts, Mapping):
        for label, group in evaluation_artifacts.items():
            if not isinstance(group, Mapping):
                continue
            for name, metadata in group.items():
                evidence_records[f"{label}_{name}"] = metadata
    errors = []
    for name, metadata in evidence_records.items():
        if not isinstance(metadata, Mapping):
            errors.append(f"validated artifact evidence is malformed: {name}")
            continue
        path_value = metadata.get("path")
        if not isinstance(path_value, str):
            errors.append(f"validated artifact path is malformed: {name}")
            continue
        path = Path(path_value)
        try:
            stat = path.stat()
        except OSError:
            errors.append(f"validated artifact disappeared: {name}={path}")
            continue
        if (
            stat.st_size != metadata.get("size")
            or stat.st_mtime_ns != metadata.get("mtime_ns")
        ):
            errors.append(f"validated artifact changed after success: {name}={path}")
            continue
        recorded_sha = metadata.get("sha256")
        if not _sha256_text(recorded_sha):
            errors.append(f"validated artifact lacks SHA-256 evidence: {name}={path}")
        elif sha256_file(path) != recorded_sha:
            errors.append(f"validated artifact hash changed after success: {name}={path}")
    return errors


def _recorded_plan_file_errors(plan: Mapping[str, Any]) -> list[str]:
    provenance = plan.get("provenance")
    files = provenance.get("files") if isinstance(provenance, Mapping) else None
    if not isinstance(files, Mapping):
        return []
    errors = []
    for label, record in files.items():
        if not isinstance(record, Mapping):
            errors.append(f"recorded provenance for {label} is malformed")
            continue
        path_value = record.get("path")
        if not isinstance(path_value, str):
            errors.append(f"recorded provenance path for {label} is malformed")
            continue
        path = Path(path_value)
        try:
            stat = path.stat()
        except OSError as exc:
            errors.append(f"recorded provenance file unavailable: {label}: {exc}")
            continue
        if stat.st_size != record.get("size"):
            errors.append(f"recorded provenance size changed: {label}={path}")
        recorded_sha = record.get("sha256")
        if not _sha256_text(recorded_sha):
            errors.append(f"recorded provenance SHA-256 is malformed: {label}")
        elif sha256_file(path) != recorded_sha:
            errors.append(f"recorded provenance hash changed: {label}={path}")
    return errors


def recover_terminal_race(
    run_dir: Path,
    requested: Sequence[str],
    audit_path: Path,
    *,
    apply: bool = False,
) -> dict[str, Any]:
    """Archive an exact false-orphan race and make its cells retryable.

    This administrative path intentionally does not compare the active source
    checkout with the frozen worker snapshot. It verifies the source and input
    hashes recorded by the immutable plan instead, because recovery is run
    from the maintained checkout while the campaign runs from its snapshot.
    """

    plan = load_plan(run_dir)
    cell_ids = list(dict.fromkeys(requested))
    if not cell_ids:
        raise ValueError("at least one explicit cell is required")
    unknown = set(cell_ids) - {cell["id"] for cell in plan["cells"]}
    if unknown:
        raise ValueError(f"unknown cell ids: {sorted(unknown)}")
    audit_path = audit_path.resolve()
    if not audit_path.is_relative_to(run_dir.resolve()):
        raise ValueError("audit record must be below the run directory")
    if audit_path.exists():
        raise FileExistsError(f"audit record already exists: {audit_path}")

    expected_errors = [
        "worker tmux session disappeared before a final status was written",
        "orphan cleanup: process identity no longer matches; signal refused",
    ]
    cells = [_cell_for(plan, cell_id) for cell_id in cell_ids]
    preflight = []
    errors = _recorded_plan_file_errors(plan)
    for cell in cells:
        cell_dir = Path(cell["cell_dir"])
        status = _read_status(cell)
        success_path = cell_dir / ".success"
        failed_path = cell_dir / ".failed.json"
        try:
            success = read_json(success_path)
            failed = read_json(failed_path)
        except (OSError, TypeError, ValueError) as exc:
            errors.append(f"{cell['id']}: terminal marker is unreadable: {exc}")
            continue
        cell_errors = []
        if status.get("state") != "success":
            cell_errors.append(f"status state is {status.get('state')!r}")
        if status.get("attempts") != 1 or success.get("attempt") != 1:
            cell_errors.append("terminal race recovery requires attempt 1")
        if success.get("fingerprint") != cell["fingerprint"]:
            cell_errors.append("success marker fingerprint does not match plan")
        if failed.get("fingerprint") != cell["fingerprint"]:
            cell_errors.append("failed marker fingerprint does not match plan")
        if status.get("errors") != expected_errors:
            cell_errors.append("status errors do not match the known orphan race")
        if failed.get("errors") != expected_errors:
            cell_errors.append("failed marker errors do not match the known orphan race")
        if failed.get("watchdog", {}).get("reason") != "orphaned_worker":
            cell_errors.append("failed marker is not an orphaned-worker record")
        if success.get("exit", {}).get("returncode") != 0:
            cell_errors.append("success marker does not record return code zero")
        if not (
            isinstance(failed.get("failed_at"), str)
            and isinstance(success.get("completed_at"), str)
            and failed["failed_at"] < success["completed_at"]
        ):
            cell_errors.append("orphan failure did not precede successful completion")
        if tmux_has_session(cell_session_name(plan, cell)):
            cell_errors.append("cell tmux session is still active")
        if _worker_identity_matches(status):
            cell_errors.append("cell worker process is still active")
        cell_errors.extend(_success_evidence_errors(cell, success))
        if cell_errors:
            errors.extend(f"{cell['id']}: {error}" for error in cell_errors)
        preflight.append({
            "cell_id": cell["id"],
            "status": status,
            "success_sha256": sha256_file(success_path),
            "failed_sha256": sha256_file(failed_path),
            "success_path": str(success_path),
            "failed_path": str(failed_path),
        })
    if errors:
        raise RuntimeError("terminal race recovery preflight failed: " + "; ".join(errors))

    recovery_id = "terminal-race-" + dt.datetime.now(dt.timezone.utc).strftime(
        "%Y%m%dT%H%M%SZ"
    )
    audit = {
        "schema_version": CONTROLLER_SCHEMA_VERSION,
        "kind": "administrative_terminal_race_recovery",
        "recovery_id": recovery_id,
        "run_id": plan["run_id"],
        "run_dir": str(run_dir.resolve()),
        "plan_fingerprint": plan["fingerprint"],
        "cells": [],
        "applied": apply,
        "created_at": utc_now(),
    }
    if not apply:
        audit["preflight"] = preflight
        return audit

    moved: list[tuple[Path, Path]] = []
    try:
        for item, cell in zip(preflight, cells):
            cell_dir = Path(cell["cell_dir"])
            archive_dir = cell_dir / "administrative-recovery" / recovery_id
            archive_dir.mkdir(parents=True, exist_ok=False)
            archived = {}
            for name in ("status.json", ".success", ".failed.json"):
                source = cell_dir / name
                destination = archive_dir / name
                os.replace(source, destination)
                moved.append((destination, source))
                archived[name] = str(destination)
            _write_status(
                cell,
                "failed",
                attempts=1,
                failed_at=utc_now(),
                errors=["administrative rerun requested after terminal race"],
                returncode=None,
                isolated_retry=False,
                validation={},
                performance=[],
                watchdog={"reason": "administrative_terminal_race_recovery"},
            )
            audit["cells"].append({
                "cell_id": cell["id"],
                "before": item,
                "archive": archived,
                "state_after": "failed",
            })
        atomic_write_json(audit_path, audit)
    except Exception:
        try:
            for destination, source in reversed(moved):
                if destination.exists():
                    os.replace(destination, source)
        finally:
            raise
    return audit


def _mark_failed(cell: Mapping[str, Any], errors: Sequence[str], *,
                 returncode: int | None, attempt: int,
                 validation: Mapping[str, Any] | None = None,
                 performance: Sequence[Mapping[str, Any]] | None = None,
                 watchdog: Mapping[str, Any] | None = None,
                 expected_states: Sequence[str] | None = None) -> bool:
    document = {
        "schema_version": CONTROLLER_SCHEMA_VERSION,
        "cell_id": cell["id"],
        "fingerprint": cell["fingerprint"],
        "attempt": attempt,
        "failed_at": utc_now(),
        "returncode": returncode,
        "errors": list(errors),
        "validation": dict(validation or {}),
        "performance": list(performance or []),
        "watchdog": dict(watchdog or {}),
    }
    with _terminal_transition(cell):
        if (
            expected_states is not None
            and _read_status(cell).get("state") not in set(expected_states)
        ):
            return False
        atomic_write_json(Path(cell["cell_dir"]) / ".failed.json", document)
        _write_status(
            cell, "failed", attempts=attempt, failed_at=document["failed_at"],
            errors=list(errors), returncode=returncode, isolated_retry=False,
            validation=document["validation"], performance=document["performance"],
            watchdog=document["watchdog"],
        )
    return True


def mem_available_gib(meminfo: Path = Path("/proc/meminfo")) -> float:
    text = Path(meminfo).read_text(encoding="utf-8")
    match = re.search(r"^MemAvailable:\s+(\d+)\s+kB$", text, re.MULTILINE)
    if not match:
        raise RuntimeError(f"MemAvailable is absent from {meminfo}")
    return int(match.group(1)) / (1024.0 * 1024.0)


def resource_snapshot() -> dict[str, float]:
    return {"load1": os.getloadavg()[0], "available_gib": mem_available_gib()}


def launch_allowed(cell: Mapping[str, Any], active: Sequence[Mapping[str, Any]],
                   resources: Mapping[str, float], policy: Policy) -> tuple[bool, str]:
    if resources["load1"] >= policy.load1_limit:
        return False, f"load1 {resources['load1']:.2f} >= {policy.load1_limit:.2f}"
    if resources["available_gib"] < policy.min_available_gib:
        return False, (
            f"MemAvailable {resources['available_gib']:.1f} GiB "
            f"< {policy.min_available_gib:.1f} GiB"
        )
    if any(item.get("exclusive") for item in active):
        return False, "an exclusive paired cell is active"
    if cell.get("exclusive") and active:
        return False, "paired cell must run exclusively"
    if len(active) >= policy.max_active:
        return False, "active-cell limit reached"
    active_cost = sum(
        int(item.get("scheduler_thread_cost", item["threads"]))
        for item in active
    )
    cell_cost = int(cell.get("scheduler_thread_cost", cell["threads"]))
    if active_cost + cell_cost > policy.max_worker_threads:
        return False, "worker-thread limit reached"
    if cell.get("full_job") and sum(bool(item.get("full_job")) for item in active) >= policy.max_full:
        return False, "full-job limit reached"
    if any(item.get("isolated_retry") for item in active):
        return False, "an isolated performance retry is active"
    if cell.get("isolated_retry") and active:
        return False, "performance retry must run in isolation"
    return True, "ok"


def controller_session_name(plan: Mapping[str, Any]) -> str:
    return safe_name(f"lifton-{plan['run_id']}-ctl", limit=75)


def cell_session_name(plan: Mapping[str, Any], cell: Mapping[str, Any]) -> str:
    suffix = cell["fingerprint"][:8]
    return safe_name(f"lifton-{plan['run_id']}-{cell['id']}-{suffix}", limit=75)


def tmux_has_session(name: str) -> bool:
    return subprocess.run(
        ["tmux", "has-session", "-t", f"={name}"],
        stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=False,
    ).returncode == 0


def _tmux_new_session(name: str, command: Sequence[str]) -> None:
    result = subprocess.run(
        ["tmux", "new-session", "-d", "-s", name, shlex.join(command)],
        text=True, capture_output=True, check=False,
    )
    if result.returncode != 0:
        raise RuntimeError(f"tmux could not create {name!r}: {result.stderr.strip()}")


def _python_worker_command(*arguments: str) -> list[str]:
    """Launch repository workers without importing bytecode into snapshots."""

    return [
        "/usr/bin/env",
        "PYTHONDONTWRITEBYTECODE=1",
        "PYTHONNOUSERSITE=1",
        sys.executable,
        "-B",
        *arguments,
    ]


class RunLogger:
    def __init__(self, path: Path):
        self.path = path

    def __call__(self, message: str) -> None:
        line = f"{utc_now()} {message}"
        print(line, flush=True)
        self.path.parent.mkdir(parents=True, exist_ok=True)
        with self.path.open("a", encoding="utf-8") as handle:
            handle.write(line + "\n")
            handle.flush()


class ControllerIsolationError(RuntimeError):
    pass


def _foreign_active_cells(plan: Mapping[str, Any]) -> list[str]:
    """Find live cell tmux sessions belonging to another run in this root."""

    current = Path(plan["run_dir"]).resolve()
    runs_root = current.parent
    now_ns = time.time_ns()
    active = []
    for plan_path in runs_root.glob("*/plan.json"):
        if plan_path.parent.resolve() == current:
            continue
        try:
            other = load_plan(plan_path.parent)
        except (OSError, TypeError, ValueError):
            continue
        for cell in other["cells"]:
            status = _read_status(cell)
            if status.get("state") not in ACTIVE_STATES:
                continue
            session_alive = tmux_has_session(cell_session_name(other, cell))
            process_alive = False
            try:
                process_alive = _pid_identity_matches(
                    int(status["command_pid"]),
                    int(status["command_pgid"]),
                    int(status["command_start_ticks"]),
                )
            except (KeyError, TypeError, ValueError):
                pass
            started_ns = int(status.get("started_ns", 0))
            launch_is_recent = bool(
                started_ns and now_ns - started_ns < int(5.0 * 1e9)
            )
            worker_alive = _worker_identity_matches(status)
            if session_alive or process_alive or worker_alive or launch_is_recent:
                active.append(f"{other['run_id']}:{cell['id']}")
    return active


@contextlib.contextmanager
def _controller_lease(plan: Mapping[str, Any]):
    """Serialize schedulers and reject stale foreign workers in one runs root."""

    if fcntl is None:
        raise ControllerIsolationError("controller isolation requires POSIX fcntl.flock")
    lock_path = Path(plan["run_dir"]).resolve().parent / ".controller.lock"
    lock_path.parent.mkdir(parents=True, exist_ok=True)
    try:
        handle = lock_path.open("a+", encoding="utf-8")
    except OSError as exc:
        raise ControllerIsolationError(
            f"could not open benchmark controller lock {lock_path}: {exc}"
        ) from exc
    acquired = False
    try:
        try:
            fcntl.flock(handle.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
            acquired = True
        except BlockingIOError as exc:
            handle.seek(0)
            owner = handle.read().strip() or "unknown owner"
            raise ControllerIsolationError(
                f"another benchmark controller holds {lock_path}: {owner}"
            ) from exc
        except OSError as exc:
            raise ControllerIsolationError(
                f"could not acquire benchmark controller lock {lock_path}: {exc}"
            ) from exc

        foreign = _foreign_active_cells(plan)
        if foreign:
            raise ControllerIsolationError(
                "foreign benchmark cells are still active: " + ", ".join(foreign)
            )
        owner = {
            "schema_version": CONTROLLER_SCHEMA_VERSION,
            "state": "running",
            "run_id": plan["run_id"],
            "run_dir": plan["run_dir"],
            "pid": os.getpid(),
            "session": controller_session_name(plan),
            "acquired_at": utc_now(),
        }
        handle.seek(0)
        handle.truncate()
        json.dump(owner, handle, sort_keys=True)
        handle.write("\n")
        handle.flush()
        os.fsync(handle.fileno())
        try:
            yield
        finally:
            owner.update({"state": "released", "released_at": utc_now()})
            handle.seek(0)
            handle.truncate()
            json.dump(owner, handle, sort_keys=True)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
    finally:
        if acquired:
            fcntl.flock(handle.fileno(), fcntl.LOCK_UN)
        handle.close()


def _within_launch_grace(
    status: Mapping[str, Any], grace_seconds: float = LAUNCH_GRACE_SECONDS,
) -> bool:
    """True while a just-launched worker may not be observable yet.

    ``tmux new-session -d`` can return before the session is listable, and the
    worker's identity is only recorded once it starts. During that window a
    cell is genuinely consuming the host even though neither liveness probe
    can see it, so admission control and orphan detection must both treat it
    as active. Counting it as neither is what allowed seven whole-genome cells
    to run against ``max_active: 2``.
    """

    try:
        started_ns = int(status.get("started_ns") or 0)
    except (TypeError, ValueError):
        return False
    if not started_ns:
        return False
    return time.time_ns() - started_ns < int(grace_seconds * 1e9)


def _cell_is_live(
    plan: Mapping[str, Any],
    cell: Mapping[str, Any],
    status: Mapping[str, Any],
) -> bool:
    return (
        tmux_has_session(cell_session_name(plan, cell))
        or _worker_identity_matches(status)
        or _within_launch_grace(status)
    )


def _active_cells(plan: Mapping[str, Any]) -> list[dict[str, Any]]:
    active = []
    for cell in plan["cells"]:
        status = _read_status(cell)
        if status.get("state") in ACTIVE_STATES and _cell_is_live(
            plan, cell, status,
        ):
            active.append({**cell, "isolated_retry": bool(status.get("isolated_retry"))})
    return active


def _mark_orphans(
    plan: Mapping[str, Any], *, grace_seconds: float = LAUNCH_GRACE_SECONDS,
) -> None:
    policy = Policy(**plan["policy"])
    now_ns = time.time_ns()
    for cell in plan["cells"]:
        status = _read_status(cell)
        if status.get("state") not in ACTIVE_STATES:
            continue
        if (
            tmux_has_session(cell_session_name(plan, cell))
            or _worker_identity_matches(status)
        ):
            continue
        updated_ns = int(status.get("started_ns", 0))
        if updated_ns and now_ns - updated_ns < int(grace_seconds * 1e9):
            continue
        cleanup: dict[str, Any] = {}
        errors = ["worker tmux session disappeared before a final status was written"]
        try:
            pid = int(status["command_pid"])
            pgid = int(status["command_pgid"])
            start_ticks = int(status["command_start_ticks"])
        except (KeyError, TypeError, ValueError):
            pass
        else:
            cleanup = _terminate_process_group(
                pid=pid, pgid=pgid, start_ticks=start_ticks,
                grace_seconds=policy.terminate_grace_seconds,
            )
            if cleanup.get("error"):
                errors.append(f"orphan cleanup: {cleanup['error']}")
        _mark_failed(
            cell, errors,
            returncode=None, attempt=int(status.get("attempts", 0)),
            watchdog={"reason": "orphaned_worker", "cleanup": cleanup},
            expected_states=ACTIVE_STATES,
        )


def _launch_cell(plan: Mapping[str, Any], cell: Mapping[str, Any]) -> None:
    status = _read_status(cell)
    session = cell_session_name(plan, cell)
    _write_status(
        cell, "launching", attempts=int(status.get("attempts", 0)),
        session=session, launch_requested_at=utc_now(),
        isolated_retry=bool(status.get("isolated_retry")),
        started_ns=time.time_ns(),
        command_pid=None, command_pgid=None, command_start_ticks=None,
    )
    command = _python_worker_command(
        str(Path(__file__).resolve()), "worker",
        "--run-dir", plan["run_dir"], "--cell", cell["id"],
    )
    try:
        _tmux_new_session(session, command)
    except Exception as exc:
        _mark_failed(
            cell, [str(exc)], returncode=None,
            attempt=int(status.get("attempts", 0)),
        )
        raise


def scheduler_loop(run_dir: Path, *, sleep=time.sleep) -> int:
    plan = load_plan(run_dir)
    log = RunLogger(Path(run_dir) / "controller.log")
    try:
        with _controller_lease(plan):
            return _scheduler_loop_locked(run_dir, plan=plan, log=log, sleep=sleep)
    except ControllerIsolationError as exc:
        controller_path = Path(run_dir) / "controller.json"
        try:
            current = read_json(controller_path)
        except (OSError, TypeError, ValueError):
            current = {}
        if current.get("state") != "running":
            atomic_write_json(controller_path, {
                "schema_version": CONTROLLER_SCHEMA_VERSION,
                "run_id": plan["run_id"],
                "fingerprint": plan["fingerprint"],
                "state": "blocked",
                "reason": str(exc),
                "updated_at": utc_now(),
            })
        log(f"BLOCKED: {exc}")
        return 2


def _scheduler_loop_locked(
    run_dir: Path,
    *,
    plan: Mapping[str, Any],
    log: RunLogger,
    sleep=time.sleep,
) -> int:
    policy = Policy(**plan["policy"])
    try:
        assert_matching_provenance(plan)
    except Exception as exc:
        atomic_write_json(Path(run_dir) / "controller.json", {
            "schema_version": CONTROLLER_SCHEMA_VERSION,
            "run_id": plan["run_id"], "state": "blocked",
            "reason": str(exc), "updated_at": utc_now(),
        })
        log(f"BLOCKED: {exc}")
        return 2
    atomic_write_json(Path(run_dir) / "controller.json", {
        "schema_version": CONTROLLER_SCHEMA_VERSION,
        "run_id": plan["run_id"], "fingerprint": plan["fingerprint"],
        "state": "running", "session": controller_session_name(plan),
        "started_at": utc_now(),
    })
    log(f"scheduler start: {len(plan['cells'])} cells policy={plan['policy']}")
    while True:
        _mark_orphans(plan)
        active = _active_cells(plan)
        pending = []
        for cell in plan["cells"]:
            status = _read_status(cell)
            if status.get("state") in PENDING_STATES:
                pending.append({**cell, "isolated_retry": bool(status.get("isolated_retry"))})
        if not active and not pending:
            break

        resources = resource_snapshot()
        launched = 0
        isolated = [cell for cell in pending if cell.get("isolated_retry")]
        candidates = isolated[:1] if isolated else pending
        for cell in candidates:
            allowed, reason = launch_allowed(cell, active, resources, policy)
            if not allowed:
                continue
            _launch_cell(plan, cell)
            active.append(cell)
            launched += 1
            log(
                f"launch {cell['id']} active={len(active)}/{policy.max_active} "
                f"full={sum(bool(c.get('full_job')) for c in active)}/{policy.max_full} "
                f"load1={resources['load1']:.2f} available={resources['available_gib']:.1f}GiB"
            )
            if isolated:
                break
            sleep(policy.stagger_seconds)
            resources = resource_snapshot()
        if not launched:
            counts = Counter(_read_status(cell).get("state", "unknown") for cell in plan["cells"])
            log(
                f"wait states={dict(counts)} load1={resources['load1']:.2f} "
                f"available={resources['available_gib']:.1f}GiB"
            )
        sleep(policy.poll_seconds)

    summary = summarize_run(plan)
    state = "success" if summary["counts"].get("failed", 0) == 0 else "failed"
    summary["state"] = state
    summary["finished_at"] = utc_now()
    atomic_write_json(Path(run_dir) / "summary.json", summary)
    atomic_write_json(Path(run_dir) / "controller.json", {
        "schema_version": CONTROLLER_SCHEMA_VERSION,
        "run_id": plan["run_id"], "fingerprint": plan["fingerprint"],
        "state": state, "finished_at": summary["finished_at"],
        "counts": summary["counts"],
    })
    log(f"scheduler exit state={state} counts={summary['counts']}")
    return 0 if state == "success" else 1


def summarize_run(plan: Mapping[str, Any]) -> dict[str, Any]:
    rows = []
    for cell in plan["cells"]:
        status = _read_status(cell)
        rows.append({
            "cell_id": cell["id"], "benchmark": cell["benchmark"],
            "kind": cell["kind"], "state": status.get("state", "unknown"),
            "attempts": status.get("attempts", 0),
            "errors": status.get("errors", []),
            "session": cell_session_name(plan, cell),
        })
    counts = dict(Counter(row["state"] for row in rows))
    try:
        controller_state = read_json(
            Path(plan["run_dir"]) / "controller.json"
        ).get("state", "unknown")
    except (OSError, TypeError, ValueError):
        controller_state = "unknown"
    return {
        "schema_version": CONTROLLER_SCHEMA_VERSION,
        "run_id": plan["run_id"], "fingerprint": plan["fingerprint"],
        "stage": plan["stage"], "counts": counts, "cells": rows,
        "controller_state": controller_state,
        "updated_at": utc_now(),
    }


def status_exit_code(summary: Mapping[str, Any]) -> int:
    """Return a stable code: 0 success, 1 failed, 2 invalid, 3 incomplete."""

    if summary.get("provenance_error"):
        return STATUS_EXIT_INVALID
    counts = summary.get("counts")
    if not isinstance(counts, Mapping):
        return STATUS_EXIT_INVALID
    if (
        any(int(counts.get(state, 0)) > 0 for state in ("failed", "blocked"))
        or summary.get("controller_state") in {"failed", "blocked"}
    ):
        return STATUS_EXIT_FAILED
    total = sum(
        int(count) for count in counts.values()
        if isinstance(count, int) and not isinstance(count, bool)
    )
    if (
        total > 0
        and int(counts.get("success", 0)) == total
        and summary.get("controller_state") == "success"
    ):
        return STATUS_EXIT_SUCCESS
    return STATUS_EXIT_INCOMPLETE


def start_controller_tmux(plan: Mapping[str, Any]) -> str:
    session = controller_session_name(plan)
    if tmux_has_session(session):
        return session
    command = _python_worker_command(
        str(Path(__file__).resolve()), "worker",
        "--run-dir", plan["run_dir"],
    )
    _tmux_new_session(session, command)
    return session


def retry_failed(run_dir: Path, requested: Sequence[str] | None = None) -> list[str]:
    plan = load_plan(run_dir)
    assert_matching_provenance(plan)
    wanted = set(requested or [])
    unknown = wanted - {cell["id"] for cell in plan["cells"]}
    if unknown:
        raise ValueError(f"unknown cell ids: {sorted(unknown)}")
    reset = []
    for cell in plan["cells"]:
        status = _read_status(cell)
        if status.get("state") != "failed" or (wanted and cell["id"] not in wanted):
            continue
        cell_dir = Path(cell["cell_dir"])
        try:
            (cell_dir / ".failed.json").unlink()
        except FileNotFoundError:
            pass
        _write_status(
            cell, "pending", attempts=int(status.get("attempts", 0)),
            manual_retries=int(status.get("manual_retries", 0)) + 1,
            errors=[], isolated_retry=False,
        )
        reset.append(cell["id"])
    return reset


def retry_candidates(run_dir: Path, requested: Sequence[str] | None = None) -> list[str]:
    plan = load_plan(run_dir)
    assert_matching_provenance(plan)
    wanted = set(requested or [])
    unknown = wanted - {cell["id"] for cell in plan["cells"]}
    if unknown:
        raise ValueError(f"unknown cell ids: {sorted(unknown)}")
    return [
        cell["id"] for cell in plan["cells"]
        if _read_status(cell).get("state") == "failed"
        and (not wanted or cell["id"] in wanted)
    ]


def reconcile_run(run_dir: Path, *, deep: bool = False) -> dict[str, Any]:
    plan = load_plan(run_dir)
    provenance_error = None
    try:
        assert_matching_provenance(plan)
    except Exception as exc:
        provenance_error = str(exc)
    _mark_orphans(plan, grace_seconds=0)
    current_results: dict[str, Any] = {}
    for cell in plan["cells"]:
        status = _read_status(cell)
        if status.get("state") != "success":
            continue
        cell_dir = Path(cell["cell_dir"])
        try:
            success = read_json(cell_dir / ".success")
        except (OSError, ValueError, TypeError) as exc:
            _mark_failed(
                cell, [f"success marker is unreadable: {exc}"], returncode=None,
                attempt=int(status.get("attempts", 0)),
            )
            continue
        errors = []
        if success.get("fingerprint") != cell["fingerprint"]:
            errors.append("success marker fingerprint does not match the plan")
        recorded = success.get("validation", {}).get("artifacts", {})
        evidence_records = {
            str(name): metadata for name, metadata in recorded.items()
        } if isinstance(recorded, Mapping) else {}
        evaluation_records = success.get("validation", {}).get(
            "evaluation_artifacts", {}
        )
        if isinstance(evaluation_records, Mapping):
            for label, group in evaluation_records.items():
                if not isinstance(group, Mapping):
                    continue
                for name, metadata in group.items():
                    evidence_records[f"{label}_{name}"] = metadata
        for name, metadata in evidence_records.items():
            if not isinstance(metadata, Mapping):
                errors.append(f"validated artifact evidence is malformed: {name}")
                continue
            path = Path(metadata.get("path", ""))
            try:
                stat = path.stat()
            except OSError:
                errors.append(f"validated artifact disappeared: {name}={path}")
                continue
            if stat.st_size != metadata.get("size") or stat.st_mtime_ns != metadata.get("mtime_ns"):
                errors.append(f"validated artifact changed after success: {name}={path}")
                continue
            recorded_sha = metadata.get("sha256")
            if not _sha256_text(recorded_sha):
                errors.append(
                    f"validated artifact lacks SHA-256 evidence: {name}={path}"
                )
            elif sha256_file(path) != recorded_sha:
                errors.append(
                    f"validated artifact hash changed after success: {name}={path}"
                )
        if deep and not errors and cell["kind"] != "gate":
            started_ns = int(success.get("exit", {}).get("started_ns", 0))
            deep_errors, _ = validate_artifacts(
                cell,
                started_ns,
                persist_validator_reports=False,
            )
            errors.extend(deep_errors)
        if errors:
            _mark_failed(
                cell, errors, returncode=None, attempt=int(status.get("attempts", 0)),
            )
            continue
        if cell["kind"] != "gate":
            result = read_json(Path(cell["artifacts"]["result_json"]))
            key = cell["artifacts"].get("result_key", cell["id"])
            if cell["kind"] == "subset":
                current_results[key] = result[key]
            else:
                current_results[key] = result
    summary = summarize_run(plan)
    summary["provenance_error"] = provenance_error
    summary["deep_validation"] = deep
    atomic_write_json(Path(run_dir) / "reconciled_results.json", {
        "schema_version": CONTROLLER_SCHEMA_VERSION,
        "run_id": plan["run_id"],
        "note": "Review-only results; the canonical baseline is never modified.",
        "results": current_results,
    })
    atomic_write_json(Path(run_dir) / "summary.json", summary)
    return summary


def _run_refresh(benchmark: str, result_dir: Path, threads: int) -> int:
    """Read canonical cached alignments but isolate every newly written output."""

    from benchmarks.compare import devel_refresh

    result_dir.mkdir(parents=True, exist_ok=True)
    try:
        cached_inputs = devel_refresh._cached_inputs(benchmark)
    except Exception as exc:
        print(f"[{benchmark}] cached-input preflight failed: {exc}", file=sys.stderr)
        return 1
    devel_refresh._cached_inputs = lambda _benchmark: cached_inputs
    devel_refresh.WORK = result_dir / "work"
    devel_refresh.OUT_DIR = result_dir
    try:
        devel_refresh.run_devel_refresh(
            benchmark, threads=threads, threads_eval=threads, log=print,
        )
    except Exception as exc:
        print(f"[{benchmark}] refresh failed: {exc}", file=sys.stderr)
        return 1
    return 0


def _run_subset(benchmark: str, cell_dir: Path, threads: int) -> int:
    """Run canonical four-way subset scoring with isolated current outputs."""

    from benchmarks.compare import fourway_compare

    cell_dir.mkdir(parents=True, exist_ok=True)
    try:
        cached_inputs = fourway_compare._ensure_inputs(
            benchmark, threads=threads, log=print,
        )
        fourway_compare.vc.provenance_gate(["stable", "devel"], log=print)
    except Exception as exc:
        print(f"[{benchmark}] subset input/provenance preflight failed: {exc}", file=sys.stderr)
        return 1
    fourway_compare._ensure_inputs = lambda _benchmark, threads=8, log=print: cached_inputs
    fourway_compare.WORK = cell_dir / "work"
    fourway_compare.RESULTS = cell_dir / "result.json"
    try:
        fourway_compare.run_subset(benchmark, threads_eval=threads, log=print)
    except Exception as exc:
        print(f"[{benchmark}] subset comparison failed: {exc}", file=sys.stderr)
        return 1
    return 0


def _run_e2e(dataset: str, cell_dir: Path, registry_path: Path, threads: int) -> int:
    """Run Phase 16 with an isolated registry/result tree and strict exit status.

    ``datasets.json`` contains useful ``_note`` metadata alongside real dataset
    fields.  The legacy loader passes those keys directly to ``Dataset`` and
    raises ``TypeError``.  Sanitizing only underscore-prefixed metadata here
    avoids changing that historical harness during a benchmark build.
    """

    from benchmarks import run_benchmarks

    raw = read_json(registry_path)
    entries = raw.get("datasets", []) if isinstance(raw, dict) else []
    sanitized_entries = []
    for entry in entries:
        if not isinstance(entry, dict) or not entry.get("id"):
            continue
        sanitized_entries.append({
            key: value for key, value in entry.items() if not key.startswith("_")
        })
    if dataset not in {entry["id"] for entry in sanitized_entries}:
        print(f"unknown end-to-end dataset: {dataset}", file=sys.stderr)
        return 2
    flags = list(raw.get("lifton_flags", []))
    try:
        thread_index = flags.index("-t") + 1
        flags[thread_index] = str(threads)
    except (ValueError, IndexError):
        flags.extend(["-t", str(threads)])
    sanitized = {
        "_comment": "Generated by build_controller; source registry is fingerprinted in plan.json.",
        "datasets": sanitized_entries,
        "lifton_flags": flags,
        "evaluation_flags": list(raw.get("evaluation_flags", [])),
    }
    cell_dir.mkdir(parents=True, exist_ok=True)
    sanitized_registry = cell_dir / "datasets.sanitized.json"
    atomic_write_json(sanitized_registry, sanitized)
    registry = run_benchmarks.load_registry(sanitized_registry)
    selected = next(item for item in registry.datasets if item.id == dataset)
    data_dir = HERE.parent / "data"
    results_dir = cell_dir / "artifacts"
    dataset_dir = data_dir / dataset
    # LiftOn's candidate-side ``-E`` evaluation is self-contained.  A
    # published target annotation is an optional comparison input and must not
    # block an otherwise complete cached E2E cell.
    required_urls = [selected.reference_fa, selected.target_fa, selected.reference_gff]
    required_paths = [
        dataset_dir / run_benchmarks._filename_for(url) for url in required_urls
    ]
    missing = [path for path in required_paths if not path.is_file() or path.stat().st_size < 1024]
    if missing:
        print(
            "cached-only end-to-end preflight failed; downloads are disabled: "
            + ", ".join(str(path) for path in missing),
            file=sys.stderr,
        )
        return 1
    try:
        row = run_benchmarks.run_dataset(
            selected, data_dir=data_dir, results_dir=results_dir,
            lifton_flags=registry.lifton_flags,
            evaluation_flags=registry.evaluation_flags,
            do_download=False, do_lift=True, do_evaluate=True, force=True,
        )
    except Exception as exc:
        row = {"dataset": dataset, "species": selected.species, "error": str(exc)}
    payload = {"runtime": run_benchmarks._ensure_runtime(), "rows": [row]}
    atomic_write_json(cell_dir / "result.json", payload)
    if row.get("error"):
        return 1
    lift_profile = row.get("lift_profile") or {}
    eval_profile = row.get("eval_profile") or {}
    return 0 if lift_profile.get("exit_code") == 0 and eval_profile.get("exit_code") == 0 else 1


def _policy_from_args(args: argparse.Namespace) -> Policy:
    defaults = Policy()

    def selected(name: str) -> Any:
        value = getattr(args, name, None)
        return getattr(defaults, name) if value is None else value

    return Policy(
        threads_per_cell=selected("threads_per_cell"),
        scheduler_threads_per_cell=selected("scheduler_threads_per_cell"),
        max_active=selected("max_active"),
        max_full=selected("max_full"),
        max_worker_threads=selected("max_worker_threads"),
        load1_limit=selected("load1_limit"),
        min_available_gib=selected("min_available_gib"),
        stagger_seconds=selected("stagger_seconds"),
        poll_seconds=selected("poll_seconds"),
        subset_timeout_seconds=selected("subset_timeout_seconds"),
        full_timeout_seconds=selected("full_timeout_seconds"),
        e2e_timeout_seconds=selected("e2e_timeout_seconds"),
        stall_timeout_seconds=selected("stall_timeout_seconds"),
        watchdog_poll_seconds=selected("watchdog_poll_seconds"),
        terminate_grace_seconds=selected("terminate_grace_seconds"),
    )


def _run_dir_argument(args: argparse.Namespace) -> Path:
    if getattr(args, "run_dir", None):
        return Path(args.run_dir).resolve()
    return Path(args.runs_root).resolve() / args.run_id


def _print_summary(summary: Mapping[str, Any]) -> None:
    print(f"run: {summary['run_id']}  stage: {summary['stage']}")
    print("states: " + " ".join(
        f"{key}={value}" for key, value in sorted(summary["counts"].items())
    ))
    for row in summary["cells"]:
        detail = f" attempts={row['attempts']}" if row.get("attempts") else ""
        print(f"  {row['state']:<13} {row['cell_id']}{detail}")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Guarded tmux controller for LiftOn benchmark builds.",
    )
    subparsers = parser.add_subparsers(dest="action", required=True)

    start = subparsers.add_parser("start", help="create/resume a run and start its controller")
    start.add_argument("--run-id")
    start.add_argument(
        "--stage",
        choices=(
            "gates", "subset-canary", "subset", "full-canary", "full",
            "e2e-canary", "e2e", *PAIRED_STAGES,
        ),
        default=None,
    )
    start.add_argument("--ids", nargs="+")
    start.add_argument("--runs-root", default=str(DEFAULT_RUNS_ROOT))
    start.add_argument("--registry")
    start.add_argument("--dataset-registry")
    start.add_argument("--baseline")
    start.add_argument("--campaign-profile")
    start.add_argument("--campaign-id")
    start.add_argument("--profile-registry")
    start.add_argument("--candidate-root")
    start.add_argument("--candidate-sha")
    start.add_argument("--reference-root")
    start.add_argument("--reference-sha")
    start.add_argument(
        "--paired-repetitions",
        type=int,
        help=(
            "paired AB/BA repetitions (default: 4; non-canary stages require "
            "an even count; canaries allow 1-5)"
        ),
    )
    start.add_argument("--paired-lifton-executable")
    start.add_argument(
        "--candidate-e2e-mode", choices=PAIRED_E2E_MODES,
    )
    start.add_argument(
        "--reference-e2e-mode", choices=PAIRED_E2E_MODES,
    )
    start.add_argument("--dry-run", action="store_true")
    start.add_argument("--foreground", action="store_true")
    _add_policy_arguments(start)

    worker = subparsers.add_parser(
        "worker", help="run the scheduler or one cell (normally invoked by tmux)",
    )
    worker.add_argument("--run-dir", required=True)
    worker.add_argument("--cell")

    status = subparsers.add_parser("status", help="show immutable-plan cell status")
    status.add_argument("run_id", nargs="?")
    status.add_argument("--run-dir")
    status.add_argument("--runs-root", default=str(DEFAULT_RUNS_ROOT))
    status.add_argument("--json", action="store_true")

    retry = subparsers.add_parser("retry", help="reset failed cells and resume the controller")
    retry.add_argument("run_id", nargs="?")
    retry.add_argument("--run-dir")
    retry.add_argument("--runs-root", default=str(DEFAULT_RUNS_ROOT))
    retry.add_argument("--cells", nargs="+")
    retry.add_argument("--dry-run", action="store_true")
    retry.add_argument("--foreground", action="store_true")

    reconcile = subparsers.add_parser("reconcile", help="audit markers without changing baselines")
    reconcile.add_argument("run_id", nargs="?")
    reconcile.add_argument("--run-dir")
    reconcile.add_argument("--runs-root", default=str(DEFAULT_RUNS_ROOT))
    reconcile.add_argument("--deep", action="store_true")
    reconcile.add_argument("--json", action="store_true")

    recover = subparsers.add_parser(
        "recover-terminal-race",
        help="archive exact false-orphan records and make cells retryable",
    )
    recover.add_argument("--run-dir", required=True)
    recover.add_argument("--cells", nargs="+", required=True)
    recover.add_argument("--audit-record", required=True)
    recover.add_argument("--apply", action="store_true")

    refresh = subparsers.add_parser(
        "_run-refresh", help="internal isolated devel-refresh runner",
    )
    refresh.add_argument("--benchmark", required=True)
    refresh.add_argument("--result-dir", required=True)
    refresh.add_argument("--threads", type=int, default=8)
    subset = subparsers.add_parser(
        "_run-subset", help="internal isolated four-way subset runner",
    )
    subset.add_argument("--benchmark", required=True)
    subset.add_argument("--cell-dir", required=True)
    subset.add_argument("--threads", type=int, default=8)
    e2e = subparsers.add_parser(
        "_run-e2e", help="internal isolated Phase 16 runner",
    )
    e2e.add_argument("--dataset", required=True)
    e2e.add_argument("--cell-dir", required=True)
    e2e.add_argument("--registry", required=True)
    e2e.add_argument("--threads", type=int, default=8)
    return parser


def _add_policy_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--threads-per-cell", type=int)
    parser.add_argument("--scheduler-threads-per-cell", type=int)
    parser.add_argument("--max-active", type=int)
    parser.add_argument("--max-full", type=int)
    parser.add_argument("--max-worker-threads", type=int)
    parser.add_argument("--load1-limit", type=float)
    parser.add_argument("--min-available-gib", type=float)
    parser.add_argument("--stagger-seconds", type=float)
    parser.add_argument("--poll-seconds", type=float)
    parser.add_argument("--subset-timeout-seconds", type=float)
    parser.add_argument("--full-timeout-seconds", type=float)
    parser.add_argument("--e2e-timeout-seconds", type=float)
    parser.add_argument("--stall-timeout-seconds", type=float)
    parser.add_argument("--watchdog-poll-seconds", type=float)
    parser.add_argument("--terminate-grace-seconds", type=float)


def _campaign_case_from_args(
    args: argparse.Namespace,
) -> dict[str, Any] | None:
    supplied = {
        "--campaign-profile": getattr(args, "campaign_profile", None),
        "--campaign-id": getattr(args, "campaign_id", None),
    }
    if not any(supplied.values()):
        if getattr(args, "profile_registry", None) is not None:
            raise ValueError(
                "--profile-registry requires --campaign-profile and --campaign-id"
            )
        return None
    missing = [flag for flag, value in supplied.items() if value is None]
    if missing:
        raise ValueError(
            "profile selection requires both --campaign-profile and "
            f"--campaign-id; missing={missing}"
        )
    from benchmarks.compare import campaign_profiles

    registry = Path(
        args.profile_registry or DEFAULT_PROFILE_REGISTRY
    ).resolve()
    profile, case = campaign_profiles.select_profile_case(
        args.campaign_profile,
        args.campaign_id,
        registry=registry,
    )
    if case["stage"] in campaign_profiles.PROTOCOL_STAGES:
        raise ValueError(
            f"campaign {case['id']!r} is a single-source protocol schedule; "
            "run it through campaign_orchestrator"
        )
    conflicts = []
    for option, current, expected in (
        ("--stage", args.stage, case["stage"]),
        ("--ids", args.ids, case["ids"]),
        (
            "--paired-repetitions",
            args.paired_repetitions,
            case["repetitions"],
        ),
        ("--threads-per-cell", args.threads_per_cell, case["threads"]),
        (
            "--candidate-e2e-mode",
            args.candidate_e2e_mode,
            case["candidate_mode"],
        ),
        (
            "--reference-e2e-mode",
            args.reference_e2e_mode,
            case["reference_mode"],
        ),
    ):
        if current is not None and current != expected:
            conflicts.append(f"{option}={current!r} (profile={expected!r})")
    if conflicts:
        raise ValueError(
            "explicit options conflict with the selected campaign profile: "
            + "; ".join(conflicts)
        )
    args.stage = case["stage"]
    args.ids = list(case["ids"])
    args.threads_per_cell = case["threads"]
    if case["stage"].startswith(PAIRED_STAGE_PREFIX):
        args.paired_repetitions = case["repetitions"]
        args.candidate_e2e_mode = case["candidate_mode"]
        args.reference_e2e_mode = case["reference_mode"]
    return campaign_profiles.case_identity(profile, case)


def _paired_cli_supplied(args: argparse.Namespace) -> bool:
    return any(
        getattr(args, name, None) is not None
        for name in (
            "candidate_root", "candidate_sha",
            "reference_root", "reference_sha",
            "paired_repetitions", "paired_lifton_executable",
            "candidate_e2e_mode", "reference_e2e_mode",
        )
    )


def _paired_configuration_from_args(
    args: argparse.Namespace,
    stage: str,
    *,
    fallback: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    def selected(name: str, default: Any) -> Any:
        value = getattr(args, name, None)
        return default if value is None else value

    candidate_fallback = fallback.get("candidate", {}) if fallback else {}
    reference_fallback = fallback.get("reference", {}) if fallback else {}
    registry_fallback = fallback.get("registries", {}) if fallback else {}
    candidate_root = selected("candidate_root", candidate_fallback.get("root"))
    candidate_sha = selected("candidate_sha", candidate_fallback.get("sha"))
    reference_root = selected("reference_root", reference_fallback.get("root"))
    reference_sha = selected("reference_sha", reference_fallback.get("sha"))
    missing = [
        flag for flag, value in (
            ("--candidate-root", candidate_root),
            ("--candidate-sha", candidate_sha),
            ("--reference-root", reference_root),
            ("--reference-sha", reference_sha),
        )
        if value is None
    ]
    if missing:
        raise ValueError(
            "paired stages require exact source inputs: " + ", ".join(missing)
        )
    executable_default = (
        fallback.get("lifton_executable")
        if fallback else Path(sys.executable).with_name("lifton")
    )
    repetitions_default = (
        fallback.get("repetitions", DEFAULT_PAIRED_REPETITIONS)
        if fallback else DEFAULT_PAIRED_REPETITIONS
    )
    return paired_configuration(
        stage=stage,
        candidate_root=Path(candidate_root),
        candidate_sha=str(candidate_sha),
        reference_root=Path(reference_root),
        reference_sha=str(reference_sha),
        repetitions=int(selected("paired_repetitions", repetitions_default)),
        lifton_executable=Path(selected(
            "paired_lifton_executable", executable_default,
        )),
        candidate_e2e_mode=str(selected(
            "candidate_e2e_mode", candidate_fallback.get("e2e_mode", "fast"),
        )),
        reference_e2e_mode=str(selected(
            "reference_e2e_mode", reference_fallback.get("e2e_mode", "fast"),
        )),
        benchmark_registry=Path(selected(
            "registry", registry_fallback.get("benchmark", DEFAULT_REGISTRY),
        )),
        dataset_registry=Path(selected(
            "dataset_registry",
            registry_fallback.get("dataset", DEFAULT_DATASET_REGISTRY),
        )),
    )


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    if args.action == "start":
        try:
            campaign_case = _campaign_case_from_args(args)
        except (OSError, ValueError) as exc:
            parser.error(str(exc))
        runs_root = Path(args.runs_root).resolve()
        registry = Path(args.registry or DEFAULT_REGISTRY).resolve()
        dataset_registry = Path(
            args.dataset_registry or DEFAULT_DATASET_REGISTRY
        ).resolve()
        baseline = Path(args.baseline or DEFAULT_BASELINE).resolve()
        if args.run_id and (runs_root / safe_name(args.run_id, limit=100) / "plan.json").exists():
            run_dir = runs_root / safe_name(args.run_id, limit=100)
            plan = load_plan(run_dir)
            if (
                campaign_case is not None
                and campaign_case != plan.get("campaign_case")
            ):
                parser.error(
                    "existing run campaign profile/case does not match the "
                    "requested selection"
                )
            if (
                campaign_case is None
                and (
                    args.campaign_profile is not None
                    or args.campaign_id is not None
                )
            ):
                parser.error("existing run has no matching campaign profile")
            if args.stage and args.stage != plan["stage"]:
                parser.error(f"existing run stage is {plan['stage']!r}, not {args.stage!r}")
            if args.ids and list(args.ids) != plan["ids"]:
                parser.error("existing run ids do not match the requested ids")
            for option, supplied, requested_path, planned_path in (
                ("--registry", args.registry, registry, plan["inputs"]["registry"]),
                (
                    "--dataset-registry", args.dataset_registry,
                    dataset_registry, plan["inputs"]["dataset_registry"],
                ),
                ("--baseline", args.baseline, baseline, plan["inputs"]["baseline"]),
            ):
                if (
                    supplied is not None
                    and requested_path != Path(planned_path).resolve()
                ):
                    parser.error(
                        f"existing run {option} does not match the requested path"
                    )
            if plan.get("paired") is not None:
                if _paired_cli_supplied(args):
                    try:
                        requested_paired = _paired_configuration_from_args(
                            args, plan["stage"], fallback=plan["paired"],
                        )
                    except ValueError as exc:
                        parser.error(str(exc))
                    if requested_paired != plan["paired"]:
                        parser.error(
                            "existing paired run source/mode configuration "
                            "does not match the requested configuration"
                        )
            elif _paired_cli_supplied(args):
                parser.error(
                    "paired source options cannot be applied to an existing "
                    "non-paired run"
                )
            try:
                assert_matching_provenance(plan)
            except RuntimeError as exc:
                parser.error(str(exc))
        else:
            stage = args.stage or "subset-canary"
            try:
                paired = None
                if _paired_panel(stage) is not None:
                    paired = _paired_configuration_from_args(args, stage)
                elif _paired_cli_supplied(args):
                    raise ValueError(
                        "paired source options require a paired-* stage"
                    )
                run_dir, plan = create_plan(
                    run_id=args.run_id, stage=stage, requested_ids=args.ids,
                    runs_root=runs_root, repo_root=REPO_ROOT,
                    registry=registry,
                    dataset_registry=dataset_registry,
                    baseline=baseline, policy=_policy_from_args(args),
                    paired=paired,
                    campaign_case=campaign_case,
                    profile_registry=(
                        Path(args.profile_registry or DEFAULT_PROFILE_REGISTRY)
                        if campaign_case is not None else None
                    ),
                )
                initialize_run(run_dir, plan)
            except (OSError, RuntimeError, ValueError) as exc:
                parser.error(str(exc))
        print(f"run directory: {run_dir}")
        print(f"plan fingerprint: {plan['fingerprint']}")
        for cell in plan["cells"]:
            print(f"  {cell['id']}: {shlex.join(cell['command'])}")
        if args.dry_run:
            return 0
        if args.foreground:
            return scheduler_loop(run_dir)
        try:
            session = start_controller_tmux(plan)
        except (OSError, RuntimeError) as exc:
            parser.error(str(exc))
        print(f"controller tmux session: {session}")
        return 0

    if args.action == "worker":
        run_dir = Path(args.run_dir).resolve()
        return execute_cell(run_dir, args.cell) if args.cell else scheduler_loop(run_dir)

    if args.action == "recover-terminal-race":
        try:
            result = recover_terminal_race(
                Path(args.run_dir).resolve(),
                args.cells,
                Path(args.audit_record),
                apply=args.apply,
            )
        except (OSError, RuntimeError, ValueError) as exc:
            parser.error(str(exc))
        print(json.dumps(result, indent=2, sort_keys=True))
        return 0

    if args.action in {"status", "retry", "reconcile"}:
        if not args.run_dir and not args.run_id:
            parser.error("provide run_id or --run-dir")
        run_dir = _run_dir_argument(args)
        if args.action == "status":
            try:
                summary = summarize_run(load_plan(run_dir))
            except (OSError, RuntimeError, ValueError) as exc:
                parser.error(str(exc))
        elif args.action == "retry":
            try:
                reset = (
                    retry_candidates(run_dir, args.cells)
                    if args.dry_run else retry_failed(run_dir, args.cells)
                )
            except (OSError, RuntimeError, ValueError) as exc:
                parser.error(str(exc))
            verb = "would reset" if args.dry_run else "reset"
            print(f"{verb}: " + (", ".join(reset) if reset else "no failed cells"))
            if reset and not args.dry_run:
                if args.foreground:
                    return scheduler_loop(run_dir)
                try:
                    print(f"controller tmux session: {start_controller_tmux(load_plan(run_dir))}")
                except (OSError, RuntimeError) as exc:
                    parser.error(str(exc))
            return 0
        else:
            try:
                summary = reconcile_run(run_dir, deep=args.deep)
            except (OSError, RuntimeError, ValueError) as exc:
                parser.error(str(exc))
        if args.json:
            print(json.dumps(summary, indent=2, sort_keys=True))
        else:
            _print_summary(summary)
            if summary.get("provenance_error"):
                print(f"provenance: MISMATCH ({summary['provenance_error']})")
        return status_exit_code(summary)

    if args.action == "_run-refresh":
        return _run_refresh(args.benchmark, Path(args.result_dir), args.threads)
    if args.action == "_run-subset":
        return _run_subset(args.benchmark, Path(args.cell_dir), args.threads)
    if args.action == "_run-e2e":
        return _run_e2e(
            args.dataset, Path(args.cell_dir), Path(args.registry), args.threads,
        )
    parser.error(f"unsupported action: {args.action}")
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
