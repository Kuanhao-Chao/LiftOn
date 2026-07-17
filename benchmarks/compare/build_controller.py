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
import datetime as dt
import hashlib
import json
import os
import re
import shlex
import shutil
import subprocess
import sys
import time
from collections import Counter
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Mapping, Sequence


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
CONTROLLER_SCHEMA_VERSION = 1

SUBSET_CANARIES = ("human_mane", "human_to_zebrafish")
FULL_CANARIES = ("arabidopsis", "human_to_zebrafish")
E2E_DATASETS = ("bee", "human", "arabidopsis", "rice", "mouse")
E2E_CANARIES = ("bee", "human")

ACTIVE_STATES = {"launching", "running"}
PENDING_STATES = {"pending", "retry_pending"}


@dataclass(frozen=True)
class Policy:
    """Resource policy applied before every cell launch."""

    threads_per_cell: int = 8
    max_active: int = 4
    max_full: int = 2
    max_worker_threads: int = 32
    load1_limit: float = 32.0
    min_available_gib: float = 256.0
    stagger_seconds: float = 15.0
    poll_seconds: float = 30.0

    def validate(self) -> None:
        numeric_positive = {
            "threads_per_cell": self.threads_per_cell,
            "max_active": self.max_active,
            "max_full": self.max_full,
            "max_worker_threads": self.max_worker_threads,
            "load1_limit": self.load1_limit,
            "min_available_gib": self.min_available_gib,
        }
        bad = [name for name, value in numeric_positive.items() if value <= 0]
        if bad:
            raise ValueError(f"policy values must be positive: {', '.join(bad)}")
        if self.max_full > self.max_active:
            raise ValueError("max_full cannot exceed max_active")
        if self.threads_per_cell > self.max_worker_threads:
            raise ValueError("threads_per_cell cannot exceed max_worker_threads")


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
        text = ((result.stdout or "") + (result.stderr or "")).strip()
        record["version"] = text[:1000]
        record["version_exit_code"] = result.returncode
    except (OSError, subprocess.TimeoutExpired) as exc:
        record["version_error"] = str(exc)
    return record


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
) -> dict[str, Any]:
    files: dict[str, Any] = {}
    for label, path in (
        ("benchmark_registry", registry),
        ("dataset_registry", dataset_registry),
        ("baseline", baseline),
    ):
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


def select_ids(stage: str, *, baseline: Path, dataset_registry: Path,
               requested: Sequence[str] | None = None) -> list[str]:
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
        available = set(_baseline_keys(baseline, "subset"))
    elif stage.startswith("full"):
        available = set(_baseline_keys(baseline, "full"))
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


def build_cells(stage: str, ids: Sequence[str], *, run_dir: Path,
                policy: Policy, dataset_registry: Path) -> list[dict[str, Any]]:
    cells: list[dict[str, Any]] = []
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
    return cells


def _cell_fingerprint(cell: Mapping[str, Any], provenance_fingerprint: str) -> str:
    material = {
        "provenance": provenance_fingerprint,
        "kind": cell["kind"],
        "benchmark": cell["benchmark"],
        "mode": cell["mode"],
        "threads": cell["threads"],
        "command": cell["command"],
        "environment": cell.get("environment", {}),
        "artifacts": cell.get("artifacts", {}),
    }
    return canonical_hash(material)


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
) -> tuple[Path, dict[str, Any]]:
    policy.validate()
    provenance = collect_provenance(
        repo_root=repo_root, registry=registry,
        dataset_registry=dataset_registry, baseline=baseline,
    )
    if run_id is None:
        stamp = dt.datetime.now(dt.timezone.utc).strftime("%Y%m%dT%H%M%SZ")
        run_id = f"{stamp}-{stage}-{provenance['fingerprint'][:8]}"
    run_id = safe_name(run_id, limit=100)
    run_dir = Path(runs_root).resolve() / run_id
    ids = select_ids(stage, baseline=baseline,
                     dataset_registry=dataset_registry, requested=requested_ids)
    cells = build_cells(stage, ids, run_dir=run_dir, policy=policy,
                        dataset_registry=dataset_registry)
    for cell in cells:
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
    plan["fingerprint"] = canonical_hash({
        "schema_version": plan["schema_version"],
        "stage": stage,
        "ids": ids,
        "policy": plan["policy"],
        "provenance": provenance["fingerprint"],
        "cells": [cell["fingerprint"] for cell in cells],
    })
    return run_dir, plan


def initialize_run(run_dir: Path, plan: Mapping[str, Any]) -> None:
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
        for name in ("result_json", "gff", "manifest"):
            output = Path(cell["artifacts"][name]).resolve()
            if not output.is_relative_to(cell_dir):
                raise ValueError(
                    f"cell {cell['id']} artifact {name} escapes isolation: {output}"
                )


def collect_current_provenance(plan: Mapping[str, Any]) -> dict[str, Any]:
    inputs = plan["inputs"]
    return collect_provenance(
        repo_root=Path(plan["repo_root"]),
        registry=Path(inputs["registry"]),
        dataset_registry=Path(inputs["dataset_registry"]),
        baseline=Path(inputs["baseline"]),
    )


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


def _write_status(cell: Mapping[str, Any], state: str, **fields: Any) -> None:
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
    atomic_write_json(Path(cell["cell_dir"]) / "status.json", payload)


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


def _number(value: Any) -> bool:
    return isinstance(value, (int, float)) and not isinstance(value, bool)


def validate_result_schema(cell: Mapping[str, Any], path: Path) -> list[str]:
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
        matches = [row for row in (rows or []) if row.get("dataset") == benchmark]
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


def validate_artifacts(cell: Mapping[str, Any], started_ns: int) -> tuple[list[str], dict[str, Any]]:
    if cell["kind"] == "gate":
        return [], {"gate": "exit-code only"}
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
    atomic_write_json(Path(cell["cell_dir"]) / "gff_validation.json", validator_report)
    stats = {
        name: {
            "path": str(path),
            "size": path.stat().st_size,
            "mtime_ns": path.stat().st_mtime_ns,
        }
        for name, path in required.items()
    }
    return errors, {"artifacts": stats, "gff_validation": validator_report}


def performance_regressions(cell: Mapping[str, Any], baseline_path: Path) -> list[dict[str, Any]]:
    if cell["kind"] not in {"subset", "full_refresh"}:
        return []
    baseline = read_json(baseline_path)
    base = baseline.get(cell["artifacts"]["result_key"], {})
    result = read_json(Path(cell["artifacts"]["result_json"]))
    if cell["kind"] == "subset":
        current_record = result.get(cell["artifacts"]["result_key"], {})
        current_wall = current_record.get("wall_s", {}).get("lifton_devel")
        current_rss = current_record.get("peak_rss_mb", {}).get("lifton_devel")
    else:
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
        ratio = float(new) / float(old)
        if ratio > 1.25:
            regressions.append({
                "metric": metric, "baseline": old, "current": new,
                "ratio": ratio, "threshold": 1.25,
            })
    return regressions


def _unlink_markers(cell_dir: Path) -> None:
    for name in (".success", ".failed.json"):
        try:
            (cell_dir / name).unlink()
        except FileNotFoundError:
            pass


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
    started_ns = time.time_ns()
    started_at = utc_now()
    _unlink_markers(cell_dir)
    _write_status(
        cell, "running", attempts=attempt, started_at=started_at,
        started_ns=started_ns, session=cell_session_name(plan, cell),
        isolated_retry=bool(old_status.get("isolated_retry")),
    )
    environment = os.environ.copy()
    environment.update({
        "PYTHONNOUSERSITE": "1",
        "PYTHONHASHSEED": "0",
        "LIFTON_BENCH_THREADS": str(cell["threads"]),
    })
    environment.update({str(key): str(value) for key, value in cell.get("environment", {}).items()})
    stdout_path = cell_dir / f"attempt-{attempt:02d}.stdout.log"
    stderr_path = cell_dir / f"attempt-{attempt:02d}.stderr.log"
    started = time.monotonic()
    returncode: int
    launch_error: str | None = None
    try:
        with stdout_path.open("wb") as stdout, stderr_path.open("wb") as stderr:
            result = subprocess.run(
                cell["command"], cwd=plan["repo_root"], env=environment,
                stdout=stdout, stderr=stderr, check=False,
            )
        returncode = result.returncode
    except OSError as exc:
        returncode = 127
        launch_error = str(exc)
    elapsed = time.monotonic() - started
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
    }
    atomic_write_json(cell_dir / "exit.json", exit_document)

    errors = []
    validation: dict[str, Any] = {}
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
    atomic_write_json(success_path, success_document)
    _write_status(
        cell, "success", attempts=attempt, completed_at=success_document["completed_at"],
        isolated_retry=False, validation=validation, performance=regressions,
    )
    return 0


def _mark_failed(cell: Mapping[str, Any], errors: Sequence[str], *,
                 returncode: int | None, attempt: int,
                 validation: Mapping[str, Any] | None = None,
                 performance: Sequence[Mapping[str, Any]] | None = None) -> None:
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
    }
    atomic_write_json(Path(cell["cell_dir"]) / ".failed.json", document)
    _write_status(
        cell, "failed", attempts=attempt, failed_at=document["failed_at"],
        errors=list(errors), returncode=returncode, isolated_retry=False,
    )


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
    if len(active) >= policy.max_active:
        return False, "active-cell limit reached"
    if sum(int(item["threads"]) for item in active) + int(cell["threads"]) > policy.max_worker_threads:
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


def _active_cells(plan: Mapping[str, Any]) -> list[dict[str, Any]]:
    active = []
    for cell in plan["cells"]:
        status = _read_status(cell)
        if status.get("state") in ACTIVE_STATES and tmux_has_session(cell_session_name(plan, cell)):
            active.append({**cell, "isolated_retry": bool(status.get("isolated_retry"))})
    return active


def _mark_orphans(plan: Mapping[str, Any], *, grace_seconds: float = 5.0) -> None:
    now_ns = time.time_ns()
    for cell in plan["cells"]:
        status = _read_status(cell)
        if status.get("state") not in ACTIVE_STATES:
            continue
        if tmux_has_session(cell_session_name(plan, cell)):
            continue
        updated_ns = int(status.get("started_ns", 0))
        if updated_ns and now_ns - updated_ns < int(grace_seconds * 1e9):
            continue
        _mark_failed(
            cell, ["worker tmux session disappeared before a final status was written"],
            returncode=None, attempt=int(status.get("attempts", 0)),
        )


def _launch_cell(plan: Mapping[str, Any], cell: Mapping[str, Any]) -> None:
    status = _read_status(cell)
    session = cell_session_name(plan, cell)
    _write_status(
        cell, "launching", attempts=int(status.get("attempts", 0)),
        session=session, launch_requested_at=utc_now(),
        isolated_retry=bool(status.get("isolated_retry")),
        started_ns=time.time_ns(),
    )
    command = [
        sys.executable, str(Path(__file__).resolve()), "worker",
        "--run-dir", plan["run_dir"], "--cell", cell["id"],
    ]
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
    policy = Policy(**plan["policy"])
    log = RunLogger(Path(run_dir) / "controller.log")
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
    return {
        "schema_version": CONTROLLER_SCHEMA_VERSION,
        "run_id": plan["run_id"], "fingerprint": plan["fingerprint"],
        "stage": plan["stage"], "counts": counts, "cells": rows,
        "updated_at": utc_now(),
    }


def start_controller_tmux(plan: Mapping[str, Any]) -> str:
    session = controller_session_name(plan)
    if tmux_has_session(session):
        return session
    command = [
        sys.executable, str(Path(__file__).resolve()), "worker",
        "--run-dir", plan["run_dir"],
    ]
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
        for name, metadata in recorded.items():
            path = Path(metadata.get("path", ""))
            try:
                stat = path.stat()
            except OSError:
                errors.append(f"validated artifact disappeared: {name}={path}")
                continue
            if stat.st_size != metadata.get("size") or stat.st_mtime_ns != metadata.get("mtime_ns"):
                errors.append(f"validated artifact changed after success: {name}={path}")
        if deep and not errors and cell["kind"] != "gate":
            started_ns = int(success.get("exit", {}).get("started_ns", 0))
            deep_errors, _ = validate_artifacts(cell, started_ns)
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
    required_urls = [selected.reference_fa, selected.target_fa, selected.reference_gff]
    if selected.target_gff:
        required_urls.append(selected.target_gff)
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
    return Policy(
        threads_per_cell=args.threads_per_cell,
        max_active=args.max_active,
        max_full=args.max_full,
        max_worker_threads=args.max_worker_threads,
        load1_limit=args.load1_limit,
        min_available_gib=args.min_available_gib,
        stagger_seconds=args.stagger_seconds,
        poll_seconds=args.poll_seconds,
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
        "--stage", choices=("gates", "subset-canary", "subset", "full-canary", "full",
                            "e2e-canary", "e2e"), default=None,
    )
    start.add_argument("--ids", nargs="+")
    start.add_argument("--runs-root", default=str(DEFAULT_RUNS_ROOT))
    start.add_argument("--registry", default=str(DEFAULT_REGISTRY))
    start.add_argument("--dataset-registry", default=str(DEFAULT_DATASET_REGISTRY))
    start.add_argument("--baseline", default=str(DEFAULT_BASELINE))
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
    parser.add_argument("--threads-per-cell", type=int, default=8)
    parser.add_argument("--max-active", type=int, default=4)
    parser.add_argument("--max-full", type=int, default=2)
    parser.add_argument("--max-worker-threads", type=int, default=32)
    parser.add_argument("--load1-limit", type=float, default=32.0)
    parser.add_argument("--min-available-gib", type=float, default=256.0)
    parser.add_argument("--stagger-seconds", type=float, default=15.0)
    parser.add_argument("--poll-seconds", type=float, default=30.0)


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    if args.action == "start":
        runs_root = Path(args.runs_root).resolve()
        if args.run_id and (runs_root / safe_name(args.run_id, limit=100) / "plan.json").exists():
            run_dir = runs_root / safe_name(args.run_id, limit=100)
            plan = load_plan(run_dir)
            if args.stage and args.stage != plan["stage"]:
                parser.error(f"existing run stage is {plan['stage']!r}, not {args.stage!r}")
            if args.ids and list(args.ids) != plan["ids"]:
                parser.error("existing run ids do not match the requested ids")
            try:
                assert_matching_provenance(plan)
            except RuntimeError as exc:
                parser.error(str(exc))
        else:
            stage = args.stage or "subset-canary"
            try:
                run_dir, plan = create_plan(
                    run_id=args.run_id, stage=stage, requested_ids=args.ids,
                    runs_root=runs_root, repo_root=REPO_ROOT,
                    registry=Path(args.registry),
                    dataset_registry=Path(args.dataset_registry),
                    baseline=Path(args.baseline), policy=_policy_from_args(args),
                )
                initialize_run(run_dir, plan)
            except (OSError, ValueError) as exc:
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

    if args.action in {"status", "retry", "reconcile"}:
        if not args.run_dir and not args.run_id:
            parser.error("provide run_id or --run-dir")
        run_dir = _run_dir_argument(args)
        if args.action == "status":
            summary = summarize_run(load_plan(run_dir))
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
            summary = reconcile_run(run_dir, deep=args.deep)
        if args.json:
            print(json.dumps(summary, indent=2, sort_keys=True))
        else:
            _print_summary(summary)
            if summary.get("provenance_error"):
                print(f"provenance: MISMATCH ({summary['provenance_error']})")
        if summary.get("provenance_error"):
            return 2
        return 0 if summary["counts"].get("failed", 0) == 0 else 1

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
