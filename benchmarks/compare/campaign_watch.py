"""Monitor multiple detached LiftOn benchmark controller runs."""

from __future__ import annotations

import argparse
import datetime as dt
import json
import os
import re
import shutil
import socket
import time
from collections import Counter
from pathlib import Path
from typing import Any, Mapping, Sequence

from . import build_controller


def _utc_now() -> str:
    return dt.datetime.now(dt.timezone.utc).isoformat().replace("+00:00", "Z")


def _memory_available_gib() -> float:
    for line in Path("/proc/meminfo").read_text(encoding="utf-8").splitlines():
        if line.startswith("MemAvailable:"):
            return float(line.split()[1]) / (1024.0 * 1024.0)
    raise RuntimeError("MemAvailable is absent from /proc/meminfo")


def _cpu_iowait_percent_from_stat(stat_text: str) -> float:
    line = next(
        line for line in stat_text.splitlines()
        if line.startswith("cpu ")
    )
    fields = [int(value) for value in line.split()[1:]]
    # guest and guest_nice are already included in user and nice.
    total = sum(fields[:8])
    iowait = fields[4] if len(fields) > 4 else 0
    return 100.0 * iowait / total if total else 0.0


def _cpu_iowait_percent() -> float:
    return _cpu_iowait_percent_from_stat(
        Path("/proc/stat").read_text()
    )


def _io_pressure_avg10() -> float | None:
    try:
        text = Path("/proc/pressure/io").read_text(encoding="utf-8")
    except OSError:
        return None
    match = re.search(r"^some\s+avg10=([0-9.]+)", text, re.MULTILINE)
    return float(match.group(1)) if match else None


def _process_state_counts() -> dict[str, int]:
    counts = Counter()
    for entry in Path("/proc").iterdir():
        if not entry.name.isdigit():
            continue
        try:
            stat = (entry / "stat").read_text(encoding="utf-8")
            fields = stat[stat.rfind(")") + 2:].split()
        except (OSError, ValueError):
            continue
        if fields:
            counts[fields[0]] += 1
    return {
        "d_state_processes": counts.get("D", 0),
        "zombie_processes": counts.get("Z", 0),
    }


def _resource_snapshot() -> dict[str, Any]:
    return {
        "load1": os.getloadavg()[0],
        "available_gib": _memory_available_gib(),
        "iowait_percent_since_boot": _cpu_iowait_percent(),
        "io_pressure_some_avg10": _io_pressure_avg10(),
        **_process_state_counts(),
    }


def _parse_time(value: Any) -> dt.datetime | None:
    if not isinstance(value, str):
        return None
    try:
        parsed = dt.datetime.fromisoformat(value.replace("Z", "+00:00"))
    except ValueError:
        return None
    return parsed if parsed.tzinfo else parsed.replace(tzinfo=dt.timezone.utc)


def _elapsed_seconds(status: Mapping[str, Any]) -> float | None:
    started = _parse_time(
        status.get("started_at") or status.get("launch_requested_at")
    )
    if started is None:
        started_ns = status.get("started_ns")
        if isinstance(started_ns, int) and started_ns > 0:
            return max(0.0, time.time() - started_ns / 1e9)
        return None
    return max(0.0, (dt.datetime.now(dt.timezone.utc) - started).total_seconds())


def _progress_record(cell: Mapping[str, Any]) -> dict[str, Any]:
    cell_dir = Path(str(cell["cell_dir"]))
    progress = cell_dir / "progress.log"
    record: dict[str, Any] = {"path": str(progress), "exists": progress.is_file()}
    if progress.is_file():
        stat = progress.stat()
        record.update({"size": stat.st_size, "mtime": stat.st_mtime})
    return record


def _session_name(
    plan: Mapping[str, Any],
    cell: Mapping[str, Any],
    status: Mapping[str, Any],
) -> str:
    persisted = status.get("session")
    if isinstance(persisted, str) and persisted:
        return persisted
    return build_controller.cell_session_name(plan, cell)


def _tmux_session_alive(name: str) -> bool:
    candidates = [name]
    normalized = name.replace(".", "_")
    if normalized != name:
        candidates.append(normalized)
    return any(
        build_controller.tmux_has_session(candidate)
        for candidate in candidates
    )


def _read_run(run_dir: Path) -> dict[str, Any]:
    plan = build_controller.load_plan(run_dir)
    rows = []
    for cell in plan["cells"]:
        status = build_controller.read_json(Path(cell["cell_dir"]) / "status.json")
        state = status.get("state", "unknown")
        row = {
            "cell_id": cell["id"],
            "benchmark": cell.get("benchmark"),
            "repetition": cell.get("repetition"),
            "state": state,
            "threads": int(cell.get("threads", 0)),
            "full_job": bool(cell.get("full_job")),
            "attempts": status.get("attempts", 0),
            "updated_at": status.get("updated_at"),
            "elapsed_seconds": _elapsed_seconds(status),
        }
        if state in build_controller.ACTIVE_STATES:
            session = _session_name(plan, cell, status)
            row["tmux_session"] = session
            row["session_alive"] = _tmux_session_alive(session)
            row["worker_alive"] = bool(
                build_controller._worker_identity_matches(status)
            )
            row["progress"] = _progress_record(cell)
        rows.append(row)
    counts = Counter(row["state"] for row in rows)
    return {
        "run_id": plan["run_id"],
        "run_dir": str(run_dir),
        "stage": plan.get("stage"),
        "plan_fingerprint": plan.get("fingerprint"),
        "policy": dict(plan.get("policy", {})),
        "counts": dict(sorted(counts.items())),
        "active": [row for row in rows if row["state"] in build_controller.ACTIVE_STATES],
        "errors": [],
    }


def _free_space_gib(path: Path) -> float:
    return shutil.disk_usage(path).free / (1024.0 * 1024.0 * 1024.0)


def snapshot(
    run_dirs: Sequence[Path],
    *,
    max_active: int = 12,
    max_worker_threads: int = 96,
    load1_limit: float = 112.0,
    min_available_gib: float = 384.0,
) -> dict[str, Any]:
    runs = []
    for raw_dir in run_dirs:
        run_dir = Path(raw_dir).resolve()
        try:
            runs.append(_read_run(run_dir))
        except (OSError, RuntimeError, TypeError, ValueError) as exc:
            runs.append({
                "run_id": run_dir.name,
                "run_dir": str(run_dir),
                "counts": {},
                "active": [],
                "errors": [str(exc)],
            })
    active = [row for run in runs for row in run["active"]]
    active_threads = sum(int(row["threads"]) for row in active)
    resources = _resource_snapshot()
    alerts = []
    if len(active) > max_active:
        alerts.append(f"active cells {len(active)} exceed limit {max_active}")
    if active_threads > max_worker_threads:
        alerts.append(
            f"worker threads {active_threads} exceed limit {max_worker_threads}"
        )
    if resources["load1"] >= load1_limit:
        alerts.append(f"load1 {resources['load1']:.2f} reaches limit {load1_limit:.2f}")
    if resources["available_gib"] < min_available_gib:
        alerts.append(
            f"available memory {resources['available_gib']:.1f} GiB below "
            f"floor {min_available_gib:.1f} GiB"
        )
    if resources.get("d_state_processes", 0) > 0:
        alerts.append(
            f"{resources['d_state_processes']} process(es) in uninterruptible I/O wait"
        )
    run_counts = Counter()
    for run in runs:
        run_counts.update(run.get("counts", {}))
        alerts.extend(f"{run['run_id']}: {error}" for error in run["errors"])
    return {
        "schema_version": 1,
        "observed_at": _utc_now(),
        "host": socket.gethostname(),
        "resources": {
            **resources,
            "free_space_gib": _free_space_gib(Path.cwd()),
        },
        "limits": {
            "max_active": max_active,
            "max_worker_threads": max_worker_threads,
            "load1_limit": load1_limit,
            "min_available_gib": min_available_gib,
        },
        "totals": {
            "counts": dict(sorted(run_counts.items())),
            "active_cells": len(active),
            "active_threads": active_threads,
        },
        "alerts": alerts,
        "runs": runs,
    }


def _atomic_write_json(path: Path, payload: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp-{os.getpid()}")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    os.replace(temporary, path)


def _write_evidence_manifest(
    path: Path,
    payload: Mapping[str, Any],
    *,
    latest: Path,
    history: Path | None,
) -> None:
    previous: Mapping[str, Any] = {}
    if path.is_file():
        try:
            loaded = json.loads(path.read_text(encoding="utf-8"))
            if isinstance(loaded, Mapping):
                previous = loaded
        except (OSError, TypeError, ValueError):
            previous = {}
    plans = []
    for run in payload.get("runs", []):
        if not isinstance(run, Mapping):
            continue
        plans.append({
            "run_id": run.get("run_id"),
            "run_dir": run.get("run_dir"),
            "stage": run.get("stage"),
            "plan_fingerprint": run.get("plan_fingerprint"),
            "policy": run.get("policy", {}),
        })
    manifest = {
        "schema_version": 1,
        "kind": "campaign_evidence_manifest",
        "host": payload.get("host"),
        "run_plans": plans,
        "limits": payload.get("limits", {}),
        "first_observed_at": previous.get(
            "first_observed_at", payload.get("observed_at")
        ),
        "last_observed_at": payload.get("observed_at"),
        "snapshot_count": int(previous.get("snapshot_count", 0)) + 1,
        "latest_snapshot": str(latest.resolve()),
        "history": str(history.resolve()) if history else None,
        "latest_counts": payload.get("totals", {}),
        "latest_alerts": payload.get("alerts", []),
    }
    _atomic_write_json(path, manifest)


def write_snapshot(
    run_dirs: Sequence[Path],
    *,
    latest: Path,
    history: Path | None = None,
    evidence: Path | None = None,
    **limits: Any,
) -> dict[str, Any]:
    payload = snapshot(run_dirs, **limits)
    _atomic_write_json(latest, payload)
    if history is not None:
        history.parent.mkdir(parents=True, exist_ok=True)
        with history.open("a", encoding="utf-8") as handle:
            handle.write(json.dumps(payload, sort_keys=True) + "\n")
    if evidence is not None:
        _write_evidence_manifest(
            evidence, payload, latest=latest, history=history,
        )
        payload = {**payload, "evidence_manifest": str(evidence.resolve())}
        _atomic_write_json(latest, payload)
    return payload


def _positive_float(value: str) -> float:
    parsed = float(value)
    if parsed <= 0:
        raise argparse.ArgumentTypeError("must be positive")
    return parsed


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", action="append", required=True)
    parser.add_argument("--latest", required=True)
    parser.add_argument("--history")
    parser.add_argument("--evidence")
    parser.add_argument("--interval", type=_positive_float, default=60.0)
    parser.add_argument("--once", action="store_true")
    parser.add_argument("--max-active", type=int, default=12)
    parser.add_argument("--max-worker-threads", type=int, default=96)
    parser.add_argument("--load1-limit", type=_positive_float, default=112.0)
    parser.add_argument("--min-available-gib", type=_positive_float, default=384.0)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    run_dirs = [Path(value).resolve() for value in args.run_dir]
    limits = {
        "max_active": args.max_active,
        "max_worker_threads": args.max_worker_threads,
        "load1_limit": args.load1_limit,
        "min_available_gib": args.min_available_gib,
    }
    while True:
        payload = write_snapshot(
            run_dirs,
            latest=Path(args.latest).resolve(),
            history=Path(args.history).resolve() if args.history else None,
            evidence=Path(args.evidence).resolve() if args.evidence else None,
            **limits,
        )
        print(
            f"{payload['observed_at']} active={payload['totals']['active_cells']} "
            f"threads={payload['totals']['active_threads']} "
            f"states={payload['totals']['counts']} alerts={payload['alerts']}",
            flush=True,
        )
        if args.once:
            return 0
        time.sleep(args.interval)


if __name__ == "__main__":
    raise SystemExit(main())
