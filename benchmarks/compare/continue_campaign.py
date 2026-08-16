"""Create and optionally launch a disjoint continuation of a benchmark run."""

from __future__ import annotations

import argparse
import copy
import datetime as dt
import json
import os
from pathlib import Path
import subprocess
import sys
from typing import Any, Mapping, Sequence

from . import build_controller


DEFAULT_STATES = ("pending", "failed", "retry_pending")
TERMINAL_STATES = frozenset(DEFAULT_STATES)


def _replace_paths(value: Any, old_root: str, new_root: str) -> Any:
    if isinstance(value, str):
        return value.replace(old_root, new_root)
    if isinstance(value, list):
        return [_replace_paths(item, old_root, new_root) for item in value]
    if isinstance(value, dict):
        return {
            key: _replace_paths(item, old_root, new_root)
            for key, item in value.items()
        }
    return value


def _replace_threads(command: Sequence[str], threads: int) -> list[str]:
    output = list(command)
    for index, token in enumerate(output[:-1]):
        if token in {"--threads", "-t"}:
            output[index + 1] = str(threads)
    return output


def _selected_cell_ids(
    plan: Mapping[str, Any],
    requested: Sequence[str] | None,
    states: Sequence[str],
) -> list[str]:
    cell_by_id = {cell["id"]: cell for cell in plan["cells"]}
    status_by_id = {
        cell_id: build_controller.read_json(
            Path(cell["cell_dir"]) / "status.json"
        )
        for cell_id, cell in cell_by_id.items()
    }
    allowed = set(states)
    unknown_states = allowed - TERMINAL_STATES
    if unknown_states:
        raise ValueError(
            f"continuation states are not terminal: {sorted(unknown_states)}"
        )
    if requested:
        cell_ids = list(dict.fromkeys(requested))
        unknown = set(cell_ids) - set(cell_by_id)
        if unknown:
            raise ValueError(f"unknown cell ids: {sorted(unknown)}")
    else:
        cell_ids = [
            cell_id for cell_id, status in status_by_id.items()
            if status.get("state") in allowed
        ]
    active = [
        cell_id for cell_id in cell_ids
        if status_by_id[cell_id].get("state") in build_controller.ACTIVE_STATES
    ]
    if active:
        raise RuntimeError(
            "refusing to continue active cells: " + ", ".join(sorted(active))
        )
    invalid = [
        cell_id for cell_id in cell_ids
        if status_by_id[cell_id].get("state") not in TERMINAL_STATES
    ]
    if invalid:
        raise RuntimeError(
            "continuation cells must be pending, failed, or retry_pending: "
            + ", ".join(sorted(invalid))
        )
    if not cell_ids:
        raise RuntimeError("no eligible terminal cells remain")
    return cell_ids


def _continuation_policy(
    original: Mapping[str, Any],
    *,
    max_active: int,
    max_full: int,
    max_worker_threads: int,
    load1_limit: float,
    min_available_gib: float,
) -> build_controller.Policy:
    values = dict(original)
    values.update({
        "max_active": max_active,
        "max_full": max_full,
        "max_worker_threads": max_worker_threads,
        "load1_limit": load1_limit,
        "min_available_gib": min_available_gib,
    })
    policy = build_controller.Policy(**values)
    policy.validate()
    return policy


def _start_controller(plan: Mapping[str, Any], tooling_root: Path | None) -> str:
    if tooling_root is None:
        return build_controller.start_controller_tmux(plan)
    session = build_controller.controller_session_name(plan)
    if build_controller.tmux_has_session(session):
        return session
    controller = Path(tooling_root).resolve() / "benchmarks" / "compare" / "build_controller.py"
    command = build_controller._python_worker_command(
        str(controller), "worker", "--run-dir", plan["run_dir"],
    )
    build_controller._tmux_new_session(session, command)
    return session


def _collect_provenance(
    plan: Mapping[str, Any], tooling_root: Path | None,
) -> dict[str, Any]:
    if tooling_root is None:
        return build_controller.collect_current_provenance(plan)
    code = (
        "import json, sys; "
        "from benchmarks.compare import build_controller; "
        "plan = json.load(sys.stdin); "
        "print(json.dumps(build_controller.collect_current_provenance(plan)))"
    )
    environment = dict(os.environ)
    environment.pop("PYTHONPATH", None)
    environment.update({
        "PYTHONDONTWRITEBYTECODE": "1",
        "PYTHONNOUSERSITE": "1",
    })
    result = subprocess.run(
        [sys.executable, "-B", "-c", code],
        cwd=str(Path(tooling_root).resolve()),
        env=environment,
        input=json.dumps(plan),
        capture_output=True,
        text=True,
        check=False,
    )
    if result.returncode != 0:
        raise RuntimeError(
            "snapshot provenance collection failed: "
            + (result.stderr or result.stdout).strip()
        )
    try:
        return json.loads(result.stdout)
    except (TypeError, ValueError) as exc:
        raise RuntimeError("snapshot provenance output was not JSON") from exc


def build_continuation_plan(
    original: Mapping[str, Any],
    cell_ids: Sequence[str],
    *,
    run_dir: Path,
    policy: build_controller.Policy,
    tooling_root: Path | None = None,
) -> dict[str, Any]:
    original_run_dir = Path(original["run_dir"]).resolve()
    new_run_dir = Path(run_dir).resolve()
    original_tooling_root = Path(original["repo_root"]).resolve()
    new_tooling_root = (
        Path(tooling_root).resolve()
        if tooling_root is not None else original_tooling_root
    )
    selected = [
        copy.deepcopy(build_controller._cell_for(original, cell_id))
        for cell_id in cell_ids
    ]
    provenance = copy.deepcopy(original["provenance"])
    provenance.pop("campaign_case", None)
    provenance["fingerprint"] = build_controller.canonical_hash({
        key: value for key, value in provenance.items()
        if key != "fingerprint"
    })
    old_root = str(original_run_dir)
    new_root = str(new_run_dir)
    for cell in selected:
        cell["cell_dir"] = _replace_paths(cell["cell_dir"], old_root, new_root)
        cell["command"] = _replace_threads(
            _replace_paths(cell["command"], old_root, new_root),
            policy.threads_per_cell,
        )
        cell["command"] = _replace_paths(
            cell["command"], str(original_tooling_root), str(new_tooling_root),
        )
        cell["threads"] = policy.threads_per_cell
        for key, value in list(cell.get("artifacts", {}).items()):
            cell["artifacts"][key] = _replace_paths(value, old_root, new_root)
        if "environment" in cell:
            cell["environment"] = _replace_paths(
                cell["environment"],
                str(original_tooling_root),
                str(new_tooling_root),
            )
        if "paired" in cell:
            cell["paired"] = _replace_paths(
                cell["paired"],
                str(original_tooling_root),
                str(new_tooling_root),
            )
    plan = copy.deepcopy(dict(original))
    plan.pop("campaign_case", None)
    plan["run_id"] = new_run_dir.name
    plan["run_dir"] = str(new_run_dir)
    plan["created_at"] = build_controller.utc_now()
    plan["ids"] = list(dict.fromkeys(cell["benchmark"] for cell in selected))
    plan["policy"] = build_controller.asdict(policy)
    plan["repo_root"] = str(new_tooling_root)
    plan["inputs"] = _replace_paths(
        plan["inputs"], str(original_tooling_root), str(new_tooling_root),
    )
    if "paired" in plan:
        plan["paired"] = _replace_paths(
            plan["paired"],
            str(original_tooling_root),
            str(new_tooling_root),
        )
    plan["provenance"] = _replace_paths(
        provenance,
        str(original_tooling_root),
        str(new_tooling_root),
    )
    plan["cells"] = selected
    if tooling_root is not None:
        plan["provenance"] = _collect_provenance(plan, tooling_root)
    for cell in plan["cells"]:
        cell["fingerprint"] = build_controller._cell_fingerprint(
            cell, plan["provenance"]["fingerprint"]
        )
    plan["continuation_of"] = {
        "run_id": original["run_id"],
        "run_dir": str(original_run_dir),
        "plan_fingerprint": original["fingerprint"],
        "cell_ids": list(cell_ids),
        "created_at": build_controller.utc_now(),
    }
    plan["fingerprint"] = build_controller._plan_fingerprint(plan)
    build_controller.validate_plan_integrity(plan)
    build_controller.validate_plan_layout(plan)
    return plan


def continue_run(
    run_dir: Path,
    *,
    cell_ids: Sequence[str] | None = None,
    states: Sequence[str] = DEFAULT_STATES,
    output_run_dir: Path | None = None,
    tooling_root: Path | None = None,
    max_active: int = 16,
    max_full: int = 4,
    max_worker_threads: int = 128,
    load1_limit: float = 112.0,
    min_available_gib: float = 384.0,
    dry_run: bool = False,
    foreground: bool = False,
) -> dict[str, Any]:
    original = build_controller.load_plan(Path(run_dir).resolve())
    selected = _selected_cell_ids(original, cell_ids, states)
    policy = _continuation_policy(
        original["policy"],
        max_active=max_active,
        max_full=max_full,
        max_worker_threads=max_worker_threads,
        load1_limit=load1_limit,
        min_available_gib=min_available_gib,
    )
    if output_run_dir is None:
        stamp = dt.datetime.now(dt.timezone.utc).strftime("%Y%m%dT%H%M%SZ")
        output_run_dir = (
            Path(original["run_dir"]).resolve().parent
            / "continuations"
            / f"{original['run_id']}-continuation-{stamp}"
        )
    plan = build_continuation_plan(
        original,
        selected,
        run_dir=output_run_dir,
        policy=policy,
        tooling_root=tooling_root,
    )
    result = {
        "run_dir": str(Path(plan["run_dir"]).resolve()),
        "run_id": plan["run_id"],
        "plan_fingerprint": plan["fingerprint"],
        "source_plan_fingerprint": original["fingerprint"],
        "cell_ids": selected,
        "policy": dict(plan["policy"]),
        "dry_run": dry_run,
    }
    if dry_run:
        return result
    build_controller.initialize_run(Path(plan["run_dir"]), plan)
    if foreground:
        result["exit_code"] = build_controller.scheduler_loop(Path(plan["run_dir"]))
    else:
        result["tmux_session"] = _start_controller(plan, tooling_root)
    return result


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", required=True)
    parser.add_argument("--cells", nargs="+")
    parser.add_argument("--states", nargs="+", default=list(DEFAULT_STATES))
    parser.add_argument("--output-run-dir")
    parser.add_argument("--tooling-root")
    parser.add_argument("--max-active", type=int, default=16)
    parser.add_argument("--max-full", type=int, default=4)
    parser.add_argument("--max-worker-threads", type=int, default=128)
    parser.add_argument("--load1-limit", type=float, default=112.0)
    parser.add_argument("--min-available-gib", type=float, default=384.0)
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--foreground", action="store_true")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    result = continue_run(
        Path(args.run_dir),
        cell_ids=args.cells,
        states=args.states,
        output_run_dir=(Path(args.output_run_dir) if args.output_run_dir else None),
        tooling_root=(Path(args.tooling_root) if args.tooling_root else None),
        max_active=args.max_active,
        max_full=args.max_full,
        max_worker_threads=args.max_worker_threads,
        load1_limit=args.load1_limit,
        min_available_gib=args.min_available_gib,
        dry_run=args.dry_run,
        foreground=args.foreground,
    )
    print(json.dumps(result, indent=2, sort_keys=True))
    return int(result.get("exit_code", 0))


if __name__ == "__main__":
    raise SystemExit(main())
