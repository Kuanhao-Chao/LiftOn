"""Deterministic schedules and summaries for LiftOn protocol benchmarks.

The release comparison uses paired AB/BA cells.  Thread scaling and fast-path
mode coverage instead exercise one immutable candidate SHA repeatedly.  This
module defines those order-balanced schedules and validates their completed
records without depending on the controller implementation.
"""
from __future__ import annotations

import math
import statistics
from collections import defaultdict
from typing import Any, Iterable, Mapping

from .release_evaluation import E2E_MODES


SCHEMA_VERSION = 1
THREAD_LEVELS = (1, 2, 4, 8, 16, 32)
THREAD_WILLIAMS_ORDERS = (
    (1, 2, 32, 4, 16, 8),
    (2, 4, 1, 8, 32, 16),
    (4, 8, 2, 16, 1, 32),
    (8, 16, 4, 32, 2, 1),
    (16, 32, 8, 1, 4, 2),
    (32, 1, 16, 2, 8, 4),
)
IO_MODE_ORDER = tuple(E2E_MODES)
IO_ROTATIONS = (0, 4, 2, 6)


def _rotated(values: tuple[Any, ...], offset: int) -> tuple[Any, ...]:
    return values[offset:] + values[:offset]


def thread_scaling_cases(dataset: str = "bee") -> list[dict[str, Any]]:
    """Return the fixed 36-arm Williams schedule."""

    return [
        {
            "case_id": (
                f"thread_scaling__{dataset}__repetition_{repetition:02d}"
                f"__position_{position:02d}__threads_{threads:02d}"
            ),
            "kind": "thread_scaling",
            "dataset": dataset,
            "repetition": repetition,
            "position": position,
            "threads": threads,
            "mode": "fast",
        }
        for repetition, order in enumerate(THREAD_WILLIAMS_ORDERS, start=1)
        for position, threads in enumerate(order, start=1)
    ]


def io_mode_cases(dataset: str = "arabidopsis") -> list[dict[str, Any]]:
    """Return the fixed 32-arm rotation-balanced mode schedule."""

    cases = []
    for repetition, offset in enumerate(IO_ROTATIONS, start=1):
        for position, mode in enumerate(_rotated(IO_MODE_ORDER, offset), start=1):
            cases.append({
                "case_id": (
                    f"io_modes__{dataset}__repetition_{repetition:02d}"
                    f"__position_{position:02d}__{mode}"
                ),
                "kind": "io_modes",
                "dataset": dataset,
                "repetition": repetition,
                "position": position,
                "threads": 8,
                "mode": mode,
            })
    return cases


def protocol_cases(kind: str) -> list[dict[str, Any]]:
    if kind == "thread_scaling":
        return thread_scaling_cases()
    if kind == "io_modes":
        return io_mode_cases()
    raise ValueError(f"unknown protocol campaign: {kind}")


def _positive_number(value: Any, *, field: str, case_id: str) -> float:
    if (
        not isinstance(value, (int, float))
        or isinstance(value, bool)
        or not math.isfinite(float(value))
        or float(value) <= 0
    ):
        raise ValueError(f"{case_id}: {field} must be positive and finite")
    return float(value)


def _optional_counter(value: Any, *, field: str, case_id: str) -> int | None:
    if value is None:
        return None
    if not isinstance(value, int) or isinstance(value, bool) or value < 0:
        raise ValueError(f"{case_id}: {field} must be a non-negative integer or null")
    return value


def _validate_records(
    records: Iterable[Mapping[str, Any]],
    *,
    kind: str,
) -> list[dict[str, Any]]:
    expected = {case["case_id"]: case for case in protocol_cases(kind)}
    observed: dict[str, dict[str, Any]] = {}
    for item in records:
        row = dict(item)
        case_id = row.get("case_id")
        if not isinstance(case_id, str) or case_id not in expected:
            raise ValueError(f"unexpected {kind} case ID: {case_id!r}")
        if case_id in observed:
            raise ValueError(f"duplicate {kind} case ID: {case_id}")
        for field in ("kind", "dataset", "repetition", "position", "threads", "mode"):
            if row.get(field) != expected[case_id][field]:
                raise ValueError(
                    f"{case_id}: {field} does not match the frozen protocol schedule"
                )
        profile = row.get("profile")
        if not isinstance(profile, Mapping):
            raise ValueError(f"{case_id}: profile is missing")
        _positive_number(
            profile.get("wall_clock_seconds"),
            field="wall_clock_seconds",
            case_id=case_id,
        )
        _positive_number(
            profile.get("peak_rss_mb"),
            field="peak_rss_mb",
            case_id=case_id,
        )
        for field in (
            "filesystem_inputs",
            "filesystem_outputs",
            "major_page_faults",
            "minor_page_faults",
        ):
            _optional_counter(profile.get(field), field=field, case_id=case_id)
        fingerprints = row.get("fingerprints")
        if (
            not isinstance(fingerprints, Mapping)
            or not isinstance(fingerprints.get("byte_sha256"), str)
            or not isinstance(fingerprints.get("semantic_sha256"), str)
        ):
            raise ValueError(f"{case_id}: output fingerprints are missing")
        candidate_sha = row.get("candidate_sha")
        if not isinstance(candidate_sha, str) or len(candidate_sha) != 40:
            raise ValueError(f"{case_id}: candidate SHA is malformed")
        observed[case_id] = row
    missing = sorted(set(expected) - set(observed))
    if missing:
        raise ValueError(
            f"{kind} campaign is incomplete: {len(missing)} missing cases; "
            f"first={missing[0]}"
        )
    return [observed[case["case_id"]] for case in protocol_cases(kind)]


def _median_optional(rows: list[Mapping[str, Any]], field: str) -> float | None:
    values = [
        row["profile"][field]
        for row in rows
        if row["profile"].get(field) is not None
    ]
    if len(values) != len(rows):
        return None
    return float(statistics.median(values))


def summarize_protocol(
    records: Iterable[Mapping[str, Any]],
    *,
    kind: str,
) -> dict[str, Any]:
    """Validate a complete protocol campaign and calculate robust summaries."""

    rows = _validate_records(records, kind=kind)
    candidate_shas = {row["candidate_sha"] for row in rows}
    if len(candidate_shas) != 1:
        raise ValueError("protocol campaign mixes candidate SHAs")
    byte_hashes = {row["fingerprints"]["byte_sha256"] for row in rows}
    semantic_hashes = {row["fingerprints"]["semantic_sha256"] for row in rows}
    if len(byte_hashes) != 1 or len(semantic_hashes) != 1:
        raise ValueError("protocol outputs are not byte- and semantic-identical")

    axis = "threads" if kind == "thread_scaling" else "mode"
    grouped: dict[Any, list[Mapping[str, Any]]] = defaultdict(list)
    for row in rows:
        grouped[row[axis]].append(row)
    summaries = []
    baseline_wall = None
    if kind == "thread_scaling":
        baseline_wall = statistics.median(
            row["profile"]["wall_clock_seconds"] for row in grouped[1]
        )
    order = THREAD_LEVELS if kind == "thread_scaling" else IO_MODE_ORDER
    for value in order:
        group = grouped[value]
        wall = float(statistics.median(
            row["profile"]["wall_clock_seconds"] for row in group
        ))
        rss = float(statistics.median(
            row["profile"]["peak_rss_mb"] for row in group
        ))
        cpu_utilization = float(statistics.median(
            (
                float(row["profile"].get("user_cpu_seconds", 0.0))
                + float(row["profile"].get("sys_cpu_seconds", 0.0))
            )
            / float(row["profile"]["wall_clock_seconds"])
            for row in group
        ))
        summary = {
            axis: value,
            "repetitions": len(group),
            "median_wall_clock_seconds": wall,
            "median_peak_rss_mb": rss,
            "median_cpu_utilization": cpu_utilization,
            "median_filesystem_inputs": _median_optional(group, "filesystem_inputs"),
            "median_filesystem_outputs": _median_optional(group, "filesystem_outputs"),
            "median_major_page_faults": _median_optional(group, "major_page_faults"),
            "median_minor_page_faults": _median_optional(group, "minor_page_faults"),
        }
        if kind == "thread_scaling":
            speedup = float(baseline_wall) / wall
            summary["speedup"] = speedup
            summary["parallel_efficiency"] = speedup / int(value)
        summaries.append(summary)
    return {
        "schema_version": SCHEMA_VERSION,
        "kind": kind,
        "candidate_sha": next(iter(candidate_shas)),
        "n_arms": len(rows),
        "byte_sha256": next(iter(byte_hashes)),
        "semantic_sha256": next(iter(semantic_hashes)),
        "outputs_identical": True,
        "summaries": summaries,
    }
