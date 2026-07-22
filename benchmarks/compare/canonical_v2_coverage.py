#!/usr/bin/env python3
"""Generate the deterministic canonical-v2 coverage and gap matrix."""
from __future__ import annotations

import argparse
import json
from collections import Counter
from pathlib import Path
from typing import Any, Mapping, Sequence

from benchmarks import manifest_tools
from benchmarks.compare import campaign_profiles


HERE = Path(__file__).resolve().parent
DEFAULT_MANIFEST = manifest_tools.DEFAULT_MANIFEST
DEFAULT_PROFILES = campaign_profiles.DEFAULT_PROFILE_REGISTRY
DEFAULT_BACKLOG = HERE / "canonical_v3_backlog.json"
DEFAULT_JSON = HERE / "canonical_v2_coverage.json"
DEFAULT_MARKDOWN = HERE / "canonical_v2_coverage.md"


def _load_json(path: Path) -> dict[str, Any]:
    document = json.loads(Path(path).read_text(encoding="utf-8"))
    if not isinstance(document, dict):
        raise ValueError(f"expected an object: {path}")
    return document


def _count(values) -> dict[str, int]:
    return dict(sorted(Counter(str(value) for value in values).items()))


def _validate_backlog(document: Mapping[str, Any]) -> None:
    if set(document) != {
        "schema_version", "profile_id", "promotion_requirements", "items",
    } or document["schema_version"] != 1:
        raise ValueError("unsupported canonical-v3 backlog schema")
    items = document["items"]
    if not isinstance(items, list) or not items:
        raise ValueError("canonical-v3 backlog must contain items")
    ids = []
    for item in items:
        if not isinstance(item, Mapping) or set(item) != {
            "id", "priority", "dimension", "description",
        }:
            raise ValueError("canonical-v3 backlog item is malformed")
        if item["priority"] not in {"P1", "P2", "P3"}:
            raise ValueError(f"invalid backlog priority: {item['priority']}")
        ids.append(item["id"])
    if len(ids) != len(set(ids)):
        raise ValueError("canonical-v3 backlog IDs are duplicated")


def build_coverage(
    manifest_path: Path = DEFAULT_MANIFEST,
    profile_registry: Path = DEFAULT_PROFILES,
    backlog_path: Path = DEFAULT_BACKLOG,
) -> dict[str, Any]:
    manifest = manifest_tools.validate_manifest(_load_json(manifest_path))
    profile = campaign_profiles.load_profile(
        "canonical-v2", registry=profile_registry,
    )
    backlog = _load_json(backlog_path)
    _validate_backlog(backlog)
    scenarios = manifest["scenarios"]
    cases = profile["campaigns"]
    case_counts = campaign_profiles.campaign_counts(profile)
    scenario_rows = []
    for scenario in scenarios:
        design = scenario["design"]
        benchmark_id = design.get("benchmark_id")
        campaigns = sorted(
            case["id"] for case in cases
            if (
                scenario["id"] in case["ids"]
                or scenario["id"] == case["id"]
                or benchmark_id in case["ids"]
            )
        )
        scenario_rows.append({
            "id": scenario["id"],
            "kind": scenario["kind"],
            "panels": sorted(scenario["panels"]),
            "input_format": design.get("input_format"),
            "truth_policy": design.get("truth_policy"),
            "campaigns": campaigns,
        })
    release_cases = [
        case for case in cases
        if case["baseline_policy"] == "paired_required"
    ]
    document = {
        "schema_version": 1,
        "profile_id": "canonical-v2",
        "manifest_sha256": manifest_tools.canonical_sha256(manifest),
        "profile_digest": profile["digest"],
        "campaign_summary": case_counts,
        "release_matrix": {
            "subset_ids": sum(
                len(case["ids"]) for case in release_cases
                if case["stage"] == "paired-subset"
                and case["truth_policy"] != "synthetic_exact_required"
            ),
            "full_ids": sum(
                len(case["ids"]) for case in release_cases
                if case["stage"] == "paired-full"
            ),
            "e2e_ids": sum(
                len(case["ids"]) for case in release_cases
                if case["stage"] == "paired-e2e"
            ),
            "synthetic_ids": sum(
                len(case["ids"]) for case in release_cases
                if case["truth_policy"] == "synthetic_exact_required"
            ),
            "repetitions": sorted({case["repetitions"] for case in release_cases}),
        },
        "dimension_counts": {
            "kind": _count(row["kind"] for row in scenario_rows),
            "panel": _count(
                panel for row in scenario_rows for panel in row["panels"]
            ),
            "input_format": _count(
                row["input_format"] or "not_applicable" for row in scenario_rows
            ),
            "truth_policy": _count(
                row["truth_policy"] or "not_applicable" for row in scenario_rows
            ),
        },
        "scenarios": scenario_rows,
        "future_backlog": backlog,
    }
    document["coverage_sha256"] = manifest_tools.canonical_sha256(document)
    return document


def render_json(document: Mapping[str, Any]) -> str:
    return json.dumps(document, indent=2, sort_keys=True) + "\n"


def render_markdown(document: Mapping[str, Any]) -> str:
    summary = document["campaign_summary"]
    release = document["release_matrix"]
    lines = [
        "# Canonical-v2 Benchmark Coverage",
        "",
        (
            f"Canonical-v2 contains **{summary['campaigns']} campaigns**, "
            f"**{summary['cells']} cells**, and **{summary['arms']} arms**. "
            f"The release matrix has {release['subset_ids']} subset, "
            f"{release['full_ids']} full-genome, and {release['e2e_ids']} "
            f"end-to-end dataset IDs, plus {release['synthetic_ids']} exact "
            "synthetic cases."
        ),
        "",
        "| Scenario | Kind | Panels | Format | Campaigns |",
        "|---|---|---|---|---|",
    ]
    for row in document["scenarios"]:
        lines.append(
            f"| `{row['id']}` | {row['kind']} | "
            f"{', '.join(row['panels'])} | "
            f"{row['input_format'] or 'n/a'} | "
            f"{', '.join(row['campaigns']) or 'manifest-only'} |"
        )
    lines.extend(["", "## Prioritized expansion backlog", ""])
    for item in document["future_backlog"]["items"]:
        lines.append(
            f"- **{item['priority']} — {item['dimension']} / "
            f"`{item['id']}`:** {item['description']}"
        )
    lines.extend([
        "",
        "Backlog cases remain diagnostic until every promotion requirement in "
        "`canonical_v3_backlog.json` is satisfied.",
        "",
    ])
    return "\n".join(lines)


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--profile-registry", type=Path, default=DEFAULT_PROFILES)
    parser.add_argument("--backlog", type=Path, default=DEFAULT_BACKLOG)
    parser.add_argument("--json-output", type=Path, default=DEFAULT_JSON)
    parser.add_argument("--markdown-output", type=Path, default=DEFAULT_MARKDOWN)
    args = parser.parse_args(argv)
    document = build_coverage(
        args.manifest, args.profile_registry, args.backlog,
    )
    outputs = {
        args.json_output: render_json(document),
        args.markdown_output: render_markdown(document),
    }
    if args.check:
        stale = [
            str(path) for path, content in outputs.items()
            if not path.is_file() or path.read_text(encoding="utf-8") != content
        ]
        if stale:
            parser.error("coverage outputs are stale: " + ", ".join(stale))
        return 0
    for path, content in outputs.items():
        path.write_text(content, encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
