#!/usr/bin/env python3
"""Verify exact-v1.0.11 fast/safe semantic equivalence canaries."""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
from typing import Any, Mapping, Sequence

from . import release_evaluation, whole_genome_study


SCHEMA_VERSION = 1
METHOD = "lifton-v1.0.11-fast-safe-equivalence-v1"
EXPECTED_SHA = "c623f0bddc5b5051a3670a3b0064c52cc61bb719"
EXPECTED_IDS = ("bee", "t2_human_to_gorilla")
SUMMARY_FIELDS = (
    "n_reference_total",
    "n_reference_coding",
    "n_tool_features",
    "n_unmapped_id",
    "n_recovered_any",
    "n_recovered_coding",
    "completeness_all",
    "completeness_coding",
    "completeness_feature_total",
    "n_extra_copies",
    "protein_identity",
    "dna_identity",
    "structural",
    "stable_id_preservation",
)
TARGET_FIELDS = (
    "parameters",
    "scope",
    "gene",
    "transcript",
    "structure",
    "target_sequence",
)


class CanaryError(RuntimeError):
    """Raised when fast/safe evidence is incomplete or inequivalent."""


def _pair_paths(root: str | Path) -> list[Path]:
    paths = sorted(Path(root).resolve().rglob("pair_result.json"))
    if len(paths) != len(EXPECTED_IDS):
        raise CanaryError(
            f"expected {len(EXPECTED_IDS)} canary pairs, observed {len(paths)}"
        )
    return paths


def _selected(document: Mapping[str, Any], fields: Sequence[str]) -> dict[str, Any]:
    return {field: document.get(field) for field in fields}


def audit_canary(root: str | Path) -> dict[str, Any]:
    rows = []
    observed_ids = set()
    for path in _pair_paths(root):
        pair = json.loads(path.read_text(encoding="utf-8"))
        benchmark = pair.get("benchmark")
        if (
            pair.get("schema_version") != release_evaluation.SCHEMA_VERSION
            or pair.get("panel") != "e2e"
            or pair.get("repetition") != 1
            or benchmark not in EXPECTED_IDS
            or benchmark in observed_ids
        ):
            raise CanaryError(f"malformed or duplicate canary result: {path}")
        observed_ids.add(benchmark)
        if pair.get("modes") != {"candidate": "fast", "reference": "safe"}:
            raise CanaryError(f"{benchmark}: canary modes are not fast/safe")
        versions = pair.get("versions")
        provenance = pair.get("provenance")
        if not isinstance(versions, Mapping) or not isinstance(provenance, Mapping):
            raise CanaryError(f"{benchmark}: canary arm evidence is missing")
        if any(
            provenance.get(label, {}).get("sha") != EXPECTED_SHA
            for label in ("candidate", "reference")
        ):
            raise CanaryError(f"{benchmark}: canary is not exact v1.0.11")
        fast = versions["candidate"]
        safe = versions["reference"]
        semantic_equal = (
            fast.get("fingerprints", {}).get("semantic_sha256")
            == safe.get("fingerprints", {}).get("semantic_sha256")
            and fast.get("fingerprints", {}).get("feature_counts")
            == safe.get("fingerprints", {}).get("feature_counts")
        )
        summary_equal = _selected(
            fast.get("summary", {}), SUMMARY_FIELDS,
        ) == _selected(safe.get("summary", {}), SUMMARY_FIELDS)
        fast_truth = fast.get("summary", {}).get("target_truth", {})
        safe_truth = safe.get("summary", {}).get("target_truth", {})
        truth_equal = _selected(fast_truth, TARGET_FIELDS) == _selected(
            safe_truth, TARGET_FIELDS,
        )
        validation_equal = {
            "is_valid": fast.get("validation", {}).get("is_valid"),
            "issue_totals": {
                name: value.get("issue_totals")
                for name, value in fast.get("validation", {}).get(
                    "passes", {}
                ).items()
            },
        } == {
            "is_valid": safe.get("validation", {}).get("is_valid"),
            "issue_totals": {
                name: value.get("issue_totals")
                for name, value in safe.get("validation", {}).get(
                    "passes", {}
                ).items()
            },
        }
        equivalent = all((
            semantic_equal,
            summary_equal,
            truth_equal,
            validation_equal,
        ))
        rows.append({
            "benchmark": benchmark,
            "equivalent": equivalent,
            "semantic_equal": semantic_equal,
            "summary_equal": summary_equal,
            "target_truth_equal": truth_equal,
            "validation_equal": validation_equal,
            "fast_byte_sha256": fast["fingerprints"]["byte_sha256"],
            "safe_byte_sha256": safe["fingerprints"]["byte_sha256"],
            "pair_result": whole_genome_study.file_record(path),
        })
    rows.sort(key=lambda row: EXPECTED_IDS.index(row["benchmark"]))
    document = {
        "schema_version": SCHEMA_VERSION,
        "method": METHOD,
        "exact_sha": EXPECTED_SHA,
        "fast_supported": all(row["equivalent"] for row in rows),
        "rows": rows,
    }
    document["fingerprint"] = whole_genome_study.canonical_sha256(document)
    return document


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", type=Path)
    parser.add_argument("--output", type=Path)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    arguments = build_parser().parse_args(argv)
    try:
        result = audit_canary(arguments.root)
        if arguments.output:
            whole_genome_study._atomic_json(arguments.output, result)
    except (CanaryError, OSError, ValueError) as exc:
        print(f"biology-canary: {exc}", file=os.sys.stderr)
        return 2
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0 if result["fast_supported"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
