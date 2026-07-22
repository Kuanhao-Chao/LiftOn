#!/usr/bin/env python3
"""Build a deterministic inventory of curated LiftOn benchmark artifacts.

Only Git-tracked files and the explicitly declared ``managed_paths`` are
eligible for the committed inventory.  Local run directories are intentionally
ignored: their patterns are documented as ineligible in the rules file, but
their machine-specific contents can never leak into release evidence.

The generated JSON contains a byte fingerprint and a compact JSON-shape
fingerprint for every curated file.  The Markdown view is a concise index of
the same family classifications.
"""
from __future__ import annotations

import argparse
import fnmatch
import hashlib
import json
import subprocess
import sys
from collections import Counter
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence


HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parent
DEFAULT_RULES = HERE / "inventory_rules.json"
DEFAULT_JSON = HERE / "BENCHMARK_INVENTORY.json"
DEFAULT_MARKDOWN = HERE / "BENCHMARK_INVENTORY.md"
OUTPUT_PATHS = {
    "benchmarks/BENCHMARK_INVENTORY.json",
    "benchmarks/BENCHMARK_INVENTORY.md",
}
VALID_CLASSIFICATIONS = {
    "canonical", "diagnostic", "historical", "obsolete", "ineligible",
}
VALID_COMPLETION_STATES = {
    "complete", "incomplete", "superseded", "supporting", "retired",
}
SHA256_LENGTH = 64


class InventoryError(ValueError):
    """Raised when inventory rules or curated artifacts are inconsistent."""


def canonical_json(value: Any) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"),
                      ensure_ascii=True)


def canonical_sha256(value: Any) -> str:
    return hashlib.sha256(canonical_json(value).encode("utf-8")).hexdigest()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _tracked_paths(repo_root: Path) -> list[str]:
    result = subprocess.run(
        ["git", "ls-files", "--", "benchmarks"],
        cwd=repo_root,
        check=False,
        text=True,
        capture_output=True,
    )
    if result.returncode != 0:
        raise InventoryError(
            f"git ls-files failed: {(result.stderr or result.stdout).strip()}"
        )
    return sorted(line for line in result.stdout.splitlines() if line)


def load_rules(path: Path = DEFAULT_RULES) -> dict[str, Any]:
    try:
        rules = json.loads(Path(path).read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise InventoryError(f"cannot load inventory rules {path}: {exc}") from exc
    if not isinstance(rules, dict):
        raise InventoryError("inventory rules must be a JSON object")
    return rules


def validate_rules(rules: Mapping[str, Any]) -> None:
    required = {
        "schema_version", "scope", "managed_paths", "output_paths",
        "ineligible_local_patterns", "frozen_artifacts", "families",
    }
    if set(rules) != required:
        raise InventoryError(
            f"inventory rules keys must be exactly {sorted(required)}"
        )
    if rules["schema_version"] != 1 or rules["scope"] != "benchmarks":
        raise InventoryError("unsupported inventory rules schema or scope")
    if set(rules["output_paths"]) != OUTPUT_PATHS:
        raise InventoryError("rules output_paths do not match inventory outputs")
    if not isinstance(rules["managed_paths"], list):
        raise InventoryError("managed_paths must be a list")
    if not isinstance(rules["ineligible_local_patterns"], list):
        raise InventoryError("ineligible_local_patterns must be a list")
    if not isinstance(rules["frozen_artifacts"], dict):
        raise InventoryError("frozen_artifacts must be an object")

    family_ids: set[str] = set()
    for family in rules["families"]:
        keys = {
            "id", "classification", "claim_eligible", "completion_state",
            "description", "superseded_by", "patterns",
        }
        if not isinstance(family, dict) or set(family) != keys:
            raise InventoryError(
                f"each family must have exactly {sorted(keys)}"
            )
        family_id = family["id"]
        if not isinstance(family_id, str) or not family_id:
            raise InventoryError("family id must be a non-empty string")
        if family_id in family_ids:
            raise InventoryError(f"duplicate family id: {family_id}")
        family_ids.add(family_id)
        if family["classification"] not in VALID_CLASSIFICATIONS:
            raise InventoryError(
                f"{family_id}: invalid classification "
                f"{family['classification']!r}"
            )
        if family["completion_state"] not in VALID_COMPLETION_STATES:
            raise InventoryError(
                f"{family_id}: invalid completion_state "
                f"{family['completion_state']!r}"
            )
        if not isinstance(family["claim_eligible"], bool):
            raise InventoryError(f"{family_id}: claim_eligible must be boolean")
        if not isinstance(family["description"], str) or not family["description"]:
            raise InventoryError(f"{family_id}: description must be non-empty")
        if (
            family["superseded_by"] is not None
            and not isinstance(family["superseded_by"], str)
        ):
            raise InventoryError(f"{family_id}: superseded_by must be string/null")
        if not isinstance(family["patterns"], list) or not family["patterns"]:
            raise InventoryError(f"{family_id}: patterns must be non-empty")

    for family in rules["families"]:
        successor = family["superseded_by"]
        if successor is not None and successor not in family_ids:
            raise InventoryError(
                f"{family['id']}: unknown superseded_by family {successor!r}"
            )

    for path, digest in rules["frozen_artifacts"].items():
        if (
            not isinstance(path, str)
            or not isinstance(digest, str)
            or len(digest) != SHA256_LENGTH
            or any(char not in "0123456789abcdef" for char in digest)
        ):
            raise InventoryError(f"invalid frozen artifact entry: {path!r}")


def _compact_json_shape(value: Any, depth: int = 0) -> Any:
    """Return a bounded structural signature rather than copying JSON values."""

    if depth >= 4:
        return type(value).__name__
    if isinstance(value, dict):
        keys = sorted(str(key) for key in value)
        if len(keys) > 20:
            child_shapes = {
                canonical_sha256(_compact_json_shape(child, depth + 1))
                for child in value.values()
            }
            return {
                "type": "object-map",
                "key_count": len(keys),
                "value_shape_sha256": sorted(child_shapes),
            }
        return {
            "type": "object",
            "fields": {
                str(key): _compact_json_shape(value[key], depth + 1)
                for key in sorted(value, key=str)
            },
        }
    if isinstance(value, list):
        item_shapes = {
            canonical_sha256(_compact_json_shape(item, depth + 1))
            for item in value
        }
        return {
            "type": "array",
            "length": len(value),
            "item_shape_sha256": sorted(item_shapes),
        }
    if value is None:
        return "null"
    if isinstance(value, bool):
        return "boolean"
    if isinstance(value, int):
        return "integer"
    if isinstance(value, float):
        return "number"
    if isinstance(value, str):
        return "string"
    return type(value).__name__


def _json_metadata(path: Path) -> dict[str, Any] | None:
    if path.suffix.lower() != ".json":
        return None
    try:
        document = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise InventoryError(f"invalid curated JSON {path}: {exc}") from exc
    shape = _compact_json_shape(document)
    return {
        "document_type": (
            "object" if isinstance(document, dict)
            else "array" if isinstance(document, list)
            else type(document).__name__
        ),
        "top_level_count": (
            len(document) if isinstance(document, (dict, list)) else None
        ),
        "shape_sha256": canonical_sha256(shape),
    }


def _family_for_path(
    path: str,
    families: Sequence[Mapping[str, Any]],
) -> Mapping[str, Any]:
    matches = [
        family for family in families
        if any(fnmatch.fnmatchcase(path, pattern) for pattern in family["patterns"])
    ]
    if not matches:
        raise InventoryError(f"curated path is unclassified: {path}")
    if len(matches) > 1:
        raise InventoryError(
            f"curated path matches multiple families: {path}: "
            f"{[family['id'] for family in matches]}"
        )
    return matches[0]


def _registry_summary(repo_root: Path) -> dict[str, Any]:
    registry_path = repo_root / "benchmarks" / "compare" / "benchmarks.json"
    registry = json.loads(registry_path.read_text(encoding="utf-8"))
    entries = registry.get("benchmarks")
    if not isinstance(entries, list):
        raise InventoryError(f"{registry_path}: benchmarks must be a list")
    required = {"id", "ref_genome", "ref_gff", "tgt_genome"}
    ids = [entry.get("id") for entry in entries if isinstance(entry, dict)]
    duplicates = sorted(
        item for item, count in Counter(ids).items()
        if item is not None and count > 1
    )
    incomplete = sorted(
        str(entry.get("id"))
        for entry in entries
        if not isinstance(entry, dict)
        or any(not entry.get(field) for field in required)
    )
    return {
        "path": "benchmarks/compare/benchmarks.json",
        "entry_count": len(entries),
        "ids": ids,
        "duplicate_ids": duplicates,
        "schema_incomplete_ids": incomplete,
    }


def _canonical_v2_summary(repo_root: Path) -> dict[str, Any]:
    from benchmarks.manifest_tools import validate_manifest

    path = repo_root / "benchmarks" / "manifests" / "canonical_v2_datasets.json"
    try:
        document = validate_manifest(json.loads(path.read_text(encoding="utf-8")))
    except (OSError, json.JSONDecodeError, ValueError) as exc:
        raise InventoryError(f"invalid canonical-v2 dataset manifest: {exc}") from exc
    return {
        "path": "benchmarks/manifests/canonical_v2_datasets.json",
        "sha256": sha256_file(path),
        "source_count": len(document["sources"]),
        "remote_source_count": sum(
            source["transport"] != "existing_registry"
            for source in document["sources"]
        ),
        "scenario_count": len(document["scenarios"]),
        "scenarios": [
            {
                "id": scenario["id"],
                "kind": scenario["kind"],
                "panels": scenario["panels"],
            }
            for scenario in document["scenarios"]
        ],
    }


def build_inventory(
    repo_root: Path,
    rules: Mapping[str, Any],
    *,
    tracked_paths: Iterable[str] | None = None,
) -> dict[str, Any]:
    validate_rules(rules)
    repo_root = Path(repo_root).resolve()
    tracked = set(
        tracked_paths if tracked_paths is not None else _tracked_paths(repo_root)
    )
    tracked.update(rules["managed_paths"])
    tracked.difference_update(rules["output_paths"])
    curated_paths = sorted(
        path for path in tracked
        if path == "benchmarks" or path.startswith("benchmarks/")
    )

    missing = [path for path in curated_paths if not (repo_root / path).is_file()]
    if missing:
        raise InventoryError(f"curated paths are missing: {missing}")

    records: list[dict[str, Any]] = []
    for relative in curated_paths:
        absolute = repo_root / relative
        family = _family_for_path(relative, rules["families"])
        record = {
            "path": relative,
            "bytes": absolute.stat().st_size,
            "sha256": sha256_file(absolute),
            "family": family["id"],
            "classification": family["classification"],
            "claim_eligible": family["claim_eligible"],
            "completion_state": family["completion_state"],
            "json_schema": _json_metadata(absolute),
        }
        records.append(record)

    frozen: list[dict[str, Any]] = []
    by_path = {record["path"]: record for record in records}
    for path, expected in sorted(rules["frozen_artifacts"].items()):
        record = by_path.get(path)
        if record is None:
            raise InventoryError(f"frozen artifact is not curated: {path}")
        observed = record["sha256"]
        if observed != expected:
            raise InventoryError(
                f"frozen artifact changed: {path}: expected {expected}, "
                f"observed {observed}"
            )
        frozen.append({"path": path, "sha256": observed, "verified": True})

    family_rows = []
    for family in rules["families"]:
        members = [
            record for record in records if record["family"] == family["id"]
        ]
        if not members:
            raise InventoryError(
                f"family has no curated members: {family['id']}"
            )
        family_rows.append({
            key: family[key]
            for key in (
                "id", "classification", "claim_eligible", "completion_state",
                "description", "superseded_by",
            )
        } | {
            "file_count": len(members),
            "bytes": sum(record["bytes"] for record in members),
            "aggregate_sha256": canonical_sha256([
                {
                    "path": record["path"],
                    "bytes": record["bytes"],
                    "sha256": record["sha256"],
                }
                for record in members
            ]),
        })

    classification_counts = Counter(
        record["classification"] for record in records
    )
    inventory = {
        "schema_version": 1,
        "scope": "git-tracked benchmark artifacts plus declared managed paths",
        "determinism": {
            "timestamps_omitted": True,
            "local_run_contents_omitted": True,
            "ordering": "UTF-8 repository-relative path",
            "digest": "sha256",
        },
        "rules_sha256": canonical_sha256(rules),
        "source_summary": {
            "file_count": len(records),
            "bytes": sum(record["bytes"] for record in records),
            "aggregate_sha256": canonical_sha256([
                {
                    "path": record["path"],
                    "bytes": record["bytes"],
                    "sha256": record["sha256"],
                }
                for record in records
            ]),
            "classification_counts": {
                name: classification_counts.get(name, 0)
                for name in sorted(VALID_CLASSIFICATIONS)
            },
        },
        "registry_summary": _registry_summary(repo_root),
        "canonical_v2_dataset_summary": _canonical_v2_summary(repo_root),
        "frozen_artifacts": frozen,
        "families": family_rows,
        "files": records,
        "ineligible_local_patterns": rules["ineligible_local_patterns"],
        "generated_outputs_excluded_from_self_hash": sorted(OUTPUT_PATHS),
    }
    inventory["inventory_sha256"] = canonical_sha256(inventory)
    return inventory


def render_json(inventory: Mapping[str, Any]) -> str:
    return json.dumps(inventory, indent=2, sort_keys=True) + "\n"


def render_markdown(inventory: Mapping[str, Any]) -> str:
    source = inventory["source_summary"]
    registry = inventory["registry_summary"]
    canonical_v2 = inventory["canonical_v2_dataset_summary"]
    lines = [
        "# Curated Benchmark Inventory",
        "",
        "This inventory is generated deterministically from Git-tracked benchmark "
        "files and explicitly managed provenance files. Local run trees are "
        "excluded from release evidence.",
        "",
        f"- Curated files: **{source['file_count']}** "
        f"({source['bytes']:,} bytes)",
        f"- Registered benchmark IDs: **{registry['entry_count']}**",
        f"- Frozen artifacts verified: **{len(inventory['frozen_artifacts'])}**",
        f"- Inventory digest: `{inventory['inventory_sha256']}`",
        "",
        "## Classifications",
        "",
        "| Family | Class | State | Claim eligible | Files | Superseded by |",
        "| --- | --- | --- | --- | ---: | --- |",
    ]
    for family in inventory["families"]:
        lines.append(
            f"| `{family['id']}` | {family['classification']} | "
            f"{family['completion_state']} | "
            f"{'yes' if family['claim_eligible'] else 'no'} | "
            f"{family['file_count']} | "
            f"{family['superseded_by'] or '—'} |"
        )
    lines.extend([
        "",
        "## Canonical-v2 Expansion",
        "",
        f"The expansion fixes **{canonical_v2['source_count']}** source packages "
        f"({canonical_v2['remote_source_count']} remote and "
        f"{canonical_v2['source_count'] - canonical_v2['remote_source_count']} "
        "existing-registry packages) across "
        f"**{canonical_v2['scenario_count']}** approved scenario families.",
        "",
        "| Scenario | Kind | Panels |",
        "| --- | --- | --- |",
    ])
    for scenario in canonical_v2["scenarios"]:
        lines.append(
            f"| `{scenario['id']}` | {scenario['kind']} | "
            f"{', '.join(scenario['panels'])} |"
        )
    lines.extend([
        "",
        "## Evidence Policy",
        "",
        "- `canonical` files may support a claim only when the family is also "
        "marked claim-eligible.",
        "- `diagnostic` and `historical` results explain design decisions but "
        "must not be pooled into the v1.0.10.1 release campaign.",
        "- `obsolete` files are retained solely for reproducibility of older "
        "reviews.",
        "- `ineligible` files are incomplete, mutable, local, or otherwise "
        "unsuitable as publication evidence.",
        "- `benchmarks/compare/fourway_results.json` is byte-frozen; the scanner "
        "fails if its SHA-256 changes.",
        "",
        "## Regeneration",
        "",
        "```bash",
        "python -m benchmarks.inventory --check",
        "python -m benchmarks.inventory --stdout json > /tmp/inventory.json",
        "python -m benchmarks.inventory --stdout markdown > /tmp/inventory.md",
        "```",
        "",
        "Review and apply regenerated files deliberately; do not mix local "
        "`_runs/`, `work/`, rerun figures, or untracked outputs into this index.",
        "",
    ])
    return "\n".join(lines)


def _check_output(path: Path, expected: str) -> None:
    try:
        observed = path.read_text(encoding="utf-8")
    except OSError as exc:
        raise InventoryError(f"cannot read generated output {path}: {exc}") from exc
    if observed != expected:
        raise InventoryError(f"generated inventory is stale: {path}")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rules", type=Path, default=DEFAULT_RULES)
    parser.add_argument("--json-output", type=Path, default=DEFAULT_JSON)
    parser.add_argument("--markdown-output", type=Path, default=DEFAULT_MARKDOWN)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--check", action="store_true",
                      help="verify committed outputs without writing")
    mode.add_argument("--stdout", choices=("json", "markdown"),
                      help="print one generated representation without writing")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        rules = load_rules(args.rules)
        inventory = build_inventory(REPO_ROOT, rules)
        json_text = render_json(inventory)
        markdown_text = render_markdown(inventory)
        if args.check:
            _check_output(args.json_output, json_text)
            _check_output(args.markdown_output, markdown_text)
            print(
                f"benchmark inventory verified: "
                f"{inventory['source_summary']['file_count']} files "
                f"{inventory['inventory_sha256']}"
            )
        elif args.stdout:
            sys.stdout.write(
                json_text if args.stdout == "json" else markdown_text
            )
        else:
            args.json_output.write_text(json_text, encoding="utf-8")
            args.markdown_output.write_text(markdown_text, encoding="utf-8")
            print(
                f"wrote {args.json_output} and {args.markdown_output}: "
                f"{inventory['source_summary']['file_count']} files"
            )
    except InventoryError as exc:
        print(f"inventory error: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
