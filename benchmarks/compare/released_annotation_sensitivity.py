#!/usr/bin/env python3
"""Rescore exact-release outputs against released target annotations.

The primary campaign uses reciprocal-best-hit protein groups.  This reducer
adds two coordinate-overlap sensitivity thresholds for that fixed scope and,
where NCBI Gene supplies enough one-to-one orthologs, a separately labelled
provider-relationship gene analysis.  It never treats either annotation as
independent target-locus truth.
"""
from __future__ import annotations

import argparse
import concurrent.futures
import hashlib
import json
import os
from pathlib import Path
from typing import Any, Mapping, Sequence

from . import (
    provider_orthologs,
    release_evaluation,
    target_truth,
    whole_genome_study,
)


SCHEMA_VERSION = 1
METHOD = "released-target-annotation-sensitivity-v1"
CANDIDATE_SHA = "c623f0bddc5b5051a3670a3b0064c52cc61bb719"
REFERENCE_SHA = "e503643d8346c600fedabcd3a4dff5c0873a4a37"
PRIMARY_THRESHOLD = 0.5
SENSITIVITY_THRESHOLDS = (0.25, 0.8)


class SensitivityError(RuntimeError):
    """Raised when input evidence is incomplete or inconsistent."""


def _fingerprint(document: Mapping[str, Any], label: str) -> None:
    material = dict(document)
    observed = material.pop("fingerprint", None)
    if observed != whole_genome_study.canonical_sha256(material):
        raise SensitivityError(f"{label} fingerprint is invalid")


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _resolved_object(root: Path, record: Any, label: str) -> Path:
    if not isinstance(record, Mapping):
        raise SensitivityError(f"{label} object record is missing")
    try:
        return whole_genome_study._resolved_object(root, record)
    except (OSError, ValueError, whole_genome_study.StudyError) as exc:
        raise SensitivityError(f"{label} object record changed: {exc}") from exc


def _load_inputs(
    study_path: Path,
    preflight_path: Path,
    provider_lock_path: Path,
) -> tuple[dict[str, Any], dict[str, Any], dict[str, Any]]:
    study = whole_genome_study.load_study(study_path)
    preflight = json.loads(preflight_path.read_text(encoding="utf-8"))
    provider = json.loads(provider_lock_path.read_text(encoding="utf-8"))
    if (
        preflight.get("kind")
        != "lifton-v1.0.11-biology-study-preflight"
        or not preflight.get("campaign_ready")
        or preflight.get("study") != whole_genome_study.file_record(study_path)
    ):
        raise SensitivityError("preflight is not bound to the study")
    _fingerprint(preflight, "preflight")
    if (
        provider.get("kind") != "lifton-v1.0.11-provider-ortholog-lock"
        or provider.get("method", provider_orthologs.METHOD)
        != provider_orthologs.METHOD
        or provider.get("study") != whole_genome_study.file_record(study_path)
        or provider.get("preflight")
        != whole_genome_study.file_record(preflight_path)
        or set(provider.get("maps", {}))
        != set(whole_genome_study.EXPECTED_PAIR_IDS)
    ):
        raise SensitivityError("provider lock is not bound to the study")
    _fingerprint(provider, "provider lock")
    return study, preflight, provider


def _repetition_one_pairs(campaign_root: Path) -> dict[str, tuple[Path, dict]]:
    result = {}
    for path in sorted(campaign_root.resolve().rglob("pair_result.json")):
        document = json.loads(path.read_text(encoding="utf-8"))
        if document.get("repetition") != 1:
            continue
        pair_id = document.get("benchmark")
        if (
            document.get("schema_version") != release_evaluation.SCHEMA_VERSION
            or document.get("panel") != "e2e"
            or pair_id not in whole_genome_study.EXPECTED_PAIR_IDS
            or pair_id in result
            or document.get("modes")
            != {"candidate": "fast", "reference": "fast"}
            or document.get("provenance", {}).get("candidate", {}).get("sha")
            != CANDIDATE_SHA
            or document.get("provenance", {}).get("reference", {}).get("sha")
            != REFERENCE_SHA
        ):
            raise SensitivityError(f"malformed repetition-one pair: {path}")
        result[pair_id] = (path, document)
    if set(result) != set(whole_genome_study.EXPECTED_PAIR_IDS):
        raise SensitivityError("campaign lacks one repetition-one result per pair")
    return result


def _output_path(version: Mapping[str, Any], label: str) -> Path:
    path = Path(str(version.get("output_gff", ""))).resolve()
    fingerprints = version.get("fingerprints")
    if (
        not path.is_file()
        or not isinstance(fingerprints, Mapping)
        or fingerprints.get("byte_sha256") != _sha256(path)
    ):
        raise SensitivityError(f"{label}: output GFF changed after execution")
    return path


def _project(document: Mapping[str, Any]) -> dict[str, Any]:
    return {
        "parameters": document["parameters"],
        "scope": document["scope"],
        "gene": document["gene"],
        "transcript": document["transcript"],
        "structure": document["structure"],
    }


def _score(
    prediction: Path,
    truth: Path,
    source: Path,
    mapping: Path,
    threshold: float,
) -> dict[str, Any]:
    return _project(target_truth.score_target_truth(
        prediction,
        truth,
        ortholog_map=mapping,
        source_gff=source,
        id_policy="ortholog-map",
        minimum_reciprocal_overlap=threshold,
    ))


def _primary_from_pair(version: Mapping[str, Any], label: str) -> dict[str, Any]:
    document = version.get("summary", {}).get("target_truth")
    if (
        not isinstance(document, Mapping)
        or document.get("method") != target_truth.METHOD
        or document.get("parameters", {}).get("minimum_reciprocal_overlap")
        != PRIMARY_THRESHOLD
    ):
        raise SensitivityError(f"{label}: primary concordance evidence is missing")
    record = version.get("evaluation_artifacts", {}).get("target_truth")
    if not isinstance(record, Mapping):
        raise SensitivityError(f"{label}: primary concordance record is missing")
    artifact = Path(str(record.get("path", ""))).resolve()
    if (
        whole_genome_study.file_record(artifact).get("sha256")
        != record.get("sha256")
        or json.loads(artifact.read_text(encoding="utf-8")) != document
    ):
        raise SensitivityError(f"{label}: primary concordance artifact changed")
    return _project(document)


def _pair_sensitivity(
    pair: Mapping[str, Any],
    pair_result: tuple[Path, dict],
    preflight: Mapping[str, Any],
    provider: Mapping[str, Any],
    cache_root: Path,
) -> dict[str, Any]:
    pair_id = pair["id"]
    pair_path, document = pair_result
    views = preflight["model_views"][pair_id]
    source = _resolved_object(cache_root, views["source"]["output"], pair_id)
    truth = _resolved_object(cache_root, views["target"]["output"], pair_id)
    primary_map = Path(
        preflight["ortholog_scopes"][pair_id]["mapping"]["path"]
    ).resolve()
    if whole_genome_study.file_record(primary_map) != (
        preflight["ortholog_scopes"][pair_id]["mapping"]
    ):
        raise SensitivityError(f"{pair_id}: primary ortholog map changed")
    versions = document["versions"]
    predictions = {
        label: _output_path(versions[label], f"{pair_id}.{label}")
        for label in ("candidate", "reference")
    }
    primary = {
        str(PRIMARY_THRESHOLD): {
            label: _primary_from_pair(
                versions[label], f"{pair_id}.{label}",
            )
            for label in ("candidate", "reference")
        }
    }
    for threshold in SENSITIVITY_THRESHOLDS:
        primary[str(threshold)] = {
            label: _score(
                predictions[label], truth, source, primary_map, threshold,
            )
            for label in ("candidate", "reference")
        }
    provider_row = provider["maps"][pair_id]
    provider_result: dict[str, Any] = {
        "available": bool(provider_row.get("available")),
        "reason": provider_row.get("reason"),
        "groups": int(provider_row.get("groups", 0)),
    }
    if provider_result["available"]:
        mapping = _resolved_object(
            cache_root, provider_row.get("mapping"), f"{pair_id}.provider",
        )
        provider_result["threshold"] = PRIMARY_THRESHOLD
        provider_result["arms"] = {
            label: _score(
                predictions[label], truth, source, mapping, PRIMARY_THRESHOLD,
            )
            for label in ("candidate", "reference")
        }
    return {
        "id": pair_id,
        "public_label": pair["public_label"],
        "pair_result": whole_genome_study.file_record(pair_path),
        "primary_protein_rbh": primary,
        "provider_gene_relationships": provider_result,
    }


def build_sensitivity(
    study_path: str | Path,
    preflight_path: str | Path,
    provider_lock_path: str | Path,
    campaign_root: str | Path,
    *,
    max_active: int = 2,
) -> dict[str, Any]:
    if max_active < 1:
        raise SensitivityError("max_active must be positive")
    study_path = Path(study_path).resolve()
    preflight_path = Path(preflight_path).resolve()
    provider_lock_path = Path(provider_lock_path).resolve()
    campaign_root = Path(campaign_root).resolve()
    study, preflight, provider = _load_inputs(
        study_path, preflight_path, provider_lock_path,
    )
    pairs = _repetition_one_pairs(campaign_root)
    cache_root = preflight_path.parent
    results = {}
    with concurrent.futures.ThreadPoolExecutor(max_workers=max_active) as pool:
        futures = {
            pool.submit(
                _pair_sensitivity,
                pair,
                pairs[pair["id"]],
                preflight,
                provider,
                cache_root,
            ): pair["id"]
            for pair in study["pairs"]
        }
        for future in concurrent.futures.as_completed(futures):
            pair_id = futures[future]
            try:
                results[pair_id] = future.result()
            except Exception as exc:
                for pending in futures:
                    pending.cancel()
                raise SensitivityError(f"{pair_id}: sensitivity failed: {exc}") from exc
    document = {
        "schema_version": SCHEMA_VERSION,
        "method": METHOD,
        "created_at": whole_genome_study.utc_now(),
        "interpretation": (
            "concordance with released target annotations; not independent "
            "target-locus truth"
        ),
        "study": whole_genome_study.file_record(study_path),
        "preflight": whole_genome_study.file_record(preflight_path),
        "provider_ortholog_lock": whole_genome_study.file_record(
            provider_lock_path
        ),
        "parameters": {
            "primary_scope": "reciprocal-best-hit protein groups",
            "primary_threshold": PRIMARY_THRESHOLD,
            "overlap_sensitivity_thresholds": list(SENSITIVITY_THRESHOLDS),
            "provider_scope": "NCBI Gene one-to-one ortholog relationships",
            "repetition": 1,
        },
        "pairs": [results[pair["id"]] for pair in study["pairs"]],
    }
    material = dict(document)
    material.pop("created_at")
    document["fingerprint"] = whole_genome_study.canonical_sha256(material)
    return document


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--study", type=Path, default=whole_genome_study.DEFAULT_STUDY)
    parser.add_argument("--preflight", type=Path, required=True)
    parser.add_argument("--provider-lock", type=Path, required=True)
    parser.add_argument("--campaign-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--max-active", type=int, default=2)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    arguments = build_parser().parse_args(argv)
    try:
        document = build_sensitivity(
            arguments.study,
            arguments.preflight,
            arguments.provider_lock,
            arguments.campaign_root,
            max_active=arguments.max_active,
        )
        whole_genome_study._atomic_json(arguments.output, document)
    except (OSError, ValueError, SensitivityError) as exc:
        print(f"released-annotation-sensitivity: {exc}", file=os.sys.stderr)
        return 2
    print(json.dumps({
        "pairs": len(document["pairs"]),
        "fingerprint": document["fingerprint"],
    }, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
