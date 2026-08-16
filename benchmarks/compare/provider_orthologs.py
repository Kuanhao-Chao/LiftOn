#!/usr/bin/env python3
"""Build separately labeled NCBI Gene ortholog sensitivity scopes."""
from __future__ import annotations

import argparse
import gzip
import json
import os
import urllib.parse
from collections import defaultdict
from pathlib import Path
from typing import Any, Mapping, Sequence

from lifton.gff3_validator import GENE_TYPES

from . import whole_genome_study


SCHEMA_VERSION = 1
METHOD = "ncbi-gene-ortholog-sensitivity-v1"
LOCK_NAME = "provider-ortholog-lock.json"
SOURCE_URL = "https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_orthologs.gz"
MINIMUM_GROUPS = 100


class ProviderOrthologError(RuntimeError):
    """Raised when provider scope evidence is malformed or inconsistent."""


def _attributes(text: str) -> dict[str, tuple[str, ...]]:
    result = {}
    for item in text.split(";"):
        if not item or "=" not in item:
            continue
        name, value = item.split("=", 1)
        result[name] = tuple(
            urllib.parse.unquote(part, errors="strict")
            for part in value.split(",")
            if part
        )
    return result


def gene_id_index(path: str | Path) -> dict[str, str]:
    """Return unambiguous NCBI GeneID-to-GFF3-gene-ID mappings."""

    annotation = Path(path).resolve()
    by_gene_id: dict[str, set[str]] = defaultdict(set)
    by_model_id: dict[str, set[str]] = defaultdict(set)
    with annotation.open("r", encoding="utf-8", errors="strict") as handle:
        for line_number, raw in enumerate(handle, start=1):
            if not raw.strip() or raw.startswith("#"):
                continue
            columns = raw.rstrip("\r\n").split("\t")
            if len(columns) != 9:
                raise ProviderOrthologError(
                    f"{annotation}: line {line_number} has {len(columns)} columns"
                )
            if columns[2] not in GENE_TYPES:
                continue
            attributes = _attributes(columns[8])
            model_ids = attributes.get("ID", ())
            if len(model_ids) != 1:
                raise ProviderOrthologError(
                    f"{annotation}: line {line_number} gene lacks one ID"
                )
            gene_ids = {
                value.split(":", 1)[1]
                for value in attributes.get("Dbxref", ())
                if value.startswith("GeneID:")
                and value.split(":", 1)[1].isdigit()
            }
            for gene_id in gene_ids:
                by_gene_id[gene_id].add(model_ids[0])
                by_model_id[model_ids[0]].add(gene_id)
    return {
        gene_id: next(iter(model_ids))
        for gene_id, model_ids in by_gene_id.items()
        if len(model_ids) == 1
        and len(by_model_id[next(iter(model_ids))]) == 1
    }


def relevant_relationships(
    path: str | Path,
    taxon_pairs: Mapping[str, tuple[int, int]],
) -> dict[str, dict[str, set[str]]]:
    requested = {
        (source_tax, target_tax): pair_id
        for pair_id, (source_tax, target_tax) in taxon_pairs.items()
        if source_tax != target_tax
    }
    relationships: dict[str, dict[str, set[str]]] = {
        pair_id: defaultdict(set) for pair_id in requested.values()
    }
    with gzip.open(path, "rt", encoding="utf-8", errors="strict") as handle:
        header = handle.readline().rstrip("\r\n").lstrip("#").split("\t")
        if header != [
            "tax_id", "GeneID", "relationship", "Other_tax_id",
            "Other_GeneID",
        ]:
            raise ProviderOrthologError("NCBI ortholog table header is unsupported")
        for line_number, raw in enumerate(handle, start=2):
            fields = raw.rstrip("\r\n").split("\t")
            if len(fields) != 5:
                raise ProviderOrthologError(
                    f"NCBI ortholog table line {line_number} is malformed"
                )
            tax_id, gene_id, relationship, other_tax_id, other_gene_id = fields
            if relationship != "Ortholog":
                continue
            forward = requested.get((int(tax_id), int(other_tax_id)))
            reverse = requested.get((int(other_tax_id), int(tax_id)))
            if forward is not None:
                relationships[forward][gene_id].add(other_gene_id)
            if reverse is not None:
                relationships[reverse][other_gene_id].add(gene_id)
    return relationships


def build_mapping(
    pair_id: str,
    relationships: Mapping[str, set[str]],
    source_index: Mapping[str, str],
    target_index: Mapping[str, str],
) -> dict[str, Any]:
    candidates = {
        source_gene: next(iter(target_genes))
        for source_gene, target_genes in relationships.items()
        if source_gene in source_index
        and len(target_genes) == 1
        and next(iter(target_genes)) in target_index
    }
    reverse: dict[str, list[str]] = defaultdict(list)
    for source_gene, target_gene in candidates.items():
        reverse[target_gene].append(source_gene)
    retained = [
        (source_gene, target_gene)
        for source_gene, target_gene in candidates.items()
        if len(reverse[target_gene]) == 1
    ]
    mappings = [
        {
            "source_id": source_index[source_gene],
            "truth_ids": [target_index[target_gene]],
            "feature_type": "gene",
            "status": "retained",
            "evidence": {
                "source_gene_id": source_gene,
                "target_gene_id": target_gene,
                "relationship": "Ortholog",
            },
        }
        for source_gene, target_gene in sorted(retained)
    ]
    return {
        "schema_version": SCHEMA_VERSION,
        "method": METHOD,
        "metadata": {
            "pair_id": pair_id,
            "scope": "one-to-one provider-curated genes",
            "groups": len(mappings),
            "source_annotation_gene_ids": len(source_index),
            "target_annotation_gene_ids": len(target_index),
        },
        "mappings": mappings,
    }


def _preflight(root: Path) -> dict[str, Any]:
    path = root / whole_genome_study.PREFLIGHT_NAME
    try:
        document = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise ProviderOrthologError(f"cannot read study preflight: {exc}") from exc
    material = dict(document)
    fingerprint = material.pop("fingerprint", None)
    if (
        document.get("kind")
        != "lifton-v1.0.11-biology-study-preflight"
        or not document.get("campaign_ready")
        or fingerprint != whole_genome_study.canonical_sha256(material)
    ):
        raise ProviderOrthologError("study preflight is not valid and complete")
    return document


def build_provider_scopes(
    study_path: str | Path,
    cache_root: str | Path,
    table_path: str | Path,
) -> dict[str, Any]:
    study = whole_genome_study.load_study(study_path)
    root = whole_genome_study._safe_cache_root(cache_root)
    lock_path = root / LOCK_NAME
    if lock_path.exists():
        return json.loads(lock_path.read_text(encoding="utf-8"))
    preflight = _preflight(root)
    table = Path(table_path).resolve()
    table_record = whole_genome_study._publish_file(
        root, table, "ncbi-gene-orthologs.gz",
    )
    table_object = whole_genome_study._resolved_object(root, table_record)
    taxon_pairs = {
        pair["id"]: (pair["source_tax_id"], pair["target_tax_id"])
        for pair in study["pairs"]
    }
    relationships = relevant_relationships(table_object, taxon_pairs)
    maps = {}
    for pair in study["pairs"]:
        pair_id = pair["id"]
        if pair["source_tax_id"] == pair["target_tax_id"]:
            maps[pair_id] = {
                "available": False,
                "reason": "same_taxon_not_applicable",
                "groups": 0,
            }
            continue
        views = preflight["model_views"][pair_id]
        source = whole_genome_study._resolved_object(
            root, views["source"]["output"],
        )
        target = whole_genome_study._resolved_object(
            root, views["target"]["output"],
        )
        mapping = build_mapping(
            pair_id,
            relationships[pair_id],
            gene_id_index(source),
            gene_id_index(target),
        )
        temporary = root / "staging" / f"{pair_id}-provider-map.json"
        temporary.unlink(missing_ok=True)
        whole_genome_study._atomic_json(temporary, mapping)
        record = whole_genome_study._publish_file(
            root, temporary, f"{pair_id}-provider-map.json",
        )
        groups = mapping["metadata"]["groups"]
        maps[pair_id] = {
            "available": groups >= MINIMUM_GROUPS,
            "reason": (
                None if groups >= MINIMUM_GROUPS
                else "insufficient_provider_coverage"
            ),
            "groups": groups,
            "mapping": record,
        }
    lock = {
        "schema_version": SCHEMA_VERSION,
        "kind": "lifton-v1.0.11-provider-ortholog-lock",
        "created_at": whole_genome_study.utc_now(),
        "study": whole_genome_study.file_record(study_path),
        "preflight": whole_genome_study.file_record(
            root / whole_genome_study.PREFLIGHT_NAME
        ),
        "source": {
            "url": SOURCE_URL,
            "retrieved_at": "2026-08-13T09:44:43Z",
            "upstream_last_modified": "2026-08-12T07:35:03Z",
            "artifact": table_record,
        },
        "minimum_groups": MINIMUM_GROUPS,
        "maps": maps,
    }
    lock["fingerprint"] = whole_genome_study.canonical_sha256(lock)
    whole_genome_study._atomic_json(lock_path, lock)
    return lock


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--study", type=Path, default=whole_genome_study.DEFAULT_STUDY)
    parser.add_argument("--cache-root", type=Path, required=True)
    parser.add_argument("--table", type=Path, required=True)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    arguments = build_parser().parse_args(argv)
    try:
        result = build_provider_scopes(
            arguments.study, arguments.cache_root, arguments.table,
        )
    except (OSError, ValueError, ProviderOrthologError) as exc:
        print(f"provider-orthologs: {exc}", file=os.sys.stderr)
        return 2
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
