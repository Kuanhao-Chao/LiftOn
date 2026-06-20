#!/usr/bin/env python
"""Register the tiered pairs (``benchmarks/tiers.json``) into the 4-way
comparison registry (``benchmarks/compare/benchmarks.json``) so each becomes a
live cell scored by Liftoff / miniprot / LiftOn-stable / LiftOn-devel.

Idempotent: an entry whose ``id`` already exists is replaced (not duplicated);
all other existing entries keep their order and content. Re-serialises with the
registry's 2-space indent — existing entries are unchanged apart from the few
cosmetic blank separator lines that JSON re-serialisation drops; the diff is
otherwise purely the appended tier entries.

Usage (repo root)::

    python -m benchmarks.register_tiers              # merge + write
    python -m benchmarks.register_tiers --dry-run    # print summary, write nothing
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

REGISTRY = Path(__file__).resolve().parent / "tiers.json"
BENCHMARKS = Path(__file__).resolve().parent / "compare" / "benchmarks.json"

# Standard benchmarks.json key order for a new pair (mirrors the existing
# divergence-ladder entries so the registry stays visually consistent).
def build_entry(pair: dict, data_root: str) -> dict:
    pid = pair["id"]
    dest = f"{data_root}/{pair['tier']}/{pid}"
    return {
        "id": pid,
        "species": f"{pair['ref_species']} -> {pair['tgt_species']}",
        "cross_species": pair["cross_species"],
        "dimension": pair["dimension"],
        "divergence_class": pair["tier"],
        "ref_genome": f"{dest}/ref.fna",
        "ref_gff": f"{dest}/ref.gff",
        "ref_proteins": None,
        "tgt_genome": f"{dest}/tgt.fna",
        "ref_chrom": pair["ref_chrom"],
        "tgt_chrom": pair["tgt_chrom"],
        "miniprot_target_space": "transcript",
        "annotation_database": "RefSeq",
    }


def build_entries(reg: dict) -> list:
    """The list of benchmarks.json entries for every tiered pair."""
    return [build_entry(p, reg["data_root"]) for p in reg["pairs"]]


def merge(existing: list, new_entries: list) -> list:
    """Existing entries (order preserved) minus any our ids, then our entries.

    Idempotent: re-running replaces our entries in place at the tail rather than
    duplicating them.
    """
    new_ids = {e["id"] for e in new_entries}
    kept = [e for e in existing if e.get("id") not in new_ids]
    return kept + list(new_entries)


def load_registry(path: Path = REGISTRY) -> dict:
    with open(path) as fh:
        return json.load(fh)


def main(argv=None) -> int:
    argv = list(argv) if argv is not None else []
    dry = "--dry-run" in argv

    reg = load_registry()
    with open(BENCHMARKS) as fh:
        bench = json.load(fh)

    new_entries = build_entries(reg)
    before_ids = {e.get("id") for e in bench["benchmarks"] if e.get("id")}
    added = [e["id"] for e in new_entries if e["id"] not in before_ids]
    updated = [e["id"] for e in new_entries if e["id"] in before_ids]

    bench["benchmarks"] = merge(bench["benchmarks"], new_entries)

    print(f"tiered pairs: {len(new_entries)}  "
          f"(added={len(added)}, updated/replaced={len(updated)})", flush=True)
    for e in new_entries:
        print(f"  {e['id']:32s} {e['divergence_class']:26s} "
              f"ref_chrom={e['ref_chrom']} tgt_chrom={e['tgt_chrom']}", flush=True)

    if dry:
        print("\n--dry-run: benchmarks.json NOT written", flush=True)
        return 0

    with open(BENCHMARKS, "w") as fh:
        json.dump(bench, fh, indent=2)
        fh.write("\n")
    print(f"\nwrote {BENCHMARKS}", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:] or None))
