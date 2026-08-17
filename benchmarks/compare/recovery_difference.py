#!/usr/bin/env python3
"""Characterise which coding transcripts only one release recovers.

A paired comparison can report a small completeness difference without saying
anything about the models behind it. One release recovering fewer coding
transcripts is only interpretable once you know what those transcripts look
like: dropping weak, low-identity models is a different result from dropping
well-supported ones, and the two justify opposite conclusions.

This reducer reads the neutral evaluator's per-transcript tables for both arms
of a paired cell and partitions the coding transcripts into the set only the
reference recovered, the set only the candidate recovered, and reports the
quality distribution of each. It reads evaluation output only; it never runs
LiftOn and never re-scores.
"""
from __future__ import annotations

import argparse
import csv
import json
import statistics
import sys
from pathlib import Path
from typing import Any, Iterable, Mapping

SCHEMA_VERSION = 1
METHOD = "paired-recovery-difference-v1"
DEFAULT_WEAK_THRESHOLD = 0.5
ARMS = ("candidate", "reference")


def _is_true(value: str | None) -> bool:
    return (value or "").strip().lower() in {"true", "1", "yes"}


def _read_coding(path: Path) -> dict[str, dict[str, str]]:
    """Index one arm's coding transcripts by reference transcript id."""

    rows: dict[str, dict[str, str]] = {}
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if _is_true(row.get("is_coding")):
                rows[row["ref_mrna_id"]] = row
    if not rows:
        raise ValueError(f"no coding transcripts in {path}")
    return rows


def _float(value: str | None) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return 0.0


def _describe(
    identifiers: Iterable[str],
    rows: Mapping[str, Mapping[str, str]],
    *,
    weak_threshold: float,
) -> dict[str, Any]:
    identifiers = sorted(identifiers)
    if not identifiers:
        return {
            "n": 0,
            "mean_protein_identity": None,
            "median_protein_identity": None,
            "n_below_threshold": 0,
            "fraction_below_threshold": None,
            "n_orf_valid": 0,
            "fraction_orf_valid": None,
            "transcripts": [],
        }
    identities = [_float(rows[i].get("protein_identity")) for i in identifiers]
    below = sum(1 for value in identities if value < weak_threshold)
    orf_valid = sum(1 for i in identifiers if _is_true(rows[i].get("orf_valid")))
    return {
        "n": len(identifiers),
        "mean_protein_identity": round(statistics.fmean(identities), 6),
        "median_protein_identity": round(statistics.median(identities), 6),
        "n_below_threshold": below,
        "fraction_below_threshold": round(below / len(identifiers), 6),
        "n_orf_valid": orf_valid,
        "fraction_orf_valid": round(orf_valid / len(identifiers), 6),
        # Bounded so the document stays readable; the counts above are complete.
        "transcripts": identifiers[:50],
    }


def compare_arms(
    candidate_tsv: Path,
    reference_tsv: Path,
    *,
    weak_threshold: float = DEFAULT_WEAK_THRESHOLD,
) -> dict[str, Any]:
    """Partition coding transcripts by which arm recovered them."""

    candidate = _read_coding(Path(candidate_tsv))
    reference = _read_coding(Path(reference_tsv))
    shared = set(candidate) & set(reference)
    reference_only = {
        key for key in shared
        if _is_true(reference[key].get("recovered"))
        and not _is_true(candidate[key].get("recovered"))
    }
    candidate_only = {
        key for key in shared
        if _is_true(candidate[key].get("recovered"))
        and not _is_true(reference[key].get("recovered"))
    }
    return {
        "denominator_coding": len(shared),
        "weak_threshold": weak_threshold,
        "reference_only": _describe(
            reference_only, reference, weak_threshold=weak_threshold,
        ),
        "candidate_only": _describe(
            candidate_only, candidate, weak_threshold=weak_threshold,
        ),
        "net_candidate_minus_reference": len(candidate_only) - len(reference_only),
    }


def _repetition_one_directories(campaign_root: Path) -> dict[str, Path]:
    """Locate each transfer's first repetition, which carries full scoring."""

    found: dict[str, Path] = {}
    for path in sorted(Path(campaign_root).rglob("*__repetition_01/pair/evaluation")):
        cell = path.parent.parent.name
        benchmark = cell.split("__")[1]
        found.setdefault(benchmark, path)
    return found


def build_document(
    campaign_root: Path, *, weak_threshold: float = DEFAULT_WEAK_THRESHOLD,
) -> dict[str, Any]:
    transfers = {}
    for benchmark, directory in sorted(_repetition_one_directories(campaign_root).items()):
        candidate_tsv = directory / "candidate.transcripts.tsv"
        reference_tsv = directory / "reference.transcripts.tsv"
        if not (candidate_tsv.is_file() and reference_tsv.is_file()):
            continue
        transfers[benchmark] = compare_arms(
            candidate_tsv, reference_tsv, weak_threshold=weak_threshold,
        )
    if not transfers:
        raise ValueError(f"no scored repetitions under {campaign_root}")
    return {
        "schema_version": SCHEMA_VERSION,
        "method": METHOD,
        "campaign_root": str(Path(campaign_root).resolve()),
        "weak_threshold": weak_threshold,
        "transfers": transfers,
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--campaign-root", type=Path, required=True)
    parser.add_argument("--output", type=Path)
    parser.add_argument(
        "--weak-threshold", type=float, default=DEFAULT_WEAK_THRESHOLD,
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    document = build_document(
        args.campaign_root, weak_threshold=args.weak_threshold,
    )
    text = json.dumps(document, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text)
    else:
        sys.stdout.write(text)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
