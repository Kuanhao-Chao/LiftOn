#!/usr/bin/env python
"""start_lost headroom diagnostic — how many transcripts does the check miss?

`variants.find_variants` decides `start_lost` from four AND-ed clauses, and the
first two read `align_dna.query_aln[0:3]`: the first three columns of the FULL
transcript alignment. For any mRNA with a 5'UTR those columns are UTR sequence,
so when the UTR matches between query and reference -- the common case on close
pairs -- the chain short-circuits and a genuine start loss is never flagged.
`start_lost` is one of four mutations that gate ORF rescue, so the miss silently
suppresses re-search.

`cds_span` (added for GH #46) already locates the real start codon. This runs
the byte-neutral `LIFTON_START_LOST_DIAG` probe, which evaluates the SAME four
clauses against the three alignment columns AT the CDS start and records both
verdicts per transcript, then reports:

  n_scored          transcripts the probe saw (coding, with a locatable CDS)
  n_current         start_lost under the shipped test
  n_scoped          start_lost under the CDS-scoped test
  n_missed          scoped says yes, shipped says no
  n_newly_rescued   THE NUMBER THAT MATTERS -- missed AND not already entering
                    ORF rescue via frameshift / stop_missing / stop_codon_gain
  n_false_positive  shipped says yes, scoped says no (the UTR happened to differ)

A transcript already re-searching for another reason gains nothing from a
start_lost verdict, so `n_missed` overstates the opportunity and
`n_newly_rescued` is the honest figure.

DECISION GATE: build the fix only if n_newly_rescued is material. An upside-only
count is not sufficient on its own -- if it proceeds, the change needs the full
ladder judged on ORF validity and structural metrics, because an ORF re-search
can move a model without moving its identity score (the Iteration-9/13/15
lesson: the headroom gates cheaply, the A/B catches the ripple).

Usage (repo root, lifton_devel env):
    python -m benchmarks.compare.start_lost_headroom [IDS...]
"""
from __future__ import annotations

import argparse
import collections
import csv
import json
import os
import sys
from pathlib import Path

from .profiling import run_profiled
from .tool_runners import _clean_input_dbs, _compose_env

HERE = Path(__file__).resolve().parent
WORK = HERE / "work"
REG = json.loads((HERE / "benchmarks.json").read_text())
TOOLS = REG["tools"]

DEFAULT_IDS = [
    "celegans_to_briggsae",       # distant
    "drosophila_to_anopheles",    # very-distant
    "zebrafish_to_medaka",        # distant
    "rice_to_sorghum",            # distant
    "t4_human_to_xenopus",        # very-distant
    "t4_human_to_chicken",        # very-distant
    "human_to_mouse",             # distant
    "drosophila",                 # same-species control
]


def _ann_db(bid: str) -> str:
    for b in REG["benchmarks"]:
        if b["id"] == bid:
            return b.get("annotation_database", "RefSeq")
    return "RefSeq"


def _divergence(bid: str) -> str:
    for b in REG["benchmarks"]:
        if b["id"] == bid:
            return b.get("divergence_class") or b.get("tier") or ""
    return ""


def _run(bid, paths, root):
    liftoff = WORK / bid / "tools" / "liftoff" / "liftoff.gff3"
    miniprot = WORK / bid / "tools" / "miniprot" / "miniprot.gff3"
    root.mkdir(parents=True, exist_ok=True)
    out = root / "lifton.gff3"
    diag = root / "start_lost.tsv"
    if diag.exists():
        diag.unlink()
    _clean_input_dbs(paths["ref_gff"], liftoff, miniprot)

    env = dict(_compose_env(TOOLS))
    env["PYTHONHASHSEED"] = "0"
    env["LIFTON_START_LOST_DIAG"] = str(diag)
    ab_pythonpath = os.environ.get("LIFTON_AB_PYTHONPATH")
    if ab_pythonpath:
        env["PYTHONPATH"] = ab_pythonpath

    argv = [TOOLS["lifton_bin"], "-t", "1", "-ad", _ann_db(bid),
            "-g", paths["ref_gff"], "-L", str(liftoff), "-M", str(miniprot),
            "-o", str(out), paths["tgt_fa"], paths["ref_fa"]]
    pr = run_profiled(argv, label=f"start_lost_{bid}", log_dir=root / "logs",
                      env=env, log=print)
    if pr.exit_code != 0 or not diag.exists():
        raise RuntimeError(f"{bid}: run failed (exit {pr.exit_code}); "
                           f"see {pr.stderr_path}")
    return diag


def _summarise(diag: Path) -> dict:
    counts = collections.Counter()
    utr_lengths = []
    examples = []
    with diag.open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            counts["n_scored"] += 1
            current = row["current"] == "1"
            scoped = row["scoped"] == "1"
            counts["n_current"] += current
            counts["n_scoped"] += scoped
            if scoped and not current:
                counts["n_missed"] += 1
                if row["newly_rescued"] == "1":
                    counts["n_newly_rescued"] += 1
                    if len(examples) < 5:
                        examples.append({
                            "transcript": row["transcript_id"],
                            "utr_len": int(row["utr_len"]),
                            "first3_full": row["dna_first3"],
                            "first3_at_cds": row["scoped_first3"],
                            "ref_at_cds": row["ref_first3"],
                            "protein_first": row["protein_first"],
                        })
            elif current and not scoped:
                counts["n_false_positive"] += 1
            if row["utr_len"]:
                utr_lengths.append(int(row["utr_len"]))
    scored = counts["n_scored"] or 1
    return {
        "n_scored": counts["n_scored"],
        "n_current": counts["n_current"],
        "n_scoped": counts["n_scoped"],
        "n_missed": counts["n_missed"],
        "n_newly_rescued": counts["n_newly_rescued"],
        "n_false_positive": counts["n_false_positive"],
        "frac_newly_rescued": round(counts["n_newly_rescued"] / scored, 6),
        "transcripts_with_5utr": sum(1 for length in utr_lengths if length > 0),
        "examples": examples,
    }


def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("ids", nargs="*")
    args = parser.parse_args(argv if argv is not None else sys.argv[1:])
    ids = args.ids or DEFAULT_IDS

    results = []
    for bid in ids:
        print(f"=== {bid}: start_lost headroom ===", flush=True)
        manifest = json.loads(
            (WORK / bid / "subset" / "subset.manifest.json").read_text())
        diag = _run(bid, manifest["paths"],
                    WORK / bid / "_start_lost_headroom")
        record = {"benchmark": bid,
                  "divergence": _divergence(bid) or manifest.get("species", "")}
        record.update(_summarise(diag))
        results.append(record)
        print(f"  [{bid}] scored={record['n_scored']} "
              f"current={record['n_current']} scoped={record['n_scoped']} "
              f"missed={record['n_missed']} "
              f"NEWLY_RESCUED={record['n_newly_rescued']} "
              f"({record['frac_newly_rescued']:.4%}) "
              f"false_pos={record['n_false_positive']}", flush=True)
        (HERE / "start_lost_headroom.json").write_text(
            json.dumps(results, indent=2))
        _write_md(results, HERE / "start_lost_headroom.md")

    total_new = sum(r["n_newly_rescued"] for r in results)
    total_scored = sum(r["n_scored"] for r in results)
    print(f"\nTOTAL newly-rescued: {total_new} of {total_scored} scored "
          f"({total_new / (total_scored or 1):.4%})", flush=True)
    print("DONE: wrote start_lost_headroom.json + .md", flush=True)
    return 0


def _write_md(results, path: Path):
    lines = [
        "## `start_lost` headroom (byte-neutral diagnostic)\n",
        "`variants.find_variants` decides `start_lost` from four AND-ed clauses "
        "whose first two read `align_dna.query_aln[0:3]` — the first three "
        "columns of the FULL transcript alignment, which for an mRNA with a "
        "5'UTR is UTR sequence. When that UTR matches, the chain "
        "short-circuits and a genuine start loss is never flagged, and "
        "`start_lost` gates ORF rescue.\n",
        "`LIFTON_START_LOST_DIAG` re-evaluates the same four clauses against "
        "the three columns AT the CDS start (`cds_span`, from GH #46) and "
        "records both verdicts. **newly rescued** — a start loss the shipped "
        "test misses AND that no other ORF trigger would have caught anyway — "
        "is the honest figure; **missed** overstates it.\n",
        "| Dataset | divergence | scored | current | scoped | missed | "
        "**newly rescued** | % | false pos |",
        "|---|---|---|---|---|---|---|---|---|",
    ]
    for r in results:
        lines.append(
            f"| {r['benchmark']} | {r.get('divergence','')} | {r['n_scored']} | "
            f"{r['n_current']} | {r['n_scoped']} | {r['n_missed']} | "
            f"**{r['n_newly_rescued']}** | {r['frac_newly_rescued']:.4%} | "
            f"{r['n_false_positive']} |"
        )
    total_new = sum(r["n_newly_rescued"] for r in results)
    total_scored = sum(r["n_scored"] for r in results)
    lines.append(
        f"\n**TOTAL newly-rescued: {total_new} of {total_scored} scored "
        f"({total_new / (total_scored or 1):.4%}).**\n"
    )
    path.write_text("\n".join(lines))


if __name__ == "__main__":
    raise SystemExit(main())
