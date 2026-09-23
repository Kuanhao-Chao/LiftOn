#!/usr/bin/env python
"""Second-locus rescue A/B: default vs ``--rescue-second-locus``.

A whole-genome duplication gives the target two genes where the reference has
one. The rescue's reference-gene dedup refuses the second; measured against
zebrafish's own GRCz11 annotation, that refusal hides 2,051 real target genes
on human to zebrafish (notes/coortholog_recall_measurement_2026-09.md), and the
single-cell gate passes (notes/second_locus_rescue_gate.md).

Human to zebrafish is also the most favourable cell imaginable -- the teleost
duplication is the reason the idea exists -- so this asks the ladder whether
anything generalises.

**What this harness can and cannot prove.** The neutral evaluator keys on
``ref_mrna_id`` and keeps the best model per reference transcript. A second
target gene for a reference gene that is ALREADY scored therefore does not
appear as "added"; it only competes for the same maximum. Reference-keyed
recall is precisely the instrument that cannot see this feature, which is why
the finding needed the target's own annotation in the first place.

So the gate here is a SAFETY gate, and it is stated as one:

  n_lost == 0                 no reference transcript stops being scored
  common_regressed == 0       no existing model gets worse
  validity not worse
  n_placed > 0                the flag did something
  mean_pi_placed >= floor     what it placed is not garbage
  n_overlapping == 0          nothing placed on top of an emitted model,
                              measured against the pipeline's own -overlap gate

Whether the placed genes are REAL needs the target's own annotation, which
exists on disk for human -> zebrafish only. For every other cell this reports
how many were placed and how good they look, and says nothing about whether
they should exist. That is the honest limit of the available data.

Usage (repo root, lifton_devel env):
    python -m benchmarks.compare.second_locus_ab [IDS...]
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import sys
from collections import defaultdict
from pathlib import Path

from . import evaluator
from .profiling import run_profiled
from .tool_runners import _clean_input_dbs, _compose_env

HERE = Path(__file__).resolve().parent
#: The subset trees are large and live with the main checkout, not with a
#: detached worktree running a candidate build, so they are addressable
#: separately. Both arms still share one `work/`, which is the point.
WORK = Path(os.environ.get("LIFTON_AB_WORK") or (HERE / "work"))
REG = json.loads((HERE / "benchmarks.json").read_text())
TOOLS = REG["tools"]

DEFAULT_IDS = [
    "drosophila",                 # same-species control, expected inert
    "human_to_mouse",             # distant
    "celegans_to_briggsae",       # distant
    "rice_to_sorghum",            # distant, a plant with its own duplications
    "drosophila_to_anopheles",    # very distant
    "zebrafish_to_medaka",        # very distant, teleost on both sides
    "t4_human_to_chicken",        # very distant
    "t4_human_to_xenopus",        # very distant
]
#: The rescue's own identity floor; the adaptive floor may lower it per cell.
FLOOR = 0.5


def _ann_db(bid):
    for b in REG["benchmarks"]:
        if b["id"] == bid:
            return b.get("annotation_database", "RefSeq")
    return "RefSeq"


def _run_state(bid, state, flags, p, root, pythonpath):
    liftoff = WORK / bid / "tools" / "liftoff" / "liftoff.gff3"
    miniprot = WORK / bid / "tools" / "miniprot" / "miniprot.gff3"
    statedir = root / state
    statedir.mkdir(parents=True, exist_ok=True)
    out = statedir / f"{state}.gff3"
    _clean_input_dbs(p["ref_gff"], liftoff, miniprot)
    # DETERMINISTIC: -t1, no -copies, cached -L/-M, so the two arms differ only
    # by the flag. (A -t8 -copies A/B on repeat-rich genomes injects run-to-run
    # noise that reads as n_lost -- the CLAUDE.md caveat.)
    argv = [TOOLS["lifton_bin"], "-t", "1", "-ad", _ann_db(bid),
            "-g", p["ref_gff"], "-L", str(liftoff), "-M", str(miniprot),
            "-o", str(out), *flags, p["tgt_fa"], p["ref_fa"]]
    env = dict(_compose_env(TOOLS))
    if pythonpath:
        # Prepend: the tool environment may need its own entries, and a
        # candidate build has to win without removing them.
        existing = env.get("PYTHONPATH", "")
        env["PYTHONPATH"] = (f"{pythonpath}{os.pathsep}{existing}"
                             if existing else pythonpath)
    for name in ("LIFTON_RESCUE_SECOND_LOCUS", "LIFTON_RESCUE_SECOND_LOCUS_MAX"):
        env.pop(name, None)
    pr = run_profiled(argv, label=f"second_locus_ab_{bid}_{state}",
                      log_dir=root / "logs", env=env, log=print)
    if pr.exit_code != 0 or not out.exists() or out.stat().st_size == 0:
        raise RuntimeError(f"{bid}: state {state} failed (exit {pr.exit_code}); "
                           f"see {pr.stderr_path}")
    return out


def _pi_by_ref(tsv):
    out = {}
    if not tsv.exists():
        return out
    with tsv.open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if str(row.get("is_coding")).strip().lower() not in ("1", "true", "yes"):
                continue
            value = row.get("protein_identity")
            if value in (None, "", "None"):
                continue
            try:
                value = float(value)
            except ValueError:
                continue
            key = row["ref_mrna_id"]
            out[key] = value if key not in out else max(out[key], value)
    return out


def _genes(path):
    out = {}
    with open(path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "gene":
                continue
            a = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)
            out[a.get("ID", "")] = (f[0], int(f[3]), int(f[4]))
    return out


def _placed(path):
    """Protein identities of the transcripts tagged as second-locus models."""
    values = []
    with open(path) as handle:
        for line in handle:
            if line.startswith("#") or "lifton_rescue_second_locus=true" not in line:
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] not in ("mRNA", "transcript"):
                continue
            a = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)
            if "protein_identity" in a:
                values.append(float(a["protein_identity"]))
    return values


def _overlapping(off_path, on_path, overlap_gate=0.10):
    """Added genes covering more than the pipeline's gate of an emitted gene."""
    off, on = _genes(off_path), _genes(on_path)
    index = defaultdict(list)
    for key, (seqid, start, end) in off.items():
        index[seqid].append((start, end, key))
    for seqid in index:
        index[seqid].sort()
    worst = []
    for key in set(on) - set(off):
        seqid, start, end = on[key]
        length = end - start + 1
        for other_start, other_end, other in index.get(seqid, ()):
            if other_start > end:
                break
            if other_end < start:
                continue
            shared = min(end, other_end) - max(start, other_start) + 1
            if shared / length >= overlap_gate:
                worst.append((round(shared / length, 3), key, other))
    return sorted(worst, reverse=True)


def _validate(gff, pythonpath=""):
    """Validator error count, from the structured summary the other A/B
    harnesses share -- a line grep for "error" counts the clean summary line
    `Errors : 0` itself -- and run on the build under test: without the
    PYTHONPATH pin the validator came from whatever `lifton` the tool
    environment has installed."""
    env = dict(_compose_env(TOOLS))
    if pythonpath:
        existing = env.get("PYTHONPATH", "")
        env["PYTHONPATH"] = (f"{pythonpath}{os.pathsep}{existing}"
                             if existing else pythonpath)
    return evaluator.count_gff3_validator_errors(gff, TOOLS["lifton_python"], env)


def _mean(values):
    values = [v for v in values if v is not None]
    return round(sum(values) / len(values), 5) if values else None


def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("--pythonpath", default=os.environ.get("LIFTON_AB_PYTHONPATH", ""),
                        help="pin both arms to one build (a detached worktree)")
    parser.add_argument("ids", nargs="*")
    args = parser.parse_args(argv if argv is not None else sys.argv[1:])
    ids = args.ids or DEFAULT_IDS

    results = []
    for bid in ids:
        print(f"=== {bid}: second-locus rescue A/B ===", flush=True)
        work = WORK / bid
        manifest = json.loads((work / "subset" / "subset.manifest.json").read_text())
        paths = manifest["paths"]
        root = work / "_second_locus_ab"
        root.mkdir(parents=True, exist_ok=True)

        # The feature is default-on since the promotion, so the baseline arm
        # has to say so explicitly. Passing [] here would compare on with on.
        off = _run_state(bid, "default", ["--no-rescue-second-locus"],
                         paths, root, args.pythonpath)
        reference, ref_index = evaluator.build_reference(
            paths["ref_gff"], paths["ref_fa"], log=print)
        eval_dir = root / "eval"
        eval_dir.mkdir(parents=True, exist_ok=True)
        evaluator.evaluate_tool("default", str(off), paths["tgt_fa"], reference,
                                manifest, eval_dir, None, log=print,
                                ref_index=ref_index, threads=8)
        pi_off = _pi_by_ref(eval_dir / "default.transcripts.tsv")
        validity_off = _validate(off, args.pythonpath)

        on = _run_state(bid, "second_locus", ["--rescue-second-locus"],
                        paths, root, args.pythonpath)
        evaluator.evaluate_tool("second_locus", str(on), paths["tgt_fa"],
                                reference, manifest, eval_dir, None, log=print,
                                ref_index=ref_index, threads=8)
        pi_on = _pi_by_ref(eval_dir / "second_locus.transcripts.tsv")
        validity_on = _validate(on, args.pythonpath)

        lost = set(pi_off) - set(pi_on)
        common = set(pi_off) & set(pi_on)
        regressed = sum(1 for k in common if pi_on[k] < pi_off[k] - 1e-9)
        improved = sum(1 for k in common if pi_on[k] > pi_off[k] + 1e-9)
        placed = _placed(on)
        overlapping = _overlapping(off, on)

        gate = bool(not lost and regressed == 0 and validity_on <= validity_off
                    and not overlapping
                    and (not placed or _mean(placed) >= FLOOR))
        record = {
            "benchmark": bid,
            "divergence": manifest.get("divergence_class") or manifest.get("species", ""),
            "n_placed": len(placed),
            "mean_pi_placed": _mean(placed),
            "n_lost": len(lost),
            "common_improved": improved,
            "common_regressed": regressed,
            "n_overlapping_above_gate": len(overlapping),
            "overlapping_examples": overlapping[:5],
            "validity": {"default": validity_off, "second_locus": validity_on},
            "safety_gate_pass": gate,
        }
        results.append(record)
        print(f"  [{bid}] placed={len(placed)} mean_pi={record['mean_pi_placed']} "
              f"lost={len(lost)} regressed={regressed} "
              f"overlapping={len(overlapping)} "
              f"val {validity_off}->{validity_on} "
              f"safety {'PASS' if gate else 'FAIL'}", flush=True)

    (HERE / "second_locus_ab.json").write_text(json.dumps(results, indent=2) + "\n")
    passed = sum(1 for r in results if r["safety_gate_pass"])
    print(f"\nSAFETY GATE: {passed}/{len(results)} cells pass.", flush=True)
    print("Reference-keyed recall cannot see what this feature adds; whether "
          "the placed genes are real needs the target's own annotation.",
          flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
