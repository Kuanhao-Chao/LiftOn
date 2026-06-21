#!/usr/bin/env python
"""Iteration 22 A/B — regime-gated miniprot-only rescue: default vs --miniprot-rescue.

Tests the NEW (Iteration-22) rescue on the regime the audit identified as the
real deficit (distant + very-distant), where the parent's Iter-15 NO-GO does NOT
apply. The Iter-15 trap was rescuing at *weak-Liftoff* loci (overlapping models
-> duplicates); this rescue fires ONLY for genes the DNA lift missed ENTIRELY
(the candidate already cleared Step 8's overlap gate, so there is no lifted gene
at the locus -> 0-redundant by construction), gated by a protein-identity floor
instead of the tight length band. Mirrors lifton2's shipped Iter-11 (+0.047 mean
completeness, 0 redundant on all 10).

States (all reuse cached -L/-M; rescue toggled by the flag + env floor):
  "default"          : no flag (rescue OFF)
  "rescue_<floor>"   : --miniprot-rescue, LIFTON_MINIPROT_RESCUE_MIN_ID=<floor>

Scored by the independent re-alignment evaluator (keyed by ref_mrna_id):
  - n_added        : coding ref transcripts rescue emits but default does not
                     (genuine completeness — the DNA lift missed them).
  - mean_pi_added  : their protein identity (legit high, garbage low).
  - n_lost         : ref transcripts default has but rescue does not (off ⊆ on
                     gate — MUST be 0).
  - common regressed/improved : per-transcript change on the shared set (rescue
                     must not perturb existing models — expect 0 regressed).
  - n_rescued_tagged : features tagged lifton_rescue=miniprot_only in the output.
  - n_redundant    : n_rescued_tagged − n_added (rescued models that did NOT add
                     a new ref id = duplicates; the Iter-15 failure mode). MUST
                     be ~0 — this is the load-bearing gate.
  - validity (gff3_validator error-line delta).

Promotion gate (per floor): n_added > 0 AND n_lost == 0 AND common regressed == 0
AND n_redundant == 0 AND mean_pi_added >= floor AND validity not worse, on the
distant/very-distant tier. A clean pass -> promote (opt-in first, then default-on
after the full-genome confirm). Duplicate-heavy or regressing -> NO-GO (revert).

Usage (repo root, lifton_devel env):
    python -m benchmarks.compare.miniprot_rescue_ab [--floors 0.3,0.5,0.7] [IDS...]
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import subprocess
import sys
from pathlib import Path

from . import evaluator
from .profiling import run_profiled
from .tool_runners import _clean_input_dbs, _compose_env

HERE = Path(__file__).resolve().parent
WORK = HERE / "work"
REG = json.loads((HERE / "benchmarks.json").read_text())
TOOLS = REG["tools"]

# Distant + very-distant datasets with strong headroom (genuine_new@0.5) and
# cached -L/-M, plus a same-species control (drosophila, expected ~flat).
DEFAULT_IDS = [
    "drosophila",                 # same-species control
    "human_to_mouse",             # distant
    "celegans_to_briggsae",       # distant (strong)
    "rice_to_sorghum",            # distant
    "drosophila_to_anopheles",    # very-distant (strong)
    "zebrafish_to_medaka",        # very-distant (strong)
    "t4_human_to_chicken",        # very-distant (strong)
    "t4_human_to_xenopus",        # very-distant (strong)
]
DEFAULT_FLOORS = [0.5]


def _ann_db(bid: str) -> str:
    for b in REG["benchmarks"]:
        if b["id"] == bid:
            return b.get("annotation_database", "RefSeq")
    return "RefSeq"


def _run_state(bid, state, floor, p, root):
    liftoff = WORK / bid / "tools" / "liftoff" / "liftoff.gff3"
    miniprot = WORK / bid / "tools" / "miniprot" / "miniprot.gff3"
    statedir = root / state
    statedir.mkdir(parents=True, exist_ok=True)
    out = statedir / f"{state}.gff3"
    _clean_input_dbs(p["ref_gff"], liftoff, miniprot)
    flags = [] if floor is None else ["--miniprot-rescue"]
    # DETERMINISTIC: -t1, NO -copies, NO --native. The rescue is a purely
    # additive elif in Step 8, so deterministically default ⊆ rescue and n_lost
    # MUST be 0. Running two independent -t8 -copies invocations on repeat-rich
    # distant genomes injects -copies run-to-run noise that masquerades as
    # n_lost/n_redundant (the CLAUDE.md caveat: byte-identity authority is the
    # fixed-input run, not a fresh -copies A/B). -t1 no-copies + cached -L/-M is
    # deterministic (proven by the 24-cell matrix + the drosophila cmp), so
    # default vs rescue then differ ONLY by the flag.
    argv = [TOOLS["lifton_bin"], "-t", "1", "-ad", _ann_db(bid),
            "-g", p["ref_gff"], "-L", str(liftoff), "-M", str(miniprot),
            "-o", str(out), *flags,
            p["tgt_fa"], p["ref_fa"]]
    env = dict(_compose_env(TOOLS))
    if floor is not None:
        env["LIFTON_MINIPROT_RESCUE"] = "1"
        env["LIFTON_MINIPROT_RESCUE_MIN_ID"] = str(floor)
    else:
        env.pop("LIFTON_MINIPROT_RESCUE", None)
    pr = run_profiled(argv, label=f"mp_rescue_ab_{state}",
                      log_dir=root / "logs", env=env, log=print)
    if pr.exit_code != 0 or not out.exists() or out.stat().st_size == 0:
        raise RuntimeError(f"{bid}: state {state} failed (exit {pr.exit_code}); "
                           f"see {pr.stderr_path}")
    return out, pr


def _pi_by_ref(tsv: Path) -> dict:
    out = {}
    if not tsv.exists():
        return out
    with tsv.open() as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if str(row.get("is_coding")).strip().lower() not in ("1", "true", "yes"):
                continue
            pi = row.get("protein_identity")
            if pi in (None, "", "None"):
                continue
            try:
                v = float(pi)
            except ValueError:
                continue
            rid = row["ref_mrna_id"]
            out[rid] = v if rid not in out else max(out[rid], v)
    return out


def _n_tagged_rescues(gff: Path) -> int:
    n = 0
    with gff.open() as fh:
        for ln in fh:
            if ln.startswith("#"):
                continue
            c = ln.split("\t")
            if len(c) > 8 and c[2] == "mRNA" and "lifton_rescue=miniprot_only" in c[8]:
                n += 1
    return n


def _validate(gff: Path):
    r = subprocess.run([TOOLS["lifton_python"], "-m", "lifton.gff3_validator", str(gff)],
                       env=_compose_env(TOOLS), capture_output=True, text=True)
    out = (r.stdout or "") + (r.stderr or "")
    return r.returncode, sum(1 for ln in out.splitlines() if "error" in ln.lower())


def _mean(xs):
    xs = [x for x in xs if x is not None]
    return round(sum(xs) / len(xs), 5) if xs else None


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--floors", default=",".join(str(f) for f in DEFAULT_FLOORS))
    ap.add_argument("ids", nargs="*")
    args = ap.parse_args(argv if argv is not None else sys.argv[1:])
    ids = args.ids or DEFAULT_IDS
    floors = [float(x) for x in args.floors.split(",") if x]

    results = []
    for bid in ids:
        print(f"=== {bid}: regime-gated miniprot-only rescue A/B (Iter 22) ===", flush=True)
        W = WORK / bid
        man = json.loads((W / "subset" / "subset.manifest.json").read_text())
        p = man["paths"]
        root = W / "_mp_rescue22_ab"
        root.mkdir(parents=True, exist_ok=True)

        out_def, pr_def = _run_state(bid, "default", None, p, root)
        ref, ref_index = evaluator.build_reference(p["ref_gff"], p["ref_fa"], log=print)
        eval_dir = root / "eval"
        eval_dir.mkdir(parents=True, exist_ok=True)
        evaluator.evaluate_tool("default", str(out_def), p["tgt_fa"], ref, man,
                                eval_dir, None, log=print, ref_index=ref_index, threads=8)
        pi_def = _pi_by_ref(eval_dir / "default.transcripts.tsv")
        val_def = _validate(out_def)

        for floor in floors:
            state = f"rescue_{floor}"
            print(f"--- {bid} {state} ---", flush=True)
            out_on, pr_on = _run_state(bid, state, floor, p, root)
            evaluator.evaluate_tool(state, str(out_on), p["tgt_fa"], ref, man,
                                    eval_dir, None, log=print, ref_index=ref_index, threads=8)
            pi_on = _pi_by_ref(eval_dir / f"{state}.transcripts.tsv")

            added = set(pi_on) - set(pi_def)
            lost = set(pi_def) - set(pi_on)
            common = set(pi_on) & set(pi_def)
            regressed = sum(1 for k in common if pi_on[k] < pi_def[k] - 1e-9)
            improved = sum(1 for k in common if pi_on[k] > pi_def[k] + 1e-9)
            added_pis = [pi_on[k] for k in added]
            mean_pi_added = _mean(added_pis)
            n_tagged = _n_tagged_rescues(out_on)
            n_redundant = max(0, n_tagged - len(added))
            val_on = _validate(out_on)

            gate = bool(len(added) > 0 and len(lost) == 0 and regressed == 0
                        and n_redundant == 0 and mean_pi_added is not None
                        and mean_pi_added >= floor and val_on[1] <= val_def[1])
            rec = {
                "benchmark": bid, "floor": floor,
                "divergence": man.get("divergence_class") or man.get("species", ""),
                "n_scored_default": len(pi_def), "n_scored_rescue": len(pi_on),
                "mean_pi_default": _mean(pi_def.values()),
                "mean_pi_rescue": _mean(pi_on.values()),
                "n_added": len(added), "mean_pi_added": mean_pi_added,
                "frac_added_ge_floor": (round(sum(1 for v in added_pis if v >= floor) / len(added_pis), 3)
                                        if added_pis else None),
                "n_lost": len(lost),
                "common_improved": improved, "common_regressed": regressed,
                "n_rescued_tagged": n_tagged, "n_redundant": n_redundant,
                "validity": {"default": val_def[1], "rescue": val_on[1]},
                "completeness_gain": round((len(pi_on) - len(pi_def)) / (man.get("n_subset_proteins") or len(pi_def) or 1), 5),
                "gate_pass": gate,
            }
            results.append(rec)
            print(f"  [{bid}@{floor}] n_added={len(added)} mean_pi_added={mean_pi_added} "
                  f"n_lost={len(lost)} regressed={regressed} tagged={n_tagged} "
                  f"redundant={n_redundant} val {val_def[1]}->{val_on[1]} "
                  f"gate {'PASS' if gate else 'FAIL'}", flush=True)

    (HERE / "miniprot_rescue_ab.json").write_text(json.dumps(results, indent=2))
    _write_md(results, HERE / "miniprot_rescue_ab.md")
    n_pass = sum(1 for r in results if r["gate_pass"])
    print(f"\nGATE: {n_pass}/{len(results)} (dataset,floor) cells pass.", flush=True)
    print("DONE: wrote miniprot_rescue_ab.json + miniprot_rescue_ab.md", flush=True)
    return 0


def _write_md(results, path: Path):
    lines = [
        "## Regime-gated miniprot-only rescue A/B — Iteration 22\n",
        "`rescue_<floor>` = `--miniprot-rescue` with PI floor `<floor>`; cached "
        "`-L`/`-M`, scored by the independent re-alignment evaluator (keyed by "
        "ref_mrna_id). **n_added** = coding ref transcripts rescue emits that the "
        "DNA lift missed (genuine completeness); **n_redundant** = rescued models "
        "that did NOT add a new ref id (the Iter-15 duplicate failure mode — MUST "
        "be 0); **n_lost** = off ⊆ on gate (MUST be 0). Gate: n_added>0 AND "
        "n_lost=0 AND 0 regressed AND n_redundant=0 AND mean PI added ≥ floor AND "
        "validity not worse.\n",
        "| Dataset | floor | n_added | mean PI added | Δcompl | n_lost | regr | tagged | redundant | val d→r | gate |",
        "|---|---|---|---|---|---|---|---|---|---|---|",
    ]
    for r in results:
        v = r["validity"]
        lines.append("| {} | {} | {} | {} | {} | {} | {} | {} | {} | {}→{} | {} |".format(
            r["benchmark"], r["floor"], r["n_added"], r["mean_pi_added"],
            r["completeness_gain"], r["n_lost"], r["common_regressed"],
            r["n_rescued_tagged"], r["n_redundant"], v["default"], v["rescue"],
            "**PASS**" if r["gate_pass"] else "NO-GO"))
    path.write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    sys.exit(main())
