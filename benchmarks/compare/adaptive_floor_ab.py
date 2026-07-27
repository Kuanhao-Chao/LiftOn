#!/usr/bin/env python
"""Divergence-adaptive rescue-floor A/B (Phase 2) — adaptive vs fixed-0.5.

The miniprot-only rescue (Iteration 23, default-ON) uses a FIXED protein-identity
floor of 0.5. This A/B tests the divergence-adaptive floor (mirrors lifton2's
shipped design): lower the floor toward 0.30 as DNA-lift gene recall drops, so
more genuinely-missing genes are recovered on the hardest pairs while
same/close-species behaviour is unchanged.

Both arms run the rescue (--miniprot-rescue); only the floor policy differs:
  "fixed"    : LIFTON_RESCUE_ADAPTIVE_FLOOR=0  (fixed floor 0.5 — precision baseline)
  "adaptive" : LIFTON_RESCUE_ADAPTIVE_FLOOR=1  (floor in [0.30, 0.50] by recall)

The adaptive floor only ever LOWERS the floor within an already-additive
separate pass, so adaptive is a strict SUPERSET of fixed: n_lost (fixed − adaptive)
MUST be 0. Scored by the independent re-alignment evaluator (keyed by ref_mrna_id):
  - n_new        : coding ref transcripts adaptive emits that fixed does not
                   (the extra recall the lower floor buys).
  - mean_pi_new  : their protein identity (>= the adaptive floor by construction).
  - n_lost       : fixed − adaptive (MUST be 0 — adaptive only lowers the floor).
  - n_redundant  : tagged rescues that added no new ref id (MUST be 0).
  - inert cells  : high-recall pairs (recall >= R_HIGH) keep floor 0.5 -> n_new=0.
  - validity     : gff3_validator error-line delta (must not worsen).

Per-cell gate: n_lost==0 AND n_redundant==0 AND (n_new>0 on divergent OR inert on
high-recall) AND common regressed==0 AND validity not worse. Report the honest
precision/recall trade (mean PI over the recovered set dips as lower-identity but
genuine genes are added).

Usage (repo root, lifton_devel env):
    python -m benchmarks.compare.adaptive_floor_ab [IDS...]
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

from . import evaluator
from .profiling import run_profiled
from .tool_runners import _clean_input_dbs, _compose_env

HERE = Path(__file__).resolve().parent
WORK = HERE / "work"
REG = json.loads((HERE / "benchmarks.json").read_text())
TOOLS = REG["tools"]

# Same-species/close control (inert, recall high) + distant + very-distant (where
# the adaptive floor should lower toward 0.30 and recover genuinely-missing genes).
DEFAULT_IDS = [
    "celegans_to_briggsae",       # distant
    "drosophila_to_anopheles",    # very-distant
    "zebrafish_to_medaka",        # very-distant
    "rice_to_sorghum",            # distant (higher recall -> may be ~inert)
    "t4_human_to_xenopus",        # very-distant
    "t4_human_to_chicken",        # very-distant
    "human_to_mouse",             # distant (high recall -> inert control)
    "drosophila",                 # same-species control (inert)
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


def _run_arm(bid, arm, p, root):
    liftoff = WORK / bid / "tools" / "liftoff" / "liftoff.gff3"
    miniprot = WORK / bid / "tools" / "miniprot" / "miniprot.gff3"
    statedir = root / arm
    statedir.mkdir(parents=True, exist_ok=True)
    out = statedir / f"{arm}.gff3"
    _clean_input_dbs(p["ref_gff"], liftoff, miniprot)
    # DETERMINISTIC: -t1, NO -copies. Both arms run --miniprot-rescue; only the
    # adaptive-floor env differs, so fixed vs adaptive differ ONLY by the floor.
    argv = [TOOLS["lifton_bin"], "-t", "1", "-ad", _ann_db(bid),
            "-g", p["ref_gff"], "-L", str(liftoff), "-M", str(miniprot),
            "-o", str(out), "--miniprot-rescue", p["tgt_fa"], p["ref_fa"]]
    env = dict(_compose_env(TOOLS))
    env["LIFTON_MINIPROT_RESCUE"] = "1"
    env["LIFTON_RESCUE_ADAPTIVE_FLOOR"] = "0" if arm == "fixed" else "1"
    pr = run_profiled(argv, label=f"adaptive_floor_{bid}_{arm}",
                      log_dir=root / "logs", env=env, log=print)
    if pr.exit_code != 0 or not out.exists() or out.stat().st_size == 0:
        raise RuntimeError(f"{bid}: arm {arm} failed (exit {pr.exit_code}); "
                           f"see {pr.stderr_path}")
    return out


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
    # Structured validator error count (the single "Errors : N" summary line),
    # shared with the candidate-3 ladder A/B — NOT a crude grep of "error"
    # substrings (which undercounted strand regressions 1->2 vs the true 0->4).
    return evaluator.count_gff3_validator_errors(
        gff, TOOLS["lifton_python"], _compose_env(TOOLS))


def _mean(xs):
    xs = [x for x in xs if x is not None]
    return round(sum(xs) / len(xs), 5) if xs else None


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("ids", nargs="*")
    args = ap.parse_args(argv if argv is not None else sys.argv[1:])
    ids = args.ids or DEFAULT_IDS

    results = []
    for bid in ids:
        print(f"=== {bid}: adaptive-floor A/B ===", flush=True)
        W = WORK / bid
        man = json.loads((W / "subset" / "subset.manifest.json").read_text())
        p = man["paths"]
        root = W / "_adaptive_floor_ab"
        root.mkdir(parents=True, exist_ok=True)

        out_fixed = _run_arm(bid, "fixed", p, root)
        out_adapt = _run_arm(bid, "adaptive", p, root)

        ref, ref_index = evaluator.build_reference(p["ref_gff"], p["ref_fa"], log=print)
        eval_dir = root / "eval"
        eval_dir.mkdir(parents=True, exist_ok=True)
        evaluator.evaluate_tool("fixed", str(out_fixed), p["tgt_fa"], ref, man,
                                eval_dir, None, log=print, ref_index=ref_index, threads=8)
        evaluator.evaluate_tool("adaptive", str(out_adapt), p["tgt_fa"], ref, man,
                                eval_dir, None, log=print, ref_index=ref_index, threads=8)
        pi_fixed = _pi_by_ref(eval_dir / "fixed.transcripts.tsv")
        pi_adapt = _pi_by_ref(eval_dir / "adaptive.transcripts.tsv")

        new = set(pi_adapt) - set(pi_fixed)
        lost = set(pi_fixed) - set(pi_adapt)
        common = set(pi_fixed) & set(pi_adapt)
        regressed = sum(1 for k in common if pi_adapt[k] < pi_fixed[k] - 1e-9)
        new_pis = [pi_adapt[k] for k in new]
        mean_pi_new = _mean(new_pis)
        n_tagged_fixed = _n_tagged_rescues(out_fixed)
        n_tagged_adapt = _n_tagged_rescues(out_adapt)
        # extra tagged rescues adaptive emitted beyond fixed, minus the genuinely
        # new ref ids they added = redundant (should be 0).
        n_redundant = max(0, (n_tagged_adapt - n_tagged_fixed) - len(new))
        val_fixed = _validate(out_fixed)
        val_adapt = _validate(out_adapt)

        gate = bool(len(lost) == 0 and n_redundant == 0 and regressed == 0
                    and val_adapt <= val_fixed)
        rec = {
            "benchmark": bid, "divergence": _divergence(bid) or man.get("species", ""),
            "n_fixed": len(pi_fixed), "n_adaptive": len(pi_adapt),
            "n_new": len(new), "mean_pi_new": mean_pi_new,
            "n_lost": len(lost), "common_regressed": regressed,
            "mean_pi_fixed": _mean(pi_fixed.values()),
            "mean_pi_adaptive": _mean(pi_adapt.values()),
            "n_tagged_fixed": n_tagged_fixed, "n_tagged_adaptive": n_tagged_adapt,
            "n_redundant": n_redundant,
            "completeness_gain": round(len(new) / (man.get("n_subset_proteins") or len(pi_fixed) or 1), 5),
            "validity": {"fixed": val_fixed, "adaptive": val_adapt},
            "gate_pass": gate,
        }
        results.append(rec)
        print(f"  [{bid}] n_new={len(new)} mean_pi_new={mean_pi_new} n_lost={len(lost)} "
              f"regressed={regressed} tagged {n_tagged_fixed}->{n_tagged_adapt} "
              f"redundant={n_redundant} val {val_fixed}->{val_adapt} "
              f"gate {'PASS' if gate else 'FAIL'}", flush=True)
        (HERE / "adaptive_floor_ab.json").write_text(json.dumps(results, indent=2))
        _write_md(results, HERE / "adaptive_floor_ab.md")

    (HERE / "adaptive_floor_ab.json").write_text(json.dumps(results, indent=2))
    _write_md(results, HERE / "adaptive_floor_ab.md")
    n_pass = sum(1 for r in results if r["gate_pass"])
    print(f"\nGATE: {n_pass}/{len(results)} cells pass.", flush=True)
    print("DONE: wrote adaptive_floor_ab.json + .md", flush=True)
    return 0


def _write_md(results, path: Path):
    lines = [
        "## Divergence-adaptive rescue-floor A/B (Phase 2) — adaptive vs fixed-0.5\n",
        "Both arms run `--miniprot-rescue`; `fixed` = `LIFTON_RESCUE_ADAPTIVE_FLOOR=0` "
        "(fixed 0.5), `adaptive` = floor in [0.30, 0.50] by DNA-lift recall. Cached "
        "`-L`/`-M`, `-t1` no-`copies`, scored by the neutral re-alignment evaluator. "
        "Adaptive only LOWERS the floor within an already-additive pass, so it is a "
        "strict superset of fixed: **n_lost** and **n_redundant** MUST be 0. "
        "**n_new** = the extra genuinely-missing genes the lower floor recovers "
        "(divergent cells); high-recall cells stay at 0.5 -> **inert** (n_new=0). "
        "Gate: n_lost=0 AND n_redundant=0 AND 0 regressed AND validity not worse.\n",
        "| Dataset | divergence | n_new | mean PI new | Δcompl | n_lost | regr | tagged f→a | redundant | val f→a | gate |",
        "|---|---|---|---|---|---|---|---|---|---|---|",
    ]
    for r in results:
        v = r["validity"]
        lines.append("| {} | {} | {} | {} | {} | {} | {} | {}→{} | {} | {}→{} | {} |".format(
            r["benchmark"], r["divergence"], r["n_new"], r["mean_pi_new"],
            r["completeness_gain"], r["n_lost"], r["common_regressed"],
            r["n_tagged_fixed"], r["n_tagged_adaptive"], r["n_redundant"],
            v["fixed"], v["adaptive"],
            "**PASS**" if r["gate_pass"] else "NO-GO"))
    path.write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    sys.exit(main())
