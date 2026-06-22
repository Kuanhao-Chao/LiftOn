#!/usr/bin/env python
"""Self-contained 2-state A/B for --orf-best-match (Iteration 9).

AUDIT RECORD (Iteration-9 NO-GO). The `--orf-best-match` / `LIFTON_ORF_BEST_MATCH`
instrumentation in `lifton_class.__find_orfs` that this A/B drove was implemented,
measured, and **REVERTED** as net-negative (drosophila −0.00003, 4 improved / 59
regressed; mouse_to_rat +0.00008, 3 / 23 — the strip/re-score ripples into the
best-of-outcome merge). This script + its committed `orf_best_match_ab.{json,md}`
are retained as the audit trail; re-running requires restoring the reverted
instrumentation. See the CLAUDE.md Iteration-9 scope clarification + memory.


ORF-rescue (`lifton_class.__find_orfs`) legacy behaviour keeps only the LONGEST
ORF per reading frame, then identity-scores those <=3 survivors. `--orf-best-match`
keeps the top-K(=3) longest per frame and selects by reference-protein identity,
so a shorter but better-matching ORF is no longer discarded. Output-changing
ONLY on the divergent transcripts where ORF-rescue fires (stop-gain / frameshift
/ start-lost).

  state "default"   (no flag)                  legacy longest-per-frame  ← default
  state "orf_best"  LIFTON_ORF_BEST_MATCH=1     best-match top-K          ← candidate

Both states run on the same cached `-L`/`-M` (so Step 4 is skipped and the only
difference is the ORF-rescue selection) and are scored by the SAME independent
re-alignment evaluator. Per-transcript protein identity is compared.

Decision gate (consistent with the best-of-outcome promotion, +0.0037 / +0.0052):
  PROMOTE iff mean-PI delta >= +0.003 on BOTH datasets AND completeness identical
  AND regressed <= 0.1 * improved. Else keep opt-in (default OFF), document NO-GO.

Usage (repo root, lifton_devel env):
    python -m benchmarks.compare.orf_best_match_ab [IDS...]   # default below
"""
from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

from . import evaluator
from .profiling import run_profiled
from .tool_runners import _clean_input_dbs, _compose_env
from .align_window_ab import _TRANSCRIPT_TYPES

HERE = Path(__file__).resolve().parent
WORK = HERE / "work"
REG = json.loads((HERE / "benchmarks.json").read_text())
TOOLS = REG["tools"]

DEFAULT_IDS = ["drosophila", "mouse_to_rat"]
STATES = ("default", "orf_best")
PROMOTE_BAR = 0.003


def _ann_db(bid: str) -> str:
    for b in REG["benchmarks"]:
        if b["id"] == bid:
            return b.get("annotation_database", "RefSeq")
    return "RefSeq"


def _run_state(bid: str, p: dict, state: str, ab_root: Path, log=print):
    liftoff_gff = WORK / bid / "tools" / "liftoff" / "liftoff.gff3"
    miniprot_gff = WORK / bid / "tools" / "miniprot" / "miniprot.gff3"
    statedir = ab_root / state
    statedir.mkdir(parents=True, exist_ok=True)
    out = statedir / f"{state}.gff3"
    argv = [TOOLS["lifton_bin"], "-t", "8", "-copies", "-ad", _ann_db(bid),
            "-g", p["ref_gff"], "-L", str(liftoff_gff), "-M", str(miniprot_gff),
            "-o", str(out), "--native", "--locus-pipeline",
            p["tgt_fa"], p["ref_fa"]]
    env = _compose_env(TOOLS)
    if state == "orf_best":
        env["LIFTON_ORF_BEST_MATCH"] = "1"   # env-flip, like LIFTON_FULL_DP_ALIGN
    pr = run_profiled(argv, label=f"orf_best_ab_{state}",
                      log_dir=ab_root / "logs", env=env, log=log)
    if pr.exit_code != 0 or not out.exists() or out.stat().st_size == 0:
        raise RuntimeError(
            f"{bid}: state {state} lifton failed (exit {pr.exit_code}); "
            f"see {pr.stderr_path}")
    return out, pr


def _pi_by_ref(tsv: Path) -> dict:
    out = {}
    if not tsv.exists():
        return out
    with tsv.open() as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row.get("is_coding") not in ("1", "True", "true"):
                continue
            pi = row.get("protein_identity")
            if pi in (None, "", "None"):
                continue
            try:
                out[row["ref_mrna_id"]] = float(pi)
            except ValueError:
                pass
    return out


def _feature_counts(gff_path: Path) -> dict:
    g = t = c = tot = 0
    for ln in Path(gff_path).read_text().splitlines():
        if not ln or ln.startswith("#"):
            continue
        cols = ln.split("\t")
        if len(cols) < 9:
            continue
        tot += 1
        ft = cols[2]
        if ft == "gene" or ft.endswith("_gene") or ft == "pseudogene":
            g += 1
        elif ft in _TRANSCRIPT_TYPES:
            t += 1
        elif ft == "CDS":
            c += 1
    return {"gene": g, "transcript": t, "CDS": c, "total": tot}


def _mean(vals):
    vals = [v for v in vals if v is not None]
    return round(sum(vals) / len(vals), 5) if vals else None


def _n_complete(pi: dict) -> int:
    return sum(1 for v in pi.values() if v >= 1.0 - 1e-9)


def _delta_counts(base: dict, new: dict):
    common = set(base) & set(new)
    improved = sum(1 for k in common if new[k] > base[k] + 1e-9)
    regressed = sum(1 for k in common if new[k] < base[k] - 1e-9)
    net = round(sum(new[k] - base[k] for k in common) / len(common), 6) if common else 0.0
    return {"n_common": len(common), "improved": improved, "regressed": regressed,
            "net_per_transcript": net}


def main(argv=None):
    ids = (argv if argv is not None else sys.argv[1:]) or DEFAULT_IDS
    results = []
    for bid in ids:
        print(f"=== {bid} (default vs --orf-best-match) ===", flush=True)
        W = WORK / bid
        man = json.loads((W / "subset" / "subset.manifest.json").read_text())
        p = man["paths"]
        ab_root = W / "_orf_best_ab"
        ab_root.mkdir(parents=True, exist_ok=True)

        liftoff_gff = W / "tools" / "liftoff" / "liftoff.gff3"
        miniprot_gff = W / "tools" / "miniprot" / "miniprot.gff3"
        _clean_input_dbs(p["ref_gff"], liftoff_gff, miniprot_gff)

        gffs, prs = {}, {}
        for state in STATES:
            print(f"--- {bid} state {state} ---", flush=True)
            gffs[state], prs[state] = _run_state(bid, p, state, ab_root, log=print)

        ref, ref_index = evaluator.build_reference(p["ref_gff"], p["ref_fa"], log=print)
        eval_dir = ab_root / "eval"
        eval_dir.mkdir(parents=True, exist_ok=True)
        for state in STATES:
            evaluator.evaluate_tool(state, str(gffs[state]), p["tgt_fa"],
                                    ref, man, eval_dir, None, log=print,
                                    ref_index=ref_index, threads=8)

        pi = {s: _pi_by_ref(eval_dir / f"{s}.transcripts.tsv") for s in STATES}
        counts = {s: _feature_counts(gffs[s]) for s in STATES}
        delta = _delta_counts(pi["default"], pi["orf_best"])
        mean_default = _mean(pi["default"].values())
        mean_orf = _mean(pi["orf_best"].values())
        mean_delta = (round(mean_orf - mean_default, 6)
                      if (mean_orf is not None and mean_default is not None) else None)
        completeness_ok = counts["default"] == counts["orf_best"]
        gate_pass = (mean_delta is not None and mean_delta >= PROMOTE_BAR
                     and completeness_ok
                     and delta["regressed"] <= 0.1 * max(delta["improved"], 1))

        rec = {
            "benchmark": bid,
            "n_scored": {s: len(pi[s]) for s in STATES},
            "n_complete": {s: _n_complete(pi[s]) for s in STATES},
            "mean_pi": {"default": mean_default, "orf_best": mean_orf},
            "mean_pi_delta": mean_delta,
            "delta_counts": delta,
            "completeness": {"counts": counts, "unchanged": completeness_ok},
            "wall_s": {s: round(prs[s].wall_clock_seconds, 2) for s in STATES},
            "gate": {"bar": PROMOTE_BAR, "mean_delta_ge_bar":
                     bool(mean_delta is not None and mean_delta >= PROMOTE_BAR),
                     "completeness_unchanged": completeness_ok,
                     "regressed_le_10pct_improved":
                     bool(delta["regressed"] <= 0.1 * max(delta["improved"], 1)),
                     "pass": bool(gate_pass)},
        }
        results.append(rec)
        print(f"  [{bid}] mean PI: default={mean_default} orf_best={mean_orf} "
              f"(Δ {mean_delta:+.6f}); {delta['improved']} improved / "
              f"{delta['regressed']} regressed; completeness "
              f"{'OK' if completeness_ok else 'CHANGED'}; gate "
              f"{'PASS' if gate_pass else 'FAIL'}", flush=True)

    (HERE / "orf_best_match_ab.json").write_text(json.dumps(results, indent=2))
    _write_md(results, HERE / "orf_best_match_ab.md")
    n_pass = sum(1 for r in results if r["gate"]["pass"])
    print(f"\nGATE: {n_pass}/{len(results)} datasets pass (need ALL for promotion).",
          flush=True)
    print("DONE: wrote orf_best_match_ab.json + orf_best_match_ab.md", flush=True)
    return 0


def _write_md(results, path: Path):
    lines = [
        "## --orf-best-match A/B (ORF-rescue best-matching ORF) — Iteration 9\n",
        "Both states run on the same cached `-L`/`-M` (Step 4 skipped → only the "
        "ORF-rescue selection differs), scored by the independent re-alignment "
        "evaluator. **default** = legacy longest-per-frame; **orf_best** = "
        "`LIFTON_ORF_BEST_MATCH=1` top-K best-match. Gate: mean-PI Δ ≥ "
        f"+{PROMOTE_BAR} on BOTH datasets AND completeness identical AND "
        "regressed ≤ 0.1×improved.\n",
        "| Dataset | default mean PI | orf_best mean PI | mean Δ | improved/regressed | complete | wall d/o | gate |",
        "|---|---|---|---|---|---|---|---|",
    ]
    for r in results:
        m, d, w = r["mean_pi"], r["delta_counts"], r["wall_s"]
        lines.append("| {} | {} | {} | {:+.6f} | {}/{} | {} | {}/{} | {} |".format(
            r["benchmark"], m["default"], m["orf_best"],
            r["mean_pi_delta"] if r["mean_pi_delta"] is not None else 0.0,
            d["improved"], d["regressed"],
            "yes" if r["completeness"]["unchanged"] else "**NO**",
            w["default"], w["orf_best"],
            "**PASS**" if r["gate"]["pass"] else "NO-GO"))
    path.write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    sys.exit(main())
