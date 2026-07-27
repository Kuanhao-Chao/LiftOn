#!/usr/bin/env python
"""Phase-3 structural re-score: re-evaluate existing A/B output GFF3s with the new
structural metrics (evaluator.py: orf_valid / intron_chain_exact / exon Sn·Sp) to
confirm the Phase-1 (candidate-3) and Phase-2 (adaptive-floor) accuracy/recall
wins carry NO hidden structural regression that mean protein-identity masks
(exactly the failure mode the Phase-2 n_lost=1 swap showed mean-PI can hide).

Reuses the existing A/B output GFF3s ON DISK (no re-lift); only re-runs the
neutral evaluator (parasail) with the new structural columns.

The correct comparison is COMMON-SET, not a mean over all recovered: an ADDITIVE
change (candidate-3 swaps a model; adaptive-floor adds genes) grows or perturbs
the recovered set, so a mean over all recovered drops purely from denominator
DILUTION by structurally-weaker newly-recovered divergent genes -- which is not a
regression. So per cell we report:
  * COMMON set (ref transcripts recovered by BOTH arms) -- the real regression
    check: structure here must be unchanged (adaptive never swaps) or net
    non-negative (candidate-3 swaps a model only when strictly better).
  * ADDED set (changed − base) -- informational: the honest structural quality of
    the newly-recovered genes (expected weaker; that is the recall/precision trade).

Gate per cell: on the COMMON set every structural metric delta (changed − base)
>= −EPS. The added-set structure is reported, never gated.

Usage (repo root, lifton_devel env):
    python -m benchmarks.compare.structural_rescore cand3 [IDS...]
    python -m benchmarks.compare.structural_rescore adaptive [IDS...]
"""
from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

from . import evaluator

HERE = Path(__file__).resolve().parent
WORK = HERE / "work"
EPS = 0.005   # structural fracs are noisy at the tail; allow a small tolerance

AB = {
    "cand3":    {"dir": "_cand3_ladder_ab", "base": "off", "changed": "on",
                 "ids": ["t4_human_to_chicken", "human_to_mouse", "rice_to_sorghum",
                         "celegans_to_briggsae"]},
    "adaptive": {"dir": "_adaptive_floor_ab", "base": "fixed", "changed": "adaptive",
                 "ids": ["drosophila_to_anopheles", "zebrafish_to_medaka",
                         "t4_human_to_xenopus", "celegans_to_briggsae"]},
}

KEYS = ["orf_valid", "intron_chain_exact", "exon_sn", "exon_sp"]


def _load_structural(tsv: Path) -> dict:
    """{ref_mrna_id: {key: float}} for recovered coding rows with structural data."""
    out = {}
    with tsv.open() as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if str(row.get("is_coding")).strip() not in ("1", "True", "true"):
                continue
            if str(row.get("recovered")).strip() not in ("1", "True", "true"):
                continue
            if row.get("orf_valid") in ("", None):
                continue
            try:
                out[row["ref_mrna_id"]] = {k: float(row[k]) for k in KEYS}
            except (ValueError, KeyError):
                continue
    return out


def _mean(vals):
    vals = [v for v in vals if v is not None]
    return round(sum(vals) / len(vals), 5) if vals else None


def _rescore_arm(bid, cfg, arm, p, man, ref, ref_index):
    out = WORK / bid / cfg["dir"] / arm / f"{arm}.gff3"
    if not out.exists():
        raise FileNotFoundError(f"{bid}: missing arm output {out}")
    eval_dir = WORK / bid / cfg["dir"] / "eval_structural"
    eval_dir.mkdir(parents=True, exist_ok=True)
    summ = evaluator.evaluate_tool(arm, str(out), p["tgt_fa"], ref, man, eval_dir,
                                   None, log=lambda *a, **k: None,
                                   ref_index=ref_index, threads=8)
    return (_load_structural(eval_dir / f"{arm}.transcripts.tsv"),
            summ.get("protein_identity", {}).get("mean"))


def main(argv=None):
    argv = argv if argv is not None else sys.argv[1:]
    which = argv[0] if argv else "cand3"
    cfg = AB[which]
    ids = argv[1:] or cfg["ids"]

    results = []
    for bid in ids:
        print(f"=== {bid}: structural re-score ({cfg['base']} vs {cfg['changed']}) ===",
              flush=True)
        W = WORK / bid
        man = json.loads((W / "subset" / "subset.manifest.json").read_text())
        p = man["paths"]
        ref, ref_index = evaluator.build_reference(p["ref_gff"], p["ref_fa"],
                                                   log=lambda *a, **k: None)
        base, pi_base = _rescore_arm(bid, cfg, cfg["base"], p, man, ref, ref_index)
        chg, pi_chg = _rescore_arm(bid, cfg, cfg["changed"], p, man, ref, ref_index)

        common = set(base) & set(chg)
        added = set(chg) - set(base)
        # COMMON-set structural deltas (the real regression check).
        common_delta = {k: (None if not common else round(
            _mean([chg[r][k] for r in common]) - _mean([base[r][k] for r in common]), 5))
            for k in KEYS}
        # ADDED-set structure (informational).
        added_struct = {k: (_mean([chg[r][k] for r in added]) if added else None)
                        for k in KEYS}
        ok = all(d is None or d >= -EPS for d in common_delta.values())
        rec = {"benchmark": bid, "which": which,
               "n_common": len(common), "n_added": len(added),
               "mean_pi": {"base": pi_base, "changed": pi_chg},
               "common_structural_delta": common_delta,
               "added_structural": added_struct,
               "no_common_regression": ok}
        results.append(rec)
        print(f"  [{bid}] common={len(common)} added={len(added)} | "
              f"COMMON Δ {common_delta} | ADDED {added_struct} "
              f"-> {'OK' if ok else 'COMMON REGRESSION'}", flush=True)
        (HERE / f"structural_rescore_{which}.json").write_text(json.dumps(results, indent=2))
        _write_md(results, which, cfg, HERE / f"structural_rescore_{which}.md")

    n_ok = sum(1 for r in results if r["no_common_regression"])
    print(f"\nSTRUCTURAL (common-set): {n_ok}/{len(results)} cells clean.", flush=True)
    return 0


def _write_md(results, which, cfg, path):
    lines = [
        f"## Phase-3 structural re-score — {which} ({cfg['base']} vs {cfg['changed']})\n",
        "Re-evaluates the existing A/B output GFF3s with the coordinate-independent "
        "structural metrics (ORF-validity, intron-chain exactness, exon Sn·Sp). "
        "**COMMON Δ** (changed − base on ref transcripts recovered by BOTH arms) is "
        "the real regression check — it must be >= −{:.3f}. **ADDED** = the mean "
        "structure of the newly-recovered genes (informational: expected weaker — "
        "that is the honest recall/precision trade, NOT a regression).\n".format(EPS),
        "| Dataset | common | added | COMMON Δ orf/intron/sn/sp | ADDED orf/intron/sn/sp | verdict |",
        "|---|---|---|---|---|---|",
    ]
    for r in results:
        d = r["common_structural_delta"]; a = r["added_structural"]
        cd = "/".join(str(d[k]) for k in KEYS)
        ad = "/".join(str(a[k]) for k in KEYS)
        lines.append("| {} | {} | {} | {} | {} | {} |".format(
            r["benchmark"], r["n_common"], r["n_added"], cd, ad,
            "**OK**" if r["no_common_regression"] else "REGRESSION"))
    path.write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    sys.exit(main())
