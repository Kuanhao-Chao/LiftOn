#!/usr/bin/env python
"""figure4_outliers.py — read-only deep-dive into the Figure-4 low performers.

Figure 4 plots LiftOn v1.0.9's apples-to-apples per-transcript protein-identity
lead over miniprot on the COMMON recovered set. Four very-distant pairs are
negative (LiftOn < miniprot). This ranks the per-transcript outliers driving the
deficit and quantifies the additive headroom a miniprot-only merge candidate
would deliver (= mean over common of max(0, miniprot_pi - lifton_pi)).

Pure analysis over the existing per-transcript eval TSVs
(work/<pair>/_fourway_full/eval/<tool>.transcripts.tsv) — no lifton run.

Run (repo root):  python -m benchmarks.compare.figure4_outliers
Writes benchmarks/compare/figure4_outliers.{json,md} + per-pair ranked
benchmarks/compare/figure4_outliers/<pair>.tsv.
"""
from __future__ import annotations
import csv, json, os, statistics as st
from pathlib import Path

HERE = Path(__file__).resolve().parent
WORK = HERE / "work"
OUT = HERE / "figure4_outliers"
EPS = 1e-9

# negative very-distant pairs + a positive very-distant control + a close control
PAIRS = ["human_to_zebrafish", "t4_human_to_xenopus", "t4_human_to_chicken",
         "arabidopsis_to_rice", "drosophila_to_bee", "drosophila"]


def load_rows(eval_dir: Path, tool: str):
    """{ref_id: best-PI row dict} over recovered coding rows (keep the highest-PI copy)."""
    path = eval_dir / f"{tool}.transcripts.tsv"
    if not path.exists():
        return None
    best: dict[str, dict] = {}
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if str(r.get("is_coding")).strip() not in ("1", "true", "True"):
                continue
            if str(r.get("recovered")).strip() not in ("1", "true", "True"):
                continue
            try:
                pi = float(r["protein_identity"])
            except (ValueError, TypeError, KeyError):
                continue
            rid = r["ref_mrna_id"]
            if rid not in best or pi > float(best[rid]["protein_identity"] or 0):
                best[rid] = r
    return best


def fnum(r, k):
    try:
        return float(r[k])
    except (ValueError, TypeError, KeyError):
        return None


def analyze(pair: str):
    eval_dir = WORK / pair / "_fourway_full" / "eval"
    if not eval_dir.exists():
        eval_dir = WORK / pair / "_fourway" / "eval"
    dev = load_rows(eval_dir, "lifton_devel")
    mini = load_rows(eval_dir, "miniprot")
    if not dev or not mini:
        return None
    common = set(dev) & set(mini)
    rows = []
    for rid in common:
        dpi = fnum(dev[rid], "protein_identity") or 0.0
        mpi = fnum(mini[rid], "protein_identity") or 0.0
        rows.append({
            "ref_id": rid, "lifton_pi": round(dpi, 4), "miniprot_pi": round(mpi, 4),
            "delta": round(dpi - mpi, 4),
            "lifton_dna_id": fnum(dev[rid], "dna_identity"),
            "lifton_plen": fnum(dev[rid], "lifted_prot_len"),
            "ref_plen": fnum(dev[rid], "ref_prot_len"),
            "lifton_fid": dev[rid].get("tool_feature_id", ""),
        })
    rows.sort(key=lambda x: x["delta"])
    deltas = [r["delta"] for r in rows]
    reg = [r for r in rows if r["delta"] < -EPS]
    imp = [r for r in rows if r["delta"] > EPS]
    headroom = st.mean(max(0.0, -r["delta"]) for r in rows) if rows else 0.0
    # catastrophic-Liftoff signature: LiftOn near-garbage where miniprot is good
    catastrophic = [r for r in rows if r["lifton_pi"] < 0.3 and r["miniprot_pi"] > 0.7]
    # of those, how many have decent lifton DNA identity (frameshift signature)?
    cat_frameshift = [r for r in catastrophic if (r["lifton_dna_id"] or 0) >= 0.6]
    summary = {
        "pair": pair, "n_common": len(rows),
        "mean_delta": round(st.mean(deltas), 5) if deltas else None,
        "n_regressed": len(reg), "n_improved": len(imp),
        "mean_deficit_over_regressed": round(st.mean(-r["delta"] for r in reg), 4) if reg else 0,
        "fix_headroom_meanpi": round(headroom, 5),
        "n_catastrophic(LiftOn<0.3,mini>0.7)": len(catastrophic),
        "n_catastrophic_with_dna>=0.6(frameshift)": len(cat_frameshift),
        "catastrophic_share_of_regressed": round(len(catastrophic) / len(reg), 3) if reg else 0,
    }
    return summary, rows


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    summaries = []
    md = ["# Figure-4 outlier deep-dive (LiftOn v1.0.9 vs miniprot, common recovered coding set)\n"]
    md.append("`fix_headroom` = mean over common of max(0, miniprot_pi − lifton_pi) — the additive upside a "
              "miniprot-only merge candidate would deliver.\n")
    md.append("| pair | n_common | mean Δ | regressed | improved | headroom | catastrophic (LiftOn<0.3, mini>0.7) | of which frameshift (dna≥0.6) |")
    md.append("|---|---|---|---|---|---|---|---|")
    for pair in PAIRS:
        res = analyze(pair)
        if not res:
            md.append(f"| {pair} | — (no eval TSV) | | | | | | |")
            continue
        s, rows = res
        summaries.append(s)
        md.append(f"| {pair} | {s['n_common']} | {s['mean_delta']} | {s['n_regressed']} | {s['n_improved']} | "
                  f"{s['fix_headroom_meanpi']} | {s['n_catastrophic(LiftOn<0.3,mini>0.7)']} "
                  f"({s['catastrophic_share_of_regressed']:.0%} of reg) | {s['n_catastrophic_with_dna>=0.6(frameshift)']} |")
        # per-pair ranked outliers (top 40 most negative)
        with open(OUT / f"{pair}.tsv", "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
            w.writeheader()
            for r in rows[:40]:
                w.writerow(r)
    (HERE / "figure4_outliers.json").write_text(json.dumps(summaries, indent=2))
    (HERE / "figure4_outliers.md").write_text("\n".join(md) + "\n")
    print("\n".join(md))
    print(f"\nwrote figure4_outliers.{{json,md}} + per-pair ranked TSVs in {OUT}/")


if __name__ == "__main__":
    main()
