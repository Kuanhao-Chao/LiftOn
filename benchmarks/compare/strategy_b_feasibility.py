#!/usr/bin/env python
"""Phase-4 feasibility probe: incremental headroom of a Strategy-B splice-junction
rebuild OVER what candidate-3 (Phase 1) already captures. FEASIBILITY-FIRST: this
is a byte-neutral measurement (no lifting, no lifton/ change) that gates whether
Strategy B is worth building at all.

Strategy B (lifton2-proven) rebuilds a transcript's exon/CDS structure from
miniprot's OWN splice junctions -- a HYBRID of miniprot CDS junctions + DNA-lift
UTR -- keeping it only when strictly better vs the reference protein. candidate-3
(already shipped, v1.0.10 P1) emits miniprot's NATIVE CDS-only model when its
ORF-rescued identity strictly beats the 2-way merge max(Liftoff+ORF, merge+ORF).

The key question: does Strategy B's hybrid add PROTEIN identity over candidate-3?
Both use miniprot's CDS junctions, so they translate the SAME protein -- the only
structural difference (DNA-lift UTR) is untranslated. And candidate-3's decision
bar (the 2-way merge) is STRONGER than Strategy B's (a pure DNA chain), so
candidate-3 already captures >= what Strategy B would. So the incremental protein
headroom should be ~0 -> candidate-3 subsumes Strategy B on the accuracy axis.

This probe MEASURES that (reusing the candidate-3 ladder A/B's `on` + standalone
`miniprot` TSVs, keyed by ref_mrna_id, scored by the neutral evaluator):
  * subsumption check -- for the transcripts candidate-3 FIRED on (tagged
    LiftOn_miniprot), is pi(lifton_on) already ~= pi(miniprot)?  If yes,
    candidate-3 emits exactly miniprot's model -> Strategy B's hybrid (same CDS)
    adds 0 protein there.
  * oracle best-of headroom -- mean over the common coding set of
    max(0, pi_miniprot - pi_lifton_on): the UPPER bound a perfect
    best-of-{lifton,miniprot} selection could reach. This is context; it is a
    candidate-3 FIRING-CRITERION opportunity (internal-ORF vs neutral metric),
    NOT a Strategy-B structural one (Strategy B gives the same CDS).

Gate: Strategy-B incremental (the tagged-set mean pi_miniprot - pi_lifton_on)
>= +0.003 => headroom, worth building; < +0.003 => NO-GO, candidate-3 subsumes it.

Usage (repo root, lifton_devel env):
    python -m benchmarks.compare.strategy_b_feasibility [IDS...]
"""
from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
WORK = HERE / "work"
DEFAULT_IDS = ["celegans_to_briggsae", "drosophila_to_anopheles",
               "zebrafish_to_medaka", "rice_to_sorghum", "t4_human_to_chicken"]
GO_BAR = 0.003


def _pi_by_ref(tsv: Path) -> dict:
    out = {}
    if not tsv.exists():
        return out
    with tsv.open() as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if str(row.get("is_coding")).strip() not in ("1", "True", "true"):
                continue
            pi = row.get("protein_identity")
            if pi in (None, "", "None"):
                continue
            try:
                v = float(pi)
            except ValueError:
                continue
            rid = row["ref_mrna_id"]
            out[rid] = max(out.get(rid, v), v)
    return out


def _tagged_ref_ids(on_gff: Path, on_tsv: Path) -> set:
    """ref_mrna_ids whose emitted model is a candidate-3 win (tool mRNA tagged
    status=LiftOn_miniprot). Map tagged tool-feature-id -> ref via on.tsv."""
    tagged_tool = set()
    with on_gff.open() as fh:
        for ln in fh:
            if ln.startswith("#"):
                continue
            c = ln.split("\t")
            if len(c) > 8 and c[2] == "mRNA" and "status=LiftOn_miniprot" in c[8]:
                for kv in c[8].split(";"):
                    if kv.startswith("ID="):
                        tagged_tool.add(kv[3:].strip())
    ref_ids = set()
    with on_tsv.open() as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row.get("tool_feature_id") in tagged_tool:
                ref_ids.add(row["ref_mrna_id"])
    return ref_ids


def _mean(xs):
    xs = [x for x in xs if x is not None]
    return round(sum(xs) / len(xs), 5) if xs else None


def main(argv=None):
    ids = (argv if argv is not None else sys.argv[1:]) or DEFAULT_IDS
    results = []
    for bid in ids:
        d = WORK / bid / "_cand3_ladder_ab"
        pi_on = _pi_by_ref(d / "eval" / "on.transcripts.tsv")
        pi_mp = _pi_by_ref(d / "eval" / "miniprot.transcripts.tsv")
        tagged = _tagged_ref_ids(d / "on" / "on.gff3", d / "eval" / "on.transcripts.tsv")

        common = set(pi_on) & set(pi_mp)
        # Strategy-B incremental over candidate-3: on the transcripts candidate-3
        # fired (tagged), miniprot vs the emitted lifton model. Strategy B uses the
        # SAME miniprot CDS -> same protein, so this measures what (if anything) a
        # hybrid rebuild would add beyond candidate-3's already-emitted model.
        tagged_common = tagged & common
        strat_b_incr = _mean([pi_mp[r] - pi_on[r] for r in tagged_common])
        # oracle best-of upper bound (context; a firing-criterion lever, not Strat-B).
        oracle_headroom = _mean([max(0.0, pi_mp[r] - pi_on[r]) for r in common])
        # how close candidate-3's tagged emissions already are to miniprot.
        tagged_abs_gap = _mean([abs(pi_on[r] - pi_mp[r]) for r in tagged_common])

        go = bool(strat_b_incr is not None and strat_b_incr >= GO_BAR)
        rec = {
            "benchmark": bid, "n_common": len(common), "n_tagged": len(tagged),
            "n_tagged_common": len(tagged_common),
            "strategy_b_incremental_pi": strat_b_incr,
            "tagged_mean_abs_gap_vs_miniprot": tagged_abs_gap,
            "oracle_bestof_headroom": oracle_headroom,
            "go": go,
        }
        results.append(rec)
        print(f"  [{bid}] tagged_common={len(tagged_common)} "
              f"strat_b_incr={strat_b_incr} (bar +{GO_BAR}) "
              f"tagged|gap|vs_mp={tagged_abs_gap} oracle_headroom={oracle_headroom} "
              f"-> {'GO' if go else 'NO-GO'}", flush=True)
        (HERE / "strategy_b_feasibility.json").write_text(json.dumps(results, indent=2))
        _write_md(results, HERE / "strategy_b_feasibility.md")

    n_go = sum(1 for r in results if r["go"])
    print(f"\nFEASIBILITY: Strategy-B incremental >= +{GO_BAR} on {n_go}/{len(results)} cells.",
          flush=True)
    if n_go == 0:
        print("VERDICT: NO-GO -- candidate-3 subsumes Strategy B on the protein axis "
              "(same miniprot CDS -> same protein; candidate-3's 2-way-merge bar is "
              "stronger than Strategy B's DNA-chain bar). Residual oracle headroom is a "
              "candidate-3 firing-criterion lever (internal-ORF vs neutral metric), NOT "
              "a Strategy-B structural one.", flush=True)
    return 0


def _write_md(results, path: Path):
    lines = [
        "## Phase-4 Strategy-B feasibility probe — incremental headroom OVER candidate-3\n",
        "Byte-neutral (reuses the candidate-3 ladder A/B `on` + standalone `miniprot` "
        "TSVs, neutral re-scorer). **strategy_b_incremental_pi** = mean(pi_miniprot − "
        "pi_lifton_on) on the transcripts candidate-3 FIRED (tagged LiftOn_miniprot): "
        "Strategy B's hybrid uses the SAME miniprot CDS, so this is the protein it would "
        "add beyond candidate-3's emitted model. **tagged |gap|** = how close "
        "candidate-3's tagged emissions already are to miniprot (≈0 => candidate-3 emits "
        "miniprot's model). **oracle headroom** = mean(max(0, pi_miniprot − pi_lifton_on)) "
        "over the common set — the upper bound of a perfect best-of selection (a "
        "firing-criterion lever, NOT Strategy B). Gate: strategy_b_incremental_pi >= +0.003.\n",
        "| Dataset | tagged∩common | Strat-B incr PI | tagged \\|gap\\| vs mp | oracle headroom | verdict |",
        "|---|---|---|---|---|---|",
    ]
    for r in results:
        lines.append("| {} | {} | {} | {} | {} | {} |".format(
            r["benchmark"], r["n_tagged_common"], r["strategy_b_incremental_pi"],
            r["tagged_mean_abs_gap_vs_miniprot"], r["oracle_bestof_headroom"],
            "GO" if r["go"] else "**NO-GO**"))
    path.write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    sys.exit(main())
