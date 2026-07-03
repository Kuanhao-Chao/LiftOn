#!/usr/bin/env python
"""candidate-3 multi-cell A/B (Phase 1 gate) — divergence-ladder, cached -L/-M.

Extends the single-cell full-genome human->chicken deep-dive
(miniprot_candidate_ab.py) to the whole divergence ladder using the cached
subset -L/-M inputs (no external aligner). Candidate-3 is the 3rd
best-of-outcome candidate in run_liftoff.process_liftoff_with_protein: for a
transcript where the DNA lift is imperfect and miniprot overlaps the SAME locus,
LiftOn also scores miniprot's NATIVE model and keeps it ONLY when strictly
better than max(Liftoff+ORF, chained-merge+ORF). It therefore SWAPS a
transcript's model (never adds/drops a transcript), so it is additive
per-transcript (identity non-decreasing).

Arms (both cached -L/-M, -t1, NO -copies = deterministic; only the env differs):
  "off" : LIFTON_MINIPROT_CANDIDATE=0  (the 2-way best-of-outcome merge)
  "on"  : default (3-way)

Scored by the independent re-alignment evaluator (keyed by ref_mrna_id):
  - n_lost            : ref transcripts off has but on does not (off subset of on
                        gate — MUST be 0; candidate-3 never drops a transcript).
  - common_regressed  : per-transcript PI drop off->on (MUST be 0 — the additive
                        "strictly better only" guarantee, re-verified neutrally).
  - common_improved   : per-transcript PI gain (the win; >0 on divergent cells).
  - n_tagged          : mRNA tagged LiftOn_miniprot in the on output.
  - a2a vs miniprot    : mean(LiftOn - miniprot) on the common-with-miniprot set;
                        candidate-3 should SHRINK the deficit.
  - validity           : structured gff3_validator error count delta (the single
                        "Errors : N" summary line; must not worsen).

Per-cell gate: n_lost==0 AND common_regressed==0 AND validity not worse AND
(deficit vs miniprot not worse) AND (common_improved>0 OR the cell is INERT —
0 improved / 0 regressed / 0 lost / validity flat, a legitimate no-op). A clean
pass across the ladder confirms the promotion (candidate-3 is already default-ON).

Usage (repo root, lifton_devel env):
    python -m benchmarks.compare.miniprot_candidate_ladder_ab [IDS...]
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

# Same ladder as miniprot_rescue_ab.py: same-species control + distant +
# very-distant (where the DNA lift collapses and miniprot's same-locus model is
# strictly better -> candidate-3 fires).
# Ordered smaller/faster cells first (early harness validation via the
# incremental JSON write); the large human chromosomes + drosophila last.
DEFAULT_IDS = [
    "celegans_to_briggsae",       # distant
    "drosophila_to_anopheles",    # very-distant
    "zebrafish_to_medaka",        # very-distant
    "rice_to_sorghum",            # distant
    "t4_human_to_xenopus",        # very-distant
    "t4_human_to_chicken",        # very-distant (83% same-locus deficit)
    "human_to_mouse",             # distant
    "drosophila",                 # same-species control (slowest cell)
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
    # DETERMINISTIC: -t1, NO -copies. candidate-3 is a per-transcript swap in
    # Step 7, so off vs on differ ONLY by the env flag (byte-identity of the
    # rest is pinned by the 24-cell matrix + the drosophila cmp).
    argv = [TOOLS["lifton_bin"], "-t", "1", "-ad", _ann_db(bid),
            "-g", p["ref_gff"], "-L", str(liftoff), "-M", str(miniprot),
            "-o", str(out), p["tgt_fa"], p["ref_fa"]]
    env = dict(_compose_env(TOOLS))
    if arm == "off":
        env["LIFTON_MINIPROT_CANDIDATE"] = "0"
    else:
        env.pop("LIFTON_MINIPROT_CANDIDATE", None)
    pr = run_profiled(argv, label=f"cand3_ladder_{bid}_{arm}",
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


def _n_tagged(gff: Path) -> int:
    n = 0
    with gff.open() as fh:
        for ln in fh:
            if ln.startswith("#"):
                continue
            c = ln.split("\t")
            if len(c) > 8 and c[2] == "mRNA" and "LiftOn_miniprot" in c[8]:
                n += 1
    return n


def _validate(gff: Path):
    # Structured validator error count (the single "Errors : N" summary line),
    # NOT a crude grep of "error" substrings — the latter undercounted the
    # candidate-3 drosophila strand regression (read 1->2; true delta 0->4).
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
        print(f"=== {bid}: candidate-3 ladder A/B ===", flush=True)
        W = WORK / bid
        man = json.loads((W / "subset" / "subset.manifest.json").read_text())
        p = man["paths"]
        root = W / "_cand3_ladder_ab"
        root.mkdir(parents=True, exist_ok=True)

        out_off = _run_arm(bid, "off", p, root)
        out_on = _run_arm(bid, "on", p, root)

        ref, ref_index = evaluator.build_reference(p["ref_gff"], p["ref_fa"], log=print)
        eval_dir = root / "eval"
        eval_dir.mkdir(parents=True, exist_ok=True)
        evaluator.evaluate_tool("off", str(out_off), p["tgt_fa"], ref, man,
                                eval_dir, None, log=print, ref_index=ref_index, threads=8)
        evaluator.evaluate_tool("on", str(out_on), p["tgt_fa"], ref, man,
                                eval_dir, None, log=print, ref_index=ref_index, threads=8)
        pi_off = _pi_by_ref(eval_dir / "off.transcripts.tsv")
        pi_on = _pi_by_ref(eval_dir / "on.transcripts.tsv")

        # standalone miniprot for the apples-to-apples deficit (best-effort).
        pi_mini = {}
        try:
            mgff = WORK / bid / "tools" / "miniprot" / "miniprot.gff3"
            evaluator.evaluate_tool("miniprot", str(mgff), p["tgt_fa"], ref, man,
                                    eval_dir, None, log=print, ref_index=ref_index, threads=8)
            pi_mini = _pi_by_ref(eval_dir / "miniprot.transcripts.tsv")
        except Exception as e:
            print(f"  [warn] miniprot eval skipped: {e}", flush=True)

        common = set(pi_off) & set(pi_on)
        improved = sum(1 for k in common if pi_on[k] > pi_off[k] + 1e-9)
        regressed = sum(1 for k in common if pi_on[k] < pi_off[k] - 1e-9)
        lost = set(pi_off) - set(pi_on)
        cm_off = set(pi_off) & set(pi_mini)
        cm_on = set(pi_on) & set(pi_mini)
        a2a_off = _mean([pi_off[k] - pi_mini[k] for k in cm_off]) if pi_mini else None
        a2a_on = _mean([pi_on[k] - pi_mini[k] for k in cm_on]) if pi_mini else None
        val_off = _validate(out_off)
        val_on = _validate(out_on)

        deficit_ok = (a2a_off is None or a2a_on is None or a2a_on >= a2a_off - 1e-9)
        # An INERT cell (candidate-3 found nothing strictly better: 0 improved /
        # 0 regressed / 0 lost / validity flat) is a legitimate PASS — the flag is
        # a no-op there, not a failure. Only a real regression (regressed>0, a lost
        # transcript, worse validity, or a widened deficit) fails.
        inert = (improved == 0 and regressed == 0 and len(lost) == 0
                 and val_on == val_off)
        gate = bool(len(lost) == 0 and regressed == 0
                    and val_on <= val_off and deficit_ok
                    and (improved > 0 or inert))
        rec = {
            "benchmark": bid,
            "divergence": _divergence(bid) or man.get("species", ""),
            "n_off": len(pi_off), "n_on": len(pi_on), "n_common": len(common),
            "common_improved": improved, "common_regressed": regressed,
            "n_lost": len(lost),
            "mean_pi_off": _mean(pi_off.values()), "mean_pi_on": _mean(pi_on.values()),
            "n_tagged_LiftOn_miniprot": _n_tagged(out_on),
            "a2a_vs_miniprot_off": a2a_off, "a2a_vs_miniprot_on": a2a_on,
            "deficit_shrunk_or_flat": deficit_ok,
            "validity": {"off": val_off, "on": val_on},
            "gate_pass": gate,
        }
        results.append(rec)
        print(f"  [{bid}] improved={improved} regressed={regressed} lost={len(lost)} "
              f"tagged={rec['n_tagged_LiftOn_miniprot']} a2a {a2a_off}->{a2a_on} "
              f"val {val_off}->{val_on} gate {'PASS' if gate else 'FAIL'}", flush=True)
        # Incremental write so a mid-run failure preserves completed cells.
        (HERE / "miniprot_candidate_ladder_ab.json").write_text(json.dumps(results, indent=2))
        _write_md(results, HERE / "miniprot_candidate_ladder_ab.md")

    (HERE / "miniprot_candidate_ladder_ab.json").write_text(json.dumps(results, indent=2))
    _write_md(results, HERE / "miniprot_candidate_ladder_ab.md")
    n_pass = sum(1 for r in results if r["gate_pass"])
    print(f"\nGATE: {n_pass}/{len(results)} cells pass.", flush=True)
    print("DONE: wrote miniprot_candidate_ladder_ab.json + .md", flush=True)
    return 0


def _write_md(results, path: Path):
    lines = [
        "## candidate-3 divergence-ladder A/B (Phase 1)\n",
        "`off` = `LIFTON_MINIPROT_CANDIDATE=0` (2-way merge); `on` = default "
        "(3-way, miniprot-only candidate). Cached `-L`/`-M`, `-t1` no-`copies`, "
        "scored by the independent re-alignment evaluator (keyed by ref_mrna_id). "
        "candidate-3 SWAPS a transcript's model (never adds/drops), so **n_lost** "
        "and **common_regressed** MUST be 0 (the additive guarantee, re-verified "
        "neutrally); **common_improved** is the win; **a2a vs miniprot** is the "
        "mean(LiftOn − miniprot) deficit, which should shrink. **validity** is the "
        "structured `gff3_validator` error count (the `Errors : N` summary line). "
        "Gate: n_lost=0 AND 0 regressed AND validity not worse AND deficit not worse "
        "AND (improved>0 OR the cell is inert — a legitimate no-op).\n",
        "| Dataset | divergence | improved | regressed | n_lost | tagged | mean PI off→on | a2a off→on | val off→on | gate |",
        "|---|---|---|---|---|---|---|---|---|---|",
    ]
    for r in results:
        v = r["validity"]
        lines.append("| {} | {} | {} | {} | {} | {} | {}→{} | {}→{} | {}→{} | {} |".format(
            r["benchmark"], r["divergence"], r["common_improved"],
            r["common_regressed"], r["n_lost"], r["n_tagged_LiftOn_miniprot"],
            r["mean_pi_off"], r["mean_pi_on"], r["a2a_vs_miniprot_off"],
            r["a2a_vs_miniprot_on"], v["off"], v["on"],
            "**PASS**" if r["gate_pass"] else "NO-GO"))
    path.write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    sys.exit(main())
