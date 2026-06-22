#!/usr/bin/env python
"""Iteration-13 feasibility diagnostic — Step-8 miniprot-rescue headroom.

AUDIT RECORD (Iteration-13 NO-GO). The byte-neutral `LIFTON_STEP8_DIAG`
instrumentation this driver relies on (in `run_miniprot.process_miniprot` +
`lifton_utils.check_ovps_ratio(return_max=True)`) was implemented, measured, and
**REVERTED** — the Step-8 overlap/length filters are well-tuned. This diagnostic
found negligible near-miss headroom (drosophila 5 absolute candidates,
mouse_to_rat 0), and the follow-on A/B (`step8_rescue_ab.py`) confirmed on
emitted output that relaxing the thresholds adds only low-identity garbage
(drosophila +2 transcripts at mean PI 0.35; mouse_to_rat +0). This script + its
committed `step8_rescue_headroom.{json,md}` are the audit trail; re-running
requires restoring the reverted instrumentation. See the CLAUDE.md Iteration-13
note + memory.

Read-only / BYTE-NEUTRAL. Step 8 (`lifton/run_miniprot.process_miniprot`) emits a
miniprot-only mRNA only if it passes the overlap (`-overlap` default 0.1) and
length-ratio (`-min_miniprot`/`-max_miniprot` 0.9/1.5) filters. Those thresholds
are already CLI args, so the open question is purely empirical: do the filters
drop genuinely-missing transcripts that relaxing the defaults would recover?

This runs LiftOn ONCE per dataset with `LIFTON_STEP8_DIAG=<path>` set, DEFAULT
thresholds (output byte-frozen), cached `-L`/`-M`. `process_miniprot` then writes
one TSV line per candidate: `mtid, reason, max_overlap_ratio, length_ratio,
ref_gene, ref_trans, n_cds`. We aggregate the drop-reason breakdown and — the
heart of the gate — the **near-miss headroom**: drops whose tunable metric sits
just past the threshold (overlap in (0.1,0.2]; length in [0.8,0.9) ∪ [1.5,1.6)).
Near-miss counts are an UPPER BOUND on what relaxing could recover (the real
yield, after the identity check, is lower — measured by the follow-on A/B).

GO/NO-GO (the Iter-4/9 feasibility-first recipe): if near-miss is negligible vs
emitted on BOTH datasets → NO-GO (relaxing won't materially help; the defaults
are well-tuned). If material → GO to `step8_rescue_ab.py` with relaxed thresholds
chosen from the near-miss bucket distribution.

Usage (repo root, lifton_devel env):
    python -m benchmarks.compare.step8_rescue_headroom [IDS...]   # default below
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

from .profiling import run_profiled
from .tool_runners import _clean_input_dbs, _compose_env

HERE = Path(__file__).resolve().parent
WORK = HERE / "work"
REG = json.loads((HERE / "benchmarks.json").read_text())
TOOLS = REG["tools"]

DEFAULT_IDS = ["drosophila", "mouse_to_rat"]
# Default Step-8 thresholds (for the near-miss windows).
OVERLAP_THR, MIN_MP, MAX_MP = 0.1, 0.9, 1.5


def _ann_db(bid: str) -> str:
    for b in REG["benchmarks"]:
        if b["id"] == bid:
            return b.get("annotation_database", "RefSeq")
    return "RefSeq"


def _run_diag(bid: str, p: dict, root: Path, log=print):
    liftoff = WORK / bid / "tools" / "liftoff" / "liftoff.gff3"
    miniprot = WORK / bid / "tools" / "miniprot" / "miniprot.gff3"
    out = root / "default.gff3"
    diag = root / "step8_diag.tsv"
    if diag.exists():
        diag.unlink()
    _clean_input_dbs(p["ref_gff"], liftoff, miniprot)
    env = _compose_env(TOOLS)
    env["LIFTON_STEP8_DIAG"] = str(diag)
    # DEFAULT thresholds (no -overlap/-min_miniprot/-max_miniprot) → output
    # byte-frozen; the diag only MEASURES the drop population.
    argv = [TOOLS["lifton_bin"], "-t", "8", "-copies", "-ad", _ann_db(bid),
            "-g", p["ref_gff"], "-L", str(liftoff), "-M", str(miniprot),
            "-o", str(out), "--native", "--locus-pipeline",
            p["tgt_fa"], p["ref_fa"]]
    pr = run_profiled(argv, label=f"step8_headroom_{bid}",
                      log_dir=root / "logs", env=env, log=log)
    if pr.exit_code != 0 or not out.exists() or out.stat().st_size == 0:
        raise RuntimeError(f"{bid}: lifton failed (exit {pr.exit_code}); "
                           f"see {pr.stderr_path}")
    return diag


def _parse_diag(diag: Path):
    rows = []
    if not diag.exists():
        return rows
    for line in diag.read_text().splitlines():
        c = line.split("\t")
        if len(c) < 7:
            continue

        def _f(x):
            try:
                return float(x)
            except (ValueError, TypeError):
                return None
        rows.append({"mtid": c[0], "reason": c[1], "max_ovp": _f(c[2]),
                     "len_ratio": _f(c[3]), "ref_gene": c[4], "ref_trans": c[5],
                     "n_cds": c[6]})
    return rows


def _bucket(vals, edges):
    """Count vals into half-open buckets [edges[i], edges[i+1])."""
    out = []
    for lo, hi in zip(edges, edges[1:]):
        out.append(sum(1 for v in vals if v is not None and lo <= v < hi))
    return out


def main(argv=None):
    ids = (argv if argv is not None else sys.argv[1:]) or DEFAULT_IDS
    results = []
    for bid in ids:
        print(f"=== {bid}: Step-8 rescue headroom diagnostic ===", flush=True)
        W = WORK / bid
        man = json.loads((W / "subset" / "subset.manifest.json").read_text())
        p = man["paths"]
        root = W / "_step8_headroom"
        root.mkdir(parents=True, exist_ok=True)

        diag = _run_diag(bid, p, root)
        rows = _parse_diag(diag)

        n_candidates = len(rows)
        reasons = ["overlap", "no_ref_gene", "processed_pseudogene",
                   "length_ratio", "no_ref_protein", "no_ref_trans_id", "EMITTED"]
        by_reason = {r: sum(1 for x in rows if x["reason"] == r) for r in reasons}
        n_emitted = by_reason["EMITTED"]

        # Overlap near-miss: overlap-drops with max_ovp in (0.1, 0.2].
        ovp_drops = [x["max_ovp"] for x in rows if x["reason"] == "overlap"]
        n_overlap_nearmiss = sum(1 for v in ovp_drops
                                 if v is not None and OVERLAP_THR < v <= 0.2)
        ovp_buckets = dict(zip(
            ["(0.1,0.15]", "(0.15,0.2]", "(0.2,0.3]", "(0.3,0.5]", "(0.5,1.0]"],
            _bucket([v for v in ovp_drops if v is not None and v > OVERLAP_THR],
                    [0.1, 0.15, 0.2, 0.3, 0.5, 1.0001])))

        # Length near-miss: length-drops with ratio in [0.8,0.9) ∪ [1.5,1.6).
        len_drops = [x["len_ratio"] for x in rows if x["reason"] == "length_ratio"]
        n_len_nearmiss = sum(1 for v in len_drops if v is not None and (
            0.8 <= v < MIN_MP or MAX_MP <= v < 1.6))
        len_buckets = dict(zip(
            ["<0.8", "[0.8,0.9)", "[1.5,1.6)", ">=1.6"],
            [sum(1 for v in len_drops if v is not None and v < 0.8),
             sum(1 for v in len_drops if v is not None and 0.8 <= v < MIN_MP),
             sum(1 for v in len_drops if v is not None and MAX_MP <= v < 1.6),
             sum(1 for v in len_drops if v is not None and v >= 1.6)]))

        nearmiss_total = n_overlap_nearmiss + n_len_nearmiss
        nearmiss_over_emitted = round(nearmiss_total / max(n_emitted, 1), 4)
        verdict = ("GO (worth an A/B)" if nearmiss_over_emitted >= 0.05
                   else "NO-GO / negligible headroom")

        rec = {
            "benchmark": bid,
            "n_candidates": n_candidates,
            "by_reason": by_reason,
            "n_emitted": n_emitted,
            "n_overlap_nearmiss_0.1_0.2": n_overlap_nearmiss,
            "overlap_drop_buckets": ovp_buckets,
            "n_length_nearmiss": n_len_nearmiss,
            "length_drop_buckets": len_buckets,
            "nearmiss_total": nearmiss_total,
            "nearmiss_over_emitted": nearmiss_over_emitted,
            "verdict": verdict,
        }
        results.append(rec)
        print(f"  [{bid}] candidates={n_candidates} emitted={n_emitted} "
              f"by_reason={by_reason}", flush=True)
        print(f"    overlap near-miss (0.1,0.2]={n_overlap_nearmiss} "
              f"buckets={ovp_buckets}", flush=True)
        print(f"    length near-miss={n_len_nearmiss} buckets={len_buckets}", flush=True)
        print(f"    near-miss/emitted = {nearmiss_over_emitted}  -->  {verdict}",
              flush=True)

    (HERE / "step8_rescue_headroom.json").write_text(json.dumps(results, indent=2))
    _write_md(results, HERE / "step8_rescue_headroom.md")
    any_go = any(r["nearmiss_over_emitted"] >= 0.05 for r in results)
    print(f"\nOVERALL: {'GO' if any_go else 'NO-GO'} "
          f"(GO if near-miss/emitted >= 0.05 on any dataset).", flush=True)
    print("DONE: wrote step8_rescue_headroom.json + step8_rescue_headroom.md",
          flush=True)
    return 0


def _write_md(results, path: Path):
    lines = [
        "## Step-8 miniprot-rescue headroom diagnostic (Iteration 13, read-only)\n",
        "Default run on cached `-L`/`-M`, `LIFTON_STEP8_DIAG` set (byte-frozen "
        "output). Per candidate miniprot mRNA: drop reason + tunable-filter "
        "metrics. **Near-miss** = drops just past the threshold (overlap in "
        "(0.1,0.2]; length in [0.8,0.9) ∪ [1.5,1.6)) — the UPPER BOUND on what "
        "relaxing could recover. GO if near-miss/emitted ≥ 0.05 on any dataset; "
        "else NO-GO (defaults well-tuned).\n",
        "| Dataset | candidates | emitted | overlap-drop | length-drop | overlap near-miss | length near-miss | near-miss/emitted | verdict |",
        "|---|---|---|---|---|---|---|---|---|",
    ]
    for r in results:
        br = r["by_reason"]
        lines.append("| {} | {} | {} | {} | {} | {} | {} | {} | {} |".format(
            r["benchmark"], r["n_candidates"], r["n_emitted"], br["overlap"],
            br["length_ratio"], r["n_overlap_nearmiss_0.1_0.2"],
            r["n_length_nearmiss"], r["nearmiss_over_emitted"], r["verdict"]))
    path.write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    sys.exit(main())
