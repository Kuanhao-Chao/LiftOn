#!/usr/bin/env python
"""Iteration-9 feasibility diagnostic — measure ORF-rescue best-match headroom.

AUDIT RECORD (Iteration-9 NO-GO). The `lifton/` instrumentation this drives
(`LIFTON_ORF_DIAG` in `lifton_class.__find_orfs`) was implemented, measured, and
**REVERTED** — the broadened best-match selection was net-negative (see
`orf_best_match_ab.py` + the CLAUDE.md Iteration-9 scope clarification). This
script + its committed `orf_best_match_headroom.{json,md}` are retained as the
audit trail; re-running requires restoring the reverted instrumentation.


Read-only / BYTE-NEUTRAL. Runs LiftOn ONCE per benchmark on the cached
`-L`/`-M` with the `--orf-best-match` flag **OFF** (output stays byte-frozen)
but `LIFTON_ORF_DIAG=<path>` set, which makes `lifton_class.__find_orfs` emit
one TSV line per rescued transcript:

  trans_id  n_total_orfs  n_multi_orf  id_longest  id_bestK  id_bestAll
            delta_K(=id_bestK-id_longest)  delta_all  crossed_gate  cur_lifton_aa

`id_longest` is the APPLIED legacy choice (single longest ORF per frame);
`id_bestK` / `id_bestAll` are what the broadened top-K (K=3) / unbounded
best-match selection WOULD achieve (no-interior-skip scan, trailing-`*`
stripped — the shipped ON behaviour). `crossed_gate` = the broadened pick
differs from the longest AND clears the 1 % rescue gate (i.e. it would actually
change the emitted CDS).

This is the go/no-go BEFORE building the A/B (mirrors Iteration 4's
chaining_headroom.py). Signals:
  - n_gate_crossed ~0 on both datasets        -> NO-GO (nothing changes).
  - projected corpus mean-PI delta < ~+0.002  -> below the bar; NO-GO/marginal.
  - delta_all >> delta_K frequently           -> K=3 too small; revisit K.

Usage (repo root, lifton_devel env):
    python -m benchmarks.compare.orf_best_match_headroom [IDS...]   # default below
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
_EPS = 1e-6


def _ann_db(bid: str) -> str:
    for b in REG["benchmarks"]:
        if b["id"] == bid:
            return b.get("annotation_database", "RefSeq")
    return "RefSeq"


def _run_diag(bid: str, p: dict, root: Path, log=print):
    liftoff = WORK / bid / "tools" / "liftoff" / "liftoff.gff3"
    miniprot = WORK / bid / "tools" / "miniprot" / "miniprot.gff3"
    out = root / "default.gff3"
    diag = root / "orf_diag.tsv"
    if diag.exists():
        diag.unlink()
    # Flag OFF (no LIFTON_ORF_BEST_MATCH) so the emitted GFF3 is byte-frozen; the
    # diagnostic only MEASURES the broadened-selection ceiling.
    env = _compose_env(TOOLS)
    env["LIFTON_ORF_DIAG"] = str(diag)
    argv = [TOOLS["lifton_bin"], "-t", "8", "-copies", "-ad", _ann_db(bid),
            "-g", p["ref_gff"], "-L", str(liftoff), "-M", str(miniprot),
            "-o", str(out), "--native", "--locus-pipeline",
            p["tgt_fa"], p["ref_fa"]]
    pr = run_profiled(argv, label=f"orf_headroom_{bid}",
                      log_dir=root / "logs", env=env, log=log)
    if pr.exit_code != 0 or not out.exists() or out.stat().st_size == 0:
        raise RuntimeError(f"{bid}: lifton failed (exit {pr.exit_code}); "
                           f"see {pr.stderr_path}")
    return diag, root / "lifton_output" / "score.txt"


def _parse_diag(diag: Path):
    rows = []
    if not diag.exists():
        return rows
    for line in diag.read_text().splitlines():
        c = line.split("\t")
        if len(c) < 10:
            continue
        try:
            rows.append({
                "tid": c[0],
                "n_total": int(c[1]),
                "n_multi": int(c[2]),
                "id_longest": float(c[3]),
                "id_bestK": float(c[4]),
                "id_bestAll": float(c[5]),
                "delta_K": float(c[6]),
                "delta_all": float(c[7]),
                "crossed": c[8] == "1",
            })
        except ValueError:
            continue
    return rows


def _n_transcripts(score_txt: Path) -> int:
    if not score_txt.exists():
        return 0
    return sum(1 for ln in score_txt.read_text().splitlines() if ln.strip())


def _mean(xs):
    return round(sum(xs) / len(xs), 5) if xs else 0.0


def main(argv=None):
    ids = (argv if argv is not None else sys.argv[1:]) or DEFAULT_IDS
    results = []
    for bid in ids:
        print(f"=== {bid}: ORF best-match headroom diagnostic ===", flush=True)
        W = WORK / bid
        man = json.loads((W / "subset" / "subset.manifest.json").read_text())
        p = man["paths"]
        root = W / "_orf_headroom"
        root.mkdir(parents=True, exist_ok=True)
        _clean_input_dbs(p["ref_gff"],
                         W / "tools" / "liftoff" / "liftoff.gff3",
                         W / "tools" / "miniprot" / "miniprot.gff3")

        diag, score_txt = _run_diag(bid, p, root)
        rows = _parse_diag(diag)
        n_corpus = _n_transcripts(score_txt)

        n_rescued = len(rows)
        crossed = [r for r in rows if r["crossed"]]
        n_multi = sum(1 for r in rows if r["n_multi"] > 0)
        n_improved = sum(1 for r in rows if r["delta_K"] > _EPS)
        # K=3 sufficiency: how often a candidate BEYOND the top-3 would do better.
        n_K_insufficient = sum(
            1 for r in rows if r["delta_all"] > r["delta_K"] + _EPS)
        # Projected corpus mean-PI delta = sum of per-transcript gains over the
        # population that actually changes output, divided by ALL transcripts.
        sum_gain_crossed = sum(r["delta_K"] for r in crossed)
        proj_corpus_delta = (round(sum_gain_crossed / n_corpus, 5)
                             if n_corpus else 0.0)

        rec = {
            "benchmark": bid,
            "n_corpus_transcripts": n_corpus,
            "n_rescued": n_rescued,
            "n_multi_orf": n_multi,
            "n_strict_improved": n_improved,
            "n_gate_crossed": len(crossed),
            "mean_delta_K_over_rescued": _mean([r["delta_K"] for r in rows]),
            "mean_delta_over_crossed": _mean([r["delta_K"] for r in crossed]),
            "max_delta_K": round(max((r["delta_K"] for r in rows), default=0.0), 5),
            "projected_corpus_pi_delta": proj_corpus_delta,
            "n_K3_insufficient": n_K_insufficient,
        }
        results.append(rec)
        print(f"  [{bid}] corpus={n_corpus} rescued={n_rescued} "
              f"multi-ORF={n_multi}", flush=True)
        print(f"    gate-crossed (would change CDS) = {rec['n_gate_crossed']}; "
              f"mean delta over crossed = {rec['mean_delta_over_crossed']}", flush=True)
        print(f"    PROJECTED corpus mean-PI delta = "
              f"{rec['projected_corpus_pi_delta']:+.5f}  "
              f"(K=3 insufficient on {n_K_insufficient} txpts)", flush=True)
        verdict = ("GO (worth an A/B)"
                   if (rec["n_gate_crossed"] > 0
                       and rec["projected_corpus_pi_delta"] >= 0.002)
                   else "NO-GO / marginal")
        print(f"    --> {verdict}", flush=True)

    (HERE / "orf_best_match_headroom.json").write_text(json.dumps(results, indent=2))
    _write_md(results, HERE / "orf_best_match_headroom.md")
    print("\nDONE: wrote orf_best_match_headroom.json + orf_best_match_headroom.md",
          flush=True)
    return 0


def _write_md(results, path: Path):
    lines = [
        "## ORF best-match headroom diagnostic (Iteration 9, read-only)\n",
        "Default run on cached `-L`/`-M`, flag OFF, `LIFTON_ORF_DIAG` set "
        "(byte-frozen output). Per rescued transcript: applied longest-per-frame "
        "identity vs the top-K(=3)/unbounded best-match identity. "
        "**gate-crossed** = the broadened pick differs from the longest AND "
        "clears the 1 % rescue gate (would change the emitted CDS). "
        "**projected corpus mean-PI delta** = sum of gains over gate-crossed "
        "transcripts / all transcripts (the number the A/B must beat the "
        "+0.003 bar on). Decision: gate-crossed ~0 or projected < ~+0.002 on "
        "both -> NO-GO (keep opt-in, document); otherwise build the A/B.\n",
        "| Dataset | corpus | rescued | multi-ORF | gate-crossed | mean Δ over crossed | projected corpus PI Δ | K=3 insufficient |",
        "|---|---|---|---|---|---|---|---|",
    ]
    for r in results:
        lines.append("| {} | {} | {} | {} | {} | {} | {:+.5f} | {} |".format(
            r["benchmark"], r["n_corpus_transcripts"], r["n_rescued"],
            r["n_multi_orf"], r["n_gate_crossed"], r["mean_delta_over_crossed"],
            r["projected_corpus_pi_delta"], r["n_K3_insufficient"]))
    path.write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    sys.exit(main())
