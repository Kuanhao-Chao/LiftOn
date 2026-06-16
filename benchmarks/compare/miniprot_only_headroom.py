#!/usr/bin/env python
"""Iteration-15 feasibility diagnostic — miniprot-only merge candidate headroom (Gap #1).

AUDIT RECORD (Iteration-15 NO-GO). The byte-neutral `LIFTON_MINIPROT_ONLY_DIAG`
instrumentation this driver relies on (in `run_liftoff.process_liftoff_with_protein`)
was implemented, measured, and **REVERTED** — adding a third "miniprot-only"
candidate to the 2-way best-of-outcome merge is negligible: projected corpus
mean-PI delta **+0.00005 (drosophila) / +0.00041 (mouse_to_rat)**, both far below
the +0.003 bar (only 26/6699 and 20/1983 merge-fired transcripts would change).
The 2-way merge already captures the miniprot signal at the chaining sync points.
This script + its committed `miniprot_only_headroom.{json,md}` are the audit
trail; re-running requires restoring the reverted instrumentation. See the
CLAUDE.md Iteration-15 note + memory.

Read-only / BYTE-NEUTRAL. The best-of-outcome merge
(`run_liftoff.process_liftoff_with_protein`) picks the better EMITTED protein of
only {chained-merge+ORF, Liftoff+ORF}; it never considers a **miniprot-only**
candidate. This driver runs LiftOn ONCE per dataset on the cached `-L`/`-M` with
the feature OFF (output byte-frozen) but `LIFTON_MINIPROT_ONLY_DIAG=<path>` set,
which makes `process_liftoff_with_protein` emit one TSV line per merge-fired
transcript:

  trans_id  ref_trans_id  applied_annotation  applied_aa  merge_aa  liftoff_aa
            mp_only_aa  delta(=mp_only_aa-applied_aa)  would_change  same_locus  locus

`applied_aa` is the APPLIED best-of-two identity; `mp_only_aa` is what a
miniprot-only CDS (all miniprot blocks, faithful strand handling) WOULD score
after ORF-rescue, trial-applied then reverted (the function returns byte-exact).
`would_change` = mp_only_aa strictly beats the applied result. `same_locus` =
the miniprot CDS span overlaps the Liftoff transcript span — only same-locus
candidates are coordinate-compatible with an in-place swap (the Gap #1 path);
the rest belong to the Step-8 separate-gene path (Gap #2,
`miniprot_only_rescue_headroom.py`).

GO/NO-GO (the Iter-4/9/13 feasibility-first recipe), over same_locus=1 rows only:
  - n_would_change ~0 on both datasets        -> NO-GO (the merge already
                                                  captures miniprot at sync points).
  - projected corpus mean-PI delta >= +0.003  -> GO to the A/B (the orf_best_match
                                                  bar); < +0.002 on both -> NO-GO.
The projection is the UPPER BOUND (it ignores downstream ripple); the A/B is
mandatory on a GO.

Usage (repo root, lifton_devel env):
    python -m benchmarks.compare.miniprot_only_headroom [IDS...]   # default below
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
_EPS = 1e-9
GO_BAR = 0.003


def _ann_db(bid: str) -> str:
    for b in REG["benchmarks"]:
        if b["id"] == bid:
            return b.get("annotation_database", "RefSeq")
    return "RefSeq"


def _run_diag(bid: str, p: dict, root: Path, log=print):
    liftoff = WORK / bid / "tools" / "liftoff" / "liftoff.gff3"
    miniprot = WORK / bid / "tools" / "miniprot" / "miniprot.gff3"
    out = root / "default.gff3"
    diag = root / "mp_only_diag.tsv"
    if diag.exists():
        diag.unlink()
    # Feature OFF -> emitted GFF3 byte-frozen; the diagnostic only MEASURES.
    env = _compose_env(TOOLS)
    env["LIFTON_MINIPROT_ONLY_DIAG"] = str(diag)
    argv = [TOOLS["lifton_bin"], "-t", "8", "-copies", "-ad", _ann_db(bid),
            "-g", p["ref_gff"], "-L", str(liftoff), "-M", str(miniprot),
            "-o", str(out), "--native", "--locus-pipeline",
            p["tgt_fa"], p["ref_fa"]]
    pr = run_profiled(argv, label=f"mp_only_headroom_{bid}",
                      log_dir=root / "logs", env=env, log=log)
    if pr.exit_code != 0 or not out.exists() or out.stat().st_size == 0:
        raise RuntimeError(f"{bid}: lifton failed (exit {pr.exit_code}); "
                           f"see {pr.stderr_path}")
    return diag, root / "lifton_output" / "score.txt"


def _parse_diag(diag: Path):
    """Dedupe by trans_id (defensive); keep the first occurrence."""
    rows = {}
    if not diag.exists():
        return []
    for line in diag.read_text().splitlines():
        c = line.split("\t")
        if len(c) < 11:
            continue
        tid = c[0]
        if tid in rows:
            continue

        def _f(x):
            try:
                return float(x)
            except (ValueError, TypeError):
                return None
        rows[tid] = {
            "tid": tid,
            "applied_ann": c[2],
            "applied_aa": _f(c[3]),
            "merge_aa": _f(c[4]),
            "liftoff_aa": _f(c[5]),
            "mp_only_aa": _f(c[6]),
            "delta": _f(c[7]),
            "would_change": c[8] == "1",
            "same_locus": c[9] == "1",
        }
    return list(rows.values())


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
        print(f"=== {bid}: miniprot-only merge-candidate headroom (Gap #1) ===", flush=True)
        W = WORK / bid
        man = json.loads((W / "subset" / "subset.manifest.json").read_text())
        p = man["paths"]
        root = W / "_mp_only_headroom"
        root.mkdir(parents=True, exist_ok=True)
        _clean_input_dbs(p["ref_gff"],
                         W / "tools" / "liftoff" / "liftoff.gff3",
                         W / "tools" / "miniprot" / "miniprot.gff3")

        diag, score_txt = _run_diag(bid, p, root)
        rows = _parse_diag(diag)
        n_corpus = _n_transcripts(score_txt)

        n_merge_fired = len(rows)
        # Gap #1 is the in-place swap; only same-locus candidates are valid here.
        sl = [r for r in rows if r["same_locus"] and r["mp_only_aa"] is not None]
        changed = [r for r in sl if r["would_change"] and (r["delta"] or 0) > _EPS]
        sum_gain = sum((r["delta"] or 0) for r in changed)
        proj_corpus_delta = round(sum_gain / n_corpus, 5) if n_corpus else 0.0
        n_diff_locus_change = sum(
            1 for r in rows if (not r["same_locus"]) and r["would_change"])

        rec = {
            "benchmark": bid,
            "n_corpus_transcripts": n_corpus,
            "n_merge_fired": n_merge_fired,
            "n_same_locus": len(sl),
            "n_would_change_same_locus": len(changed),
            "mean_delta_over_changed": _mean([(r["delta"] or 0) for r in changed]),
            "max_delta_same_locus": round(
                max(((r["delta"] or 0) for r in sl), default=0.0), 5),
            "projected_corpus_pi_delta": proj_corpus_delta,
            "n_diff_locus_would_change": n_diff_locus_change,
        }
        results.append(rec)
        print(f"  [{bid}] corpus={n_corpus} merge-fired={n_merge_fired} "
              f"same-locus={len(sl)}", flush=True)
        print(f"    would-change (same-locus, mp-only beats applied) = "
              f"{len(changed)}; mean Δ over changed = "
              f"{rec['mean_delta_over_changed']}", flush=True)
        print(f"    PROJECTED corpus mean-PI delta = "
              f"{proj_corpus_delta:+.5f}  (diff-locus would-change = "
              f"{n_diff_locus_change}, deferred to Gap #2)", flush=True)
        verdict = ("GO (worth an A/B)"
                   if (len(changed) > 0 and proj_corpus_delta >= GO_BAR)
                   else "NO-GO / marginal")
        print(f"    --> {verdict}", flush=True)

    (HERE / "miniprot_only_headroom.json").write_text(json.dumps(results, indent=2))
    _write_md(results, HERE / "miniprot_only_headroom.md")
    any_go = any(r["n_would_change_same_locus"] > 0
                 and r["projected_corpus_pi_delta"] >= GO_BAR for r in results)
    print(f"\nOVERALL: {'GO' if any_go else 'NO-GO'} "
          f"(GO if projected corpus PI Δ >= +{GO_BAR} on any dataset).", flush=True)
    print("DONE: wrote miniprot_only_headroom.json + miniprot_only_headroom.md",
          flush=True)
    return 0


def _write_md(results, path: Path):
    lines = [
        "## miniprot-only merge-candidate headroom (Iteration 15, Gap #1, read-only)\n",
        "Default run on cached `-L`/`-M`, feature OFF, `LIFTON_MINIPROT_ONLY_DIAG` "
        "set (byte-frozen output). Per merge-fired transcript: the APPLIED "
        "best-of-{merge, Liftoff} identity vs what a **miniprot-only** CDS would "
        "score after ORF-rescue. **would-change (same-locus)** = miniprot-only "
        "strictly beats the applied result AND its CDS span overlaps the Liftoff "
        "transcript (coordinate-compatible with an in-place swap). **projected "
        "corpus mean-PI delta** = sum of same-locus gains / all transcripts (the "
        f"number the A/B must beat the +{GO_BAR} bar on). Decision: would-change "
        f"~0 or projected < +0.002 on both -> NO-GO (document); else build the A/B.\n",
        "| Dataset | corpus | merge-fired | same-locus | would-change | mean Δ over changed | projected corpus PI Δ | diff-locus would-change |",
        "|---|---|---|---|---|---|---|---|",
    ]
    for r in results:
        lines.append("| {} | {} | {} | {} | {} | {} | {:+.5f} | {} |".format(
            r["benchmark"], r["n_corpus_transcripts"], r["n_merge_fired"],
            r["n_same_locus"], r["n_would_change_same_locus"],
            r["mean_delta_over_changed"], r["projected_corpus_pi_delta"],
            r["n_diff_locus_would_change"]))
    path.write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    sys.exit(main())
