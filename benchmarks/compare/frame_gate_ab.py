#!/usr/bin/env python
"""3-state A/B for Iteration-1's frame-consistency gate (Step 7 merge accuracy).

For a divergent transcript-space dataset (default: drosophila), attribute the
frame gate's marginal effect by comparing three merge states, all on the same
cached -L/-M inputs (genes mode):

  state 1  default                                  unconditional merge (byte-frozen)
  state 2  --optimize  LIFTON_MERGE_FRAME_GATE=0     best-of-outcome alone
  state 3  --optimize  LIFTON_MERGE_FRAME_GATE=1     best-of-outcome + frame gate

States 1 and 3 are produced by the main 4-variant harness (`lifton` and
`lifton_optimize` tools — run `run_compare ... --feature-modes genes` first);
this script runs only the missing state 2 and reads 1/3 from `work/<id>/eval/`.
It scores every state with the SAME independent re-alignment evaluator
(`evaluator.evaluate_tool`) and counts best-of-outcome fallbacks from each run's
`lifton_output/score.txt` (column 6 = annotation: "Liftoff" = fell back,
"LiftOn_chaining_algorithm" = merge kept).

Acceptance: state3 mean protein identity >= state2 >= state1; fewer fallbacks
(more "LiftOn_chaining_algorithm" kept) under state 3 than state 2.

Usage (repo root, lifton_devel env, AFTER the harness genes run):
    python -m benchmarks.compare.frame_gate_ab [IDS...]   # default = drosophila
"""
from __future__ import annotations

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

DEFAULT_IDS = ["drosophila"]


def _ann_db(bid: str) -> str:
    for b in REG["benchmarks"]:
        if b["id"] == bid:
            return b.get("annotation_database", "RefSeq")
    return "RefSeq"


def _run_state2(bid: str, p: dict):
    """Run state 2: --optimize with the frame gate explicitly OFF."""
    W = WORK / bid
    fg = W / "_framegate"
    log_dir = fg / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    liftoff_gff = W / "tools" / "liftoff" / "liftoff.gff3"
    miniprot_gff = W / "tools" / "miniprot" / "miniprot.gff3"
    out = fg / "state2_optimize_nogate.gff3"
    _clean_input_dbs(p["ref_gff"], liftoff_gff, miniprot_gff)
    argv = [TOOLS["lifton_bin"], "-t", "8", "-copies", "-ad", _ann_db(bid),
            "-g", p["ref_gff"], "-L", str(liftoff_gff), "-M", str(miniprot_gff),
            "-o", str(out), "--native", "--locus-pipeline", "--optimize",
            p["tgt_fa"], p["ref_fa"]]
    env = _compose_env(TOOLS)
    env["LIFTON_MERGE_FRAME_GATE"] = "0"   # state 2: best-of-outcome alone
    pr = run_profiled(argv, label="framegate_state2", log_dir=log_dir, env=env, log=print)
    if pr.exit_code != 0 or not out.exists() or out.stat().st_size == 0:
        raise RuntimeError(f"{bid}: state2 lifton failed (exit {pr.exit_code}); see {pr.stderr_path}")
    return out, pr


def _score_annotation_counts(score_txt: Path) -> dict:
    """Count per-transcript merge provenance from score.txt col 6 (annotation)."""
    counts = {"LiftOn_chaining_algorithm": 0, "Liftoff": 0, "other": 0, "total": 0}
    if not score_txt.exists():
        return counts
    for line in score_txt.read_text().splitlines():
        cols = line.split("\t")
        if len(cols) < 6:
            continue
        ann = cols[5]
        counts["total"] += 1
        counts[ann if ann in ("LiftOn_chaining_algorithm", "Liftoff") else "other"] += 1
    return counts


def _pi_by_ref(tsv: Path) -> dict:
    """ref_mrna_id -> protein_identity for coding, scored rows."""
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


def _mean(vals):
    vals = [v for v in vals if v is not None]
    return round(sum(vals) / len(vals), 5) if vals else None


def _delta_counts(base: dict, new: dict):
    """Per-transcript improved/regressed of `new` vs `base` on common ref ids."""
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
        print(f"=== {bid} ===", flush=True)
        W = WORK / bid
        man = json.loads((W / "subset" / "subset.manifest.json").read_text())
        p = man["paths"]

        # state 1 / 3 from the main harness (must be fresh: re-run the genes
        # harness AFTER the frame-gate code change before running this).
        s1_tsv = W / "eval" / "lifton.transcripts.tsv"
        s3_tsv = W / "eval" / "lifton_optimize.transcripts.tsv"
        for t in (s1_tsv, s3_tsv):
            if not t.exists():
                raise RuntimeError(f"{bid}: missing harness eval {t}. Run "
                                   f"run_compare --feature-modes genes first.")

        # state 2: run + evaluate with the same scorer the harness uses
        s2_gff, s2_pr = _run_state2(bid, p)
        ref, ref_index = evaluator.build_reference(p["ref_gff"], p["ref_fa"], log=print)
        eval_dir = W / "_framegate" / "eval"
        eval_dir.mkdir(parents=True, exist_ok=True)
        evaluator.evaluate_tool("state2_nogate", str(s2_gff), p["tgt_fa"], ref, man,
                                eval_dir, None, log=print, ref_index=ref_index, threads=8)
        s2_tsv = eval_dir / "state2_nogate.transcripts.tsv"

        pi1, pi2, pi3 = _pi_by_ref(s1_tsv), _pi_by_ref(s2_tsv), _pi_by_ref(s3_tsv)
        fb1 = _score_annotation_counts(W / "tools" / "lifton" / "lifton_output" / "score.txt")
        fb2 = _score_annotation_counts(W / "_framegate" / "lifton_output" / "score.txt")
        fb3 = _score_annotation_counts(W / "tools" / "lifton_optimize" / "lifton_output" / "score.txt")

        rec = {
            "benchmark": bid,
            "mean_pi": {"state1_default": _mean(pi1.values()),
                        "state2_bestof": _mean(pi2.values()),
                        "state3_bestof_framegate": _mean(pi3.values())},
            "delta_vs_state1": {"state2": _delta_counts(pi1, pi2),
                                "state3": _delta_counts(pi1, pi3)},
            "delta_state3_vs_state2": _delta_counts(pi2, pi3),
            "merge_kept": {"state1": fb1["LiftOn_chaining_algorithm"],
                           "state2": fb2["LiftOn_chaining_algorithm"],
                           "state3": fb3["LiftOn_chaining_algorithm"]},
            "fallback_to_liftoff": {"state1": fb1["Liftoff"], "state2": fb2["Liftoff"],
                                    "state3": fb3["Liftoff"]},
            "state2_wall_s": round(s2_pr.wall_clock_seconds, 2),
        }
        results.append(rec)
        m = rec["mean_pi"]
        print(f"  [{bid}] mean PI: s1={m['state1_default']} s2={m['state2_bestof']} "
              f"s3={m['state3_bestof_framegate']}  merge_kept s2={rec['merge_kept']['state2']} "
              f"s3={rec['merge_kept']['state3']}", flush=True)

    (HERE / "frame_gate_ab.json").write_text(json.dumps(results, indent=2))
    _write_md(results, HERE / "frame_gate_ab.md")
    print(f"\nDONE: wrote frame_gate_ab.json + frame_gate_ab.md", flush=True)
    return 0


def _write_md(results, path: Path):
    lines = [
        "## Frame-consistency gate: 3-state A/B (Step 7 merge)\n",
        "All three states on the same cached `-L`/`-M` (genes mode), scored by the "
        "independent re-alignment evaluator. **state 1** = default (unconditional "
        "merge); **state 2** = `--optimize` best-of-outcome alone "
        "(`LIFTON_MERGE_FRAME_GATE=0`); **state 3** = best-of-outcome + frame gate "
        "(`=1`). Acceptance: state3 ≥ state2 ≥ state1 mean protein identity; more "
        "merges kept (fewer fallbacks) under state 3.\n",
        "| Dataset | s1 mean PI | s2 mean PI | s3 mean PI | s3 vs s2 (impr/regr, net) "
        "| merge-kept s2→s3 | fallback s2→s3 |",
        "|---|---|---|---|---|---|---|",
    ]
    for r in results:
        m, d32 = r["mean_pi"], r["delta_state3_vs_state2"]
        mk, fb = r["merge_kept"], r["fallback_to_liftoff"]
        lines.append("| {} | {} | {} | {} | {}/{}, {:+.5f} | {}→{} | {}→{} |".format(
            r["benchmark"], m["state1_default"], m["state2_bestof"],
            m["state3_bestof_framegate"], d32["improved"], d32["regressed"],
            d32["net_per_transcript"], mk["state2"], mk["state3"],
            fb["state2"], fb["state3"]))
    path.write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    sys.exit(main())
