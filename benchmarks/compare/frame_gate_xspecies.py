#!/usr/bin/env python
"""Self-contained 3-state A/B for the Step-7 merge on a *mammalian cross-species*
dataset (default: mouse_to_rat), where the harness artifacts are absent or stale.

Unlike ``frame_gate_ab.py`` (which reads states 1 & 3 from a prior 4-variant
harness run), this driver runs ALL THREE states fresh on the same cached
``-L``/``-M`` inputs and scores every one with the SAME independent re-alignment
evaluator (``evaluator.evaluate_tool``). That removes any dependence on
possibly-stale ``eval/*.transcripts.tsv`` and makes the experiment fully
reproducible from the cached Liftoff/miniprot GFFs alone.

  state 1  default                                  unconditional merge (byte-frozen)
  state 2  --optimize  LIFTON_MERGE_FRAME_GATE=0     best-of-outcome alone
  state 3  --optimize  LIFTON_MERGE_FRAME_GATE=1     best-of-outcome + frame gate

The three lifton runs are SERIAL (they reuse the same in-place gffutils DBs next
to the cached GFFs; serial execution avoids any SQLite read/write contention).
Input DBs are cleaned ONCE up front so state 1 rebuilds them cleanly and states
2/3 reuse them.

This is the mammalian counterpart to the drosophila A/B documented in
``lifton-frame-gate-inert``: it answers (a) does best-of-outcome (state 2 > 1)
pay off on mammalian cross-species, and (b) does the frame gate add anything on
top (state 3 > 2)?

Usage (repo root, lifton_devel env; miniprot must already be transcript-space):
    python -m benchmarks.compare.frame_gate_xspecies [IDS...]   # default = mouse_to_rat
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

DEFAULT_IDS = ["mouse_to_rat"]


def _ann_db(bid: str) -> str:
    for b in REG["benchmarks"]:
        if b["id"] == bid:
            return b.get("annotation_database", "RefSeq")
    return "RefSeq"


def _state_argv(bid: str, p: dict, state: int, out: Path, liftoff_gff: Path,
                miniprot_gff: Path):
    """Build (argv, env) for one of the three merge states."""
    argv = [TOOLS["lifton_bin"], "-t", "8", "-copies", "-ad", _ann_db(bid),
            "-g", p["ref_gff"], "-L", str(liftoff_gff), "-M", str(miniprot_gff),
            "-o", str(out), "--native", "--locus-pipeline"]
    env = _compose_env(TOOLS)
    if state == 1:
        pass                                   # default path, byte-frozen
    elif state == 2:
        argv.append("--optimize")
        env["LIFTON_MERGE_FRAME_GATE"] = "0"   # best-of-outcome alone
    elif state == 3:
        argv.append("--optimize")
        env["LIFTON_MERGE_FRAME_GATE"] = "1"   # best-of-outcome + frame gate
    else:
        raise ValueError(state)
    argv += [p["tgt_fa"], p["ref_fa"]]
    return argv, env


def _run_state(bid: str, p: dict, state: int, fg_root: Path, log=print):
    liftoff_gff = WORK / bid / "tools" / "liftoff" / "liftoff.gff3"
    miniprot_gff = WORK / bid / "tools" / "miniprot" / "miniprot.gff3"
    statedir = fg_root / f"state{state}"
    statedir.mkdir(parents=True, exist_ok=True)
    out = statedir / f"state{state}.gff3"
    argv, env = _state_argv(bid, p, state, out, liftoff_gff, miniprot_gff)
    pr = run_profiled(argv, label=f"framegate_xs_state{state}",
                      log_dir=fg_root / "logs", env=env, log=log)
    if pr.exit_code != 0 or not out.exists() or out.stat().st_size == 0:
        raise RuntimeError(
            f"{bid}: state{state} lifton failed (exit {pr.exit_code}); "
            f"see {pr.stderr_path}")
    score_txt = statedir / "lifton_output" / "score.txt"
    return out, pr, score_txt


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


def _n_complete(pi: dict) -> int:
    return sum(1 for v in pi.values() if v >= 1.0 - 1e-9)


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
        print(f"=== {bid} (self-contained 3-state) ===", flush=True)
        W = WORK / bid
        man = json.loads((W / "subset" / "subset.manifest.json").read_text())
        p = man["paths"]
        fg_root = W / "_framegate_xs"
        fg_root.mkdir(parents=True, exist_ok=True)

        liftoff_gff = W / "tools" / "liftoff" / "liftoff.gff3"
        miniprot_gff = W / "tools" / "miniprot" / "miniprot.gff3"
        # Clean ONCE so state 1 rebuilds the input DBs cleanly; states 2/3 reuse.
        _clean_input_dbs(p["ref_gff"], liftoff_gff, miniprot_gff)

        gffs, prs, score_txts = {}, {}, {}
        for state in (1, 2, 3):
            print(f"--- {bid} state {state} ---", flush=True)
            g, pr, st = _run_state(bid, p, state, fg_root, log=print)
            gffs[state], prs[state], score_txts[state] = g, pr, st

        # Score all three with the SAME evaluator.
        ref, ref_index = evaluator.build_reference(p["ref_gff"], p["ref_fa"], log=print)
        eval_dir = fg_root / "eval"
        eval_dir.mkdir(parents=True, exist_ok=True)
        for state in (1, 2, 3):
            evaluator.evaluate_tool(f"state{state}", str(gffs[state]), p["tgt_fa"],
                                    ref, man, eval_dir, None, log=print,
                                    ref_index=ref_index, threads=8)

        pi = {s: _pi_by_ref(eval_dir / f"state{s}.transcripts.tsv") for s in (1, 2, 3)}
        fb = {s: _score_annotation_counts(score_txts[s]) for s in (1, 2, 3)}

        rec = {
            "benchmark": bid,
            "n_scored": {f"state{s}": len(pi[s]) for s in (1, 2, 3)},
            "n_complete": {f"state{s}": _n_complete(pi[s]) for s in (1, 2, 3)},
            "mean_pi": {"state1_default": _mean(pi[1].values()),
                        "state2_bestof": _mean(pi[2].values()),
                        "state3_bestof_framegate": _mean(pi[3].values())},
            "delta_vs_state1": {"state2": _delta_counts(pi[1], pi[2]),
                                "state3": _delta_counts(pi[1], pi[3])},
            "delta_state3_vs_state2": _delta_counts(pi[2], pi[3]),
            "merge_kept": {f"state{s}": fb[s]["LiftOn_chaining_algorithm"] for s in (1, 2, 3)},
            "fallback_to_liftoff": {f"state{s}": fb[s]["Liftoff"] for s in (1, 2, 3)},
            "wall_s": {f"state{s}": round(prs[s].wall_clock_seconds, 2) for s in (1, 2, 3)},
        }
        results.append(rec)
        m = rec["mean_pi"]
        print(f"  [{bid}] mean PI: s1={m['state1_default']} s2={m['state2_bestof']} "
              f"s3={m['state3_bestof_framegate']}  merge_kept "
              f"s1={rec['merge_kept']['state1']} s2={rec['merge_kept']['state2']} "
              f"s3={rec['merge_kept']['state3']}", flush=True)

    (HERE / "frame_gate_xspecies.json").write_text(json.dumps(results, indent=2))
    _write_md(results, HERE / "frame_gate_xspecies.md")
    print("\nDONE: wrote frame_gate_xspecies.json + frame_gate_xspecies.md", flush=True)
    return 0


def _write_md(results, path: Path):
    lines = [
        "## Frame-consistency gate: mammalian cross-species 3-state A/B (Step 7 merge)\n",
        "Self-contained — all three states run fresh on the same cached `-L`/`-M` "
        "(genes mode, transcript-space miniprot), scored by the independent "
        "re-alignment evaluator. **state 1** = default (unconditional merge); "
        "**state 2** = `--optimize` best-of-outcome alone (`LIFTON_MERGE_FRAME_GATE=0`); "
        "**state 3** = best-of-outcome + frame gate (`=1`). Acceptance: state3 ≥ "
        "state2 ≥ state1 mean protein identity; completeness not reduced.\n",
        "| Dataset | s1 mean PI | s2 mean PI | s3 mean PI | s2 vs s1 (impr/regr, net) "
        "| s3 vs s2 (impr/regr, net) | merge-kept s1→s2→s3 | complete s1→s2→s3 | wall s1/s2/s3 |",
        "|---|---|---|---|---|---|---|---|---|",
    ]
    for r in results:
        m = r["mean_pi"]
        d21, d32 = r["delta_vs_state1"]["state2"], r["delta_state3_vs_state2"]
        mk, nc, w = r["merge_kept"], r["n_complete"], r["wall_s"]
        lines.append(
            "| {} | {} | {} | {} | {}/{}, {:+.5f} | {}/{}, {:+.5f} | {}→{}→{} | {}→{}→{} | {}/{}/{} |".format(
                r["benchmark"], m["state1_default"], m["state2_bestof"],
                m["state3_bestof_framegate"],
                d21["improved"], d21["regressed"], d21["net_per_transcript"],
                d32["improved"], d32["regressed"], d32["net_per_transcript"],
                mk["state1"], mk["state2"], mk["state3"],
                nc["state1"], nc["state2"], nc["state3"],
                w["state1"], w["state2"], w["state3"]))
    path.write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    sys.exit(main())
