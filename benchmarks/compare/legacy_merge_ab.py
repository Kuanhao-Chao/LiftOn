#!/usr/bin/env python
"""Self-contained 2-state A/B for the PROMOTED Step-7 merge.

As of the best-of-outcome promotion, the default LiftOn path runs the verified
best-of-outcome merge, and the pre-promotion *unconditional* Liftoff/miniprot
merge is reachable only via ``--legacy-merge`` (``--optimize`` is now a kept
no-op alias). This driver runs BOTH states fresh on the same cached
``-L``/``-M`` inputs and scores each with the SAME independent re-alignment
evaluator, so the comparison needs no harness eval TSVs.

  state "default"  (no flag)        best-of-outcome merge   ← the new default
  state "legacy"   --legacy-merge   unconditional merge     ← pre-promotion default

This is the inverse of the historical ``frame_gate_xspecies.py`` mapping (whose
``--optimize`` / ``LIFTON_MERGE_FRAME_GATE`` toggles are now no-ops): the new
``default`` reproduces that driver's old *state 2* (best-of-outcome) and the new
``--legacy-merge`` reproduces its old *state 1* (unconditional merge). It is the
integration-level regression guard for the promotion: ``default`` must score
>= ``legacy`` mean protein identity with a favorable improved/regressed ratio,
confirming on real divergent data what ``TestMergePromotion`` pins synthetically.

The two lifton runs are SERIAL (they reuse the same in-place gffutils DBs next
to the cached GFFs; serial execution avoids SQLite read/write contention).
Input DBs are cleaned ONCE up front so the first run rebuilds them cleanly and
the second reuses them.

Usage (repo root, lifton_devel env; miniprot must already be transcript-space):
    python -m benchmarks.compare.legacy_merge_ab [IDS...]   # default = drosophila mouse_to_rat
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

DEFAULT_IDS = ["drosophila", "mouse_to_rat"]
STATES = ("default", "legacy")


def _ann_db(bid: str) -> str:
    for b in REG["benchmarks"]:
        if b["id"] == bid:
            return b.get("annotation_database", "RefSeq")
    return "RefSeq"


def _state_argv(bid: str, p: dict, state: str, out: Path, liftoff_gff: Path,
                miniprot_gff: Path):
    """Build (argv, env) for one of the two merge states."""
    argv = [TOOLS["lifton_bin"], "-t", "8", "-copies", "-ad", _ann_db(bid),
            "-g", p["ref_gff"], "-L", str(liftoff_gff), "-M", str(miniprot_gff),
            "-o", str(out), "--native", "--locus-pipeline"]
    env = _compose_env(TOOLS)
    if state == "default":
        pass                          # best-of-outcome (the promoted default)
    elif state == "legacy":
        argv.append("--legacy-merge")  # pre-promotion unconditional merge
    else:
        raise ValueError(state)
    argv += [p["tgt_fa"], p["ref_fa"]]
    return argv, env


def _run_state(bid: str, p: dict, state: str, ab_root: Path, log=print):
    liftoff_gff = WORK / bid / "tools" / "liftoff" / "liftoff.gff3"
    miniprot_gff = WORK / bid / "tools" / "miniprot" / "miniprot.gff3"
    statedir = ab_root / state
    statedir.mkdir(parents=True, exist_ok=True)
    out = statedir / f"{state}.gff3"
    argv, env = _state_argv(bid, p, state, out, liftoff_gff, miniprot_gff)
    pr = run_profiled(argv, label=f"legacy_ab_{state}",
                      log_dir=ab_root / "logs", env=env, log=log)
    if pr.exit_code != 0 or not out.exists() or out.stat().st_size == 0:
        raise RuntimeError(
            f"{bid}: state {state} lifton failed (exit {pr.exit_code}); "
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
        print(f"=== {bid} (default best-of-outcome vs --legacy-merge) ===", flush=True)
        W = WORK / bid
        man = json.loads((W / "subset" / "subset.manifest.json").read_text())
        p = man["paths"]
        ab_root = W / "_legacy_ab"
        ab_root.mkdir(parents=True, exist_ok=True)

        liftoff_gff = W / "tools" / "liftoff" / "liftoff.gff3"
        miniprot_gff = W / "tools" / "miniprot" / "miniprot.gff3"
        # Clean ONCE so the first run rebuilds the input DBs cleanly; reuse after.
        _clean_input_dbs(p["ref_gff"], liftoff_gff, miniprot_gff)

        gffs, prs, score_txts = {}, {}, {}
        for state in STATES:
            print(f"--- {bid} state {state} ---", flush=True)
            g, pr, st = _run_state(bid, p, state, ab_root, log=print)
            gffs[state], prs[state], score_txts[state] = g, pr, st

        # Score both with the SAME evaluator.
        ref, ref_index = evaluator.build_reference(p["ref_gff"], p["ref_fa"], log=print)
        eval_dir = ab_root / "eval"
        eval_dir.mkdir(parents=True, exist_ok=True)
        for state in STATES:
            evaluator.evaluate_tool(state, str(gffs[state]), p["tgt_fa"],
                                    ref, man, eval_dir, None, log=print,
                                    ref_index=ref_index, threads=8)

        pi = {s: _pi_by_ref(eval_dir / f"{s}.transcripts.tsv") for s in STATES}
        fb = {s: _score_annotation_counts(score_txts[s]) for s in STATES}

        rec = {
            "benchmark": bid,
            "n_scored": {s: len(pi[s]) for s in STATES},
            "n_complete": {s: _n_complete(pi[s]) for s in STATES},
            "mean_pi": {"default_bestof": _mean(pi["default"].values()),
                        "legacy_uncond": _mean(pi["legacy"].values())},
            # default vs legacy: the promotion win (best-of-outcome vs unconditional)
            "delta_default_vs_legacy": _delta_counts(pi["legacy"], pi["default"]),
            "merge_kept": {s: fb[s]["LiftOn_chaining_algorithm"] for s in STATES},
            "fallback_to_liftoff": {s: fb[s]["Liftoff"] for s in STATES},
            "wall_s": {s: round(prs[s].wall_clock_seconds, 2) for s in STATES},
        }
        results.append(rec)
        m, d = rec["mean_pi"], rec["delta_default_vs_legacy"]
        print(f"  [{bid}] mean PI: default={m['default_bestof']} legacy={m['legacy_uncond']} "
              f"(default vs legacy: {d['improved']} improved / {d['regressed']} regressed, "
              f"net {d['net_per_transcript']:+.6f}); merge-kept "
              f"default={rec['merge_kept']['default']} legacy={rec['merge_kept']['legacy']}",
              flush=True)

    (HERE / "legacy_merge_ab.json").write_text(json.dumps(results, indent=2))
    _write_md(results, HERE / "legacy_merge_ab.md")
    print("\nDONE: wrote legacy_merge_ab.json + legacy_merge_ab.md", flush=True)
    return 0


def _write_md(results, path: Path):
    lines = [
        "## Promoted merge: default best-of-outcome vs `--legacy-merge` A/B (Step 7)\n",
        "Self-contained — both states run fresh on the same cached `-L`/`-M` "
        "(genes mode, transcript-space miniprot), scored by the independent "
        "re-alignment evaluator. **default** = the promoted best-of-outcome merge "
        "(no flag); **legacy** = the pre-promotion unconditional merge "
        "(`--legacy-merge`). Acceptance: default >= legacy mean protein identity "
        "with a favorable improved/regressed ratio; completeness not reduced.\n",
        "| Dataset | default mean PI | legacy mean PI | default vs legacy (impr/regr, net) "
        "| merge-kept default/legacy | complete default/legacy | wall default/legacy |",
        "|---|---|---|---|---|---|---|",
    ]
    for r in results:
        m, d = r["mean_pi"], r["delta_default_vs_legacy"]
        mk, nc, w = r["merge_kept"], r["n_complete"], r["wall_s"]
        lines.append(
            "| {} | {} | {} | {}/{}, {:+.6f} | {}/{} | {}/{} | {}/{} |".format(
                r["benchmark"], m["default_bestof"], m["legacy_uncond"],
                d["improved"], d["regressed"], d["net_per_transcript"],
                mk["default"], mk["legacy"],
                nc["default"], nc["legacy"],
                w["default"], w["legacy"]))
    path.write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    sys.exit(main())
