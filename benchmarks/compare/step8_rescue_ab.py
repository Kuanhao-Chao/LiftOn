#!/usr/bin/env python
"""Iteration-13 A/B — Step-8 rescue: default vs RELAXED thresholds.

AUDIT RECORD (Iteration-13 NO-GO). This A/B uses only the existing
`-overlap`/`-min_miniprot`/`-max_miniprot` CLI args (no reverted instrumentation
needed), so it remains runnable. Its committed result is the decisive NO-GO:
relaxing the Step-8 thresholds adds only **2 transcripts on drosophila at mean
protein identity 0.35** (both <0.7 — garbage) and **0 on mouse_to_rat**; the
overall mean PI even ticks down. The Step-8 overlap/length filters are well-tuned
(they correctly remove redundant duplicates + truncated fragments). See the
CLAUDE.md Iteration-13 note + memory.

Confirms (on emitted output) what the headroom diagnostic
(`step8_rescue_headroom.py`) projects: does relaxing the Step-8 miniprot-rescue
filters capture more genuinely-missing transcripts, or just admit garbage?

The thresholds are already CLI args, so NO `lifton/` code change — the `relaxed`
state just passes looser `-overlap`/`-min_miniprot`/`-max_miniprot` (values chosen
from the diagnostic's near-miss buckets: overlap near-misses in (0.1,0.2], length
in [0.8,0.9)). Both states reuse cached `-L`/`-M`.

  state "default"  (no threshold flags)              0.1 / 0.9 / 1.5
  state "relaxed"  -overlap 0.2 -min_miniprot 0.8 -max_miniprot 1.6

The CRITICAL false-rescue measure: how many transcripts does `relaxed` ADD, and
what is the **mean protein identity of the ADDED transcripts** (legit rescues are
high-identity; garbage is low)? Plus validity (no new gff3_validator errors) and
no perturbation of the COMMON set (a relaxed -overlap could shadow a Liftoff
transcript). Decision gate: relaxing is worth promoting as the default iff
n_added>0 AND mean_pi_added >= 0.75 AND regressed==0 AND no new validity errors,
on BOTH datasets.

Usage (repo root, lifton_devel env):
    python -m benchmarks.compare.step8_rescue_ab [IDS...]
"""
from __future__ import annotations

import csv
import json
import subprocess
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
STATES = {
    "default": [],
    "relaxed": ["-overlap", "0.2", "-min_miniprot", "0.8", "-max_miniprot", "1.6"],
}
PI_BAR = 0.75


def _ann_db(bid: str) -> str:
    for b in REG["benchmarks"]:
        if b["id"] == bid:
            return b.get("annotation_database", "RefSeq")
    return "RefSeq"


def _run_state(bid, state, flags, p, root):
    liftoff = WORK / bid / "tools" / "liftoff" / "liftoff.gff3"
    miniprot = WORK / bid / "tools" / "miniprot" / "miniprot.gff3"
    statedir = root / state
    statedir.mkdir(parents=True, exist_ok=True)
    out = statedir / f"{state}.gff3"
    _clean_input_dbs(p["ref_gff"], liftoff, miniprot)
    argv = [TOOLS["lifton_bin"], "-t", "8", "-copies", "-ad", _ann_db(bid),
            "-g", p["ref_gff"], "-L", str(liftoff), "-M", str(miniprot),
            "-o", str(out), "--native", "--locus-pipeline", *flags,
            p["tgt_fa"], p["ref_fa"]]
    pr = run_profiled(argv, label=f"step8_ab_{state}",
                      log_dir=root / "logs", env=_compose_env(TOOLS), log=print)
    if pr.exit_code != 0 or not out.exists() or out.stat().st_size == 0:
        raise RuntimeError(f"{bid}: state {state} failed (exit {pr.exit_code}); "
                           f"see {pr.stderr_path}")
    return out, pr


def _pi_by_ref(tsv: Path) -> dict:
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


def _validate(gff: Path):
    r = subprocess.run([TOOLS["lifton_python"], "-m", "lifton.gff3_validator", str(gff)],
                       env=_compose_env(TOOLS), capture_output=True, text=True)
    out = (r.stdout or "") + (r.stderr or "")
    return r.returncode, sum(1 for ln in out.splitlines() if "error" in ln.lower())


def _mean(xs):
    xs = [x for x in xs if x is not None]
    return round(sum(xs) / len(xs), 5) if xs else None


def main(argv=None):
    ids = (argv if argv is not None else sys.argv[1:]) or DEFAULT_IDS
    results = []
    for bid in ids:
        print(f"=== {bid}: Step-8 rescue default vs relaxed A/B ===", flush=True)
        W = WORK / bid
        man = json.loads((W / "subset" / "subset.manifest.json").read_text())
        p = man["paths"]
        root = W / "_step8_rescue_ab"
        root.mkdir(parents=True, exist_ok=True)

        outs, prs = {}, {}
        for state, flags in STATES.items():
            print(f"--- {bid} state {state} ---", flush=True)
            outs[state], prs[state] = _run_state(bid, state, flags, p, root)

        ref, ref_index = evaluator.build_reference(p["ref_gff"], p["ref_fa"], log=print)
        eval_dir = root / "eval"
        eval_dir.mkdir(parents=True, exist_ok=True)
        for state in STATES:
            evaluator.evaluate_tool(state, str(outs[state]), p["tgt_fa"], ref, man,
                                    eval_dir, None, log=print, ref_index=ref_index, threads=8)

        pi = {s: _pi_by_ref(eval_dir / f"{s}.transcripts.tsv") for s in STATES}
        added = set(pi["relaxed"]) - set(pi["default"])
        common = set(pi["relaxed"]) & set(pi["default"])
        regressed = sum(1 for k in common if pi["relaxed"][k] < pi["default"][k] - 1e-9)
        improved = sum(1 for k in common if pi["relaxed"][k] > pi["default"][k] + 1e-9)
        added_pis = [pi["relaxed"][k] for k in added]
        mean_pi_added = _mean(added_pis)
        frac_added_ge_07 = (round(sum(1 for v in added_pis if v >= 0.7) / len(added_pis), 3)
                            if added_pis else None)
        val_def = _validate(outs["default"])
        val_rel = _validate(outs["relaxed"])

        gate_pass = bool(
            len(added) > 0 and mean_pi_added is not None and mean_pi_added >= PI_BAR
            and regressed == 0 and val_rel[1] <= val_def[1])

        rec = {
            "benchmark": bid,
            "n_scored": {s: len(pi[s]) for s in STATES},
            "mean_pi": {s: _mean(pi[s].values()) for s in STATES},
            "n_added": len(added),
            "mean_pi_added": mean_pi_added,
            "frac_added_ge_0.7": frac_added_ge_07,
            "added_pis": sorted(round(v, 3) for v in added_pis),
            "common_improved": improved,
            "common_regressed": regressed,
            "validate": {"default_errs": val_def[1], "relaxed_errs": val_rel[1]},
            "wall_s": {s: round(prs[s].wall_clock_seconds, 1) for s in STATES},
            "gate_pass": gate_pass,
        }
        results.append(rec)
        print(f"\n  [{bid}] n_added={len(added)} mean_pi_added={mean_pi_added} "
              f"(frac≥0.7={frac_added_ge_07}); added PIs={rec['added_pis']}", flush=True)
        print(f"  common improved/regressed = {improved}/{regressed}; "
              f"validity def/relax errs = {val_def[1]}/{val_rel[1]}; "
              f"gate {'PASS' if gate_pass else 'FAIL'}", flush=True)

    (HERE / "step8_rescue_ab.json").write_text(json.dumps(results, indent=2))
    _write_md(results, HERE / "step8_rescue_ab.md")
    n_pass = sum(1 for r in results if r["gate_pass"])
    print(f"\nGATE: {n_pass}/{len(results)} datasets pass (need ALL to promote "
          f"relaxed defaults).", flush=True)
    print("DONE: wrote step8_rescue_ab.json + step8_rescue_ab.md", flush=True)
    return 0


def _write_md(results, path: Path):
    lines = [
        "## Step-8 rescue: default vs relaxed thresholds A/B — Iteration 13\n",
        "`relaxed` = `-overlap 0.2 -min_miniprot 0.8 -max_miniprot 1.6` (values "
        "from the headroom near-miss buckets); cached `-L`/`-M`, scored by the "
        "independent re-alignment evaluator. **n_added** = transcripts the relaxed "
        "filters newly emit; **mean PI added** = their protein identity (legit "
        "rescues high, garbage low — the false-rescue discriminator). Gate: "
        f"n_added>0 AND mean PI added ≥ {PI_BAR} AND 0 regressed AND no new "
        "validity errors, on BOTH datasets.\n",
        "| Dataset | n_added | mean PI added | frac ≥0.7 | added PIs | common impr/regr | validity def/relax | gate |",
        "|---|---|---|---|---|---|---|---|",
    ]
    for r in results:
        v = r["validate"]
        lines.append("| {} | {} | {} | {} | {} | {}/{} | {}/{} | {} |".format(
            r["benchmark"], r["n_added"], r["mean_pi_added"], r["frac_added_ge_0.7"],
            r["added_pis"], r["common_improved"], r["common_regressed"],
            v["default_errs"], v["relaxed_errs"],
            "**PASS**" if r["gate_pass"] else "NO-GO"))
    path.write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    sys.exit(main())
