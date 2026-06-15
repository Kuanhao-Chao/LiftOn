#!/usr/bin/env python
"""Iteration-15 A/B — weak-Liftoff miniprot rescue: default vs --miniprot-rescue-weak-liftoff (Gap #2).

AUDIT RECORD (Iteration-15 NO-GO, decisive). This A/B is what turned the
promising Gap-#2 headroom into a NO-GO. Result: drosophila **n_added=5** new
genes (mean PI 0.905, 0 regressed, 0 lost, no new validity errors) but
**44 extra gene models** (mostly DUPLICATES of weakly-lifted ref genes — the
"36 common improved" is an artifact of the max-per-ref evaluator crediting an
additive redundant model); mouse_to_rat **n_added=0** (all 8 rescues duplicate
already-lifted genes). mean-PI deltas +0.0009 / +0.0025, below the +0.003 bar.
The genuine missing-gene recovery is negligible and the dominant effect is
redundant overlapping models, so the `--miniprot-rescue-weak-liftoff` feature
was **REVERTED**. The flag is gone; this script needs it restored to re-run. Its
committed `miniprot_only_rescue_ab.{json,md}` are the audit trail. Jointly with
the Gap-#1 headroom and Iter-13, this confirms Step-8 suppression correctly drops
redundant duplicates. See the CLAUDE.md Iteration-15 note + memory.

Confirms on EMITTED output what the headroom diagnostic
(`miniprot_only_rescue_headroom.py`) projects: does recovering an overlapping
miniprot mRNA at a genuinely WEAK Liftoff locus add real, high-identity gene
models — or just duplicate/garbage? Distinct from the Iter-13 threshold A/B
(`step8_rescue_ab.py`): this is the Liftoff-QUALITY-aware gate, toggled by the
new `--miniprot-rescue-weak-liftoff` flag. Both states reuse cached `-L`/`-M`.

  state "default"  (no flag)
  state "weak_rescue"  --miniprot-rescue-weak-liftoff

Measures, scored by the independent re-alignment evaluator (keyed by
ref_mrna_id, NOT LiftOn's own score):
  - n_added       : coding ref transcripts emitted by weak_rescue but NOT default
                    (genuine completeness — Liftoff never produced them).
  - mean_pi_added : their protein identity (legit high, garbage low).
  - n_lost        : ref transcripts in default but NOT weak_rescue (the off ⊆ on
                    completeness gate — MUST be 0).
  - common improved/regressed (a rescue can raise the best PI for a ref Liftoff
                    lifted weakly — captured here if the evaluator keeps the better).
  - n_extra_gene_models : raw `\tgene\t` line delta (catches duplicate/overlapping
                    models when a rescue's ref id was already lifted: extra models
                    >> n_added means duplicates, not new genes).
  - validity (gff3_validator error-line delta).

Decision gate (promote — as opt-in, given mouse_to_rat headroom was marginal):
  n_added > 0 AND mean_pi_added >= PI_BAR AND n_lost == 0 AND common regressed == 0
  AND no new validity errors, on at least one dataset; a clean both-dataset pass
  argues for default-on, a drosophila-only pass for opt-in. Net duplicate-heavy or
  validity-degrading -> NO-GO (revert, keep audit).

Usage (repo root, lifton_devel env):
    python -m benchmarks.compare.miniprot_only_rescue_ab [IDS...]
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
    "weak_rescue": ["--miniprot-rescue-weak-liftoff"],
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
    pr = run_profiled(argv, label=f"mp_rescue_ab_{state}",
                      log_dir=root / "logs", env=_compose_env(TOOLS), log=print)
    if pr.exit_code != 0 or not out.exists() or out.stat().st_size == 0:
        raise RuntimeError(f"{bid}: state {state} failed (exit {pr.exit_code}); "
                           f"see {pr.stderr_path}")
    return out, pr


def _pi_by_ref(tsv: Path) -> dict:
    """ref_mrna_id -> best coding protein identity (max if duplicated)."""
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
                v = float(pi)
            except ValueError:
                continue
            rid = row["ref_mrna_id"]
            out[rid] = v if rid not in out else max(out[rid], v)
    return out


def _n_gene_lines(gff: Path) -> int:
    n = 0
    with gff.open() as fh:
        for ln in fh:
            if ln.startswith("#"):
                continue
            c = ln.split("\t")
            if len(c) > 2 and c[2] == "gene":
                n += 1
    return n


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
        print(f"=== {bid}: weak-Liftoff miniprot-rescue A/B (Gap #2) ===", flush=True)
        W = WORK / bid
        man = json.loads((W / "subset" / "subset.manifest.json").read_text())
        p = man["paths"]
        root = W / "_mp_rescue_ab"
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
        added = set(pi["weak_rescue"]) - set(pi["default"])
        lost = set(pi["default"]) - set(pi["weak_rescue"])
        common = set(pi["weak_rescue"]) & set(pi["default"])
        regressed = sum(1 for k in common if pi["weak_rescue"][k] < pi["default"][k] - 1e-9)
        improved = sum(1 for k in common if pi["weak_rescue"][k] > pi["default"][k] + 1e-9)
        added_pis = [pi["weak_rescue"][k] for k in added]
        mean_pi_added = _mean(added_pis)
        frac_added_ge_07 = (round(sum(1 for v in added_pis if v >= 0.7) / len(added_pis), 3)
                            if added_pis else None)
        n_genes = {s: _n_gene_lines(outs[s]) for s in STATES}
        n_extra_gene_models = n_genes["weak_rescue"] - n_genes["default"]
        val_def = _validate(outs["default"])
        val_weak = _validate(outs["weak_rescue"])

        gate_pass = bool(
            len(added) > 0 and mean_pi_added is not None and mean_pi_added >= PI_BAR
            and len(lost) == 0 and regressed == 0 and val_weak[1] <= val_def[1])

        rec = {
            "benchmark": bid,
            "n_scored": {s: len(pi[s]) for s in STATES},
            "mean_pi": {s: _mean(pi[s].values()) for s in STATES},
            "n_added": len(added),
            "mean_pi_added": mean_pi_added,
            "frac_added_ge_0.7": frac_added_ge_07,
            "added_pis": sorted(round(v, 3) for v in added_pis),
            "n_lost": len(lost),
            "common_improved": improved,
            "common_regressed": regressed,
            "n_gene_lines": n_genes,
            "n_extra_gene_models": n_extra_gene_models,
            "validate": {"default_errs": val_def[1], "weak_errs": val_weak[1]},
            "wall_s": {s: round(prs[s].wall_clock_seconds, 1) for s in STATES},
            "gate_pass": gate_pass,
        }
        results.append(rec)
        print(f"\n  [{bid}] n_added={len(added)} mean_pi_added={mean_pi_added} "
              f"(frac>=0.7={frac_added_ge_07}); added PIs={rec['added_pis']}", flush=True)
        print(f"  n_lost={len(lost)} common impr/regr={improved}/{regressed}; "
              f"extra gene models={n_extra_gene_models} "
              f"(genes {n_genes['default']}->{n_genes['weak_rescue']})", flush=True)
        print(f"  validity def/weak errs = {val_def[1]}/{val_weak[1]}; "
              f"gate {'PASS' if gate_pass else 'FAIL'}", flush=True)

    (HERE / "miniprot_only_rescue_ab.json").write_text(json.dumps(results, indent=2))
    _write_md(results, HERE / "miniprot_only_rescue_ab.md")
    n_pass = sum(1 for r in results if r["gate_pass"])
    print(f"\nGATE: {n_pass}/{len(results)} datasets pass. ALL pass -> default-on; "
          f"drosophila-only -> opt-in; none -> NO-GO.", flush=True)
    print("DONE: wrote miniprot_only_rescue_ab.json + miniprot_only_rescue_ab.md",
          flush=True)
    return 0


def _write_md(results, path: Path):
    lines = [
        "## Weak-Liftoff miniprot-rescue A/B — Iteration 15 (Gap #2)\n",
        "`weak_rescue` = `--miniprot-rescue-weak-liftoff`; cached `-L`/`-M`, scored "
        "by the independent re-alignment evaluator (keyed by ref_mrna_id). "
        "**n_added** = coding ref transcripts weak_rescue emits but default does "
        "not (genuine completeness); **mean PI added** = their identity. "
        "**n_lost** = the off ⊆ on gate (must be 0). **extra gene models** = raw "
        "gene-line delta (>> n_added ⇒ duplicate/overlapping models, a quality "
        f"concern). Gate: n_added>0 AND mean PI added ≥ {PI_BAR} AND n_lost=0 AND "
        "0 regressed AND no new validity errors.\n",
        "| Dataset | n_added | mean PI added | frac ≥0.7 | n_lost | common impr/regr | extra gene models | validity def/weak | gate |",
        "|---|---|---|---|---|---|---|---|---|",
    ]
    for r in results:
        v = r["validate"]
        lines.append("| {} | {} | {} | {} | {} | {}/{} | {} | {}/{} | {} |".format(
            r["benchmark"], r["n_added"], r["mean_pi_added"], r["frac_added_ge_0.7"],
            r["n_lost"], r["common_improved"], r["common_regressed"],
            r["n_extra_gene_models"], v["default_errs"], v["weak_errs"],
            "**PASS**" if r["gate_pass"] else "NO-GO"))
    path.write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    sys.exit(main())
