#!/usr/bin/env python
"""candidate-3 A/B (Phase C) — the 3-way miniprot-only best-of-outcome candidate.

Tests `run_liftoff.process_liftoff_with_protein`'s new 3rd candidate (emit
miniprot's native model when STRICTLY better than max(Liftoff+ORF, merge+ORF)),
on the FULL human->chicken pair where 83% of the Figure-4 deficit is SAME-LOCUS
(so has_valid_miniprot=True and the candidate can fire). Step-7-only via the
-T/-P cached-FASTA bypass (sidesteps the Step-3 env blocker) + cached -L/-M.

Arms (both -t8 --locus-pipeline, NO -copies = deterministic):
  off : LIFTON_MINIPROT_CANDIDATE=0 (the 2-way best-of-outcome — pre-fix)
  on  : default (3-way)

Gate (the candidate is additive by construction): off subset of on (n_lost==0),
n_improved>0, n_regressed==0, completeness/validity not worse, and the
apples-to-apples deficit vs miniprot shrinks (re-scored common-set delta).

Run (repo root, lifton_devel env):  python -m benchmarks.compare.miniprot_candidate_ab
"""
from __future__ import annotations
import csv, json, os, subprocess, sys
from pathlib import Path
from . import evaluator
from .profiling import run_profiled
from .tool_runners import _compose_env

HERE = Path(__file__).resolve().parent
WORK = HERE / "work"
TOOLS = json.loads((HERE / "benchmarks.json").read_text())["tools"]

# FULL-run inputs (the outliers live in _fourway_full, not the subset).
# Config comes from benchmarks.json so ref/target match the cached -L/-M/-T/-P
# EXACTLY (an earlier version hardcoded a LiftOn-OUTPUT gff as the reference and
# an NC_-accession genome that mismatched the chr-style RefSeq seqids -> the
# ref-db build crashed; the real ref.gff is a symlink to NCBI_RefSeq_no_rRNA.gff
# with a valid cached .gff_db, so lifton opens it directly and never rebuilds).
PAIR = "t4_human_to_chicken"
_CFG = next(c for c in
            __import__("json").loads((HERE / "benchmarks.json").read_text())["benchmarks"]
            if c.get("id") == PAIR)
REF_GFF = Path(_CFG["ref_gff"])
REF_FA = Path(_CFG["ref_genome"])
TGT_FA = Path(_CFG["tgt_genome"])
FULL = WORK / PAIR / "_fourway_full" / "devel" / "lifton_output"
L = FULL / "liftoff" / "liftoff.gff3"
M = FULL / "miniprot" / "miniprot.gff3"
T = FULL / "intermediate_files" / "transcripts.fa"
P = FULL / "intermediate_files" / "proteins.fa"
MINI_TSV = WORK / PAIR / "_fourway_full" / "eval" / "miniprot.transcripts.tsv"


def _run(arm, root):
    out = root / arm / f"{arm}.gff3"
    out.parent.mkdir(parents=True, exist_ok=True)
    # Reuse a completed arm (the lift is the slow part; the harness can crash
    # later in scoring). A fully-written gff3 has its trailing report; treat a
    # >1 MB file as complete for these full-genome runs.
    if out.exists() and out.stat().st_size > 1_000_000:
        print(f"[bench] reusing existing {arm}.gff3 ({out.stat().st_size} bytes)", flush=True)
        return out
    argv = [TOOLS["lifton_bin"], "-t", "8", "--locus-pipeline", "-ad", "RefSeq",
            "-g", str(REF_GFF), "-L", str(L), "-M", str(M),
            "-T", str(T), "-P", str(P), "-o", str(out),
            str(TGT_FA), str(REF_FA)]
    env = dict(_compose_env(TOOLS))
    if arm == "off":
        env["LIFTON_MINIPROT_CANDIDATE"] = "0"
    else:
        env.pop("LIFTON_MINIPROT_CANDIDATE", None)
    pr = run_profiled(argv, label=f"cand3_ab_{arm}", log_dir=root / "logs", env=env, log=print)
    if pr.exit_code != 0 or not out.exists() or out.stat().st_size == 0:
        raise RuntimeError(f"arm {arm} failed (exit {pr.exit_code}); see {pr.stderr_path}")
    return out


def _pi_by_ref(tsv: Path) -> dict:
    out = {}
    if not tsv.exists():
        return out
    for r in csv.DictReader(tsv.open(), delimiter="\t"):
        if str(r.get("is_coding")).strip().lower() not in ("1", "true", "yes"):
            continue
        pi = r.get("protein_identity")
        if pi in (None, "", "None"):
            continue
        try:
            v = float(pi)
        except ValueError:
            continue
        rid = r["ref_mrna_id"]
        out[rid] = v if rid not in out else max(out[rid], v)
    return out


def _ntag(gff: Path) -> int:
    n = 0
    for ln in gff.open():
        if not ln.startswith("#"):
            c = ln.split("\t")
            if len(c) > 8 and c[2] == "mRNA" and "LiftOn_miniprot" in c[8]:
                n += 1
    return n


def _mean(xs):
    xs = [x for x in xs if x is not None]
    return round(sum(xs) / len(xs), 5) if xs else None


def main():
    root = WORK / PAIR / "_cand3_ab"
    root.mkdir(parents=True, exist_ok=True)
    for f in (REF_GFF, REF_FA, TGT_FA, L, M, T, P):
        if not Path(f).exists():
            sys.exit(f"missing input: {f}")
    man = {"id": PAIR, "species": _CFG.get("species", PAIR),
           "cross_species": _CFG.get("cross_species", True),
           "paths": {"ref_gff": str(REF_GFF), "ref_fa": str(REF_FA), "tgt_fa": str(TGT_FA)},
           "miniprot_target_space": "transcript", "protein_acc_to_mrna": {}}
    print("=== build reference (full human) ===", flush=True)
    ref, ref_index = evaluator.build_reference(str(REF_GFF), str(REF_FA), log=print)
    eval_dir = root / "eval"; eval_dir.mkdir(parents=True, exist_ok=True)
    pis = {}
    for arm in ("off", "on"):
        print(f"=== run arm {arm} ===", flush=True)
        out = _run(arm, root)
        evaluator.evaluate_tool(arm, str(out), str(TGT_FA), ref, man, eval_dir, None,
                                log=print, ref_index=ref_index, threads=8)
        pis[arm] = _pi_by_ref(eval_dir / f"{arm}.transcripts.tsv")
    mini = _pi_by_ref(MINI_TSV)

    off, on = pis["off"], pis["on"]
    common = set(off) & set(on)
    improved = sum(1 for k in common if on[k] > off[k] + 1e-9)
    regressed = sum(1 for k in common if on[k] < off[k] - 1e-9)
    lost = set(off) - set(on)
    # apples-to-apples vs miniprot on the common-with-miniprot set
    cm_off = set(off) & set(mini); cm_on = set(on) & set(mini)
    a2a_off = _mean([off[k] - mini[k] for k in cm_off])
    a2a_on = _mean([on[k] - mini[k] for k in cm_on])
    rec = {
        "pair": PAIR, "n_off": len(off), "n_on": len(on),
        "n_common": len(common), "n_improved": improved, "n_regressed": regressed,
        "n_lost(off-minus-on)": len(lost),
        "mean_pi_off": _mean(list(off.values())), "mean_pi_on": _mean(list(on.values())),
        "apples_to_apples_vs_miniprot_off": a2a_off,
        "apples_to_apples_vs_miniprot_on": a2a_on,
        "deficit_shrunk": (a2a_on is not None and a2a_off is not None and a2a_on > a2a_off),
        "n_tagged_LiftOn_miniprot": _ntag(root / "on" / "on.gff3"),
        "GATE_PASS": bool(improved > 0 and regressed == 0 and len(lost) == 0
                          and a2a_on is not None and a2a_off is not None and a2a_on > a2a_off),
    }
    (HERE / "miniprot_candidate_ab.json").write_text(json.dumps(rec, indent=2))
    print("\n=== candidate-3 A/B (human->chicken) ===")
    print(json.dumps(rec, indent=2))
    print("done:GATE_PASS=" + str(rec["GATE_PASS"]))


if __name__ == "__main__":
    main()
