#!/usr/bin/env python
"""A/B for the `--native` fresh-Liftoff fix (Iteration 7).

Before the fix, a FRESH `--native` run mapped NOTHING: the mappy facade built
its `Aligner` without translating Liftoff's `mm2_options`, so `--eqx` was
dropped, mappy emitted plain `M` CIGAR, and the downstream Liftoff parser
(which only counts `=`/`X` ops 7/8) saw zero aligned bases → empty
`liftoff.gff3` → exit 1. The fix (`minimap_facade._translate_mm2_options`)
passes `extra_flags |= MM_F_EQX` (+ `best_n` from `-N`).

This A/B proves the native path is now a valid drop-in for the subprocess
minimap2 path. Each state runs LiftOn FRESH (no cached `-L`/`-M`):

  state "subprocess"  (no flag)   Liftoff drives minimap2 as a subprocess  ← baseline
  state "native"      --native    Liftoff drives minimap2 via mappy in-process

Reports per dataset: gene/mRNA/CDS completeness each, shared-mRNA count, and
mean protein-identity Δ on shared transcripts. Byte-identity is NOT the gate
(mappy ≠ subprocess minimap2, and Liftoff's `-copies` is non-deterministic
across fresh runs anyway). The fix is good iff native maps a comparable gene
set (completeness ≈ subprocess) at ≈equal identity — vs the pre-fix ZERO.

Usage (repo root, lifton_devel env):
    python -m benchmarks.compare.native_liftoff_ab [IDS...]   # default below
"""
from __future__ import annotations

import json
import os
import sys
from pathlib import Path

from .profiling import run_profiled
from .tool_runners import _clean_input_dbs, _compose_env
from .align_window_ab import _ann_db, _mrna_blocks, _TRANSCRIPT_TYPES

HERE = Path(__file__).resolve().parent
WORK = HERE / "work"
REG = json.loads((HERE / "benchmarks.json").read_text())
TOOLS = REG["tools"]

DEFAULT_IDS = ["bee", "drosophila", "mouse_to_rat"]
THREADS = os.environ.get("LIFTON_AB_THREADS", "8")

STATES = {
    "subprocess": [],
    "native": ["--native"],
}


def _run(bid, state, flag, p, root):
    out_dir = root / state
    out_dir.mkdir(parents=True, exist_ok=True)
    out = out_dir / f"{state}.gff3"
    _clean_input_dbs(p["ref_gff"])
    argv = [TOOLS["lifton_bin"], "-t", THREADS, "-copies", "-ad", _ann_db(bid),
            "-g", p["ref_gff"],
            "-o", str(out), *flag, p["tgt_fa"], p["ref_fa"]]
    pr = run_profiled(argv, label=f"native_liftoff_{state}",
                      log_dir=root / "logs", env=_compose_env(TOOLS), log=print)
    if pr.exit_code != 0 or not out.exists() or out.stat().st_size == 0:
        raise RuntimeError(f"{bid}: state {state} failed (exit {pr.exit_code}); "
                           f"see {pr.stderr_path}")
    return out, pr


def _feature_counts(gff_path):
    n_gene = n_tx = n_cds = 0
    for ln in Path(gff_path).read_text().splitlines():
        if not ln or ln.startswith("#"):
            continue
        c = ln.split("\t")
        if len(c) < 9:
            continue
        t = c[2]
        if t == "gene" or t.endswith("_gene") or t == "pseudogene":
            n_gene += 1
        elif t in _TRANSCRIPT_TYPES:
            n_tx += 1
        elif t == "CDS":
            n_cds += 1
    return {"gene": n_gene, "transcript": n_tx, "CDS": n_cds}


def _identity_summary(spid, npid):
    shared = sorted(set(spid) & set(npid))
    deltas = [npid[t] - spid[t] for t in shared]
    mean = (sum(deltas) / len(deltas)) if deltas else 0.0
    regressed = sum(1 for d in deltas if d < -1e-3)
    improved = sum(1 for d in deltas if d > 1e-3)
    return {
        "n_sub": len(spid), "n_native": len(npid),
        "n_shared": len(shared),
        "native_only": len(set(npid) - set(spid)),
        "sub_only": len(set(spid) - set(npid)),
        "mean_id_delta": round(mean, 6),
        "improved": improved, "regressed": regressed,
    }


def main(argv=None):
    ids = (argv if argv is not None else sys.argv[1:]) or DEFAULT_IDS
    results = []
    for bid in ids:
        print(f"=== {bid}: --native vs subprocess Liftoff A/B "
              f"(-t {THREADS}, FRESH Step 4) ===", flush=True)
        man = json.loads((WORK / bid / "subset" / "subset.manifest.json").read_text())
        p = man["paths"]
        root = WORK / bid / "_native_liftoff_ab"
        root.mkdir(parents=True, exist_ok=True)

        outs, prs = {}, {}
        for state, flag in STATES.items():
            print(f"--- {bid} state {state} ---", flush=True)
            outs[state], prs[state] = _run(bid, state, flag, p, root)

        cnt_s = _feature_counts(outs["subprocess"])
        cnt_n = _feature_counts(outs["native"])
        _, spid = _mrna_blocks(outs["subprocess"])
        _, npid = _mrna_blocks(outs["native"])
        idsum = _identity_summary(spid, npid)

        gene_ratio = round(cnt_n["gene"] / max(cnt_s["gene"], 1), 4)
        rec = {
            "benchmark": bid,
            "counts": {"subprocess": cnt_s, "native": cnt_n},
            "native_gene_ratio": gene_ratio,
            "identity": idsum,
            "wall_s": {s: round(prs[s].wall_clock_seconds, 1) for s in STATES},
        }
        results.append(rec)
        print(f"\n  [{bid}] genes sub={cnt_s['gene']} native={cnt_n['gene']} "
              f"(ratio {gene_ratio}) | mRNA sub={cnt_s['transcript']} "
              f"native={cnt_n['transcript']} | CDS sub={cnt_s['CDS']} native={cnt_n['CDS']}",
              flush=True)
        print(f"  shared mRNA={idsum['n_shared']} (native-only={idsum['native_only']}, "
              f"sub-only={idsum['sub_only']}) | mean protein-id Δ(native-sub)="
              f"{idsum['mean_id_delta']:+f} (improved={idsum['improved']}, "
              f"regressed={idsum['regressed']})", flush=True)
        print(f"  wall: {rec['wall_s']}", flush=True)

    (HERE / "native_liftoff_ab.json").write_text(json.dumps(results, indent=2))
    _write_md(results, HERE / "native_liftoff_ab.md")
    print("\nDONE: wrote native_liftoff_ab.json + native_liftoff_ab.md", flush=True)
    return 0


def _write_md(results, path):
    lines = [
        "## --native vs subprocess Liftoff A/B — Iteration 7 (fresh-Liftoff fix)\n",
        "Each state runs LiftOn **fresh** (no cached `-L`/`-M`). **subprocess** "
        "= Liftoff drives minimap2 as a subprocess (baseline); **native** = "
        "`--native` drives minimap2 via mappy in-process. Before the "
        "`mm2_options`/`--eqx` fix, native mapped **0** genes; the gate is that "
        "it now maps a comparable gene set at ≈equal protein identity. "
        "Byte-identity is NOT expected (mappy ≠ subprocess minimap2; Liftoff "
        "`-copies` is non-deterministic across runs).\n",
        f"Threads: `-t {THREADS}`.\n",
        "| Dataset | genes sub→native (ratio) | mRNA sub→native | CDS sub→native | shared mRNA | mean protein-id Δ | improved/regressed | wall sub→native |",
        "|---|---|---|---|---|---|---|---|",
    ]
    for r in results:
        c, i, w = r["counts"], r["identity"], r["wall_s"]
        lines.append("| {} | {}→{} ({}) | {}→{} | {}→{} | {} | {:+f} | {}/{} | {}s→{}s |".format(
            r["benchmark"], c["subprocess"]["gene"], c["native"]["gene"],
            r["native_gene_ratio"], c["subprocess"]["transcript"],
            c["native"]["transcript"], c["subprocess"]["CDS"], c["native"]["CDS"],
            i["n_shared"], i["mean_id_delta"], i["improved"], i["regressed"],
            w["subprocess"], w["native"]))
    path.write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    sys.exit(main())
