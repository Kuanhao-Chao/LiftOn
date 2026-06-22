#!/usr/bin/env python
"""A/B for the "band everything" fast-align kernel (Iteration 3).

Runs LiftOn on a benchmark TWICE on the same cached `-L`/`-M`, differing only by
the alignment-kernel mode:
  state "default"  giant-only windowing (Iteration 2)            byte-frozen path
  state "fast"     LIFTON_FAST_ALIGN=1                           band everything

Unlike `align_window_ab.py` (which isolates the giant subset), this compares
EVERY transcript — fast-align deliberately changes the mid-size tail too — and
reports the decision-gate signals:
  - Δwall, Δpeak-RSS (speed/memory win),
  - mean & worst per-transcript Δprotein-identity (fast − default), with
    improved/regressed counts at the 1e-3 material threshold (accuracy),
  - completeness: gene/transcript/CDS feature counts must be unchanged.

Promotion criteria (the loop's decision gate): speedup ≥ ~5% wall AND mean
identity Δ ≥ −1e-3 AND completeness unchanged. This script only REPORTS the Δ;
promotion is a human decision after reading it.

Usage (repo root, lifton_devel env):
    python -m benchmarks.compare.fast_align_ab [IDS...]   # default = all three
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

from .profiling import run_profiled
from .tool_runners import _clean_input_dbs, _compose_env
from .align_window_ab import _mrna_blocks, _ann_db, _TRANSCRIPT_TYPES

HERE = Path(__file__).resolve().parent
WORK = HERE / "work"
REG = json.loads((HERE / "benchmarks.json").read_text())
TOOLS = REG["tools"]

DEFAULT_IDS = ["mouse", "drosophila", "mouse_to_rat"]
MATERIAL = 1e-3          # per-transcript identity Δ that counts as a real change

STATES = {
    "default": {},
    "fast": {"LIFTON_FAST_ALIGN": "1"},
}


def _run(bid, state, env_extra, p, root):
    out_dir = root / state
    out_dir.mkdir(parents=True, exist_ok=True)
    out = out_dir / f"{state}.gff3"
    liftoff = WORK / bid / "tools" / "liftoff" / "liftoff.gff3"
    miniprot = WORK / bid / "tools" / "miniprot" / "miniprot.gff3"
    argv = [TOOLS["lifton_bin"], "-t", "1", "-copies", "-ad", _ann_db(bid),
            "-g", p["ref_gff"], "-L", str(liftoff), "-M", str(miniprot),
            "-o", str(out), p["tgt_fa"], p["ref_fa"]]
    env = _compose_env(TOOLS)
    env.update(env_extra)
    pr = run_profiled(argv, label=f"fast_align_{state}",
                      log_dir=root / "logs", env=env, log=print)
    if pr.exit_code != 0 or not out.exists() or out.stat().st_size == 0:
        raise RuntimeError(f"{bid}: state {state} failed (exit {pr.exit_code}); "
                           f"see {pr.stderr_path}")
    return out, pr


def _feature_counts(gff_path):
    """Top-line completeness: count gene / transcript / CDS feature lines."""
    n_gene = n_tx = n_cds = n_total = 0
    for ln in Path(gff_path).read_text().splitlines():
        if not ln or ln.startswith("#"):
            continue
        c = ln.split("\t")
        if len(c) < 9:
            continue
        n_total += 1
        t = c[2]
        if t == "gene" or t.endswith("_gene") or t == "pseudogene":
            n_gene += 1
        elif t in _TRANSCRIPT_TYPES:
            n_tx += 1
        elif t == "CDS":
            n_cds += 1
    return {"gene": n_gene, "transcript": n_tx, "CDS": n_cds, "total": n_total}


def _identity_delta(dpid, fpid):
    """Per-transcript protein_identity delta (fast − default) over transcripts
    scored in BOTH states."""
    shared = sorted(set(dpid) & set(fpid))
    deltas = [(t, round(fpid[t] - dpid[t], 6)) for t in shared]
    nz = [d for _, d in deltas if abs(d) >= 1e-9]
    improved = [d for d in nz if d > MATERIAL]
    regressed = [(t, d) for t, d in deltas if d < -MATERIAL]
    mean = (sum(d for _, d in deltas) / len(deltas)) if deltas else 0.0
    worst = min(deltas, key=lambda kv: kv[1]) if deltas else None
    best = max(deltas, key=lambda kv: kv[1]) if deltas else None
    return {
        "n_scored_both": len(deltas),
        "n_changed": len(nz),
        "n_improved_gt_1e-3": len(improved),
        "n_regressed_gt_1e-3": len(regressed),
        "mean_delta": round(mean, 6),
        "worst": worst,
        "best": best,
        "regressed_ids_top": sorted(regressed, key=lambda kv: kv[1])[:10],
    }


def main(argv=None):
    ids = (argv if argv is not None else sys.argv[1:]) or DEFAULT_IDS
    results = []
    for bid in ids:
        print(f"=== {bid}: fast-align (band everything) vs default A/B ===", flush=True)
        man = json.loads((WORK / bid / "subset" / "subset.manifest.json").read_text())
        p = man["paths"]
        root = WORK / bid / "_fast_align_ab"
        root.mkdir(parents=True, exist_ok=True)
        _clean_input_dbs(p["ref_gff"],
                         WORK / bid / "tools" / "liftoff" / "liftoff.gff3",
                         WORK / bid / "tools" / "miniprot" / "miniprot.gff3")

        outs, prs = {}, {}
        for state, env_extra in STATES.items():
            print(f"--- {bid} state {state} ---", flush=True)
            outs[state], prs[state] = _run(bid, state, env_extra, p, root)

        _, dpid = _mrna_blocks(outs["default"])
        _, fpid = _mrna_blocks(outs["fast"])
        idelta = _identity_delta(dpid, fpid)
        counts_d = _feature_counts(outs["default"])
        counts_f = _feature_counts(outs["fast"])
        completeness_ok = (counts_d == counts_f)

        wall_d, wall_f = prs["default"].wall_clock_seconds, prs["fast"].wall_clock_seconds
        rss_d = getattr(prs["default"], "peak_rss_mb", 0) or 0
        rss_f = getattr(prs["fast"], "peak_rss_mb", 0) or 0
        speedup = round(wall_d / max(wall_f, 1e-9), 3)
        wall_pct = round(100.0 * (wall_d - wall_f) / max(wall_d, 1e-9), 1)

        # decision gate
        gate_speed = wall_pct >= 5.0
        gate_acc = idelta["mean_delta"] >= -MATERIAL and idelta["n_regressed_gt_1e-3"] == 0
        promote = bool(gate_speed and gate_acc and completeness_ok)

        rec = {
            "benchmark": bid,
            "wall_s": {"default": round(wall_d, 1), "fast": round(wall_f, 1)},
            "wall_saved_pct": wall_pct,
            "speedup": speedup,
            "peak_rss_mb": {"default": round(rss_d, 0), "fast": round(rss_f, 0)},
            "identity": idelta,
            "completeness": {"default": counts_d, "fast": counts_f,
                             "unchanged": completeness_ok},
            "gate": {"speed_ge_5pct": gate_speed, "accuracy_safe": gate_acc,
                     "completeness_unchanged": completeness_ok, "promote": promote},
        }
        results.append(rec)
        print(f"\n  [{bid}] wall default={rec['wall_s']['default']}s "
              f"fast={rec['wall_s']['fast']}s ({wall_pct:+.1f}% / {speedup}x) | "
              f"RSS {rec['peak_rss_mb']['default']}→{rec['peak_rss_mb']['fast']}MB",
              flush=True)
        print(f"  identity: scored={idelta['n_scored_both']} changed={idelta['n_changed']} "
              f"meanΔ={idelta['mean_delta']:+f} improved={idelta['n_improved_gt_1e-3']} "
              f"regressed={idelta['n_regressed_gt_1e-3']} worst={idelta['worst']}",
              flush=True)
        print(f"  completeness unchanged={completeness_ok} | GATE promote={promote} "
              f"(speed{'✓' if gate_speed else '✗'} acc{'✓' if gate_acc else '✗'})",
              flush=True)

    (HERE / "fast_align_ab.json").write_text(json.dumps(results, indent=2))
    _write_md(results, HERE / "fast_align_ab.md")
    print("\nDONE: wrote fast_align_ab.json + fast_align_ab.md", flush=True)
    return 0


def _write_md(results, path):
    lines = [
        "## Fast-align (band everything) vs default A/B — Iteration 3\n",
        "Same cached `-L`/`-M`, single-threaded, differing only by "
        "`LIFTON_FAST_ALIGN`. Promotion gate: **speedup ≥ 5% wall AND mean "
        "identity Δ ≥ −1e-3 AND completeness unchanged**.\n",
        "| Dataset | wall default→fast | Δwall | speedup | meanΔ id | improved | regressed | worst id Δ | complete | PROMOTE |",
        "|---|---|---|---|---|---|---|---|---|---|",
    ]
    for r in results:
        w, i, g = r["wall_s"], r["identity"], r["gate"]
        lines.append(
            "| {} | {}s→{}s | {:+.1f}% | {}× | {:+f} | {} | {} | {} | {} | {} |".format(
                r["benchmark"], w["default"], w["fast"], r["wall_saved_pct"],
                r["speedup"], i["mean_delta"], i["n_improved_gt_1e-3"],
                i["n_regressed_gt_1e-3"], i["worst"],
                "yes" if r["completeness"]["unchanged"] else "NO",
                "yes" if g["promote"] else "no"))
    path.write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    sys.exit(main())
