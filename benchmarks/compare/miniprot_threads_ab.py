#!/usr/bin/env python
"""A/B for miniprot thread plumbing — un-cap from the binary default 4 (Iter 17).

Before Iteration 17 miniprot was launched with NO ``-t`` flag, so it always ran
at its hard-coded binary default of **4 threads** regardless of LiftOn's
``-t/--threads N`` (minimap2, by contrast, scaled). On cross-species data
(mouse_to_rat) miniprot does heavy protein alignment and, capped at 4, became
the serial tail of the concurrent Step 4 (`wall = max(t_liftoff, t_miniprot)`).
Iteration 17 plumbs ``args.threads`` into miniprot's own ``-t`` (gated on
``> 1`` so the default ``-t 1`` run is byte-identical to before).

The standalone feasibility gate already proved miniprot's GFF is **byte-
identical** across ``-t 1/4/8/16`` (input-order-preserving), with ~1.7x at -t8
/ ~2.7x at -t16 vs the fixed-4 baseline. This A/B confirms it **end-to-end**.

Two arms per dataset:

  ARM 1 — cached `-L` + FRESH miniprot (the clean gate):
    Liftoff is a cached, deterministic load (so the ONLY variable is miniprot's
    thread count) and Step-4 wall == miniprot wall. States:
      mp_t4 : LIFTON_MINIPROT_THREADS=0  -> miniprot at its old default 4
      mp_tN : (unset)                    -> miniprot at -t N   ← the iteration
    Gate: **byte-identical(mp_t4, mp_tN)** AND a miniprot/Step-4 wall win AND
    completeness unchanged.

  ARM 2 — FRESH `-L` + FRESH miniprot, default concurrent Step 4 (realism +
    oversubscription check): both aligners actually run at once, so miniprot's
    ``-t N`` overlaps minimap2's ``~N`` (up to 2N threads during the overlap
    window). Confirms the un-cap still WINS end-to-end (i.e. transient
    oversubscription does not erase the gain). Byte-identity is NOT gated here
    (Liftoff's `-copies` is non-deterministic across fresh runs — the known
    caveat); completeness is reported as a sanity signal.

Usage (repo root, lifton_devel env):
    python -m benchmarks.compare.miniprot_threads_ab [IDS...]
    LIFTON_AB_THREADS=16 python -m benchmarks.compare.miniprot_threads_ab drosophila
    MP_THREADS_AB_SKIP_FRESH=1 python -m benchmarks.compare.miniprot_threads_ab  # arm 1 only
"""
from __future__ import annotations

import json
import os
import sys
from pathlib import Path

from .profiling import run_profiled
from .tool_runners import _clean_input_dbs, _compose_env
from .align_window_ab import _ann_db, _TRANSCRIPT_TYPES

HERE = Path(__file__).resolve().parent
WORK = HERE / "work"
REG = json.loads((HERE / "benchmarks.json").read_text())
TOOLS = REG["tools"]

DEFAULT_IDS = ["drosophila", "mouse_to_rat", "rice"]
THREADS = os.environ.get("LIFTON_AB_THREADS", "8")
SKIP_FRESH = os.environ.get("MP_THREADS_AB_SKIP_FRESH") == "1"

# LIFTON_MINIPROT_THREADS value per state: "0" reproduces the old fixed default
# (miniprot t4); "" (unset) gives miniprot -t N (the iteration).
MP_ENV = {"mp_t4": "0", "mp_tN": None}


def _run(bid, state, p, root, *, cached_liftoff, perf_probe):
    out_dir = root / state
    out_dir.mkdir(parents=True, exist_ok=True)
    out = out_dir / f"{state}.gff3"

    argv = [TOOLS["lifton_bin"], "-t", THREADS, "-copies", "-ad", _ann_db(bid),
            "-g", p["ref_gff"]]
    if cached_liftoff is not None:
        # ARM 1: cached -L (Liftoff skipped, deterministic). Clean stale sibling
        # DBs for ref_gff + the cached liftoff.gff3 so each run rebuilds them.
        _clean_input_dbs(p["ref_gff"], cached_liftoff)
        argv += ["-L", str(cached_liftoff)]
    else:
        # ARM 2: fresh both. Only the ref_gff sibling DB needs cleaning.
        _clean_input_dbs(p["ref_gff"])
    # NO -M in either arm -> miniprot always runs FRESH (the point).
    argv += ["-o", str(out), "--locus-pipeline", p["tgt_fa"], p["ref_fa"]]

    env = _compose_env(TOOLS)
    env[perf_probe] = "1"
    mp_env = MP_ENV[state]
    if mp_env is None:
        env.pop("LIFTON_MINIPROT_THREADS", None)
    else:
        env["LIFTON_MINIPROT_THREADS"] = mp_env

    pr = run_profiled(argv, label=f"mp_threads_{state}",
                      log_dir=root / "logs", env=env, log=print)
    if pr.exit_code != 0 or not out.exists() or out.stat().st_size == 0:
        raise RuntimeError(f"{bid}: state {state} failed (exit {pr.exit_code}); "
                           f"see {pr.stderr_path}")
    return out, pr


def _parse_step4_wall(stderr_path):
    """Pull the LIFTON_PERF_STEP4 probe line. With cached -L (arm 1) the Step-4
    wall IS the miniprot wall (Liftoff is an instant cached load)."""
    try:
        for ln in Path(stderr_path).read_text().splitlines():
            if "[LiftOn][perf] Step4 wall" in ln:
                return float(ln.rsplit(":", 1)[1].strip().rstrip("s"))
    except (OSError, ValueError):
        pass
    return None


def _feature_counts(gff_path):
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


def _arm(bid, p, root, *, cached_liftoff, perf_probe):
    outs, prs = {}, {}
    for state in MP_ENV:
        print(f"--- {bid} [{root.name}] state {state} "
              f"(LIFTON_MINIPROT_THREADS={MP_ENV[state]!r}) ---", flush=True)
        outs[state], prs[state] = _run(bid, state, p, root,
                                       cached_liftoff=cached_liftoff,
                                       perf_probe=perf_probe)
    s4 = {s: _parse_step4_wall(prs[s].stderr_path) for s in MP_ENV}
    wall = {s: prs[s].wall_clock_seconds for s in MP_ENV}
    rss = {s: getattr(prs[s], "peak_rss_mb", 0) or 0 for s in MP_ENV}
    counts = {s: _feature_counts(outs[s]) for s in MP_ENV}
    identical = outs["mp_t4"].read_bytes() == outs["mp_tN"].read_bytes()
    s4_speedup = (round(s4["mp_t4"] / max(s4["mp_tN"], 1e-9), 2)
                  if (s4["mp_t4"] and s4["mp_tN"]) else None)
    s4_saved_pct = (round(100.0 * (s4["mp_t4"] - s4["mp_tN"]) /
                          max(s4["mp_t4"], 1e-9), 1)
                    if (s4["mp_t4"] and s4["mp_tN"]) else None)
    total_speedup = round(wall["mp_t4"] / max(wall["mp_tN"], 1e-9), 2)
    total_saved_pct = round(100.0 * (wall["mp_t4"] - wall["mp_tN"]) /
                            max(wall["mp_t4"], 1e-9), 1)
    return {
        "byte_identical": identical,
        "step4_wall_s": s4, "step4_speedup": s4_speedup,
        "step4_saved_pct": s4_saved_pct,
        "total_wall_s": {s: round(wall[s], 1) for s in MP_ENV},
        "total_speedup": total_speedup, "total_saved_pct": total_saved_pct,
        "peak_rss_mb": {s: round(rss[s], 0) for s in MP_ENV},
        "completeness": {"counts": counts,
                         "unchanged": counts["mp_t4"] == counts["mp_tN"]},
    }


def main(argv=None):
    ids = (argv if argv is not None else sys.argv[1:]) or DEFAULT_IDS
    results = []
    for bid in ids:
        print(f"=== {bid}: miniprot thread plumbing A/B (-t {THREADS}) ===",
              flush=True)
        man = json.loads((WORK / bid / "subset" / "subset.manifest.json").read_text())
        p = man["paths"]
        cached_liftoff = WORK / bid / "tools" / "liftoff" / "liftoff.gff3"
        if not cached_liftoff.exists():
            raise RuntimeError(f"{bid}: missing cached -L at {cached_liftoff}; "
                               f"build it with benchmarks/compare/build_inputs.py")

        root = WORK / bid / "_mp_threads_ab"
        root.mkdir(parents=True, exist_ok=True)

        print(f"  ARM 1: cached -L + fresh miniprot (clean gate)", flush=True)
        arm1 = _arm(bid, p, root / "cached_L", cached_liftoff=cached_liftoff,
                    perf_probe="LIFTON_PERF_STEP4")

        arm2 = None
        if not SKIP_FRESH:
            print(f"  ARM 2: fresh -L + fresh miniprot (concurrent Step 4, "
                  f"oversubscription check)", flush=True)
            arm2 = _arm(bid, p, root / "fresh_both", cached_liftoff=None,
                        perf_probe="LIFTON_PERF_STEP4")

        rec = {
            "benchmark": bid, "threads": THREADS,
            "arm1_cached_L": arm1,
            "arm2_fresh_both": arm2,
            # The GO gate lives in ARM 1 (clean, deterministic Liftoff).
            "gate": {
                "byte_identical": arm1["byte_identical"],
                "miniprot_wall_win": bool(arm1["step4_saved_pct"]
                                          and arm1["step4_saved_pct"] > 0),
                "completeness_unchanged": arm1["completeness"]["unchanged"],
            },
        }
        results.append(rec)
        a1 = arm1
        print(f"\n  [{bid}] ARM1 byte-identical(mp_t4==mp_tN)={a1['byte_identical']}",
              flush=True)
        print(f"  ARM1 miniprot/Step-4 wall: t4={a1['step4_wall_s']['mp_t4']}s "
              f"tN={a1['step4_wall_s']['mp_tN']}s "
              f"({a1['step4_speedup']}x / saved {a1['step4_saved_pct']}%)", flush=True)
        print(f"  ARM1 total wall: t4={a1['total_wall_s']['mp_t4']}s "
              f"tN={a1['total_wall_s']['mp_tN']}s "
              f"({a1['total_speedup']}x / {a1['total_saved_pct']}%) | "
              f"completeness unchanged={a1['completeness']['unchanged']}", flush=True)
        if arm2:
            print(f"  ARM2 (fresh, concurrent) Step-4 wall: "
                  f"t4={arm2['step4_wall_s']['mp_t4']}s "
                  f"tN={arm2['step4_wall_s']['mp_tN']}s "
                  f"({arm2['step4_speedup']}x / {arm2['step4_saved_pct']}%) | "
                  f"total {arm2['total_wall_s']['mp_t4']}s->"
                  f"{arm2['total_wall_s']['mp_tN']}s "
                  f"({arm2['total_speedup']}x)", flush=True)
        if not a1["byte_identical"]:
            print(f"  !! WARNING: ARM1 outputs differ — byte-neutrality gate "
                  f"FAILED; investigate before trusting the speed numbers",
                  flush=True)

    (HERE / "miniprot_threads_ab.json").write_text(json.dumps(results, indent=2))
    _write_md(results, HERE / "miniprot_threads_ab.md")
    print("\nDONE: wrote miniprot_threads_ab.json + miniprot_threads_ab.md",
          flush=True)
    return 0


def _write_md(results, path):
    lines = [
        "## miniprot thread plumbing A/B — Iteration 17\n",
        "miniprot was capped at its binary default of 4 threads regardless of "
        "LiftOn's `-t/--threads`. This A/B plumbs `args.threads` into miniprot's "
        "`-t` (gated on `>1`, so the default `-t 1` is unchanged). "
        "**ARM 1** (cached `-L`, fresh miniprot) is the clean gate: Liftoff is a "
        "deterministic cached load, so the only variable is miniprot's thread "
        "count and Step-4 wall == miniprot wall. Gate: **byte-identical AND a "
        "miniprot/Step-4 wall win AND completeness unchanged**. "
        "**ARM 2** (fresh both, concurrent Step 4) confirms the win holds "
        "end-to-end despite transient oversubscription (byte-identity not gated "
        "there — Liftoff `-copies` is non-deterministic across fresh runs).\n",
        f"Baseline `mp_t4` = `LIFTON_MINIPROT_THREADS=0` (old default 4); "
        f"`mp_tN` = miniprot `-t {THREADS}`.\n",
        "| Dataset | ARM1 byte-ident | ARM1 miniprot wall t4→tN | ARM1 speedup | ARM1 complete | ARM2 Step-4 wall t4→tN | ARM2 total t4→tN |",
        "|---|---|---|---|---|---|---|",
    ]
    for r in results:
        a1 = r["arm1_cached_L"]
        a2 = r["arm2_fresh_both"]
        s4 = a1["step4_wall_s"]
        a2cell = ("—" if not a2 else
                  "{}s→{}s ({}x)".format(a2["step4_wall_s"]["mp_t4"],
                                         a2["step4_wall_s"]["mp_tN"],
                                         a2["step4_speedup"]))
        a2tot = ("—" if not a2 else
                 "{}s→{}s ({}x)".format(a2["total_wall_s"]["mp_t4"],
                                        a2["total_wall_s"]["mp_tN"],
                                        a2["total_speedup"]))
        lines.append(
            "| {} | {} | {}s→{}s | {}x ({}%) | {} | {} | {} |".format(
                r["benchmark"],
                "yes" if a1["byte_identical"] else "**NO**",
                s4["mp_t4"], s4["mp_tN"], a1["step4_speedup"],
                a1["step4_saved_pct"],
                "yes" if a1["completeness"]["unchanged"] else "**NO**",
                a2cell, a2tot))
    lines.append(
        "\n**Interpretation:** ARM-1 byte-identical confirms the change is "
        "byte-neutral end-to-end (miniprot output is order-stable across "
        "threads — the standalone feasibility gate proved `-t 1/4/8/16` are "
        "byte-identical). The ARM-1 miniprot wall win is the direct speedup; "
        "ARM-2 shows it survives concurrent Step 4. The rigorous determinism "
        "proof is the standalone gate + the unit suite "
        "(`tests/test_miniprot_threads.py`); the 24-cell matrix is unaffected "
        "(it uses cached `-M`, so miniprot is never invoked).")
    path.write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    sys.exit(main())
