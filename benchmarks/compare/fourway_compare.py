#!/usr/bin/env python
"""4-way per-benchmark comparison: standalone **Liftoff** vs standalone **miniprot**
vs **LiftOn v1.0.8** (stable) vs **LiftOn devel** — all scored by the SAME
version-agnostic neutral evaluator (``benchmarks.compare.evaluator``).

The design is deliberately "LiftOn vs exactly the inputs it combined": the two
LiftOn columns run on the *same* cached ``tools/liftoff/liftoff.gff3`` +
``tools/miniprot/miniprot.gff3`` (via ``-L``/``-M``) that the standalone Liftoff
and miniprot columns are scored on. So the four columns are mutually consistent
and the deltas (``devel - best(liftoff,miniprot)``, ``devel - v1.0.8``) are fair.

Both LiftOn versions run with the **neutral common flag set only**
(``-t 1 -copies -ad <db> -g -L -M -o``) — no devel-only flags — so the v1.0.8
binary (which lacks ``--native``/``--lift-gene-like``/…) is valid and the two
versions differ only in their in-process logic. (devel ``-t 1`` is byte-identical
to its threaded path per the 24-cell matrix.)

Reuses ``version_compare`` for VERSIONS / provenance_gate / run_lift / validate_gff,
``build_inputs`` to bootstrap the cached ``-L``/``-M`` for new pairs, and
``evaluator.build_reference`` + ``evaluator.evaluate_tool`` for the 4× scoring.

Writes a dedicated ``fourway_results.json`` (NOT ``version_compare.results.json``)
so it never races the in-flight ``lifton_vc`` writer.

Usage (repo root, lifton_devel env):
    python -m benchmarks.compare.fourway_compare <id> [<id> ...]      # subset 4-way
    python -m benchmarks.compare.fourway_compare --full bee drosophila # full-genome 4-way
"""
from __future__ import annotations

import argparse
import dataclasses
import json
import os
import sys
from pathlib import Path

from . import build_inputs, evaluator, version_compare as vc
from .profiling import run_profiled
from .tool_runners import _clean_input_dbs

HERE = Path(__file__).resolve().parent
WORK = HERE / "work"
DATA = HERE.parent / "data"
REG = json.loads((HERE / "benchmarks.json").read_text())
# Per-run results file (parallel-safe): FOURWAY_RESULTS_JSON lets concurrent
# `--full` runs each write their OWN file (no race on the shared, unlocked _save
# read-modify-write); merge them into fourway_results.json afterward.
RESULTS = Path(os.environ.get("FOURWAY_RESULTS_JSON", str(HERE / "fourway_results.json")))

TOOLS = ["liftoff", "miniprot", "lifton_stable", "lifton_devel"]
# tool label -> version_compare VERSIONS key (only the two LiftOn columns map)
_LIFTON_VERSION = {"lifton_stable": "stable", "lifton_devel": "devel"}


def _bench(bid):
    for b in REG["benchmarks"]:
        if b["id"] == bid:
            return b
    raise KeyError(bid)


def _load_profile(path: Path) -> dict:
    try:
        return json.loads(Path(path).read_text())
    except Exception:
        return {}


def _best_baseline(d: dict, keys=("liftoff", "miniprot")):
    """max over the two baseline tools, ignoring None."""
    vals = [d.get(k) for k in keys if d.get(k) is not None]
    return max(vals) if vals else None


def _sub(a, b):
    return round(a - b, 6) if (a is not None and b is not None) else None


# ---------------------------------------------------------------------------
# subset 4-way
# ---------------------------------------------------------------------------

def _ensure_inputs(bid, threads=8, log=print):
    """Guarantee subset manifest + cached standalone liftoff/miniprot exist."""
    work = WORK / bid
    man_p = work / "subset" / "subset.manifest.json"
    lf = work / "tools" / "liftoff" / "liftoff.gff3"
    mp = work / "tools" / "miniprot" / "miniprot.gff3"
    if not (man_p.exists() and lf.exists() and lf.stat().st_size > 0
            and mp.exists() and mp.stat().st_size > 0):
        log(f"  [{bid}] inputs missing -> build_inputs (subset+liftoff+miniprot)")
        rc = build_inputs.main([bid, "-t", str(threads)])
        if rc != 0:
            raise RuntimeError(f"build_inputs({bid}) failed rc={rc}")
    return man_p, lf, mp


def run_subset(bid, threads_eval=8, log=print):
    bench = _bench(bid)
    anndb = bench.get("annotation_database", "RefSeq")
    man_p, liftoff_gff, miniprot_gff = _ensure_inputs(bid, log=log)
    man = json.loads(man_p.read_text())
    paths = man["paths"]
    work = WORK / bid
    ab = work / "_fourway"
    ab.mkdir(parents=True, exist_ok=True)

    # the two LiftOn versions on the SAME cached -L/-M, -t1, neutral flags.
    # Clean stale _db siblings before EACH version (each rebuilds its own in its
    # own conda env); tolerate one version crashing (e.g. a v1.0.8 robustness bug)
    # so the surviving columns still score.
    gffs = {"liftoff": liftoff_gff, "miniprot": miniprot_gff}
    profs = {"liftoff": _load_profile(work / "tools" / "logs" / "liftoff.profile.json"),
             "miniprot": _load_profile(work / "tools" / "logs" / "miniprot.profile.json")}
    for label, version in _LIFTON_VERSION.items():
        _clean_input_dbs(paths["ref_gff"], liftoff_gff, miniprot_gff)
        try:
            out_gff, pr = vc.run_lift(version, paths, anndb, 1, liftoff_gff, miniprot_gff,
                                      ab / label, devel_fast=False, log=log)
            gffs[label] = out_gff
            profs[label] = dataclasses.asdict(pr)
        except Exception as e:   # noqa: BLE001
            log(f"  !! {label} CRASHED on {bid} subset ({e}); skipping that column")

    return _score_and_record(bid, "subset", bench, man, paths, gffs, profs,
                             ab / "eval", threads_eval, log)


# ---------------------------------------------------------------------------
# full-genome 4-way (headline pairs)
# ---------------------------------------------------------------------------

def _full_paths(bid):
    """Full-genome (ref_fa, ref_gff, tgt_fa) for a headline pair. Reuses
    version_compare.FULL_INPUTS when available (bee/rice/arabidopsis), else the
    on-disk benchmarks.json paths (cross-species headlines)."""
    if bid in vc.FULL_INPUTS:
        ref_fa, ref_gff, tgt_fa = vc.FULL_INPUTS[bid]
        d = DATA / bid
        return {"ref_fa": d / ref_fa, "ref_gff": d / ref_gff, "tgt_fa": d / tgt_fa}
    b = _bench(bid)
    return {"ref_fa": Path(b["ref_genome"]), "ref_gff": Path(b["ref_gff"]),
            "tgt_fa": Path(b["tgt_genome"])}


def run_full(bid, threads_eval=8, log=print):
    bench = _bench(bid)
    anndb = bench.get("annotation_database", "RefSeq")
    paths = _full_paths(bid)
    for k, p in paths.items():
        if not Path(p).exists():
            raise RuntimeError(f"full input missing: {k}={p}")
    work = WORK / bid
    ab = work / "_fourway_full"
    ab.mkdir(parents=True, exist_ok=True)

    man = {"id": bid, "species": bench["species"],
           "cross_species": bench["cross_species"],
           "miniprot_target_space": "transcript", "protein_acc_to_mrna": {},
           "paths": {k: str(v) for k, v in paths.items()}}

    # Source the two LiftOn outputs + their INTERNAL standalone liftoff/miniprot.
    # Prefer reusing the lifton_vc "full" arm outputs (work/<id>/_version_ab/full/);
    # else run both versions end-to-end here (fresh, no -L/-M).
    vc_full = work / "_version_ab" / "full"

    def _statedir(label, tool):  # label in {stable, devel}
        reused = vc_full / label / f"{label}.gff3"
        if reused.exists() and reused.stat().st_size > 0:
            return vc_full / label, reused, None
        sd = ab / label
        _clean_input_dbs(paths["ref_gff"])   # both fresh runs share the full ref_gff
        try:
            out_gff, pr = vc.run_lift(label, paths, anndb, 8 if label == "devel" else 1,
                                      None, None, sd, devel_fast=(label == "devel"), log=log)
            return sd, out_gff, dataclasses.asdict(pr)
        except Exception as e:   # noqa: BLE001 — v1.0.8 crashes on some full genomes
            log(f"  !! {tool} CRASHED on {bid} full ({e}); column skipped, "
                f"statedir kept for its Step-4 liftoff/miniprot intermediates")
            return sd, None, None

    gffs, profs, statedirs = {}, {}, {}
    for label, tool in (("stable", "lifton_stable"), ("devel", "lifton_devel")):
        sd, out_gff, pr = _statedir(label, tool)
        statedirs[label] = sd
        if out_gff is not None and Path(out_gff).exists() and Path(out_gff).stat().st_size > 0:
            gffs[tool] = out_gff
            profs[tool] = pr or _load_profile(sd / "logs" / f"{label}.profile.json")
    # standalone liftoff/miniprot = the STABLE run's INTERNAL intermediates: stable
    # (v1.0.8) runs end-to-end with NO --stream/--inmemory-liftoff, so it writes
    # lifton_output/{liftoff,miniprot}/*.gff3 to disk in Step 4 — present EVEN if it
    # later crashes in Step 7 (the FeatureNotFoundError bug). devel with devel_fast
    # pipes them in-memory. The internal miniprot is transcript-space.
    src_sd = statedirs["stable"]
    lo = src_sd / "lifton_output" / "liftoff" / "liftoff.gff3"
    mp = src_sd / "lifton_output" / "miniprot" / "miniprot.gff3"
    if lo.exists() and lo.stat().st_size > 0:
        gffs["liftoff"] = lo; profs["liftoff"] = {}
    if mp.exists() and mp.stat().st_size > 0:
        gffs["miniprot"] = mp; profs["miniprot"] = {}
    for t in ("liftoff", "miniprot", "lifton_stable", "lifton_devel"):
        if t not in gffs:
            log(f"  [{bid}/full] note: {t} column n/a (missing/crashed)")

    return _score_and_record(bid, "full", bench, man, paths, gffs, profs,
                             ab / "eval", threads_eval, log)


# ---------------------------------------------------------------------------
# shared scoring + record assembly
# ---------------------------------------------------------------------------

def _score_and_record(bid, mode, bench, man, paths, gffs, profs, eval_dir,
                      threads_eval, log):
    ref, ref_index = evaluator.build_reference(str(paths["ref_gff"]),
                                               str(paths["ref_fa"]), log=log)
    summaries, validity = {}, {}
    for tool in TOOLS:
        gff = gffs.get(tool)
        if gff is None or not Path(gff).exists() or Path(gff).stat().st_size == 0:
            log(f"  [{bid}/{mode}] skip {tool}: missing/empty gff")
            continue
        summaries[tool] = evaluator.evaluate_tool(
            tool, str(gff), str(paths["tgt_fa"]), ref, man, eval_dir,
            profs.get(tool), log=log, ref_index=ref_index, threads=threads_eval)
        validity[tool] = vc.validate_gff(gff, log=log)

    present = [t for t in TOOLS if t in summaries]
    comp = {t: summaries[t]["completeness_coding"] for t in present}
    comp_feat = {t: summaries[t]["completeness_feature_total"] for t in present}
    mean_pi = {t: summaries[t]["protein_identity"]["mean"] for t in present}
    median_pi = {t: summaries[t]["protein_identity"]["median"] for t in present}
    pct_id = {t: summaries[t]["protein_identity"]["pct_identical"] for t in present}
    nrec = {t: summaries[t]["n_recovered_coding"] for t in present}
    nrec_any = {t: summaries[t]["n_recovered_any"] for t in present}
    wall = {t: (profs.get(t) or {}).get("wall_clock_seconds") for t in present}
    rss = {t: (profs.get(t) or {}).get("peak_rss_mb") for t in present}

    dev, sta = "lifton_devel", "lifton_stable"
    rec = {
        "benchmark": bid, "mode": mode, "key": f"{mode}:{bid}",
        "species": bench["species"], "cross_species": bench["cross_species"],
        "dimension": bench.get("dimension"),
        "divergence_class": ("same_species" if not bench["cross_species"]
                             else bench.get("divergence_class", "cross_species")),
        "annotation_database": bench.get("annotation_database", "RefSeq"),
        "n_reference_coding": summaries[present[0]]["n_reference_coding"] if present else None,
        "n_reference_total": summaries[present[0]]["n_reference_total"] if present else None,
        "tools": present,
        "completeness_coding": comp, "completeness_feature_total": comp_feat,
        "mean_pi": mean_pi, "median_pi": median_pi, "pct_identical": pct_id,
        "n_recovered_coding": nrec, "n_recovered_any": nrec_any,
        "wall_s": wall, "peak_rss_mb": rss, "validity": validity,
        "feature_census": {t: summaries[t].get("completeness_by_type", {}) for t in present},
        "devel_vs_stable": {
            "meanpi": _sub(mean_pi.get(dev), mean_pi.get(sta)),
            "completeness": _sub(comp.get(dev), comp.get(sta)),
            "n_recovered": _sub(nrec.get(dev), nrec.get(sta)),
        },
        "devel_vs_best_baseline": {
            "meanpi": _sub(mean_pi.get(dev), _best_baseline(mean_pi)),
            "completeness": _sub(comp.get(dev), _best_baseline(comp)),
        },
    }
    # Audit finding #1: the joint recall-vs-identity block (apples-to-apples
    # common-set PI + covPI + recall@PI), computed from the per-transcript TSVs
    # the evaluator just wrote into eval_dir. The headline set-mean metrics above
    # are over each tool's OWN recovered set; this is the fair head-to-head on
    # the transcripts both tools recover. See joint_metrics.py.
    try:
        from . import joint_metrics as _jm
        rec["joint"] = _jm.compute_joint_metrics(str(eval_dir), rec["n_reference_coding"])
    except Exception as _e:  # never let metric enrichment break a cell
        log(f"  [{bid}/{mode}] joint-metric computation skipped: {_e}")
    _save(rec)
    log(f"  [{bid}/{mode}] completeness_coding: " +
        " ".join(f"{t}={comp.get(t)}" for t in present))
    log(f"  [{bid}/{mode}] mean_pi: " +
        " ".join(f"{t}={mean_pi.get(t)}" for t in present))
    log(f"      devel vs v1.0.8: Δmeanpi={rec['devel_vs_stable']['meanpi']} "
        f"Δcompleteness={rec['devel_vs_stable']['completeness']}; "
        f"devel vs best(liftoff,miniprot): Δmeanpi={rec['devel_vs_best_baseline']['meanpi']}")
    return rec


def _save(rec):
    db = json.loads(RESULTS.read_text()) if RESULTS.exists() else {}
    db[rec["key"]] = rec
    RESULTS.write_text(json.dumps(db, indent=2))


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--full", action="store_true",
                    help="full-genome 4-way (reuse lifton_vc full outputs / run fresh)")
    ap.add_argument("ids", nargs="+")
    a = ap.parse_args(argv)
    # one provenance gate up front (stable SHA e503643d + no --native; devel HEAD)
    vc.provenance_gate(["stable", "devel"], log=print)
    for bid in a.ids:
        print(f"\n=== {bid} [{'full' if a.full else 'subset'}] : liftoff | miniprot | "
              f"LiftOn v1.0.8 | LiftOn devel ===", flush=True)
        try:
            (run_full if a.full else run_subset)(bid, log=print)
        except Exception as e:
            import traceback
            print(f"!! {bid} FAILED: {e}\n{traceback.format_exc()}", flush=True)
    print(f"\nDONE -> {RESULTS}", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
