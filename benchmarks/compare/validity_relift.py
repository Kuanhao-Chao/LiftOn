#!/usr/bin/env python
"""Track-7 validity-focused re-lift for the LiftOn v1.0.9 pre-release refresh.

Re-lifts ONLY the ``lifton_devel`` column for the 17 full-genome cells with the
CURRENT devel code (Iteration-24 containment normalization default-ON), then
re-validates every available tool GFF (liftoff / miniprot / lifton_stable / the
new lifton_devel) with the CURRENT ``gff3-validate``, and writes back ONLY the
``validity[*]`` blocks of the ``full:*`` cells in ``fourway_results.json``.
Accuracy / completeness / performance are left exactly as committed (provably
unchanged by Iter-24, which never touches CDS coordinates).

Fast path: when cached ``-L``/``-M`` (a prior run's ``liftoff.gff3`` /
``miniprot.gff3``) exist, the devel re-lift is merge-only (skips the external
aligners) and isolates the Iter-24 change; otherwise a fresh full devel lift runs
(arabidopsis/bee/rice, which have no cached aligner outputs).

Concurrency: a process pool runs several genomes at once; only the PARENT writes
the JSON (workers return validity dicts), so there is no race. ``--threads N`` is
byte-identical to ``-t1`` (the 24-cell matrix), so concurrency does not change any
validity number. Run from the repo root:

    PYTHONNOUSERSITE=1 PYTHONHASHSEED=0 python -m benchmarks.compare.validity_relift \
        --workers 4 --threads 8
"""
import argparse
import json
import shutil
import sys
import time
import traceback
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

sys.path.insert(0, "/ccb/salz3/kh.chao/LiftOn")
from benchmarks.compare import fourway_compare as fw  # noqa: E402
from benchmarks.compare import version_compare as vc  # noqa: E402

ROOT = Path("/ccb/salz3/kh.chao/LiftOn")
WORK = ROOT / "benchmarks/compare/work"
RESULTS = ROOT / "benchmarks/compare/fourway_results.json"

# Roughly small -> large so the quick wins land (and surface any bug) first.
FULL_IDS = [
    "t4_drosophila_to_bee", "drosophila", "arabidopsis_to_rice",
    "human_to_zebrafish", "t4_human_to_xenopus", "t4_human_to_chicken",
    "t2_tomato_to_potato", "t1_tomato_microtom_to_heinz", "t1_maize_b73_to_mo17",
    "t3_dog_to_cat", "t2_mouse_to_caroli", "t3_human_to_macaque",
    "t3_human_to_marmoset", "t2_human_to_gorilla",
    "bee", "arabidopsis", "rice",
]


def _cached_LM(bid):
    """Cached liftoff/miniprot GFFs from a prior full run (prefer devel's own)."""
    ff = WORK / bid / "_fourway_full"
    for sub in ("devel", "stable"):
        L = ff / sub / "lifton_output" / "liftoff" / "liftoff.gff3"
        M = ff / sub / "lifton_output" / "miniprot" / "miniprot.gff3"
        if L.exists() and L.stat().st_size > 0 and M.exists() and M.stat().st_size > 0:
            return L, M
    return None, None


def _tool_gffs(bid, new_devel):
    """Map the 4 tool columns to GFF paths to validate (current validator)."""
    ff = WORK / bid / "_fourway_full"
    g = {}
    for col in ("liftoff", "miniprot"):
        for c in (new_devel.parent / "lifton_output" / col / f"{col}.gff3",
                  ff / "devel" / "lifton_output" / col / f"{col}.gff3",
                  ff / "stable" / "lifton_output" / col / f"{col}.gff3"):
            if c.exists() and c.stat().st_size > 0:
                g[col] = c
                break
    sg = ff / "stable" / "stable.gff3"
    if sg.exists() and sg.stat().st_size > 0:
        g["lifton_stable"] = sg
    g["lifton_devel"] = new_devel
    return g


def process_one(bid, threads):
    t0 = time.time()
    try:
        paths = fw._full_paths(bid)
        for k, p in paths.items():
            if not Path(p).exists():
                return {"bid": bid, "error": f"missing input {k}={p}"}
        L, M = _cached_LM(bid)
        mode = "merge-only" if L else "fresh"
        statedir = WORK / bid / "_fourway_full" / "devel_v109"
        if statedir.exists():
            shutil.rmtree(statedir)
        out_gff, _pr = vc.run_lift("devel", paths, "RefSeq", threads, L, M,
                                   statedir, devel_fast=True, log=lambda *_: None)
        validity = {col: vc.validate_gff(str(gff))
                    for col, gff in _tool_gffs(bid, out_gff).items()}
        return {"bid": bid, "mode": mode, "secs": round(time.time() - t0),
                "validity": validity, "out": str(out_gff),
                "size_mb": out_gff.stat().st_size // 1048576}
    except Exception as e:  # noqa: BLE001
        return {"bid": bid, "error": str(e), "tb": traceback.format_exc(),
                "secs": round(time.time() - t0)}


RELIFT_DIR = ROOT / "benchmarks/compare/_validity_relift"


def _log(m):
    print(f"[{time.strftime('%H:%M:%S')}] {m}", flush=True)


def _single(bid, threads, out_dir):
    """Process ONE genome and write its result to ``<out_dir>/<bid>.result.json``
    plus a ``<bid>.done`` sentinel. Lets each genome run in its OWN tmux session
    in parallel with no shared-JSON race (the parallel runs only touch per-genome
    files; a later ``--merge`` folds them into fourway_results.json)."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    _log(f"single {bid}: starting (threads={threads})")
    r = process_one(bid, threads)
    (out_dir / f"{bid}.result.json").write_text(json.dumps(r, indent=1))
    (out_dir / f"{bid}.done").write_text(time.strftime("%Y-%m-%d %H:%M:%S"))
    if r.get("error"):
        _log(f"single {bid}: ERROR ({r.get('secs', '?')}s): {r['error']}")
        return 1
    nerr = {k: v.get("n_errors") for k, v in r["validity"].items()}
    _log(f"single {bid}: DONE [{r['mode']}, {r['secs']}s, {r['size_mb']}M] errors={nerr}")
    return 0


def _merge(out_dir, base_ref="HEAD"):
    """Fold all per-genome ``<bid>.result.json`` files into fourway_results.json,
    updating ONLY the ``validity[*]`` blocks of the ``full:<bid>`` cells.

    Starts from the committed (``base_ref``) JSON so a re-merge is idempotent, and
    updates validity PER TOOL with two guards: (1) never overwrite a committed
    value that was null / had ``n_errors=None`` -- that marks a tool that crashed
    or produced no scorable output (e.g. v1.0.8 on maize/tomato), so a leftover
    partial GFF must not "resurrect" it; (2) never DELETE a committed tool value
    just because a fresh re-lift lacked that tool's cached GFF (e.g. arabidopsis/
    bee/rice have no cached stable.gff3 -> keep the committed lifton_stable). This
    keeps Table 8 consistent with the report's crash markers and accuracy table."""
    import subprocess
    out_dir = Path(out_dir)
    base = subprocess.run(
        ["git", "-C", str(ROOT), "show",
         f"{base_ref}:benchmarks/compare/fourway_results.json"],
        capture_output=True, text=True)
    res = json.loads(base.stdout)
    updated, skipped = 0, 0
    for f in sorted(out_dir.glob("*.result.json")):
        r = json.loads(f.read_text())
        bid = r.get("bid")
        key = f"full:{bid}"
        if r.get("error"):
            _log(f"MERGE skip {bid}: error result ({r.get('error')})")
            skipped += 1
            continue
        if key not in res:
            _log(f"MERGE skip {bid}: no {key} cell")
            skipped += 1
            continue
        old = res[key].setdefault("validity", {})
        changes = []
        for tool, nv in r["validity"].items():
            ov = old.get(tool)
            if ov is None or (isinstance(ov, dict) and ov.get("n_errors") is None):
                continue  # preserve crashed/absent committed tool (don't resurrect)
            oe = ov.get("n_errors") if isinstance(ov, dict) else None
            old[tool] = nv
            changes.append(f"{tool}:{oe}->{nv.get('n_errors')}")
        res[key]["_v109_validity_refresh"] = {
            "mode": r["mode"], "secs": r["secs"], "devel_gff": r["out"]}
        updated += 1
        _log(f"MERGE {bid}: " + ", ".join(changes))
    tmp = RESULTS.with_suffix(".json.tmp")
    tmp.write_text(json.dumps(res, indent=2))
    tmp.replace(RESULTS)
    _log(f"MERGE DONE: {updated} full cells updated, {skipped} skipped, "
         f"base={base_ref}, in {RESULTS.name}")


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--ids", nargs="*", default=FULL_IDS)
    ap.add_argument("--single", default=None,
                    help="process ONE genome and write a per-genome result file")
    ap.add_argument("--out-dir", default=str(RELIFT_DIR))
    ap.add_argument("--merge", action="store_true",
                    help="fold per-genome result files into fourway_results.json")
    a = ap.parse_args(argv)

    if a.single:
        return _single(a.single, a.threads, a.out_dir)
    if a.merge:
        _merge(a.out_dir)
        return 0

    # default: in-process pool over a.ids; the parent writes the JSON (no race).
    res = json.load(open(RESULTS))
    _log(f"validity re-lift: {len(a.ids)} genomes, workers={a.workers}, threads={a.threads}")
    done = 0
    with ProcessPoolExecutor(max_workers=a.workers) as ex:
        futs = {ex.submit(process_one, bid, a.threads): bid for bid in a.ids}
        for fut in as_completed(futs):
            r = fut.result()
            bid = r["bid"]
            key = f"full:{bid}"
            if r.get("error"):
                _log(f"ERROR {bid} ({r.get('secs', '?')}s): {r['error']}")
                continue
            if key not in res:
                _log(f"WARN {bid}: no {key} cell; skip JSON update")
                continue
            old = {k: v.get("n_errors") for k, v in res[key].get("validity", {}).items()}
            res[key]["validity"] = r["validity"]
            res[key]["_v109_validity_refresh"] = {
                "mode": r["mode"], "secs": r["secs"], "devel_gff": r["out"]}
            tmp = RESULTS.with_suffix(".json.tmp")
            tmp.write_text(json.dumps(res, indent=2))
            tmp.replace(RESULTS)
            new = {k: v.get("n_errors") for k, v in r["validity"].items()}
            done += 1
            _log(f"DONE {bid} [{r['mode']}, {r['secs']}s] ({done}/{len(a.ids)}) "
                 + ", ".join(f"{k}:{old.get(k)}->{new.get(k)}" for k in new))
    _log(f"ALL DONE ({done}/{len(a.ids)} updated)")
    return 0


if __name__ == "__main__":
    sys.exit(main() or 0)
