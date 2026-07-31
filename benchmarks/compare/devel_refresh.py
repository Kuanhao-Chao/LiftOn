#!/usr/bin/env python
"""devel_refresh.py — re-run ONLY `lifton_devel` on a full-genome pair with the
CURRENT working-tree code (this session's candidate-3 fix + the 3 already-
committed correctness fixes since the report's cited SHA) and re-score it,
producing a per-pair JSON patch for `fourway_results.json`.

Liftoff, miniprot, and `lifton_stable` (v1.0.8) are UNCHANGED by any of this
session's work, so their cached fourway_results.json fields stay as-is; only
`lifton_devel`'s fields (mean_pi, completeness, validity, n_recovered_*,
feature_census, joint metrics, devel_vs_best_baseline, devel_vs_stable) need
recomputing.

Protocol fidelity: mirrors `version_compare._build_argv`'s devel-fast full-genome
invocation EXACTLY (-t 8 -copies -ad <db> --locus-pipeline --no-miniprot-rescue,
target then reference) with ONE deliberate addition: cached -L/-M (Liftoff/
miniprot GFF3, unaffected by this session's changes) and -T/-P (extracted
transcripts.fa/proteins.fa, unaffected — Step 3 depends only on the reference)
from the existing `_fourway_full/devel/lifton_output/` directory, bypassing a
fresh Step-3 extraction. Path resolution reuses `fourway_compare._full_paths`
directly so cross-species pairs and the three special-cased same-species pairs
(bee/rice/arabidopsis, resolved via `version_compare.FULL_INPUTS`) are handled
identically to the canonical harness.

Run (repo root, lifton_devel env):
  python -m benchmarks.compare.devel_refresh <pair_id> [<pair_id> ...]
Writes benchmarks/compare/devel_refresh/<pair_id>.json
"""
from __future__ import annotations
import datetime as _dt
import json, os, subprocess, sys
from pathlib import Path

from . import evaluator, fourway_compare as fc, version_compare as vc
from .profiling import run_profiled

HERE = Path(__file__).resolve().parent
WORK = HERE / "work"
OUT_DIR = HERE / "devel_refresh"
OUT_DIR_DEFAULT = HERE / "devel_refresh_default"

ALL_PAIRS = [
    "arabidopsis", "bee", "rice", "t1_maize_b73_to_mo17", "t1_tomato_microtom_to_heinz",
    "drosophila", "t2_human_to_gorilla", "t2_mouse_to_caroli", "t2_tomato_to_potato",
    "t3_dog_to_cat", "t3_human_to_macaque", "t3_human_to_marmoset",
    "arabidopsis_to_rice", "human_to_zebrafish", "t4_human_to_chicken",
    "t4_human_to_xenopus", "t4_drosophila_to_bee",
]


def _cached_inputs(bid):
    """The existing _fourway_full/devel/lifton_output/ cache for this pair —
    Liftoff/miniprot GFF3 (unaffected by this session) + extracted transcripts/
    proteins FASTA (Step 3 depends only on the reference, also unaffected)."""
    lo_dir = WORK / bid / "_fourway_full" / "devel" / "lifton_output"
    L = lo_dir / "liftoff" / "liftoff.gff3"
    M = lo_dir / "miniprot" / "miniprot.gff3"
    T = lo_dir / "intermediate_files" / "transcripts.fa"
    P = lo_dir / "intermediate_files" / "proteins.fa"
    for f in (L, M, T, P):
        if not f.exists() or f.stat().st_size == 0:
            raise RuntimeError(f"{bid}: missing cached input {f}")
    return L, M, T, P


def campaign_provenance():
    """Host and build identity for a measurement, so campaigns can be compared.

    The previously frozen `fourway_results.json` recorded none of this, so when
    the report's wall-clock and peak-memory tables needed refreshing there was no
    way to tell whether the old numbers came from this machine or another one --
    which forced re-measuring *both* arms rather than just the changed one.
    Recording it costs nothing and removes that ambiguity next time.
    """
    import platform
    import socket

    def _git(worktree):
        try:
            return subprocess.run(
                ["git", "-C", str(worktree), "rev-parse", "HEAD"],
                capture_output=True, text=True, timeout=15,
            ).stdout.strip()[:12] or None
        except Exception:
            return None

    mem_gib = None
    try:
        mem_gib = round(os.sysconf("SC_PAGE_SIZE") * os.sysconf("SC_PHYS_PAGES")
                        / (1024 ** 3), 1)
    except Exception:
        pass

    return {
        "host": socket.gethostname(),
        "cpu_count": os.cpu_count(),
        "total_mem_gib": mem_gib,
        "platform": platform.platform(),
        "python": platform.python_version(),
        "sha": {
            "devel": _git(vc.VERSIONS["devel"]["worktree"]),
            "stable": _git(vc.VERSIONS["stable"]["worktree"]),
        },
        "measured_at": _dt.datetime.now(_dt.timezone.utc).isoformat(timespec="seconds"),
    }


def run_devel_refresh(bid, threads=8, threads_eval=8, log=print, rescue_on=False):
    bench = fc._bench(bid)
    anndb = bench.get("annotation_database", "RefSeq")
    paths = fc._full_paths(bid)
    for k, p in paths.items():
        if not Path(p).exists():
            raise RuntimeError(f"full input missing: {k}={p}")
    L, M, T, P = _cached_inputs(bid)

    # The shipped-default arm gets its own work + output trees so it can never
    # overwrite the pinned arm that feeds the frozen four-way column.
    suffix = "_default" if rescue_on else ""
    out_dir = WORK / bid / f"_devel_refresh{suffix}"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_gff = out_dir / f"devel_refresh{suffix}.gff3"

    # Mirror version_compare._build_argv's devel-fast full-genome protocol
    # EXACTLY (order + flags), adding only -L/-M/-T/-P for the cached bypass.
    argv = [
        vc.VERSIONS["devel"]["bin"], "-t", str(threads), "-copies",
        "-ad", anndb, "-g", str(paths["ref_gff"]),
        "-L", str(L), "-M", str(M), "-T", str(T), "-P", str(P),
        "--locus-pipeline",
    ]
    if not rescue_on:
        # Pinned arm: matches the frozen four-way protocol, comparable to v1.0.8.
        # The rescue-on arm omits this so it measures the SHIPPED default.
        argv += ["--no-miniprot-rescue"]
    argv += [
        "-o", str(out_gff), str(paths["tgt_fa"]), str(paths["ref_fa"]),
    ]
    log(f"[{bid}] invoking: {' '.join(argv)}")
    pr = run_profiled(argv, label=f"devel_refresh_{bid}", log_dir=out_dir / "logs",
                      env=vc._env(), cwd=out_dir, log=log)
    if pr.exit_code != 0 or not out_gff.exists() or out_gff.stat().st_size == 0:
        raise RuntimeError(f"{bid}: devel refresh lift failed (exit {pr.exit_code}); "
                          f"see {pr.stderr_path}")

    # Re-score exactly like _score_and_record does for a single tool.
    ref, ref_index = evaluator.build_reference(str(paths["ref_gff"]), str(paths["ref_fa"]), log=log)
    man = {"id": bid, "species": bench["species"], "cross_species": bench["cross_species"],
          "miniprot_target_space": "transcript", "protein_acc_to_mrna": {},
          "paths": {k: str(v) for k, v in paths.items()}}
    eval_dir = out_dir / "eval"; eval_dir.mkdir(parents=True, exist_ok=True)
    summary = evaluator.evaluate_tool(
        "lifton_devel", str(out_gff), str(paths["tgt_fa"]), ref, man, eval_dir,
        dataclasses_asdict(pr), log=log, ref_index=ref_index, threads=threads_eval)
    validity = vc.validate_gff(out_gff, log=log)

    result = {"bid": bid, "out_gff": str(out_gff), "eval_dir": str(eval_dir),
             "summary": summary, "validity": validity,
             "arm": "shipped_default" if rescue_on else "pinned_no_rescue",
             "argv": argv, "provenance": campaign_provenance()}
    out_root = OUT_DIR_DEFAULT if rescue_on else OUT_DIR
    out_root.mkdir(parents=True, exist_ok=True)
    (out_root / f"{bid}.json").write_text(json.dumps(result, indent=2, default=str))
    log(f"[{bid}] DONE: completeness_coding={summary.get('completeness_coding')} "
        f"mean_pi={summary.get('mean_pi')} validity_errors={validity.get('n_errors')}")
    return result


def dataclasses_asdict(pr):
    import dataclasses
    return dataclasses.asdict(pr) if dataclasses.is_dataclass(pr) else (pr or {})


def main(argv=None):
    import argparse
    ap = argparse.ArgumentParser(description=__doc__ and __doc__.splitlines()[0])
    ap.add_argument("pairs", nargs="*", help="pair ids (default: all 17)")
    ap.add_argument("--rescue-on", action="store_true",
                    help="measure the SHIPPED default (miniprot-only rescue ON). "
                         "Writes to devel_refresh_default/ so the pinned arm that "
                         "feeds the frozen four-way column is never overwritten.")
    ap.add_argument("-t", "--threads", type=int, default=8)
    ap.add_argument("--threads-eval", type=int, default=8)
    args = ap.parse_args(argv if argv is not None else sys.argv[1:])

    pairs = args.pairs or ALL_PAIRS
    out_root = OUT_DIR_DEFAULT if args.rescue_on else OUT_DIR
    failures = []
    for bid in pairs:
        try:
            run_devel_refresh(bid, threads=args.threads,
                              threads_eval=args.threads_eval,
                              rescue_on=args.rescue_on)
        except Exception as e:
            print(f"[{bid}] FAILED: {e}", file=sys.stderr)
            out_root.mkdir(parents=True, exist_ok=True)
            (out_root / f"{bid}.FAILED.txt").write_text(str(e))
            failures.append(bid)
    if failures:
        print(f"FAILED pairs ({len(failures)}): {failures}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
