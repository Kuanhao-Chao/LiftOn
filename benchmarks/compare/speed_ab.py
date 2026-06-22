#!/usr/bin/env python
"""Serial-vs-threaded speed A/B for LiftOn's byte-identical fast paths.

For each fast-subset dataset, run LiftOn TWICE on the same cached -L/-M inputs
(genes mode, NON-optimize — byte-identity only holds for the default path):

  * serial   : -t 1                              (no --native / --locus-pipeline)
  * threaded : -t 8 --native --locus-pipeline    (Phase-17b/17c materialised path)

Captures wall-clock + peak RSS via the audited /usr/bin/time profiler, ASSERTS
the two outputs are byte-identical (the harness-side echo of the 24-cell matrix
gate), and writes speed_ab.json + speed_ab.md. This confirms the Phase-17 speed
work is both faster AND output-preserving.

Usage (repo root, lifton_devel env):
    python -m benchmarks.compare.speed_ab [IDS...]   # default = fast subset
"""
from __future__ import annotations

import filecmp
import json
import sys
from pathlib import Path

from .profiling import run_profiled
from .tool_runners import _clean_input_dbs, _compose_env

HERE = Path(__file__).resolve().parent
WORK = HERE / "work"
REG = json.loads((HERE / "benchmarks.json").read_text())
TOOLS = REG["tools"]

FAST_DEFAULT = ["human_mane", "drosophila", "arabidopsis", "rice"]
SERIAL_THREADS = 1
THREADED_THREADS = 8


def _ann_db(bid: str) -> str:
    for b in REG["benchmarks"]:
        if b["id"] == bid:
            return b.get("annotation_database", "RefSeq")
    return "RefSeq"


def _run(bid: str, p: dict, threaded: bool):
    W = WORK / bid
    ab = W / "_speedab"
    log_dir = ab / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    liftoff_gff = W / "tools" / "liftoff" / "liftoff.gff3"
    miniprot_gff = W / "tools" / "miniprot" / "miniprot.gff3"
    tag = "threaded" if threaded else "serial"
    out = ab / f"{tag}.gff3"
    threads = THREADED_THREADS if threaded else SERIAL_THREADS
    # fresh input DBs (guards against a stale miniprot db suppressing the merge)
    _clean_input_dbs(p["ref_gff"], liftoff_gff, miniprot_gff)
    argv = [TOOLS["lifton_bin"], "-t", str(threads), "-copies",
            "-ad", _ann_db(bid),
            "-g", p["ref_gff"],
            "-L", str(liftoff_gff),
            "-M", str(miniprot_gff),
            "-o", str(out)]
    if threaded:
        argv += ["--native", "--locus-pipeline"]
    argv += [p["tgt_fa"], p["ref_fa"]]
    pr = run_profiled(argv, label=f"speed_{tag}", log_dir=log_dir,
                      env=_compose_env(TOOLS), log=print)
    if pr.exit_code != 0:
        raise RuntimeError(f"{bid}:{tag} lifton failed (exit {pr.exit_code}); see {pr.stderr_path}")
    if not out.exists() or out.stat().st_size == 0:
        raise RuntimeError(f"{bid}:{tag} produced no/empty output {out}")
    return out, pr


def _write_md(results, path: Path):
    lines = [
        "## Speed: serial vs threaded (byte-identical)\n",
        "LiftOn run twice on the same cached `-L`/`-M` inputs (genes mode, "
        f"non-`--optimize`): **serial** `-t {SERIAL_THREADS}` (no fast path) vs "
        f"**threaded** `-t {THREADED_THREADS} --native --locus-pipeline` "
        "(Phase-17b/17c materialised). `Byte-identical` MUST be ✓ — the "
        "harness-side echo of the 24-cell matrix gate.\n",
        "| Dataset | Serial wall (s) | Threaded wall (s) | Speedup | Serial RSS (MB) "
        "| Threaded RSS (MB) | Byte-identical |",
        "|---|---|---|---|---|---|---|",
    ]
    for r in results:
        lines.append("| {} | {} | {} | {} | {} | {} | {} |".format(
            r["benchmark"], r["serial_wall_s"], r["threaded_wall_s"],
            f"{r['speedup']}×" if r["speedup"] else "—",
            r["serial_rss_mb"], r["threaded_rss_mb"],
            "✓" if r["byte_identical"] else "✗ DIFFERS"))
    path.write_text("\n".join(lines) + "\n")


def main(argv=None):
    ids = (argv if argv is not None else sys.argv[1:]) or FAST_DEFAULT
    results = []
    for bid in ids:
        print(f"=== {bid} ===", flush=True)
        man = json.loads((WORK / bid / "subset" / "subset.manifest.json").read_text())
        p = man["paths"]
        s_out, s_pr = _run(bid, p, threaded=False)
        t_out, t_pr = _run(bid, p, threaded=True)
        identical = filecmp.cmp(str(s_out), str(t_out), shallow=False)
        speedup = (s_pr.wall_clock_seconds / t_pr.wall_clock_seconds
                   if t_pr.wall_clock_seconds else None)
        rec = {
            "benchmark": bid,
            "serial_threads": SERIAL_THREADS,
            "threaded_threads": THREADED_THREADS,
            "serial_wall_s": round(s_pr.wall_clock_seconds, 2),
            "threaded_wall_s": round(t_pr.wall_clock_seconds, 2),
            "speedup": round(speedup, 2) if speedup else None,
            "serial_rss_mb": round(s_pr.peak_rss_mb, 0),
            "threaded_rss_mb": round(t_pr.peak_rss_mb, 0),
            "byte_identical": identical,
        }
        results.append(rec)
        print(f"  [{bid}] serial={rec['serial_wall_s']}s threaded={rec['threaded_wall_s']}s "
              f"speedup={rec['speedup']}x byte_identical={identical}", flush=True)
        if not identical:
            print(f"  [{bid}] !! WARNING: serial vs threaded output DIFFERS — the fast "
                  f"path is no longer byte-frozen!", flush=True)

    (HERE / "speed_ab.json").write_text(json.dumps(results, indent=2))
    _write_md(results, HERE / "speed_ab.md")
    n_ident = sum(1 for r in results if r["byte_identical"])
    print(f"\nDONE: {n_ident}/{len(results)} byte-identical; wrote speed_ab.json + speed_ab.md",
          flush=True)
    return 0 if n_ident == len(results) else 1


if __name__ == "__main__":
    sys.exit(main())
