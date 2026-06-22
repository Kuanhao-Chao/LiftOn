#!/usr/bin/env python
"""Build just the A/B inputs for one benchmark: subset + Liftoff + miniprot.

The standalone A/B drivers (``legacy_merge_ab.py``, ``frame_gate_xspecies.py``)
consume cached ``tools/liftoff/liftoff.gff3`` + ``tools/miniprot/miniprot.gff3``
(genes mode, transcript-space miniprot). This helper regenerates exactly those
for a single benchmark id WITHOUT running the full ``run_compare`` pipeline (no
lifton variants, no eval/report), which is what you want after changing a
benchmark's ``ref_chrom`` in ``benchmarks.json`` (e.g. pinning mouse_to_rat to
the smaller chr18 for a faster mammalian A/B).

Usage (repo root, lifton_devel env):
    python -m benchmarks.compare.build_inputs mouse_to_rat [--threads N]
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from . import subset_builder, tool_runners

HERE = Path(__file__).resolve().parent
WORK = HERE / "work"
REG = json.loads((HERE / "benchmarks.json").read_text())


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("benchmark", help="benchmark id from benchmarks.json")
    ap.add_argument("-t", "--threads", type=int, default=8)
    ap.add_argument("--no-force", action="store_true",
                    help="reuse cached stages (.done sentinels) instead of rebuilding")
    args = ap.parse_args(argv)

    bench = next((b for b in REG["benchmarks"] if b["id"] == args.benchmark), None)
    if bench is None:
        ap.error(f"unknown benchmark id {args.benchmark!r}")
    tools = REG["tools"]
    work = WORK / args.benchmark
    force = not args.no_force

    print(f"=== building A/B inputs for {args.benchmark} "
          f"(ref_chrom={bench.get('ref_chrom')}, force={force}) ===", flush=True)

    manifest = subset_builder.build_subset(bench, work, tools, threads=args.threads,
                                           force=force, log=print)
    print(f"  [subset] ref_chrom={manifest.get('ref_chrom')} "
          f"tgt_chroms={manifest.get('tgt_chroms')} "
          f"n_proteins={manifest.get('n_subset_proteins')}", flush=True)

    tool_runners.run_liftoff(bench, manifest, work, tools, threads=args.threads,
                             force=force, log=print, mode="genes")
    tool_runners.run_miniprot(bench, manifest, work, tools, threads=args.threads,
                              force=force, log=print)

    lf = work / "tools" / "liftoff" / "liftoff.gff3"
    mp = work / "tools" / "miniprot" / "miniprot.gff3"
    ok = lf.exists() and lf.stat().st_size > 0 and mp.exists() and mp.stat().st_size > 0
    print(f"\n{'DONE' if ok else 'INCOMPLETE'}: liftoff={lf} ({lf.stat().st_size if lf.exists() else 0} B), "
          f"miniprot={mp} ({mp.stat().st_size if mp.exists() else 0} B)", flush=True)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
