#!/usr/bin/env python
"""Fold the per-cell parallel-safe ``_full_runs/<id>.json`` files into the
shared ``fourway_results.json`` (the documented post-step from
``fourway_compare.py``: concurrent ``--full`` runs each write their OWN file via
``FOURWAY_RESULTS_JSON`` to avoid racing the unlocked read-modify-write; this
merges them afterward).

Single-writer by construction (run once, serially). Skips degenerate cells whose
``tools`` list is empty (e.g. the soybean GCA that resolves to metadata-only).
Idempotent: re-running overwrites each ``full:<id>`` key with the per-cell file's
record. The subset cells already in ``fourway_results.json`` are untouched.

Usage:
    python -m benchmarks.compare.merge_full_runs            # dry-run (report)
    python -m benchmarks.compare.merge_full_runs --apply    # write
"""
import argparse, json, sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
RUNS = HERE / "_full_runs"
RESULTS = HERE / "fourway_results.json"


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--apply", action="store_true",
                    help="write the merge (default: dry-run report only)")
    a = ap.parse_args(argv)

    db = json.loads(RESULTS.read_text()) if RESULTS.exists() else {}
    before_full = sorted(k for k in db if k.startswith("full:"))

    add, skip = [], []
    for p in sorted(RUNS.glob("*.json")):
        cell = json.loads(p.read_text())
        for key, rec in cell.items():
            tools = rec.get("tools") or []
            if not tools:
                skip.append((key, "empty tools (degenerate cell)"))
                continue
            status = "update" if key in db else "add"
            add.append((key, status, tools))
            if a.apply:
                db[key] = rec

    print(f"=== merge_full_runs ({'APPLY' if a.apply else 'DRY-RUN'}) ===")
    print(f"fourway_results.json full: cells BEFORE ({len(before_full)}): "
          f"{[k.split(':',1)[1] for k in before_full]}")
    for key, status, tools in add:
        print(f"  [{status:6}] {key}  tools={tools}")
    for key, why in skip:
        print(f"  [SKIP  ] {key}  ({why})")

    if a.apply:
        RESULTS.write_text(json.dumps(db, indent=2))
        after = sorted(k for k in db if k.startswith("full:"))
        print(f"WROTE {RESULTS}")
        print(f"full: cells AFTER ({len(after)}): "
              f"{[k.split(':',1)[1] for k in after]}")
    else:
        print("(dry-run — pass --apply to write)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
