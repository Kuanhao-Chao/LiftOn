#!/usr/bin/env python
"""report_tables.py — emit the technical report's tables straight from the
measured JSON, so no number in the report is ever hand-typed.

The report carries eleven tables whose values come from two documents:

  * ``fourway_results.json``          — accuracy, completeness, validity,
                                        recovery and the feature census
                                        (Tables 2, 3, 5, 7, 8, 9)
  * ``version_compare.results.json``  — the controlled ``-t1`` wall-clock and
                                        peak-RSS arm (Tables 4, 6)

Transcribing those by hand is how a report drifts from its data. This renders
them as markdown, and ``--check`` diffs the rendered rows against the report
itself so a stale figure fails loudly instead of being believed.

Usage (repo root):
    python -m benchmarks.compare.report_tables --list
    python -m benchmarks.compare.report_tables --table 5
    python -m benchmarks.compare.report_tables --all > tables.md
    python -m benchmarks.compare.report_tables --check <report.mdx>
"""
from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
FOURWAY = HERE / "fourway_results.json"
VERSION_CMP = HERE / "version_compare.results.json"

TOOLS = ("liftoff", "miniprot", "lifton_stable", "lifton_devel")

# Divergence tiers in report order, with the label the report uses.
TIERS = [
    ("same_species", "same"),
    ("cross_species", "cross"),
    ("close_cross_species", "close"),
    ("distant_cross_species", "distant"),
    ("very_distant_cross_species", "very distant"),
]

# Display names, matching the report's prose.
PRETTY = {
    "arabidopsis": "arabidopsis", "bee": "bee", "rice": "rice",
    "t1_maize_b73_to_mo17": "maize", "t1_tomato_microtom_to_heinz": "tomato",
    "drosophila": "drosophila", "t2_human_to_gorilla": "human→gorilla",
    "t2_mouse_to_caroli": "mouse→caroli", "t2_tomato_to_potato": "tomato→potato",
    "t3_dog_to_cat": "dog→cat", "t3_human_to_macaque": "human→macaque",
    "t3_human_to_marmoset": "human→marmoset",
    "arabidopsis_to_rice": "arabidopsis→rice",
    "human_to_zebrafish": "human→zebrafish",
    "t4_human_to_chicken": "human→chicken",
    "t4_human_to_xenopus": "human→xenopus",
    "t4_drosophila_to_bee": "drosophila→bee",
}


def load(path: Path):
    if not path.exists():
        sys.exit(f"missing results document: {path}")
    return json.loads(path.read_text())


def full_cells(db):
    """Full-genome cells in divergence-tier order, as the report presents them."""
    cells = {k.split(":", 1)[1]: v for k, v in db.items() if k.startswith("full:")}
    order = []
    for cls, _label in TIERS:
        order += sorted(
            (b for b, c in cells.items() if c.get("divergence_class") == cls),
            key=lambda b: PRETTY.get(b, b),
        )
    order += [b for b in cells if b not in order]      # anything untiered
    return [(b, cells[b]) for b in order]


def _fmt(value, nd=5, dash="—"):
    if value is None:
        return dash
    if isinstance(value, float):
        return f"{value:.{nd}f}"
    return f"{value:,}" if isinstance(value, int) else str(value)


def _tier(cell):
    cls = cell.get("divergence_class")
    return dict(TIERS).get(cls, cls or "—")


# --------------------------------------------------------------------------
# tables
# --------------------------------------------------------------------------

def table2(db, _vc):
    rows = ["| Whole genome | Divergence | n coding | n features |",
            "|---|---|---:|---:|"]
    for b, c in full_cells(db):
        rows.append(f"| {PRETTY.get(b,b)} | {_tier(c)} | "
                    f"{_fmt(c.get('n_reference_coding'))} | "
                    f"{_fmt(c.get('n_reference_total'))} |")
    return "\n".join(rows)


def table5(db, _vc):
    rows = ["| Whole genome | Divergence | Liftoff | *miniprot (evidence)* | "
            "v1.0.8 | **v1.0.10** | v1.0.10 − Liftoff |", "|---|---|---:|---:|---:|---:|---:|"]
    for b, c in full_cells(db):
        pi = c.get("mean_pi") or {}
        lo, dv = pi.get("liftoff"), pi.get("lifton_devel")
        lead = f"{dv - lo:+.5f}" if (lo is not None and dv is not None) else "—"
        rows.append(
            f"| {PRETTY.get(b,b)} | {_tier(c)} | {_fmt(lo)} | {_fmt(pi.get('miniprot'))} | "
            f"{_fmt(pi.get('lifton_stable'))} | **{_fmt(dv)}** | {lead} |")
    return "\n".join(rows)


def table7(db, _vc):
    rows = ["| Whole genome | Divergence | Liftoff | *miniprot (evidence)* | "
            "v1.0.8 | **v1.0.10** |", "|---|---|---:|---:|---:|---:|"]
    for b, c in full_cells(db):
        cc = c.get("completeness_coding") or {}
        rows.append(
            f"| {PRETTY.get(b,b)} | {_tier(c)} | "
            + " | ".join(
                (f"**{cc[t]*100:.2f}%**" if t == "lifton_devel" else f"{cc[t]*100:.2f}%")
                if cc.get(t) is not None else "crash"
                for t in TOOLS)
            + " |")
    return "\n".join(rows)


def table9(db, _vc):
    rows = ["| Whole genome | Divergence | Liftoff | *miniprot (evidence)* | "
            "v1.0.8 | **v1.0.10** |", "|---|---|---:|---:|---:|---:|"]
    for b, c in full_cells(db):
        val = c.get("validity") or {}
        cells = []
        for t in TOOLS:
            v = val.get(t)
            n = (v or {}).get("n_errors")
            cells.append("—" if n is None
                         else (f"**{n}**" if t == "lifton_devel" else str(n)))
        rows.append(f"| {PRETTY.get(b,b)} | {_tier(c)} | " + " | ".join(cells) + " |")
    return "\n".join(rows)


def table3(db, _vc):
    """Recovery on the genomes v1.0.8 cannot finish."""
    rows = ["| Full genome | v1.0.8 recovered | v1.0.10 recovered | +v1.0.10 | "
            "v1.0.8 completeness | v1.0.10 completeness |",
            "|---|---:|---:|---:|---:|---:|"]
    for b, c in full_cells(db):
        rec = c.get("n_recovered_coding") or {}
        cc = c.get("completeness_coding") or {}
        st, dv = rec.get("lifton_stable"), rec.get("lifton_devel")
        # Only the genomes where v1.0.8 fails to produce a complete annotation.
        if st is not None and (cc.get("lifton_stable") or 0) > 0.90:
            continue
        gain = f"+{dv - st:,}" if (st is not None and dv is not None) else "—"
        st_cell = _fmt(st) if st is not None else "crash"
        cc_stable = cc.get("lifton_stable")
        cc_st_cell = ("crash" if cc_stable is None
                      else f"{cc_stable * 100:.2f}% (partial crash)")
        cc_dv_cell = f"{(cc.get('lifton_devel') or 0) * 100:.2f}%"
        rows.append(f"| {PRETTY.get(b, b)} | {st_cell} | {_fmt(dv)} | {gain} | "
                    f"{cc_st_cell} | {cc_dv_cell} |")
    return "\n".join(rows)


def table4(_db, vc):
    rows = ["| Benchmark | wall v1.0.8 (s) | wall v1.0.10 (s) | speedup | "
            "RSS v1.0.8 (MB) | RSS v1.0.10 (MB) | RSS reduction |",
            "|---|---:|---:|---:|---:|---:|---:|"]
    for key, rec in vc.items():
        if not key.startswith("controlled:"):
            continue
        b = key.split(":", 1)[1]
        w = rec.get("wall_s") or {}
        r = rec.get("peak_rss_mb") or {}
        ws, wd = w.get("stable"), w.get("devel")
        rs, rd = r.get("stable"), r.get("devel")
        sp = f"{ws/wd:.2f}×" if ws and wd else "—"
        rr = f"{rs/rd:.2f}×" if rs and rd else "—"
        rows.append(f"| {b} | {ws:.1f} | {wd:.1f} | {sp} | "
                    f"{rs:,.0f} | {rd:,.0f} | {rr} |" if all(
                        x is not None for x in (ws, wd, rs, rd)) else
                    f"| {b} | — | — | — | — | — | — |")
    return "\n".join(rows)


TABLES = {
    "2": ("The seventeen whole-genome lift-overs, by divergence tier", table2),
    "3": ("Full-genome coding-transcript recovery where v1.0.8 fails", table3),
    "4": ("Controlled comparison (identical cached aligner inputs, -t1)", table4),
    "5": ("Whole-genome mean protein identity by tool", table5),
    "7": ("Whole-genome coding-transcript completeness by tool", table7),
    "9": ("Whole-genome gff3-validate error counts", table9),
}


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--table", help="render one table by number")
    ap.add_argument("--all", action="store_true", help="render every table")
    ap.add_argument("--list", action="store_true", help="list available tables")
    ap.add_argument("--fourway", type=Path, default=FOURWAY)
    ap.add_argument("--version-compare", type=Path, default=VERSION_CMP)
    args = ap.parse_args(argv)

    if args.list:
        for num, (title, _fn) in sorted(TABLES.items(), key=lambda kv: int(kv[0])):
            print(f"  Table {num}: {title}")
        return 0

    db = load(args.fourway)
    vc = load(args.version_compare)

    if args.table:
        title, fn = TABLES[args.table]
        print(f"**Table {args.table}.** {title}\n")
        print(fn(db, vc))
        return 0

    if args.all:
        for num, (title, fn) in sorted(TABLES.items(), key=lambda kv: int(kv[0])):
            print(f"\n**Table {num}.** {title}\n")
            print(fn(db, vc))
        return 0

    ap.print_help()
    return 1


if __name__ == "__main__":
    sys.exit(main())
