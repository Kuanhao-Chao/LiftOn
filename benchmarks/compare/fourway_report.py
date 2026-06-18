#!/usr/bin/env python
"""Render benchmarks/compare/fourway_results.json into a comprehensive markdown
report (+ a flattened JSON) comparing FOUR tools per benchmark:
standalone **Liftoff** | standalone **miniprot** | **LiftOn v1.0.8** | **LiftOn devel**.

Usage (repo root): python -m benchmarks.compare.fourway_report
Writes fourway_report.md and fourway_report.json next to this.
"""
from __future__ import annotations

import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
RESULTS = HERE / "fourway_results.json"
MD = HERE / "fourway_report.md"
JSON_OUT = HERE / "fourway_report.json"

TOOLS = ["liftoff", "miniprot", "lifton_stable", "lifton_devel"]
TOOL_LABEL = {"liftoff": "Liftoff", "miniprot": "miniprot",
              "lifton_stable": "LiftOn v1.0.8", "lifton_devel": "LiftOn devel"}
# gene-like feature types that distinguish devel's gene-like lift from v1.0.8
GENE_LIKE = ("pseudogene", "ncRNA", "lnc_RNA", "lncRNA", "tRNA", "rRNA",
             "snoRNA", "snRNA", "miRNA", "ncRNA_gene", "transposable_element_gene")


def _fmt(x, nd=5):
    if x is None:
        return "n/a"
    if isinstance(x, float):
        s = f"{x:.{nd}f}".rstrip("0").rstrip(".")
        return s if s else "0"
    return str(x)


def _signed(x, nd=5):
    return "n/a" if x is None else f"{x:+.{nd}f}"


def _speedup(slow, fast):
    if not slow or not fast:
        return "n/a"
    return f"{slow / fast:.2f}x"


def _recs(db, mode):
    rs = [r for r in db.values() if r.get("mode") == mode]
    # same-species first, then cross-species; stable order within
    rs.sort(key=lambda r: (r.get("cross_species", False), r["benchmark"]))
    return rs


def _present(r):
    return [t for t in TOOLS if t in r.get("tools", [])]


def _section(L, r, mode):
    present = _present(r)
    name = f"{r['benchmark']} — {r['species']}"
    tag = "same-species" if not r["cross_species"] else r.get("divergence_class", "cross-species")
    dim = f" · _{r['dimension']}_" if r.get("dimension") else ""
    L.append(f"### {name}  ({tag}{dim}; {r.get('annotation_database','RefSeq')}; "
             f"n_coding={r.get('n_reference_coding')})\n")

    # accuracy + completeness
    L.append("| tool | completeness_coding | feature_total | mean PI | median PI | %identical | n recovered | wall (s) | peak RSS (MB) | errors/warn |")
    L.append("|---|---|---|---|---|---|---|---|---|---|")
    for t in present:
        val = r.get("validity", {}).get(t, {})
        ew = f"{val.get('n_errors','?')}/{val.get('n_warnings','?')}"
        L.append("| {} | {} | {} | {} | {} | {} | {} | {} | {} | {} |".format(
            TOOL_LABEL[t],
            _fmt(r["completeness_coding"].get(t)),
            _fmt(r["completeness_feature_total"].get(t)),
            _fmt(r["mean_pi"].get(t)), _fmt(r["median_pi"].get(t)),
            _fmt(r["pct_identical"].get(t)), r["n_recovered_coding"].get(t),
            _fmt(r["wall_s"].get(t), 1), _fmt(r["peak_rss_mb"].get(t), 0), ew))
    L.append("")

    # headline deltas
    dvs = r.get("devel_vs_stable", {})
    dvb = r.get("devel_vs_best_baseline", {})
    L.append(f"- **LiftOn devel − v1.0.8:** mean PI {_signed(dvs.get('meanpi'))}, "
             f"completeness {_signed(dvs.get('completeness'))}, "
             f"n_recovered {_signed(dvs.get('n_recovered'), 0) if dvs.get('n_recovered') is not None else 'n/a'}")
    L.append(f"- **LiftOn devel − best(Liftoff, miniprot):** mean PI {_signed(dvb.get('meanpi'))}, "
             f"completeness {_signed(dvb.get('completeness'))}")
    if present and r["wall_s"].get("lifton_stable") and r["wall_s"].get("lifton_devel"):
        L.append(f"- **LiftOn devel speedup vs v1.0.8:** "
                 f"{_speedup(r['wall_s']['lifton_stable'], r['wall_s']['lifton_devel'])}")
    L.append("")

    # gene-like feature-type census where devel recovers more than v1.0.8
    census = r.get("feature_census", {})
    sd = census.get("lifton_devel", {}) or {}
    ss = census.get("lifton_stable", {}) or {}
    rows = []
    for ftype in sorted(set(sd) | set(ss)):
        if ftype == "_overall_":
            continue
        d_rec = (sd.get(ftype) or {}).get("n_recovered", 0)
        s_rec = (ss.get(ftype) or {}).get("n_recovered", 0)
        n_ref = (sd.get(ftype) or ss.get(ftype) or {}).get("n_reference", 0)
        if d_rec != s_rec and n_ref:
            star = " ⬅ gene-like" if any(g.lower() in ftype.lower() for g in GENE_LIKE) else ""
            rows.append((ftype, n_ref, s_rec, d_rec, star))
    if rows:
        L.append("Feature-type recovery where LiftOn devel ≠ v1.0.8 (gene-like lift):\n")
        L.append("| feature type | n_reference | v1.0.8 recovered | devel recovered | |")
        L.append("|---|---|---|---|---|")
        for ftype, n_ref, s_rec, d_rec, star in rows:
            L.append(f"| {ftype} | {n_ref} | {s_rec} | {d_rec} |{star} |")
        L.append("")


def main():
    if not RESULTS.exists():
        print(f"no results at {RESULTS}")
        return 1
    db = json.loads(RESULTS.read_text())

    L = []
    L.append("# 4-way tool comparison — Liftoff vs miniprot vs LiftOn v1.0.8 vs LiftOn devel\n")
    L.append("Each tool's predicted annotation is scored by the **same version-agnostic neutral "
             "evaluator** (`benchmarks.compare.evaluator`): every predicted protein is re-aligned to "
             "the reference protein with LiftOn's own parasail kernel, and completeness is measured by "
             "reference-id recovery — identically for all four tools. The two LiftOn columns run on the "
             "**same** standalone Liftoff (`-L`) + miniprot (`-M`) outputs that the Liftoff/miniprot "
             "columns are scored on, with the **neutral common flag set only** (`-t 1 -copies -ad -g -L "
             "-M -o`; no devel-only flags, so the v1.0.8 binary is valid). So LiftOn is measured against "
             "exactly the two inputs it combines.\n")
    L.append("Provenance: LiftOn v1.0.8 = the real `v1.0.8` tag binary (`e503643d`, no `--native`) in an "
             "isolated `lifton_stable` env; LiftOn devel = current `devel` HEAD. miniprot's "
             "`completeness_feature_total` is **n/a** by construction (its `MP*` ids never match "
             "reference feature ids; its protein identity and coding-completeness are still scored).\n")

    n_sub = len(_recs(db, "subset"))
    n_full = len(_recs(db, "full"))
    L.append(f"_Records: {n_sub} subset benchmark(s), {n_full} full-genome headline(s)._\n")

    # ---- compact executive summary (deltas only) ----
    L.append("## Executive summary — the two headline deltas\n")
    L.append("| benchmark | mode | class | devel−v1.0.8 meanPI | devel−v1.0.8 compl | devel−best(LO,MP) meanPI | devel−best(LO,MP) compl |")
    L.append("|---|---|---|---|---|---|---|")
    for mode in ("subset", "full"):
        for r in _recs(db, mode):
            dvs, dvb = r.get("devel_vs_stable", {}), r.get("devel_vs_best_baseline", {})
            cls = "same" if not r["cross_species"] else r.get("divergence_class", "cross")
            L.append("| {} | {} | {} | {} | {} | {} | {} |".format(
                r["benchmark"], mode, cls,
                _signed(dvs.get("meanpi")), _signed(dvs.get("completeness")),
                _signed(dvb.get("meanpi")), _signed(dvb.get("completeness"))))
    L.append("")

    for mode, title in (("subset", "Subset 4-way matrix"), ("full", "Full-genome headlines")):
        recs = _recs(db, mode)
        if not recs:
            continue
        L.append(f"## {title}\n")
        for r in recs:
            _section(L, r, mode)

    MD.write_text("\n".join(L) + "\n")
    JSON_OUT.write_text(json.dumps(db, indent=2))
    print(f"wrote {MD}\nwrote {JSON_OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
