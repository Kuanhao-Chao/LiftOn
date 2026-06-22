#!/usr/bin/env python
"""Render benchmarks/compare/version_compare.results.json into a comprehensive
markdown report (+ a flattened JSON) comparing LiftOn v1.0.8 (stable) vs devel.

Usage (repo root): python -m benchmarks.compare.version_report
Writes v1_0_8_vs_devel_report.md and v1_0_8_vs_devel_report.json next to this.
"""
from __future__ import annotations

import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
RESULTS = HERE / "version_compare.results.json"
MD = HERE / "v1_0_8_vs_devel_report.md"
JSON_OUT = HERE / "v1_0_8_vs_devel_report.json"

ARM_ORDER = ["controlled", "fresh", "full"]
DS_ORDER = ["drosophila", "mouse_to_rat", "rice", "arabidopsis", "bee"]


def _fmt(x, nd=5):
    if x is None:
        return "n/a"
    if isinstance(x, float):
        return f"{x:.{nd}f}".rstrip("0").rstrip(".") if nd else str(x)
    return str(x)


def _signed(x, nd=6):
    if x is None:
        return "n/a"
    return f"{x:+.{nd}f}"


def _sorted_recs(db, arm):
    recs = [r for k, r in db.items() if r["arm"] == arm]
    recs.sort(key=lambda r: DS_ORDER.index(r["benchmark"]) if r["benchmark"] in DS_ORDER else 99)
    return recs


def _speedup(stable, devel):
    if not stable or not devel:
        return "n/a"
    return f"{stable / devel:.2f}x" if devel else "n/a"


def _revalidate(db):
    """Recompute validity for every record from the on-disk output GFFs with the
    authoritative parser (records written before the parser fix carry crude
    counts; the GFFs persist, so re-derive)."""
    from . import version_compare as vc
    for r in db.values():
        arm, bid = r["arm"], r["benchmark"]
        new = {}
        for v in r["versions"]:
            gff = HERE / "work" / bid / "_version_ab" / arm / v / f"{v}.gff3"
            if gff.exists():
                try:
                    new[v] = vc.validate_gff(gff)
                except Exception:
                    new[v] = r.get("validity", {}).get(v, {})
        if new:
            r["validity"] = new
    return db


def main():
    if not RESULTS.exists():
        print(f"no results at {RESULTS}")
        return 1
    db = json.loads(RESULTS.read_text())
    db = _revalidate(db)

    L = []
    L.append("# LiftOn v1.0.8 (previous stable release) vs `devel` HEAD — comprehensive comparison\n")
    L.append("**stable** = the real LiftOn **v1.0.8** binary (git tag `v1.0.8` = `e503643d`, "
             "what `pip install lifton` delivers), built into an isolated `lifton_stable` conda env. "
             "**devel** = current `devel` HEAD (Iteration 19, `4496dd5`), **41 commits ahead**. "
             "**devel_legacy** = the devel binary run with `--legacy-merge --full-dp-align --gene-only` "
             "(emulates v1.0.8's *algorithm* defaults) — used only on the controlled arm to attribute "
             "deltas.\n")
    L.append("Accuracy is the **version-agnostic neutral re-score** (`benchmarks.compare.evaluator`): "
             "each predicted protein is re-aligned to the reference protein with LiftOn's own parasail "
             "kernel, identically for both versions — it never trusts either tool's self-reported "
             "`protein_identity`. Validity uses the **devel** `gff3-validate` on both outputs (one "
             "consistent yardstick). `-V` is NOT a version discriminator (devel's `__version__` still "
             "reads `v1.0.8`); provenance is pinned on each install's git SHA + `lifton.__file__`.\n")

    # ---- executive summary ----
    L.append("## Executive summary\n")
    ctrl = _sorted_recs(db, "controlled")
    if ctrl:
        L.append("Controlled arm (cached `-L`/`-M`, `-t 1` both — isolates LiftOn's in-process logic):\n")
        L.append("| dataset | neutral mean PI stable→devel (Δ) | impr/regr | recovered +dev/−dev | wall stable→devel | RSS stable→devel |")
        L.append("|---|---|---|---|---|---|")
        for r in ctrl:
            nm, d = r["neutral_mean_pi"], r.get("delta_devel_vs_stable", {})
            L.append("| {} | {} → {} ({}) | {}/{} | +{}/−{} | {}→{}s | {}→{} MB |".format(
                r["benchmark"], _fmt(nm.get("stable")), _fmt(nm.get("devel")),
                _signed(r.get("neutral_mean_pi_delta")),
                d.get("improved", "?"), d.get("regressed", "?"),
                r.get("recovered_gained_by_devel", "?"), r.get("recovered_lost_by_devel", "?"),
                _fmt(r["wall_s"].get("stable"), 1), _fmt(r["wall_s"].get("devel"), 1),
                _fmt(r["peak_rss_mb"].get("stable"), 0), _fmt(r["peak_rss_mb"].get("devel"), 0)))
        L.append("")

    # ---- per-arm detail ----
    for arm in ARM_ORDER:
        recs = _sorted_recs(db, arm)
        if not recs:
            continue
        L.append(f"## Arm: {arm}\n")

        crashed_recs = [(r["benchmark"], list(r.get("crashed", {}))) for r in recs if r.get("crashed")]
        if crashed_recs:
            L.append("> **Recorded crashes:** " + "; ".join(
                f"`{bid}` — {', '.join(cv)} crashed" for bid, cv in crashed_recs) +
                ". A crashed version's cells show `n/a`; the surviving version still scored. "
                "(v1.0.8 crashes on some full genomes with an unhandled "
                "`gffutils.FeatureNotFoundError` whose error handler itself raises "
                "`TypeError: __str__ returned non-string` — a robustness bug devel fixed.)\n")

        L.append("### Accuracy (neutral re-score)\n")
        L.append("| dataset | n coding | mean PI stable | mean PI devel | Δ mean | median s→d | %identical s→d | improved | regressed | net/txn | self−neutral bias (s/d) |")
        L.append("|---|---|---|---|---|---|---|---|---|---|---|")
        for r in recs:
            nm, md, pid = r["neutral_mean_pi"], r["neutral_median_pi"], r["neutral_pct_identical"]
            d = r.get("delta_devel_vs_stable", {})
            sb = r.get("self_vs_neutral_bias", {})
            L.append("| {} | {} | {} | {} | {} | {}→{} | {}→{} | {} | {} | {} | {}/{} |".format(
                r["benchmark"], r.get("n_reference_coding", "?"),
                _fmt(nm.get("stable")), _fmt(nm.get("devel")),
                _signed(r.get("neutral_mean_pi_delta")),
                _fmt(md.get("stable")), _fmt(md.get("devel")),
                _fmt(pid.get("stable")), _fmt(pid.get("devel")),
                d.get("improved", "?"), d.get("regressed", "?"),
                _signed(d.get("net_per_transcript")),
                _signed(sb.get("stable"), 4), _signed(sb.get("devel"), 4)))
        L.append("")

        L.append("### Completeness\n")
        L.append("| dataset | recovered coding stable→devel | completeness_coding s→d | feature_total s→d | recovered +dev/−dev |")
        L.append("|---|---|---|---|---|")
        for r in recs:
            nr, cc, ft = r["n_recovered_coding"], r["completeness_coding"], r["completeness_feature_total"]
            L.append("| {} | {}→{} | {}→{} | {}→{} | +{}/−{} |".format(
                r["benchmark"], nr.get("stable"), nr.get("devel"),
                _fmt(cc.get("stable")), _fmt(cc.get("devel")),
                _fmt(ft.get("stable")), _fmt(ft.get("devel")),
                r.get("recovered_gained_by_devel", "?"), r.get("recovered_lost_by_devel", "?")))
        L.append("")
        # gene-like feature-type recovery (where it differs) — fresh/full arms
        if arm in ("fresh", "full"):
            L.append("Per-type reference-feature recovery (stable→devel `n_recovered`), types where they differ:\n")
            L.append("| dataset | feature type | n_reference | stable recovered | devel recovered |")
            L.append("|---|---|---|---|---|")
            for r in recs:
                fc = r.get("feature_census", {})
                s_by = fc.get("stable", {}) or {}
                d_by = fc.get("devel", {}) or {}
                types = sorted(set(s_by) | set(d_by))
                for t in types:
                    if t == "_overall_":
                        continue
                    s_rec = (s_by.get(t) or {}).get("n_recovered", 0)
                    d_rec = (d_by.get(t) or {}).get("n_recovered", 0)
                    n_ref = (d_by.get(t) or s_by.get(t) or {}).get("n_reference", 0)
                    if s_rec != d_rec and n_ref:
                        L.append(f"| {r['benchmark']} | {t} | {n_ref} | {s_rec} | {d_rec} |")
            L.append("")

        L.append("### Performance\n")
        L.append("| dataset | wall stable | wall devel | wall speedup | RSS stable | RSS devel | CPU stable | CPU devel |")
        L.append("|---|---|---|---|---|---|---|---|")
        for r in recs:
            w, rss, cpu = r["wall_s"], r["peak_rss_mb"], r["cpu_s"]
            L.append("| {} | {}s | {}s | {} | {} MB | {} MB | {}s | {}s |".format(
                r["benchmark"], w.get("stable"), w.get("devel"),
                _speedup(w.get("stable"), w.get("devel")),
                rss.get("stable"), rss.get("devel"), cpu.get("stable"), cpu.get("devel")))
        L.append("")

        L.append("### Validity (devel `gff3-validate`, one yardstick on both outputs)\n")
        L.append("| dataset | stable errors/warnings | devel errors/warnings |")
        L.append("|---|---|---|")
        for r in recs:
            val = r.get("validity", {})

            def _vc(d):
                e = d.get("n_errors", d.get("n_error_lines", "?"))
                w = d.get("n_warnings", d.get("n_warning_lines", "?"))
                return f"{e}/{w}"
            L.append("| {} | {} | {} |".format(
                r["benchmark"], _vc(val.get("stable", {})), _vc(val.get("devel", {}))))
        L.append("")

        # attribution (controlled only)
        if arm == "controlled" and any("delta_devel_legacy_vs_stable" in r for r in recs):
            L.append("### Attribution (controlled): devel_legacy emulates v1.0.8's algorithm defaults\n")
            L.append("`devel_legacy` = devel + `--legacy-merge --full-dp-align --gene-only`. "
                     "**devel_legacy vs stable** = the residual from the NON-flag-gated changes "
                     "(Phase-4 gene-ID fix + Phase-13.5C canonical writer). "
                     "**devel vs devel_legacy** = the flag-gated algorithm promotions "
                     "(best-of-outcome merge + band-everything alignment + gene-like lift).\n")
            L.append("| dataset | mean PI stable | mean PI devel_legacy | mean PI devel | residual (legacy−stable) impr/regr | promotions (devel−legacy) impr/regr |")
            L.append("|---|---|---|---|---|---|")
            for r in recs:
                nm = r["neutral_mean_pi"]
                rs = r.get("delta_devel_legacy_vs_stable", {})
                pr = r.get("delta_devel_vs_devel_legacy", {})
                L.append("| {} | {} | {} | {} | {}/{} (net {}) | {}/{} (net {}) |".format(
                    r["benchmark"], _fmt(nm.get("stable")), _fmt(nm.get("devel_legacy")),
                    _fmt(nm.get("devel")),
                    rs.get("improved", "?"), rs.get("regressed", "?"), _signed(rs.get("net_per_transcript")),
                    pr.get("improved", "?"), pr.get("regressed", "?"), _signed(pr.get("net_per_transcript"))))
            L.append("")

    MD.write_text("\n".join(L) + "\n")
    JSON_OUT.write_text(json.dumps(db, indent=2))
    print(f"wrote {MD}\nwrote {JSON_OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
