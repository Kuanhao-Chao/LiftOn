#!/usr/bin/env python
"""report_figures.py — figures for the LiftOn v1.0.9 *technical report*.

Distinct from blog_figures.py: the report (a) EXCLUDES the two full genomes
v1.0.8 crashed on (arabidopsis*full, rice*full) from every head-to-head
*comparison* figure — their partial-run metrics are not comparable — while
keeping a dedicated robustness/recovery figure that shows the crash itself; and
(b) carries a DEEPER, four-way performance analysis (Liftoff, miniprot, LiftOn
v1.0.8, LiftOn v1.0.9) plus the v1.0.8->v1.0.9 improvement with wall / RSS / CPU.

Reuses master_report's data loaders + palette/constants so numbers never drift
from the report tables. Touches no lifton/ source.

Usage (repo root):
  PYTHONNOUSERSITE=1 python benchmarks/compare/report_figures.py [OUTDIR]
Writes PNGs into benchmarks/compare/figures/report/ by default.
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt          # noqa: E402
from matplotlib.patches import Patch      # noqa: E402
import numpy as np                        # noqa: E402

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import master_report as mr                # noqa: E402

OUTDIR = Path(sys.argv[1]) if len(sys.argv) > 1 else (HERE / "figures" / "report")

# a full-genome run "completed" for both versions when neither collapsed below
# this coding-completeness floor (v1.0.8 crashed arabidopsis at 0.28, rice 0.77)
_COMPLETE = 0.90

# report-facing tool labels
LABEL = dict(mr.TOOL_LABEL)
LABEL["lifton_devel"] = "LiftOn v1.0.9"
LABEL["lifton_stable"] = "LiftOn v1.0.8"

# Divergence ramp — a SEQUENTIAL single-hue (purple) ladder so "darker = more
# divergent" reads at a glance, and deliberately OUTSIDE the categorical tool
# palette (blue/orange/brown/green) so a divergence band is never confused with
# a tool bar.
DIVC_COLOR = {
    "same_species": "#dadaeb",
    "cross_species": "#bcbddc",
    "close_cross_species": "#9e9ac8",
    "distant_cross_species": "#807dba",
    "very_distant_cross_species": "#6a51a3",
}

# Reader-friendly whole-genome labels (the JSON keys are pipeline ids).
_NICE = {
    "arabidopsis": "Arabidopsis", "bee": "honeybee", "rice": "rice",
    "t1_maize_b73_to_mo17": "maize", "t1_tomato_microtom_to_heinz": "tomato",
    "drosophila": "Drosophila→erecta",
    "t2_human_to_gorilla": "human→gorilla", "t2_mouse_to_caroli": "mouse→caroli",
    "t2_tomato_to_potato": "tomato→potato",
    "t3_dog_to_cat": "dog→cat", "t3_human_to_macaque": "human→macaque",
    "t3_human_to_marmoset": "human→marmoset",
    "arabidopsis_to_rice": "Arabidopsis→rice", "human_to_zebrafish": "human→zebrafish",
    "t4_drosophila_to_bee": "Drosophila→bee", "t4_human_to_chicken": "human→chicken",
    "t4_human_to_xenopus": "human→Xenopus",
}


def _flabel(r):
    sk = mr._short_key(r["key"])
    return _NICE.get(sk, sk)


# whole genomes whose residual gff3-validate errors are reference-inherited
# (organellar CDS-under-gene + duplicate organellar IDs the reference carries).
_REF_INHERITED = {"arabidopsis", "rice", "t1_maize_b73_to_mo17",
                  "t1_tomato_microtom_to_heinz", "t2_tomato_to_potato",
                  "arabidopsis_to_rice"}

plt.rcParams.update({
    "font.size": 11,
    "axes.titlesize": 11,
    "axes.titleweight": "bold",
    "axes.spines.top": False,
    "axes.spines.right": False,
    "figure.facecolor": "white",
    "savefig.facecolor": "white",
})


def _panel_title(ax, letter, text):
    ax.set_title(f"{letter}.  {text}", loc="left", fontsize=11, fontweight="bold")


def _tag(r):
    return mr._short_key(r["key"]) + ("·full" if r["mode"] == "full" else "")


DPI = 300  # print-quality source (Astro downsamples to responsive widths)


def _save(fig, name):
    OUTDIR.mkdir(parents=True, exist_ok=True)
    p = OUTDIR / name
    fig.savefig(p, format="png", dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {p}")
    return p


def _tier_bands(ax, recs, legend=False):
    """Shade light horizontal background spans, one per divergence tier, behind
    the bars of a barh figure whose rows are in DIV_ORDER (as _fw_recs returns).
    Makes 'which regime' legible at a glance — the spine of right-metric-per-tier
    reading. `recs` are the per-row records, top→bottom (y = 0..n-1, axis inverted)."""
    i, n = 0, len(recs)
    seen = []
    while i < n:
        c = recs[i].get("divergence_class")
        j = i
        while j < n and recs[j].get("divergence_class") == c:
            j += 1
        ax.axhspan(i - 0.5, j - 0.5, color=DIVC_COLOR.get(c, "#eeeeee"),
                   alpha=0.16, zorder=0, lw=0)
        if c not in seen:
            seen.append(c)
        i = j
    if legend:
        from matplotlib.patches import Patch as _P
        handles = [_P(facecolor=DIVC_COLOR[c], alpha=0.5, label=mr.DIV_LABEL.get(c, c))
                   for c in mr.DIV_ORDER if c in seen]
        ax.legend(handles=handles, title="divergence tier", fontsize=7.5,
                  title_fontsize=7.5, loc="lower left", framealpha=0.95)


# standardized crash marker: hatched stub + label in the v1.0.8 colour, so a
# v1.0.8 crash (partial or no-output) reads identically wherever it appears.
CRASH_HATCH = "////"
CRASH_EDGE = mr.TOOL_COLORS["lifton_stable"]


def _completed(r):
    """True unless this is a full-genome record v1.0.8 *crashed* on — i.e. devel
    finished but v1.0.8 aborted mid-run, emitting a partial annotation whose
    numbers are not comparable. Genomes where BOTH versions completed are kept
    even at low recall (e.g. the extreme-distance flagships), because the
    head-to-head comparison is then valid."""
    if r["mode"] != "full":
        return True
    cc = r.get("completeness_coding") or {}
    dev = cc.get("lifton_devel") or 0
    sta = cc.get("lifton_stable") or 0
    return not (dev >= _COMPLETE and sta < _COMPLETE)


def _eligible(fw):
    """All comparison-eligible records (subsets + completed fulls), in a stable
    order: subsets first (by key), then completed fulls."""
    recs = [r for r in mr._fw_recs(fw, "subset") + mr._fw_recs(fw, "full")
            if _completed(r)]
    return recs


# --------------------------------------------------------------------------- #
# ACCURACY — four-way mean protein identity (comparison-eligible)
# --------------------------------------------------------------------------- #
def fig_accuracy_fourway(fw):
    recs = [r for r in _eligible(fw)
            if any(isinstance(r["mean_pi"].get(t), float) for t in mr.TOOLS)]
    labels = [_tag(r) for r in recs]
    y = np.arange(len(recs))
    h = 0.2
    fig, ax = plt.subplots(figsize=(10, max(5, 0.5 * len(recs))))
    for i, t in enumerate(mr.TOOLS):
        vals = [r["mean_pi"].get(t) if isinstance(r["mean_pi"].get(t), float) else 0
                for r in recs]
        ax.barh(y + (1.5 - i) * h, vals, height=h, color=mr.TOOL_COLORS[t],
                label=LABEL[t])
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=8)
    ax.invert_yaxis()
    ax.set_xlabel("mean protein identity (tool-neutral parasail re-score)")
    ax.set_xlim(0.4, 1.005)
    ax.legend(fontsize=8, loc="lower left", ncol=2, framealpha=0.9)
    ax.grid(axis="x", alpha=0.3)
    ax.set_title("Mean protein identity — four tools, every completed dataset",
                 fontsize=12, fontweight="bold", loc="left")
    return _save(fig, "rfig_accuracy_fourway.png")


# --------------------------------------------------------------------------- #
# ACCURACY — divergence ladder + deltas (devel vs field)
# --------------------------------------------------------------------------- #
def fig_divergence_ladder(fw):
    classes = [c for c in mr.DIV_ORDER
               if any(r.get("divergence_class") == c and r.get("mode") == "subset"
                      for r in fw.values())]
    x = np.arange(len(classes))
    fig, ax = plt.subplots(figsize=(7.5, 5.0))
    for t in mr.TOOLS:
        ys = []
        for c in classes:
            recs = [r for r in fw.values() if r.get("divergence_class") == c
                    and r.get("mode") == "subset"]
            vals = [r["mean_pi"].get(t) for r in recs
                    if isinstance(r["mean_pi"].get(t), float)]
            ys.append(sum(vals) / len(vals) if vals else np.nan)
        # v1.0.8 ≈ v1.0.9 on this axis, so dash it (drawn on top) to stay visible
        st = dict(lw=2.3, ms=8, ls="-", zorder=3)
        if t == "lifton_stable":
            st.update(ls=(0, (5, 2)), ms=6, zorder=5)
        elif t == "lifton_devel":
            st.update(zorder=4)
        ax.plot(x, ys, marker="o", color=mr.TOOL_COLORS[t], label=LABEL[t], **st)
    ax.set_xticks(x)
    ax.set_xticklabels([mr.DIV_LABEL.get(c, c) for c in classes],
                       rotation=15, ha="right", fontsize=9)
    ax.set_ylabel("mean protein identity")
    ax.legend(fontsize=9, loc="lower left", framealpha=0.9)
    ax.grid(alpha=0.3)
    ax.set_title("The lead over the baselines widens as species diverge",
                 fontsize=12, fontweight="bold", loc="left")
    return _save(fig, "rfig_divergence_ladder.png")


def fig_devel_vs_field(fw, vc):
    fig = plt.figure(figsize=(13.5, 7.5))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.35, 0.8], wspace=0.26)
    axB = fig.add_subplot(gs[0, 0])
    axC = fig.add_subplot(gs[0, 1])

    # (A) devel - best baseline, comparison-eligible, grouped by divergence.
    #     HORIZONTAL bars so all ~27 dataset names read cleanly (no rotation).
    recs = [r for r in _eligible(fw)]
    pts = []
    for r in recs:
        d = (r.get("devel_vs_best_baseline") or {}).get("meanpi")
        if isinstance(d, float):
            pts.append((r, d))
    pts.sort(key=lambda rd: (
        mr.DIV_ORDER.index(rd[0].get("divergence_class"))
        if rd[0].get("divergence_class") in mr.DIV_ORDER else 99, rd[1]))
    labels = [_tag(r) for r, _ in pts]
    vals = [d for _, d in pts]
    colors = [DIVC_COLOR.get(r.get("divergence_class"), "#999999") for r, _ in pts]
    by = np.arange(len(pts))
    axB.barh(by, vals, color=colors)
    axB.axvline(0, color="k", lw=0.8)
    axB.set_yticks(by)
    axB.set_yticklabels(labels, fontsize=8)
    axB.invert_yaxis()
    axB.margins(x=0.10)
    axB.set_xlabel("Δ mean protein identity vs best of (Liftoff, miniprot)")
    axB.grid(axis="x", alpha=0.3)
    for i, (r, d) in enumerate(pts):
        if d < 0:
            axB.annotate(f"{d:+.4f}", xy=(d, i), xytext=(-4, 0),
                         textcoords="offset points", ha="right", va="center",
                         fontsize=7.5, color="#e45756", fontweight="bold")
    legend_handles = [Patch(facecolor=DIVC_COLOR[c], label=mr.DIV_LABEL[c])
                      for c in mr.DIV_ORDER if c in DIVC_COLOR]
    axB.legend(handles=legend_handles, fontsize=8, loc="upper left",
               framealpha=0.9, title="divergence class", title_fontsize=8)
    _nwin = sum(1 for _, d in pts if d >= 0)
    _panel_title(axB, "A",
                 f"LiftOn v1.0.9 beats the best single method on {_nwin}/{len(pts)}")

    # (C) controlled per-transcript improved vs regressed vs v1.0.8
    ctrl = sorted([r for r in vc.values() if r.get("arm") == "controlled"
                   and ((r.get("delta_devel_vs_stable") or {}).get("improved", 0)
                        + (r.get("delta_devel_vs_stable") or {}).get("regressed", 0)) > 0],
                  key=lambda r: r["benchmark"])
    cl = [r["benchmark"] for r in ctrl]
    imp = [(r.get("delta_devel_vs_stable") or {}).get("improved", 0) for r in ctrl]
    reg = [(r.get("delta_devel_vs_stable") or {}).get("regressed", 0) for r in ctrl]
    cy = np.arange(len(cl))
    w = 0.4
    bi = axC.bar(cy - w / 2, imp, w, color="#54a24b", label="improved")
    br = axC.bar(cy + w / 2, reg, w, color="#e45756", label="regressed")
    axC.bar_label(bi, fmt="%d", fontsize=8, padding=2)
    axC.bar_label(br, fmt="%d", fontsize=8, padding=2)
    axC.set_xticks(cy)
    axC.set_xticklabels(cl, rotation=20, ha="right", fontsize=8.5)
    axC.set_ylabel("transcripts vs v1.0.8")
    axC.legend(fontsize=8, loc="upper right")
    axC.grid(axis="y", alpha=0.3)
    axC.margins(y=0.18)
    _panel_title(axC, "B", "Same inputs: many improved, ~none regressed")

    return _save(fig, "rfig_devel_vs_field.png")


# --------------------------------------------------------------------------- #
# COMPLETENESS (comparison-eligible)
# --------------------------------------------------------------------------- #
def fig_completeness(fw):
    recs = [r for r in _eligible(fw)
            if any(isinstance(r["completeness_coding"].get(t), float) for t in mr.TOOLS)]
    labels = [_tag(r) for r in recs]
    y = np.arange(len(recs))
    h = 0.2
    fig, ax = plt.subplots(figsize=(10, max(5, 0.5 * len(recs))))
    for i, t in enumerate(mr.TOOLS):
        vals = [(r["completeness_coding"].get(t) or 0) * 100
                if isinstance(r["completeness_coding"].get(t), float) else 0
                for r in recs]
        ax.barh(y + (1.5 - i) * h, vals, height=h, color=mr.TOOL_COLORS[t],
                label=LABEL[t])
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=8)
    ax.invert_yaxis()
    ax.set_xlabel("coding-transcript completeness (% of reference)")
    ax.set_xlim(0, 105)
    ax.axvline(100, color="k", lw=0.6, ls=":", alpha=0.5)
    ax.legend(fontsize=8, loc="lower left", ncol=2, framealpha=0.9)
    ax.grid(axis="x", alpha=0.3)
    ax.set_title("Coding completeness — four tools, every completed dataset",
                 fontsize=12, fontweight="bold", loc="left")
    return _save(fig, "rfig_completeness.png")


# --------------------------------------------------------------------------- #
# ROBUSTNESS — full-genome crash recovery + gene-like breadth (the crashed
# genomes are the POINT here, so they are shown deliberately)
# --------------------------------------------------------------------------- #
def fig_robustness(fw):
    """HEADLINE robustness figure. (A) the FIVE full RefSeq genomes v1.0.8 cannot
    finish — two it aborts partway (Arabidopsis 28%, rice 77%) and three it
    produces no scorable output on at all (maize, two tomato pairs) — every one
    completed by v1.0.9. (B) the gene-like feature types v1.0.8 drops by design
    that v1.0.9 lifts (full Arabidopsis)."""
    fig = plt.figure(figsize=(12.6, 5.6))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.05, 1.0], wspace=0.32)
    axA = fig.add_subplot(gs[0, 0])
    axB = fig.add_subplot(gs[0, 1])

    # (A) ALL five genomes v1.0.8 crashes on. A partial crash leaves a scorable
    # n_recovered (Arabidopsis/rice); a no-output crash leaves None — shown as a
    # standardized hatched "crash (no output)" stub so the count of five reads.
    full = [r for r in mr._fw_recs(fw, "full") if _is_stable_crash(r)]
    labels = [_flabel(r) for r in full]
    dev = [r["n_recovered_coding"]["lifton_devel"] for r in full]
    sta = [r["n_recovered_coding"].get("lifton_stable") for r in full]
    x = np.arange(len(full))
    w = 0.38
    ymax = max(dev) if dev else 1
    axA.bar(x - w / 2, [s if isinstance(s, int) else 0 for s in sta], w,
            color=mr.TOOL_COLORS["lifton_stable"], label="LiftOn v1.0.8")
    axA.bar(x + w / 2, dev, w, color=mr.TOOL_COLORS["lifton_devel"],
            label="LiftOn v1.0.9")
    for xi, r, s in zip(x, full, sta):
        dpc = r["completeness_coding"].get("lifton_devel")
        if isinstance(dpc, float):
            axA.annotate(f"{dpc*100:.0f}%",
                         (xi + w / 2, r["n_recovered_coding"]["lifton_devel"]),
                         ha="center", va="bottom", fontsize=8, fontweight="bold")
        if isinstance(s, int):                       # partial crash → show the % it reached
            spc = r["completeness_coding"].get("lifton_stable")
            if isinstance(spc, float):
                axA.annotate(f"crash\n{spc*100:.0f}%", (xi - w / 2, s), ha="center",
                             va="bottom", fontsize=7, color=CRASH_EDGE,
                             fontweight="bold", linespacing=0.95)
        else:                                        # no-output crash → hatched stub
            axA.bar(xi - w / 2, ymax * 0.02, w, color="none", edgecolor=CRASH_EDGE,
                    hatch=CRASH_HATCH, lw=0.8, zorder=3)
            axA.annotate("crash\n(no\noutput)", (xi - w / 2, ymax * 0.025),
                         ha="center", va="bottom", fontsize=6.6, color=CRASH_EDGE,
                         fontweight="bold", linespacing=0.92)
    axA.set_xticks(x)
    axA.set_xticklabels(labels, rotation=18, ha="right", fontsize=8.5)
    axA.set_ylabel("coding transcripts recovered")
    axA.legend(fontsize=8, loc="upper right")
    axA.grid(axis="y", alpha=0.3)
    axA.margins(y=0.18)
    _panel_title(axA, "A", "v1.0.9 finishes five full RefSeq genomes v1.0.8 cannot")

    arab = next((r for r in mr._fw_recs(fw, "full")
                 if mr._short_key(r["key"]) == "arabidopsis"), None)
    feats_order = ["pseudogene", "lnc_RNA", "tRNA", "ncRNA", "snoRNA",
                   "antisense_RNA", "rRNA"]
    rows = []
    if arab:
        fc = arab.get("feature_census", {})
        sta_c, dev_c = fc.get("lifton_stable", {}), fc.get("lifton_devel", {})
        for ft in feats_order:
            ns = (sta_c.get(ft) or {}).get("n_recovered")
            nd = (dev_c.get(ft) or {}).get("n_recovered")
            if isinstance(ns, int) and isinstance(nd, int) and nd - ns > 0:
                rows.append((ft, ns, nd))
    if rows:
        y = np.arange(len(rows))
        h = 0.38
        axB.barh(y + h / 2, [r[1] for r in rows], h,
                 color=mr.TOOL_COLORS["lifton_stable"], label="LiftOn v1.0.8")
        axB.barh(y - h / 2, [r[2] for r in rows], h,
                 color=mr.TOOL_COLORS["lifton_devel"], label="LiftOn v1.0.9")
        for yi, (ft, ns, nd) in zip(y, rows):
            axB.annotate(f"+{nd-ns:,}", (nd, yi - h / 2), ha="left", va="center",
                         fontsize=8, xytext=(3, 0), textcoords="offset points")
        axB.set_yticks(y)
        axB.set_yticklabels([r[0] for r in rows], fontsize=9)
        axB.invert_yaxis()
        axB.set_xlabel("features recovered (Arabidopsis, full genome)")
        axB.legend(fontsize=8, loc="lower right")
        axB.grid(axis="x", alpha=0.3)
        axB.margins(x=0.12)
    _panel_title(axB, "B", "and lifts the gene-like types v1.0.8 drops by design")
    fig.suptitle("Whole-genome robustness and feature breadth",
                 fontsize=12.5, fontweight="bold", y=1.02)
    return _save(fig, "rfig_robustness.png")


# --------------------------------------------------------------------------- #
# VALIDITY (comparison-eligible)
# --------------------------------------------------------------------------- #
def fig_validity(fw):
    recs = [r for r in _eligible(fw)
            if isinstance(mr._val_errs((r.get("validity") or {}).get("lifton_stable")), int)
            and isinstance(mr._val_errs((r.get("validity") or {}).get("lifton_devel")), int)]
    labels = [_tag(r) for r in recs]
    sta = [mr._val_errs(r["validity"]["lifton_stable"]) for r in recs]
    dev = [mr._val_errs(r["validity"]["lifton_devel"]) for r in recs]
    y = np.arange(len(recs))
    h = 0.38
    fig, ax = plt.subplots(figsize=(9, max(5, 0.5 * len(recs))))
    ax.barh(y - h / 2, sta, h, color=mr.TOOL_COLORS["lifton_stable"], label="LiftOn v1.0.8")
    ax.barh(y + h / 2, dev, h, color=mr.TOOL_COLORS["lifton_devel"], label="LiftOn v1.0.9")
    for yi, s, d in zip(y, sta, dev):
        if d > s:
            ax.annotate(f"{s}→{d}", (max(s, d), yi), ha="left", va="center",
                        fontsize=7.5, color="#e45756", fontweight="bold",
                        xytext=(4, 0), textcoords="offset points")
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=8)
    ax.invert_yaxis()
    ax.set_xlabel("gff3-validate error count (lower is better)")
    ax.legend(fontsize=8, loc="lower right")
    ax.grid(axis="x", alpha=0.3)
    ax.margins(x=0.12)
    ax.set_title("Output validity — LiftOn v1.0.9 vs v1.0.8",
                 fontsize=12, fontweight="bold", loc="left")
    return _save(fig, "rfig_validity.png")


# --------------------------------------------------------------------------- #
# PERFORMANCE 1 — four-way resource footprint (fresh end-to-end, subsets)
# --------------------------------------------------------------------------- #
def fig_perf_fourway(fw):
    recs = [r for r in mr._fw_recs(fw, "subset")
            if all(isinstance(r["wall_s"].get(t), (int, float)) for t in mr.TOOLS)
            and all(isinstance(r["peak_rss_mb"].get(t), (int, float)) for t in mr.TOOLS)]
    # one shared order for both panels: by v1.0.8 peak RSS (giant-gene blowups on top)
    recs.sort(key=lambda r: r["peak_rss_mb"]["lifton_stable"], reverse=True)
    labels = [mr._short_key(r["key"]) for r in recs]
    y = np.arange(len(recs))
    h = 0.2

    fig = plt.figure(figsize=(12.5, max(6, 0.62 * len(recs))))
    gs = fig.add_gridspec(1, 2, wspace=0.30)
    axW = fig.add_subplot(gs[0, 0])
    axR = fig.add_subplot(gs[0, 1])

    def panel(ax, key, xlabel, title, gb=False):
        for i, t in enumerate(mr.TOOLS):
            vals = [r[key][t] / (1000 if gb else 1) for r in recs]
            ax.barh(y + (1.5 - i) * h, vals, height=h, color=mr.TOOL_COLORS[t],
                    label=LABEL[t])
        ax.set_yticks(y)
        ax.set_yticklabels(labels, fontsize=8)
        ax.invert_yaxis()
        ax.set_xscale("log")
        ax.set_xlabel(xlabel)
        ax.grid(axis="x", alpha=0.3, which="both")
        ax.set_title(title, fontsize=11, fontweight="bold", loc="left")

    panel(axW, "wall_s", "wall-clock (s, log scale)",
          "A.  End-to-end wall-clock — four tools")
    axW.legend(fontsize=8, loc="lower right", ncol=2, framealpha=0.9)
    panel(axR, "peak_rss_mb", "peak resident memory (GB, log scale)",
          "B.  Peak memory — LiftOn v1.0.9 is lowest of all four", gb=True)
    fig.suptitle("Fresh end-to-end runs (each tool does its own alignment) — "
                 f"{len(recs)} subset benchmarks",
                 fontsize=11.5, fontweight="bold", y=1.005)
    return _save(fig, "rfig_perf_fourway.png")


# --------------------------------------------------------------------------- #
# PERFORMANCE 2 — the v1.0.8 -> v1.0.9 improvement (controlled arm: cached aligner
# inputs at -t1, isolates LiftOn's own code). wall / RSS / CPU + gain factors.
# --------------------------------------------------------------------------- #
def fig_perf_improvement(vc):
    recs = [r for r in vc.values() if r.get("arm") == "controlled"
            and isinstance((r.get("wall_s") or {}).get("stable"), (int, float))
            and isinstance((r.get("peak_rss_mb") or {}).get("stable"), (int, float))]
    recs.sort(key=lambda r: r["peak_rss_mb"]["stable"] / r["peak_rss_mb"]["devel"],
              reverse=True)
    labels = [r["benchmark"] for r in recs]
    x = np.arange(len(recs))
    w = 0.38
    cs, cd = mr.TOOL_COLORS["lifton_stable"], mr.TOOL_COLORS["lifton_devel"]

    fig = plt.figure(figsize=(13.5, 4.6))
    gs = fig.add_gridspec(1, 3, wspace=0.30)
    axW, axR, axF = (fig.add_subplot(gs[0, i]) for i in range(3))

    def grouped(ax, key, ylabel, title, log, fmt):
        sta = [r[key]["stable"] for r in recs]
        dev = [r[key]["devel"] for r in recs]
        ax.bar(x - w / 2, sta, w, color=cs, label="LiftOn v1.0.8")
        ax.bar(x + w / 2, dev, w, color=cd, label="LiftOn v1.0.9")
        # ratio label centered over each pair, above the taller bar (no clipping)
        for xi, s, d in zip(x, sta, dev):
            if d:
                ax.annotate(fmt(s, d), (xi, max(s, d)), ha="center", va="bottom",
                            fontsize=8.5, fontweight="bold", color=cd,
                            xytext=(0, 2), textcoords="offset points")
        ax.set_xticks(x)
        ax.set_xticklabels(labels, rotation=22, ha="right", fontsize=8.5)
        ax.set_ylabel(ylabel)
        if log:
            ax.set_yscale("log")
        ax.legend(fontsize=8, loc="upper right")
        ax.grid(axis="y", alpha=0.3)
        ax.margins(y=0.2)
        ax.set_title(title, fontsize=11, fontweight="bold", loc="left")

    grouped(axW, "wall_s", "wall-clock (s)", "A.  Wall-clock", False,
            lambda s, d: f"{s/d:.2f}×")
    grouped(axR, "peak_rss_mb", "peak memory (MB, log)", "B.  Peak memory", True,
            lambda s, d: f"{s/d:.1f}×")

    # (C) gain factors: speedup and memory-reduction side by side
    spd = [r["wall_s"]["stable"] / r["wall_s"]["devel"] for r in recs]
    mem = [r["peak_rss_mb"]["stable"] / r["peak_rss_mb"]["devel"] for r in recs]
    b1 = axF.bar(x - w / 2, spd, w, color="#9e6ebd", label="wall-clock speedup")
    b2 = axF.bar(x + w / 2, mem, w, color="#3a8f5a", label="memory reduction")
    axF.bar_label(b1, fmt="%.1f×", fontsize=7.5, padding=2)
    axF.bar_label(b2, fmt="%.1f×", fontsize=7.5, padding=2)
    axF.axhline(1, color="k", lw=0.7, ls=":")
    axF.set_xticks(x)
    axF.set_xticklabels(labels, rotation=22, ha="right", fontsize=8.5)
    axF.set_ylabel("improvement factor (×)")
    axF.legend(fontsize=8, loc="upper right")
    axF.grid(axis="y", alpha=0.3)
    axF.margins(y=0.2)
    axF.set_title("C.  Gain factors (v1.0.8 → v1.0.9)", fontsize=11,
                  fontweight="bold", loc="left")

    fig.suptitle("Bounded-memory engine: up to 23.9× less RAM and ~2× faster on "
                 "identical cached inputs (-t1, controlled arm)", fontsize=11.5,
                 fontweight="bold", y=1.02)
    return _save(fig, "rfig_perf_improvement.png")


# --------------------------------------------------------------------------- #
# FULL-GENOME FIGURES — real-world whole-genome end-to-end runs (the 6 fulls).
# These power the full-genome-focused website report; the mixed/subset figures
# above remain generated but are no longer imported by the report.
# --------------------------------------------------------------------------- #
def _full_recs(fw):
    """The full-genome records, ordered easy→hard by divergence class."""
    return mr._fw_recs(fw, "full")


def _is_stable_crash(r):
    cc = r.get("completeness_coding") or {}
    return (cc.get("lifton_devel") or 0) >= 0.90 and (cc.get("lifton_stable") or 0) < 0.90


def fig_full_accuracy(fw):
    """Whole-genome accuracy. (A) four-way mean protein identity, every genome,
    rows tier-banded. (B) Δ vs Liftoff — the DNA baseline LiftOn extends: an
    artifact-free, same-metric/same-tool-family comparison that is positive on
    all 17 and grows with divergence (the miniprot head-to-head is Fig 4)."""
    recs = _full_recs(fw)
    labels = [_flabel(r) for r in recs]
    y = np.arange(len(recs))
    h = 0.2
    fig = plt.figure(figsize=(13.6, 5.8))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.25, 0.95], wspace=0.30)
    axA = fig.add_subplot(gs[0, 0])
    axB = fig.add_subplot(gs[0, 1])

    # (A) four-way mean protein identity
    _tier_bands(axA, recs)
    for i, t in enumerate(mr.TOOLS):
        vals = [r["mean_pi"].get(t) if isinstance(r["mean_pi"].get(t), float) else 0
                for r in recs]
        axA.barh(y + (1.5 - i) * h, vals, height=h, color=mr.TOOL_COLORS[t],
                 label=LABEL[t])
    axA.set_yticks(y)
    axA.set_yticklabels(labels, fontsize=9)
    axA.invert_yaxis()
    axA.set_xlim(0, 1.02)
    axA.set_xlabel("mean protein identity (tool-neutral parasail re-score)")
    axA.grid(axis="x", alpha=0.3)
    axA.legend(fontsize=8.5, loc="lower right", ncol=2, framealpha=0.9)
    _panel_title(axA, "A", "Four tools, every whole genome (rows grouped by divergence tier)")

    # (B) devel − Liftoff: LiftOn's value-add over the pure-DNA model, which grows
    # from negligible (same-species, DNA already near-perfect) to enormous
    # (very-distant, where DNA collapses and the protein evidence carries it).
    deltas = []
    for r in recs:
        lo, dv = r["mean_pi"].get("liftoff"), r["mean_pi"].get("lifton_devel")
        deltas.append(dv - lo if None not in (lo, dv) else None)
    _tier_bands(axB, recs)
    colors = [DIVC_COLOR.get(r.get("divergence_class"), "#999999") for r in recs]
    axB.barh(y, [d or 0 for d in deltas], color=colors)
    axB.axvline(0, color="k", lw=0.8)
    axB.set_yticks(y)
    axB.set_yticklabels(labels, fontsize=9)
    axB.invert_yaxis()
    axB.margins(x=0.24)
    axB.set_xlabel("Δ mean protein identity\nvs Liftoff (the DNA baseline LiftOn extends)")
    axB.grid(axis="x", alpha=0.3)
    for yi, d in zip(y, deltas):
        if d is not None:
            axB.annotate(f"{d:+.3f}", xy=(d, yi), xytext=(4, 0),
                         textcoords="offset points", ha="left", va="center",
                         fontsize=7.5, fontweight="bold")
    _nwin = sum(1 for d in deltas if d is not None and d >= 0)
    _panel_title(axB, "B", f"Above the DNA baseline on every genome ({_nwin}/{len(recs)})")
    fig.suptitle("Whole-genome accuracy: LiftOn v1.0.9 vs the DNA and protein baselines",
                 fontsize=12.5, fontweight="bold", y=1.02)
    return _save(fig, "rfig_full_accuracy.png")


def fig_full_completeness(fw):
    recs = _full_recs(fw)
    labels = [_flabel(r) for r in recs]
    y = np.arange(len(recs))
    h = 0.2
    fig, ax = plt.subplots(figsize=(11.6, 6.6))
    _tier_bands(ax, recs)
    for i, t in enumerate(mr.TOOLS):
        vals = [(r["completeness_coding"].get(t) or 0) * 100
                if isinstance(r["completeness_coding"].get(t), float) else 0
                for r in recs]
        ax.barh(y + (1.5 - i) * h, vals, height=h, color=mr.TOOL_COLORS[t],
                label=LABEL[t])
    # mark the genomes v1.0.8 crashed on. Where it left a scorable partial
    # annotation (arabidopsis 28%, rice 77%) show that %; where it aborted with
    # no scorable whole-genome output (maize, tomato pairs) show just "crash"
    # at the axis origin rather than a misleading "0%".
    stable_y = y + (1.5 - 2) * h
    for yi, r in zip(stable_y, recs):
        if _is_stable_crash(r):
            scc = r["completeness_coding"].get("lifton_stable")
            if isinstance(scc, float):
                txt, xpos = f"v1.0.8 crash ({scc*100:.0f}%)", scc * 100
            else:
                txt, xpos = "v1.0.8 crash", 0
            ax.annotate(txt, xy=(xpos, yi), xytext=(5, 0),
                        textcoords="offset points", ha="left", va="center",
                        fontsize=7.5, color=CRASH_EDGE, fontweight="bold")
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=9)
    ax.invert_yaxis()
    ax.set_xlim(0, 112)
    ax.axvline(100, color="k", lw=0.6, ls=":", alpha=0.5)
    ax.set_xlabel("coding-transcript completeness (% of reference)")
    ax.legend(fontsize=8.5, loc="lower right", ncol=2, framealpha=0.9)
    ax.grid(axis="x", alpha=0.3)
    ax.set_title("Whole-genome completeness — at extreme distance miniprot trades "
                 "identity for recall", fontsize=11.5, fontweight="bold", loc="left")
    return _save(fig, "rfig_full_completeness.png")


def fig_full_validity(fw):
    """gff3-validate error counts, v1.0.8 vs the released v1.0.9, on the 17 whole
    genomes. v1.0.9's containment normalization drives the count to zero on 9 of
    17 and far below v1.0.8 on every genome both produced; the † plant rows carry
    only reference-inherited organellar errors the strict validator flags."""
    recs = _full_recs(fw)
    labels = [_flabel(r) + (" †" if mr._short_key(r["key"]) in _REF_INHERITED else "")
              for r in recs]
    # raw values keep None for cells where v1.0.8 crashed (no output to validate);
    # coerce to 0 only for the bar height (→ no bar), keep raw for the annotation.
    sta_raw = [mr._val_errs((r.get("validity") or {}).get("lifton_stable")) for r in recs]
    dev_raw = [mr._val_errs((r.get("validity") or {}).get("lifton_devel")) for r in recs]
    sta = [s if isinstance(s, int) else 0 for s in sta_raw]
    dev = [d if isinstance(d, int) else 0 for d in dev_raw]
    n_clean = sum(1 for d in dev_raw if d == 0)
    y = np.arange(len(recs))
    h = 0.38
    fig, ax = plt.subplots(figsize=(10.2, 6.0))
    _tier_bands(ax, recs)
    ax.barh(y - h / 2, sta, h, color=mr.TOOL_COLORS["lifton_stable"], label="LiftOn v1.0.8")
    ax.barh(y + h / 2, dev, h, color=mr.TOOL_COLORS["lifton_devel"], label="LiftOn v1.0.9")
    for yi, s, d in zip(y, sta_raw, dev_raw):
        if isinstance(s, int) and isinstance(d, int):
            ax.annotate(f"{s}→{d}", (max(s, d), yi), ha="left", va="center",
                        fontsize=8, fontweight="bold",
                        color="#3a8f5a" if d < s else CRASH_EDGE,
                        xytext=(4, 0), textcoords="offset points")
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=9)
    ax.invert_yaxis()
    ax.set_xlabel("gff3-validate error count (lower is better)")
    ax.legend(fontsize=8.5, loc="lower right")
    ax.grid(axis="x", alpha=0.3)
    ax.margins(x=0.18)
    ax.annotate("† residual errors reference-inherited (organellar CDS-under-gene + "
                "duplicate organellar IDs)", xy=(0.0, -0.155),
                xycoords="axes fraction", fontsize=7.5, color="#666666")
    ax.set_title(f"Cleaner output on every whole genome — 0 errors on {n_clean} of "
                 f"{len(recs)} (v1.0.9 vs v1.0.8)",
                 fontsize=12, fontweight="bold", loc="left")
    return _save(fig, "rfig_full_validity.png")


def fig_full_apples_to_apples(fw):
    """The honest summary. (A) per genome, the grey bar is v1.0.9's set-mean lead
    over the better single baseline; the coloured bar is its apples-to-apples lead
    over miniprot on the transcripts BOTH recover (common-set size n labelled) —
    green through the distant tier, red on four of five very-distant genomes,
    exposing the very-distant set-mean lead as a small-denominator effect.
    (B) coverage-weighted identity (recall × accuracy): miniprot leads wherever it
    recovers more transcripts."""
    recs = [r for r in _full_recs(fw)
            if (r.get("joint") or {}).get("devel_vs_miniprot_common")]
    labels = [_flabel(r) for r in recs]
    y = np.arange(len(recs))
    h = 0.38
    fig = plt.figure(figsize=(13.8, 5.8))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.08, 1.0], wspace=0.34)
    axA = fig.add_subplot(gs[0, 0])
    axB = fig.add_subplot(gs[0, 1])

    _tier_bands(axA, recs)
    setmean = [(r.get("devel_vs_best_baseline") or {}).get("meanpi") for r in recs]
    common = [r["joint"]["devel_vs_miniprot_common"]["meanpi_delta"] for r in recs]
    ncom = [r["joint"]["devel_vs_miniprot_common"].get("n_common") for r in recs]
    axA.barh(y - h / 2, [d or 0 for d in setmean], height=h, color="#9e9e9e",
             label="set-mean lead vs best(Liftoff, miniprot)")
    axA.barh(y + h / 2, [d or 0 for d in common], height=h,
             color=["#2ca02c" if (d or 0) >= 0 else "#d62728" for d in common],
             label="apples-to-apples vs miniprot (common recovered set)")
    for yi, c, n in zip(y, common, ncom):
        if isinstance(n, int):
            axA.annotate(f"n={n:,}", (c or 0, yi + h / 2),
                         xytext=(5 if (c or 0) >= 0 else -5, 0),
                         textcoords="offset points",
                         ha="left" if (c or 0) >= 0 else "right", va="center",
                         fontsize=6.6, color="#444444")
    axA.axvline(0, color="k", lw=0.8)
    axA.set_yticks(y)
    axA.set_yticklabels(labels, fontsize=9)
    axA.invert_yaxis()
    axA.margins(x=0.22)
    axA.set_xlabel("Δ mean protein identity")
    axA.grid(axis="x", alpha=0.3)
    axA.legend(fontsize=7.6, loc="lower left", framealpha=0.95)
    _panel_title(axA, "A", "Set-mean vs apples-to-apples common-set lead")

    _tier_bands(axB, recs)
    cov_dev = [(r["joint"].get("covpi") or {}).get("lifton_devel", 0) for r in recs]
    cov_mp = [(r["joint"].get("covpi") or {}).get("miniprot", 0) for r in recs]
    axB.barh(y - h / 2, cov_dev, height=h, color=mr.TOOL_COLORS["lifton_devel"], label="LiftOn v1.0.9")
    axB.barh(y + h / 2, cov_mp, height=h, color=mr.TOOL_COLORS["miniprot"], label="miniprot")
    axB.set_yticks(y)
    axB.set_yticklabels(labels, fontsize=9)
    axB.invert_yaxis()
    axB.set_xlim(0, 1.02)
    axB.set_xlabel("coverage-weighted protein identity (recall × accuracy)")
    axB.grid(axis="x", alpha=0.3)
    axB.legend(fontsize=8.5, loc="lower right", framealpha=0.9)
    _panel_title(axB, "B", "Coverage-weighted PI (recall × accuracy)")
    fig.suptitle("The joint recall-vs-identity view — locating the genuine per-transcript win",
                 fontsize=12.5, fontweight="bold", y=1.02)
    return _save(fig, "rfig_full_apples_to_apples.png")


def main():
    vc, fw, _ = mr._load()
    if not fw or not vc:
        sys.exit("ERROR: source JSON not found / empty — run the comparison harness first.")
    print(f"OUTDIR = {OUTDIR}")
    # full-genome-focused figures (imported by the website report)
    fig_full_accuracy(fw)
    fig_full_completeness(fw)
    fig_full_apples_to_apples(fw)
    fig_full_validity(fw)
    fig_robustness(fw)
    fig_perf_improvement(vc)
    # legacy mixed/subset figures (still generated; no longer imported)
    fig_accuracy_fourway(fw)
    fig_divergence_ladder(fw)
    fig_devel_vs_field(fw, vc)
    fig_completeness(fw)
    fig_validity(fw)
    fig_perf_fourway(fw)
    print("done.")


if __name__ == "__main__":
    main()
