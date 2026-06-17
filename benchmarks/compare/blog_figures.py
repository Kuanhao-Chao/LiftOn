#!/usr/bin/env python
"""blog_figures.py — render the LiftOn 2.0 announcement's composite figures.

Six publication-style figures for the website blog post, built from the SAME
source JSON as master_report.py (so the numbers stay in lockstep):

  * fig1_accuracy.png           (A) divergence ladder  (B) devel − best baseline
                                (C) controlled per-transcript improved/regressed
                                (only datasets with a non-zero change)
  * fig4_fourway_pi.png         four-way mean protein identity, every dataset that
                                BOTH versions completed (the two genomes v1.0.8
                                crashed on are excluded — its partial-run accuracy
                                there is not comparable; they appear in the
                                completeness figures instead)
  * fig6_completeness_bytool.png  four-way coding completeness, every benchmark
  * fig2_completeness.png       (A) full-genome crash recovery (v1.0.8 vs devel)
                                (B) gene-like feature recovery (Arabidopsis full)
  * fig5_validity.png           gff3-validate error counts, devel vs v1.0.8
                                (honest: the lone yeast regression is visible)
  * fig3_performance.png        (A) end-to-end wall-clock  (B) peak RSS (log)

Reuses master_report's data loaders + palette/constants; touches no `lifton/`
source and is orthogonal to the 24-cell byte-identity contract.

Usage (repo root):
  PYTHONNOUSERSITE=1 python benchmarks/compare/blog_figures.py [OUTDIR]

Writes PNGs into benchmarks/compare/figures/blog/ by default (or OUTDIR if given).
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt        # noqa: E402
from matplotlib.patches import Patch    # noqa: E402
import numpy as np                      # noqa: E402

# reuse the report's loaders + constants so figures and report never drift
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import master_report as mr              # noqa: E402

OUTDIR = Path(sys.argv[1]) if len(sys.argv) > 1 else (HERE / "figures" / "blog")

# a full-genome run "completed" for both versions when neither collapsed below
# this coding-completeness floor (v1.0.8 crashed arabidopsis at 0.28, rice at 0.77)
_COMPLETE = 0.90

# blog-facing tool labels: the devel build is announced as "LiftOn 2.0"
LABEL = dict(mr.TOOL_LABEL)
LABEL["lifton_devel"] = "LiftOn 2.0"

# divergence-class palette (easy → hard), used to color the per-benchmark bars
DIVC_COLOR = {
    "same_species": "#4c78a8",
    "cross_species": "#72b7b2",
    "close_cross_species": "#f58518",
    "distant_cross_species": "#e45756",
}

plt.rcParams.update({
    "font.size": 10,
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


def _save(fig, name):
    OUTDIR.mkdir(parents=True, exist_ok=True)
    p = OUTDIR / name
    fig.savefig(p, format="png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {p}")
    return p


def _both_completed(r):
    """A subset record always counts; a full record only if NEITHER version's
    coding completeness collapsed (i.e. v1.0.8 didn't crash partway)."""
    if r["mode"] != "full":
        return True
    cc = r.get("completeness_coding") or {}
    return (cc.get("lifton_stable") or 0) >= _COMPLETE and \
           (cc.get("lifton_devel") or 0) >= _COMPLETE


# --------------------------------------------------------------------------- #
# FIGURE 1 — accuracy (A ladder, B devel−best baseline, C per-transcript)
# --------------------------------------------------------------------------- #
def fig_accuracy(fw, vc):
    fig = plt.figure(figsize=(11.5, 9.0))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 1.25],
                          hspace=0.42, wspace=0.26)
    axA = fig.add_subplot(gs[0, 0])
    axC = fig.add_subplot(gs[0, 1])
    axB = fig.add_subplot(gs[1, :])

    # ---- (A) divergence ladder: mean PI vs class, per tool ----
    classes = [c for c in mr.DIV_ORDER
               if any(r.get("divergence_class") == c and r.get("mode") == "subset"
                      for r in fw.values())]
    x = np.arange(len(classes))
    for t in mr.TOOLS:
        ys = []
        for c in classes:
            recs = [r for r in fw.values() if r.get("divergence_class") == c
                    and r.get("mode") == "subset"]
            vals = [r["mean_pi"].get(t) for r in recs
                    if isinstance(r["mean_pi"].get(t), float)]
            ys.append(sum(vals) / len(vals) if vals else np.nan)
        axA.plot(x, ys, "o-", color=mr.TOOL_COLORS[t], label=LABEL[t],
                 lw=2.2, ms=7)
    axA.set_xticks(x)
    axA.set_xticklabels([mr.DIV_LABEL.get(c, c) for c in classes],
                        rotation=18, ha="right", fontsize=8.5)
    axA.set_ylabel("mean protein identity")
    axA.legend(fontsize=7.5, loc="lower left", framealpha=0.9)
    axA.grid(alpha=0.3)
    _panel_title(axA, "A", "The lead widens as species diverge")

    # ---- (C) controlled per-transcript improved vs regressed vs v1.0.8 ----
    # only datasets with a non-zero change (same-species bee/rice are flat 0/0)
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
    _panel_title(axC, "C", "Same inputs: many improved, ~none regressed")

    # ---- (B) devel − best baseline, per benchmark, grouped by divergence ----
    recs = mr._fw_recs(fw, "subset") + mr._fw_recs(fw, "full")
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
    bx = np.arange(len(pts))
    axB.bar(bx, vals, color=colors)
    axB.axhline(0, color="k", lw=0.8)
    axB.set_xticks(bx)
    axB.set_xticklabels(labels, rotation=55, ha="right", fontsize=7.5)
    axB.set_ylabel("Δ mean protein identity\nvs best of (Liftoff, miniprot)")
    axB.grid(axis="y", alpha=0.3)
    for i, (r, d) in enumerate(pts):
        if d < 0:
            axB.annotate(f"{mr._short_key(r['key'])}\n{d:+.4f}",
                         xy=(i, d), xytext=(i, d - 0.012),
                         ha="center", va="top", fontsize=7.5,
                         arrowprops=dict(arrowstyle="-", lw=0.6))
    legend_handles = [Patch(facecolor=DIVC_COLOR[c], label=mr.DIV_LABEL[c])
                      for c in mr.DIV_ORDER if c in DIVC_COLOR]
    axB.legend(handles=legend_handles, fontsize=8, loc="upper left", ncol=2,
               title="divergence class", title_fontsize=8)
    _panel_title(axB, "B", "LiftOn 2.0 beats the best single method on 19/20 datasets")

    return _save(fig, "fig1_accuracy.png")


# --------------------------------------------------------------------------- #
# FIGURE 4 — four-way mean PI, every dataset BOTH versions completed
# --------------------------------------------------------------------------- #
def fig_fourway_pi(fw):
    recs = [r for r in mr._fw_recs(fw, "subset") + mr._fw_recs(fw, "full")
            if _both_completed(r)
            and any(isinstance(r["mean_pi"].get(t), float) for t in mr.TOOLS)]
    if not recs:
        return None
    labels = [_tag(r) for r in recs]
    y = np.arange(len(recs))
    h = 0.2
    fig, ax = plt.subplots(figsize=(10, max(5, 0.46 * len(recs))))
    for i, t in enumerate(mr.TOOLS):
        vals = [r["mean_pi"].get(t) if isinstance(r["mean_pi"].get(t), float) else 0
                for r in recs]
        ax.barh(y + (1.5 - i) * h, vals, height=h, color=mr.TOOL_COLORS[t],
                label=LABEL[t])
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=8)
    ax.invert_yaxis()
    ax.set_xlabel("mean protein identity (tool-neutral re-score)")
    ax.set_xlim(0.4, 1.005)
    ax.legend(fontsize=8, loc="lower left", ncol=2, framealpha=0.9)
    ax.grid(axis="x", alpha=0.3)
    ax.set_title("Mean protein identity — four tools, every completed dataset",
                 fontsize=12, fontweight="bold", loc="left")
    return _save(fig, "fig4_fourway_pi.png")


# --------------------------------------------------------------------------- #
# FIGURE 6 — four-way coding completeness, every benchmark
# --------------------------------------------------------------------------- #
def fig_completeness_bytool(fw):
    recs = [r for r in mr._fw_recs(fw, "subset") + mr._fw_recs(fw, "full")
            if any(isinstance(r["completeness_coding"].get(t), float) for t in mr.TOOLS)]
    if not recs:
        return None
    labels = [_tag(r) for r in recs]
    y = np.arange(len(recs))
    h = 0.2
    fig, ax = plt.subplots(figsize=(10, max(5, 0.46 * len(recs))))
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
    ax.set_title("Coding completeness — four tools, every benchmark "
                 "(v1.0.8's short bars on arabidopsis·full / rice·full are its crashes)",
                 fontsize=10.5, fontweight="bold", loc="left")
    return _save(fig, "fig6_completeness_bytool.png")


# --------------------------------------------------------------------------- #
# FIGURE 2 — completeness & robustness (A crash recovery, B gene-like)
# --------------------------------------------------------------------------- #
def fig_completeness(fw):
    fig = plt.figure(figsize=(11.5, 5.4))
    gs = fig.add_gridspec(1, 2, width_ratios=[0.85, 1.15], wspace=0.32)
    axA = fig.add_subplot(gs[0, 0])
    axB = fig.add_subplot(gs[0, 1])

    # ---- (A) full-genome crash recovery ----
    full = [r for r in mr._fw_recs(fw, "full")
            if isinstance(r["n_recovered_coding"].get("lifton_stable"), int)
            and isinstance(r["n_recovered_coding"].get("lifton_devel"), int)
            and r["n_recovered_coding"]["lifton_devel"]
            - r["n_recovered_coding"]["lifton_stable"] > 0]
    labels = [mr._short_key(r["key"]) for r in full]
    sta = [r["n_recovered_coding"]["lifton_stable"] for r in full]
    dev = [r["n_recovered_coding"]["lifton_devel"] for r in full]
    sta_pct = [r["completeness_coding"].get("lifton_stable") for r in full]
    dev_pct = [r["completeness_coding"].get("lifton_devel") for r in full]
    x = np.arange(len(full))
    w = 0.38
    b1 = axA.bar(x - w / 2, sta, w, color=mr.TOOL_COLORS["lifton_stable"],
                 label="LiftOn v1.0.8")
    b2 = axA.bar(x + w / 2, dev, w, color=mr.TOOL_COLORS["lifton_devel"],
                 label="LiftOn 2.0")
    for rect, p in zip(b1, sta_pct):
        if isinstance(p, float):
            axA.annotate(f"{p*100:.0f}%", (rect.get_x() + rect.get_width() / 2,
                         rect.get_height()), ha="center", va="bottom", fontsize=8)
    for rect, p in zip(b2, dev_pct):
        if isinstance(p, float):
            axA.annotate(f"{p*100:.0f}%", (rect.get_x() + rect.get_width() / 2,
                         rect.get_height()), ha="center", va="bottom", fontsize=8)
    axA.set_xticks(x)
    axA.set_xticklabels(labels, fontsize=9)
    axA.set_ylabel("coding transcripts recovered")
    axA.legend(fontsize=8, loc="upper left")
    axA.grid(axis="y", alpha=0.3)
    axA.margins(y=0.16)
    _panel_title(axA, "A", "v1.0.8 crashes partway; 2.0 finishes")

    # ---- (B) gene-like feature recovery (Arabidopsis full) ----
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
        ns_v = [r[1] for r in rows]
        nd_v = [r[2] for r in rows]
        axB.barh(y + h / 2, ns_v, h, color=mr.TOOL_COLORS["lifton_stable"],
                 label="LiftOn v1.0.8")
        axB.barh(y - h / 2, nd_v, h, color=mr.TOOL_COLORS["lifton_devel"],
                 label="LiftOn 2.0")
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
    _panel_title(axB, "B", "2.0 lifts gene-like types v1.0.8 drops")

    return _save(fig, "fig2_completeness.png")


# --------------------------------------------------------------------------- #
# FIGURE 5 — output validity, devel vs v1.0.8 (honest)
# --------------------------------------------------------------------------- #
def fig_validity(fw):
    recs = [r for r in mr._fw_recs(fw, "subset") + mr._fw_recs(fw, "full")
            if isinstance(mr._val_errs((r.get("validity") or {}).get("lifton_stable")), int)
            and isinstance(mr._val_errs((r.get("validity") or {}).get("lifton_devel")), int)]
    if not recs:
        return None
    labels = [_tag(r) for r in recs]
    sta = [mr._val_errs(r["validity"]["lifton_stable"]) for r in recs]
    dev = [mr._val_errs(r["validity"]["lifton_devel"]) for r in recs]
    y = np.arange(len(recs))
    h = 0.38
    fig, ax = plt.subplots(figsize=(9, max(5, 0.46 * len(recs))))
    ax.barh(y - h / 2, sta, h, color=mr.TOOL_COLORS["lifton_stable"], label="LiftOn v1.0.8")
    ax.barh(y + h / 2, dev, h, color=mr.TOOL_COLORS["lifton_devel"], label="LiftOn 2.0")
    # mark the one dataset where 2.0 has MORE errors than v1.0.8 (honest)
    for yi, s, d, lab in zip(y, sta, dev, labels):
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
    ax.set_title("Output validity — LiftOn 2.0 vs v1.0.8",
                 fontsize=12, fontweight="bold", loc="left")
    return _save(fig, "fig5_validity.png")


# --------------------------------------------------------------------------- #
# FIGURE 3 — performance (A wall-clock, B peak RSS log)
# --------------------------------------------------------------------------- #
def fig_performance(vc):
    ctrl = sorted([r for r in vc.values() if r.get("arm") == "controlled"],
                  key=lambda r: r["benchmark"])
    fig = plt.figure(figsize=(11.5, 4.9))
    gs = fig.add_gridspec(1, 2, wspace=0.26)
    axA = fig.add_subplot(gs[0, 0])
    axB = fig.add_subplot(gs[0, 1])

    def grouped(ax, metric, ylabel, log, ratio_fmt):
        recs = [r for r in ctrl
                if isinstance((r.get(metric) or {}).get("stable"), (int, float))
                and isinstance((r.get(metric) or {}).get("devel"), (int, float))]
        labels = [r["benchmark"] for r in recs]
        sta = [r[metric]["stable"] for r in recs]
        dev = [r[metric]["devel"] for r in recs]
        x = np.arange(len(recs))
        w = 0.38
        ax.bar(x - w / 2, sta, w, color=mr.TOOL_COLORS["lifton_stable"],
               label="LiftOn v1.0.8")
        b2 = ax.bar(x + w / 2, dev, w, color=mr.TOOL_COLORS["lifton_devel"],
                    label="LiftOn 2.0")
        for rect, s, d in zip(b2, sta, dev):
            if d:
                ax.annotate(ratio_fmt(s, d),
                            (rect.get_x() + rect.get_width() / 2, rect.get_height()),
                            ha="center", va="bottom", fontsize=8, fontweight="bold",
                            color=mr.TOOL_COLORS["lifton_devel"])
        ax.set_xticks(x)
        ax.set_xticklabels(labels, rotation=20, ha="right", fontsize=8.5)
        ax.set_ylabel(ylabel)
        if log:
            ax.set_yscale("log")
        ax.legend(fontsize=8, loc="upper right")
        ax.grid(axis="y", alpha=0.3)
        ax.margins(y=0.18)

    grouped(axA, "wall_s", "wall-clock (s)", False,
            lambda s, d: f"{s/d:.2f}×")
    _panel_title(axA, "A", "Faster end-to-end (same inputs)")
    grouped(axB, "peak_rss_mb", "peak memory (MB, log)", True,
            lambda s, d: f"{s/d:.1f}×")
    _panel_title(axB, "B", "Far lighter on RAM (windowed aligner)")

    return _save(fig, "fig3_performance.png")


def main():
    vc, fw, _ = mr._load()
    if not fw or not vc:
        sys.exit("ERROR: source JSON not found / empty — run the comparison harness first.")
    print(f"OUTDIR = {OUTDIR}")
    fig_accuracy(fw, vc)
    fig_fourway_pi(fw)
    fig_completeness_bytool(fw)
    fig_completeness(fw)
    fig_validity(fw)
    fig_performance(vc)
    print("done.")


if __name__ == "__main__":
    main()
