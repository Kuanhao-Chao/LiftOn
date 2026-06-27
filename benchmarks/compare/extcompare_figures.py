#!/usr/bin/env python
"""extcompare_figures.py — figures for the broader-comparison sub-study.

Feeds on benchmarks/compare/extcompare_results.json (the NEW file written by
extcompare.py; the committed fourway_results.json is never touched). The
sub-study compares the three homology lift-over tools (Liftoff, LiftOn v1.0.8,
LiftOn v1.0.9) against coordinate **LiftOver** (CrossMap over a per-pair Cactus
chain) on a representative 5-pair chromosome-subset, with miniprot kept as
protein evidence. CAT (v1/v2) was attempted but its toil sub-workflows proved
intractable on the RefSeq subset — recorded, no column.

LiftOver is reported with a DUAL metric, the honest out-of-class framing:
  (A) re-scored protein identity — LOW, because coordinate lift does NOT correct
      frame; and (C) feature-lift rate — the fraction of reference coding
      transcripts whose CDS lifted at all (its FAIR metric, typically high).

Reuses master_report's palette/labels so colours/names never drift from the
report. Renders ONLY the pairs present in the JSON, so it is safe to run while
the remaining Cactus alignments are still recovering.

Usage (repo root):
  PYTHONNOUSERSITE=1 python benchmarks/compare/extcompare_figures.py [OUTDIR]
Writes PNGs into benchmarks/compare/figures/extcompare/ by default.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt           # noqa: E402
from matplotlib.patches import Patch       # noqa: E402
import numpy as np                         # noqa: E402

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import master_report as mr                 # noqa: E402

RESULTS = HERE / "extcompare_results.json"
OUTDIR = Path(sys.argv[1]) if len(sys.argv) > 1 else (HERE / "figures" / "extcompare")

# headline lift-over tools (bars) + LiftOver (out-of-class bar); miniprot overlaid
# as an evidence marker, never a headline bar.
BAR_TOOLS = ["liftoff", "lifton_stable", "lifton_devel", "liftover"]
EVIDENCE = "miniprot"

# canonical pair order: same-species -> close -> distant -> very-distant
PAIR_ORDER = ["t1_maize_b73_to_mo17", "t2_human_to_gorilla", "t3_dog_to_cat",
              "arabidopsis_to_rice", "t4_human_to_chicken"]
PAIR_LABEL = {
    "t1_maize_b73_to_mo17": "maize B73→Mo17\n(same-species)",
    "t2_human_to_gorilla": "human→gorilla\n(close)",
    "t3_dog_to_cat": "dog→cat\n(distant)",
    "arabidopsis_to_rice": "Arabidopsis→rice\n(very-distant)",
    "t4_human_to_chicken": "human→chicken\n(very-distant)",
}
LABEL = dict(mr.TOOL_LABEL)
LABEL["lifton_stable"] = "LiftOn v1.0.8"
LABEL["lifton_devel"] = "LiftOn v1.0.9"
LABEL["liftover"] = "LiftOver (coord-only)"
COLOR = dict(mr.TOOL_COLORS)


def _load():
    if not RESULTS.exists():
        sys.exit(f"no results file: {RESULTS}")
    data = json.loads(RESULTS.read_text())
    cells = {}
    for pid in PAIR_ORDER:
        c = data.get(f"extcompare:{pid}")
        if c:
            cells[pid] = c
    if not cells:
        sys.exit("extcompare_results.json has no recognised pair cells")
    return cells


def _grouped_bars(ax, cells, field, title, ylabel, ylim=None):
    """Grouped bars: one group per present pair, one bar per BAR_TOOL, with the
    miniprot evidence value drawn as a diamond marker over each group."""
    pids = [p for p in PAIR_ORDER if p in cells]
    n = len(BAR_TOOLS)
    width = 0.8 / n
    x = np.arange(len(pids))
    for i, tool in enumerate(BAR_TOOLS):
        vals = [cells[p][field].get(tool) for p in pids]
        xs = x + (i - (n - 1) / 2) * width
        hatch = "//" if tool == "liftover" else None
        bars = ax.bar(xs, [v if v is not None else 0 for v in vals], width,
                      label=LABEL.get(tool, tool), color=COLOR.get(tool, "#888"),
                      edgecolor="white", linewidth=0.4, hatch=hatch)
        for b, v in zip(bars, vals):
            if v is not None:
                ax.text(b.get_x() + b.get_width() / 2, v + 0.005, f"{v:.3f}",
                        ha="center", va="bottom", fontsize=5.5, rotation=90)
    # miniprot evidence marker
    ev = [cells[p][field].get(EVIDENCE) for p in pids]
    ax.scatter(x, [v if v is not None else np.nan for v in ev], marker="D", s=34,
               color=mr.TOOL_COLORS["miniprot"], edgecolor="black", linewidth=0.5,
               zorder=5, label="miniprot (protein evidence)")
    ax.set_xticks(x)
    ax.set_xticklabels([PAIR_LABEL[p] for p in pids], fontsize=7.5)
    ax.set_ylabel(ylabel, fontsize=9)
    ax.set_title(title, fontsize=10, fontweight="bold")
    if ylim:
        ax.set_ylim(*ylim)
    ax.grid(axis="y", alpha=0.25, linewidth=0.5)
    ax.set_axisbelow(True)


def _liftover_dual(ax, cells):
    """LiftOver's honest dual metric: feature-lift rate (high) vs re-scored
    protein identity (low, no frame correction)."""
    pids = [p for p in PAIR_ORDER if p in cells]
    x = np.arange(len(pids))
    w = 0.38
    lift = [(cells[p].get("liftover_feature_lift_rate") or {}).get("cds_lift_rate") for p in pids]
    pi = [cells[p]["mean_pi"].get("liftover") for p in pids]
    b1 = ax.bar(x - w / 2, [v or 0 for v in lift], w, color="#7aa6c2",
                edgecolor="white", label="feature-lift rate (fraction of ref CDS lifted)")
    b2 = ax.bar(x + w / 2, [v or 0 for v in pi], w, color="#9e9e9e", hatch="//",
                edgecolor="white", label="re-scored protein identity (no frame correction)")
    for bars, vals in ((b1, lift), (b2, pi)):
        for b, v in zip(bars, vals):
            if v is not None:
                ax.text(b.get_x() + b.get_width() / 2, v + 0.01, f"{v:.3f}",
                        ha="center", va="bottom", fontsize=6, rotation=90)
    ax.set_xticks(x)
    ax.set_xticklabels([PAIR_LABEL[p] for p in pids], fontsize=7.5)
    ax.set_ylabel("rate / identity", fontsize=9)
    ax.set_ylim(0, 1.08)
    ax.set_title("C. LiftOver dual metric — lifts most features, but coordinate-only\n"
                 "(no frame correction → re-translated CDS scores low)",
                 fontsize=10, fontweight="bold")
    ax.legend(fontsize=7.5, loc="lower left")
    ax.grid(axis="y", alpha=0.25, linewidth=0.5)
    ax.set_axisbelow(True)


def main():
    cells = _load()
    OUTDIR.mkdir(parents=True, exist_ok=True)
    n_pairs = sum(p in cells for p in PAIR_ORDER)

    fig, axes = plt.subplots(3, 1, figsize=(max(8, 1.9 * n_pairs + 2), 13))
    _grouped_bars(axes[0], cells, "mean_pi",
                  "A. Mean protein identity (tool-neutral parasail re-score)",
                  "mean protein identity", ylim=(0, 1.08))
    axes[0].legend(fontsize=7.5, ncol=3, loc="lower left")
    _grouped_bars(axes[1], cells, "completeness_coding",
                  "B. Coding completeness (recall of reference coding transcripts)",
                  "coding completeness", ylim=(0, 1.08))
    axes[1].legend(fontsize=7.5, ncol=3, loc="lower left")
    _liftover_dual(axes[2], cells)

    fig.suptitle(
        "Broader comparison: homology lift-over tools vs coordinate LiftOver\n"
        f"({n_pairs} tractable chromosome-subset pairs; miniprot = protein evidence)\n"
        "distant/very-distant pairs (dog→cat, Arabidopsis→rice, human→chicken) omitted — the "
        "Cactus alignment LiftOver+CAT require was intractable (>700 GB cons; see text)",
        fontsize=11, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    out = OUTDIR / "rfig_extcompare.png"
    fig.savefig(out, dpi=160, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out}  ({n_pairs} pairs: {', '.join(p for p in PAIR_ORDER if p in cells)})")


if __name__ == "__main__":
    main()
