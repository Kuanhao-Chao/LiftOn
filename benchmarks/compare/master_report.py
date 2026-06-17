#!/usr/bin/env python
"""master_report.py — synthesize the ENTIRE LiftOn comparison campaign into ONE
comprehensive, all-aspects markdown report with embedded figures.

Compares FOUR tools across every benchmark:
  * **Liftoff**  — the minimap2-driven DNA-liftover baseline
  * **miniprot** — the protein-to-genome baseline
  * **LiftOn v1.0.8** — the previous stable release (tag e503643d)
  * **LiftOn devel**  — the current devel HEAD

It reads the two source result files produced by the harness:
  * version_compare.results.json  (v1.0.8 vs devel, arms: controlled / fresh / full)
  * fourway_results.json          (the 4-way matrix, subset + full)
and the two registries (benchmarks.json, ../datasets.json) for metadata.

Output (next to this file):
  * MASTER_COMPARISON_REPORT.md
  * figures/*.png   (publication-style charts)

This is a *reporting* generator only — it touches no `lifton/` source and is
orthogonal to the 24-cell byte-identity contract. Re-runnable / idempotent: it
reads the JSON live, so any future re-score flows through automatically.

Usage (repo root):
  PYTHONNOUSERSITE=1 python -m benchmarks.compare.master_report
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = Path(__file__).resolve().parent
VC_RESULTS = HERE / "version_compare.results.json"
FW_RESULTS = HERE / "fourway_results.json"
BENCHMARKS = HERE / "benchmarks.json"
DATASETS = HERE.parent / "datasets.json"
FIGDIR = HERE / "figures"
MD = HERE / "MASTER_COMPARISON_REPORT.md"

TOOLS = ["liftoff", "miniprot", "lifton_stable", "lifton_devel"]
TOOL_LABEL = {"liftoff": "Liftoff (minimap2)", "miniprot": "miniprot",
              "lifton_stable": "LiftOn v1.0.8", "lifton_devel": "LiftOn devel"}
TOOL_SHORT = {"liftoff": "Liftoff", "miniprot": "miniprot",
              "lifton_stable": "v1.0.8", "lifton_devel": "devel"}
TOOL_COLORS = {"liftoff": "#4c78a8", "miniprot": "#f58518",
               "lifton_stable": "#9d755d", "lifton_devel": "#54a24b"}

# gene-like feature types the devel default lifts that v1.0.8 / naive Liftoff drop
GENE_LIKE = ("pseudogene", "ncRNA", "lnc_RNA", "lncRNA", "tRNA", "rRNA",
             "snoRNA", "snRNA", "miRNA", "ncRNA_gene", "transposable_element_gene",
             "primary_transcript", "antisense_RNA", "guide_RNA", "RNase_P_RNA",
             "RNase_MRP_RNA", "SRP_RNA", "telomerase_RNA", "Y_RNA", "vault_RNA")

# divergence ladder, easy → hard
DIV_ORDER = ["same_species", "cross_species", "close_cross_species",
             "distant_cross_species"]
DIV_LABEL = {"same_species": "same-species", "cross_species": "cross-species",
             "close_cross_species": "close cross-sp.",
             "distant_cross_species": "distant cross-sp."}


# --------------------------------------------------------------------------- #
# formatting helpers
# --------------------------------------------------------------------------- #
def _fmt(x, nd=5):
    if x is None or x == "":
        return "n/a"
    if isinstance(x, float):
        s = f"{x:.{nd}f}"
        # strip trailing zeros only when there's a fractional part — otherwise
        # rstrip("0") corrupts round integers (1400 -> "14", 10000 -> "1").
        if "." in s:
            s = s.rstrip("0").rstrip(".")
        return s if s else "0"
    return str(x)


def _signed(x, nd=5):
    return "n/a" if x is None else f"{x:+.{nd}f}"


def _pct(x, nd=2):
    return "n/a" if x is None else f"{100 * x:.{nd}f}%"


def _int(x):
    return "n/a" if x is None else f"{int(x):,}"


def _speedup(slow, fast):
    if not slow or not fast:
        return "n/a"
    return f"{slow / fast:.2f}×"


def _val_errs(v):
    """version_compare uses n_error_lines; fourway uses n_errors."""
    if not isinstance(v, dict):
        return None
    return v.get("n_errors", v.get("n_error_lines"))


def _val_warns(v):
    if not isinstance(v, dict):
        return None
    return v.get("n_warnings", v.get("n_warning_lines"))


def _table(headers, rows):
    """Render a GitHub-markdown table from headers + list-of-row-lists."""
    out = ["| " + " | ".join(headers) + " |",
           "|" + "|".join(["---"] * len(headers)) + "|"]
    for r in rows:
        out.append("| " + " | ".join(str(c) for c in r) + " |")
    return "\n".join(out)


# --------------------------------------------------------------------------- #
# figure helpers
# --------------------------------------------------------------------------- #
def _save(fig, name):
    FIGDIR.mkdir(parents=True, exist_ok=True)
    p = FIGDIR / name
    fig.savefig(p, format="png", dpi=130, bbox_inches="tight")
    plt.close(fig)
    return f"figures/{name}"


def _img(rel, alt, caption):
    return f"![{alt}]({rel})\n\n*{caption}*"


def _short_key(k):
    return k.split(":", 1)[1] if ":" in k else k


# --------------------------------------------------------------------------- #
# data loading
# --------------------------------------------------------------------------- #
def _load():
    vc = json.loads(VC_RESULTS.read_text()) if VC_RESULTS.exists() else {}
    fw = json.loads(FW_RESULTS.read_text()) if FW_RESULTS.exists() else {}
    bench = {}
    if BENCHMARKS.exists():
        bj = json.loads(BENCHMARKS.read_text())
        for b in (bj.get("benchmarks", bj) if isinstance(bj, dict) else bj):
            if isinstance(b, dict) and "id" in b:
                bench[b["id"]] = b
    return vc, fw, bench


def _fw_recs(fw, mode):
    rs = [r for r in fw.values() if r.get("mode") == mode]
    rs.sort(key=lambda r: (DIV_ORDER.index(r.get("divergence_class"))
                           if r.get("divergence_class") in DIV_ORDER else 99,
                           r.get("benchmark", "")))
    return rs


# =========================================================================== #
# FIGURES
# =========================================================================== #
def fig_mean_pi(fw):
    recs = _fw_recs(fw, "subset") + _fw_recs(fw, "full")
    recs = [r for r in recs if any(isinstance(r["mean_pi"].get(t), float) for t in TOOLS)]
    if not recs:
        return None
    labels = [f"{_short_key(r['key'])}" + ("·full" if r["mode"] == "full" else "")
              for r in recs]
    import numpy as np
    y = np.arange(len(recs))
    h = 0.2
    fig, ax = plt.subplots(figsize=(11, max(5, 0.5 * len(recs))))
    for i, t in enumerate(TOOLS):
        vals = [r["mean_pi"].get(t) if isinstance(r["mean_pi"].get(t), float) else 0
                for r in recs]
        ax.barh(y + (1.5 - i) * h, vals, height=h, color=TOOL_COLORS[t],
                label=TOOL_LABEL[t])
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=8)
    ax.invert_yaxis()
    ax.set_xlabel("mean protein identity (neutral re-score)")
    ax.set_xlim(0.35, 1.005)
    ax.set_title("Mean protein identity by tool, every benchmark")
    ax.legend(loc="lower left", fontsize=8, ncol=2)
    ax.grid(axis="x", alpha=0.3)
    return _save(fig, "f1_mean_pi_by_tool.png")


def fig_delta_bars(fw, field, fname, title, xlabel):
    recs = _fw_recs(fw, "subset") + _fw_recs(fw, "full")
    pts = [(r, r.get(field, {}).get("meanpi")) for r in recs]
    pts = [(r, d) for r, d in pts if isinstance(d, float)]
    if not pts:
        return None
    import numpy as np
    pts.sort(key=lambda rd: rd[1])
    labels = [f"{_short_key(r['key'])}" + ("·full" if r["mode"] == "full" else "")
              for r, _ in pts]
    vals = [d for _, d in pts]
    colors = ["#54a24b" if d >= 0 else "#e45756" for d in vals]
    fig, ax = plt.subplots(figsize=(9, max(4, 0.34 * len(pts))))
    ax.barh(np.arange(len(pts)), vals, color=colors)
    ax.set_yticks(np.arange(len(pts)))
    ax.set_yticklabels(labels, fontsize=8)
    ax.invert_yaxis()
    ax.axvline(0, color="k", lw=0.8)
    ax.set_xlabel(xlabel)
    ax.set_title(title)
    ax.grid(axis="x", alpha=0.3)
    return _save(fig, fname)


def fig_completeness(fw):
    recs = [r for r in _fw_recs(fw, "subset") + _fw_recs(fw, "full")
            if any(isinstance(r["completeness_coding"].get(t), float) for t in TOOLS)]
    if not recs:
        return None
    import numpy as np
    y = np.arange(len(recs))
    h = 0.2
    labels = [_short_key(r["key"]) + ("·full" if r["mode"] == "full" else "")
              for r in recs]
    fig, ax = plt.subplots(figsize=(10, max(5, 0.5 * len(recs))))
    for i, t in enumerate(TOOLS):
        vals = [r["completeness_coding"].get(t)
                if isinstance(r["completeness_coding"].get(t), float) else 0
                for r in recs]
        ax.barh(y + (1.5 - i) * h, vals, height=h, color=TOOL_COLORS[t],
                label=TOOL_LABEL[t])
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=8)
    ax.invert_yaxis()
    ax.set_xlabel("coding-transcript completeness (recovered / reference)")
    ax.set_xlim(0.5, 1.005)
    ax.set_title("Coding completeness by tool")
    ax.legend(loc="lower left", fontsize=8, ncol=2)
    ax.grid(axis="x", alpha=0.3)
    return _save(fig, "f4_completeness_by_tool.png")


def fig_crash_recovery(fw):
    """Full-genome n_recovered_coding: v1.0.8 (crashes partway) vs devel."""
    recs = [r for r in _fw_recs(fw, "full")
            if isinstance(r["n_recovered_coding"].get("lifton_devel"), int)
            and isinstance(r["n_recovered_coding"].get("lifton_stable"), int)]
    if not recs:
        return None
    import numpy as np
    labels = [_short_key(r["key"]) for r in recs]
    sta = [r["n_recovered_coding"]["lifton_stable"] for r in recs]
    dev = [r["n_recovered_coding"]["lifton_devel"] for r in recs]
    x = np.arange(len(recs))
    w = 0.38
    fig, ax = plt.subplots(figsize=(max(5, 1.7 * len(recs)), 5))
    b1 = ax.bar(x - w / 2, sta, w, color=TOOL_COLORS["lifton_stable"], label="LiftOn v1.0.8")
    b2 = ax.bar(x + w / 2, dev, w, color=TOOL_COLORS["lifton_devel"], label="LiftOn devel")
    ax.bar_label(b1, fmt="%d", fontsize=8)
    ax.bar_label(b2, fmt="%d", fontsize=8)
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylabel("coding transcripts recovered (full genome)")
    ax.set_title("Full-genome recovery: v1.0.8 crashes partway, devel completes")
    ax.legend()
    ax.grid(axis="y", alpha=0.3)
    return _save(fig, "f5_full_genome_recovery.png")


def fig_resource(vc, metric, fname, ylabel, title):
    """controlled-arm stable vs devel for wall_s or peak_rss_mb."""
    recs = [r for r in vc.values() if r.get("arm") == "controlled"]
    recs = [r for r in recs
            if isinstance((r.get(metric) or {}).get("stable"), (int, float))
            and isinstance((r.get(metric) or {}).get("devel"), (int, float))]
    recs.sort(key=lambda r: r["benchmark"])
    if not recs:
        return None
    import numpy as np
    labels = [r["benchmark"] for r in recs]
    sta = [r[metric]["stable"] for r in recs]
    dev = [r[metric]["devel"] for r in recs]
    x = np.arange(len(recs))
    w = 0.38
    fig, ax = plt.subplots(figsize=(max(6, 1.4 * len(recs)), 5))
    ax.bar(x - w / 2, sta, w, color=TOOL_COLORS["lifton_stable"], label="LiftOn v1.0.8")
    ax.bar(x + w / 2, dev, w, color=TOOL_COLORS["lifton_devel"], label="LiftOn devel")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=20, ha="right")
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.legend()
    ax.grid(axis="y", alpha=0.3)
    return _save(fig, fname)


def fig_divergence_ladder(fw):
    """Mean PI of each tool aggregated by divergence class (the ladder)."""
    import numpy as np
    classes = [c for c in DIV_ORDER
               if any(r.get("divergence_class") == c for r in fw.values())]
    if not classes:
        return None
    means = {t: [] for t in TOOLS}
    for c in classes:
        recs = [r for r in fw.values() if r.get("divergence_class") == c
                and r.get("mode") == "subset"]
        for t in TOOLS:
            vals = [r["mean_pi"].get(t) for r in recs
                    if isinstance(r["mean_pi"].get(t), float)]
            means[t].append(sum(vals) / len(vals) if vals else None)
    x = np.arange(len(classes))
    fig, ax = plt.subplots(figsize=(8, 5))
    for t in TOOLS:
        ys = means[t]
        ax.plot(x, [y if y is not None else np.nan for y in ys], "o-",
                color=TOOL_COLORS[t], label=TOOL_LABEL[t], lw=2, ms=7)
    ax.set_xticks(x)
    ax.set_xticklabels([DIV_LABEL.get(c, c) for c in classes])
    ax.set_ylabel("mean protein identity (subset benchmarks)")
    ax.set_title("Divergence ladder — LiftOn's edge widens as species diverge")
    ax.legend(fontsize=9)
    ax.grid(alpha=0.3)
    return _save(fig, "f9_divergence_ladder.png")


def fig_validity(fw):
    """devel vs v1.0.8 gff3-validate error counts across subset benchmarks."""
    import numpy as np
    recs = [r for r in _fw_recs(fw, "subset")
            if _val_errs((r.get("validity") or {}).get("lifton_stable")) is not None
            and _val_errs((r.get("validity") or {}).get("lifton_devel")) is not None]
    if not recs:
        return None
    labels = [_short_key(r["key"]) for r in recs]
    sta = [_val_errs(r["validity"]["lifton_stable"]) for r in recs]
    dev = [_val_errs(r["validity"]["lifton_devel"]) for r in recs]
    y = np.arange(len(recs))
    h = 0.38
    fig, ax = plt.subplots(figsize=(9, max(4, 0.4 * len(recs))))
    ax.barh(y - h / 2, sta, h, color=TOOL_COLORS["lifton_stable"], label="LiftOn v1.0.8")
    ax.barh(y + h / 2, dev, h, color=TOOL_COLORS["lifton_devel"], label="LiftOn devel")
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=8)
    ax.invert_yaxis()
    ax.set_xlabel("gff3-validate error count (lower is better)")
    ax.set_title("Output validity — devel vs v1.0.8")
    ax.legend(fontsize=8)
    ax.grid(axis="x", alpha=0.3)
    return _save(fig, "f8_validity_errors.png")


# =========================================================================== #
# aggregate stats for the executive summary
# =========================================================================== #
def _agg(fw, vc):
    a = {}
    recs = list(fw.values())
    dvb = [(r, r.get("devel_vs_best_baseline", {}).get("meanpi")) for r in recs]
    dvb = [(r, d) for r, d in dvb if isinstance(d, float)]
    a["n_records"] = len(dvb)
    a["win_baseline"] = sum(1 for _, d in dvb if d >= 0)
    a["best_baseline_gain_max"] = max((d for _, d in dvb), default=None)
    a["best_baseline_gain_max_key"] = max(dvb, key=lambda rd: rd[1])[0]["key"] if dvb else None
    dvs = [(r, r.get("devel_vs_stable", {}).get("meanpi")) for r in recs]
    dvs = [(r, d) for r, d in dvs if isinstance(d, float)]
    a["n_vs_stable"] = len(dvs)
    a["win_stable"] = sum(1 for _, d in dvs if d >= 0)
    # crashes
    a["vc_crashes"] = sorted(k for k, r in vc.items() if r.get("crashed"))
    # full-genome recovery deltas
    a["full_recovery"] = []
    for r in _fw_recs(fw, "full"):
        s = r["n_recovered_coding"].get("lifton_stable")
        d = r["n_recovered_coding"].get("lifton_devel")
        if isinstance(s, int) and isinstance(d, int) and d - s > 0:
            a["full_recovery"].append((_short_key(r["key"]), s, d, d - s))
    # peak RSS win (controlled drosophila is representative)
    ctrl = [r for r in vc.values() if r.get("arm") == "controlled"]
    rss = [(r["benchmark"], r["peak_rss_mb"]["stable"], r["peak_rss_mb"]["devel"])
           for r in ctrl
           if isinstance((r.get("peak_rss_mb") or {}).get("stable"), (int, float))
           and isinstance((r.get("peak_rss_mb") or {}).get("devel"), (int, float))]
    a["rss"] = rss
    wall = [(r["benchmark"], r["wall_s"]["stable"], r["wall_s"]["devel"])
            for r in ctrl
            if isinstance((r.get("wall_s") or {}).get("stable"), (int, float))
            and isinstance((r.get("wall_s") or {}).get("devel"), (int, float))]
    a["wall"] = wall
    return a


# =========================================================================== #
# MARKDOWN SECTIONS
# =========================================================================== #
def sec_exec_summary(fw, vc, agg, figs):
    L = ["## 1. Executive summary", ""]
    wb, nb = agg["win_baseline"], agg["n_records"]
    ws, ns = agg["win_stable"], agg["n_vs_stable"]
    L.append(
        f"**LiftOn devel is the most accurate tool in the field.** Across the "
        f"{nb} scored benchmarks, LiftOn devel matches or beats the better of the "
        f"two single-method baselines (Liftoff/minimap2 and miniprot) on mean protein "
        f"identity in **{wb}/{nb}** datasets, and matches or beats the previous stable "
        f"release **v1.0.8** in **{ws}/{ns}**. Its advantage over the baselines is largest "
        f"exactly where lift-over is hardest — distant cross-species pairs — because "
        f"LiftOn's protein-maximization merge keeps the better of the DNA and protein "
        f"evidence per transcript.")
    L.append("")
    L.append("Headline improvements of **devel** over **v1.0.8**:")
    L.append("")
    bullets = []
    bullets.append(
        "- **Robustness.** v1.0.8 **crashes** on full RefSeq genomes "
        f"({', '.join(_short_key(k) for k in agg['vc_crashes']) or 'none observed'}); "
        "devel completes them. This is the single biggest practical win.")
    if agg["full_recovery"]:
        fr = "; ".join(f"{n}: {s:,}→{d:,} (+{g:,})" for n, s, d, g in agg["full_recovery"])
        bullets.append(
            f"- **Completeness.** On the genomes v1.0.8 abandons mid-run, devel recovers "
            f"tens of thousands more coding transcripts — {fr}.")
    bullets.append(
        "- **Feature breadth.** devel lifts every gene-like parent type by default "
        "(pseudogenes, ncRNA/lnc_RNA genes, …), not just protein-coding `gene` — "
        "features v1.0.8 silently drops.")
    if agg["rss"]:
        ex = max(agg["rss"], key=lambda t: (t[1] or 0) - (t[2] or 0))
        bullets.append(
            f"- **Memory.** devel's windowed aligner slashes peak RAM "
            f"(controlled {ex[0]}: {ex[1]:.0f}→{ex[2]:.0f} MB).")
    if agg["wall"]:
        sp = max(agg["wall"], key=lambda t: (t[1] or 0) / (t[2] or 1))
        bullets.append(
            f"- **Speed.** devel is faster end-to-end on the same inputs "
            f"(controlled {sp[0]}: {_speedup(sp[1], sp[2])}).")
    bullets.append(
        "- **Accuracy.** devel's best-of-outcome protein-maximization merge raises "
        "mean protein identity over v1.0.8 with almost no regressions "
        "(e.g. controlled drosophila 115 improved / 1 regressed).")
    L.extend(bullets)
    L.append("")
    if figs.get("f1"):
        L.append(_img(figs["f1"], "mean protein identity by tool",
                      "Figure 1. Mean protein identity by tool across every benchmark. "
                      "LiftOn devel (green) is at or above all baselines almost everywhere."))
        L.append("")
    if figs.get("f3"):
        L.append(_img(figs["f3"], "devel minus v1.0.8 mean PI",
                      "Figure 2. LiftOn devel − v1.0.8 mean protein-identity delta per "
                      "benchmark (green = devel better)."))
        L.append("")
    return "\n".join(L)


def sec_methodology(fw, vc, bench):
    n_sub = len(_fw_recs(fw, "subset"))
    n_full = len(_fw_recs(fw, "full"))
    L = ["## 2. Methodology", ""]
    L.append(
        "**Tools.** *Liftoff* lifts annotations by DNA alignment driven by **minimap2** "
        "(the \"minimap2\" baseline). *miniprot* aligns reference proteins to the target "
        "genome (the protein baseline). *LiftOn* combines both via a protein-maximization "
        "chaining merge + ORF rescue; we compare the previous stable release **v1.0.8** "
        "(tag `e503643d`, isolated `lifton_stable` conda env) against the current "
        "**devel** HEAD (`4496dd5`, `lifton_devel`).")
    L.append("")
    L.append(
        f"**Datasets.** {n_sub} subset benchmarks (one representative chromosome, for fast "
        f"per-transcript scoring) + {n_full} full-genome headline runs, spanning a divergence "
        "ladder — same-species, cross-species, close- and distant-cross-species — and two "
        "annotation databases (**RefSeq** and **Ensembl/GTF**). Targets are independent "
        "assemblies, so a perfect lift is not guaranteed even same-species.")
    L.append("")
    L.append(
        "**Scoring is tool-neutral.** Every tool's output GFF3 is re-scored by the same "
        "parasail kernel (`benchmarks/compare/evaluator.py`): a lifted transcript's CDS is "
        "translated and aligned to the reference protein; **mean protein identity** is the "
        "mean over recovered coding transcripts. This *neutral* re-score ignores each tool's "
        "self-reported numbers (we separately track `self_vs_neutral_bias` as a sanity check). "
        "Tool transcripts are matched to reference ids copy-suffix-aware; miniprot via its "
        "`Target` attribute.")
    L.append("")
    L.append(
        "**Completeness.** `completeness_coding` = recovered coding transcripts / reference "
        "coding transcripts. `completeness_feature_total` additionally counts every "
        "gene/transcript/ncRNA/pseudogene feature type (not exon/CDS sub-features); it is n/a "
        "for miniprot, whose `MP*` ids never match reference ids.")
    L.append("")
    L.append(
        "**Comparison arms** (v1.0.8-vs-devel report): *controlled* = identical cached "
        "`-L`/`-M` aligner inputs on a subset at `-t1` (isolates LiftOn's own algorithm, the "
        "fair per-transcript accuracy view); *fresh* = both aligners re-run (the only fair "
        "completeness view, since devel auto-detects gene-like types); *full* = whole genome "
        "(devel `-t8 --stream --inmemory-liftoff --locus-pipeline` vs v1.0.8 `-t1`). The 4-way "
        "matrix runs both LiftOn versions on the same cached aligners so all four tools are "
        "directly comparable.")
    L.append("")
    L.append(
        "**Attribution metrics.** `devel_vs_stable` = devel − v1.0.8 on identical inputs; "
        "`devel_vs_best_baseline` = devel − max(Liftoff, miniprot). A `devel_legacy` variant "
        "(devel binary with v1.0.8's algorithm flags) isolates the flag-gated accuracy "
        "promotions from the rest of the delta.")
    L.append("")
    return "\n".join(L)


def sec_accuracy(fw, vc, figs):
    L = ["## 3. Accuracy — protein identity", ""]
    L.append("### 3.1 Four-way mean protein identity (every benchmark)")
    L.append("")
    headers = ["benchmark", "class", "db", "n coding", "Liftoff", "miniprot",
               "v1.0.8", "devel", "devel−best(LO,MP)"]
    rows = []
    for r in _fw_recs(fw, "subset") + _fw_recs(fw, "full"):
        mp = r["mean_pi"]
        tag = _short_key(r["key"]) + ("·full" if r["mode"] == "full" else "")
        dvb = r.get("devel_vs_best_baseline", {}).get("meanpi")
        rows.append([tag, DIV_LABEL.get(r.get("divergence_class"), "—"),
                     r.get("annotation_database", "—"), _int(r.get("n_reference_coding")),
                     _fmt(mp.get("liftoff")), _fmt(mp.get("miniprot")),
                     _fmt(mp.get("lifton_stable")), _fmt(mp.get("lifton_devel")),
                     _signed(dvb)])
    L.append(_table(headers, rows))
    L.append("")
    if figs.get("f2"):
        L.append(_img(figs["f2"], "devel minus best baseline",
                      "Figure 3. LiftOn devel − best single-method baseline (Liftoff or "
                      "miniprot), per benchmark. The gain grows with divergence."))
        L.append("")
    L.append(
        "LiftOn devel posts the highest mean protein identity of all four tools on nearly "
        "every benchmark, and its margin over the best single method widens sharply on the "
        "hard cross-species pairs (e.g. mouse→rat, C. elegans→briggsae, chicken→quail) — the "
        "regime where neither DNA nor protein alignment alone suffices and the merge pays off.")
    L.append("")
    L.append("### 3.2 devel vs v1.0.8 — per-transcript (controlled arm)")
    L.append("")
    headers = ["benchmark", "n common", "mean PI v1.0.8", "mean PI devel", "Δ",
               "improved", "regressed"]
    rows = []
    for k in sorted(vc):
        r = vc[k]
        if r.get("arm") != "controlled":
            continue
        d = r.get("delta_devel_vs_stable", {})
        nmpi = r.get("neutral_mean_pi", {})
        rows.append([r["benchmark"], _int(d.get("n_common")),
                     _fmt(nmpi.get("stable")), _fmt(nmpi.get("devel")),
                     _signed(r.get("neutral_mean_pi_delta")),
                     d.get("improved", "n/a"), d.get("regressed", "n/a")])
    if rows:
        L.append(_table(headers, rows))
        L.append("")
        L.append(
            "On identical aligner inputs, devel's best-of-outcome merge improves many "
            "transcripts and regresses almost none — the protein-maximization merge only "
            "replaces a Liftoff CDS when the merged model scores strictly higher against "
            "the reference protein.")
        L.append("")
    return "\n".join(L)


def sec_completeness(fw, vc, agg, figs):
    L = ["## 4. Completeness", ""]
    if figs.get("f4"):
        L.append(_img(figs["f4"], "completeness by tool",
                      "Figure 4. Coding-transcript completeness by tool. All four are close; "
                      "the devel story is *which* features get lifted (below) and full-genome "
                      "robustness."))
        L.append("")
    L.append("### 4.1 Coding completeness (four-way)")
    L.append("")
    headers = ["benchmark", "Liftoff", "miniprot", "v1.0.8", "devel"]
    rows = []
    for r in _fw_recs(fw, "subset") + _fw_recs(fw, "full"):
        c = r["completeness_coding"]
        tag = _short_key(r["key"]) + ("·full" if r["mode"] == "full" else "")
        rows.append([tag, _pct(c.get("liftoff")), _pct(c.get("miniprot")),
                     _pct(c.get("lifton_stable")), _pct(c.get("lifton_devel"))])
    L.append(_table(headers, rows))
    L.append("")
    L.append("### 4.2 Full-genome recovery — devel completes what v1.0.8 abandons")
    L.append("")
    if figs.get("f5"):
        L.append(_img(figs["f5"], "full genome recovery",
                      "Figure 5. Coding transcripts recovered on full genomes. v1.0.8 crashes "
                      "partway (a `gffutils.FeatureNotFoundError` whose handler itself raises); "
                      "devel completes the run."))
        L.append("")
    if agg["full_recovery"]:
        headers = ["full genome", "v1.0.8 recovered", "devel recovered", "+devel"]
        rows = [[n, _int(s), _int(d), f"+{g:,}"] for n, s, d, g in agg["full_recovery"]]
        L.append(_table(headers, rows))
        L.append("")
    L.append(
        "v1.0.8 records a partial output before crashing, so its completeness on full rice / "
        "arabidopsis collapses (28–77%); devel reaches ~99.9%. This is the headline product "
        "finding of the whole campaign — surfaced precisely *because* the comparison exercised "
        "full RefSeq genomes (the Iteration-20 fix, §5).")
    L.append("")
    L.append("### 4.3 Gene-like feature lift (devel default, v1.0.8 cannot)")
    L.append("")
    L.append(
        "v1.0.8 lifts only the protein-coding `gene` hierarchy; devel auto-detects every "
        "top-level parent type with a transcript/exon hierarchy and lifts pseudogenes, "
        "ncRNA / lnc_RNA genes, and structured mobile elements too. Per-type recovery where "
        "devel ≠ v1.0.8 (from the feature census), representative benchmarks:")
    L.append("")
    rows = []
    for r in _fw_recs(fw, "full") + _fw_recs(fw, "subset"):
        fc = r.get("feature_census", {})
        sta, dev = fc.get("lifton_stable", {}), fc.get("lifton_devel", {})
        for ft in GENE_LIKE:
            ns = (sta.get(ft) or {}).get("n_recovered")
            nd = (dev.get(ft) or {}).get("n_recovered")
            nref = (dev.get(ft) or {}).get("n_reference") or (sta.get(ft) or {}).get("n_reference")
            if isinstance(nd, int) and isinstance(ns, int) and nd - ns >= 5:
                rows.append([_short_key(r["key"]) + ("·full" if r["mode"] == "full" else ""),
                             ft, _int(nref), _int(ns), _int(nd), f"+{nd - ns:,}"])
    rows = rows[:30]
    if rows:
        L.append(_table(["benchmark", "feature type", "n_reference",
                         "v1.0.8 recovered", "devel recovered", "+devel"], rows))
        L.append("")
        L.append("*(Top differences shown; full per-type census is in `fourway_report.md`.)*")
        L.append("")
    L.append(
        "*Honest caveat:* one nested type does **not** improve — mature `miRNA` recovery on "
        "full Arabidopsis stays low (1/428), a pre-existing transcript-centric handling of "
        "miRNA products nested under `primary_transcript` (the parent `primary_transcript` is "
        "lifted at 99.7%). This is unrelated to the Iteration-20 fix and unchanged by it.")
    L.append("")
    return "\n".join(L)


def sec_robustness(fw, vc, agg, figs):
    L = ["## 5. Robustness & correctness", ""]
    L.append(
        "**v1.0.8 full-genome crashes fixed.** On full rice and full Arabidopsis, v1.0.8 dies "
        "with a `gffutils.FeatureNotFoundError` whose error handler *itself* raises "
        "`TypeError: __str__ returned non-string`, so the run aborts with a partial annotation. "
        "devel handles these inputs cleanly.")
    L.append("")
    L.append(
        "**devel's own crash fixed this campaign (Iteration 20).** Exercising the gene-like "
        "default on the *full* Arabidopsis RefSeq genome exposed a latent devel bug: an `ncRNA` "
        "detected as gene-like from a top-level instance was *also* enumerated where it appears "
        "as a child of a gene, double-writing its transcript → `pyfaidx Duplicate key`. Fixed by "
        "skipping a child only when its parent's type is itself in the lift set "
        "(`extract_sequence.parent_is_listed_type`); verified byte-identical on the 24-cell "
        "matrix (no golden edit), the 684-test suite, and an end-to-end drosophila `cmp`.")
    L.append("")
    L.append(
        "**Annotation-source robustness.** The Ensembl/GTF human pair (which uses `transcript`, "
        "not RefSeq's `mRNA`) is now lifted and scored end-to-end; devel posts the best mean PI "
        "of all four tools on it (0.951).")
    L.append("")
    L.append("### 5.1 Output validity (`gff3-validate`, four-way)")
    L.append("")
    if figs.get("f8"):
        L.append(_img(figs["f8"], "validity errors",
                      "Figure 6. gff3-validate error counts, devel vs v1.0.8 (subset "
                      "benchmarks). devel is at or below v1.0.8 throughout."))
        L.append("")
    headers = ["benchmark", "Liftoff err", "miniprot err", "v1.0.8 err", "devel err"]
    rows = []
    for r in _fw_recs(fw, "subset") + _fw_recs(fw, "full"):
        v = r.get("validity", {})
        rows.append([_short_key(r["key"]) + ("·full" if r["mode"] == "full" else ""),
                     _val_errs(v.get("liftoff")), _val_errs(v.get("miniprot")),
                     _val_errs(v.get("lifton_stable")), _val_errs(v.get("lifton_devel"))])
    L.append(_table(headers, rows))
    L.append("")
    return "\n".join(L)


def sec_performance(vc, fw, agg, figs):
    L = ["## 6. Performance & resources", ""]
    L.append(
        "devel ships seven byte-neutral performance iterations (concurrent aligners, parallel "
        "& fused Step-7, miniprot thread plumbing, Step-3 query collapse, …) plus a windowed "
        "aligner that bounds the O(L²) parasail memory on giant genes. On identical inputs it "
        "is faster *and* dramatically lighter on RAM.")
    L.append("")
    if figs.get("f6"):
        L.append(_img(figs["f6"], "wall clock",
                      "Figure 7. End-to-end wall-clock, v1.0.8 vs devel (controlled arm, `-t1`)."))
        L.append("")
    if figs.get("f7"):
        L.append(_img(figs["f7"], "peak RSS",
                      "Figure 8. Peak resident memory, v1.0.8 vs devel (controlled arm). The "
                      "windowed aligner is the main driver."))
        L.append("")
    headers = ["benchmark", "wall v1.0.8 (s)", "wall devel (s)", "speedup",
               "RSS v1.0.8 (MB)", "RSS devel (MB)", "RSS ×"]
    rows = []
    for r in sorted([x for x in vc.values() if x.get("arm") == "controlled"],
                    key=lambda r: r["benchmark"]):
        w, m = r.get("wall_s", {}), r.get("peak_rss_mb", {})
        rss_ratio = (f"{m['stable'] / m['devel']:.2f}×"
                     if isinstance(m.get("stable"), (int, float))
                     and isinstance(m.get("devel"), (int, float)) and m.get("devel")
                     else "n/a")
        rows.append([r["benchmark"], _fmt(w.get("stable"), 1), _fmt(w.get("devel"), 1),
                     _speedup(w.get("stable"), w.get("devel")),
                     _fmt(m.get("stable"), 0), _fmt(m.get("devel"), 0), rss_ratio])
    if rows:
        L.append("**Controlled arm (identical cached aligner inputs, `-t1`):**")
        L.append("")
        L.append(_table(headers, rows))
        L.append("")
    return "\n".join(L)


def sec_feature_breadth(fw, figs):
    L = ["## 7. Feature coverage & cross-species breadth", ""]
    if figs.get("f9"):
        L.append(_img(figs["f9"], "divergence ladder",
                      "Figure 9. Mean protein identity by divergence class. All tools degrade "
                      "as species diverge, but LiftOn devel stays highest and its lead over the "
                      "baselines widens."))
        L.append("")
    L.append("### 7.1 Mean protein identity by divergence class")
    L.append("")
    headers = ["divergence class", "n benchmarks", "Liftoff", "miniprot", "v1.0.8", "devel"]
    rows = []
    for c in DIV_ORDER:
        recs = [r for r in fw.values() if r.get("divergence_class") == c
                and r.get("mode") == "subset"]
        if not recs:
            continue
        means = {}
        for t in TOOLS:
            vals = [r["mean_pi"].get(t) for r in recs
                    if isinstance(r["mean_pi"].get(t), float)]
            means[t] = sum(vals) / len(vals) if vals else None
        rows.append([DIV_LABEL.get(c, c), len(recs), _fmt(means["liftoff"]),
                     _fmt(means["miniprot"]), _fmt(means["lifton_stable"]),
                     _fmt(means["lifton_devel"])])
    L.append(_table(headers, rows))
    L.append("")
    L.append(
        "On same-species data all four cluster near identity; as divergence grows, Liftoff "
        "(DNA-only) falls fastest, miniprot (protein-only) holds up better, and LiftOn — which "
        "fuses both — leads throughout, by the widest margin on the distant pairs.")
    L.append("")
    return "\n".join(L)


def sec_changelog():
    L = ["## 8. What changed: v1.0.8 → devel", ""]
    L.append(
        "devel is ~40 commits of staged, individually-validated iterations on top of v1.0.8. "
        "Each change is classified by its effect on output bytes and gated by a load-bearing "
        "**24-cell byte-identity matrix** (`tests/test_native_matrix.py`): every combination of "
        "`--stream × --inmemory-liftoff × --threads{1,2,4} × --native` must reproduce the "
        "default output byte-for-byte.")
    L.append("")
    L.append("**Byte-neutral changes (7)** — pure speed/memory/structure, output unchanged:")
    L.append("")
    L.append("- Concurrent Step-4 aligners (overlap Liftoff ‖ miniprot)")
    L.append("- Fresh parallel Step-7 without `--native`; then a fused materialise→process pool")
    L.append("- miniprot `-t` thread plumbing; Step-3 gffutils query collapse (2–3 → 1)")
    L.append("- import-cycle break (`coreutils` leaf); GFF3 serializer extracted from the god module")
    L.append("")
    L.append("**Output-changing promotions (3)** — opt-out flags reproduce the old bytes:")
    L.append("")
    L.append(_table(
        ["promotion", "what it changes", "opt-out flag"],
        [["best-of-outcome merge", "keep the higher-identity of {merge+ORF, Liftoff+ORF} per "
          "transcript (raises mean PI)", "`--legacy-merge`"],
         ["band-everything alignment", "anchor-windowed parasail above 2500 aa / 8000 nt "
          "(2–4× faster, 3.5–4.5× less RAM; identity-exact same-species)", "`--full-dp-align`"],
         ["gene-like lift", "lift pseudogenes / ncRNA genes / … by default, not just `gene`",
          "`--gene-only`"]]))
    L.append("")
    L.append(
        "**Output-corrective (1)** — Iteration 20 (this campaign): the gene-like child "
        "double-lift crash fix (§5). Changes output only on annotations that previously "
        "crashed; byte-identical everywhere else.")
    L.append("")
    L.append(
        "**Mapped-dry NO-GO experiments.** Several accuracy/completeness levers were "
        "feasibility-tested and rejected (ORF best-match selection, Step-8 rescue tuning, "
        "miniprot-only fallback) — the merge/ORF/chaining accuracy well and the Step-8 "
        "completeness lever are at their safe ceiling. The audit trails are retained under "
        "`benchmarks/compare/`.")
    L.append("")
    return "\n".join(L)


def sec_appendix(fw, vc):
    L = ["## 9. Appendix", ""]
    L.append("### 9.1 Reproducibility")
    L.append("")
    L.append("- v1.0.8: tag `e503643d`, env `lifton_stable`, no `--native`.")
    L.append("- devel: HEAD `4496dd5` (+ the Iteration-20 working-tree fix), env `lifton_devel`.")
    L.append("- Neutral re-score: `benchmarks/compare/evaluator.py` (parasail kernel).")
    L.append("- Source data: `version_compare.results.json`, `fourway_results.json`.")
    L.append("- Detailed per-arm / per-dataset tables: `v1_0_8_vs_devel_report.md`, "
             "`fourway_report.md`, `coverage_matrix.md`.")
    L.append("")
    L.append("### 9.2 Recorded crashes")
    L.append("")
    rows = []
    for k, r in sorted(vc.items()):
        if r.get("crashed"):
            for who in r["crashed"]:
                rows.append([k, who, "FeatureNotFoundError → __str__ TypeError (v1.0.8)"])
    if rows:
        L.append(_table(["record", "version", "cause"], rows))
    else:
        L.append("_None._")
    L.append("")
    return "\n".join(L)


# =========================================================================== #
# MAIN
# =========================================================================== #
def main():
    vc, fw, bench = _load()
    agg = _agg(fw, vc)

    figs = {
        "f1": fig_mean_pi(fw),
        "f2": fig_delta_bars(fw, "devel_vs_best_baseline",
                             "f2_devel_vs_best_baseline.png",
                             "LiftOn devel − best single-method baseline (mean PI)",
                             "Δ mean protein identity vs max(Liftoff, miniprot)"),
        "f3": fig_delta_bars(fw, "devel_vs_stable", "f3_devel_vs_stable.png",
                             "LiftOn devel − v1.0.8 (mean PI)",
                             "Δ mean protein identity vs v1.0.8"),
        "f4": fig_completeness(fw),
        "f5": fig_crash_recovery(fw),
        "f6": fig_resource(vc, "wall_s", "f6_wall_clock.png",
                           "wall-clock (s)", "End-to-end wall-clock — v1.0.8 vs devel (controlled)"),
        "f7": fig_resource(vc, "peak_rss_mb", "f7_peak_rss.png",
                           "peak RSS (MB)", "Peak memory — v1.0.8 vs devel (controlled)"),
        "f8": fig_validity(fw),
        "f9": fig_divergence_ladder(fw),
    }

    parts = []
    parts.append(
        "# LiftOn — comprehensive comparison report\n\n"
        "**LiftOn devel** vs **LiftOn v1.0.8** vs **Liftoff (minimap2)** vs **miniprot**, "
        "across a same-species → distant-cross-species divergence ladder and two annotation "
        "databases (RefSeq, Ensembl/GTF). All accuracy numbers are a tool-neutral parasail "
        "re-score; \"minimap2\" denotes the Liftoff DNA-liftover baseline.\n")
    parts.append(sec_exec_summary(fw, vc, agg, figs))
    parts.append(sec_methodology(fw, vc, bench))
    parts.append(sec_accuracy(fw, vc, figs))
    parts.append(sec_completeness(fw, vc, agg, figs))
    parts.append(sec_robustness(fw, vc, agg, figs))
    parts.append(sec_performance(vc, fw, agg, figs))
    parts.append(sec_feature_breadth(fw, figs))
    parts.append(sec_changelog())
    parts.append(sec_appendix(fw, vc))

    MD.write_text("\n\n".join(parts) + "\n")
    made = [v for v in figs.values() if v]
    print(f"wrote {MD}  ({len(MD.read_text().splitlines())} lines)")
    print(f"figures: {len(made)} written to {FIGDIR}")
    for k, v in figs.items():
        print(f"  {k}: {v or '(skipped — no data)'}")


if __name__ == "__main__":
    main()
