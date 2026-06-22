"""Reporter: roll up per-benchmark/per-tool summaries into comparison tables,
matplotlib plots (base64-inlined for a self-contained HTML), and a markdown
report rendered to HTML via the repo-root ``build_spec_html.py``.
"""
from __future__ import annotations

import base64
import json
import subprocess
from io import BytesIO
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

# TOOLS drives the 3-way grouped-bar/scatter/distribution plots (kept uncluttered
# at 3 series). MASTER_TOOLS additionally includes the optional lifton_optimize
# 4th variant for the tabular views (master table + flat tsv); the optimize story
# is told in detail by the dedicated default-vs-optimize delta section.
TOOLS = ["liftoff", "miniprot", "lifton"]
MASTER_TOOLS = ["liftoff", "miniprot", "lifton", "lifton_optimize"]
TOOL_COLORS = {"liftoff": "#4c78a8", "miniprot": "#f58518", "lifton": "#54a24b",
               "lifton_optimize": "#b279a2"}
REPO_ROOT = Path(__file__).resolve().parents[2]

# sub-feature types where id-match recovery is unreliable (LiftOn re-ids CDS/exon
# during ORF rescue), so the per-type recovery view focuses on id-stable
# gene/transcript-level types. Mirrors evaluator._SUBFEATURE_TYPES.
_SUBFEATURE_TYPES = frozenset({
    "exon", "CDS", "start_codon", "stop_codon",
    "five_prime_UTR", "three_prime_UTR", "intron",
})


def _fig_to_b64(fig) -> str:
    buf = BytesIO()
    fig.savefig(buf, format="png", dpi=120, bbox_inches="tight")
    plt.close(fig)
    return base64.b64encode(buf.getvalue()).decode("ascii")


def _img_md(b64: str, alt: str) -> str:
    return f'<img alt="{alt}" src="data:image/png;base64,{b64}" style="max-width:100%;">'


def _fmt(x, pct=False, nd=3):
    if x is None or x == "":
        return "—"
    try:
        return f"{100*float(x):.1f}%" if pct else f"{float(x):.{nd}f}"
    except (TypeError, ValueError):
        return str(x)


def load(bench_ids, work_root: Path, modes=("genes",)) -> dict:
    data = {}
    for bid in bench_ids:
        wd = work_root / bid
        sm = wd / "eval" / "summaries.json"
        mf = wd / "subset" / "subset.manifest.json"
        if not sm.exists():
            continue
        summaries_by_mode = {}
        for m in modes:
            eval_sub = "eval" if m == "genes" else f"eval_{m}"
            mp = wd / eval_sub / "summaries.json"
            if mp.exists():
                summaries_by_mode[m] = json.loads(mp.read_text())
        data[bid] = {
            # genes view — every existing table/plot reads this, so they are unaffected
            "summaries": json.loads(sm.read_text()),
            "summaries_by_mode": summaries_by_mode,
            "manifest": json.loads(mf.read_text()) if mf.exists() else {},
            "crosscheck": _load_json(wd / "eval" / "liftofftools" / "crosscheck.json"),
            "work_dir": wd,
        }
    return data


def _load_json(p: Path):
    return json.loads(p.read_text()) if p.exists() else {}


# ---------------------------------------------------------------------------
# plots
# ---------------------------------------------------------------------------

def _grouped_bar(data, metric_fn, ylabel, title):
    bids = list(data.keys())
    fig, ax = plt.subplots(figsize=(max(7, 1.6 * len(bids)), 4.2))
    width = 0.26
    xs = range(len(bids))
    for i, tool in enumerate(TOOLS):
        vals = [metric_fn(data[b]["summaries"].get(tool)) for b in bids]
        vals = [v if v is not None else 0 for v in vals]
        ax.bar([x + (i - 1) * width for x in xs], vals, width,
               label=tool, color=TOOL_COLORS[tool])
    ax.set_xticks(list(xs))
    ax.set_xticklabels([data[b]["manifest"].get("id", b) for b in bids], rotation=30, ha="right")
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.legend()
    ax.grid(axis="y", alpha=0.3)
    return fig


def _identity_distribution(main_data, work_dir: Path):
    fig, ax = plt.subplots(figsize=(7, 4.2))
    for tool in TOOLS:
        tsv = work_dir / "eval" / f"{tool}.transcripts.tsv"
        if not tsv.exists():
            continue
        vals = []
        for line in tsv.read_text().splitlines()[1:]:
            c = line.split("\t")
            if len(c) >= 6 and c[5] not in ("", "—"):
                try:
                    vals.append(float(c[5]))
                except ValueError:
                    pass
        if vals:
            ax.hist(vals, bins=20, range=(0, 1), histtype="step", linewidth=2,
                    label=f"{tool} (n={len(vals)})", color=TOOL_COLORS[tool])
    ax.set_xlabel("protein-sequence identity")
    ax.set_ylabel("transcripts")
    ax.set_title("MAIN benchmark — protein identity distribution")
    ax.legend()
    ax.grid(alpha=0.3)
    return fig


def _scatter_lifton_vs_liftoff(work_dir: Path):
    def _load(tool):
        tsv = work_dir / "eval" / f"{tool}.transcripts.tsv"
        d = {}
        if tsv.exists():
            for line in tsv.read_text().splitlines()[1:]:
                c = line.split("\t")
                if len(c) >= 6 and c[5] not in ("",):
                    try:
                        d[c[0]] = float(c[5])
                    except ValueError:
                        pass
        return d
    lo, lt = _load("liftoff"), _load("lifton")
    common = sorted(set(lo) & set(lt))
    if not common:
        return None
    xs = [lo[k] for k in common]
    ys = [lt[k] for k in common]
    fig, ax = plt.subplots(figsize=(5.2, 5))
    ax.scatter(xs, ys, s=8, alpha=0.4, color="#54a24b")
    ax.plot([0, 1], [0, 1], "k--", linewidth=1)
    ax.set_xlabel("Liftoff protein identity")
    ax.set_ylabel("LiftOn protein identity")
    above = sum(1 for x, y in zip(xs, ys) if y > x + 1e-9)
    below = sum(1 for x, y in zip(xs, ys) if y < x - 1e-9)
    ax.set_title(f"MAIN — LiftOn vs Liftoff per transcript\n{above} improved, {below} worse, {len(common)} shared")
    ax.grid(alpha=0.3)
    return fig


def _resource_bar(data):
    return _grouped_bar(
        data,
        lambda s: (s or {}).get("profile", {}).get("wall_clock_seconds"),
        "wall-clock (s)", "Runtime per tool")


# ---------------------------------------------------------------------------
# markdown
# ---------------------------------------------------------------------------

def _master_table(data) -> str:
    hdr = ("| Benchmark | Tool | Completeness (coding) | Mean prot id | Median prot id "
           "| % identical | Mean DNA id | Wall (s) | Peak RSS (MB) | DNA basis |\n"
           "|---|---|---|---|---|---|---|---|---|---|\n")
    rows = []
    for bid, d in data.items():
        for tool in MASTER_TOOLS:
            s = d["summaries"].get(tool)
            if not s:
                continue
            pi, di = s["protein_identity"], s["dna_identity"]
            prof = s.get("profile", {})
            rows.append("| {} | {} | {} | {} | {} | {} | {} | {} | {} | {} |".format(
                d["manifest"].get("id", bid), tool,
                _fmt(s["completeness_coding"], pct=True),
                _fmt(pi["mean"]), _fmt(pi["median"]), _fmt(pi["pct_identical"], pct=True),
                _fmt(di["mean"]),
                _fmt(prof.get("wall_clock_seconds"), nd=1),
                _fmt(prof.get("peak_rss_mb"), nd=0),
                s["dna_basis"]))
    return hdr + "\n".join(rows) + "\n"


def _headtohead_table(data) -> str:
    hdr = ("| Benchmark | LiftOn compl. | Δ vs Liftoff | LiftOn mean prot id | Δ vs Liftoff "
           "| Δ vs miniprot |\n|---|---|---|---|---|---|\n")
    rows = []
    for bid, d in data.items():
        s = d["summaries"]
        lt, lo, mp = s.get("lifton"), s.get("liftoff"), s.get("miniprot")
        if not lt:
            continue
        def comp(x): return (x or {}).get("completeness_coding")
        def pid(x): return ((x or {}).get("protein_identity") or {}).get("mean")
        dc = (comp(lt) - comp(lo)) if comp(lt) is not None and comp(lo) is not None else None
        dp_lo = (pid(lt) - pid(lo)) if pid(lt) is not None and pid(lo) is not None else None
        dp_mp = (pid(lt) - pid(mp)) if pid(lt) is not None and pid(mp) is not None else None
        rows.append("| {} | {} | {} | {} | {} | {} |".format(
            d["manifest"].get("id", bid),
            _fmt(comp(lt), pct=True), _signed(dc, pct=True),
            _fmt(pid(lt)), _signed(dp_lo), _signed(dp_mp)))
    return hdr + "\n".join(rows) + "\n"


def _signed(x, pct=False):
    if x is None:
        return "—"
    s = f"{100*x:+.1f}%" if pct else f"{x:+.3f}"
    return s


def _load_pi_tsv(work_dir: Path, tool: str) -> dict:
    """ref_mrna_id -> protein_identity for coding, recovered, scored transcripts,
    read from the genes-mode eval TSV (mirrors score_abc.load_pi)."""
    tsv = work_dir / "eval" / f"{tool}.transcripts.tsv"
    out = {}
    if not tsv.exists():
        return out
    lines = tsv.read_text().splitlines()
    if not lines:
        return out
    hdr = lines[0].split("\t")
    try:
        i_id, i_pi, i_cod = (hdr.index("ref_mrna_id"), hdr.index("protein_identity"),
                             hdr.index("is_coding"))
    except ValueError:
        return out
    for ln in lines[1:]:
        c = ln.split("\t")
        if len(c) <= max(i_id, i_pi, i_cod) or c[i_cod] != "1" or c[i_pi] == "":
            continue
        try:
            out[c[i_id]] = float(c[i_pi])
        except ValueError:
            pass
    return out


def _optimize_delta_table(data) -> str:
    """LiftOn default vs --optimize, per benchmark: aggregate means + per-transcript
    deltas (n improved / regressed / mean Δ), then top improvements and EVERY
    regression listed explicitly. The merge is double-edged, so any regression must
    be visible, not hidden inside a favorable aggregate."""
    if not any(d["summaries"].get("lifton_optimize") for d in data.values()):
        return ""
    hdr = ("| Benchmark | Mean prot id (default) | Mean (optimize) | Δ mean | n common "
           "| improved | regressed | mean Δ/transcript |\n"
           "|---|---|---|---|---|---|---|---|\n")
    rows, details = [], []
    for bid, d in data.items():
        s = d["summaries"]
        lt, op = s.get("lifton"), s.get("lifton_optimize")
        if not lt or not op:
            continue
        m_def = (lt.get("protein_identity") or {}).get("mean")
        m_opt = (op.get("protein_identity") or {}).get("mean")
        dmean = (m_opt - m_def) if (m_def is not None and m_opt is not None) else None
        pd_, po = _load_pi_tsv(d["work_dir"], "lifton"), _load_pi_tsv(d["work_dir"], "lifton_optimize")
        common = sorted(set(pd_) & set(po))
        improved = [t for t in common if po[t] > pd_[t] + 1e-6]
        regressed = [t for t in common if po[t] < pd_[t] - 1e-6]
        mdelta = (sum(po[t] - pd_[t] for t in common) / len(common)) if common else None
        rows.append("| {} | {} | {} | {} | {} | {} | {} | {} |".format(
            d["manifest"].get("id", bid), _fmt(m_def), _fmt(m_opt), _signed(dmean),
            len(common), len(improved), len(regressed),
            f"{mdelta:+.5f}" if mdelta is not None else "—"))
        bullets = [f"  - +{po[t]-pd_[t]:.4f}  `{t}`  ({pd_[t]:.4f} → {po[t]:.4f})"
                   for t in sorted(improved, key=lambda t: po[t] - pd_[t], reverse=True)[:5]]
        bullets += [f"  - **{po[t]-pd_[t]:.4f}  `{t}`  ({pd_[t]:.4f} → {po[t]:.4f})  ← regression**"
                    for t in sorted(regressed, key=lambda t: po[t] - pd_[t])]
        if bullets:
            details.append(f"\n**{d['manifest'].get('id', bid)}** — top improvements"
                           + (" + all regressions" if regressed else "") + ":\n"
                           + "\n".join(bullets) + "\n")
    if not rows:
        return ""
    return hdr + "\n".join(rows) + "\n" + "".join(details)


def _crosscheck_table(data) -> str:
    hdr = ("| Benchmark | Tool | LiftoffTools identical | LiftoffTools n | Custom %identical×recovered | available |\n"
           "|---|---|---|---|---|---|\n")
    rows = []
    for bid, d in data.items():
        cc = d.get("crosscheck", {})
        for tool in ("liftoff", "lifton"):
            r = cc.get(tool)
            s = d["summaries"].get(tool)
            if not r:
                continue
            custom = None
            if s:
                pi = s["protein_identity"]
                if pi["pct_identical"] is not None:
                    custom = int(round(pi["pct_identical"] * s["n_recovered_coding"]))
            rows.append("| {} | {} | {} | {} | {} | {} |".format(
                d["manifest"].get("id", bid), tool,
                r.get("n_identical", "—"), r.get("n", "—"),
                custom if custom is not None else "—",
                "yes" if r.get("available") else "no"))
    return hdr + ("\n".join(rows) + "\n" if rows else "_(LiftoffTools cross-check unavailable)_\n")


def _feature_completeness_table(data) -> str:
    """Cross-benchmark overall all-feature-type completeness (Liftoff & LiftOn).

    'Overall' = recovered / reference across all gene/transcript/ncRNA/pseudogene
    feature types (exon/CDS/UTR sub-parts excluded). miniprot is protein-only and
    emits its own ids, so per-type id-recovery does not apply to it.
    """
    hdr = ("| Benchmark | Liftoff all-feature compl. | LiftOn all-feature compl. "
           "| Δ (LiftOn−Liftoff) | Ref features |\n|---|---|---|---|---|\n")
    rows = []
    for bid, d in data.items():
        s = d["summaries"]

        def overall(tool):
            cbt = (s.get(tool) or {}).get("completeness_by_type") or {}
            return cbt.get("_overall_", {}) or {}
        lo, lt = overall("liftoff"), overall("lifton")
        lo_pct, lt_pct = lo.get("pct_recovered"), lt.get("pct_recovered")
        if lo_pct is None and lt_pct is None:
            continue
        delta = (lt_pct - lo_pct) if (lo_pct is not None and lt_pct is not None) else None
        n_ref = lt.get("n_reference") or lo.get("n_reference")
        rows.append("| {} | {} | {} | {} | {} |".format(
            d["manifest"].get("id", bid),
            _fmt(lo_pct, pct=True), _fmt(lt_pct, pct=True),
            _signed(delta, pct=True), n_ref if n_ref is not None else "—"))
    if not rows:
        return ""
    return hdr + "\n".join(rows) + "\n"


def _feature_type_detail_table(data, main_id) -> str:
    """Per-feature-type recovery for the MAIN benchmark (Liftoff vs LiftOn)."""
    d = data.get(main_id)
    if not d:
        return ""
    s = d["summaries"]
    lt = (s.get("lifton") or {}).get("completeness_by_type") or {}
    lo = (s.get("liftoff") or {}).get("completeness_by_type") or {}
    # id-stable gene/transcript-level types only (sub-features are re-id'd by
    # LiftOn, so id-recovery there is not comparable — see census via TSV).
    types = [t for t in sorted(set(lt) | set(lo))
             if t != "_overall_" and t not in _SUBFEATURE_TYPES]
    if not types:
        return ""
    hdr = ("| Feature type | Ref n | Liftoff recovered | Liftoff % | LiftOn recovered "
           "| LiftOn % |\n|---|---|---|---|---|---|\n")
    rows = []
    for t in types:
        rt, ro = lt.get(t, {}), lo.get(t, {})
        n_ref = rt.get("n_reference") or ro.get("n_reference")
        rows.append("| {} | {} | {} | {} | {} | {} |".format(
            t, n_ref if n_ref is not None else "—",
            ro.get("n_recovered", "—"), _fmt(ro.get("pct_recovered"), pct=True),
            rt.get("n_recovered", "—"), _fmt(rt.get("pct_recovered"), pct=True)))
    return hdr + "\n".join(rows) + "\n"


def _feature_completeness_plot(data):
    """Grouped bar of overall all-feature completeness for the two full-annotation
    tools (Liftoff, LiftOn); miniprot omitted (protein-only)."""
    bids = list(data.keys())
    fig, ax = plt.subplots(figsize=(max(7, 1.6 * len(bids)), 4.2))
    width = 0.38
    xs = range(len(bids))
    for i, tool in enumerate(["liftoff", "lifton"]):
        vals = []
        for b in bids:
            cbt = (data[b]["summaries"].get(tool) or {}).get("completeness_by_type") or {}
            vals.append((cbt.get("_overall_", {}) or {}).get("pct_recovered") or 0)
        ax.bar([x + (i - 0.5) * width for x in xs], vals, width,
               label=tool, color=TOOL_COLORS[tool])
    ax.set_xticks(list(xs))
    ax.set_xticklabels([data[b]["manifest"].get("id", b) for b in bids],
                       rotation=30, ha="right")
    ax.set_ylabel("all-feature completeness")
    ax.set_title("Annotation completeness across ALL feature types (Liftoff vs LiftOn)")
    ax.legend()
    ax.grid(axis="y", alpha=0.3)
    return fig


def _has_allfeat(data) -> bool:
    return any(d.get("summaries_by_mode", {}).get("allfeat") for d in data.values())


def _newly_recovered_types(g_tool, a_tool):
    """Top-level types recovered in the allfeat set but NOT in the genes set
    (genes n_recovered==0, allfeat n_recovered>0). Returns [(type, n_ref, n_rec)]."""
    gcbt = (g_tool or {}).get("completeness_by_type", {}) or {}
    acbt = (a_tool or {}).get("completeness_by_type", {}) or {}
    out = []
    for t, ae in acbt.items():
        if t == "_overall_" or t in _SUBFEATURE_TYPES:
            continue
        g_rec = (gcbt.get(t, {}) or {}).get("n_recovered", 0)
        a_rec = (ae or {}).get("n_recovered", 0)
        if a_rec > 0 and g_rec == 0:
            out.append((t, ae.get("n_reference", 0), a_rec))
    return sorted(out)


def _genes_vs_allfeat_table(data) -> str:
    """Per benchmark + tool: all-feature completeness in the genes set vs the
    allfeat set, and the extra feature types/instances the allfeat set recovers."""
    hdr = ("| Benchmark | Tool | Genes all-feat compl. | AllFeat all-feat compl. | Δ "
           "| Newly-lifted types | Extra recovered |\n|---|---|---|---|---|---|---|\n")
    rows = []
    for bid, d in data.items():
        sbm = d.get("summaries_by_mode", {})
        g, a = sbm.get("genes"), sbm.get("allfeat")
        if not a:
            continue
        for tool in ("liftoff", "lifton"):
            gt, at = (g or {}).get(tool, {}), a.get(tool, {})
            gov = (gt.get("completeness_by_type", {}) or {}).get("_overall_", {}) or {}
            aov = (at.get("completeness_by_type", {}) or {}).get("_overall_", {}) or {}
            gp, ap_ = gov.get("pct_recovered"), aov.get("pct_recovered")
            if gp is None and ap_ is None:
                continue
            delta = (ap_ - gp) if (gp is not None and ap_ is not None) else None
            newly = _newly_recovered_types(gt, at)
            extra_rec = (aov.get("n_recovered") or 0) - (gov.get("n_recovered") or 0)
            rows.append("| {} | {} | {} | {} | {} | {} | {} |".format(
                d["manifest"].get("id", bid), tool,
                _fmt(gp, pct=True), _fmt(ap_, pct=True), _signed(delta, pct=True),
                ", ".join(t for t, _, _ in newly) if newly else "—",
                f"+{extra_rec}" if extra_rec else "0"))
    if not rows:
        return ""
    return hdr + "\n".join(rows) + "\n"


def _allfeat_new_types_table(data) -> str:
    """Per benchmark, the feature types LiftOn lifts ONLY in the allfeat set, with
    reference count and how many LiftOn recovered."""
    hdr = ("| Benchmark | Newly-lifted type | Ref n | LiftOn recovered (allfeat) | % |\n"
           "|---|---|---|---|---|\n")
    rows = []
    for bid, d in data.items():
        sbm = d.get("summaries_by_mode", {})
        g, a = sbm.get("genes"), sbm.get("allfeat")
        if not a:
            continue
        newly = _newly_recovered_types((g or {}).get("lifton", {}), a.get("lifton", {}))
        for t, n_ref, n_rec in newly:
            pct = (n_rec / n_ref) if n_ref else None
            rows.append("| {} | {} | {} | {} | {} |".format(
                d["manifest"].get("id", bid), t, n_ref, n_rec, _fmt(pct, pct=True)))
    if not rows:
        return ""
    return hdr + "\n".join(rows) + "\n"


def _genes_vs_allfeat_plot(data):
    """Grouped bars per benchmark: LiftOn all-feature completeness, genes vs allfeat."""
    bids = [b for b, d in data.items() if d.get("summaries_by_mode", {}).get("allfeat")]
    fig, ax = plt.subplots(figsize=(max(7, 1.6 * len(bids)), 4.2))
    width = 0.38
    xs = range(len(bids))
    for i, mode in enumerate(["genes", "allfeat"]):
        vals = []
        for b in bids:
            sm = data[b]["summaries_by_mode"].get(mode, {})
            ov = (sm.get("lifton", {}).get("completeness_by_type", {}) or {}).get("_overall_", {}) or {}
            vals.append(ov.get("pct_recovered") or 0)
        ax.bar([x + (i - 0.5) * width for x in xs], vals, width,
               label=mode, color=("#4c78a8" if mode == "genes" else "#e45756"))
    ax.set_xticks(list(xs))
    ax.set_xticklabels([data[b]["manifest"].get("id", b) for b in bids], rotation=30, ha="right")
    ax.set_ylabel("LiftOn all-feature completeness")
    ax.set_title("Feature-mode comparison — LiftOn all-feature completeness (genes vs allfeat)")
    ax.legend()
    ax.grid(axis="y", alpha=0.3)
    return fig


def write_report(data, out_md: Path, main_id: str = "human_mane", log=print) -> Path:
    if not data:
        raise RuntimeError("no benchmark data to report")
    # plots
    plots = {}
    plots["completeness"] = _img_md(_fig_to_b64(_grouped_bar(
        data, lambda s: (s or {}).get("completeness_coding"),
        "completeness (coding)", "Annotation completeness (coding transcripts)")), "completeness")
    plots["protein_id"] = _img_md(_fig_to_b64(_grouped_bar(
        data, lambda s: ((s or {}).get("protein_identity") or {}).get("mean"),
        "mean protein identity", "Mean protein-sequence identity")), "protein identity")
    plots["resources"] = _img_md(_fig_to_b64(_resource_bar(data)), "runtime")
    if main_id in data:
        wd = data[main_id]["work_dir"]
        plots["dist"] = _img_md(_fig_to_b64(_identity_distribution(data[main_id], wd)), "dist")
        sc = _scatter_lifton_vs_liftoff(wd)
        if sc is not None:
            plots["scatter"] = _img_md(_fig_to_b64(sc), "scatter")

    benches = ", ".join(d["manifest"].get("species", b) for b, d in data.items())
    md = []
    md.append("# LiftOn vs Liftoff vs miniprot — Benchmark Comparison\n")
    md.append(f"Three-way comparison of standalone **Liftoff**, standalone **miniprot**, and "
              f"**LiftOn** (reusing both via `-L`/`-M`) across {len(data)} chromosome-subset "
              f"benchmarks. Primary metrics from a custom evaluator reusing LiftOn's own "
              f"parasail kernel (`align.py` + `get_id_fraction.py`); cross-checked with "
              f"`liftofftools variants`. Benchmarks: {benches}.\n")
    md.append("> **Reading the metrics.** *Completeness (coding)* = recovered coding "
              "transcripts / reference coding transcripts. *Protein identity* is CDS-translated "
              "and directly comparable across tools. *DNA identity* basis is full-transcript "
              "(exon-spliced) for Liftoff/LiftOn but CDS-spliced for miniprot (it has no UTR/"
              "exon records) — see the `DNA basis` column. Cross-species benchmarks subset both "
              "genomes to one syntenic chromosome, so non-syntenic reference genes count as lost "
              "by design.\n")
    md.append("## 1. Master comparison\n")
    md.append(_master_table(data))
    ft_overall = _feature_completeness_table(data)
    if ft_overall:
        md.append("\n## 1b. All-feature-type completeness\n")
        md.append("Beyond coding transcripts, this tracks recovery of **every** reference "
                  "feature type (genes, mRNAs, ncRNAs, tRNAs, rRNAs, lncRNAs, pseudogenes, …) "
                  "by id. LiftOn inherits Liftoff's full-annotation lift, so the two track "
                  "closely; miniprot is protein-only (emits its own ids) and is excluded from "
                  "id-recovery here — its coding recovery is in the master table above.\n")
        md.append(ft_overall)
        md.append("\n" + _img_md(_fig_to_b64(_feature_completeness_plot(data)),
                                 "all-feature completeness") + "\n")
        detail = _feature_type_detail_table(data, main_id)
        if detail:
            md.append(f"\n**Per-feature-type detail — MAIN benchmark (`{main_id}`):**\n")
            md.append(detail)
    if _has_allfeat(data):
        md.append("\n## 1c. Feature-mode comparison (genes vs allfeat)\n")
        md.append("Two maintained benchmark sets: **genes** runs Liftoff with only the gene "
                  "hierarchy (default); **allfeat** runs Liftoff with `-f` listing every top-level "
                  "annotation type (pseudogenes, mobile elements, ncRNA genes, regulatory "
                  "features, …) and passes the same `-f` to LiftOn. This shows how many more "
                  "features the all-feature lift recovers — the gap iteration-1 surfaced.\n")
        md.append(_genes_vs_allfeat_table(data))
        md.append("\n" + _img_md(_fig_to_b64(_genes_vs_allfeat_plot(data)),
                                 "genes vs allfeat") + "\n")
        new_types = _allfeat_new_types_table(data)
        if new_types:
            md.append("\n**Feature types LiftOn lifts only in the allfeat set:**\n")
            md.append(new_types)
    md.append("\n## 2. Head-to-head — does LiftOn improve on Liftoff?\n")
    md.append("Positive Δ = LiftOn better. The thesis: LiftOn ≥ Liftoff on protein identity and "
              "completeness by fusing miniprot's protein evidence.\n")
    md.append(_headtohead_table(data))
    opt_delta = _optimize_delta_table(data)
    if opt_delta:
        md.append("\n## 2b. LiftOn default vs `--optimize` (per-transcript protein-identity delta)\n")
        md.append("The opt-in `--optimize` lane runs a best-of-outcome Liftoff↔miniprot merge "
                  "(keeps, per locus, whichever of Liftoff+ORF-rescue vs merge+ORF-rescue scores "
                  "higher) plus a `update_cds_list` multi-exon corruption fix. Positive Δ = optimize "
                  "better. The merge is **double-edged** — it rescues some transcripts hugely but can "
                  "backfire on others — so every per-transcript regression is listed explicitly below, "
                  "not hidden in the aggregate. (Deltas computed by re-aligning the emitted protein, "
                  "independent of LiftOn's internal score.)\n")
        md.append(opt_delta)
    md.append("\n## 3. Plots\n")
    md.append("### Completeness\n\n" + plots["completeness"] + "\n")
    md.append("\n### Mean protein identity\n\n" + plots["protein_id"] + "\n")
    if "dist" in plots:
        md.append("\n### MAIN benchmark — protein-identity distribution\n\n" + plots["dist"] + "\n")
    if "scatter" in plots:
        md.append("\n### MAIN benchmark — LiftOn vs Liftoff per transcript\n\n" + plots["scatter"] + "\n")
    md.append("\n### Runtime\n\n" + plots["resources"] + "\n")
    md.append("\n## 4. LiftoffTools cross-check\n")
    md.append("Independent confirmation on the Liftoff/LiftOn outputs (miniprot excluded — no "
              "gene hierarchy). `identical` count should track the custom evaluator's "
              "%identical × recovered.\n")
    md.append(_crosscheck_table(data))
    md.append("\n## 5. Pipeline\n")
    md.append("Reproduce with `python -m benchmarks.compare.run_compare --all -t 8` (resumable). "
              "Per benchmark: subset both genomes to one chromosome → run Liftoff and miniprot "
              "standalone → run LiftOn reusing both (`-L`/`-M`) → evaluate each output → "
              "cross-check with LiftoffTools → aggregate. Code under `benchmarks/compare/`.\n")

    out_md.parent.mkdir(parents=True, exist_ok=True)
    out_md.write_text("\n".join(md))
    log(f"  [report] wrote {out_md}")
    return out_md


def render_html(md_path: Path, html_path: Path, log=print) -> bool:
    builder = REPO_ROOT / "build_spec_html.py"
    mermaid = REPO_ROOT / "spec_assets" / "mermaid.min.js"
    if not builder.exists():
        log(f"  [report] build_spec_html.py not found; skipping HTML")
        return False
    argv = ["python", str(builder), "--md", str(md_path), "--out", str(html_path),
            "--title", "LiftOn vs Liftoff vs miniprot — Benchmark"]
    if mermaid.exists():
        argv += ["--mermaid", str(mermaid)]
    proc = subprocess.run(argv, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    log("  [report] " + (proc.stdout or "").strip().splitlines()[-1] if proc.stdout else "")
    return proc.returncode == 0


def build_master_json_tsv(data, out_dir: Path):
    out_dir.mkdir(parents=True, exist_ok=True)
    flat = []
    for bid, d in data.items():
        by_mode = d.get("summaries_by_mode") or {"genes": d["summaries"]}
        for mode, summ in by_mode.items():
            for tool in MASTER_TOOLS:
                s = summ.get(tool)
                if not s:
                    continue
                pi, di = s["protein_identity"], s["dna_identity"]
                prof = s.get("profile", {})
                flat.append({
                    "benchmark": bid, "feature_mode": mode, "species": s["species"], "tool": tool,
                    "completeness_coding": s["completeness_coding"],
                    "completeness_all": s["completeness_all"],
                    "completeness_feature_total": s.get("completeness_feature_total"),
                    "n_recovered_coding": s["n_recovered_coding"],
                    "n_reference_coding": s["n_reference_coding"],
                    "mean_protein_identity": pi["mean"], "median_protein_identity": pi["median"],
                    "pct_identical": pi["pct_identical"], "mean_dna_identity": di["mean"],
                    "dna_basis": s["dna_basis"], "n_extra_copies": s["n_extra_copies"],
                    "wall_clock_seconds": prof.get("wall_clock_seconds"),
                    "peak_rss_mb": prof.get("peak_rss_mb"),
                })
    (out_dir / "comparison.json").write_text(json.dumps(flat, indent=2))
    if flat:
        cols = list(flat[0].keys())
        with open(out_dir / "comparison.tsv", "w") as fh:
            fh.write("\t".join(cols) + "\n")
            for r in flat:
                fh.write("\t".join("" if r[c] is None else str(r[c]) for c in cols) + "\n")
    return flat
