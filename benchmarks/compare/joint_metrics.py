"""Joint recall-vs-identity metrics for the 4-way comparison (audit finding #1).

The headline ``mean_pi`` / ``devel_vs_best_baseline`` are computed as
*mean-PI-over-each-tool's-OWN-recovered-set*. When two tools recover wildly
different fractions of the reference (as at the very-distant tier, where devel
recovers the easy 4-24% DNA-anchored core while miniprot recovers 67-86%), that
set-mean comparison is a denominator artifact: it compares means over different
populations. This module adds the honest, apples-to-apples view:

  * ``devel_vs_<tool>_common`` — mean-PI delta over the COMMON recovered-coding
    set (transcripts both tools recover), plus a paired ``n_improved:n_regressed``
    sign-test count.
  * ``covpi`` — coverage-weighted PI per tool (= Sum(PI over recovered) /
    n_reference_coding = recall x accuracy): a single number that a tool cannot
    win by sacrificing recall.
  * ``recall_at_0.5`` / ``recall_at_0.9`` — fraction of the reference coding set
    recovered at protein-identity >= threshold.

All values are derived from the per-transcript ``<tool>.transcripts.tsv`` files
the evaluator already writes; nothing here re-runs a lift. Shared by
``enrich_joint_metrics.py`` (post-hoc) and ``fourway_compare.py`` (live).
"""
import csv
import os
import statistics as st

TOOLS = ["liftoff", "miniprot", "lifton_stable", "lifton_devel"]
_TRUTHY = ("1", "true", "yes", "t")


def _truthy(x):
    return str(x).strip().lower() in _TRUTHY


def load_recovered_coding(eval_dir, tool):
    """Return {ref_mrna_id: max protein_identity} over recovered coding rows, or
    None if the tool's TSV is absent."""
    path = os.path.join(eval_dir, f"{tool}.transcripts.tsv")
    if not os.path.exists(path):
        return None
    rec = {}
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if not _truthy(row.get("is_coding")) or not _truthy(row.get("recovered")):
                continue
            try:
                pi = float(row["protein_identity"])
            except (ValueError, TypeError, KeyError):
                continue
            rid = row["ref_mrna_id"]
            if rid not in rec or pi > rec[rid]:
                rec[rid] = pi
    return rec


def _paired(devel_rec, base_rec, eps=1e-9):
    keys = set(devel_rec) & set(base_rec)
    if not keys:
        return None
    delta = st.mean(devel_rec[k] for k in keys) - st.mean(base_rec[k] for k in keys)
    imp = sum(1 for k in keys if devel_rec[k] - base_rec[k] > eps)
    reg = sum(1 for k in keys if base_rec[k] - devel_rec[k] > eps)
    return {
        "meanpi_delta": round(delta, 5),
        "n_common": len(keys),
        "n_improved": imp,
        "n_regressed": reg,
    }


def compute_joint_metrics(eval_dir, n_reference_coding):
    """Compute the joint recall-vs-identity block for one cell from its eval dir.
    Returns {} if the devel TSV is missing (e.g. a crashed cell)."""
    recs = {t: load_recovered_coding(eval_dir, t) for t in TOOLS}
    out = {}
    denom = n_reference_coding or 0

    covpi, recall05, recall09 = {}, {}, {}
    for t, r in recs.items():
        if r is None or not denom:
            continue
        covpi[t] = round(sum(r.values()) / denom, 5)
        recall05[t] = round(sum(1 for v in r.values() if v >= 0.5) / denom, 5)
        recall09[t] = round(sum(1 for v in r.values() if v >= 0.9) / denom, 5)
    if covpi:
        out["covpi"] = covpi
        out["recall_at_0.5"] = recall05
        out["recall_at_0.9"] = recall09

    dv = recs.get("lifton_devel")
    if dv:
        for base in ("miniprot", "liftoff", "lifton_stable"):
            br = recs.get(base)
            if br is None:
                continue
            paired = _paired(dv, br)
            if paired:
                out[f"devel_vs_{base}_common"] = paired
    return out


def eval_dir_for_cell(repo_root, cell):
    """Resolve the per-transcript eval dir for a fourway_results.json cell.
    Full cells live under _fourway_full/eval; subset cells under _fourway/eval."""
    bench = cell.get("benchmark")
    if not bench:
        return None
    mode = cell.get("mode", "subset")
    sub = "_fourway_full" if mode == "full" else "_fourway"
    return os.path.join(repo_root, "benchmarks", "compare", "work", bench, sub, "eval")
