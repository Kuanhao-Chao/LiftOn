"""Joint recall-vs-identity metrics for the 4-way comparison (audit finding #1).

The headline ``mean_pi`` / ``devel_vs_best_baseline`` are computed as
*mean-PI-over-each-tool's-OWN-recovered-set*. When two tools recover wildly
different fractions of the reference (as at the very-distant tier, where devel
recovers the easy 4-24% DNA-anchored core while miniprot recovers 67-86%), that
set-mean comparison is a denominator artifact: it compares means over different
populations. This module adds the honest, apples-to-apples view:

  * ``devel_vs_<tool>_common`` — mean-PI delta over the COMMON recovered-coding
    set (transcripts both tools recover), paired sign counts, and a deterministic
    paired-bootstrap 95% confidence interval.
  * ``covpi`` — coverage-weighted PI per tool (= Sum(PI over recovered) /
    n_reference_coding = recall x accuracy): a single number that a tool cannot
    win by sacrificing recall.
  * ``recall_at_<threshold>`` — fraction of the reference coding set recovered
    at PI >= 0.50, 0.75, 0.90, and 0.95.
  * ``structural`` — intron-chain/ORF recall and exon sensitivity/specificity
    when the evaluator TSV contains those truth-based columns.

All values are derived from the per-transcript ``<tool>.transcripts.tsv`` files
the evaluator already writes; nothing here re-runs a lift. Shared by
``enrich_joint_metrics.py`` (post-hoc) and ``fourway_compare.py`` (live).
"""
import csv
import hashlib
import os
import random
import statistics as st

TOOLS = ["liftoff", "miniprot", "lifton_stable", "lifton_devel"]
_TRUTHY = ("1", "true", "yes", "t")
RECALL_THRESHOLDS = (0.5, 0.75, 0.9, 0.95)
BOOTSTRAP_REPLICATES = 1000


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


def _percentile(values, quantile):
    """Return a linearly interpolated percentile without a NumPy dependency."""
    ordered = sorted(values)
    if not ordered:
        return None
    if len(ordered) == 1:
        return ordered[0]
    position = (len(ordered) - 1) * quantile
    lower = int(position)
    upper = min(lower + 1, len(ordered) - 1)
    fraction = position - lower
    return ordered[lower] * (1.0 - fraction) + ordered[upper] * fraction


def _paired_bootstrap(deltas, keys, replicates=BOOTSTRAP_REPLICATES):
    """Return a deterministic paired-bootstrap 95% CI for mean PI delta."""
    if not deltas:
        return None
    if len(deltas) == 1:
        value = round(deltas[0], 5)
        return {"low": value, "high": value, "replicates": 1}

    seed_material = "\0".join(keys).encode("utf-8", "surrogateescape")
    seed = int.from_bytes(hashlib.sha256(seed_material).digest()[:8], "big")
    count = len(deltas)
    replicates = max(1, int(replicates))
    try:
        import numpy as np

        values = np.asarray(deltas, dtype=float)
        rng = np.random.default_rng(seed)
        # Keep the temporary index matrix around 8 MiB even for chromosome-
        # scale transcript sets; scoring cost remains O(n * replicates) but
        # memory does not grow with the full bootstrap matrix.
        batch_size = max(1, min(128, 1_000_000 // count))
        means = []
        remaining = replicates
        while remaining:
            batch = min(batch_size, remaining)
            indices = rng.integers(
                0, count, size=(batch, count), dtype=np.int32,
            )
            means.extend(values[indices].mean(axis=1).tolist())
            remaining -= batch
    except ImportError:  # pragma: no cover - NumPy is a LiftOn dependency.
        rng = random.Random(seed)
        means = [
            sum(deltas[rng.randrange(count)] for _ in range(count)) / count
            for _ in range(replicates)
        ]
    return {
        "low": round(_percentile(means, 0.025), 5),
        "high": round(_percentile(means, 0.975), 5),
        "replicates": len(means),
    }


def _paired(devel_rec, base_rec, eps=1e-9):
    keys = sorted(set(devel_rec) & set(base_rec))
    if not keys:
        return None
    deltas = [devel_rec[k] - base_rec[k] for k in keys]
    delta = st.mean(deltas)
    imp = sum(1 for value in deltas if value > eps)
    reg = sum(1 for value in deltas if value < -eps)
    return {
        "meanpi_delta": round(delta, 5),
        "devel_meanpi": round(st.mean(devel_rec[k] for k in keys), 5),
        "baseline_meanpi": round(st.mean(base_rec[k] for k in keys), 5),
        "bootstrap_95ci": _paired_bootstrap(deltas, keys),
        "n_common": len(keys),
        "n_improved": imp,
        "n_regressed": reg,
        "n_tied": len(keys) - imp - reg,
    }


def _load_structural_rows(eval_dir, tool):
    """Load recovered coding rows that carry optional truth-based metrics."""
    path = os.path.join(eval_dir, f"{tool}.transcripts.tsv")
    if not os.path.exists(path):
        return None
    rows = []
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if not _truthy(row.get("is_coding")) or not _truthy(row.get("recovered")):
                continue
            rows.append(row)
    return rows


def _structural_metrics(eval_dir, n_reference_coding):
    """Summarize structural accuracy only when evaluator truth columns exist."""
    denominator = int(n_reference_coding or 0)
    if denominator <= 0:
        return {}
    out = {}
    for tool in TOOLS:
        rows = _load_structural_rows(eval_dir, tool)
        if rows is None:
            continue
        metrics = {}
        for key in ("intron_chain_exact", "orf_valid"):
            scored = [row.get(key) for row in rows if row.get(key) not in (None, "")]
            if scored:
                successes = sum(1 for value in scored if _truthy(value))
                metrics[f"{key}_recall"] = round(successes / denominator, 5)
                metrics[f"{key}_n_scored"] = len(scored)
        for key in ("exon_sn", "exon_sp"):
            values = []
            for row in rows:
                try:
                    values.append(float(row[key]))
                except (KeyError, TypeError, ValueError):
                    continue
            if values:
                metrics[f"{key}_mean"] = round(st.mean(values), 5)
                metrics[f"{key}_n_scored"] = len(values)
        if metrics:
            out[tool] = metrics
    return out


def compute_joint_metrics(eval_dir, n_reference_coding):
    """Compute the joint recall-vs-identity block for one cell from its eval dir.
    Returns {} if the devel TSV is missing (e.g. a crashed cell)."""
    recs = {t: load_recovered_coding(eval_dir, t) for t in TOOLS}
    out = {}
    denom = n_reference_coding or 0

    covpi = {}
    recalls = {threshold: {} for threshold in RECALL_THRESHOLDS}
    for t, r in recs.items():
        if r is None or not denom:
            continue
        covpi[t] = round(sum(r.values()) / denom, 5)
        for threshold in RECALL_THRESHOLDS:
            recalls[threshold][t] = round(
                sum(1 for value in r.values() if value >= threshold) / denom,
                5,
            )
    if covpi:
        out["covpi"] = covpi
        for threshold, values in recalls.items():
            out[f"recall_at_{threshold:g}"] = values

    dv = recs.get("lifton_devel")
    if dv:
        for base in ("miniprot", "liftoff", "lifton_stable"):
            br = recs.get(base)
            if br is None:
                continue
            paired = _paired(dv, br)
            if paired:
                out[f"devel_vs_{base}_common"] = paired
    structural = _structural_metrics(eval_dir, n_reference_coding)
    if structural:
        out["structural"] = structural
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
