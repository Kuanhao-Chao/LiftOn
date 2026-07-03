#!/usr/bin/env python
"""apply_devel_refresh.py — fold the per-pair `devel_refresh/<pair>.json` patches
(produced by `devel_refresh.py`) into the shared `fourway_results.json`.

ONLY the `lifton_devel` fields of the 17 `full:<pair>` cells change — Liftoff,
miniprot, and `lifton_stable` (v1.0.8) are UNCHANGED by this session's work, so
their cached fourway_results.json numbers stay byte-identical, and all 34 subset
cells are untouched. The refresh ran `--no-miniprot-rescue` so the devel column
stays apples-to-apples with the committed baselines (the Iter-23 pin).

For each refreshed pair this:
  1. overwrites the `lifton_devel` key in every per-tool scalar dict of the
     `full:<pair>` cell from the refresh `summary`;
  2. recomputes `devel_vs_stable` / `devel_vs_best_baseline` by REUSING
     `fourway_compare._sub` / `_best_baseline` (exactly as `_score_and_record`);
  3. copies the refresh's `lifton_devel.{transcripts,feature_types}.tsv` over the
     stale devel TSVs in the canonical `_fourway_full/eval/` dir (baseline TSVs
     untouched) so the `joint` block can be recomputed from the on-disk 4-tool
     TSVs.
After patching every cell it runs `enrich_joint_metrics` once to refresh the
`joint` block of every cell from those TSVs (reusing the shipped audit-#1 code).

Idempotent. Dry-run by default; `--apply` writes. A post-write guard asserts that
only `full:*` cells' `lifton_devel`/`devel_vs_*`/`joint` keys moved.

Usage (repo root, lifton_devel env):
    python -m benchmarks.compare.apply_devel_refresh            # dry-run report
    python -m benchmarks.compare.apply_devel_refresh --apply    # write
"""
from __future__ import annotations
import argparse
import copy
import json
import shutil
import subprocess
import sys
from pathlib import Path

from . import fourway_compare as fc
from . import joint_metrics as jm

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parent.parent
RESULTS = HERE / "fourway_results.json"
REFRESH_DIR = HERE / "devel_refresh"
DEV = "lifton_devel"
STA = "lifton_stable"

# The eval TSVs the joint-metric recompute reads, keyed by tool.
_TSV_SUFFIXES = ("transcripts.tsv", "feature_types.tsv")


def _summary_to_cell_fields(s: dict) -> dict:
    """Map a refresh `summary` dict to the per-tool scalar values a fourway cell
    stores, exactly as `fourway_compare._score_and_record` derives them."""
    pi = s["protein_identity"]
    prof = s.get("profile") or {}
    return {
        "completeness_coding": s["completeness_coding"],
        "completeness_feature_total": s["completeness_feature_total"],
        "mean_pi": pi["mean"],
        "median_pi": pi["median"],
        "pct_identical": pi["pct_identical"],
        "n_recovered_coding": s["n_recovered_coding"],
        "n_recovered_any": s["n_recovered_any"],
        "wall_s": prof.get("wall_clock_seconds"),
        "peak_rss_mb": prof.get("peak_rss_mb"),
        "feature_census": s.get("completeness_by_type", {}),
    }


def _patch_cell(cell: dict, refresh: dict) -> dict:
    """Return a per-pair report of the devel deltas after patching `cell` in place."""
    s = refresh["summary"]
    v = refresh["validity"]
    fields = _summary_to_cell_fields(s)

    before_mpi = cell["mean_pi"].get(DEV)
    before_comp = cell["completeness_coding"].get(DEV)
    before_nrec = cell["n_recovered_coding"].get(DEV)
    before_val = (cell.get("validity", {}).get(DEV) or {}).get("n_errors")

    # 1. per-tool scalar dicts (only the DEV key moves)
    for k in ("completeness_coding", "completeness_feature_total", "mean_pi",
              "median_pi", "pct_identical", "n_recovered_coding", "n_recovered_any",
              "wall_s", "peak_rss_mb", "feature_census"):
        cell.setdefault(k, {})[DEV] = fields[k]
    cell.setdefault("validity", {})[DEV] = v

    # 2. recompute devel deltas exactly like _score_and_record (REUSE fc helpers)
    mean_pi, comp, nrec = cell["mean_pi"], cell["completeness_coding"], cell["n_recovered_coding"]
    cell["devel_vs_stable"] = {
        "meanpi": fc._sub(mean_pi.get(DEV), mean_pi.get(STA)),
        "completeness": fc._sub(comp.get(DEV), comp.get(STA)),
        "n_recovered": fc._sub(nrec.get(DEV), nrec.get(STA)),
    }
    cell["devel_vs_best_baseline"] = {
        "meanpi": fc._sub(mean_pi.get(DEV), fc._best_baseline(mean_pi)),
        "completeness": fc._sub(comp.get(DEV), fc._best_baseline(comp)),
    }

    return {
        "mean_pi": (before_mpi, mean_pi.get(DEV)),
        "completeness_coding": (before_comp, comp.get(DEV)),
        "n_recovered_coding": (before_nrec, nrec.get(DEV)),
        "val_errors": (before_val, (v or {}).get("n_errors")),
    }


def _copy_devel_tsvs(bid: str, refresh: dict, apply: bool, log=print) -> bool:
    """Copy the refresh's devel eval TSVs over the canonical _fourway_full/eval
    devel TSVs so enrich_joint_metrics recomputes joint from refreshed data.
    Baseline TSVs are never touched. Returns True if the source TSVs exist."""
    src_dir = Path(refresh["eval_dir"])
    dst_dir = HERE / "work" / bid / "_fourway_full" / "eval"
    ok = True
    for suf in _TSV_SUFFIXES:
        src = src_dir / f"{DEV}.{suf}"
        dst = dst_dir / f"{DEV}.{suf}"
        if not src.exists() or src.stat().st_size == 0:
            log(f"    [WARN] {bid}: missing refresh TSV {src} — joint block will be stale")
            ok = False
            continue
        if apply:
            dst_dir.mkdir(parents=True, exist_ok=True)
            shutil.copy2(src, dst)
        log(f"    {'copy' if apply else 'would copy'} {src.name} -> {dst}")
    return ok


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--apply", action="store_true", help="write (default: dry-run report)")
    ap.add_argument("--skip-enrich", action="store_true",
                    help="don't run enrich_joint_metrics after patching")
    a = ap.parse_args(argv)

    db = json.loads(RESULTS.read_text())
    before = copy.deepcopy(db)

    refreshes = sorted(REFRESH_DIR.glob("*.json"))
    failed = sorted(REFRESH_DIR.glob("*.FAILED.txt"))
    if failed:
        print(f"[WARN] {len(failed)} pair(s) FAILED the refresh (skipped): "
              f"{[p.stem for p in failed]}", file=sys.stderr)

    print(f"=== apply_devel_refresh ({'APPLY' if a.apply else 'DRY-RUN'}) ===")
    print(f"{len(refreshes)} refreshed pair(s) found\n")
    print(f"{'pair':<28}{'meanpi Δ':>22}{'compl Δ':>22}{'nrec Δ':>16}{'valerr':>12}")

    patched = []
    for p in refreshes:
        bid = p.stem
        key = f"full:{bid}"
        refresh = json.loads(p.read_text())
        if key not in db:
            print(f"  [SKIP] {bid}: no `{key}` cell in fourway_results.json")
            continue
        rep = _patch_cell(db[key], refresh)
        _copy_devel_tsvs(bid, refresh, a.apply, log=print)
        patched.append(bid)

        def fmt(pair):
            o, n = pair
            return f"{o}→{n}" if o != n else f"={n}"
        print(f"{bid:<28}{fmt(rep['mean_pi']):>22}{fmt(rep['completeness_coding']):>22}"
              f"{fmt(rep['n_recovered_coding']):>16}{fmt(rep['val_errors']):>12}")

    # Guard: only full:* cells' DEV/devel_vs_*/joint keys may have moved.
    _assert_scope(before, db, patched)

    if a.apply:
        RESULTS.write_text(json.dumps(db, indent=2))
        print(f"\n[APPLY] wrote {RESULTS} ({len(patched)} full cells patched)")
        if not a.skip_enrich:
            print("[APPLY] running enrich_joint_metrics to refresh joint blocks...")
            r = subprocess.run([sys.executable, "-m", "benchmarks.compare.enrich_joint_metrics"],
                               cwd=str(REPO_ROOT))
            if r.returncode != 0:
                print("[WARN] enrich_joint_metrics exited non-zero; joint blocks may be stale",
                      file=sys.stderr)
            # enrich_joint_metrics recomputes EVERY cell's joint from on-disk TSVs.
            # That over-reaches this refresh's scope in two ways we must undo:
            #  (i) subset cells (out of scope, frozen) whose on-disk devel TSVs drifted
            #      since commit get a recomputed joint — restore them byte-for-byte;
            # (ii) a baseline tool whose eval TSV is no longer on disk (lifton_stable
            #      CRASHED partway on full arabidopsis/rice, so its TSV was never
            #      retained) gets nulled in covpi/recall_at_*/devel_vs_<tool>_common —
            #      restore just those committed baseline contributions (they can't be
            #      recomputed) while KEEPING enrich's legitimate devel + present-baseline
            #      refresh.
            _restore_scope_after_enrich(before, patched)
    else:
        print(f"\n[DRY-RUN] {len(patched)} full cells would be patched; re-run with --apply")


# joint sub-dicts keyed per-tool that enrich nulls when a tool's TSV is missing.
_JOINT_TOOL_KEYED = ("covpi", "recall_at_0.5", "recall_at_0.9")


def _restore_scope_after_enrich(before: dict, patched: list):
    """After enrich rewrote the file, re-read it and undo out-of-scope joint moves:
    freeze every subset / unpatched cell to committed, and in each patched full cell
    restore only the committed baseline joint contributions that enrich dropped to
    None (missing-TSV baselines). Idempotent; leaves devel + present-baseline joint
    refresh intact. Then re-verifies scope and rewrites."""
    final = json.loads(RESULTS.read_text())
    patched_keys = {f"full:{b}" for b in patched}
    n_frozen = 0
    n_baseline_restored = 0
    for key, cell_b in before.items():
        cell_f = final[key]
        if key not in patched_keys:
            if cell_f != cell_b:               # subset / unpatched full cell drifted
                final[key] = copy.deepcopy(cell_b)
                n_frozen += 1
            continue
        jb = cell_b.get("joint") or {}
        jf = cell_f.get("joint") or {}
        for sub in _JOINT_TOOL_KEYED:          # covpi / recall_at_* per-tool dicts
            db_sub = jb.get(sub) or {}
            df_sub = jf.get(sub)
            if df_sub is None:
                continue
            for tool in ("liftoff", "miniprot", STA):
                if db_sub.get(tool) is not None and df_sub.get(tool) is None:
                    df_sub[tool] = db_sub[tool]
                    n_baseline_restored += 1
        for tool in ("liftoff", "miniprot", STA):   # devel_vs_<baseline>_common
            fld = f"devel_vs_{tool}_common"
            if jb.get(fld) is not None and jf.get(fld) is None:
                jf[fld] = copy.deepcopy(jb[fld])
                n_baseline_restored += 1
    RESULTS.write_text(json.dumps(final, indent=2))
    print(f"[APPLY] post-enrich scope restore: froze {n_frozen} out-of-scope cell(s) "
          f"to committed; restored {n_baseline_restored} missing-TSV baseline joint "
          f"field(s) (arabidopsis/rice lifton_stable crash residual)")
    _assert_scope_final(before, final, patched)


def _assert_scope_final(before: dict, final: dict, patched: list):
    """Post-restore invariant: every subset + unpatched cell byte-identical to
    committed; in each patched full cell every committed BASELINE number is unchanged
    — scalars (liftoff/miniprot/lifton_stable) AND the baseline covpi/recall_at_*
    joint entries. The devel_vs_<tool>_common blocks are DEVEL-side comparisons and
    are EXPECTED to move (new devel vs the same baseline), so they are not frozen; we
    only assert no baseline joint contribution was left None (missing-TSV baselines
    must have been restored)."""
    patched_keys = {f"full:{b}" for b in patched}
    for key, cell_b in before.items():
        cell_f = final[key]
        if key not in patched_keys:
            if cell_f != cell_b:
                raise AssertionError(f"FINAL SCOPE VIOLATION: frozen cell {key} moved")
            continue
        # scalar baselines unchanged
        for fld in ("completeness_coding", "completeness_feature_total", "mean_pi",
                    "median_pi", "pct_identical", "n_recovered_coding",
                    "n_recovered_any", "wall_s", "peak_rss_mb", "validity",
                    "feature_census"):
            for tool in ("liftoff", "miniprot", STA):
                if (cell_b.get(fld) or {}).get(tool) != (cell_f.get(fld) or {}).get(tool):
                    raise AssertionError(f"FINAL SCOPE VIOLATION: {key}.{fld}[{tool}] moved")
        # joint per-tool BASELINE entries (covpi / recall_at_*) unchanged; and no
        # baseline contribution left None (missing-TSV baseline must be restored).
        jb, jf = cell_b.get("joint") or {}, cell_f.get("joint") or {}
        for sub in _JOINT_TOOL_KEYED:
            for tool in ("liftoff", "miniprot", STA):
                vb = (jb.get(sub) or {}).get(tool)
                vf = (jf.get(sub) or {}).get(tool)
                if vb != vf:
                    raise AssertionError(
                        f"FINAL SCOPE VIOLATION: {key}.joint.{sub}[{tool}] {vb}->{vf}")
        # a baseline devel_vs_<tool>_common that enrich nulled must have been restored
        for tool in ("liftoff", "miniprot", STA):
            fld = f"devel_vs_{tool}_common"
            if jb.get(fld) is not None and jf.get(fld) is None:
                raise AssertionError(
                    f"FINAL SCOPE VIOLATION: {key}.joint.{fld} left None (not restored)")
    print("[guard] FINAL scope OK: subset + baseline numbers (scalars AND joint "
          f"covpi/recall) byte-identical to committed; only the lifton_devel side of "
          f"{len(patched)} full cells moved (devel scalars + devel_vs_* + devel joint)")


def _assert_scope(before: dict, after: dict, patched: list):
    """Fail loudly if anything outside the intended scope changed:
    - every subset cell is byte-identical;
    - every full cell NOT in `patched` is byte-identical;
    - in each patched cell, only lifton_devel entries + devel_vs_* moved
      (baseline per-tool values are unchanged). The joint block is refreshed
      separately by enrich_joint_metrics, so it is exempted here."""
    patched_keys = {f"full:{b}" for b in patched}
    for key, cell_b in before.items():
        cell_a = after[key]
        if key not in patched_keys:
            if cell_a != cell_b:
                raise AssertionError(f"SCOPE VIOLATION: unpatched cell {key} changed")
            continue
        # patched full cell: baseline per-tool scalars must be identical
        for fld in ("completeness_coding", "completeness_feature_total", "mean_pi",
                    "median_pi", "pct_identical", "n_recovered_coding",
                    "n_recovered_any", "wall_s", "peak_rss_mb", "validity",
                    "feature_census"):
            for tool in ("liftoff", "miniprot", "lifton_stable"):
                vb = (cell_b.get(fld) or {}).get(tool)
                va = (cell_a.get(fld) or {}).get(tool)
                if vb != va:
                    raise AssertionError(
                        f"SCOPE VIOLATION: {key}.{fld}[{tool}] changed {vb}->{va}")
    print("\n[guard] scope OK: baselines + subset cells byte-identical; "
          f"only lifton_devel + devel_vs_* of {len(patched)} full cells moved")


if __name__ == "__main__":
    main()
