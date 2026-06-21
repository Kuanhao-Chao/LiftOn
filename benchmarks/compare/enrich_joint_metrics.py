"""One-time / idempotent enrichment of fourway_results.json with the joint
recall-vs-identity metrics (audit finding #1). Reads each cell's per-transcript
eval TSVs (already on disk) and writes a `joint` block per cell. Does NOT re-run
any lift. Re-runnable (overwrites the `joint` key each time).

Usage:
    python -m benchmarks.compare.enrich_joint_metrics            # enrich in place
    python -m benchmarks.compare.enrich_joint_metrics --dry-run  # print, no write
"""
import argparse
import json
import os
import sys

sys.path.insert(0, os.path.dirname(__file__))
import joint_metrics as jm  # noqa: E402

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
FW_PATH = os.path.join(REPO_ROOT, "benchmarks", "compare", "fourway_results.json")


def _iter_cells(data):
    if isinstance(data, list):
        return enumerate(data)
    return data.items()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    with open(FW_PATH) as fh:
        data = json.load(fh)

    n_enriched = n_skip = 0
    rows = []
    for _key, cell in _iter_cells(data):
        eval_dir = jm.eval_dir_for_cell(REPO_ROOT, cell)
        if not eval_dir or not os.path.isdir(eval_dir):
            n_skip += 1
            continue
        joint = jm.compute_joint_metrics(eval_dir, cell.get("n_reference_coding"))
        if not joint:
            n_skip += 1
            continue
        cell["joint"] = joint
        n_enriched += 1
        if cell.get("mode") == "full":
            cm = joint.get("devel_vs_miniprot_common", {})
            rows.append((cell.get("divergence_class", ""), cell.get("benchmark", ""),
                         cm.get("meanpi_delta"), cm.get("n_improved"), cm.get("n_regressed"),
                         (joint.get("covpi") or {}).get("lifton_devel"),
                         (joint.get("covpi") or {}).get("miniprot")))

    print(f"enriched {n_enriched} cells, skipped {n_skip}")
    DIV = ["same_species", "cross_species", "close_cross_species",
           "distant_cross_species", "very_distant_cross_species"]
    rows.sort(key=lambda r: (DIV.index(r[0]) if r[0] in DIV else 9, r[1]))
    print(f"\n{'cell':<30}{'tier':<10}{'dvMP_common':>12}{'imp:reg':>12}{'covPI_dev':>10}{'covPI_mp':>10}")
    for tier, b, d, imp, reg, cpd, cpm in rows:
        ds = f"{d:+.4f}" if d is not None else "NA"
        ir = f"{imp}:{reg}" if imp is not None else "NA"
        print(f"{b:<30}{tier[:9]:<10}{ds:>12}{ir:>12}{(cpd if cpd is not None else 0):>10.4f}{(cpm if cpm is not None else 0):>10.4f}")

    if args.dry_run:
        print("\n[dry-run] not written")
        return
    tmp = FW_PATH + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(data, fh, indent=2)
    os.replace(tmp, FW_PATH)
    print(f"\nwrote {FW_PATH}")


if __name__ == "__main__":
    main()
