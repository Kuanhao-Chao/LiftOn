#!/usr/bin/env python
"""Qualify source snapshots against the previous release, preserving evidence.

Fresh native alignment and optional copy search are explicit protocol choices.
The historical evaluator is supplemented by structural evidence; unscored models
remain unresolved. Legacy campaign rescoring never overwrites its source records.
"""
from __future__ import annotations

import argparse
import json
import math
from collections import Counter
from dataclasses import dataclass
import re
import sys
from pathlib import Path

from . import evaluator
from . import fourway_compare as fc
from . import gene_level
from . import version_compare as vc
from .profiling import run_profiled
from . import release_provenance as provenance
from lifton.run_manifest import atomic_write_json

HERE = Path(__file__).resolve().parent
OUT = HERE / "_v1012_validation"
RESULTS = OUT / "results"
FROZEN = HERE / "fourway_results.json"

#: Previous release, for the "does upgrading cost anything" question.
OLD_TREE = "/ccb/salz3/kh.chao/lifton_v1_0_11_worktree"
OLD_LABEL = "v1.0.11"
NEW_LABEL = "v1.0.12"


@dataclass(frozen=True)
class Configuration:
    output: Path
    candidate: Path
    baseline: Path
    python: str = sys.executable
    threads: int = 8
    copies: bool = False
    resume: bool = False


def _argv(python, paths, anndb, threads, out_gff, *, copies=False):
    # lifton.lifton does not have a __main__ block. Invoke the actual CLI main,
    # using the interpreter and import path whose provenance was just verified.
    argv = [python, "-c", "from lifton.lifton import main; main()",
            "-t", str(threads), "-ad", anndb, "-g", str(paths["ref_gff"]),
            "-o", str(out_gff)]
    if copies:
        argv.append("-copies")
    return argv + [str(paths["tgt_fa"]), str(paths["ref_fa"])]


def _run_arm(bid, arm, tree, paths, anndb, config, log=print):
    statedir = config.output / "cells" / bid / arm
    statedir.mkdir(parents=True, exist_ok=True)
    out_gff = statedir / f"{arm}.gff3"
    manifest = statedir / "lifton_output" / "run_manifest.json"
    receipt_path = statedir / "completion.json"
    env = provenance.isolated_env(tree)
    argv = _argv(config.python, paths, anndb, config.threads, out_gff, copies=config.copies)
    expected = provenance.run_evidence(python=config.python, root=tree, inputs=paths,
                                      argv=argv, cwd=statedir, env=env)
    receipt = provenance.read_receipt(receipt_path, expected, out_gff, manifest)
    if config.resume and receipt:
        log(f"[{bid}/{arm}] verified completion receipt; reusing output")
        return out_gff, receipt
    if out_gff.exists() or receipt_path.exists() or manifest.exists():
        raise RuntimeError(f"{statedir}: existing evidence does not authorize reuse; choose a new run ID")
    atomic_write_json(statedir / "expected.json", expected)
    profile = run_profiled(argv, label=arm, log_dir=statedir / "logs", env=env,
                           cwd=statedir, log=log)
    if profile.exit_code != 0 or not out_gff.is_file() or out_gff.stat().st_size == 0:
        raise RuntimeError(f"{bid}/{arm} failed (exit {profile.exit_code}); see {profile.stderr_path}")
    receipt = provenance.write_receipt(receipt_path, expected=expected, output=out_gff,
                                       manifest=manifest, profile=profile, validation=validate_output(out_gff))
    return out_gff, receipt


def _recovered_rows(rows):
    return [row for row in rows
            if str(row.get("is_coding", "")).strip().lower() in ("1", "true")
            and str(row.get("recovered", "")).strip().lower() in ("1", "true")]


def _finite_identity(value):
    try:
        value = float(value)
    except (TypeError, ValueError):
        return None
    return value if math.isfinite(value) and 0 <= value <= 1 else None


def _recovered_pi(rows):
    return {row["ref_mrna_id"]: value for row in _recovered_rows(rows)
            if (value := _finite_identity(row.get("protein_identity"))) is not None}


def _missing_cds_ids(rows):
    # This is a resolved biological absence, not an alignment-engine failure.
    # Never infer it from a blank score alone: require the extraction evidence.
    return {r["ref_mrna_id"] for r in rows
            if str(r.get("n_cds_lifted")) == "0"
            and str(r.get("lifted_prot_len")) == "0"
            and _finite_identity(r.get("protein_identity")) is None}


def common_set(new_rows, old_rows):
    """Track recovered IDs separately from the subset with usable scores.

    An increased count can hide lost IDs, and a non-regressing scored common
    set can hide models whose scoring failed. Both remain explicit evidence.
    """
    nr, oldr = _recovered_rows(new_rows), _recovered_rows(old_rows)
    ni, oi = Counter(r["ref_mrna_id"] for r in nr), Counter(r["ref_mrna_id"] for r in oldr)
    new, old = _recovered_pi(nr), _recovered_pi(oldr)
    shared = sorted(set(new) & set(old))
    regressed = [k for k in shared if new[k] < old[k] - 1e-9]
    mean_old = math.fsum(old[k] for k in shared) / len(shared) if shared else None
    mean_new = math.fsum(new[k] for k in shared) / len(shared) if shared else None
    new_missing, old_missing = _missing_cds_ids(nr), _missing_cds_ids(oldr)
    return {
        "n": len(shared), "shared_recovered": len(ni.keys() & oi.keys()),
        "mean_old": mean_old, "mean_new": mean_new,
        "delta": mean_new - mean_old if shared else None,
        "n_improved": sum(new[k] > old[k] + 1e-9 for k in shared),
        "n_regressed": len(regressed), "regressed_examples": regressed[:10],
        "regressed_ids": regressed,
        "lost_ids": sorted(oi.keys() - ni.keys()),
        "added_ids": sorted(ni.keys() - oi.keys()),
        "duplicate_new_ids": sorted(k for k, count in ni.items() if count > 1),
        "duplicate_old_ids": sorted(k for k, count in oi.items() if count > 1),
        "unscored_new_ids": sorted(ni.keys() - new.keys()),
        "unscored_old_ids": sorted(oi.keys() - old.keys()),
        "newly_unscored_ids": sorted((old.keys() & ni.keys()) - new.keys()),
        "newly_scored_ids": sorted((new.keys() & oi.keys()) - old.keys()),
        "missing_cds_new_ids": sorted(new_missing),
        "missing_cds_old_ids": sorted(old_missing),
        "lost_coding_model_ids": sorted(old.keys() - new.keys()),
        "unresolved_score_new_ids": sorted(ni.keys() - new.keys() - new_missing),
        "unresolved_score_old_ids": sorted(oi.keys() - old.keys() - old_missing),
        "statuses_new": dict(sorted(Counter(r.get("status", "unknown") for r in nr).items())),
        "statuses_old": dict(sorted(Counter(r.get("status", "unknown") for r in oldr).items())),
    }


def new_validity_issues(new, old):
    def counts(result):
        return Counter((i["check"], i["feature_id"]) for i in result.get("issues", [])
                       if i["severity"] == "ERROR")
    return [[check, fid, n] for (check, fid), n in sorted((counts(new) - counts(old)).items())]


def validate_output(path):
    """Use structured, uncapped issues; never interpret a crash as zero errors."""
    from lifton.gff3_validator import validate_gff3_file
    result = validate_gff3_file(str(path), max_issues_per_check=sys.maxsize)
    issues = [{"severity": i.severity, "check": i.check, "feature_id": i.feature_id,
               "line": i.lineno, "message": i.message} for i in result.issues]
    fatal = {"parse_error", "file_exists", "file_readable", "features_present"}
    return {"exit": 0 if result.is_valid else 1, "valid": result.is_valid,
            "n_errors": len(result.errors), "n_warnings": len(result.warnings),
            "issues": issues, "complete": result.data_lines > 0
            and not any(i["check"] in fatal for i in issues)}


def _summarize(rows, genes):
    summary = gene_level.gene_level_summary(gene_level.annotate_rows(rows, genes))
    coding = [r for r in rows
              if str(r.get("is_coding", "")).strip().lower() in ("1", "true")]
    recovered = [r for r in coding
                 if str(r.get("recovered", "")).strip().lower() in ("1", "true")]
    identities = [value for r in recovered
                  if (value := _finite_identity(r.get("protein_identity"))) is not None]
    summary["n_recovered_coding"] = len(recovered)
    summary["n_scored_coding"] = len(identities)
    summary["n_missing_cds"] = len(_missing_cds_ids(recovered))
    summary["mean_protein_identity"] = (
        round(sum(identities) / len(identities), 5) if identities else None)
    return summary


def _evaluation_evidence():
    return {p.name: provenance.sha256(p) for p in sorted(HERE.glob("*.py"))}


def run_cell(bid, config, log=print):
    bench = fc._bench(bid)
    paths = fc._full_paths(bid)
    anndb = bench.get("annotation_database", "RefSeq")
    root = config.output / "cells" / bid
    arms = {}
    for arm, tree in ((NEW_LABEL, config.candidate), (OLD_LABEL, config.baseline)):
        log(f"=== {bid}: {arm}, fresh alignment, -t {config.threads}, copies={config.copies} ===")
        arms[arm] = _run_arm(bid, arm, tree, paths, anndb, config, log)
    manifest = {"id": bid, "species": bench["species"], "cross_species": bench["cross_species"],
                "miniprot_target_space": "transcript", "protein_acc_to_mrna": {},
                "paths": {k: str(v) for k, v in paths.items()}}
    eval_dir = root / "eval"
    eval_dir.mkdir(parents=True, exist_ok=True)
    evaluation_evidence = _evaluation_evidence()
    ref, ref_index = evaluator.build_reference(str(paths["ref_gff"]), str(paths["ref_fa"]), log=log)
    genes = gene_level.load_transcript_genes(str(paths["ref_gff"]))
    record = {"schema_version": 2, "cell": bid, "threads": config.threads,
              "divergence": bench.get("divergence_class", ""), "arms": {},
              "protocol": {"fresh_alignment": True, "copies": config.copies},
              "evaluator": evaluation_evidence}
    rows_by_arm = {}
    for arm, (gff, receipt) in arms.items():
        evaluator.evaluate_tool(arm, str(gff), str(paths["tgt_fa"]), ref, manifest,
                                eval_dir, None, log=log, ref_index=ref_index, threads=config.threads)
        table = eval_dir / f"{arm}.transcripts.tsv"
        rows = list(gene_level.read_transcript_tsv(table))
        rows_by_arm[arm] = rows
        record["arms"][arm] = {
            "recall": _summarize(rows, genes), "validity": receipt["validation"],
            "wall_seconds": receipt["wall_seconds"], "peak_rss_mb": receipt["peak_rss_mb"],
            "output": str(gff), "completed": True, "provenance_verified": True,
            "receipt": provenance.fingerprint(gff.parent / "completion.json"),
            "transcript_table": provenance.fingerprint(table),
        }
    if _evaluation_evidence() != evaluation_evidence:
        raise RuntimeError("Evaluator source changed while scoring")
    record["reference_validity"] = validate_output(paths["ref_gff"])
    record["common_set"] = common_set(rows_by_arm[NEW_LABEL], rows_by_arm[OLD_LABEL])
    _finish(record)
    atomic_write_json(config.output / "results" / f"{bid}.json", record)
    log(f"[{bid}] lost={len(record['common_set']['lost_ids'])}; unresolved={record['unresolved_gates']}")
    return record


def _finish(record):
    """Fail closed on missing evidence, per-ID loss and new validity issues."""
    new, old = record["arms"][NEW_LABEL], record["arms"][OLD_LABEL]
    common = record["common_set"]
    nv, ov = new.get("validity", {}), old.get("validity", {})
    validity_complete = all(v.get("complete") is True and v.get("exit") in (0, 1)
                            and isinstance(v.get("n_errors"), int)
                            and isinstance(v.get("issues"), list) for v in (nv, ov))
    introduced = new_validity_issues(nv, ov) if validity_complete else None
    own_new = new["recall"].get("mean_protein_identity")
    own_old = old["recall"].get("mean_protein_identity")
    record["delta"] = {
        "recovered_coding": new["recall"]["n_recovered_coding"] - old["recall"]["n_recovered_coding"],
        "own_set_mean_identity": own_new - own_old if own_new is not None and own_old is not None else None,
        "common_set_mean_identity": common["delta"],
        "validity_errors": nv["n_errors"] - ov["n_errors"] if validity_complete else None,
        "introduced_validity_issues": introduced,
    }
    record["gate"] = {
        "both_arms_completed": all(a.get("completed") is True for a in (new, old)),
        "provenance_verified": all(a.get("provenance_verified") is True for a in (new, old)),
        "no_recovery_loss": not common["lost_ids"],
        "no_coding_model_loss": not common["lost_coding_model_ids"],
        "no_duplicate_reference_rows": not (common["duplicate_new_ids"] or common["duplicate_old_ids"]),
        "no_common_set_regression": common["n_regressed"] == 0,
        "scoring_resolved": not (common["unresolved_score_new_ids"] or common["unresolved_score_old_ids"]),
        "validity_complete": validity_complete,
        "no_new_validity_errors": introduced == [],
    }
    record["gate_pass"] = all(record["gate"].values())
    record["unresolved_gates"] = [k for k, passed in record["gate"].items() if not passed]
    return record


def rescore(bid, source, destination, log=print):
    """Re-derive legacy tables into a NEW report; never upgrade missing provenance."""
    source = Path(source)
    destination = Path(destination)
    if source.resolve() == destination.resolve():
        raise ValueError("Rescoring must preserve the original campaign directory")
    record = json.loads((source / "results" / f"{bid}.json").read_text())
    paths = fc._full_paths(bid)
    genes = gene_level.load_transcript_genes(str(paths["ref_gff"]))
    rows_by_arm = {}
    for arm in (NEW_LABEL, OLD_LABEL):
        details = record["arms"][arm]
        output = Path(details["output"])
        table = output.parent.parent / "eval" / f"{arm}.transcripts.tsv"
        rows = list(gene_level.read_transcript_tsv(table))
        rows_by_arm[arm] = rows
        details["recall"] = _summarize(rows, genes)
        details["validity"] = validate_output(output)
        manifest_path = output.parent / "lifton_output" / "run_manifest.json"
        manifest = json.loads(manifest_path.read_text()) if manifest_path.exists() else {}
        details["completed"] = manifest.get("run", {}).get("status") == "success"
        details["provenance_verified"] = False
        details["legacy_manifest"] = provenance.fingerprint(manifest_path) if manifest else None
        details["transcript_table"] = provenance.fingerprint(table)
    record["schema_version"] = 2
    record["legacy_source"] = provenance.fingerprint(source / "results" / f"{bid}.json")
    record["protocol"] = {"legacy": True, "fresh_alignment": True, "copies": True}
    record["evaluator"] = _evaluation_evidence()
    record["reference_validity"] = validate_output(paths["ref_gff"])
    record["common_set"] = common_set(rows_by_arm[NEW_LABEL], rows_by_arm[OLD_LABEL])
    _finish(record)
    atomic_write_json(destination / "results" / f"{bid}.json", record)
    log(f"[{bid}] legacy evidence rescored; unresolved={record['unresolved_gates']}")
    return record


def merge(output, expected_cells=None):
    output = Path(output)
    records = [json.loads(p.read_text()) for p in sorted((output / "results").glob("*.json"))]
    expected_cells = set(expected_cells or (r["cell"] for r in records))
    missing = sorted(expected_cells - {r["cell"] for r in records})
    lines = ["# LiftOn release qualification", "",
             "Generated from the accompanying JSON. PASS requires complete evidence; "
             "legacy records and unexplained scoring failures remain unresolved. "
             "An ID present without CDS is not a recovered protein. Own-set means are descriptive only.", "",
             "| Cell | Reference coding IDs present | Added IDs | Lost IDs | Shared scored | PI regressions | "
             "Unscored new | New validity issues | Gate |",
             "|---|---:|---:|---:|---:|---:|---:|---:|---|"]
    for r in records:
        if r.get("schema_version") != 2:
            raise ValueError("Legacy records must be rescored before merging")
        c = r["common_set"]
        issues = r["delta"]["introduced_validity_issues"]
        lines.append(f"| {r['cell']} | {r['arms'][NEW_LABEL]['recall']['n_recovered_coding']} | "
                     f"{len(c['added_ids'])} | {len(c['lost_ids'])} | {c['n']} | {c['n_regressed']} | "
                     f"{len(c['unscored_new_ids'])} | {len(issues) if issues is not None else 'unknown'} | "
                     f"{'PASS' if r['gate_pass'] else ', '.join(r['unresolved_gates'])} |")
    passed = sum(r["gate_pass"] for r in records)
    lines += ["", f"{passed}/{len(expected_cells)} expected cells pass.",
              f"Missing cells: {', '.join(missing) if missing else 'none'}.", ""]
    output.mkdir(parents=True, exist_ok=True)
    atomic_write_json(output / "release_validation.json", {
        "schema_version": 2, "expected_cells": sorted(expected_cells), "missing_cells": missing,
        "records": records, "gate_pass": bool(records) and not missing and passed == len(records)})
    (output / "release_validation.md").write_text("\n".join(lines))
    return 0 if records and not missing and passed == len(records) else 1


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("cells", nargs="*")
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--run-id", required=False)
    ap.add_argument("--output-root", type=Path, default=HERE / "_runs" / "release-qualification")
    ap.add_argument("--candidate-tree", type=Path, default=HERE.parents[1])
    ap.add_argument("--baseline-tree", type=Path, default=Path(OLD_TREE))
    ap.add_argument("--python", default=sys.executable)
    ap.add_argument("--copies", action="store_true", help="explicitly enable copy search in both arms")
    ap.add_argument("--resume", action="store_true", help="reuse only verified completion receipts")
    ap.add_argument("--merge", action="store_true")
    ap.add_argument("--rescore", type=Path, metavar="LEGACY_CAMPAIGN")
    args = ap.parse_args(argv)
    if not args.run_id or not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", args.run_id):
        ap.error("--run-id must be a nonempty path-safe identifier")
    if args.threads < 1:
        ap.error("--threads must be positive")
    if not args.cells and not (args.merge or args.rescore):
        ap.error("supply at least one benchmark cell")
    output = (args.output_root / args.run_id).resolve()
    if args.merge:
        return merge(output, args.cells)
    cells = args.cells or [p.stem for p in sorted((args.rescore / "results").glob("*.json"))]
    if not cells:
        ap.error("no cells found")
    if any(not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", cell) for cell in cells):
        ap.error("cell IDs must be path-safe")
    if output.exists() and not args.resume:
        ap.error("run directory exists; use a new --run-id or explicit --resume")
    config = Configuration(output, args.candidate_tree.resolve(), args.baseline_tree.resolve(),
                           args.python, args.threads, args.copies, args.resume)
    failures = []
    for bid in cells:
        try:
            if args.rescore:
                rescore(bid, args.rescore, output)
            else:
                run_cell(bid, config)
        except Exception as error:
            print(f"[{bid}] FAILED: {error}", file=sys.stderr, flush=True)
            failures.append(bid)
            atomic_write_json(output / "failures" / f"{bid}.json", {"cell": bid, "error": str(error)})
    result = merge(output, cells)
    return 1 if failures else result


if __name__ == "__main__":
    sys.exit(main())
