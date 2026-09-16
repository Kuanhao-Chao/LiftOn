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
import tempfile
import uuid
from itertools import chain
from pathlib import Path

from . import evaluator
from . import fourway_compare as fc
from . import gene_level
from .profiling import run_profiled
from . import release_provenance as provenance
from lifton.run_manifest import atomic_write_json

HERE = Path(__file__).resolve().parent
# Historical arm labels remain readable; fresh campaigns use roles.
OLD_LABEL = "v1.0.11"
NEW_LABEL = "v1.0.12"
CANDIDATE = "candidate"
REFERENCE = "reference"
REPORT_SCHEMA_VERSION = 3


def _arm_labels(record):
    """Read version-neutral reports and preserved historical paired reports."""
    labels = set(record["arms"])
    if labels == {CANDIDATE, REFERENCE}:
        return CANDIDATE, REFERENCE
    if labels == {NEW_LABEL, OLD_LABEL}:
        return NEW_LABEL, OLD_LABEL
    raise ValueError("Report must contain exactly candidate/reference or the historical release pair")


@dataclass(frozen=True)
class Configuration:
    output: Path
    candidate: Path
    baseline: Path
    python: str = sys.executable
    threads: int = 8
    copies: bool = False
    resume: bool = False
    paths: dict | None = None
    metadata: dict | None = None
    role_options: dict | None = None


def _argv(python, paths, anndb, threads, out_gff, *, root, copies=False, options=()):
    # lifton.lifton does not have a __main__ block. Invoke the actual CLI main,
    # using the interpreter and import path whose provenance was just verified.
    argv = [python, "-c", provenance.guarded_code(root, "from lifton.lifton import main; main()"),
            "-t", str(threads), "-ad", anndb, "-g", str(paths["ref_gff"]),
            "-o", str(out_gff)]
    if copies:
        argv.append("-copies")
    return argv + list(options) + [str(paths["tgt_fa"]), str(paths["ref_fa"])]


def _run_arm(bid, arm, tree, paths, anndb, config, log=print):
    statedir = config.output / "cells" / bid / arm
    statedir.mkdir(parents=True, exist_ok=True)
    out_gff = statedir / f"{arm}.gff3"
    manifest = statedir / "lifton_output" / "run_manifest.json"
    receipt_path = statedir / "completion.json"
    env = provenance.isolated_env(tree)
    env.update({'OMP_NUM_THREADS': '1', 'OPENBLAS_NUM_THREADS': '1', 'MKL_NUM_THREADS': '1',
                'LIFTON_MINIPROT_THREADS': str(config.threads)})
    argv = _argv(config.python, paths, anndb, config.threads, out_gff, root=tree, copies=config.copies,
                 options=(config.role_options or {}).get(arm, ()))
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
                           cwd=statedir, inherit_env=False, log=log)
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


def common_set(new_rows, old_rows, *, expected_ids=None, expected_coding_ids=None):
    """Track recovered IDs separately from the subset with usable scores.

    An increased count can hide lost IDs, and a non-regressing scored common
    set can hide models whose scoring failed. Both remain explicit evidence.
    """
    all_new = Counter(r["ref_mrna_id"] for r in new_rows)
    all_old = Counter(r["ref_mrna_id"] for r in old_rows)
    expected = set(expected_ids) if expected_ids is not None else set(all_old)
    coding_ids = (set(expected_coding_ids) if expected_coding_ids is not None else
                  {r['ref_mrna_id'] for r in old_rows if str(r.get('is_coding', '')).lower() in ('1', 'true')})
    coding_consistent = all((str(r.get('is_coding', '')).lower() in ('1', 'true')) ==
                            (r['ref_mrna_id'] in coding_ids) for r in chain(new_rows, old_rows))
    nr, oldr = _recovered_rows(new_rows), _recovered_rows(old_rows)
    ni, oi = Counter(r["ref_mrna_id"] for r in nr), Counter(r["ref_mrna_id"] for r in oldr)
    new, old = _recovered_pi(nr), _recovered_pi(oldr)
    shared = sorted(set(new) & set(old))
    regressed = [k for k in shared if new[k] < old[k] - 1e-9]
    mean_old = math.fsum(old[k] for k in shared) / len(shared) if shared else None
    mean_new = math.fsum(new[k] for k in shared) / len(shared) if shared else None
    new_missing, old_missing = _missing_cds_ids(nr), _missing_cds_ids(oldr)
    return {
        "evaluation_rows_complete": bool(expected) and all_new.keys() == all_old.keys() == expected,
        "expected_ids": sorted(expected), "expected_coding_ids": sorted(coding_ids),
        "coding_membership_consistent": coding_consistent,
        "missing_evaluation_new_ids": sorted(expected - all_new.keys()),
        "missing_evaluation_old_ids": sorted(expected - all_old.keys()),
        "unexpected_evaluation_new_ids": sorted(all_new.keys() - expected),
        "unexpected_evaluation_old_ids": sorted(all_old.keys() - expected),
        "n_reference_rows_new": sum(all_new.values()), "n_reference_rows_old": sum(all_old.values()),
        "n": len(shared), "shared_recovered": len(ni.keys() & oi.keys()),
        "mean_old": mean_old, "mean_new": mean_new,
        "delta": mean_new - mean_old if shared else None,
        "n_improved": sum(new[k] > old[k] + 1e-9 for k in shared),
        "n_regressed": len(regressed), "regressed_examples": regressed[:10],
        "regressed_ids": regressed,
        "lost_ids": sorted(oi.keys() - ni.keys()),
        "added_ids": sorted(ni.keys() - oi.keys()),
        "duplicate_new_ids": sorted(k for k, count in all_new.items() if count > 1),
        "duplicate_old_ids": sorted(k for k, count in all_old.items() if count > 1),
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


def validate_output(path, *, reference=False):
    """Retain every gate-relevant ERROR, count and sample nonerror diagnostics."""
    from lifton.gff3_validator import validate_gff3_file
    result = validate_gff3_file(str(path), max_issues_per_check=sys.maxsize,
                               check_lifton_attrs=not reference)
    counts = Counter((i.severity, i.check) for i in result.issues)
    sampled = Counter()
    issues = []
    for issue in result.issues:
        key = (issue.severity, issue.check)
        if issue.severity != "ERROR" and sampled[key] >= 5:
            continue
        sampled[key] += 1
        issues.append({"severity": issue.severity, "check": issue.check, "feature_id": issue.feature_id,
                       "line": issue.lineno, "message": issue.message})
    fatal = {"parse_error", "file_exists", "file_readable", "features_present"}
    return {"exit": 0 if result.is_valid else 1, "valid": result.is_valid,
            "n_errors": sum(count for (severity, _), count in counts.items() if severity == "ERROR"),
            "n_warnings": sum(count for (severity, _), count in counts.items() if severity == "WARNING"),
            "issue_counts": [{"severity": severity, "check": check, "count": count}
                             for (severity, check), count in sorted(counts.items())],
            "nonerror_examples_per_check": 5,
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
    # Historical name means ID presence, not CDS or scored recovery.
    summary["n_recovered_coding"] = len(recovered)
    summary["n_present_coding"] = len(recovered)
    summary["n_cds_recovered_coding"] = sum(
        str(r.get("n_cds_lifted", "")).isdigit() and int(r["n_cds_lifted"]) > 0 for r in recovered)
    summary["n_cds_evidence_unknown"] = sum(
        not str(r.get("n_cds_lifted", "")).isdigit() for r in recovered)
    summary["n_unresolved_coding"] = len(recovered) - len(identities) - len(_missing_cds_ids(recovered))
    summary["n_scored_coding"] = len(identities)
    summary["n_missing_cds"] = len(_missing_cds_ids(recovered))
    summary["mean_protein_identity"] = (
        round(sum(identities) / len(identities), 5) if identities else None)
    return summary


def _evaluation_evidence():
    return provenance.evaluation_evidence()


def reference_inventory(ref):
    """Keep annotated coding membership separate from extractable proteins."""
    coding = {identifier for identifier, value in ref.items() if value['is_coding']}
    return {'schema_version': 1, 'all_ids': sorted(ref), 'coding_ids': sorted(coding),
            'unresolved_coding_ids': sorted(identifier for identifier in coding if not ref[identifier]['prot'])}


def _isolated_inputs(paths, destination):
    """Build all evaluator indexes beside private links, never source artifacts."""
    destination.mkdir(parents=True, exist_ok=False)
    result = {}
    for key, source in paths.items():
        source = Path(source).resolve(strict=True)
        target = destination / (key + "".join(source.suffixes))
        target.symlink_to(source)
        result[key] = target
    return result


def run_cell(bid, config, log=print):
    if (config.paths is None) != (config.metadata is None):
        raise ValueError('Explicit qualification inputs require matching benchmark metadata')
    bench = config.metadata if config.metadata is not None else fc._bench(bid)
    paths = config.paths if config.paths is not None else fc._full_paths(bid)
    anndb = bench.get("annotation_database", "RefSeq")
    root = config.output / "cells" / bid
    arms = {}
    for arm, tree in ((CANDIDATE, config.candidate), (REFERENCE, config.baseline)):
        log(f"=== {bid}: {arm}, fresh alignment, -t {config.threads}, copies={config.copies} ===")
        arms[arm] = _run_arm(bid, arm, tree, paths, anndb, config, log)
    manifest = {"id": bid, "species": bench["species"], "cross_species": bench["cross_species"],
                "miniprot_target_space": "transcript", "protein_acc_to_mrna": {},
                "paths": {k: str(v) for k, v in paths.items()}}
    root.mkdir(parents=True, exist_ok=True)
    eval_dir = Path(tempfile.mkdtemp(prefix="eval-", dir=root))
    evaluation_paths = _isolated_inputs(paths, eval_dir / "inputs")
    inputs = {key: provenance.fingerprint(value) for key, value in paths.items()}
    evaluation_evidence = _evaluation_evidence()
    ref, ref_index = evaluator.build_reference(
        str(evaluation_paths["ref_gff"]), str(evaluation_paths["ref_fa"]), log=log)
    genes = gene_level.load_transcript_genes(str(paths["ref_gff"]))
    record = {"schema_version": REPORT_SCHEMA_VERSION, "cell": bid, "threads": config.threads,
              "divergence": bench.get("divergence_class", ""), "arms": {},
              "protocol": {"fresh_alignment": True, "copies": config.copies,
                           "role_options": config.role_options or {}},
              "evaluator": evaluation_evidence, "inputs": inputs}
    rows_by_arm = {}
    for arm, (gff, receipt) in arms.items():
        output_artifact = provenance.fingerprint(gff)
        isolated_output = _isolated_inputs({"output": gff}, eval_dir / arm)["output"]
        evaluator.evaluate_tool(arm, str(isolated_output), str(evaluation_paths["tgt_fa"]), ref, manifest,
                                eval_dir, None, log=log, ref_index=ref_index, threads=config.threads)
        table = eval_dir / f"{arm}.transcripts.tsv"
        table_artifact = provenance.fingerprint(table)
        rows = list(gene_level.read_transcript_tsv(table))
        rows_by_arm[arm] = rows
        record["arms"][arm] = {
            "role": arm, "version": receipt["evidence"]["runtime"]["version"],
            "commit": receipt["evidence"]["source"]["commit"],
            "recall": _summarize(rows, genes), "validity": receipt["validation"],
            "wall_seconds": receipt["wall_seconds"], "peak_rss_mb": receipt["peak_rss_mb"],
            "output": str(gff), "output_artifact": output_artifact,
            "completed": True, "provenance_verified": True,
            "receipt": provenance.fingerprint(gff.parent / "completion.json"),
            "transcript_table": table_artifact,
        }
    record["reference_validity"] = validate_output(paths["ref_gff"], reference=True)
    record["expected_reference_ids"] = sorted(ref)
    record["reference_inventory"] = reference_inventory(ref)
    record["common_set"] = common_set(*(rows_by_arm[label] for label in _arm_labels(record)), expected_ids=ref,
                                      expected_coding_ids=record['reference_inventory']['coding_ids'])
    for value in record["inputs"].values():
        provenance.verify_fingerprint(value, label="Evaluation input")
    for details in record["arms"].values():
        provenance.verify_fingerprint(details["output_artifact"], label="Output artifact")
        provenance.verify_fingerprint(details["transcript_table"], label="Transcript table")
    if _evaluation_evidence() != evaluation_evidence:
        raise RuntimeError("Evaluator source or dependencies changed while scoring")
    _finish(record)
    _write_report(config.output, bid, record)
    log(f"[{bid}] lost={len(record['common_set']['lost_ids'])}; unresolved={record['unresolved_gates']}")
    return record


def _finish(record):
    """Fail closed on missing evidence, per-ID loss and new validity issues."""
    new, old = (record["arms"][label] for label in _arm_labels(record))
    common = record["common_set"]
    inventory = record.get('reference_inventory', {})
    inventory_complete = (isinstance(inventory, dict) and inventory.get('schema_version') == 1
                          and bool(inventory.get('all_ids'))
                          and inventory.get('all_ids') == common.get('expected_ids')
                          and inventory.get('coding_ids') == common.get('expected_coding_ids')
                          and set(inventory['coding_ids']).issubset(inventory['all_ids'])
                          and isinstance(inventory.get('unresolved_coding_ids'), list)
                          and set(inventory['unresolved_coding_ids']).issubset(inventory['coding_ids']))
    nv, ov = new.get("validity", {}), old.get("validity", {})
    validity_complete = all(v.get("complete") is True and v.get("exit") in (0, 1)
                            and isinstance(v.get("n_errors"), int)
                            and isinstance(v.get("issues"), list)
                            and v["n_errors"] == sum(i.get("severity") == "ERROR" for i in v["issues"])
                            for v in (nv, ov))
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
        "reference_inventory_complete": inventory_complete,
        "reference_scoring_resolved": inventory_complete and not inventory['unresolved_coding_ids'],
        "coding_membership_consistent": common.get('coding_membership_consistent') is True,
        "evaluation_rows_complete": common.get("evaluation_rows_complete") is True,
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


def _legacy_inputs(record):
    """Only consistent, complete manifests identify the legacy inputs."""
    recorded = {}
    for arm in _arm_labels(record):
        output = Path(record["arms"][arm]["output"])
        manifest_path = output.parent / "lifton_output" / "run_manifest.json"
        try:
            manifest = json.loads(manifest_path.read_text())
            inputs = {}
            for key, name in (("ref_gff", "reference_annotation"),
                              ("ref_fa", "reference_genome"), ("tgt_fa", "target_genome")):
                value = manifest["inputs"][name]
                if (value.get("fingerprint_status") != "complete" or value.get("changed_during_hash")
                        or not value.get("sha256") or not Path(value["path"]).is_absolute()):
                    raise ValueError(f"Incomplete recorded input: {name}")
                inputs[key] = {"path": value["path"], "sha256": value["sha256"], "size": value["size_bytes"]}
        except (OSError, KeyError, TypeError, ValueError) as error:
            raise ValueError(f"Missing or invalid recorded inputs in {manifest_path}: {error}") from error
        if recorded and inputs != recorded:
            raise ValueError("Inconsistent recorded inputs between legacy arm manifests")
        recorded = inputs
    for value in recorded.values():
        provenance.verify_fingerprint(value, label="Legacy input")
    return recorded


def rescore(bid, source, destination, log=print):
    """Re-derive legacy tables into a NEW report; never upgrade missing provenance."""
    source = Path(source)
    destination = Path(destination)
    if source.resolve() == destination.resolve():
        raise ValueError("Rescoring must preserve the original campaign directory")
    legacy_source = provenance.fingerprint(source / "results" / f"{bid}.json")
    record = json.loads((source / "results" / f"{bid}.json").read_text())
    record["inputs"] = _legacy_inputs(record)
    paths = {key: Path(value["path"]) for key, value in record["inputs"].items()}
    evaluation_evidence = _evaluation_evidence()
    genes = gene_level.load_transcript_genes(str(paths["ref_gff"]))
    rows_by_arm = {}
    for arm in _arm_labels(record):
        details = record["arms"][arm]
        output = Path(details["output"])
        table = (Path(details["transcript_table"]["path"]) if details.get("transcript_table")
                 else output.parent.parent / "eval" / f"{arm}.transcripts.tsv")
        table_artifact = provenance.fingerprint(table)
        output_artifact = provenance.fingerprint(output)
        rows = list(gene_level.read_transcript_tsv(table))
        rows_by_arm[arm] = rows
        details["recall"] = _summarize(rows, genes)
        details["validity"] = validate_output(output)
        details["output_artifact"] = output_artifact
        manifest_path = output.parent / "lifton_output" / "run_manifest.json"
        manifest = json.loads(manifest_path.read_text()) if manifest_path.exists() else {}
        details["completed"] = manifest.get("run", {}).get("status") == "success"
        details["provenance_verified"] = False
        details["legacy_manifest"] = provenance.fingerprint(manifest_path) if manifest else None
        details["transcript_table"] = table_artifact
    record["schema_version"] = 2
    provenance.verify_fingerprint(legacy_source, label="Legacy source")
    record["legacy_source"] = legacy_source
    record["protocol"] = {"legacy": True, "fresh_alignment": None, "copies": None}
    record["evaluator"] = evaluation_evidence
    record["reference_validity"] = validate_output(paths["ref_gff"], reference=True)
    record["common_set"] = common_set(*(rows_by_arm[label] for label in _arm_labels(record)))
    for value in record["inputs"].values():
        provenance.verify_fingerprint(value, label="Evaluation input")
    for details in record["arms"].values():
        provenance.verify_fingerprint(details["output_artifact"], label="Output artifact")
        provenance.verify_fingerprint(details["transcript_table"], label="Transcript table")
    if _evaluation_evidence() != evaluation_evidence:
        raise RuntimeError("Evaluator source or dependencies changed while scoring")
    _finish(record)
    _write_report(destination, bid, record)
    log(f"[{bid}] legacy evidence rescored; unresolved={record['unresolved_gates']}")
    return record


def _write_report(output, bid, record):
    campaign = Path(output) / "campaign.json"
    if campaign.exists():
        record["campaign"] = provenance.fingerprint(campaign)
    path = Path(output) / "results" / f"{bid}.json"
    if path.exists():
        previous = provenance.fingerprint(path)
        archive = Path(output) / "report_history" / bid / (previous["sha256"] + ".json")
        archive.parent.mkdir(parents=True, exist_ok=True)
        archive.write_bytes(path.read_bytes())
    atomic_write_json(path, record)
    atomic_write_json(Path(output) / "report_receipts" / f"{bid}.json", provenance.fingerprint(path))


def _report_errors(output, path, record, evaluation_evidence):
    errors = []
    try:
        campaign_path = output / "campaign.json"
        if record.get("campaign") != provenance.fingerprint(campaign_path):
            raise ValueError("Pinned campaign no longer matches the report")
        seal = json.loads((output / "report_receipts" / f"{record['cell']}.json").read_text())
        if seal != provenance.fingerprint(path):
            raise ValueError("Report does not match its completion fingerprint")
        if record.get("evaluator") != evaluation_evidence:
            raise ValueError("Evaluator source or dependencies no longer match report")
        if set(record.get("inputs", {})) != {"ref_gff", "ref_fa", "tgt_fa"}:
            raise ValueError("Missing report input evidence")
        for item in record["inputs"].values():
            provenance.verify_fingerprint(item, label="Report input")
        for arm in _arm_labels(record):
            details = record["arms"][arm]
            provenance.verify_fingerprint(details["output_artifact"], label="Output")
            provenance.verify_fingerprint(details["transcript_table"], label="Transcript table")
            if record.get("legacy_source"):
                provenance.verify_fingerprint(record["legacy_source"], label="Legacy record")
                provenance.verify_fingerprint(details["legacy_manifest"], label="Legacy manifest")
            else:
                receipt_path = provenance.verify_fingerprint(details["receipt"], label="Completion receipt")["path"]
                receipt = json.loads(Path(receipt_path).read_text())
                expected = json.loads((Path(receipt_path).parent / "expected.json").read_text())
                provenance.verify_run_evidence(expected)
                if record.get("schema_version") == REPORT_SCHEMA_VERSION and (
                        arm not in (CANDIDATE, REFERENCE) or details.get("role") != arm
                        or details.get("version") != expected["runtime"]["version"]
                        or details.get("commit") != expected["source"]["commit"]):
                    raise ValueError("Arm role, version or commit does not match recorded provenance")
                if (receipt.get("output") != details["output_artifact"]
                        or receipt.get("validation") != details["validity"]
                        or not provenance.read_receipt(receipt_path, expected, details["output"],
                                                       receipt["manifest"]["path"])):
                    raise ValueError("Run receipt does not verify report")
                if receipt["evidence"]["inputs"] != record["inputs"]:
                    raise ValueError("Run receipt inputs do not match report")
        attempt = output / "attempts" / f"{record['cell']}.json"
        state = json.loads(attempt.read_text())
        if (state.get("cell") != record["cell"] or state.get("status") != "success"
                or state.get("report") != seal):
            raise ValueError("Latest cell attempt did not successfully produce this report")
        history = Path(state["attempt"]).resolve(strict=True)
        if history.parent != (output / "attempt_history" / record["cell"]).resolve():
            raise ValueError("Attempt history is outside the cell history directory")
        if json.loads(history.read_text()) != state:
            raise ValueError("Latest attempt and preserved history disagree")
    except (OSError, ValueError, TypeError, KeyError, RuntimeError) as error:
        errors.append(str(error))
    return errors


def merge(output, expected_cells=None):
    output = Path(output)
    summary = output / "release_validation.json"
    if summary.exists():
        previous = provenance.fingerprint(summary)
        archive = output / "summary_history" / (previous["sha256"] + ".json")
        archive.parent.mkdir(parents=True, exist_ok=True)
        archive.write_bytes(summary.read_bytes())
    atomic_write_json(summary, {"schema_version": 2, "gate_pass": False,
                                "campaign_errors": ["Merge has not completed"]})
    (output / "release_validation.md").write_text("# LiftOn release qualification\n\nUNRESOLVED: merge has not completed.\n")
    paths = sorted((output / "results").glob("*.json"))
    records = [json.loads(p.read_text()) for p in paths]
    expected_list = list(expected_cells or [])
    expected_cells = set(expected_list)
    observed = Counter(r["cell"] for r in records)
    missing = sorted(expected_cells - observed.keys())
    unexpected = sorted(observed.keys() - expected_cells)
    duplicates = sorted(key for key, count in observed.items() if count > 1)
    campaign_errors = []
    if not expected_list:
        campaign_errors.append("An explicit nonempty expected cell set is required")
    if len(expected_list) != len(expected_cells):
        campaign_errors.append("Expected cells contain duplicates")
    if unexpected or duplicates:
        campaign_errors.append("Unexpected or duplicate result cells")
    campaign = output / "campaign.json"
    try:
        if json.loads(campaign.read_text())["expected_cells"] != sorted(expected_cells):
            campaign_errors.append("Expected cells do not match the pinned campaign")
    except (OSError, ValueError, TypeError, KeyError):
        campaign_errors.append("Pinned campaign is missing or malformed")
    evaluation_evidence = _evaluation_evidence() if records else None
    for path, record in zip(paths, records):
        if record.get("schema_version") not in (2, REPORT_SCHEMA_VERSION):
            raise ValueError("Legacy records must be rescored before merging")
        _finish(record)
        errors = _report_errors(output, path, record, evaluation_evidence)
        if errors:
            record["gate_pass"] = False
            record["unresolved_gates"].append("report_evidence_verified")
            record["evidence_errors"] = errors
    lines = ["# LiftOn release qualification", "",
             "Generated from the accompanying JSON. PASS requires complete evidence; "
             "legacy records and unexplained scoring failures remain unresolved. "
             "An ID present without CDS is not a recovered protein. Own-set means are descriptive only.", "",
             "| Cell | Reference coding IDs present | Added IDs | Lost IDs | Shared scored | PI regressions | "
             "Unscored new | New validity issues | Gate |",
             "|---|---:|---:|---:|---:|---:|---:|---:|---|"]
    for r in records:
        if r.get("schema_version") not in (2, REPORT_SCHEMA_VERSION):
            raise ValueError("Legacy records must be rescored before merging")
        c = r["common_set"]
        issues = r["delta"]["introduced_validity_issues"]
        lines.append(f"| {r['cell']} | {r['arms'][_arm_labels(r)[0]]['recall']['n_recovered_coding']} | "
                     f"{len(c['added_ids'])} | {len(c['lost_ids'])} | {c['n']} | {c['n_regressed']} | "
                     f"{len(c['unscored_new_ids'])} | {len(issues) if issues is not None else 'unknown'} | "
                     f"{'PASS' if r['gate_pass'] else ', '.join(r['unresolved_gates'])} |")
    passed = sum(r["gate_pass"] for r in records)
    lines += ["", f"{passed}/{len(expected_cells)} expected cells pass.",
              f"Missing cells: {', '.join(missing) if missing else 'none'}.",
              f"Unexpected cells: {', '.join(unexpected) if unexpected else 'none'}.",
              f"Duplicate cells: {', '.join(duplicates) if duplicates else 'none'}.",
              *campaign_errors, ""]
    output.mkdir(parents=True, exist_ok=True)
    atomic_write_json(output / "release_validation.json", {
        "schema_version": 2, "expected_cells": sorted(expected_cells), "missing_cells": missing,
        "unexpected_cells": unexpected, "duplicate_cells": duplicates, "campaign_errors": campaign_errors,
        "records": records, "gate_pass": bool(records) and not missing and not campaign_errors and passed == len(records)})
    (output / "release_validation.md").write_text("\n".join(lines))
    return 0 if records and not missing and not campaign_errors and passed == len(records) else 1


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("cells", nargs="*")
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--run-id", required=False)
    ap.add_argument("--output-root", type=Path, default=HERE / "_runs" / "release-qualification")
    ap.add_argument("--candidate-tree", type=Path, default=HERE.parents[1])
    ap.add_argument("--reference-tree", "--baseline-tree", dest="baseline_tree", type=Path)
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
    if not args.rescore and args.baseline_tree is None:
        ap.error("fresh runs require an explicit --reference-tree (alias --baseline-tree)")
    cells = args.cells or [p.stem for p in sorted((args.rescore / "results").glob("*.json"))]
    if not cells:
        ap.error("no cells found")
    if any(not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", cell) for cell in cells):
        ap.error("cell IDs must be path-safe")
    if output.exists() and not args.resume:
        ap.error("run directory exists; use a new --run-id or explicit --resume")
    if len(cells) != len(set(cells)):
        ap.error("cell IDs must be unique")
    campaign_path = output / "campaign.json"
    if args.resume and not campaign_path.is_file():
        ap.error("resume requires the original pinned campaign; use a new --run-id")
    campaign = {"expected_cells": sorted(cells), "candidate": str(args.candidate_tree.resolve()),
                "reference": str(args.baseline_tree.resolve()) if args.baseline_tree else None, "python": args.python,
                "threads": args.threads, "copies": args.copies,
                "rescore": str(args.rescore.resolve()) if args.rescore else None}
    if campaign_path.exists() and json.loads(campaign_path.read_text()) != campaign:
        ap.error("resume configuration does not match pinned campaign")
    atomic_write_json(campaign_path, campaign)
    config = Configuration(output, args.candidate_tree.resolve(),
                           args.baseline_tree.resolve() if args.baseline_tree else None,
                           args.python, args.threads, args.copies, args.resume)
    failures = []
    for bid in cells:
        attempt_path = output / "attempt_history" / bid / f"{uuid.uuid4().hex}.json"
        active = output / "attempts" / f"{bid}.json"
        state = {"cell": bid, "status": "running", "attempt": str(attempt_path)}
        atomic_write_json(attempt_path, state)
        atomic_write_json(active, state)
        try:
            if args.rescore:
                rescore(bid, args.rescore, output)
            else:
                run_cell(bid, config)
            state.update(status="success", report=provenance.fingerprint(output / "results" / f"{bid}.json"))
        except Exception as error:
            print(f"[{bid}] FAILED: {error}", file=sys.stderr, flush=True)
            failures.append(bid)
            state.update(status="failed", error=str(error))
            atomic_write_json(output / "failures" / f"{bid}.json", state)
        atomic_write_json(attempt_path, state)
        atomic_write_json(active, state)
    result = merge(output, cells)
    return 1 if failures else result


if __name__ == "__main__":
    sys.exit(main())
