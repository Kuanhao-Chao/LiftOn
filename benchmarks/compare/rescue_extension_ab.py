#!/usr/bin/env python
"""A/B harness for the v1.0.12 miniprot-only rescue extensions.

Each experiment runs LiftOn twice per cell on cached Liftoff/miniprot inputs,
changing only that extension's environment switch, and scores both outputs
with the neutral evaluator.

Experiments:

* ``coverage_gate`` (A1): the protein-coverage sub-pass
  (``LIFTON_RESCUE_COVERAGE_GATE`` 0 vs 1).
* ``isoforms`` (A2): the isoform-aware rescue (``LIFTON_RESCUE_ISOFORMS`` 0 vs
  1), measured with the coverage sub-pass on in both arms.

Cells:

* ``ladder``: the eight subset cells of the earlier rescue A/Bs, at ``-t 1``;
* ``full``: five distant whole-genome transfers, ``-t 8 --locus-pipeline``
  (byte-identical to ``-t 1`` by the scheduling contract).

Per-cell gate, as for earlier rescue promotions:

* ``n_lost == 0``: every coding reference transcript the OFF arm recovers, the
  ON arm recovers too;
* no duplicate models: no added mRNA or gene repeats an ID (copy suffix
  removed) already in the OFF output;
* every added model explained: each added rescue mRNA is either a newly
  recovered reference transcript or one the evaluator does not score (such as
  a RefSeq ``V_gene_segment`` of an immunoglobulin gene), never neither;
* ``regressed == 0``: no common transcript loses protein identity;
* ``prefix``: when the experiment only appends, the OFF output is a byte
  prefix of the ON output;
* validity: ON has no more validator errors than OFF.

Reported, not gated: gene-level recall (all, primary-assembly, GeneID), the
added set's identity and structure, and, where a released annotation of the
target exists, how often models overlap an annotated CDS on the same strand.

Arms are pinned to a frozen source tree through ``--pythonpath`` so edits made
while a run is in flight cannot leak into it.

Usage (repository root, lifton_devel environment)::

    python -m benchmarks.compare.rescue_extension_ab --experiment coverage_gate \\
        --mode ladder --pythonpath /path/to/snapshot [IDS...]
    python -m benchmarks.compare.rescue_extension_ab --experiment coverage_gate \\
        --merge
"""
from __future__ import annotations

import argparse
import bisect
import csv
import json
import os
import re
import statistics
import sys
import threading
from collections import defaultdict
from pathlib import Path

from . import devel_refresh as dr
from . import evaluator
from . import fourway_compare as fc
from . import gene_level
from . import version_compare as vc
from .profiling import run_profiled

HERE = Path(__file__).resolve().parent
WORK = HERE / "work"
REG = json.loads((HERE / "benchmarks.json").read_text())

EXPERIMENTS = {
    "coverage_gate": {
        "title": "protein-coverage rescue sub-pass (A1)",
        "off": {"LIFTON_RESCUE_COVERAGE_GATE": "0"},
        "on": {"LIFTON_RESCUE_COVERAGE_GATE": "1"},
        "tag": "rescue_gate=protein_coverage",
        "expect_prefix": True,
    },
    # Measured on top of A1: both arms run the coverage sub-pass.
    "isoforms": {
        "title": "isoform-aware rescue (A2), on top of the coverage sub-pass",
        "off": {"LIFTON_RESCUE_COVERAGE_GATE": "1", "LIFTON_RESCUE_ISOFORMS": "0"},
        "on": {"LIFTON_RESCUE_COVERAGE_GATE": "1", "LIFTON_RESCUE_ISOFORMS": "1"},
        "tag": "rescue_isoform=true",
        "expect_prefix": False,
        "expect_same_genes": True,
    },
}

LADDER = [
    "celegans_to_briggsae", "drosophila_to_anopheles", "zebrafish_to_medaka",
    "rice_to_sorghum", "t4_human_to_xenopus", "t4_human_to_chicken",
    "human_to_mouse", "drosophila",
]

# Released annotation of each full transfer's TARGET genome, when available.
FULL = {
    "human_to_zebrafish":
        "/ccb/salz3/kh.chao/lifton_benchmark_data/zebrafish_to_medaka/ref.gff",
    "t4_human_to_chicken":
        "/ccb/salz3/kh.chao/lifton_benchmark_data/chicken_to_quail/ref.gff",
    "t4_human_to_xenopus": None,
    "arabidopsis_to_rice": "/ccb/salz2/jheinz3/shared/lifton/rice/IRGSP_genomic.gff",
    "t4_drosophila_to_bee": "/ccb/salz2/jheinz3/shared/lifton/bee/HAv3.1_genomic.gff",
}


def _bench(bid):
    for bench in REG["benchmarks"]:
        if bench["id"] == bid:
            return bench
    raise KeyError(bid)


def _link(target, link):
    link.parent.mkdir(parents=True, exist_ok=True)
    if link.is_symlink() or link.exists():
        link.unlink()
    link.symlink_to(Path(target).resolve())
    return link


def _clean_sidecars(path):
    for suffix in ("_db", "_db.lock", "_db.manifest.json", ".eval_db"):
        try:
            os.unlink(str(path) + suffix)
        except FileNotFoundError:
            pass


def cell_inputs(bid, mode):
    """Return (manifest, ref_gff, ref_fa, tgt_fa, cached inputs, threads, extra)."""
    bench = _bench(bid)
    if mode == "ladder":
        manifest = json.loads(
            (WORK / bid / "subset" / "subset.manifest.json").read_text())
        paths = manifest["paths"]
        cached = {
            "-L": WORK / bid / "tools" / "liftoff" / "liftoff.gff3",
            "-M": WORK / bid / "tools" / "miniprot" / "miniprot.gff3",
        }
        return (manifest, paths["ref_gff"], paths["ref_fa"], paths["tgt_fa"],
                cached, 1, [])
    paths = fc._full_paths(bid)
    liftoff, miniprot, transcripts, proteins = dr._cached_inputs(bid)
    manifest = {
        "id": bid, "species": bench["species"],
        "cross_species": bench["cross_species"],
        "miniprot_target_space": "transcript", "protein_acc_to_mrna": {},
        "paths": {k: str(v) for k, v in paths.items()},
    }
    cached = {"-L": liftoff, "-M": miniprot, "-T": transcripts, "-P": proteins}
    return (manifest, str(paths["ref_gff"]), str(paths["ref_fa"]),
            str(paths["tgt_fa"]), cached, 8, ["--locus-pipeline"])


def _run_arm(bid, arm, experiment, inputs, root, pythonpath):
    manifest, ref_gff, ref_fa, tgt_fa, cached, threads, extra = inputs
    arm_dir = root / arm
    (arm_dir / "inputs").mkdir(parents=True, exist_ok=True)
    out = arm_dir / f"{arm}.gff3"
    argv = [vc.VERSIONS["devel"]["bin"], "-t", str(threads), "-ad",
            _bench(bid).get("annotation_database", "RefSeq"), "-g", ref_gff]
    for flag, path in cached.items():
        # A per-arm link gives each arm its own gffutils sidecar, so the two arms
        # can run concurrently without racing on one database file.
        link = _link(path, arm_dir / "inputs" / Path(path).name)
        _clean_sidecars(link)
        argv += [flag, str(link)]
    argv += [*extra, "-o", str(out), tgt_fa, ref_fa]
    env = dict(vc._env())
    if pythonpath:
        env["PYTHONPATH"] = pythonpath
    env.update(EXPERIMENTS[experiment][arm])
    result = run_profiled(argv, label=f"{experiment}_{bid}_{arm}",
                          log_dir=arm_dir / "logs", env=env, cwd=arm_dir,
                          log=print)
    if result.exit_code != 0 or not out.exists() or out.stat().st_size == 0:
        raise RuntimeError(f"{bid}: arm {arm} failed (exit {result.exit_code}); "
                           f"see {result.stderr_path}")
    return out, result


def _scored_ids(tsv):
    with open(tsv, newline="") as handle:
        return {row["ref_mrna_id"] for row in csv.DictReader(handle, delimiter="\t")}


def _coding_rows(tsv):
    rows = {}
    with open(tsv, newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if str(row.get("is_coding")).strip().lower() in ("1", "true"):
                rows[row["ref_mrna_id"]] = row
    return rows


def _pi(row):
    try:
        return float(row["protein_identity"])
    except (TypeError, ValueError, KeyError):
        return None


def _recovered_pi(rows):
    return {tx: _pi(row) for tx, row in rows.items()
            if str(row["recovered"]).strip() in ("1", "True", "true")
            and _pi(row) is not None}


def _mean(values):
    values = [v for v in values if v is not None]
    return round(statistics.fmean(values), 5) if values else None


def _fraction(rows, key):
    values = [row.get(key) for row in rows if row.get(key) not in ("", None)]
    return round(sum(int(float(v)) for v in values) / len(values), 4) if values else None


def _quantiles(values):
    values = sorted(v for v in values if v is not None)
    if not values:
        return None

    def q(p):
        return round(values[int(p * (len(values) - 1))], 3)

    return {"n": len(values), "p10": q(0.1), "median": q(0.5), "p90": q(0.9)}


def _mrna_classes(gff, tag):
    """Map each output mRNA id to experiment/rescue/lift, and collect CDS rows."""
    classes, cds = {}, defaultdict(list)
    with open(gff) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) != 9:
                continue
            if c[2] == "mRNA":
                attributes = gene_level.parse_attributes(c[8])
                mrna = attributes["ID"][0]
                if tag and tag in c[8]:
                    classes[mrna] = "experiment"
                elif "lifton_rescue=miniprot_only" in c[8]:
                    classes[mrna] = "rescue"
                else:
                    classes[mrna] = "lift"
            elif c[2] == "CDS":
                attributes = gene_level.parse_attributes(c[8])
                for parent in attributes.get("Parent", ()):
                    cds[parent].append((c[0], c[6], int(c[3]), int(c[4])))
    return classes, cds


def _cds_index(annotation):
    """Merged same-strand CDS intervals of a released annotation."""
    raw = defaultdict(list)
    with open(annotation) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            c = line.split("\t", 8)
            if len(c) == 9 and c[2] == "CDS":
                raw[(c[0], c[6])].append((int(c[3]), int(c[4])))
    index = {}
    for key, intervals in raw.items():
        intervals.sort()
        merged = []
        for start, end in intervals:
            if merged and start <= merged[-1][1] + 1:
                merged[-1][1] = max(merged[-1][1], end)
            else:
                merged.append([start, end])
        index[key] = ([s for s, _ in merged], [e for _, e in merged])
    return index


def _overlaps(index, seqid, strand, start, end):
    starts, ends = index.get((seqid, strand), ((), ()))
    i = bisect.bisect_right(starts, end) - 1
    return i >= 0 and ends[i] >= start


def target_agreement(gff, tag, annotation):
    """Fraction of coding models per class whose CDS overlaps an annotated CDS."""
    if not annotation or not Path(annotation).exists():
        return None
    index = _cds_index(annotation)
    classes, cds = _mrna_classes(gff, tag)
    counts = defaultdict(lambda: [0, 0])
    for mrna, cls in classes.items():
        segments = cds.get(mrna)
        if not segments:
            continue
        counts[cls][1] += 1
        if any(_overlaps(index, *segment) for segment in segments):
            counts[cls][0] += 1
    return {cls: {"n": n, "agree": round(hit / n, 4) if n else None}
            for cls, (hit, n) in sorted(counts.items())}


def _gene_ids(gff):
    ids = set()
    with open(gff) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            c = line.split("\t", 8)
            if len(c) == 9 and c[2] == "gene":
                ids.add(gene_level.parse_attributes(c[8])["ID"][0])
    return ids


def _tagged_transcript_ids(gff, needle):
    ids = set()
    with open(gff) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            c = line.split("\t", 8)
            if len(c) == 9 and c[2] == "mRNA" and needle in c[8]:
                ids.add(gene_level.parse_attributes(c[8])["ID"][0])
    return ids


def _ids_by_type(gff, featuretype):
    ids = set()
    with open(gff) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            c = line.split("\t", 8)
            if len(c) == 9 and c[2] == featuretype:
                ids.add(gene_level.parse_attributes(c[8])["ID"][0])
    return ids


def redundancy(out_off, out_on, scored_ids, n_new):
    """Classify what the ON arm added beyond the OFF arm.

    A model is redundant only if it duplicates something already present: an
    added mRNA whose ID (copy suffix removed) is already an mRNA of the OFF
    output, or an added gene whose ID is already a gene of the OFF output.
    Added rescue models whose reference transcript the evaluator does not
    score (for example RefSeq ``V_gene_segment`` containers of immunoglobulin
    genes) are counted separately; they are new but invisible to the scorer.
    Whatever is neither new, redundant, nor unscored is unexplained.
    """
    def base(identifier):
        return re.sub(r"_\d+$", "", identifier)

    off_mrna = {base(i) for i in _ids_by_type(out_off, "mRNA")}
    off_genes = {base(i) for i in _ids_by_type(out_off, "gene")}
    added_mrna = (_tagged_transcript_ids(out_on, "lifton_rescue=miniprot_only")
                  - _tagged_transcript_ids(out_off, "lifton_rescue=miniprot_only"))
    added_genes = _ids_by_type(out_on, "gene") - _ids_by_type(out_off, "gene")
    duplicate_mrna = sorted(i for i in added_mrna if base(i) in off_mrna)
    duplicate_genes = sorted(i for i in added_genes if base(i) in off_genes)
    unscored = sorted(i for i in added_mrna if base(i) not in scored_ids)
    return {
        "n_added_rescue_mrna": len(added_mrna),
        "n_duplicate_mrna": len(duplicate_mrna),
        "n_duplicate_genes": len(duplicate_genes),
        "n_unscored_added": len(unscored),
        "unscored_examples": unscored[:10],
        "duplicate_examples": (duplicate_mrna + duplicate_genes)[:10],
        "n_unexplained": len(added_mrna) - n_new - len(unscored),
    }


def _structure(rows):
    """Identity and structure summary of a set of evaluator rows."""
    pis = [_pi(row) for row in rows]
    pis = [p for p in pis if p is not None]
    return {
        "n": len(rows),
        "mean_pi": _mean(pis),
        "pi_ge_0.5": round(sum(p >= 0.5 for p in pis) / len(pis), 4) if pis else None,
        "orf_valid": _fraction(rows, "orf_valid"),
        "intron_chain_exact": _fraction(rows, "intron_chain_exact"),
    }


def _count_tagged(gff, needle):
    total = 0
    with open(gff) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            c = line.split("\t", 8)
            if len(c) == 9 and c[2] == "mRNA" and needle in c[8]:
                total += 1
    return total


def run_cell(bid, mode, experiment, pythonpath, force=False):
    spec = EXPERIMENTS[experiment]
    root = WORK / bid / f"_v1012_{experiment}_ab_{mode}"
    result_path = root / "result.json"
    if result_path.exists() and not force:
        print(f"[{bid}] cached result {result_path}")
        return json.loads(result_path.read_text())
    root.mkdir(parents=True, exist_ok=True)
    inputs = cell_inputs(bid, mode)
    manifest, ref_gff, ref_fa, tgt_fa = inputs[:4]

    outputs, errors = {}, {}

    def _worker(arm):
        try:
            outputs[arm] = _run_arm(bid, arm, experiment, inputs, root, pythonpath)
        except Exception as error:  # surfaced after both arms finish
            errors[arm] = error

    threads = [threading.Thread(target=_worker, args=(arm,)) for arm in ("off", "on")]
    for thread in threads:
        thread.start()
    for thread in threads:
        thread.join()
    if errors:
        raise RuntimeError(f"{bid}: {errors}")
    out_off, prof_off = outputs["off"]
    out_on, prof_on = outputs["on"]

    # Score through a private link so the evaluator's sidecar cleanup never
    # touches the reference database cache LiftOn itself uses.
    eval_ref = _link(ref_gff, root / "eval_ref" / Path(ref_gff).name)
    ref, ref_index = evaluator.build_reference(str(eval_ref), ref_fa, log=print)
    eval_dir = root / "eval"
    eval_dir.mkdir(parents=True, exist_ok=True)
    for arm, out in (("off", out_off), ("on", out_on)):
        evaluator.evaluate_tool(arm, str(out), tgt_fa, ref, manifest, eval_dir,
                                None, log=print, ref_index=ref_index, threads=8)
    rows_off = _coding_rows(eval_dir / "off.transcripts.tsv")
    rows_on = _coding_rows(eval_dir / "on.transcripts.tsv")
    pi_off, pi_on = _recovered_pi(rows_off), _recovered_pi(rows_on)
    new = sorted(set(pi_on) - set(pi_off))
    lost = sorted(set(pi_off) - set(pi_on))
    common = set(pi_off) & set(pi_on)
    regressed = sorted(tx for tx in common if pi_on[tx] < pi_off[tx] - 1e-9)
    tagged_off = _count_tagged(out_off, "lifton_rescue=miniprot_only")
    tagged_on = _count_tagged(out_on, "lifton_rescue=miniprot_only")
    experiment_tagged = _count_tagged(out_on, spec["tag"])
    added_models = redundancy(out_off, out_on,
                              _scored_ids(eval_dir / "on.transcripts.tsv"), len(new))
    n_redundant = added_models["n_duplicate_mrna"] + added_models["n_duplicate_genes"]
    prefix = out_on.read_bytes().startswith(out_off.read_bytes())
    validity = {"off": vc.validate_gff(out_off, log=print),
                "on": vc.validate_gff(out_on, log=print)}

    genes = gene_level.load_transcript_genes(eval_ref)
    recall = {arm: gene_level.summarize_tsv(eval_dir / f"{arm}.transcripts.tsv", genes)
              for arm in ("off", "on")}
    added_rows = [rows_on[tx] for tx in new]
    added = {
        "n": len(new),
        "pi": _quantiles(pi_on[tx] for tx in new),
        "mean_pi": _mean(pi_on[tx] for tx in new),
        "pi_ge_0.5": round(sum(pi_on[tx] >= 0.5 for tx in new) / len(new), 4)
        if new else None,
        "orf_valid": _fraction(added_rows, "orf_valid"),
        "intron_chain_exact": _fraction(added_rows, "intron_chain_exact"),
        "exon_sn": _mean(float(r["exon_sn"]) for r in added_rows if r["exon_sn"] not in ("", None)),
        "exon_sp": _mean(float(r["exon_sp"]) for r in added_rows if r["exon_sp"] not in ("", None)),
    }
    # The earlier rescues (miniprot-only models already in the OFF output) are
    # the natural comparison population for anything the experiment adds.
    earlier = _tagged_transcript_ids(out_off, "lifton_rescue=miniprot_only")
    earlier_rows = [rows_off[re.sub(r"_\d+$", "", tx)] for tx in sorted(earlier)
                    if re.sub(r"_\d+$", "", tx) in rows_off]
    same_genes = _gene_ids(out_off) == _gene_ids(out_on)
    annotation = FULL.get(bid) if mode == "full" else None
    agreement = target_agreement(out_on, spec["tag"], annotation)

    errors_off = validity["off"].get("n_errors") or 0
    errors_on = validity["on"].get("n_errors") or 0
    gate = {
        "n_lost_zero": not lost,
        "no_duplicate_models": n_redundant == 0,
        "all_added_models_explained": added_models["n_unexplained"] == 0,
        "no_regression": not regressed,
        "validity_not_worse": errors_on <= errors_off,
    }
    if spec.get("expect_prefix"):
        gate["off_is_prefix_of_on"] = prefix
    if spec.get("expect_same_genes"):
        gate["gene_ids_unchanged"] = same_genes
    result = {
        "benchmark": bid, "mode": mode, "experiment": experiment,
        "divergence": _bench(bid).get("divergence_class", ""),
        "n_recovered_off": len(pi_off), "n_recovered_on": len(pi_on),
        "n_new": len(new), "n_lost": len(lost), "n_regressed": len(regressed),
        "lost_examples": lost[:10], "regressed_examples": regressed[:10],
        "tagged_rescue_off": tagged_off, "tagged_rescue_on": tagged_on,
        "experiment_tagged": experiment_tagged, "n_redundant": n_redundant,
        "added_models": added_models,
        "prefix": prefix, "gene_ids_unchanged": same_genes, "validity": validity,
        "recall": recall, "added": added,
        "earlier_rescues": _structure(earlier_rows),
        "target_agreement": agreement,
        "wall_seconds": {"off": prof_off.wall_clock_seconds,
                         "on": prof_on.wall_clock_seconds},
        "peak_rss_mb": {"off": prof_off.peak_rss_mb, "on": prof_on.peak_rss_mb},
        "gate": gate, "gate_pass": all(gate.values()),
    }
    result_path.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(f"[{bid}] new={len(new)} lost={len(lost)} regressed={len(regressed)} "
          f"redundant={n_redundant} prefix={prefix} val {errors_off}->{errors_on} "
          f"gate {'PASS' if result['gate_pass'] else 'FAIL'}", flush=True)
    return result


def _backfill(root, result):
    """Add blocks introduced after a cell ran, from its outputs and TSVs."""
    out_off, out_on = root / "off" / "off.gff3", root / "on" / "on.gff3"
    if "earlier_rescues" not in result:
        rows_off = _coding_rows(root / "eval" / "off.transcripts.tsv")
        earlier = _tagged_transcript_ids(out_off, "lifton_rescue=miniprot_only")
        result["earlier_rescues"] = _structure(
            [rows_off[re.sub(r"_\d+$", "", tx)] for tx in sorted(earlier)
             if re.sub(r"_\d+$", "", tx) in rows_off])
    if "gene_ids_unchanged" not in result:
        result["gene_ids_unchanged"] = _gene_ids(out_off) == _gene_ids(out_on)
    # Re-judge redundancy with the direct duplicate checks (cells that ran
    # before them used "extra tagged minus new", which also counted models
    # the evaluator cannot score).
    added = redundancy(out_off, out_on, _scored_ids(root / "eval" / "on.transcripts.tsv"),
                       result["n_new"])
    result["added_models"] = added
    result["n_redundant"] = added["n_duplicate_mrna"] + added["n_duplicate_genes"]
    gate = result["gate"]
    gate.pop("n_redundant_zero", None)
    gate["no_duplicate_models"] = result["n_redundant"] == 0
    gate["all_added_models_explained"] = added["n_unexplained"] == 0
    result["gate_pass"] = all(gate.values())
    return result


def merge(experiment):
    results = []
    for mode, ids in (("ladder", LADDER), ("full", list(FULL))):
        for bid in ids:
            root = WORK / bid / f"_v1012_{experiment}_ab_{mode}"
            path = root / "result.json"
            if path.exists():
                result = _backfill(root, json.loads(path.read_text()))
                path.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
                results.append(result)
    json_path = HERE / f"rescue_extension_ab.{experiment}.json"
    json_path.write_text(json.dumps(results, indent=2, sort_keys=True) + "\n")
    _write_markdown(experiment, results, json_path.with_suffix(".md"))
    passed = sum(r["gate_pass"] for r in results)
    print(f"GATE: {passed}/{len(results)} cells pass -> {json_path} (+ .md)")
    return 0 if passed == len(results) else 1


def _recall_cell(result, arm, key):
    block = result["recall"][arm]
    if key == "primary_gene":
        return block["primary"]["gene_recall_coding"]
    if key == "geneid":
        return block["geneid_collapsed"]["recall"]
    return block[key]


def _write_markdown(experiment, results, path):
    spec = EXPERIMENTS[experiment]
    lines = [
        f"## Rescue-extension A/B: {spec['title']}\n",
        "Two arms per cell on cached Liftoff/miniprot inputs; only the "
        f"experiment switch differs (`{spec['off']}` vs `{spec['on']}`). Scored "
        "by the neutral evaluator. Gate: 0 lost, 0 redundant, 0 regressed, "
        "validity not worse" + (", OFF output a byte prefix of ON" if spec.get("expect_prefix") else "")
        + ".\n",
        "| cell | mode | new tx | added mean PI | added ORF-valid | lost | regr | redundant "
        "| gene recall off→on | primary gene off→on | tx recall off→on | val off→on | gate |",
        "|---|---|---|---|---|---|---|---|---|---|---|---|---|",
    ]
    for r in results:
        v = r["validity"]
        lines.append(
            "| {} | {} | {} | {} | {} | {} | {} | {} | {}→{} | {}→{} | {}→{} | {}→{} | {} |".format(
                r["benchmark"], r["mode"], r["n_new"], r["added"]["mean_pi"],
                r["added"]["orf_valid"], r["n_lost"], r["n_regressed"],
                r["n_redundant"],
                _recall_cell(r, "off", "gene_recall_coding"),
                _recall_cell(r, "on", "gene_recall_coding"),
                _recall_cell(r, "off", "primary_gene"),
                _recall_cell(r, "on", "primary_gene"),
                _recall_cell(r, "off", "transcript_recall_coding"),
                _recall_cell(r, "on", "transcript_recall_coding"),
                v["off"].get("n_errors"), v["on"].get("n_errors"),
                "**PASS**" if r["gate_pass"] else "FAIL"))
    lines += ["", "### Added models against the earlier rescues", "",
              "Earlier rescues are the miniprot-only models already in the OFF "
              "output, the natural comparison for what the experiment adds.", "",
              "| cell | added n | added mean PI | added PI≥0.5 | added ORF-valid "
              "| earlier n | earlier mean PI | earlier PI≥0.5 | earlier ORF-valid |",
              "|---|---|---|---|---|---|---|---|---|"]
    for r in results:
        a, e = r["added"], r.get("earlier_rescues") or {}
        lines.append(f"| {r['benchmark']} ({r['mode']}) | {a['n']} | {a['mean_pi']} | "
                     f"{a['pi_ge_0.5']} | {a['orf_valid']} | {e.get('n')} | "
                     f"{e.get('mean_pi')} | {e.get('pi_ge_0.5')} | {e.get('orf_valid')} |")
    agreement_rows = [r for r in results if r.get("target_agreement")]
    if agreement_rows:
        lines += ["", "### Agreement with the released target annotation", "",
                  "Fraction of coding models whose CDS overlaps an annotated CDS on "
                  "the same strand, by model class in the ON output.", "",
                  "| cell | experiment-added | earlier rescues | DNA lift |",
                  "|---|---|---|---|"]
        for r in agreement_rows:
            a = r["target_agreement"]

            def cell(name):
                block = a.get(name)
                return f"{block['agree']} (n={block['n']})" if block else "—"

            lines.append(f"| {r['benchmark']} | {cell('experiment')} | "
                         f"{cell('rescue')} | {cell('lift')} |")
    path.write_text("\n".join(lines) + "\n")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("ids", nargs="*")
    parser.add_argument("--experiment", choices=sorted(EXPERIMENTS), required=True)
    parser.add_argument("--mode", choices=("ladder", "full"), default="ladder")
    parser.add_argument("--pythonpath", default=None,
                        help="frozen source tree the LiftOn arms import")
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--merge", action="store_true")
    args = parser.parse_args(argv)
    if args.merge:
        return merge(args.experiment)
    ids = args.ids or (LADDER if args.mode == "ladder" else list(FULL))
    failures = []
    for bid in ids:
        try:
            run_cell(bid, args.mode, args.experiment, args.pythonpath, args.force)
        except Exception as error:
            print(f"[{bid}] FAILED: {error}", file=sys.stderr, flush=True)
            failures.append(bid)
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
