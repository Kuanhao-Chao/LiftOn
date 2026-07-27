"""extcompare.py — broader-comparison sub-study (chromosome-subset scope).

Scores the added coordinate LiftOver comparator alongside the headline lift-over
tools (Liftoff, LiftOn v1.0.8, v1.0.9) on the 5 representative subset pairs,
with miniprot kept as protein evidence. CAT was attempted but its toil
sub-workflows proved intractable on the RefSeq subset (documented separately) —
no CAT column here.

LiftOver (CrossMap, coordinate-lift, no frame correction) is reported with a
DUAL metric: (a) re-scored protein identity via the same parasail kernel
(expected LOW — out-of-class) and (b) a feature-lift rate (fraction of reference
coding transcripts whose CDS lifted at all, by id). CrossMap fragments features
at chain boundaries, so we collapse the lifted GFF to one record per reference
id (the dominant target seqid) before scoring, and count the lift rate by id.

Run:  python -m benchmarks.compare.extcompare            # all 5 pairs
      python -m benchmarks.compare.extcompare t2_human_to_gorilla
Writes benchmarks/compare/extcompare_results.json (NEW file; never touches
the committed fourway_results.json).
"""
from __future__ import annotations

import collections
import json
import re
import sys
from pathlib import Path

from . import evaluator, version_compare as vc

HERE = Path(__file__).resolve().parent
WORK = HERE / "work"
EXT = Path("/ccb/salz3/kh.chao/lifton_extcompare/work_ext")
OUT_JSON = HERE / "extcompare_results.json"
REG = {b["id"]: b for b in json.loads((HERE / "benchmarks.json").read_text())["benchmarks"]}

PAIRS = ["t2_human_to_gorilla", "t1_maize_b73_to_mo17", "t3_dog_to_cat",
         "arabidopsis_to_rice", "t4_human_to_chicken"]
# scored set: headline lift-over tools + LiftOver (added) + miniprot (evidence)
TOOLS = ["liftoff", "lifton_stable", "lifton_devel", "liftover", "miniprot"]

_ATTR = re.compile(r'(?:^|;)(ID|Parent)=([^;]+)')


def _attrs(s):
    return {k: v for k, v in _ATTR.findall(s)}


def _collapse_liftover(raw_gff: Path, out_gff: Path, ref_ids: set) -> dict:
    """CrossMap fragments each feature into many lines sharing one ID (it splits a
    feature at every chain block boundary, even on the same seqid). Collapse to ONE
    clean transcript per reference id: gather all lifted CDS segments for a
    transcript, pick the target seqid carrying the most lifted CDS bases, dedupe the
    segments there, and emit a single synthetic mRNA (spanning them) + those CDS.
    The evaluator then translates the concatenated (frame-UNcorrected) CDS — which
    is the whole point: most will not re-translate to the reference protein. Returns
    the feature-lift rate (fraction of reference coding transcripts whose CDS lifted
    at all, by id)."""
    cds_by_id_seq = collections.defaultdict(lambda: collections.defaultdict(list))  # id -> seqid -> [(s,e,strand,frame)]
    for ln in open(raw_gff):
        if ln.startswith("#"):
            continue
        f = ln.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "CDS":
            continue
        for p in _attrs(f[8]).get("Parent", "").split(","):
            if p:
                cds_by_id_seq[p][f[0]].append((int(f[3]), int(f[4]), f[6], f[7]))
    lifted_ids = set()
    with open(out_gff, "w") as fh:
        fh.write("##gff-version 3\n")
        for tid, seqmap in cds_by_id_seq.items():
            best_seq = max(seqmap, key=lambda s: sum(e - st + 1 for st, e, *_ in seqmap[s]))
            segs = sorted(set(seqmap[best_seq]))  # dedupe identical fragments
            if not segs:
                continue
            lifted_ids.add(tid)
            strand = segs[0][2]
            start, end = min(s for s, *_ in segs), max(e for _, e, *_ in segs)
            gid = "gene-" + tid
            fh.write(f"{best_seq}\tliftover\tgene\t{start}\t{end}\t.\t{strand}\t.\tID={gid}\n")
            fh.write(f"{best_seq}\tliftover\tmRNA\t{start}\t{end}\t.\t{strand}\t.\tID={tid};Parent={gid}\n")
            for i, (s, e, st, fr) in enumerate(segs):
                fh.write(f"{best_seq}\tliftover\texon\t{s}\t{e}\t.\t{st}\t.\tID={tid}-e{i};Parent={tid}\n")
                fh.write(f"{best_seq}\tliftover\tCDS\t{s}\t{e}\t.\t{st}\t{fr}\tID={tid}-c{i};Parent={tid}\n")
    n_ref_coding = len(ref_ids)
    n_lifted = len(lifted_ids & ref_ids)
    return {"cds_lift_rate": round(n_lifted / n_ref_coding, 5) if n_ref_coding else None,
            "n_ref_coding_for_lift": n_ref_coding, "n_lifted_coding": n_lifted}


def _gff_paths(pid):
    w = WORK / pid
    return {
        "liftoff": w / "tools" / "liftoff" / "liftoff.gff3",
        "miniprot": w / "tools" / "miniprot" / "miniprot.gff3",
        "lifton_stable": w / "_fourway" / "lifton_stable" / "stable.gff3",
        "lifton_devel": w / "_fourway" / "lifton_devel" / "devel.gff3",
        "liftover": EXT / pid / "liftover" / "liftover.collapsed.gff3",
    }


def run_pair(pid, log=print):
    bench = REG[pid]
    man = json.loads((WORK / pid / "subset" / "subset.manifest.json").read_text())
    paths = man["paths"]
    log(f"\n=== {pid} (subset scope) ===")
    ref, ref_index = evaluator.build_reference(str(paths["ref_gff"]), str(paths["ref_fa"]), log=log)
    ref_coding_ids = {k for k, v in ref.items() if v["is_coding"]}

    # prepare the collapsed LiftOver GFF + feature-lift rate
    lift_raw = EXT / pid / "liftover" / "liftover.gff3"
    lift_info = None
    gffs = _gff_paths(pid)
    if lift_raw.exists():
        lift_info = _collapse_liftover(lift_raw, gffs["liftover"], ref_coding_ids)
        log(f"  [liftover] collapsed; cds_lift_rate={lift_info['cds_lift_rate']}")

    eval_dir = EXT / pid / "extcompare_eval"
    eval_dir.mkdir(parents=True, exist_ok=True)
    summaries, validity = {}, {}
    for tool in TOOLS:
        gff = gffs.get(tool)
        if gff is None or not Path(gff).exists() or Path(gff).stat().st_size == 0:
            log(f"  skip {tool}: missing gff ({gff})")
            continue
        try:
            summaries[tool] = evaluator.evaluate_tool(
                tool, str(gff), str(paths["tgt_fa"]), ref, man, eval_dir,
                None, log=log, ref_index=ref_index, threads=8)
            validity[tool] = vc.validate_gff(str(gff), log=log)
        except Exception as e:  # noqa: BLE001
            log(f"  !! {tool} scoring failed on {pid}: {e}")
    present = [t for t in TOOLS if t in summaries]
    rec = {
        "benchmark": pid, "mode": "subset_extcompare", "key": f"extcompare:{pid}",
        "species": bench["species"], "cross_species": bench["cross_species"],
        "divergence_class": ("same_species" if not bench["cross_species"]
                             else bench.get("divergence_class", "cross_species")),
        "scope": "chromosome_subset",
        "n_reference_coding": summaries[present[0]]["n_reference_coding"] if present else None,
        "tools": present,
        "evidence_tools": ["miniprot"],
        "completeness_coding": {t: summaries[t]["completeness_coding"] for t in present},
        "mean_pi": {t: summaries[t]["protein_identity"]["mean"] for t in present},
        "median_pi": {t: summaries[t]["protein_identity"]["median"] for t in present},
        "n_recovered_coding": {t: summaries[t]["n_recovered_coding"] for t in present},
        "validity": validity,
        "liftover_feature_lift_rate": lift_info,
        "cat": {"status": "attempted_intractable",
                "note": "CAT v2 installs + runs 18/21 tasks but its toil sub-workflows "
                        "(Chaining/GenerateHints) fail opaquely on the RefSeq subset; no column."},
    }
    return rec


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    pairs = argv if argv else PAIRS
    results = {}
    if OUT_JSON.exists():
        results = json.loads(OUT_JSON.read_text())
    for pid in pairs:
        try:
            rec = run_pair(pid)
            results[rec["key"]] = rec
            OUT_JSON.write_text(json.dumps(results, indent=2))
            print(f"  -> wrote cell extcompare:{pid}")
        except Exception as e:  # noqa: BLE001
            import traceback
            print(f"!! {pid} FAILED: {e}\n{traceback.format_exc()}")
    print(f"\nextcompare_results.json: {len(results)} cells")


if __name__ == "__main__":
    main()
