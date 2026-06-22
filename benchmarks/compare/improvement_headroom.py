#!/usr/bin/env python3
"""
improvement_headroom.py  --  READ-ONLY brainstorm/diagnostic (no engine code).

Quantifies the upper-bound headroom of four proposed LiftOn-devel improvements,
from EXISTING 4-way benchmark artifacts (per-transcript eval TSVs + tool output
GFFs). NOTHING here changes the engine or any output; it only measures.

Ideas (see plan / CLAUDE.md regime discussion):
  A  miniprot-only rescue (completeness): high-PI miniprot genes devel LOSES,
     split genuine-new vs duplicate by REAL interval overlap against devel output.
  B  3-way best-of-outcome merge (accuracy): devel-recovered transcripts where a
     standalone miniprot model would score materially higher (the broken-stub cases).
  C  map_failed recovery: devel transcripts whose region lifted but protein didn't
     map; how many miniprot recovers at high PI (a recoverability upper bound).
  D  splice-junction rebuild (lifton2 Strategy B): estimated from B's residual +
     the lifton2-measured increment (qualitative; higher-risk refactor).

Usage (repo root):
  PYTHONNOUSERSITE=1 /home/kh.chao/.../lifton_devel/bin/python \
      benchmarks/compare/improvement_headroom.py
Writes improvement_headroom.json + improvement_headroom.md next to this file.
"""
import csv, json, os, glob, statistics
from pathlib import Path

HERE = Path(__file__).resolve().parent
WORK = HERE / "work"

# divergence_class: from benchmarks.json/tiers.json, with a fallback for legacy ids.
FALLBACK_DIV = {
    "drosophila": "same_species", "mouse": "same_species", "rice": "same_species",
    "bee": "same_species", "arabidopsis": "same_species", "human_mane": "same_species",
    "mouse_to_rat": "close_cross_species", "human_to_chimp": "close_cross_species",
}
DIV_ORDER = ["same_species", "close_cross_species", "distant_cross_species",
             "very_distant_cross_species"]


def load_div():
    m = {}
    for f in ["benchmarks.json", "../tiers.json", "tiers.json"]:
        p = HERE / f
        if not p.exists():
            continue
        try:
            d = json.loads(p.read_text())
        except Exception:
            continue
        items = d if isinstance(d, list) else d.get("benchmarks", d.get("datasets", []))
        if isinstance(items, dict):
            items = list(items.values())
        for it in items:
            if isinstance(it, dict) and it.get("id"):
                dc = it.get("divergence_class")
                if dc:
                    m[it["id"]] = dc
    return m


def eval_dir(cid):
    """Prefer the FULL-genome eval dir; fall back to the subset."""
    full = WORK / cid / "_fourway_full" / "eval"
    sub = WORK / cid / "_fourway" / "eval"
    if (full / "lifton_devel.transcripts.tsv").exists() and \
       (full / "miniprot.transcripts.tsv").exists():
        return full, "full"
    if (sub / "lifton_devel.transcripts.tsv").exists() and \
       (sub / "miniprot.transcripts.tsv").exists():
        return sub, "subset"
    return None, None


def gff_paths(cid, mode):
    """Locate the devel output GFF + the miniprot GFF the eval ACTUALLY scored.

    For FULL mode the eval scores against the full miniprot run inside LiftOn
    (`_fourway_full/stable/lifton_output/miniprot/miniprot.gff3`, MP ids matching
    the full TSV) -- NOT the top-level subset `tools/miniprot/miniprot.gff3`.
    We pick the first existing candidate per mode."""
    if mode == "full":
        devel = WORK / cid / "_fourway_full" / "devel" / "devel.gff3"
        mp_candidates = [
            WORK / cid / "_fourway_full" / "stable" / "lifton_output" / "miniprot" / "miniprot.gff3",
            WORK / cid / "_fourway_full" / "devel" / "lifton_output" / "miniprot" / "miniprot.gff3",
            WORK / cid / "_fourway_full" / "tools" / "miniprot" / "miniprot.gff3",
            WORK / cid / "tools" / "miniprot" / "miniprot.gff3",
        ]
    else:
        devel = WORK / cid / "_fourway" / "lifton_devel" / "devel.gff3"
        mp_candidates = [
            WORK / cid / "_fourway" / "tools" / "miniprot" / "miniprot.gff3",
            WORK / cid / "tools" / "miniprot" / "miniprot.gff3",
            WORK / cid / "_fourway" / "stable" / "lifton_output" / "miniprot" / "miniprot.gff3",
        ]
    mp = next((p for p in mp_candidates if p.exists()), mp_candidates[-1])
    return devel, mp


def fpi(s):
    try:
        return float(s)
    except (TypeError, ValueError):
        return None


def load_tsv_byid(tsv):
    """Aggregate per ref-mRNA id: status (recovered>map_failed>lost), best coding PI,
    ref_prot_len, lifted_prot_len, lifted_dna_len, seqid, and the tool_feature_id
    (MP id for miniprot) of the best recovered row."""
    agg = {}
    with open(tsv) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            i = r["ref_mrna_id"]
            cur = agg.get(i)
            st = r["status"]
            pi = fpi(r["protein_identity"]) if r["is_coding"] == "1" and r["recovered"] == "1" else None
            if cur is None:
                agg[i] = {"status": st, "is_coding": r["is_coding"] == "1",
                          "pi": pi, "ref_prot_len": r["ref_prot_len"],
                          "lifted_prot_len": r["lifted_prot_len"],
                          "lifted_dna_len": r["lifted_dna_len"], "seqid": r["seqid"],
                          "feat_id": r["tool_feature_id"] if r["recovered"] == "1" else None}
            else:
                rank = {"recovered": 2, "map_failed": 1, "lost": 0}
                if rank.get(st, 0) > rank.get(cur["status"], 0):
                    cur["status"] = st
                    cur["seqid"] = r["seqid"]
                if pi is not None and (cur["pi"] is None or pi > cur["pi"]):
                    cur["pi"] = pi
                    cur["lifted_prot_len"] = r["lifted_prot_len"]
                    cur["feat_id"] = r["tool_feature_id"]
                elif cur.get("feat_id") is None and r["recovered"] == "1":
                    cur["feat_id"] = r["tool_feature_id"]
    return agg


def parse_attrs(field):
    d = {}
    for kv in field.strip().split(";"):
        if "=" in kv:
            k, v = kv.split("=", 1)
            d[k] = v
    return d


def load_miniprot_mpid_coords(mp_gff):
    """MP-id (ID= attr) -> (seqid,start,end). Robust join key: the eval TSV stores
    each recovered miniprot row's tool_feature_id = this MP id."""
    coords = {}
    if not mp_gff.exists():
        return coords
    with open(mp_gff) as fh:
        for ln in fh:
            if ln.startswith("#"):
                continue
            c = ln.rstrip("\n").split("\t")
            if len(c) < 9 or c[2] != "mRNA":
                continue
            a = parse_attrs(c[8])
            mpid = a.get("ID")
            if mpid:
                coords[mpid] = (c[0], int(c[3]), int(c[4]))
    return coords


def load_devel_mrna_intervals(devel_gff):
    """seqid -> sorted list of (start,end,pi) for devel mRNA features."""
    by = {}
    if not devel_gff.exists():
        return by
    with open(devel_gff) as fh:
        for ln in fh:
            if ln.startswith("#"):
                continue
            c = ln.rstrip("\n").split("\t")
            if len(c) < 9 or c[2] != "mRNA":
                continue
            a = parse_attrs(c[8])
            pi = fpi(a.get("protein_identity"))
            by.setdefault(c[0], []).append((int(c[3]), int(c[4]), pi))
    for k in by:
        by[k].sort()
    return by


def overlap_bp(a0, a1, b0, b1):
    return max(0, min(a1, b1) - max(a0, b0) + 1)


def best_devel_overlap(seqid, s, e, devel_by):
    """Return (max_overlap_ratio, overlapping_devel_pi) vs devel mRNAs on seqid."""
    ivs = devel_by.get(seqid, [])
    best_ratio, best_pi = 0.0, None
    qlen = e - s + 1
    for (ds, de, dpi) in ivs:
        if ds > e:
            break
        ov = overlap_bp(s, e, ds, de)
        if ov <= 0:
            continue
        ratio = ov / min(qlen, de - ds + 1)
        if ratio > best_ratio:
            best_ratio, best_pi = ratio, dpi
    return best_ratio, best_pi


def analyze_cell(cid, div):
    ed, mode = eval_dir(cid)
    if ed is None:
        return None
    mp = load_tsv_byid(ed / "miniprot.transcripts.tsv")
    dv = load_tsv_byid(ed / "lifton_devel.transcripts.tsv")
    devel_gff, mp_gff = gff_paths(cid, mode)
    mp_coords = load_miniprot_mpid_coords(mp_gff)
    devel_iv = load_devel_mrna_intervals(devel_gff)

    n_ref_coding = sum(1 for v in dv.values() if v["is_coding"])
    dv_rec = sum(1 for v in dv.values() if v["status"] == "recovered" and v["is_coding"])

    # ---- Idea A: miniprot-only rescue headroom ----
    gap = [i for i, v in mp.items()
           if v["status"] == "recovered" and v["is_coding"]
           and dv.get(i, {}).get("status") in ("lost", "map_failed", None)
           and dv.get(i, {}).get("is_coding", True)]
    A = {"n_gap": len(gap), "have_coords": 0, "no_overlap_genuine_new": 0,
         "overlap_dup": 0, "floor": {}, "examples_new": [], "examples_dup": [],
         "mp_gff": mp_gff.name, "coord_coverage": None}
    floors = [0.3, 0.5, 0.7]
    floor_new = {f: 0 for f in floors}
    for i in gap:
        mpid = mp[i].get("feat_id")
        sp = mp_coords.get(mpid) if mpid else None
        mppi = mp[i]["pi"]
        if sp is None or mppi is None:
            continue
        A["have_coords"] += 1
        seqid, s, e = sp
        ratio, dpi = best_devel_overlap(seqid, s, e, devel_iv)
        genuine = ratio <= 0.10   # engine-faithful: <=10% overlap would NOT be suppressed
        if genuine:
            A["no_overlap_genuine_new"] += 1
            for f in floors:
                if mppi >= f:
                    floor_new[f] += 1
            if len(A["examples_new"]) < 6 and mppi >= 0.7:
                A["examples_new"].append({"id": i, "miniprot_pi": round(mppi, 3),
                                          "ref_aa": mp[i]["ref_prot_len"], "seqid": seqid})
        else:
            A["overlap_dup"] += 1
            if len(A["examples_dup"]) < 4:
                A["examples_dup"].append({"id": i, "miniprot_pi": round(mppi, 3),
                                          "overlap_ratio": round(ratio, 2),
                                          "devel_overlap_pi": (round(dpi, 3) if dpi is not None else None)})
    for f in floors:
        A["floor"][str(f)] = {"genuine_new": floor_new[f],
                              "compl_gain": round(floor_new[f] / n_ref_coding, 5) if n_ref_coding else None}
    A["coord_coverage"] = round(A["have_coords"] / A["n_gap"], 3) if A["n_gap"] else 1.0

    # ---- Idea B: 3-way merge (accuracy) headroom ----
    B = {"thr": {}, "egregious": [], "proj_meanpi_gain": None, "n_devel_rec_coding": dv_rec}
    pairs = []  # (devel_pi, miniprot_pi) over devel recovered_coding ids
    for i, v in dv.items():
        if v["status"] != "recovered" or not v["is_coding"] or v["pi"] is None:
            continue
        mpi = mp.get(i, {}).get("pi") if mp.get(i, {}).get("status") == "recovered" else None
        pairs.append((i, v["pi"], mpi))
    for thr in [0.05, 0.10, 0.15]:
        B["thr"][str(thr)] = sum(1 for _, d, m in pairs if m is not None and (m - d) >= thr)
    if pairs:
        cur_mean = statistics.fmean([d for _, d, _ in pairs])
        new_mean = statistics.fmean([max(d, m) if m is not None else d for _, d, m in pairs])
        B["proj_meanpi_gain"] = round(new_mean - cur_mean, 5)
        B["cur_meanpi"] = round(cur_mean, 5)
    egr = sorted([(i, d, m) for i, d, m in pairs if m is not None and d < 0.10 and m > 0.5],
                 key=lambda x: x[2] - x[1], reverse=True)
    B["n_egregious_stub"] = len(egr)
    for i, d, m in egr[:5]:
        B["egregious"].append({"id": i, "devel_pi": round(d, 3), "miniprot_pi": round(m, 3)})

    # ---- Idea C: map_failed recovery ----
    mf = [i for i, v in dv.items() if v["status"] == "map_failed"]
    C = {"n_map_failed": len(mf),
         "miniprot_recoverable_pi>=0.5": sum(1 for i in mf
              if mp.get(i, {}).get("status") == "recovered" and (mp[i]["pi"] or 0) >= 0.5),
         "examples": []}
    for i in mf[:5]:
        C["examples"].append({"id": i, "ref_aa": dv[i]["ref_prot_len"],
                              "lifted_dna_len": dv[i]["lifted_dna_len"],
                              "miniprot_pi": (round(mp[i]["pi"], 3) if mp.get(i, {}).get("pi") else None)})

    # trust flag: known stale-eval full cells (see memory distant-fullgenome-bugs);
    # low coord coverage only matters when the gap is non-trivial.
    EVAL_SUSPECT = {"human_to_zebrafish"}  # full eval under-counts devel (stale gffutils DB)
    trust = "ok"
    if cid in EVAL_SUSPECT and mode == "full":
        trust = "eval_suspect"
    elif A["coord_coverage"] is not None and A["coord_coverage"] < 0.9:
        trust = "low_coverage" if A["n_gap"] >= 200 else "low_coverage_immaterial"

    return {"cell": cid, "divergence_class": div, "mode": mode, "trust": trust,
            "n_ref_coding": n_ref_coding, "devel_recovered_coding": dv_rec,
            "devel_compl": round(dv_rec / n_ref_coding, 5) if n_ref_coding else None,
            "A_rescue": A, "B_merge": B, "C_map_failed": C}


def main():
    divmap = load_div()
    cells = []
    for d in sorted(glob.glob(str(WORK / "*"))):
        cid = os.path.basename(d)
        ed, _ = eval_dir(cid)
        if ed is None:
            continue
        div = divmap.get(cid) or FALLBACK_DIV.get(cid, "?")
        cells.append((cid, div))
    results = []
    for cid, div in cells:
        try:
            r = analyze_cell(cid, div)
            if r:
                results.append(r)
                a = r["A_rescue"]
                cov = a["coord_coverage"]
                flag = "" if (cov is None or cov >= 0.9) else f" COVERAGE_LOW={cov}"
                print(f"[{div:26s}] {cid:34s} mode={r['mode']:6s} "
                      f"gap={a['n_gap']:6d} new={a['no_overlap_genuine_new']:6d} "
                      f"dup={a['overlap_dup']:5d} new@0.5={a['floor']['0.5']['genuine_new']:5d} "
                      f"(gain={a['floor']['0.5']['compl_gain']}) "
                      f"B_gain={r['B_merge'].get('proj_meanpi_gain')} "
                      f"mapfail={r['C_map_failed']['n_map_failed']}{flag}")
        except Exception as e:
            print(f"  !! {cid}: {type(e).__name__}: {e}")
    results.sort(key=lambda r: (DIV_ORDER.index(r["divergence_class"])
                                if r["divergence_class"] in DIV_ORDER else 99, r["cell"]))
    out = {"results": results}
    (HERE / "improvement_headroom.json").write_text(json.dumps(out, indent=2))
    print(f"\nwrote {HERE/'improvement_headroom.json'}  ({len(results)} cells)")
    write_md(results)
    return out


def write_md(results):
    DL = {"same_species": "Same species", "close_cross_species": "Close cross-species",
          "distant_cross_species": "Distant cross-species",
          "very_distant_cross_species": "Very-distant cross-species"}
    L = []
    L.append("# LiftOn v2.0.0 improvement headroom — brainstorm/diagnostic\n")
    L.append("Read-only headroom for four proposed engine improvements, measured from the "
             "existing 4-way benchmark artifacts (per-transcript eval TSVs + miniprot/devel "
             "output GFFs). **No engine code is changed.** Numbers are *upper bounds* — an A/B "
             "would still be required before any ship (24-cell byte-identity contract). "
             "Generated by `improvement_headroom.py`.\n")
    L.append("- **A** = miniprot-only rescue (completeness): `genuine_new@0.5` = ref genes devel "
             "loses that miniprot recovers at PI≥0.5 **and whose locus does NOT overlap (>10%) any "
             "devel feature** (0-redundancy rescues, by real interval test). `dup` = the overlapping "
             "fraction (would need suppression/replacement, not a clean add).")
    L.append("- **B** = 3-way merge (accuracy): `B_gain` = projected mean-PI lift over devel's "
             "recovered-coding set if a standalone miniprot model could win per transcript (oracle "
             "max(devel,miniprot)).")
    L.append("- **C** = map_failed recovery: `mapfail` transcripts (region lifted, protein didn't "
             "map); `recoverable` = how many miniprot recovers at PI≥0.5.\n")

    by_div = {}
    for r in results:
        by_div.setdefault(r["divergence_class"], []).append(r)
    L.append("## Per-cell headroom (grouped by divergence)\n")
    L.append("| Divergence | Cell | mode | devel compl | A: genuine_new@0.5 (compl gain) | A: dup | B_gain | C: mapfail (recov.) | trust |")
    L.append("|---|---|---|---|---|---|---|---|---|")
    for dc in DIV_ORDER + [k for k in by_div if k not in DIV_ORDER]:
        for r in sorted(by_div.get(dc, []), key=lambda r: r["cell"]):
            a = r["A_rescue"]; b = r["B_merge"]; c = r["C_map_failed"]
            f5 = a["floor"]["0.5"]
            L.append(f"| {DL.get(dc, dc)} | {r['cell']} | {r['mode']} | {r['devel_compl']} | "
                     f"{f5['genuine_new']} ({f5['compl_gain']}) | {a['overlap_dup']} | "
                     f"{b.get('proj_meanpi_gain')} | {c['n_map_failed']} "
                     f"({c['miniprot_recoverable_pi>=0.5']}) | {r['trust']} |")
    # tier aggregate over TRUSTWORTHY cells
    L.append("\n## Tier aggregate (trustworthy cells only; eval-suspect excluded)\n")
    L.append("| Divergence | cells | Σ genuine_new@0.5 | median B_gain | Σ mapfail |")
    L.append("|---|---|---|---|---|")
    for dc in DIV_ORDER:
        rs = [r for r in by_div.get(dc, []) if r["trust"] in ("ok", "low_coverage_immaterial")]
        if not rs:
            continue
        sn = sum(r["A_rescue"]["floor"]["0.5"]["genuine_new"] for r in rs)
        bgs = sorted(r["B_merge"].get("proj_meanpi_gain") or 0 for r in rs)
        medb = bgs[len(bgs) // 2]
        smf = sum(r["C_map_failed"]["n_map_failed"] for r in rs)
        L.append(f"| {DL.get(dc, dc)} | {len(rs)} | {sn} | {round(medb,5)} | {smf} |")
    L.append("\n_Note: `human_to_zebrafish` full is flagged `eval_suspect` (stale gffutils eval DB "
             "under-counts devel; see memory `lifton-distant-fullgenome-bugs`) and is excluded from "
             "aggregates. Same-species full cells with tiny gaps and low coord-coverage are "
             "`low_coverage_immaterial` (headroom ~0 regardless)._\n")
    (HERE / "improvement_headroom.md").write_text("\n".join(L) + "\n")
    print(f"wrote {HERE/'improvement_headroom.md'}")


if __name__ == "__main__":
    main()
