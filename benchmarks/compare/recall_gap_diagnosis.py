"""Why does LiftOn miss coding genes that miniprot finds?

Read-only diagnostic behind the v1.0.12 improvement analysis
(``notes/lifton_v1.0.12_improvement_analysis.md``). For one transfer it takes
the shipped LiftOn output, miniprot's own output, and the neutral evaluator's
per-transcript TSVs, and reports:

1. **Recall** per tool at transcript, gene, primary-assembly gene and
   GeneID-collapsed level (``gene_level``).
2. **Rescue-gate replay.** For every coding gene LiftOn does not recover but
   miniprot recovers at protein identity >= ``--min-pi``, take miniprot's best
   model for that gene and report the first miniprot-only rescue gate it fails,
   using the rescue's own definitions:

   * (4) overlap: more than ``--overlap`` of the model's span overlaps a gene in
     the LiftOn output (a proxy for the final suppression tree);
   * (5) processed-pseudogene filter: one miniprot CDS where the reference
     transcript has more than one CDS segment (``ref_trans_exon_num_dict``
     counts CDS rows, ``lifton_utils.get_ref_liffover_features``);
   * (6) span band: the model's genomic span divided by the reference gene's
     CDS span, i.e. first-CDS start to the end of the last-starting CDS across
     all isoforms (``lifton_utils.py`` ``ref_features_len_dict``), outside the
     open interval ``--band``.

   Candidates rejected at (6) also report the model's protein coverage from
   miniprot's ``Target=<id> <start> <end>`` attribute.
3. **Rescue hit rank.** For LiftOn models tagged ``lifton_rescue=miniprot_only``,
   the miniprot ``Rank`` of the hit each was built from.
4. **Isoform opportunity.** For each rescued gene, the other reference
   transcripts of that gene whose miniprot hits lie on the same sequence and
   strand and overlap the rescued model.
5. **Haplotype displacement.** GeneID groups with a primary-assembly member
   and an alt/fix member where LiftOn recovers only the alt/fix copy.
"""
from __future__ import annotations

import argparse
import bisect
import json
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

from benchmarks.compare import gene_level

_COPY_SUFFIX = re.compile(r"_\d+$")


def _quantiles(values):
    if not values:
        return None
    ordered = sorted(values)

    def q(p):
        return round(ordered[int(p * (len(ordered) - 1))], 3)

    return {"n": len(ordered), "p10": q(0.10), "p25": q(0.25),
            "median": q(0.50), "p75": q(0.75), "p90": q(0.90)}


def load_reference(ref_gff, alt_seqids=frozenset()):
    """One streaming pass over the reference GFF3."""
    top_level = {}
    transcripts = {}
    cds_by_transcript = defaultdict(list)
    with open(ref_gff) as handle:
        for line in handle:
            if not line or line[0] == "#":
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) != 9:
                continue
            if c[2] == "CDS":
                attributes = gene_level.parse_attributes(c[8])
                for parent in attributes.get("Parent", ()):
                    cds_by_transcript[parent].append((int(c[3]), int(c[4])))
                continue
            if c[2] in gene_level._SUBFEATURE_TYPES:
                continue
            attributes = gene_level.parse_attributes(c[8])
            ids = attributes.get("ID")
            if not ids:
                continue
            if "Parent" in attributes:
                transcripts[ids[0]] = (c[0], attributes)
            else:
                top_level[ids[0]] = attributes
    tx_genes = {}
    for tx, (seqid, attributes) in transcripts.items():
        parent = (attributes.get("Parent") or [None])[0]
        tx_genes[tx] = gene_level.transcript_gene_record(
            tx, seqid, attributes, top_level.get(parent), alt_seqids)
    gene_cds = defaultdict(list)
    for tx, rows in cds_by_transcript.items():
        record = tx_genes.get(tx)
        if record is not None:
            gene_cds[record.gene_id].extend(rows)
    gene_span = {}
    for gene, rows in gene_cds.items():
        rows.sort()
        gene_span[gene] = rows[-1][1] - rows[0][0] + 1
    tx_cds_count = {tx: len(rows) for tx, rows in cds_by_transcript.items()}
    return tx_genes, gene_span, tx_cds_count


def load_miniprot(mp_gff):
    mrna, cds_count = {}, Counter()
    with open(mp_gff) as handle:
        for line in handle:
            if not line or line[0] == "#":
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) != 9:
                continue
            if c[2] == "CDS":
                attributes = gene_level.parse_attributes(c[8])
                cds_count[(attributes.get("Parent") or [""])[0]] += 1
            elif c[2] == "mRNA":
                attributes = gene_level.parse_attributes(c[8])
                target = (attributes.get("Target") or [""])[0].split()
                mrna[attributes["ID"][0]] = {
                    "seqid": c[0], "start": int(c[3]), "end": int(c[4]),
                    "strand": c[6],
                    "rank": int((attributes.get("Rank") or ["0"])[0]),
                    "identity": float((attributes.get("Identity") or ["nan"])[0]),
                    "target": target[0] if target else None,
                    "qstart": int(target[1]) if len(target) > 2 else None,
                    "qend": int(target[2]) if len(target) > 2 else None,
                }
    return mrna, cds_count


def load_lifton(lo_gff):
    genes = defaultdict(list)
    rescued = []
    with open(lo_gff) as handle:
        for line in handle:
            if not line or line[0] == "#":
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) != 9:
                continue
            if c[2] == "gene":
                genes[c[0]].append((int(c[3]), int(c[4])))
            elif c[2] == "mRNA" and "lifton_rescue=miniprot_only" in c[8]:
                attributes = gene_level.parse_attributes(c[8])
                rescued.append({
                    "id": attributes["ID"][0],
                    "gene": (attributes.get("Parent") or [""])[0],
                    "seqid": c[0], "start": int(c[3]), "end": int(c[4]),
                    "strand": c[6],
                })
    for seqid in genes:
        genes[seqid].sort()
    starts = {seqid: [s for s, _ in rows] for seqid, rows in genes.items()}
    return genes, starts, rescued


def overlap_fraction(genes, starts, seqid, start, end):
    """Largest fraction of [start, end] covered by one gene interval."""
    rows = genes.get(seqid, [])
    index = bisect.bisect_right(starts.get(seqid, []), end)
    best = 0.0
    # Gene intervals are sorted by start; scan back far enough to reach long
    # genes that start well before this interval.
    for j in range(index - 1, max(-1, index - 200), -1):
        s, e = rows[j]
        overlap = min(e, end) - max(s, start) + 1
        if overlap > 0:
            best = max(best, overlap / (end - start + 1))
    return best


def _pi(row):
    try:
        return float(row["protein_identity"] or 0)
    except (TypeError, ValueError):
        return 0.0


def replay_gates(candidates, mp_mrna, mp_cds, gene_span, tx_cds_count,
                 lo_genes, lo_starts, overlap, band, ref_prot_len):
    first_fail = Counter()
    ratios, coverage = [], []
    lo, hi = band
    for gene, (tx, row) in candidates.items():
        model = mp_mrna.get(row["tool_feature_id"])
        if model is None or not gene_span.get(gene):
            first_fail["unresolvable"] += 1
            continue
        span = model["end"] - model["start"] + 1
        ratio = span / gene_span[gene]
        ratios.append(ratio)
        if overlap_fraction(lo_genes, lo_starts, model["seqid"],
                            model["start"], model["end"]) > overlap:
            first_fail["(4) overlaps a lifted gene"] += 1
        elif mp_cds[row["tool_feature_id"]] == 1 and tx_cds_count.get(tx, 0) > 1:
            first_fail["(5) single-CDS processed-pseudogene filter"] += 1
        elif not (lo < ratio < hi):
            first_fail["(6) span ratio outside the rescue band"] += 1
            length = ref_prot_len.get(tx)
            if length and model["qstart"] is not None:
                coverage.append(
                    min(1.0, (model["qend"] - model["qstart"] + 1) / length))
        else:
            first_fail["passes gates 4-6 (fails later: PI floor/ORF)"] += 1
    n = sum(first_fail.values())
    return {
        "n_candidates": n,
        "first_fail": {k: v for k, v in first_fail.most_common()},
        "first_fail_fraction": {
            k: round(v / n, 4) for k, v in first_fail.most_common()} if n else {},
        "span_ratio": _quantiles(ratios),
        "span_ratio_below_lo": round(sum(r <= lo for r in ratios) / len(ratios), 4)
        if ratios else None,
        "span_ratio_above_hi": round(sum(r >= hi for r in ratios) / len(ratios), 4)
        if ratios else None,
        "span_rejected_protein_coverage": _quantiles(coverage),
        "span_rejected_coverage_ge_0.8": round(
            sum(c >= 0.8 for c in coverage) / len(coverage), 4) if coverage else None,
    }


def diagnose(args):
    alt = gene_level.load_alt_seqids(args.alt_seqids)
    tx_genes, gene_span, tx_cds_count = load_reference(args.ref_gff, alt)
    mp_mrna, mp_cds = load_miniprot(args.miniprot_gff)
    lo_genes, lo_starts, rescued = load_lifton(args.lifton_gff)

    tables = {"lifton": args.lifton_tsv, "miniprot": args.miniprot_tsv}
    if args.liftoff_tsv:
        tables["liftoff"] = args.liftoff_tsv
    rows = {name: gene_level.annotate_rows(gene_level.read_transcript_tsv(path),
                                           tx_genes)
            for name, path in tables.items()}
    recall = {name: gene_level.gene_level_summary(r) for name, r in rows.items()}

    def coding(r):
        return str(r["is_coding"]).strip().lower() in ("1", "true")

    def recovered(r):
        return str(r["recovered"]).strip().lower() in ("1", "true")

    lo_rows = {r["ref_mrna_id"]: r for r in rows["lifton"] if coding(r)}
    mp_rows = {r["ref_mrna_id"]: r for r in rows["miniprot"] if coding(r)}
    ref_prot_len = {tx: int(r["ref_prot_len"] or 0) for tx, r in mp_rows.items()}
    genes_recovered = {r["ref_gene_id"] for r in lo_rows.values() if recovered(r)}

    candidates = {}
    for tx, row in mp_rows.items():
        if not recovered(row) or _pi(row) < args.min_pi:
            continue
        gene = row["ref_gene_id"]
        if gene in genes_recovered:
            continue
        if gene not in candidates or _pi(row) > _pi(candidates[gene][1]):
            candidates[gene] = (tx, row)
    primary = {g: v for g, v in candidates.items()
               if v[1]["ref_seqid_class"] == gene_level.PRIMARY}
    replay_args = (mp_mrna, mp_cds, gene_span, tx_cds_count, lo_genes, lo_starts,
                   args.overlap, tuple(args.band), ref_prot_len)
    n_primary_genes = recall["lifton"]["primary"]["n_coding_genes"]
    gates = {"all": replay_gates(candidates, *replay_args),
             "primary": replay_gates(primary, *replay_args)}
    for key in ("all", "primary"):
        gates[key]["min_pi"] = args.min_pi
    span_primary = gates["primary"]["first_fail"].get(
        "(6) span ratio outside the rescue band", 0)
    gates["primary"]["span_class_gene_recall_ceiling"] = (
        round(span_primary / n_primary_genes, 4) if n_primary_genes else None)

    # Rescue hit rank and isoform opportunity.
    hits_by_target = defaultdict(list)
    for mid, model in mp_mrna.items():
        if model["target"]:
            hits_by_target[model["target"]].append(model)
    transcripts_by_gene = defaultdict(list)
    for tx, record in tx_genes.items():
        if tx in mp_rows:
            transcripts_by_gene[record.gene_id].append(tx)
    ranks = Counter()
    rank_gt1_identity_gap = []
    isoform = Counter()
    for model in rescued:
        base = _COPY_SUFFIX.sub("", model["id"])
        hits = [h for h in hits_by_target.get(base, [])
                if h["seqid"] == model["seqid"] and h["start"] <= model["end"]
                and model["start"] <= h["end"]]
        if hits:
            used = max(hits, key=lambda h: min(h["end"], model["end"])
                       - max(h["start"], model["start"]))
            ranks[used["rank"]] += 1
            if used["rank"] > 1:
                best = max(hits_by_target[base], key=lambda h: h["identity"])
                rank_gt1_identity_gap.append(best["identity"] - used["identity"])
        record = tx_genes.get(base)
        gene = record.gene_id if record else model["gene"]
        for other in transcripts_by_gene.get(gene, []):
            if other == base:
                continue
            isoform["other_transcripts"] += 1
            colocated = any(
                h["seqid"] == model["seqid"] and h["strand"] == model["strand"]
                and h["start"] <= model["end"] and model["start"] <= h["end"]
                for h in hits_by_target.get(other, []))
            if not colocated:
                continue
            pi = _pi(mp_rows[other])
            isoform["colocated_pi_ge_0.30"] += pi >= 0.30
            isoform["colocated_pi_ge_0.50"] += pi >= 0.50
    n_rescued = len(rescued)
    rescue = {
        "n_rescued_transcripts": n_rescued,
        "hit_rank": dict(sorted(ranks.items())),
        "rank_gt1": sum(v for k, v in ranks.items() if k > 1),
        "rank_gt1_identity_gap": _quantiles(rank_gt1_identity_gap),
    }
    isoform_block = dict(isoform)
    isoform_block["per_rescued_gene_pi_ge_0.30"] = (
        round(isoform["colocated_pi_ge_0.30"] / n_rescued, 2) if n_rescued else None)

    # Haplotype displacement among GeneID groups.
    members = defaultdict(lambda: {"primary": set(), "other": set()})
    for record in tx_genes.values():
        side = "primary" if record.seqid_class == gene_level.PRIMARY else "other"
        members[record.gene_key][side].add(record.gene_id)
    both = displaced = 0
    for group in members.values():
        if not group["primary"] or not group["other"]:
            continue
        both += 1
        if (group["other"] & genes_recovered) and not (group["primary"] & genes_recovered):
            displaced += 1
    haplotype = {"geneid_groups_with_primary_and_alt_or_fix": both,
                 "recovered_only_as_alt_or_fix_copy": displaced}

    return {
        "inputs": {k: str(v) for k, v in vars(args).items() if k != "json"},
        "recall": recall,
        "gate_replay": gates,
        "rescue_hit_rank": rescue,
        "isoform_opportunity": isoform_block,
        "haplotype_displacement": haplotype,
    }


def build_parser():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--ref-gff", required=True)
    parser.add_argument("--lifton-gff", required=True)
    parser.add_argument("--lifton-tsv", required=True)
    parser.add_argument("--miniprot-gff", required=True)
    parser.add_argument("--miniprot-tsv", required=True)
    parser.add_argument("--liftoff-tsv", default=None)
    parser.add_argument("--min-pi", type=float, default=0.5)
    parser.add_argument("--overlap", type=float, default=0.1)
    parser.add_argument("--band", type=float, nargs=2, default=(0.5, 2.0),
                        metavar=("LO", "HI"))
    parser.add_argument("--alt-seqids", default=None)
    parser.add_argument("--json", default=None)
    return parser


def main(argv=None):
    args = build_parser().parse_args(argv)
    result = diagnose(args)
    text = json.dumps(result, indent=2, sort_keys=True)
    if args.json:
        Path(args.json).write_text(text + "\n")
    rec = result["recall"]
    for name, summary in rec.items():
        print(f"{name:9} tx {summary['transcript_recall_coding']}  gene "
              f"{summary['gene_recall_coding']}  primary gene "
              f"{summary['primary']['gene_recall_coding']}  GeneID "
              f"{summary['geneid_collapsed']['recall']}")
    for key in ("all", "primary"):
        block = result["gate_replay"][key]
        print(f"[{key}] candidates {block['n_candidates']}: {block['first_fail']}")
    print("rescue:", result["rescue_hit_rank"])
    print("isoform:", result["isoform_opportunity"])
    print("haplotype:", result["haplotype_displacement"])
    return 0


if __name__ == "__main__":
    sys.exit(main())
