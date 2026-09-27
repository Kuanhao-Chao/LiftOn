"""Independent target-annotation gate for default-on second-locus rescue.

Score added LiftOn CDS models against assembly-matched RefSeq transcripts,
not against the reference annotation used to build those models. The shuffled
control preserves each model's sequence, strand, CDS structure, and lengths;
only its locus is changed. It excludes coding sequence containing N bases.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import random
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from urllib.parse import unquote

from intervaltree import IntervalTree
from pyfaidx import Fasta


@dataclass
class Transcript:
    identifier: str
    gene_id: str
    seqid: str
    strand: str
    cds: list[tuple[int, int]] = field(default_factory=list)
    second_locus: bool = False

    def coding_intervals(self):
        intervals = sorted(self.cds)
        merged = []
        for start, end in intervals:
            if merged and start <= merged[-1][1] + 1:
                merged[-1] = (merged[-1][0], max(merged[-1][1], end))
            else:
                merged.append((start, end))
        return merged


def _attrs(raw):
    return {key: unquote(value) for field in raw.split(";")
            if "=" in field for key, value in [field.split("=", 1)]}


def _read_models(path, *, target=False):
    genes = {}
    transcripts = {}
    regions = {}
    assembly = None
    with Path(path).open() as handle:
        for line in handle:
            if line.startswith("##sequence-region "):
                _, seqid, start, end = line.strip().split()
                if start != "1":
                    raise ValueError(f"Unexpected sequence-region start in {path}: {line.strip()}")
                regions[seqid] = int(end)
                continue
            if line.startswith("#!genome-build-accession "):
                assembly = line.strip().split(":", 1)[-1]
                continue
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                raise ValueError(f"Malformed feature row in {path}")
            seqid, _, kind, start, end, _, strand, _, raw_attrs = fields
            attrs = _attrs(raw_attrs)
            if kind == "gene":
                genes[attrs.get("ID", "")] = attrs.get("gene_biotype") == "protein_coding"
            elif kind in ("mRNA", "transcript"):
                identifier = attrs.get("ID")
                parent = attrs.get("Parent", "").split(",")[0]
                if identifier and identifier not in transcripts:
                    transcripts[identifier] = Transcript(
                        identifier, parent, seqid, strand,
                        second_locus=attrs.get("lifton_rescue_second_locus", "").lower() == "true",
                    )
            elif kind == "CDS":
                for parent in attrs.get("Parent", "").split(","):
                    transcript = transcripts.get(parent)
                    if transcript is not None and transcript.seqid == seqid:
                        transcript.cds.append((int(start), int(end)))
    if target:
        transcripts = {key: trans for key, trans in transcripts.items()
                       if genes.get(trans.gene_id) and trans.cds}
    else:
        transcripts = {key: trans for key, trans in transcripts.items() if trans.cds}
    return transcripts, regions, assembly


def _fasta_lengths(path):
    index = Path(str(path) + ".fai")
    if not index.exists():
        raise FileNotFoundError(f"FASTA index required for assembly check: {index}")
    lengths = {}
    with index.open() as handle:
        for line in handle:
            seqid, length, *_ = line.rstrip("\n").split("\t")
            lengths[seqid] = int(length)
    return lengths


def _sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _feature_rows(path):
    with Path(path).open() as handle:
        return Counter(line for line in handle if line and not line.startswith("#"))


def _overlap_bases(left, right):
    shared = 0
    i = j = 0
    while i < len(left) and j < len(right):
        a, b = left[i], right[j]
        shared += max(0, min(a[1], b[1]) - max(a[0], b[0]) + 1)
        if a[1] < b[1]:
            i += 1
        else:
            j += 1
    return shared


class TargetIndex:
    def __init__(self, transcripts):
        self.transcripts = transcripts
        self.trees = defaultdict(IntervalTree)
        self.intervals = {}
        self.lengths = {}
        for identifier, trans in transcripts.items():
            spans = trans.coding_intervals()
            self.intervals[identifier] = spans
            self.lengths[identifier] = sum(end - start + 1 for start, end in spans)
            self.trees[(trans.seqid, trans.strand)].addi(spans[0][0], spans[-1][1] + 1, identifier)

    def matches(self, model):
        spans = model.coding_intervals()
        if not spans:
            return {}
        model_length = sum(end - start + 1 for start, end in spans)
        matches = {}
        for hit in self.trees[(model.seqid, model.strand)].overlap(spans[0][0], spans[-1][1] + 1):
            target = self.transcripts[hit.data]
            shared = _overlap_bases(spans, self.intervals[hit.data])
            model_fraction = shared / model_length
            target_fraction = shared / self.lengths[hit.data]
            if model_fraction >= 0.5 and target_fraction >= 0.5:
                score = min(model_fraction, target_fraction)
                matches[target.gene_id] = max(matches.get(target.gene_id, 0), score)
        return matches


def _distinct_matches(models, target_index, off_covered):
    candidates = []
    details = {}
    for model in models:
        matches = target_index.matches(model)
        new_matches = {gene: score for gene, score in matches.items() if gene not in off_covered}
        details[model.identifier] = {
            "source_gene": model.gene_id,
            "seqid": model.seqid,
            "strand": model.strand,
            "cds": model.coding_intervals(),
            "matching_target_genes": sorted(matches),
            "new_target_genes": sorted(new_matches),
        }
        candidates.extend((score, model.identifier, gene) for gene, score in new_matches.items())
    used_models = set()
    used_genes = set()
    for score, model_id, gene in sorted(candidates, reverse=True):
        if model_id in used_models or gene in used_genes:
            continue
        used_models.add(model_id)
        used_genes.add(gene)
        details[model_id]["assigned_target_gene"] = gene
        details[model_id]["reciprocal_cds_score"] = round(score, 5)
    return len(used_models), details


def _shuffle(model, rng, lengths, fasta, max_attempts=200):
    spans = model.coding_intervals()
    left, right = spans[0][0], spans[-1][1]
    max_start = lengths[model.seqid] - (right - left)
    if max_start < 1:
        raise ValueError(f"Model {model.identifier} exceeds sequence {model.seqid}")
    for _ in range(max_attempts):
        shift = rng.randint(1, max_start) - left
        shifted = [(start + shift, end + shift) for start, end in spans]
        if all("N" not in str(fasta[model.seqid][start - 1:end]).upper()
               for start, end in shifted):
            return Transcript(model.identifier, model.gene_id, model.seqid,
                              model.strand, shifted)
    raise RuntimeError(f"Could not draw an N-free decoy for {model.identifier}")


def evaluate(off_path, on_path, truth_path, fasta_path, assembly, *, seed=1409, replicates=1000):
    off, _, _ = _read_models(off_path)
    on, _, _ = _read_models(on_path)
    target, regions, observed_assembly = _read_models(truth_path, target=True)
    if observed_assembly != assembly:
        raise ValueError(f"Target annotation is {observed_assembly}, expected {assembly}")
    lengths = _fasta_lengths(fasta_path)
    mismatches = {seqid: (length, lengths.get(seqid)) for seqid, length in regions.items()
                  if lengths.get(seqid) != length}
    if mismatches:
        raise ValueError(f"Target annotation and FASTA sequence lengths differ: {mismatches}")
    for trans in list(off.values()) + list(on.values()):
        if trans.seqid not in lengths:
            raise ValueError(f"Output sequence {trans.seqid} is absent from target FASTA")

    target_index = TargetIndex(target)
    off_rows = _feature_rows(off_path)
    on_rows = _feature_rows(on_path)
    lost_rows = sum((off_rows - on_rows).values())
    added_rows = sum((on_rows - off_rows).values())
    off_covered = set()
    for model in off.values():
        off_covered.update(target_index.matches(model))
    added = sorted((model for model in on.values() if model.second_locus),
                   key=lambda model: model.identifier)
    observed, details = _distinct_matches(added, target_index, off_covered)

    rng = random.Random(seed)
    null_counts = []
    if replicates:
        fasta = Fasta(str(fasta_path), as_raw=True, rebuild=False)
        try:
            for _ in range(replicates):
                shuffled = [_shuffle(model, rng, lengths, fasta) for model in added]
                count, _ = _distinct_matches(shuffled, target_index, off_covered)
                null_counts.append(count)
        finally:
            fasta.close()
    p_value = ((1 + sum(count >= observed for count in null_counts)) / (replicates + 1)
               if replicates else None)
    passes = (bool(added) and lost_rows == 0 and observed * 2 >= len(added)
              and p_value is not None and p_value < 0.01)
    return {
        "assembly": assembly,
        "inputs_sha256": {str(path): _sha256(path) for path in
                          (off_path, on_path, truth_path, fasta_path)},
        "target_annotated_coding_transcripts": len(target),
        "target_covered_by_off": len(off_covered),
        "added_models": len(added),
        "off_feature_rows_lost": lost_rows,
        "on_feature_rows_added": added_rows,
        "independently_supported_new_loci": observed,
        "support_fraction": observed / len(added) if added else 0,
        "null": {"seed": seed, "replicates": replicates,
                 "mean": sum(null_counts) / replicates if replicates else None,
                 "maximum": max(null_counts) if null_counts else None,
                 "p_value": p_value},
        "gate_pass": passes,
        "placements": details,
    }


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--off", type=Path, required=True)
    parser.add_argument("--on", type=Path, required=True)
    parser.add_argument("--truth", type=Path, required=True)
    parser.add_argument("--fasta", type=Path, required=True)
    parser.add_argument("--assembly", required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--seed", type=int, default=1409)
    parser.add_argument("--replicates", type=int, default=1000)
    args = parser.parse_args(argv)
    report = evaluate(args.off, args.on, args.truth, args.fasta, args.assembly,
                      seed=args.seed, replicates=args.replicates)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(f"supported {report['independently_supported_new_loci']}/{report['added_models']}; "
          f"null p={report['null']['p_value']}; gate {'PASS' if report['gate_pass'] else 'FAIL'}")
    return 0 if report["gate_pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
