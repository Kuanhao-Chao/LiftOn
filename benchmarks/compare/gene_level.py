"""Gene-level, primary-assembly, and GeneID-collapsed recall.

The neutral evaluator (``evaluator.py``) scores one row per reference
transcript. Two properties of real references make transcript recall
misleading when it is read on its own:

* **Isoform multiplicity.** A human gene carries about six annotated
  transcripts, so a tool that recovers one model per gene moves transcript
  recall by a sixth of what it moves gene recall.
* **Duplicate haplotype copies.** GRCh38 RefSeq annotates alt-locus and
  fix-patch sequences, so 14.6 % of human coding genes are second copies of a
  primary-assembly gene. A target genome has one locus for both, so at most one
  copy can be placed and the other always reads as a miss.

This module attaches every reference transcript to its gene, to a gene key
that collapses haplotype copies (the NCBI ``GeneID`` cross-reference, which
alt and fix copies share with the primary gene), and to the class of the
sequence it sits on, then summarizes recall at those levels. It is
standard-library only, so it can re-score existing ``*.transcripts.tsv`` files
without rebuilding a feature database, and the evaluator uses it only when a
caller asks for it: the evaluator's own TSV and summary bytes are unchanged
by default, which keeps sealed evidence byte-reproducible.

Command line::

    python -m benchmarks.compare.gene_level --ref-gff REF.gff \\
        --tsv lifton=eval/lifton.transcripts.tsv \\
        --tsv miniprot=eval/miniprot.transcripts.tsv --json out.json
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
import urllib.parse
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Mapping

# Sub-feature rows never act as transcript containers (mirrors
# evaluator._SUBFEATURE_TYPES).
_SUBFEATURE_TYPES = frozenset({
    "exon", "CDS", "start_codon", "stop_codon",
    "five_prime_UTR", "three_prime_UTR", "intron",
})

PRIMARY = "primary"
ALT = "alt"
FIX = "fix"


@dataclass(frozen=True)
class TranscriptGene:
    """How one reference transcript is grouped for gene-level recall."""

    gene_id: str
    gene_key: str
    key_source: str
    seqid_class: str


def parse_attributes(text: str) -> dict[str, list[str]]:
    """Parse GFF3 column 9 into ``{key: [values]}`` with percent-decoding."""
    attributes: dict[str, list[str]] = {}
    for field in text.strip().split(";"):
        if "=" not in field:
            continue
        key, value = field.split("=", 1)
        attributes[key.strip()] = [
            urllib.parse.unquote(item) for item in value.split(",")
        ]
    return attributes


def geneid_from_dbxref(values: Iterable[str] | None) -> str | None:
    """Return ``GeneID:<n>`` from a Dbxref value list, or None."""
    for value in values or ():
        value = value.strip()
        if value.startswith("GeneID:") and len(value) > len("GeneID:"):
            return value
    return None


def seqid_class(seqid: str, alt_seqids: frozenset[str] = frozenset()) -> str:
    """Classify a reference sequence as ``primary``, ``alt`` or ``fix``.

    UCSC-style GRCh38 names mark alternate loci with ``_alt`` and patches with
    ``_fix``. Unlocalized (``_random``) and unplaced (``chrUn_``) scaffolds are
    part of the primary assembly and are not duplicates, so they stay
    ``primary``. Accession-style names cannot be classified from the name
    alone; pass those sequences explicitly through ``alt_seqids``.
    """
    if seqid in alt_seqids:
        return ALT
    if seqid.endswith("_alt"):
        return ALT
    if seqid.endswith("_fix"):
        return FIX
    return PRIMARY


def transcript_gene_record(
    transcript_id: str,
    seqid: str,
    attributes: Mapping[str, list[str]],
    parent_attributes: Mapping[str, list[str]] | None = None,
    alt_seqids: frozenset[str] = frozenset(),
) -> TranscriptGene:
    """Build the grouping record for one reference transcript.

    The gene key is the transcript's own ``GeneID`` cross-reference, then its
    parent's, and falls back to the parent ID (no collapsing) when neither is
    present, as in Ensembl-derived GFF3.
    """
    parents = attributes.get("Parent") or [transcript_id]
    gene_id = parents[0]
    key = geneid_from_dbxref(attributes.get("Dbxref"))
    source = "Dbxref:GeneID"
    if key is None and parent_attributes is not None:
        key = geneid_from_dbxref(parent_attributes.get("Dbxref"))
    if key is None:
        key, source = gene_id, "gene_id"
    return TranscriptGene(gene_id, key, source, seqid_class(seqid, alt_seqids))


def load_alt_seqids(path: str | Path | None) -> frozenset[str]:
    """Read one sequence name per line (``#`` comments allowed)."""
    if not path:
        return frozenset()
    names = set()
    with open(path) as handle:
        for line in handle:
            line = line.split("#", 1)[0].strip()
            if line:
                names.add(line)
    return frozenset(names)


def load_transcript_genes(
    ref_gff: str | Path, alt_seqids: frozenset[str] = frozenset(),
) -> dict[str, TranscriptGene]:
    """Stream a reference GFF3 into ``{transcript_id: TranscriptGene}``.

    Every non-sub-feature row with both ``ID`` and ``Parent`` is recorded, so
    the map covers every transcript-container type the evaluator scores
    (``mRNA``, ``transcript``, ncRNA types). Lookups are by the evaluator's
    ``ref_mrna_id``, so the extra entries are harmless.
    """
    top_level: dict[str, dict[str, list[str]]] = {}
    children: list[tuple[str, str, dict[str, list[str]]]] = []
    with open(ref_gff) as handle:
        for line in handle:
            if not line or line[0] == "#":
                continue
            columns = line.rstrip("\n").split("\t")
            if len(columns) != 9 or columns[2] in _SUBFEATURE_TYPES:
                continue
            attributes = parse_attributes(columns[8])
            feature_ids = attributes.get("ID")
            if not feature_ids:
                continue
            if "Parent" in attributes:
                children.append((feature_ids[0], columns[0], attributes))
            else:
                top_level[feature_ids[0]] = attributes
    records = {}
    for transcript_id, seqid, attributes in children:
        parent = (attributes.get("Parent") or [None])[0]
        records[transcript_id] = transcript_gene_record(
            transcript_id, seqid, attributes, top_level.get(parent), alt_seqids,
        )
    return records


def _flag(value) -> bool:
    return str(value).strip().lower() in ("1", "true")


def annotate_rows(
    rows: Iterable[dict], transcript_genes: Mapping[str, TranscriptGene],
) -> list[dict]:
    """Attach ``ref_gene_id``/``ref_gene_key``/``ref_seqid_class`` to rows.

    A row whose transcript is absent from the map keeps its own ID as its gene
    and is classed ``primary``, so it still counts rather than vanishing.
    """
    annotated = []
    for row in rows:
        transcript_id = row["ref_mrna_id"]
        record = transcript_genes.get(transcript_id)
        row = dict(row)
        if record is None:
            row["ref_gene_id"] = row["ref_gene_key"] = transcript_id
            row["ref_seqid_class"] = PRIMARY
        else:
            row["ref_gene_id"] = record.gene_id
            row["ref_gene_key"] = record.gene_key
            row["ref_seqid_class"] = record.seqid_class
        annotated.append(row)
    return annotated


def _ratio(numerator: int, denominator: int):
    return round(numerator / denominator, 5) if denominator else None


def _grouped_recall(rows: list[dict], key: str) -> tuple[int, int]:
    groups: dict[str, bool] = defaultdict(bool)
    for row in rows:
        groups[row[key]] |= _flag(row["recovered"])
    return len(groups), sum(groups.values())


def gene_level_summary(rows: Iterable[dict]) -> dict:
    """Summarize coding recall per transcript, gene, and GeneID group.

    ``rows`` must carry ``recovered``, ``is_coding``, ``ref_gene_id``,
    ``ref_gene_key`` and ``ref_seqid_class`` (see :func:`annotate_rows`). A
    coding gene is one with at least one coding transcript; it is recovered
    when any of its coding transcripts is.
    """
    coding = [row for row in rows if _flag(row["is_coding"])]
    primary = [row for row in coding if row["ref_seqid_class"] == PRIMARY]
    n_genes, n_genes_recovered = _grouped_recall(coding, "ref_gene_id")
    n_primary, n_primary_recovered = _grouped_recall(primary, "ref_gene_id")
    n_groups, n_groups_recovered = _grouped_recall(coding, "ref_gene_key")
    genes_by_class: dict[str, set] = defaultdict(set)
    for row in coding:
        genes_by_class[row["ref_seqid_class"]].add(row["ref_gene_id"])
    n_primary_tx_recovered = sum(_flag(row["recovered"]) for row in primary)
    return {
        "n_coding_transcripts": len(coding),
        "transcript_recall_coding": _ratio(
            sum(_flag(row["recovered"]) for row in coding), len(coding)),
        "n_coding_genes": n_genes,
        "n_coding_genes_recovered": n_genes_recovered,
        "gene_recall_coding": _ratio(n_genes_recovered, n_genes),
        "n_coding_genes_by_seqid_class": {
            name: len(genes) for name, genes in sorted(genes_by_class.items())
        },
        "primary": {
            "n_coding_transcripts": len(primary),
            "transcript_recall_coding": _ratio(
                n_primary_tx_recovered, len(primary)),
            "n_coding_genes": n_primary,
            "n_coding_genes_recovered": n_primary_recovered,
            "gene_recall_coding": _ratio(n_primary_recovered, n_primary),
        },
        "geneid_collapsed": {
            "n_groups": n_groups,
            "n_groups_recovered": n_groups_recovered,
            "recall": _ratio(n_groups_recovered, n_groups),
        },
    }


def read_transcript_tsv(path: str | Path) -> list[dict]:
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def summarize_tsv(path, transcript_genes) -> dict:
    return gene_level_summary(
        annotate_rows(read_transcript_tsv(path), transcript_genes))


def _parse_labelled(value: str) -> tuple[str, str]:
    if "=" not in value:
        raise argparse.ArgumentTypeError("expected LABEL=PATH")
    label, path = value.split("=", 1)
    if not label or not path:
        raise argparse.ArgumentTypeError("expected LABEL=PATH")
    return label, path


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--ref-gff", required=True)
    parser.add_argument("--tsv", action="append", type=_parse_labelled,
                        required=True, metavar="LABEL=PATH",
                        help="an evaluator *.transcripts.tsv; repeatable")
    parser.add_argument("--alt-seqids", default=None,
                        help="file listing extra alternate-locus sequences")
    parser.add_argument("--json", default=None, help="write the summaries here")
    return parser


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    genes = load_transcript_genes(args.ref_gff, load_alt_seqids(args.alt_seqids))
    summaries = {label: summarize_tsv(path, genes) for label, path in args.tsv}
    header = (f"{'tool':<14}{'tx recall':>10}{'gene recall':>12}"
              f"{'primary gene':>14}{'GeneID groups':>15}")
    print(header)
    for label, summary in summaries.items():
        print(f"{label:<14}{summary['transcript_recall_coding']!s:>10}"
              f"{summary['gene_recall_coding']!s:>12}"
              f"{summary['primary']['gene_recall_coding']!s:>14}"
              f"{summary['geneid_collapsed']['recall']!s:>15}")
    if args.json:
        Path(args.json).write_text(json.dumps(summaries, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    sys.exit(main())
