"""GFF3 serialization for the Lifton object model.

Iteration 19 extracted the four-level ``write_entry`` chain (the GFF3-serialisation
concern) out of the ``Lifton_TRANS`` god module (``lifton/lifton_class.py``) into
these free functions — the first concern-split of that 896-LOC class. Each function
takes a ``Lifton_*`` object and emits its GFF3 line(s) through the canonical
:func:`lifton.io.gff3_writer.format_feature`.

The functions are **duck-typed** — they read the object's ``.entry`` / ``.exons`` /
``.transcripts`` / ``.features`` attributes and never import ``lifton_class``, so the
dependency is strictly one-way (``lifton_class`` imports this module, never the
reverse — no import cycle). The bodies are a verbatim move of the prior methods, so
the emitted bytes are unchanged. The class ``write_entry`` methods remain as thin
delegations to these functions, preserving the recursive dispatch
(gene → trans → exon → CDS) exactly.
"""
from lifton import logger
from lifton.exceptions import LiftOnError
from lifton.io import gff3_writer


def write_gene(gene, fw, transcripts_stats_dict):
    if not gene.tmp:
        try:
            fw.write(gff3_writer.format_feature(gene.entry) + "\n")
        except Exception as e:
            logger.log_error(f"Failed to write GENE entry {gene.entry.id}: {e}")

    for key, trans in gene.transcripts.items():
        trans.write_entry(fw)
        TYPE = ""
        if gene.is_protein_coding and trans.entry.featuretype == "mRNA":
            TYPE = "coding"
        elif gene.is_non_coding and (trans.entry.featuretype == "ncRNA" or trans.entry.featuretype == "nc_RNA" or trans.entry.featuretype == "lncRNA" or trans.entry.featuretype == "lnc_RNA"):
            TYPE = "non-coding"
        else:
            TYPE = "other"
        if trans.ref_tran_id is None:
            continue
        if not trans.ref_tran_id in transcripts_stats_dict[TYPE].keys():
            transcripts_stats_dict[TYPE][trans.ref_tran_id] = 1
        else:
            transcripts_stats_dict[TYPE][trans.ref_tran_id] += 1


def write_feature(feature, fw):
    try:
        fw.write(gff3_writer.format_feature(feature.entry) + "\n")
    except Exception as e:
        logger.log_error(f"Failed to write FEATURE entry {feature.entry.id}: {e}")

    for key, sub_feature in feature.features.items():
        sub_feature.write_entry(fw)


def write_trans(trans, fw):
    # V1.5 + V5.4-V5.6 + V5.9: route through canonical writer that
    # percent-encodes reserved chars, sorts attributes (ID, Parent,
    # alphabetical), and validates start <= end. Skip child writes
    # if the parent fails so we don't ship orphan-Parent rows.
    #
    # Iter-21: `format_feature` raises `LiftOnInputError` (a `LiftOnError`,
    # NOT a ValueError) when a transcript has invalid coords — e.g. an
    # inverted `start > end` mRNA produced by `update_cds_list` on a
    # pathological no-valid-ORF model (full dog->cat crash). Catch the
    # `LiftOnError` base too so one malformed transcript is skip-and-logged
    # (matching the bare-`Exception` breadth of the sibling write_* funcs)
    # instead of propagating out of the parent-thread consume() write phase
    # and aborting the entire genome lift.
    try:
        fw.write(gff3_writer.format_feature(trans.entry) + "\n")
    except (OSError, ValueError, TypeError, AttributeError, LiftOnError) as e:
        logger.log_error(
            f"Failed to write TRANSCRIPT entry {trans.entry.id}: {e}"
        )
        return

    # Write out the exons first
    for exon in trans.exons:
        exon.write_entry(fw)
    # Write out the CDSs second
    for exon in trans.exons:
        if exon.cds is not None:
            exon.cds.write_entry(fw)


def write_exon(exon, fw):
    # Route through canonical writer (V5.4-V5.6, V5.9).
    try:
        fw.write(gff3_writer.format_feature(exon.entry) + "\n")
    except Exception as e:
        logger.log_error(f"Failed to write EXON entry {exon.entry.id}: {e}")


def write_cds(cds, fw):
    # Route through canonical writer (V5.4-V5.6, V5.9).
    try:
        fw.write(gff3_writer.format_feature(cds.entry) + "\n")
    except Exception as e:
        logger.log_error(f"Failed to write CDS entry {cds.entry.id}: {e}")
