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
import io
import os

from lifton import logger
from lifton.exceptions import LiftOnError
from lifton.io import gff3_writer


def _mrna_harmonize_enabled():
    """GH #28: canonicalize CODING transcripts to featuretype 'mRNA' (default ON).
    Set LIFTON_NO_MRNA_HARMONIZE=1 to preserve the reference featuretype instead."""
    return os.environ.get("LIFTON_NO_MRNA_HARMONIZE", "0").strip().lower() in (
        "", "0", "false", "no", "off")


def _render_exon(exon) -> str:
    return gff3_writer.format_feature(exon.entry) + "\n"


def _render_cds(cds) -> str:
    return gff3_writer.format_feature(cds.entry) + "\n"


def _render_trans(trans) -> str:
    """Render one transcript hierarchy without publishing partial rows."""
    if _mrna_harmonize_enabled() and any(
            getattr(e, "cds", None) is not None
            for e in getattr(trans, "exons", [])):
        trans.entry.featuretype = "mRNA"

    out = io.StringIO()
    out.write(gff3_writer.format_feature(trans.entry) + "\n")
    for exon in trans.exons:
        out.write(_render_exon(exon))
    for exon in trans.exons:
        if exon.cds is not None:
            out.write(_render_cds(exon.cds))
    return out.getvalue()


def _render_feature(feature) -> str:
    """Render a generic gene-like subtree as one atomic block."""
    out = io.StringIO()
    out.write(gff3_writer.format_feature(feature.entry) + "\n")
    for sub_feature in feature.features.values():
        # A 3-level hierarchy (gene -> primary_transcript -> miRNA -> exon) nests a
        # real Lifton_TRANS under a LiftOn_FEATURE, so dispatch on shape exactly as
        # _render_gene does -- recursing blindly into `.features` raised
        # AttributeError on the transcript and lost the whole subtree.
        if hasattr(sub_feature, "exons"):
            out.write(_render_trans(sub_feature))
        else:
            out.write(_render_feature(sub_feature))
    return out.getvalue()


def _render_gene(gene):
    """Render a gene plus every independently valid child hierarchy.

    A malformed transcript invalidates that transcript block, not its valid
    siblings. The parent and all retained children are still published in one
    write, so readers never observe an orphan or a half-written transcript.
    """
    out = io.StringIO()
    if not gene.tmp:
        out.write(gff3_writer.format_feature(gene.entry) + "\n")
    emitted_children = []
    failures = []
    for child in gene.transcripts.values():
        try:
            if hasattr(child, "exons"):
                out.write(_render_trans(child))
            else:
                out.write(_render_feature(child))
        except Exception as exc:
            failures.append((getattr(child.entry, "id", "<unknown>"), exc))
            continue
        emitted_children.append(child)
    return out.getvalue(), emitted_children, failures


def _record_transcript_stats(gene, transcripts_stats_dict, children) -> None:
    """Commit statistics only after the complete gene block was written."""
    for trans in children:
        if not hasattr(trans, "ref_tran_id"):
            continue
        if gene.is_protein_coding and trans.entry.featuretype == "mRNA":
            trans_type = "coding"
        elif gene.is_non_coding and trans.entry.featuretype in {
                "ncRNA", "nc_RNA", "lncRNA", "lnc_RNA"}:
            trans_type = "non-coding"
        else:
            trans_type = "other"
        if trans.ref_tran_id is None:
            continue
        bucket = transcripts_stats_dict[trans_type]
        bucket[trans.ref_tran_id] = bucket.get(trans.ref_tran_id, 0) + 1


def write_gene(gene, fw, transcripts_stats_dict):
    """Atomically publish a complete gene hierarchy.

    Formatting every row into a private buffer first prevents a bad parent,
    exon, or CDS from leaving an orphaned partial hierarchy in the final GFF3.
    The boolean result lets the ordered pipeline count only emitted models.
    """
    try:
        block, emitted_children, failures = _render_gene(gene)
        if gene.transcripts and not emitted_children:
            first_id, first_exc = failures[0]
            raise LiftOnError(
                f"all child hierarchies failed; first failure {first_id}: "
                f"{first_exc}"
            )
        fw.write(block)
    except Exception as exc:
        logger.log_error(
            f"Failed to write GENE hierarchy {gene.entry.id}: {exc}"
        )
        gene._serialization_failures = [
            (getattr(gene.entry, "id", "<unknown>"), str(exc))
        ]
        return False
    gene._serialization_failures = []
    for child_id, exc in failures:
        child_kind = "TRANSCRIPT" if any(
            getattr(c.entry, "id", None) == child_id and hasattr(c, "exons")
            for c in gene.transcripts.values()
        ) else "FEATURE"
        logger.log_error(
            f"Failed to write {child_kind} hierarchy {child_id} under "
            f"{gene.entry.id}: {exc}"
        )
        gene._serialization_failures.append((child_id, str(exc)))
    _record_transcript_stats(gene, transcripts_stats_dict, emitted_children)
    return True


def write_feature(feature, fw):
    try:
        fw.write(_render_feature(feature))
    except Exception as exc:
        logger.log_error(
            f"Failed to write FEATURE hierarchy {feature.entry.id}: {exc}"
        )
        return False
    return True


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
    #
    # GH #28: harmonize the transcript featuretype. LiftOn preserves the reference
    # type on the Liftoff/chaining path (often the generic "transcript") while the
    # miniprot path emits "mRNA", so a single output labels the same kind of feature
    # two ways. A transcript that HAS a CDS is, by SO definition, an mRNA -> canonicalize
    # coding transcripts to "mRNA" so the output is consistent regardless of which path
    # produced them (this also fixes the coding/other stats classification in write_gene,
    # which keys on featuretype == "mRNA"). Non-coding transcripts keep their type.
    try:
        fw.write(_render_trans(trans))
    except (OSError, ValueError, TypeError, AttributeError, LiftOnError) as e:
        logger.log_error(
            f"Failed to write TRANSCRIPT entry {trans.entry.id}: {e}"
        )
        return False
    return True


def write_exon(exon, fw):
    # Route through canonical writer (V5.4-V5.6, V5.9).
    try:
        fw.write(_render_exon(exon))
    except Exception as e:
        logger.log_error(f"Failed to write EXON entry {exon.entry.id}: {e}")
        return False
    return True


def write_cds(cds, fw):
    # Route through canonical writer (V5.4-V5.6, V5.9).
    try:
        fw.write(_render_cds(cds))
    except Exception as e:
        logger.log_error(f"Failed to write CDS entry {cds.entry.id}: {e}")
        return False
    return True
