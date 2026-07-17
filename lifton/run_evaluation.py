"""Evaluation-mode helpers.

Evaluation consumes an existing target annotation and writes per-transcript
identity/status rows.  It must not apply the copy-number or interval-tree
mutations used by the annotation-producing pipeline.
"""

import copy

import gffutils

from lifton import lifton_class, lifton_utils


_MISSING_FEATURE_ERRORS = (KeyError, gffutils.exceptions.FeatureNotFoundError)


def _chm13_mode(args):
    return bool(getattr(args, "evaluation_liftoff_chm13", False))


def _logical_id(feature_id, prefix):
    """Return the unprefixed ID used by evaluation's sequence dictionaries."""
    if prefix and feature_id.startswith(prefix):
        return feature_id[len(prefix):]
    return feature_id


def _copy_base(feature_id):
    """Return a possible LiftOn/Liftoff copy base, or ``None``.

    The base is only a candidate.  Lookup always tries the exact ID first, so
    natural reference IDs such as ``TX_1`` cannot be mistaken for copies.
    """
    base, separator, suffix = feature_id.rpartition("_")
    if separator and base and suffix.isdigit():
        return base
    return None


def _lookup_reference_feature(ref_db, feature_id, prefix="", allow_copy_suffix=True):
    """Resolve a reference feature without swallowing unexpected failures.

    CHM13 evaluation presents IDs without RefSeq's ``gene-``/``rna-``
    prefixes.  Copied target annotations may additionally append ``_<n>``.
    Exact lookup precedes copy-suffix fallback to preserve genuine numeric
    suffixes in reference IDs.
    """
    if feature_id is None:
        return None, None

    logical_id = _logical_id(str(feature_id), prefix)
    candidates = [logical_id]
    copy_base = _copy_base(logical_id) if allow_copy_suffix else None
    if copy_base is not None:
        candidates.append(copy_base)

    for candidate in candidates:
        key = f"{prefix}{candidate}" if prefix else candidate
        try:
            return candidate, ref_db[key]
        except _MISSING_FEATURE_ERRORS:
            continue
    return None, None


def _isolated_ref_features(ref_features_dict, ref_gene_id):
    """Clone only mutable copy-counter state needed by ``Lifton_GENE``."""
    isolated = dict(ref_features_dict)
    if ref_gene_id in isolated:
        isolated[ref_gene_id] = copy.copy(isolated[ref_gene_id])
    return isolated


def _feature_sort_key(feature):
    return (
        str(getattr(feature, "seqid", "")),
        int(getattr(feature, "start", 0)),
        int(getattr(feature, "end", 0)),
        str(getattr(feature, "featuretype", "")),
        str(getattr(feature, "id", "")),
    )


def _children(tgt_db, locus, *, featuretype=None, level=None):
    kwargs = {"order_by": "start"}
    if featuretype is not None:
        kwargs["featuretype"] = featuretype
    if level is not None:
        kwargs["level"] = level
    return sorted(list(tgt_db.children(locus, **kwargs)), key=_feature_sort_key)


def initialize_lifton_gene_eval(locus, ref_db, tree_dict, ref_features_dict,
                                args, with_exons=False):
    """Initialize an evaluation-only gene without mutating pipeline state."""
    ref_gene_id, ref_trans_id = lifton_utils.get_ref_ids_liftoff(
        ref_features_dict, locus.id, None,
    )
    if ref_gene_id is None:
        return None, None, None

    prefix = "gene-" if _chm13_mode(args) else ""
    ref_gene_id, ref_feature = _lookup_reference_feature(
        ref_db, ref_gene_id, prefix=prefix,
    )
    if ref_feature is None:
        return None, None, None

    # Lifton_GENE normally increments a global copy counter and inserts into a
    # shared interval tree.  Neither side effect belongs in read-only
    # evaluation, so give it isolated state for both.
    isolated_features = _isolated_ref_features(ref_features_dict, ref_gene_id)
    lifton_gene = lifton_class.Lifton_GENE(
        ref_gene_id,
        copy.deepcopy(locus),
        copy.deepcopy(ref_feature.attributes),
        {},
        isolated_features,
        args,
        tmp=with_exons,
    )
    return lifton_gene, ref_gene_id, ref_trans_id


def lifton_add_trans_exon_cds_eval(lifton_gene, locus, ref_db, tgt_db,
                                   ref_trans_id, args):
    """Add one target transcript hierarchy for evaluation.

    The target database features are deep-copied because LiftOn constructors
    rewrite IDs and Parent attributes.  If a valid target has CDS records but
    omits explicit exons, CDS spans are used as evaluation-only exon spans.
    """
    prefix = "rna-" if _chm13_mode(args) else ""
    ref_trans_id, ref_feature = _lookup_reference_feature(
        ref_db, ref_trans_id, prefix=prefix,
    )
    if ref_feature is None:
        return None, 0

    lifton_trans = lifton_gene.add_transcript(
        ref_trans_id,
        copy.deepcopy(locus),
        copy.deepcopy(ref_feature.attributes),
    )

    exons = _children(tgt_db, locus, featuretype="exon")
    cds_features = _children(
        tgt_db, locus, featuretype=("CDS", "stop_codon"),
    )
    if not exons:
        # GFF3 permits CDS children without separate exon records.  LiftOn's
        # sequence extraction operates on exon objects, so create transient
        # exon views from CDS segments without altering the target database.
        for cds in _children(tgt_db, locus, featuretype="CDS"):
            exon = copy.deepcopy(cds)
            exon.featuretype = "exon"
            exon.frame = "."
            lifton_gene.add_exon(lifton_trans.entry.id, exon)
    else:
        for exon in exons:
            lifton_gene.add_exon(lifton_trans.entry.id, copy.deepcopy(exon))

    for cds in cds_features:
        lifton_gene.add_cds(lifton_trans.entry.id, copy.deepcopy(cds))
    return lifton_trans, len(cds_features)


def evaluation(lifton_gene, locus, ref_db, tgt_db, tree_dict, tgt_fai,
               ref_proteins, ref_trans, ref_features_dict, fw_score, args,
               ENTRY_FEATURE=False):
    """Evaluate every transcript below ``locus`` and return its LiftOn gene.

    Unmapped target genes/transcripts are skipped.  Unexpected database or
    alignment errors still propagate so real input/tool failures remain
    visible to callers.
    """
    direct_exons = _children(tgt_db, locus, featuretype="exon", level=1)
    direct_cdss = _children(tgt_db, locus, featuretype="CDS", level=1)
    has_transcript_structure = bool(direct_exons or direct_cdss)

    ref_gene_id = None
    if lifton_gene is None:
        if not ENTRY_FEATURE:
            return None
        lifton_gene, ref_gene_id, _ = initialize_lifton_gene_eval(
            locus,
            ref_db,
            tree_dict,
            ref_features_dict,
            args,
            with_exons=has_transcript_structure,
        )
        if lifton_gene is None or lifton_gene.ref_gene_id is None:
            return None

    if not has_transcript_structure:
        for feature in _children(tgt_db, locus, level=1):
            if feature.featuretype in {"exon", "CDS", "stop_codon"}:
                continue
            evaluation(
                lifton_gene,
                feature,
                ref_db,
                tgt_db,
                tree_dict,
                tgt_fai,
                ref_proteins,
                ref_trans,
                ref_features_dict,
                fw_score,
                args,
            )
        return lifton_gene

    if ENTRY_FEATURE:
        ref_trans_id = ref_gene_id or lifton_gene.ref_gene_id
    else:
        _, ref_trans_id = lifton_utils.get_ref_ids_liftoff(
            ref_features_dict, lifton_gene.ref_gene_id, locus.id,
        )

    lifton_status = lifton_class.Lifton_Status()
    lifton_trans, _ = lifton_add_trans_exon_cds_eval(
        lifton_gene, locus, ref_db, tgt_db, ref_trans_id, args,
    )
    if lifton_trans is None:
        return lifton_gene

    canonical_ref_trans_id = lifton_trans.ref_tran_id
    lifton_gene.orf_search_protein(
        lifton_trans.entry.id,
        canonical_ref_trans_id,
        tgt_fai,
        ref_proteins,
        ref_trans,
        lifton_status,
        eval_only=True,
        eval_liftoff_chm13=_chm13_mode(args),
    )
    lifton_utils.print_lifton_status(
        lifton_trans.entry.id, locus, lifton_status,
        DEBUG=bool(getattr(args, "debug", False)),
    )
    lifton_utils.write_lifton_eval_status(
        fw_score, lifton_trans.entry.id, locus, lifton_status,
    )
    return lifton_gene
