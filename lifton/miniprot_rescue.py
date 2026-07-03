"""Iteration 23 — clean separate-pass miniprot-only rescue.

The default LiftOn pipeline never falls back to a miniprot-only gene model for a
coding gene the DNA lift (Liftoff / Step 7) missed ENTIRELY and that Step 8 then
dropped for falling outside the tight ``-min_miniprot``/``-max_miniprot`` length
band. The Iteration-22 in-loop rescue (an ``elif`` inside
``run_miniprot.process_miniprot``) recovered those genes but was NOT
default-ready: building a ``Lifton_GENE`` adds it to the SHARED Step-8
suppression ``tree_dict`` (``lifton_class.Lifton_GENE.__init__``), so a rescue
emitted mid-loop suppressed a LATER overlapping candidate the default would have
emitted (the ``off ⊄ on`` swap), and a multi-hit ref transcript got a redundant
second model.

This module runs the rescue as a SEPARATE POST-LIFT PASS, AFTER the Step-7 and
Step-8 loops have fully closed. Because no rescue mutates ``tree_dict`` during
any default decision, the default Step-7+8 output is byte-identical between
flag-OFF and flag-ON, so ``on = default ∪ rescues ⊇ default = off`` (off ⊆ on by
construction). A ref-gene-id dedup set (the genes Step 7 + Step 8 already
emitted) plus the now-final suppression tree give 0-redundant emission. Mirrors
the shipped lifton2 design (``lifton2/lifton2/miniprot_rescue.py``).

The whole pass is gated behind ``args.miniprot_rescue``; the caller never imports
this module when the flag is OFF, so the default path is provably inert.
"""
import os
import sys

from intervaltree import Interval

from lifton import run_miniprot, lifton_utils, logger


# Divergence-adaptive protein-identity floor (PROMOTED to default ON). The
# miniprot-only rescue (Iteration 23) used a FIXED floor (default 0.5). When the
# DNA lift placed only a small fraction of the genes miniprot ALSO found, the pair
# is divergent and the rescue is the dominant recall lever, so lower the floor
# toward ADAPT_FLOOR_MIN to admit more genuinely-missing genes; above ADAPT_R_HIGH
# recall (same/close-species) the floor stays at the user's base (rescue is
# marginal there -> keep precision). Mirrors lifton2's shipped design
# (lifton2/lifton2/miniprot_rescue.py). Thresholds env-overridable for the A/B.
ADAPT_FLOOR_MIN = float(os.environ.get("LIFTON_RESCUE_FLOOR_MIN") or 0.30)
ADAPT_R_LOW = float(os.environ.get("LIFTON_RESCUE_R_LOW") or 0.10)
ADAPT_R_HIGH = float(os.environ.get("LIFTON_RESCUE_R_HIGH") or 0.50)


def _adaptive_floor_on(args):
    """Whether the divergence-adaptive floor is active. Env
    ``LIFTON_RESCUE_ADAPTIVE_FLOOR`` wins (1/0), else ``args.adaptive_rescue_floor``
    (default True -- PROMOTED to default ON). Opt out with
    ``--no-adaptive-rescue-floor`` / ``LIFTON_RESCUE_ADAPTIVE_FLOOR=0`` (restores
    the fixed ``-miniprot_rescue_min_id`` floor)."""
    env = os.environ.get("LIFTON_RESCUE_ADAPTIVE_FLOOR")
    if env is not None:
        return env not in ("0", "", "false", "False")
    return bool(getattr(args, "adaptive_rescue_floor", True))


def _adaptive_floor(recall, base):
    """Map DNA-lift gene recall -> PI floor in ``[min(ADAPT_FLOOR_MIN, base), base]``:
    ``recall<=ADAPT_R_LOW`` -> the low floor (very-distant, recall is the priority);
    ``recall>=ADAPT_R_HIGH`` -> ``base`` (well-synteny, keep precision); linear in
    between. Never rises above ``base``."""
    lo = min(ADAPT_FLOOR_MIN, base)
    if recall <= ADAPT_R_LOW:
        return lo
    if recall >= ADAPT_R_HIGH:
        return base
    frac = (recall - ADAPT_R_LOW) / (ADAPT_R_HIGH - ADAPT_R_LOW)
    return lo + frac * (base - lo)


def _dna_lift_recall(universe, emitted):
    """Fraction of the miniprot-found ref genes (``universe``) that the DNA lift +
    default Step 8 already placed (``emitted``) -- a cheap, divergence-sensitive
    proxy: ~1 on same-species, low on very-distant. Returns 1.0 (=> base floor) if
    the universe is empty (no miniprot-found ref genes -> the rescue is marginal)."""
    universe = set(universe)
    if not universe:
        return 1.0
    return len(universe & set(emitted)) / len(universe)


def rescue_miniprot_only_pass(m_feature_db, ref_db, tree_dict, tgt_fai,
                              ref_proteins, ref_trans, ref_features_dict,
                              m_id_2_ref_id_trans_dict, ref_features_len_dict,
                              ref_trans_exon_num_dict, ref_features_reverse_dict,
                              emitted_ref_gene_ids, fw, fw_score,
                              transcripts_stats_dict, args):
    """Emit a miniprot-only gene model for each reference coding gene the DNA
    lift + default Step 8 missed but miniprot found at a non-overlapping locus,
    protein-identity >= floor, non-redundant.

    Runs strictly after the Step-8 loop closes, so it never perturbs a default
    decision. Returns the number of genes added; inert (returns 0) when
    ``m_feature_db is None``.
    """
    if m_feature_db is None:
        return 0

    floor = float(getattr(args, "miniprot_rescue_min_id", 0.5))
    added = 0

    # Deterministic iteration: sort the miniprot mRNAs by (seqid, start, end,
    # ID) so -t8 == -t1 and the cascading suppression (each emit mutates
    # tree_dict, so a later overlapping candidate is suppressed) is fully
    # reproducible run-to-run.
    try:
        mtranscripts = sorted(
            m_feature_db.features_of_type('mRNA'),
            key=lambda m: (m.seqid, m.start, m.end, m.attributes["ID"][0]),
        )
    except Exception as e:
        logger.log_error(f"miniprot-only rescue: failed to enumerate mRNAs: {e}")
        return 0

    # Divergence-adaptive PI floor (default ON): the DNA-lift gene recall is the
    # fraction of the miniprot-found ref genes already emitted by Step 7 + Step 8
    # (~1 same-species, low very-distant); lower the floor toward ADAPT_FLOOR_MIN
    # as recall drops. off subset of on is preserved: this only LOWERS the floor
    # within the already-additive separate pass, so it can only ADD rescues, never
    # remove or change a default-emitted feature.
    if _adaptive_floor_on(args):
        universe = set()
        for _m in mtranscripts:
            try:
                _gid, _tid = lifton_utils.get_ref_ids_miniprot(
                    ref_features_reverse_dict, _m.attributes["ID"][0],
                    m_id_2_ref_id_trans_dict)
            except Exception:
                _gid = None
            if _gid is not None:
                universe.add(_gid)
        recall = _dna_lift_recall(universe, emitted_ref_gene_ids)
        adapted = _adaptive_floor(recall, floor)
        if abs(adapted - floor) > 1e-9:
            sys.stderr.write(
                f"[LiftOn] miniprot-only rescue: adaptive PI floor {adapted:.2f} "
                f"(DNA-lift gene recall {recall:.3f}, base {floor:.2f}).\n")
            sys.stderr.flush()
        floor = adapted

    for mtrans in mtranscripts:
        try:
            mtrans_id = mtrans.attributes["ID"][0]

            # (1) resolve canonical ref ids (same namespace as Step 7's
            #     get_ref_ids_liftoff harvest -> dedup is comparable).
            ref_gene_id, ref_trans_id = lifton_utils.get_ref_ids_miniprot(
                ref_features_reverse_dict, mtrans_id, m_id_2_ref_id_trans_dict)
            if ref_gene_id is None or ref_trans_id is None:
                continue

            # (2) DEDUP: skip a ref gene already emitted by Step 7 (DNA lift),
            #     Step 8 (default miniprot), or an earlier rescue in this pass.
            #     This is what makes the pass 0-redundant + off ⊆ on.
            if ref_gene_id in emitted_ref_gene_ids:
                continue

            # (3) protein availability (mirror process_miniprot:381)
            if ref_trans_id not in ref_proteins or ref_trans_id not in ref_trans:
                continue

            # (4) overlap suppression vs the FINAL tree_dict: a gene already
            #     sits at this target locus (DNA lift, default Step 8, or a
            #     prior rescue). check_ovps_ratio reads tree_dict per-seqid.
            mtrans_interval = Interval(mtrans.start, mtrans.end, mtrans_id)
            if lifton_utils.check_ovps_ratio(mtrans, mtrans_interval,
                                             args.overlap, tree_dict):
                continue

            # (5) processed-pseudogene filter (mirror process_miniprot:386):
            #     1 CDS in miniprot but >1 exon in the reference.
            if (len(list(m_feature_db.children(mtrans, featuretype='CDS'))) == 1
                    and ref_trans_exon_num_dict.get(ref_trans_id, 0) > 1):
                continue

            # (6) length-ratio sanity band (the WIDE rescue band; the PI floor
            #     at step 8 is the real quality gate). Guards against a
            #     catastrophically mis-scaled miniprot hit.
            ref_len = ref_features_len_dict.get(ref_gene_id)
            if not ref_len:
                continue
            ratio = (mtrans.end - mtrans.start + 1) / ref_len
            if not run_miniprot._miniprot_rescue_band_ok(ratio, args):
                continue

            # (7) build + score. NOTE: constructing the Lifton_GENE adds it to
            #     tree_dict (lifton_class.Lifton_GENE.__init__), so a later
            #     overlapping rescue in this loop is suppressed at step (4) --
            #     cascading dedup, lifton2-style. Safe: this is post-loop, so
            #     no default Step-8 decision can be perturbed.
            lifton_gene, lifton_trans, transcript_id, lifton_status = \
                run_miniprot.lifton_miniprot_with_ref_protein(
                    mtrans, m_feature_db, ref_db.db_connection, ref_gene_id,
                    ref_trans_id, tgt_fai, ref_proteins, ref_trans, tree_dict,
                    ref_features_dict, args)

            # (8) PI floor -- the quality gate (mirror the in-loop pre-ORF
            #     check at process_miniprot:401, on the miniprot alignment
            #     identity set by lifton_miniprot_with_ref_protein).
            if lifton_status.lifton_aa < floor:
                # The build already inserted this candidate into tree_dict; the
                # stray interval only suppresses OTHER candidates at the same
                # missed locus (competing reps of the same gene), which is the
                # desired dedup -- it can never resurrect a lost default
                # feature, so off ⊆ on is unaffected.
                continue

            # (9) tag attrs BEFORE the ORF/status tail (parity with the in-loop
            #     rescue ordering at process_miniprot:403-404).
            lifton_gene.transcripts[transcript_id].entry.attributes[
                "miniprot_annotation_ratio"] = [f"{ratio:.3f}"]
            lifton_gene.transcripts[transcript_id].entry.attributes[
                "lifton_rescue"] = ["miniprot_only"]

            # (10) ORF search + status tail (mirror process_miniprot:409-413).
            lifton_gene.orf_search_protein(lifton_trans.entry.id, ref_trans_id,
                                           tgt_fai, ref_proteins, ref_trans,
                                           lifton_status)
            lifton_utils.print_lifton_status(transcript_id, mtrans,
                                             lifton_status, DEBUG=args.debug)
            lifton_gene.add_lifton_gene_status_attrs("miniprot")
            lifton_gene.add_lifton_trans_status_attrs(transcript_id, lifton_status)
            lifton_utils.write_lifton_status(fw_score, transcript_id, mtrans,
                                             lifton_status)

            # (11) emit + record dedup (so a later hit for the same ref gene
            #      and a later overlapping rescue are both skipped).
            lifton_gene.write_entry(fw, transcripts_stats_dict)
            emitted_ref_gene_ids.add(ref_gene_id)
            added += 1
        except Exception as e:
            logger.log_error(f"miniprot-only rescue error ({mtrans.id}): {e}")

    if added:
        sys.stderr.write(
            f"\n[LiftOn] miniprot-only rescue: {added} gene(s) added "
            f"(protein-identity floor {floor:.2f}).\n")
        sys.stderr.flush()
    return added
