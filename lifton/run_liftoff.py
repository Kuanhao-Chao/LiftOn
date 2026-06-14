import subprocess
import os, copy, sys
import gffutils
from lifton import align, lifton_class, logger, lifton_utils, protein_maximization, run_miniprot
from lifton.exceptions import LiftOnInputError
from lifton.liftoff import liftoff_main
from lifton.liftoff.tests import test_basic, test_advanced
from intervaltree import Interval, IntervalTree

def run_liftoff(output_dir, ref_db, args):
    """Run Liftoff and return either a path to the resulting GFF3
    (legacy default) or an in-memory GFF3 bytes blob (Phase 8).

    The bytes-blob branch fires when ``args.inmemory_liftoff`` is
    True. The output is byte-identical to the legacy file because
    :mod:`lifton.liftoff.inmemory_emitter` shares its serialisation
    helpers with :mod:`lifton.liftoff.write_new_gff`.

    Parameters
    ----------
    output_dir : str
        Output directory for intermediate Liftoff artefacts
        (``unmapped_features.txt`` is still written here in both
        modes; only the final GFF3 is held in RAM in the in-memory
        mode).
    ref_db : gffutils.FeatureDB
        Reference annotation DB Liftoff lifts from.
    args : argparse.Namespace
        Must carry ``polish``; may carry ``inmemory_liftoff``
        (default False).

    Returns
    -------
    str | bytes
        Path to the Liftoff GFF3 output, OR an in-memory bytes blob
        when ``args.inmemory_liftoff`` is True.
    """
    liftoff_args = copy.deepcopy(args)
    liftoff_outdir = output_dir + "liftoff/"
    os.makedirs(liftoff_outdir, exist_ok=True)
    liftoff_annotation = liftoff_outdir + "liftoff.gff3"
    liftoff_args.output = liftoff_annotation
    liftoff_args.u = liftoff_outdir + "unmapped_features.txt"

    inmemory = bool(getattr(args, "inmemory_liftoff", False))
    # Phase 16 Tier 3: raise the recursion limit before delegating to
    # the vendored Liftoff library. Real-world reference annotations
    # (NCBI RefSeq, GENCODE) drive Liftoff's recursive feature-hierarchy
    # traversal deep enough to exceed Python's default 1000-frame limit
    # on the human/mouse benchmark datasets. 10× headroom is the right
    # shape for bounded-but-deep traversal while staying well within
    # the OS thread stack (Linux default 8 MB ≈ ~13K Python frames).
    # Pathological cycles will still surface — Tier 2's full-traceback
    # dump in the except below pinpoints the recursive function so a
    # precise cycle-guard fix can land in a follow-up phase.
    _orig_recursion_limit = sys.getrecursionlimit()
    sys.setrecursionlimit(max(_orig_recursion_limit, 10000))
    try:
        if inmemory:
            from lifton.liftoff import inmemory_emitter
            lifted_feature_list, feature_db, _ref_parent_order, _unmapped = \
                liftoff_main.run_all_liftoff_steps_inmemory(liftoff_args, ref_db)
            gff_bytes = inmemory_emitter.lifted_features_to_gff3_bytes(
                lifted_feature_list, liftoff_args, feature_db,
            )
        else:
            liftoff_main.run_all_liftoff_steps(liftoff_args, ref_db)
    except Exception as e:
        import traceback
        logger.log_error(f"Liftoff encountered a fatal error during native execution: {e}")
        # Phase 16 Tier 2: emit the full traceback so the deepest LiftOn /
        # vendored-Liftoff frame is visible in stderr. Without this the
        # exception's str() form (e.g. "maximum recursion depth exceeded
        # while calling a Python object" for a RecursionError) gives no
        # hint as to which file is recursing.
        logger.log_error(
            "Full Python traceback (deepest frame is the failing call):\n"
            + traceback.format_exc().rstrip("\n")
        )
        logger.log_error("LiftOn cannot proceed without a valid Liftoff baseline annotation.")
        sys.exit(1)
    finally:
        sys.setrecursionlimit(_orig_recursion_limit)

    if args.polish:
        liftoff_annotation += "_polished"
    # test_basic.test_yeast(liftoff_outdir + "test_basic/")
    # test_advanced.test_yeast(liftoff_outdir + "test_advance/")

    if inmemory:
        return gff_bytes
    return liftoff_annotation


def initialize_lifton_gene(locus, ref_db, tree_dict, ref_features_dict, args, with_exons=False):
    """
        This function initializes Lifton gene instance.

        Parameters:
        - locus: gffutils feature instance
        - ref_db: reference database
        - tree_dict: intervaltree dictionary for each chromosome
        - ref_features_dict: reference features dictionary
        - with_exons: True if the gene has exons, False otherwise

        Returns:
        lifton_gene: Lifton gene instance
        ref_gene_id: reference gene
        ref_trans_id: reference transcript
    """
    ref_gene_id, ref_trans_id = lifton_utils.get_ref_ids_liftoff(ref_features_dict, locus.id, None)
    lifton_gene = lifton_class.Lifton_GENE(ref_gene_id, copy.deepcopy(locus), copy.deepcopy(ref_db[ref_gene_id].attributes), tree_dict, ref_features_dict, args, tmp=with_exons)
    return lifton_gene, ref_gene_id, ref_trans_id


def lifton_add_trans_exon_cds(lifton_gene, locus, ref_db, l_feature_db, ref_trans_id):
    """
        This function adds transcript, exons, and CDSs to the Lifton gene instance.

        Parameters:
        - lifton_gene: Lifton gene instance
        - locus: gffutils feature instance
        - ref_db: reference database
        - l_feature_db: liftoff feature database
        - ref_trans_id: reference transcript ID

        Returns:
        lifton_trans: Lifton transcript instance
        len(cdss_list): number of CDSs
    """
    lifton_trans = lifton_gene.add_transcript(ref_trans_id, copy.deepcopy(locus), copy.deepcopy(ref_db[ref_trans_id].attributes))
    exons = l_feature_db.children(locus, featuretype='exon', order_by='start')
    for exon in list(exons):
        lifton_gene.add_exon(lifton_trans.entry.id, exon)
    cdss = l_feature_db.children(locus, featuretype=('CDS', 'stop_codon'), order_by='start') 
    cdss_list = list(cdss)
    for cds in cdss_list:
        lifton_gene.add_cds(lifton_trans.entry.id, cds)
    return lifton_trans, len(cdss_list)


def _snapshot_merge_state(lifton_gene, lifton_trans, lifton_status):
    """Capture the exact state that ``update_cds_list`` + a trial re-alignment
    mutate, so the default best-of-outcome verified-merge path can revert to
    pure Liftoff when the merge does not improve protein identity.

    ``update_cds_list`` replaces ``lifton_trans.exons`` and resets the
    transcript/gene boundaries (``Lifton_TRANS.update_boundaries`` sets
    start/end from the exon list; ``Lifton_GENE.update_boundaries`` grows the
    gene span). A trial ``orf_search_protein(eval_only=True)`` additionally
    mutates ``lifton_status`` (via ``align_coding_seq``/``align_trans_seq``/
    ``find_variants``) and appends to ``entry.attributes['mutation']``. We
    snapshot all of these.
    """
    st = lifton_status
    mut = lifton_trans.entry.attributes.get("mutation")
    return {
        "exons": copy.deepcopy(lifton_trans.exons),
        "trans_start": lifton_trans.entry.start,
        "trans_end": lifton_trans.entry.end,
        "gene_start": lifton_gene.entry.start,
        "gene_end": lifton_gene.entry.end,
        "mutation": copy.deepcopy(mut) if mut is not None else None,
        "status": (st.liftoff, st.miniprot, st.lifton_dna, st.lifton_aa,
                   st.eval_dna, st.eval_aa, st.annotation, list(st.status)),
    }


def _restore_status_and_mutation(lifton_trans, lifton_status, snap):
    """Undo the trial re-alignment's side effects on ``lifton_status`` and the
    ``mutation`` attribute so the canonical ``orf_search_protein`` re-derives
    identical state. In the accepted-merge case this keeps the emitted output
    consistent with a direct merge-then-ORF pass."""
    (liftoff, miniprot, lifton_dna, lifton_aa,
     eval_dna, eval_aa, annotation, status) = snap["status"]
    lifton_status.liftoff = liftoff
    lifton_status.miniprot = miniprot
    lifton_status.lifton_dna = lifton_dna
    lifton_status.lifton_aa = lifton_aa
    lifton_status.eval_dna = eval_dna
    lifton_status.eval_aa = eval_aa
    lifton_status.annotation = annotation
    lifton_status.status = list(status)
    if snap["mutation"] is None:
        lifton_trans.entry.attributes.pop("mutation", None)
    else:
        lifton_trans.entry.attributes["mutation"] = copy.deepcopy(snap["mutation"])


def _restore_merge_structure(lifton_gene, lifton_trans, snap):
    """Revert the transcript exon/CDS structure + transcript/gene boundaries to
    the pre-merge (pure Liftoff) snapshot."""
    lifton_trans.exons = snap["exons"]
    lifton_trans.entry.start = snap["trans_start"]
    lifton_trans.entry.end = snap["trans_end"]
    lifton_gene.entry.start = snap["gene_start"]
    lifton_gene.entry.end = snap["gene_end"]


def _optimize_fast_enabled():
    """Best-of-outcome fast path: skip the redundant Liftoff+ORF candidate when
    the merge already wins regardless of it (perfect protein, or a structural
    no-op vs Liftoff). Byte-identical to computing both candidates — pinned by
    TestOptimizeFlagScaffold + the drosophila byte-diff. On by default; set
    LIFTON_OPTIMIZE_FAST=0 to force computing both candidates."""
    env = os.environ.get("LIFTON_OPTIMIZE_FAST")
    if env is None:
        return True
    return env.strip().lower() not in ("", "0", "false", "no", "off")


def process_liftoff_with_protein(locus, lifton_gene, lifton_trans,
                                 ref_id_2_m_id_trans_dict, m_feature_db, tree_dict,
                                 tgt_fai, ref_trans_id, ref_proteins, ref_trans,
                                 fw_chain, write_chains, lifton_status, DEBUG,
                                 legacy_merge=False):
    """
        This function process liftoff annotation with protein.
        (1) Only Liftoff annotation: directly apply 
        (2) Liftoff annotation with miniprot annotation: Protein-maximization algorithm

        Parameters:
        - locus: gffutils feature instance
        - lifton_gene: Lifton gene instance
        - lifton_trans: Lifton transcript instance
        - ref_id_2_m_id_trans_dict: reference id to miniprot transcript ids dictionary
        - m_feature_db: miniprot feature database
        - tree_dict: intervaltree dictionary for each chromosome
        - tgt_fai: target fasta index
        - ref_trans_id: reference transcript ID
        - ref_proteins: reference proteins dictionary
        - ref_trans: reference transcript dictionary
        - fw_chain: file writer for chains
        - write_chains: write chains or not
        - lifton_status: Lifton_Status instance
        - DEBUG: debug mode
    """
    # Liftoff alignment
    liftoff_aln = lifton_utils.LiftOn_liftoff_alignment(lifton_trans, locus, tgt_fai, ref_proteins, ref_trans_id, lifton_status)
    # miniprot alignment
    miniprot_aln, has_valid_miniprot = lifton_utils.LiftOn_miniprot_alignment(locus.seqid, locus, ref_id_2_m_id_trans_dict, m_feature_db, tree_dict, tgt_fai, ref_proteins, ref_trans_id, lifton_status)
    orf_done = False
    if liftoff_aln is None:
        lifton_status.annotation = "no_ref_protein"
    elif liftoff_aln.identity == 1:
        # Liftoff protein annotation is perfect"
        lifton_status.lifton_aa = 1
    elif liftoff_aln.identity < 1:
        # Liftoff protein annotation is not perfect
        if has_valid_miniprot:
            lifton_status.annotation = "LiftOn_chaining_algorithm"
            cds_list, chains = protein_maximization.chaining_algorithm(
                liftoff_aln, miniprot_aln, tgt_fai, DEBUG)
            if write_chains:
                lifton_utils.write_lifton_chains(fw_chain, lifton_trans.entry.id, chains)
            if legacy_merge:
                # Legacy path (--legacy-merge): apply the chained CDS
                # unconditionally — the pre-promotion default, kept for
                # published-manuscript reproduction. Can silently frameshift
                # downstream CDS (the merge double-edge the new default fixes).
                lifton_gene.update_cds_list(lifton_trans.entry.id, cds_list)
            else:
                # DEFAULT (best-of-outcome verified merge). Compare the
                # final EMITTED protein of two candidates and keep the better:
                #   (1) no-merge: pure Liftoff + ORF-rescue  (== the no-miniprot
                #       default result for this transcript), and
                #   (2) merge: chained CDS + ORF-rescue.
                # Both are scored AFTER ORF-rescue via the same path the benchmark
                # re-aligns on, so the choice is exact; enabling miniprot becomes a
                # STRICT per-transcript win (emitted identity >= no-miniprot),
                # turning the protein-maximization merge from net-negative
                # (catastrophic backfires) into a guaranteed gain. Proven on
                # drosophila (+0.0012 mean PI) and mouse_to_rat (+0.0067,
                # 294 improved / 2 regressed).
                #
                # The merge candidate (2) is evaluated FIRST so the (expensive)
                # Liftoff + ORF-rescue candidate (1) — two parasail alignments and
                # possibly a 3-frame ORF scan — can be SKIPPED whenever the merge
                # already wins regardless of it:
                #   (a) merge reaches a perfect protein (>=1.0): Liftoff is <=1.0
                #       and ties go to merge, so merge is kept; or
                #   (b) the chained CDS is structurally identical to the Liftoff
                #       CDS: the two candidates are the same computation → tie →
                #       merge kept.
                # In both, the always-both-candidates code keeps the merge
                # candidate, so skipping candidate (1) is byte-identical (pinned by
                # TestOptimizeFlagScaffold + the drosophila byte-diff). Set
                # LIFTON_OPTIMIZE_FAST=0 to force computing both candidates.
                _merge_noop = (
                    sorted((c.start, c.end) for c in liftoff_aln.cds_children)
                    == sorted((lc.entry.start, lc.entry.end) for lc in cds_list))
                _snap0 = _snapshot_merge_state(lifton_gene, lifton_trans, lifton_status)
                # Candidate 2 (merge) first, from the raw pre-merge state.
                lifton_status.lifton_aa = 0
                lifton_gene.update_cds_list(lifton_trans.entry.id, cds_list, optimize=True)
                lifton_gene.orf_search_protein(
                    lifton_trans.entry.id, ref_trans_id, tgt_fai,
                    ref_proteins, ref_trans, lifton_status)
                _merge_outcome = lifton_status.lifton_aa
                if _optimize_fast_enabled() and (_merge_outcome >= 1.0 or _merge_noop):
                    # Merge wins without computing the Liftoff candidate; keep it
                    # (structure + status are already the merge candidate's).
                    lifton_status.annotation = "LiftOn_chaining_algorithm"
                else:
                    _snapM = _snapshot_merge_state(lifton_gene, lifton_trans, lifton_status)
                    # Candidate 1: Liftoff + ORF-rescue, from the raw pre-merge state.
                    _restore_merge_structure(lifton_gene, lifton_trans, _snap0)
                    _restore_status_and_mutation(lifton_trans, lifton_status, _snap0)
                    lifton_status.lifton_aa = 0
                    lifton_gene.orf_search_protein(
                        lifton_trans.entry.id, ref_trans_id, tgt_fai,
                        ref_proteins, ref_trans, lifton_status)
                    _liftoff_outcome = lifton_status.lifton_aa
                    if _merge_outcome >= _liftoff_outcome:
                        # Keep the merge candidate: restore its snapshot.
                        _restore_merge_structure(lifton_gene, lifton_trans, _snapM)
                        _restore_status_and_mutation(lifton_trans, lifton_status, _snapM)
                        lifton_status.annotation = "LiftOn_chaining_algorithm"
                    else:
                        # Keep Liftoff (current state IS the Liftoff candidate).
                        lifton_status.annotation = "Liftoff"
                # Candidate(s) already ran ORF-rescue; tell process_liftoff to skip
                # the canonical orf_search_protein so it is not run twice.
                orf_done = True
    return orf_done


def process_liftoff(lifton_gene, locus, ref_db, l_feature_db, ref_id_2_m_id_trans_dict, m_feature_db, tree_dict, tgt_fai, ref_proteins, ref_trans, ref_features_dict, fw_score, fw_chain, args, ENTRY_FEATURE=False, _visited=None):
    """
        This function processes liftoff annotation.

        Parameters:
        - lifton_gene: Lifton gene instance
        - locus: feature instance
        - ref_db: reference database
        - l_feature_db: liftoff feature database
        - ref_id_2_m_id_trans_dict: reference id to miniprot transcript ids dictionary
        - m_feature_db: miniprot feature database
        - tree_dict: intervaltree dictionary for each chromosome
        - tgt_fai: target fasta index
        - ref_proteins: reference protein dictionary
        - ref_trans: reference transcript dictionary
        - ref_features_dict: reference features dictionary
        - fw_score: file writer for scores
        - fw_chain: file writer for chains
        - write_chains: write chains or not
        - DEBUG: debug mode
        - ENTRY_FEATURE: True if the feature is the root feature for a gene locus
        - _visited: V5.2 cycle-detection set (internal — do not pass)

        Returns:
        lifton_gene: LiftOn gene instance
    """
    # V5.2 fix: detect circular Parent= cycles before they trigger
    # RecursionError. Track every locus.id seen on this descent.
    if _visited is None:
        _visited = set()
    locus_id_for_cycle = getattr(locus, "id", None)
    if locus_id_for_cycle is not None:
        if locus_id_for_cycle in _visited:
            raise LiftOnInputError(
                f"Circular Parent reference detected at feature "
                f"{locus_id_for_cycle!r}: this row's Parent chain forms a "
                f"cycle, which is not permitted by the GFF3 specification."
            )
        _visited = _visited | {locus_id_for_cycle}

    exon_children = list(l_feature_db.children(locus, featuretype='exon', level=1, order_by='start'))
    if lifton_gene is None and ENTRY_FEATURE:
        # Gene (1st) features
        lifton_gene, ref_gene_id, ref_trans_id = initialize_lifton_gene(locus, ref_db, tree_dict, ref_features_dict, args, with_exons=len(exon_children)>0)
        if lifton_gene.ref_gene_id is None: return None
    if len(exon_children) == 0:
        parent_feature = None
        if ENTRY_FEATURE: # Gene (1st) features without direct exons
            parent_feature = lifton_gene
        else: # Middle features without exons
            parent_feature = lifton_gene.add_feature(copy.deepcopy(locus))
        features = l_feature_db.children(locus, level=1)
        for feature in list(features):
            process_liftoff(parent_feature, feature, ref_db, l_feature_db, ref_id_2_m_id_trans_dict, m_feature_db, tree_dict, tgt_fai, ref_proteins, ref_trans, ref_features_dict, fw_score, fw_chain, args, _visited=_visited)
    else:
        if ENTRY_FEATURE: # Gene (1st) features with direct exons
            ref_trans_id = ref_gene_id
        else: # Transcript features with direct exons
            ref_gene_id, ref_trans_id = lifton_utils.get_ref_ids_liftoff(ref_features_dict, lifton_gene.entry.id, locus.id)
        lifton_status = lifton_class.Lifton_Status()
        lifton_status.annotation = "Liftoff"
        # V1.1a fix: narrow bare `except:` to the actual exceptions a
        # missing reference transcript can raise. KeyError is the
        # gffbase / dict-style miss; FeatureNotFoundError is gffutils.
        # KeyboardInterrupt and SystemExit now correctly propagate.
        try: # Test if the reference transcript exists. Skip if not.
            ref_db[ref_trans_id]
        except (KeyError, gffutils.exceptions.FeatureNotFoundError):
            return None
        lifton_trans, cds_num = lifton_add_trans_exon_cds(lifton_gene, locus, ref_db, l_feature_db, ref_trans_id)
        orf_done = False
        if cds_num > 0:
            orf_done = process_liftoff_with_protein(locus, lifton_gene, lifton_trans,
                                        ref_id_2_m_id_trans_dict, m_feature_db, tree_dict,
                                        tgt_fai, ref_trans_id, ref_proteins, ref_trans,
                                        fw_chain, args.write_chains, lifton_status, args.debug,
                                        legacy_merge=getattr(args, "legacy_merge", False))
        # The default best-of-outcome merge path already ran ORF-rescue on both
        # candidates (orf_done=True); avoid running it a second time here.
        # Under --legacy-merge orf_done stays False, so the canonical ORF-rescue
        # below runs after the unconditional merge (the pre-promotion behavior).
        if not args.no_orf_search and not orf_done:
            lifton_trans_aln, lifton_aa_aln = lifton_gene.orf_search_protein(lifton_trans.entry.id, ref_trans_id, tgt_fai, ref_proteins, ref_trans, lifton_status)
        # lifton_utils.print_lifton_status(lifton_trans.entry.id, locus, lifton_status, DEBUG=args.debug)
        lifton_gene.add_lifton_gene_status_attrs("Liftoff")
        lifton_gene.add_lifton_trans_status_attrs(lifton_trans.entry.id, lifton_status)
        lifton_utils.write_lifton_status(fw_score, lifton_trans.entry.id, locus, lifton_status)
    return lifton_gene