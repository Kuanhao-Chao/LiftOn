import subprocess
import os, copy
from lifton import lifton_class, logger, fix_trans_annotation, lifton_utils, run_miniprot, align
from lifton.liftoff import liftoff_main
from lifton.liftoff.tests import test_basic, test_advanced
from intervaltree import Interval, IntervalTree

def run_liftoff(output_dir, args):
    """
        This function runs liftoff.

        Parameters:
        - output_dir: output directory
        - args: Liftoff arguments

        Returns:
        liftoff_annotation: liftoff annotation file path
    """
    liftoff_args = copy.deepcopy(args)
    liftoff_outdir = output_dir + "liftoff/"    
    os.makedirs(liftoff_outdir, exist_ok=True)
    liftoff_annotation = liftoff_outdir + "liftoff.gff3"
    liftoff_args.output = liftoff_annotation
    liftoff_main.run_all_liftoff_steps(liftoff_args)
    if args.polish:
        liftoff_annotation += "_polished"
    # test_basic.test_yeast(liftoff_outdir + "test_basic/")
    # test_advanced.test_yeast(liftoff_outdir + "test_advance/")
    return liftoff_annotation


def initialize_lifton_gene(locus, ref_db, tree_dict, ref_features_dict, with_exons=False):
    ref_gene_id, ref_trans_id = lifton_utils.get_ref_ids_liftoff(ref_features_dict, locus.id, None)
    lifton_gene = lifton_class.Lifton_GENE(ref_gene_id, copy.deepcopy(locus), copy.deepcopy(ref_db[ref_gene_id].attributes), tree_dict, ref_features_dict, tmp=with_exons)
    return lifton_gene, ref_gene_id, ref_trans_id


def initialize_lifton_miniprot_gene(miniprot_trans, ref_gene_id, ref_db, tree_dict, ref_features_dict):
    m_gene_feature = copy.deepcopy(miniprot_trans.entry)
    m_gene_feature.featuretype = "gene"
    lifton_gene = lifton_class.Lifton_GENE(ref_gene_id, m_gene_feature, copy.deepcopy(ref_db[ref_gene_id].attributes), tree_dict, ref_features_dict)    
    return lifton_gene


def lifton_add_trans_exon_cds(lifton_gene, locus, ref_db, l_feature_db, ref_trans_id):
    lifton_trans = lifton_gene.add_transcript(ref_trans_id, copy.deepcopy(locus), copy.deepcopy(ref_db[ref_trans_id].attributes))
    exons = l_feature_db.children(locus, featuretype='exon', order_by='start')
    for exon in list(exons):
        lifton_gene.add_exon(lifton_trans.entry.id, exon)
    cdss = l_feature_db.children(locus, featuretype=('CDS', 'stop_codon'), order_by='start') 
    cdss_list = list(cdss)
    for cds in cdss_list:
        lifton_gene.add_cds(lifton_trans.entry.id, cds)
    return lifton_trans, len(cdss_list)


def process_liftoff_with_protein(locus, lifton_gene, lifton_trans,
                                 ref_id_2_m_id_trans_dict, m_feature_db, tree_dict,
                                 tgt_fai, ref_trans_id, ref_proteins, ref_trans,
                                 fw_chain, write_chains, lifton_status, DEBUG):
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
        - ref_gene_id: reference gene ID
        - ref_trans_id: reference transcript ID
        - ref_proteins: reference proteins dictionary
        - ref_trans: reference transcript dictionary
        - fw_chain: file writer for chains
        - write_chains: write chains or not
        - lifton_status: Lifton_Status instance
        - DEBUG: debug mode

        Returns:
        process_locus: boolean value if the locus should be processed or not
    """
    # Liftoff alignment
    liftoff_aln = lifton_utils.LiftOn_liftoff_alignment(lifton_trans, locus, tgt_fai, ref_proteins, ref_trans_id, lifton_status)
    # miniprot alignment
    miniprot_aln, has_valid_miniprot = lifton_utils.LiftOn_miniprot_alignment(locus.seqid, locus, ref_id_2_m_id_trans_dict, m_feature_db, tree_dict, tgt_fai, ref_proteins, ref_trans_id, lifton_status)
    if liftoff_aln is None:
        # There is no reference protein -> just keep Liftoff annotation
        lifton_status.annotation = "Liftoff_no_ref_protein"
        lifton_status.status = ["no_ref_protein"]
    elif liftoff_aln.identity == 1:
        # Liftoff protein annotation is perfect
        lifton_status.lifton_aa = 1
    elif liftoff_aln.identity < 1:
        # Liftoff protein annotation is not perfect
        if has_valid_miniprot:
            lifton_status.annotation = "LiftOn_chaining_algorithm"
            cds_list, chains = fix_trans_annotation.chaining_algorithm(liftoff_aln, miniprot_aln, tgt_fai, DEBUG)
            if write_chains:
                lifton_utils.write_lifton_chains(fw_chain, lifton_trans.entry.id, chains)
            lifton_gene.update_cds_list(lifton_trans.entry.id, cds_list)
        else:
            lifton_status.annotation = "Liftoff_truncated"
        lifton_trans_aln, lifton_aa_aln = lifton_gene.fix_truncated_protein(lifton_trans.entry.id, ref_trans_id, tgt_fai, ref_proteins, ref_trans, lifton_status)


# def process_liftoff_without_protein(locus, lifton_gene, lifton_trans, ref_id_2_m_id_trans_dict, m_feature_db, tree_dict, tgt_fai, ref_gene_id, ref_trans_id, ref_proteins, ref_trans, lifton_status, DEBUG):
#     """
#         This function process liftoff annotation without protein.
#         (1) miniprot has protein: apply miniprot annotation to Liftoff
#         (2) Liftoff & miniprot don't have protein: lncRNA

#         Parameters:
#         - locus: gffutils feature instance
#         - ref_id_2_m_id_trans_dict: reference id to miniprot transcript ids dictionary
#         - m_feature_db: miniprot feature database
#         - tree_dict: intervaltree dictionary for each chromosome
#         - tgt_fai: target fasta index
#         - ref_gene_id: reference gene ID
#         - ref_trans_id: reference transcript ID
#         - ref_proteins: reference proteins dictionary
#         - ref_trans: reference transcript dictionary
#         - ref_features_dict: reference features dictionary
#         - ref_db: reference database
#         - lifton_status: Lifton_Status instance
#         - DEBUG: debug mode

#         Returns:
#         process_locus: boolean value if the locus should be processed or not
#     """
#     miniprot_aln = None
#     miniprot_trans = None
#     has_valid_miniprot = False
#     if (ref_trans_id in ref_id_2_m_id_trans_dict.keys()) and (ref_trans_id in ref_proteins.keys()):
#         m_ids = ref_id_2_m_id_trans_dict[ref_trans_id]
#         for m_id in m_ids:
#             ##################################################
#             # Check 1: Check if the miniprot transcript is overlapping the current gene locus
#             ##################################################
#             mtrans = m_feature_db[m_id]
#             mtrans_id = mtrans.attributes["ID"][0]
#             is_overlapped = lifton_utils.segments_overlap((mtrans.start, mtrans.end), (locus.start, locus.end))
#             # mtrans_interval = Interval(mtrans.start, mtrans.end, mtrans_id)
#             # is_overlapped = lifton_utils.check_ovps_ratio(mtrans, mtrans_interval, 0.7, tree_dict)
#             if not is_overlapped or mtrans.seqid != locus.seqid:
#                 # "Not overlapped"
#                 continue
#             ##################################################
#             # Check 2: reference overlapping status
#             #   1. Check it the transcript overlapping with the next gene
#             # Check the miniprot protein overlapping status
#             # The case I should not process the transcript 
#             #   1. The Liftoff does not overlap with other gene
#             #   2. The miniprot protein overlap the other gene
#             ##################################################
#             ovps_liftoff = tree_dict[locus.seqid].overlap(locus.start, locus.end)
#             ovps_miniprot = tree_dict[locus.seqid].overlap(mtrans.start, mtrans.end)
#             miniprot_cross_gene_loci = False
#             liftoff_set = set()
#             for ovp_liftoff in ovps_liftoff:
#                 liftoff_set.add(ovp_liftoff[2])
#             for ovp_miniprot in ovps_miniprot:
#                 if ovp_miniprot[2] not in liftoff_set:
#                     # Miniprot overlap to more genes
#                     miniprot_cross_gene_loci = True
#                     break
#             if miniprot_cross_gene_loci:
#                 continue
#             # Valid miniprot transcript exists => check if the miniprot transcript is valid
#             has_valid_miniprot = True

#             tmp_miniprot_trans = lifton_class.Lifton_TRANS(mtrans_id, "", "", 0, mtrans, {})
#             exons = m_feature_db.children(mtrans, featuretype=('CDS', 'stop_codon'), order_by='start')
#             for exon in list(exons):
#                 tmp_miniprot_trans.add_exon(exon)
#             cdss = m_feature_db.children(mtrans, featuretype=('CDS', 'stop_codon'), order_by='start') 
#             cds_num = 0
#             for cds in list(cdss):
#                 cds_num += 1
#                 tmp_miniprot_trans.add_cds(cds)
#             tmp_miniprot_aln = align.lifton_parasail_align("miniprot", tmp_miniprot_trans, mtrans, tgt_fai, ref_proteins, ref_trans_id)
#             if miniprot_aln == None or tmp_miniprot_aln.identity > lifton_status.miniprot:
#                 miniprot_trans = tmp_miniprot_trans
#                 miniprot_aln = tmp_miniprot_aln
#                 lifton_status.miniprot = miniprot_aln.identity
#     if has_valid_miniprot:
#         print("In")
#         # Create LifOn gene instance
#         # lifton_gene = initialize_lifton_miniprot_gene(miniprot_trans, ref_gene_id, ref_db, tree_dict, ref_features_dict)

#         # Add transcript instance
#         lifton_gene.transcripts[miniprot_trans.entry.id] = miniprot_trans
#         lifton_status.lifton_aa = miniprot_aln.identity
#         if miniprot_aln.identity == 1:
#             lifton_status.annotation =  "miniprot_identical"
#         elif miniprot_aln.identity < 1:
#             lifton_status.annotation =  "miniprot_truncated"
#         lifton_trans_aln, lifton_aa_aln = lifton_gene.fix_truncated_protein(miniprot_trans.entry.id, ref_trans_id, tgt_fai, ref_proteins, ref_trans, lifton_status)


def process_liftoff(lifton_gene, locus, ref_db, l_feature_db, ref_id_2_m_id_trans_dict, m_feature_db, tree_dict, tgt_fai, ref_proteins, ref_trans, ref_features_dict, fw_score, fw_chain, write_chains, DEBUG, ENTRY_FEATURE=False):
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

        Returns:
        lifton_gene: LiftOn gene instance
    """
    exon_children = list(l_feature_db.children(locus, featuretype='exon', level=1, order_by='start'))
    if lifton_gene is None and ENTRY_FEATURE:   
        # Gene (1st) features
        lifton_gene, ref_gene_id, ref_trans_id = initialize_lifton_gene(locus, ref_db, tree_dict, ref_features_dict, with_exons=len(exon_children)>0)
        if lifton_gene.ref_gene_id is None:
            # Cannot find the reference gene id
            return None
    
    if len(exon_children) == 0:
        parent_feature = None
        if ENTRY_FEATURE:
            # Gene (1st) features without direct exons 
            parent_feature = lifton_gene
        else:
            # Middle features without exons
            parent_feature = lifton_gene.add_feature(copy.deepcopy(locus))
        features = l_feature_db.children(locus, level=1)
        for feature in list(features):
            lifton_gene = process_liftoff(parent_feature, feature, ref_db, l_feature_db, ref_id_2_m_id_trans_dict, m_feature_db, tree_dict, tgt_fai, ref_proteins, ref_trans, ref_features_dict, fw_score, fw_chain, write_chains, DEBUG)
    else:
        if ENTRY_FEATURE:
            # Gene (1st) features with direct exons 
            ref_trans_id = ref_gene_id
        else:
            # Transcript features with direct exons 
            ref_gene_id, ref_trans_id = lifton_utils.get_ref_ids_liftoff(ref_features_dict, lifton_gene.entry.id, locus.id)
        lifton_status = lifton_class.Lifton_Status()
        lifton_trans, cds_num = lifton_add_trans_exon_cds(lifton_gene, locus, ref_db, l_feature_db, ref_trans_id)
        if cds_num > 0:
            process_liftoff_with_protein(locus, lifton_gene, lifton_trans,
                                        ref_id_2_m_id_trans_dict, m_feature_db, tree_dict,
                                        tgt_fai, ref_trans_id, ref_proteins, ref_trans,
                                        fw_chain, write_chains, lifton_status, DEBUG)
        # LiftOn DNA level alignment
        if ref_trans_id not in ref_trans.keys():
            # no transcript sequence
            return None
        _ = lifton_gene.align_trans_dna(lifton_trans.entry.id, ref_trans_id, tgt_fai, ref_trans, lifton_status, cds_num)
        lifton_utils.print_lifton_status(lifton_trans.entry.id, locus, lifton_status, DEBUG=DEBUG)
        lifton_gene.add_lifton_status_attrs(lifton_trans.entry.id, lifton_status)
        # Writing out LiftOn entries & scores
        lifton_utils.write_lifton_status(fw_score, lifton_trans.entry.id, locus, lifton_status)
    return lifton_gene