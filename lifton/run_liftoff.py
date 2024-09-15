import subprocess
import os, copy
from lifton import align, lifton_class, logger, lifton_utils, protein_maximization, run_miniprot
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
    liftoff_args.u = liftoff_outdir + "unmapped_features.txt"
    liftoff_main.run_all_liftoff_steps(liftoff_args)
    if args.polish:
        liftoff_annotation += "_polished"
    # test_basic.test_yeast(liftoff_outdir + "test_basic/")
    # test_advanced.test_yeast(liftoff_outdir + "test_advance/")
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
    if liftoff_aln is None:
        lifton_status.annotation = "no_ref_protein"
    elif liftoff_aln.identity == 1:
        # Liftoff protein annotation is perfect"
        lifton_status.lifton_aa = 1
    elif liftoff_aln.identity < 1:
        # Liftoff protein annotation is not perfect
        if has_valid_miniprot:
            lifton_status.annotation = "LiftOn_chaining_algorithm"
            cds_list, chains = protein_maximization.chaining_algorithm(liftoff_aln, miniprot_aln, tgt_fai, DEBUG)
            if write_chains:
                lifton_utils.write_lifton_chains(fw_chain, lifton_trans.entry.id, chains)
            lifton_gene.update_cds_list(lifton_trans.entry.id, cds_list)


def process_liftoff(lifton_gene, locus, ref_db, l_feature_db, ref_id_2_m_id_trans_dict, m_feature_db, tree_dict, tgt_fai, ref_proteins, ref_trans, ref_features_dict, fw_score, fw_chain, args, ENTRY_FEATURE=False):
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
            lifton_gene = process_liftoff(parent_feature, feature, ref_db, l_feature_db, ref_id_2_m_id_trans_dict, m_feature_db, tree_dict, tgt_fai, ref_proteins, ref_trans, ref_features_dict, fw_score, fw_chain, args)
    else:
        if ENTRY_FEATURE: # Gene (1st) features with direct exons 
            ref_trans_id = ref_gene_id
        else: # Transcript features with direct exons 
            ref_gene_id, ref_trans_id = lifton_utils.get_ref_ids_liftoff(ref_features_dict, lifton_gene.entry.id, locus.id)
        lifton_status = lifton_class.Lifton_Status()
        lifton_status.annotation = "Liftoff"
        try: # Test if the reference transcript exists. Skip if not.
            ref_db[ref_trans_id]
        except:
            return None
        lifton_trans, cds_num = lifton_add_trans_exon_cds(lifton_gene, locus, ref_db, l_feature_db, ref_trans_id)
        if cds_num > 0:
            process_liftoff_with_protein(locus, lifton_gene, lifton_trans,
                                        ref_id_2_m_id_trans_dict, m_feature_db, tree_dict,
                                        tgt_fai, ref_trans_id, ref_proteins, ref_trans,
                                        fw_chain, args.write_chains, lifton_status, args.debug)
        lifton_trans_aln, lifton_aa_aln = lifton_gene.orf_search_protein(lifton_trans.entry.id, ref_trans_id, tgt_fai, ref_proteins, ref_trans, lifton_status)
        
        # lifton_utils.print_lifton_status(lifton_trans.entry.id, locus, lifton_status, DEBUG=args.debug)
        lifton_gene.add_lifton_gene_status_attrs("Liftoff")
        lifton_gene.add_lifton_trans_status_attrs(lifton_trans.entry.id, lifton_status)
        lifton_utils.write_lifton_status(fw_score, lifton_trans.entry.id, locus, lifton_status)
    return lifton_gene