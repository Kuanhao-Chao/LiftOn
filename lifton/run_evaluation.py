import copy 
from lifton import lifton_class, lifton_utils, run_liftoff

def evaluation_with_protein(locus, lifton_trans, tgt_fai, ref_trans_id, ref_proteins, lifton_status):
    """
        This function process evaluated gene locus with protein.

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
    eval_aln = lifton_utils.LiftOn_eval_alignment(lifton_trans, locus, tgt_fai, ref_proteins, ref_trans_id, lifton_status)
    # miniprot alignment
    if eval_aln is None:
        lifton_status.annotation = "no_ref_protein"


def evaluation(lifton_gene, locus, ref_db, l_feature_db, tree_dict, tgt_fai, ref_proteins, ref_trans, ref_features_dict, fw_score, args, ENTRY_FEATURE=False):
    """
        This function evaluates annotation.

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
        lifton_gene, ref_gene_id, ref_trans_id = run_liftoff.initialize_lifton_gene(locus, ref_db, tree_dict, ref_features_dict, args, with_exons=len(exon_children)>0)
        if lifton_gene.ref_gene_id is None: return None            
    if len(exon_children) == 0:
        parent_feature = None
        if ENTRY_FEATURE: # Gene (1st) features without direct exons 
            parent_feature = lifton_gene
        else: # Middle features without exons
            parent_feature = lifton_gene.add_feature(copy.deepcopy(locus))
        features = l_feature_db.children(locus, level=1)
        for feature in list(features):
            lifton_gene = evaluation(parent_feature, feature, ref_db, l_feature_db, tree_dict, tgt_fai, ref_proteins, ref_trans, ref_features_dict, fw_score, args)
    else:
        if ENTRY_FEATURE: # Gene (1st) features with direct exons 
            ref_trans_id = ref_gene_id
        else: # Transcript features with direct exons 
            ref_gene_id, ref_trans_id = lifton_utils.get_ref_ids_liftoff(ref_features_dict, lifton_gene.entry.id, locus.id)
        lifton_status = lifton_class.Lifton_Status()
        lifton_status.annotation = "Evaluation"
        lifton_trans, cds_num = run_liftoff.lifton_add_trans_exon_cds(lifton_gene, locus, ref_db, l_feature_db, ref_trans_id)
        if cds_num > 0:
            evaluation_with_protein(locus, lifton_trans, tgt_fai, ref_trans_id, ref_proteins, lifton_status)
        lifton_trans_aln, lifton_aa_aln = lifton_gene.orf_search_protein(lifton_trans.entry.id, ref_trans_id, tgt_fai, ref_proteins, ref_trans, lifton_status, eval_only=True)
        lifton_utils.print_lifton_status(lifton_trans.entry.id, locus, lifton_status, DEBUG=args.debug)
    return lifton_gene