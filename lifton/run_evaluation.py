import copy 
from lifton import lifton_class, lifton_utils

def initialize_lifton_gene_eval(locus, ref_db, tree_dict, ref_features_dict, args, with_exons=False):
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
    if ref_gene_id is None: return None, None, None
    if not args.evaluation_liftoff_chm13:
        lifton_gene = lifton_class.Lifton_GENE(ref_gene_id, copy.deepcopy(locus), copy.deepcopy(ref_db[ref_gene_id].attributes), tree_dict, ref_features_dict, args, tmp=with_exons)
    else:
        lifton_gene = lifton_class.Lifton_GENE(ref_gene_id, copy.deepcopy(locus), copy.deepcopy(ref_db["gene-"+ref_gene_id].attributes), tree_dict, ref_features_dict, args, tmp=with_exons)
    return lifton_gene, ref_gene_id, ref_trans_id


def lifton_add_trans_exon_cds_eval(lifton_gene, locus, ref_db, tgt_db, ref_trans_id, args):
    """
        This function adds transcript, exons, and CDSs to the Lifton gene instance.

        Parameters:
        - lifton_gene: Lifton gene instance
        - locus: gffutils feature instance
        - ref_db: reference database
        - tgt_db: liftoff feature database
        - ref_trans_id: reference transcript ID

        Returns:
        lifton_trans: Lifton transcript instance
        len(cdss_list): number of CDSs
    """
    try:
        ref_db["rna-"+ref_trans_id]
    except:
        return None, 0
    if not args.evaluation_liftoff_chm13:
        lifton_trans = lifton_gene.add_transcript(ref_trans_id, copy.deepcopy(locus), copy.deepcopy(ref_db[ref_trans_id].attributes))
    else:
        try:
            ref_db["rna-"+ref_trans_id]
        except:
            return None, 0
        lifton_trans = lifton_gene.add_transcript(ref_trans_id, copy.deepcopy(locus), copy.deepcopy(ref_db["rna-"+ref_trans_id].attributes))
    exons = tgt_db.children(locus, featuretype='exon', order_by='start')
    for exon in list(exons):
        lifton_gene.add_exon(lifton_trans.entry.id, exon)
    cdss = tgt_db.children(locus, featuretype=('CDS', 'stop_codon'), order_by='start') 
    cdss_list = list(cdss)
    for cds in cdss_list:
        lifton_gene.add_cds(lifton_trans.entry.id, cds)
    return lifton_trans, len(cdss_list)


def evaluation(lifton_gene, locus, ref_db, tgt_db, tree_dict, tgt_fai, ref_proteins, ref_trans, ref_features_dict, fw_score, args, ENTRY_FEATURE=False):
    """
        This function evaluates annotation.

        Parameters:
        - lifton_gene: Lifton gene instance
        - locus: feature instance 
        - ref_db: reference database
        - tgt_db: liftoff feature database
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
    exon_children = list(tgt_db.children(locus, featuretype='exon', level=1, order_by='start'))
    if lifton_gene is None and ENTRY_FEATURE:   
        # Gene (1st) features
        lifton_gene, ref_gene_id, ref_trans_id = initialize_lifton_gene_eval(locus, ref_db, tree_dict, ref_features_dict, args, with_exons=len(exon_children)>0)
        if lifton_gene is None or lifton_gene.ref_gene_id is None: return None            
    if len(exon_children) == 0:
        parent_feature = None
        if ENTRY_FEATURE: # Gene (1st) features without direct exons 
            parent_feature = lifton_gene
        else: # Middle features without exons
            parent_feature = lifton_gene.add_feature(copy.deepcopy(locus))
        features = tgt_db.children(locus, level=1)
        for feature in list(features):
            lifton_gene = evaluation(parent_feature, feature, ref_db, tgt_db, tree_dict, tgt_fai, ref_proteins, ref_trans, ref_features_dict, fw_score, args)
    else:
        if ENTRY_FEATURE: # Gene (1st) features with direct exons 
            ref_trans_id = ref_gene_id
        else: # Transcript features with direct exons 
            ref_gene_id, ref_trans_id = lifton_utils.get_ref_ids_liftoff(ref_features_dict, lifton_gene.entry.id, locus.id)
        lifton_status = lifton_class.Lifton_Status()
        lifton_trans, cds_num = lifton_add_trans_exon_cds_eval(lifton_gene, locus, ref_db, tgt_db, ref_trans_id, args)
        if lifton_trans == None: return lifton_gene
        lifton_trans_aln, lifton_aa_aln = lifton_gene.orf_search_protein(lifton_trans.entry.id, ref_trans_id, tgt_fai, ref_proteins, ref_trans, lifton_status, eval_only=True, eval_liftoff_chm13=args.evaluation_liftoff_chm13)
        lifton_utils.print_lifton_status(lifton_trans.entry.id, locus, lifton_status, DEBUG=args.debug)
        lifton_utils.write_lifton_eval_status(fw_score, lifton_trans.entry.id, locus, lifton_status)
    return lifton_gene