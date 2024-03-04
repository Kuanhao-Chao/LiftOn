from lifton import logger, lifton_class, align
import subprocess, os, sys, copy

def check_miniprot_installed():
    """
        This function checks if miniprot is installed.

        Parameters:
        None

        Returns:
        True if miniprot is installed, False otherwise.
    """
    miniprot_path = "miniprot"
    command = [miniprot_path, "--version"]
    installed = False
    try:
        res = subprocess.run(command)
        installed = True
    except: 
        pass
    return installed


def run_miniprot(outdir, args, tgt_genome, ref_proteins_file):
    """
        This function runs miniprot.

        Parameters:
        - outdir: output directory
        - args: arguments
        - tgt_genome: target genome
        - ref_proteins_file: reference protein file

        Returns:
        miniprot_output: miniprot output file
    """
    miniprot_outdir = outdir + "miniprot/"
    os.makedirs(miniprot_outdir, exist_ok=True)
    miniprot_output = miniprot_outdir + "miniprot.gff3"
    miniprot_path = "miniprot"
    command = [miniprot_path, "--gff-only", tgt_genome, ref_proteins_file]
    try:
        fw = open(miniprot_output, "w")
        subprocess.run(command, stdout=fw)
        fw.close()
    except: 
        print("failed to run miniprot")
        sys.exit(1)
    return miniprot_output


def lifton_miniprot_with_ref_protein(m_feature, m_feature_db, ref_db, ref_gene_id, ref_trans_id, tgt_fai, 
ref_proteins, ref_trans, tree_dict, ref_features_dict, DEBUG):
    """
        This function create a miniprot gene entry with reference protein.

        Parameters:
        - m_feature: miniprot feature
        - m_feature_db: miniprot feature database
        - ref_db: reference database
        - ref_gene_id: reference gene ID
        - ref_trans_id: reference transcript ID
        - tgt_fai: target fasta index
        - ref_proteins: reference protein dictionary
        - ref_trans: reference transcript dictionary
        - tree_dict: tree dictionary
        - ref_features_dict: reference features dictionary
        - DEBUG: debug mode

        Returns:
        lifton_gene: LiftOn gene instance
        lifton_transcript_id: LiftOn transcript ID
        lifton_status: LiftOn status
    """
    # logger.log("\tminiprot with reference protein", debug=DEBUG)
    mtrans_id = m_feature.attributes["ID"][0]
    # Create LifOn gene instance
    m_gene_feature = copy.deepcopy(m_feature)
    m_gene_feature.featuretype = "gene"
    if ref_gene_id is None:
        # This is a place holder gene
        lifton_gene = lifton_class.Lifton_GENE(ref_gene_id, m_gene_feature, {}, tree_dict, ref_features_dict, miniprot_holder=True)
    else:
        lifton_gene = lifton_class.Lifton_GENE(ref_gene_id, m_gene_feature, copy.deepcopy(ref_db[ref_gene_id].attributes), tree_dict, ref_features_dict)
    lifton_gene.update_gene_info(m_feature.seqid, m_feature.start, m_feature.end)
    # Create LifOn transcript instance
    Lifton_trans = lifton_gene.add_miniprot_transcript(ref_trans_id, copy.deepcopy(m_feature), ref_db[ref_trans_id].attributes, ref_features_dict)
    lifton_gene.update_trans_info(Lifton_trans.entry.id, m_feature.seqid, m_feature.start, m_feature.end)
    # Create exon / CDS entries
    cdss = m_feature_db.children(m_feature, featuretype='CDS')  # Replace 'exon' with the desired child feature type
    for cds in list(cdss):
        lifton_gene.add_exon(Lifton_trans.entry.id, cds)
        cds_copy = copy.deepcopy(cds)
        lifton_gene.add_cds(Lifton_trans.entry.id, cds_copy)
    # Update LiftOn status
    lifton_status = lifton_class.Lifton_Status()                
    m_entry = m_feature_db[mtrans_id]
    m_lifton_aln = align.parasail_align("miniprot", Lifton_trans, m_entry, tgt_fai, ref_proteins, ref_trans_id, lifton_status)
    lifton_status.lifton_aa = m_lifton_aln.identity
    if m_lifton_aln.identity == 1:
        lifton_status.annotation =  "miniprot_identical"
    elif m_lifton_aln.identity < 1:
        lifton_status.annotation =  "miniprot_truncated"
    lifton_trans_aln, lifton_aa_aln = lifton_gene.fix_truncated_protein(Lifton_trans.entry.id, ref_trans_id, tgt_fai, ref_proteins, ref_trans, lifton_status)
    return lifton_gene, Lifton_trans.entry.id, lifton_status


def lifton_miniprot_no_ref_protein(m_feature, m_feature_db, ref_gene_id, ref_trans_id, ref_features_dict, tree_dict, DEBUG):
    """
        This function creates a miniprot gene entry with no reference protein.

        Parameters:
        - m_feature: miniprot feature
        - m_feature_db: miniprot feature database
        - ref_gene_id: reference gene ID
        - ref_trans _id: reference transcript ID
        - ref_features_dict: reference features dictionary
        - tree_dict: tree dictionary
        - DEBUG: debug mode

        Returns:
        lifton_gene: LiftOn gene instance
        lifton_transcript_id: LiftOn transcript ID
        lifton_status: LiftOn status
    """ 
    logger.log("\tminiprot with NO reference protein", debug=DEBUG)
    # Reason it's missing => the mRNA does not belong to gene (vdj segments) || the mRNA is not in the reference annotation
    # Create LiftOn gene entry
    m_gene_feature = copy.deepcopy(m_feature)
    m_gene_feature.featuretype = "gene"
    lifton_gene = lifton_class.Lifton_GENE(ref_gene_id, m_gene_feature, {}, tree_dict, ref_features_dict)
    lifton_gene.update_gene_info(m_feature.seqid, m_feature.start, m_feature.end)
    # Create LiftOn transcript entry
    Lifton_trans = lifton_gene.add_miniprot_transcript(ref_trans_id, copy.deepcopy(m_feature), {}, ref_features_dict)
    lifton_gene.update_trans_info(Lifton_trans.entry.id, m_feature.seqid, m_feature.start, m_feature.end)
    # Create the exon / CDS entry
    cdss = m_feature_db.children(m_feature, featuretype='CDS')  # Replace 'exon' with the desired child feature type
    for cds in list(cdss):
        lifton_gene.add_exon(Lifton_trans.entry.id, cds)
        cds_copy = copy.deepcopy(cds)
        lifton_gene.add_cds(Lifton_trans.entry.id, cds_copy)
    # Update LiftOn status
    lifton_status = lifton_class.Lifton_Status()                
    lifton_status.annotation = "miniprot_no_ref_protein"
    lifton_status.status = ["no_ref_protein"]
    return lifton_gene, Lifton_trans.entry.id, lifton_status