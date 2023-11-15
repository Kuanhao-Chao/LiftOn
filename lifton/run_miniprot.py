from lifton import logger, lifton_class, align
import subprocess, os, sys, copy

def check_miniprot_installed():
    miniprot_path = "miniprot"
    command = [miniprot_path, "--version"]
    installed = False
    try:
        res = subprocess.run(command)
        print(res)
        installed = True
    except: 
        pass
    return installed

def run_miniprot(args, tgt_genome, ref_proteins_file):
    # print(">> run_miniprot")

    miniprot_outdir = os.path.dirname(args.output) + "/miniprot/"
    os.makedirs(miniprot_outdir, exist_ok=True)
    miniprot_output = miniprot_outdir + "miniprot.gff3"
    
    miniprot_path = "miniprot"
    command = [miniprot_path, "--gff-only", tgt_genome, ref_proteins_file]
    print("miniprot command: ", command)
    try:
        fw = open(miniprot_output, "w")
        subprocess.run(command, stdout=fw)
        fw.close()
    except: 
        print("failed to run miniprot")
        sys.exit(1)
    return miniprot_output


def lifton_miniprot_with_ref_protein(m_feature, m_feature_db, ref_gene_id, ref_trans_id, ref_gene_info_dict, ref_trans_info_dict, gene_copy_num_dict, trans_copy_num_dict, tgt_fai, ref_proteins, ref_trans, tree_dict, DEBUG):
    logger.log("\tminiprot with reference protein", debug=DEBUG)

    mtrans_id = m_feature.attributes["ID"][0]

    gene_feature = copy.deepcopy(m_feature)
    trans_feature = copy.deepcopy(m_feature)
    ###########################
    # Create LifOn gene instance
    ###########################
    lifton_gene = lifton_class.Lifton_GENE(ref_gene_id, gene_feature, ref_gene_info_dict[ref_gene_id], gene_copy_num_dict, tree_dict)
    lifton_gene.update_gene_info(m_feature.seqid, m_feature.start, m_feature.end)

    ###########################
    # Create LifOn transcript instance
    ###########################
    new_m_trans_id = lifton_gene.add_miniprot_transcript(ref_trans_id, trans_feature, ref_trans_info_dict[ref_trans_id], gene_copy_num_dict, trans_copy_num_dict)
    lifton_gene.update_trans_info(new_m_trans_id, m_feature.seqid, m_feature.start, m_feature.end)

    #######################################
    # Create exon / CDS entries
    #######################################
    cdss = m_feature_db.children(m_feature, featuretype='CDS')  # Replace 'exon' with the desired child feature type
    for cds in list(cdss):
        lifton_gene.add_exon(new_m_trans_id, cds)
        cds_copy = copy.deepcopy(cds)
        lifton_gene.add_cds(new_m_trans_id, cds_copy)

    #######################################
    # Update LiftOn status
    #######################################
    lifton_status = lifton_class.Lifton_Status()                
    m_entry = m_feature_db[mtrans_id]
    m_lifton_aln = align.parasail_align("miniprot", m_feature_db, m_entry, tgt_fai, ref_proteins, ref_trans_id, lifton_status)
    lifton_status.lifton = m_lifton_aln.identity

    if m_lifton_aln.identity == 1:
        lifton_status.annotation =  "miniprot_identical"
    elif m_lifton_aln.identity < 1:
        lifton_status.annotation =  "miniprot_truncated"
    
    lifton_trans_aln, lifton_aa_aln = lifton_gene.fix_truncated_protein(new_m_trans_id, ref_trans_id, tgt_fai, ref_proteins, ref_trans, lifton_status)

    return lifton_gene, new_m_trans_id, lifton_status



def lifton_miniprot_no_ref_protein(m_feature, m_feature_db, ref_gene_id, ref_trans_id, gene_copy_num_dict, trans_copy_num_dict, tree_dict, DEBUG):
    logger.log("\tminiprot with NO reference protein", debug=DEBUG)

    gene_feature = copy.deepcopy(m_feature)
    trans_feature = copy.deepcopy(m_feature)

    # Reason it's missing => the mRNA does not belong to gene (vdj segments) || the mRNA is not in the reference annotation
    #######################################
    # Create LiftOn gene entry
    #######################################
    gene_attrs = lifton_class.Lifton_GENE_info({}, ref_gene_id)
    lifton_gene = lifton_class.Lifton_GENE(ref_gene_id, gene_feature, gene_attrs, gene_copy_num_dict, tree_dict)
    lifton_gene.update_gene_info(m_feature.seqid, m_feature.start, m_feature.end)

    #######################################
    # Create LiftOn transcript entry
    #######################################
    trans_attrs = lifton_class.Lifton_TRANS_info({}, ref_trans_id, ref_gene_id)
    new_m_trans_id = lifton_gene.add_miniprot_transcript(ref_trans_id, trans_feature, trans_attrs, gene_copy_num_dict, trans_copy_num_dict)
    lifton_gene.update_trans_info(new_m_trans_id, m_feature.seqid, m_feature.start, m_feature.end)

    #######################################
    # Create the exon / CDS entry
    #######################################
    cdss = m_feature_db.children(m_feature, featuretype='CDS')  # Replace 'exon' with the desired child feature type
    for cds in list(cdss):
        lifton_gene.add_exon(new_m_trans_id, cds)
        cds_copy = copy.deepcopy(cds)
        lifton_gene.add_cds(new_m_trans_id, cds_copy)

    #######################################
    # Update LiftOn status
    #######################################
    lifton_status = lifton_class.Lifton_Status()                
    lifton_status.annotation = "miniprot_no_ref_protein"
    lifton_status.status = ["no_ref_protein"]
    return lifton_gene, new_m_trans_id, lifton_status