from lifton import logger, lifton_class, align, lifton_utils

import subprocess, os, sys, copy
from intervaltree import Interval, IntervalTree

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


















def lifton_miniprot_with_ref_protein(chromosome, transcript, lifton_status, m_id_dict, m_feature_db, tree_dict, fai, ref_proteins, ref_trans_id):
    """
        This function checks the miniprot alignment.

        Parameters:
        - chromosome: chromosome
        - transcript: transcript gffutils feature
        - lifton_status: Lifton_Status instance
        - m_id_dict: miniprot id dictionary
        - m_feature_db: miniprot feature database
        - tree_dict: tree dictionary
        - fai: reference fasta index
        - ref_proteins: reference proteins dictionary
        - ref_trans_id: reference transcript ID

        Returns:
        m_lifton_aln: miniprot lifton alignment
        has_valid_miniprot: True if the miniprot transcript is valid, False otherwise.
    """
    m_lifton_aln = None
    has_valid_miniprot = False
    if (ref_trans_id in m_id_dict.keys()) and (ref_trans_id in ref_proteins.keys()):
        m_ids = m_id_dict[ref_trans_id]
        for m_id in m_ids:
            ##################################################
            # Check 1: Check if the miniprot transcript is overlapping the current gene locus
            ##################################################
            mtrans = m_feature_db[m_id]
            mtrans_id = mtrans.attributes["ID"][0]
            # overlap = segments_overlap((m_entry.start, m_entry.end), (transcript.start, transcript.end))

            mtrans_interval = Interval(mtrans.start, mtrans.end, mtrans_id)
            is_overlapped = lifton_utils.check_ovps_ratio(mtrans, mtrans_interval, 0.7, tree_dict)

            if not is_overlapped or mtrans.seqid != transcript.seqid:
                # "Not overlapped"
                continue
            ##################################################
            # Check 2: reference overlapping status
            #   1. Check it the transcript overlapping with the next gene
            # Check the miniprot protein overlapping status
            # The case I should not process the transcript 
            #   1. The Liftoff does not overlap with other gene
            #   2. The miniprot protein overlap the other gene
            ##################################################
            ovps_liftoff = tree_dict[chromosome].overlap(transcript.start, transcript.end)
            ovps_miniprot = tree_dict[chromosome].overlap(mtrans.start, mtrans.end)
            miniprot_cross_gene_loci = False
            liftoff_set = set()
            for ovp_liftoff in ovps_liftoff:
                liftoff_set.add(ovp_liftoff[2])
            for ovp_miniprot in ovps_miniprot:
                if ovp_miniprot[2] not in liftoff_set:
                    # Miniprot overlap to more genes
                    miniprot_cross_gene_loci = True
                    break
            if miniprot_cross_gene_loci:
                continue
            # Valid miniprot transcript exists => check if the miniprot transcript is valid
            has_valid_miniprot = True

            
            miniprot_trans = lifton_class.Lifton_TRANS(m_id, "", "", 0, mtrans, {})
            exons = m_feature_db.children(mtrans, featuretype=('CDS', 'stop_codon'), order_by='start')
            for exon in list(exons):
                miniprot_trans.add_exon(exon)
            cdss = m_feature_db.children(mtrans, featuretype=('CDS', 'stop_codon'), order_by='start') 
            cds_num = 0
            for cds in list(cdss):
                cds_num += 1
                miniprot_trans.add_cds(cds)
            tmp_m_lifton_aln = align.lifton_parasail_align("miniprot", miniprot_trans, m_entry, fai, ref_proteins, ref_trans_id)
            if m_lifton_aln == None or tmp_m_lifton_aln.identity > lifton_status.miniprot:
                m_lifton_aln = tmp_m_lifton_aln
                lifton_status.miniprot = m_lifton_aln.identity
    return m_lifton_aln, has_valid_miniprot


        # if ref_trans_id in ref_id_2_m_id_trans_dict.keys():
        #     for mtrans_id in ref_id_2_m_id_trans_dict[ref_trans_id]:
        #         mtrans = m_feature_db[mtrans_id]
        #         mtrans_interval = Interval(mtrans.start, mtrans.end, mtrans_id)
        #         is_overlapped = lifton_utils.check_ovps_ratio(mtrans, mtrans_interval, 0.7, tree_dict)
        #         if is_overlapped:
        #             lifton_gene, transcript_id, lifton_status = run_miniprot.lifton_miniprot_with_ref_protein(mtrans, m_feature_db, ref_db, ref_gene_id, ref_trans_id, tgt_fai, ref_proteins, ref_trans, tree_dict, ref_features_dict, DEBUG)
        #             lifton_trans = lifton_gene.transcripts[transcript_id]
        #         else:
        #             lifton_status.annotation = "Liftoff_no_protein"
        #             lifton_status.status = ["no_protein"]




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
    m_lifton_aln = align.lifton_parasail_align("miniprot", Lifton_trans, m_entry, tgt_fai, ref_proteins, ref_trans_id)
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