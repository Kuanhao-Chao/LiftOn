from lifton import align, logger, lifton_class, lifton_utils
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
    print("args.mp_options: ", args.mp_options)
    command = [miniprot_path, "--gff-only", tgt_genome, ref_proteins_file] + args.mp_options.split(" ")
    print("miniprot: ", " ".join(command))
    try:
        fw = open(miniprot_output, "w")
        subprocess.run(command, stdout=fw)
        fw.close()
    except: 
        print("failed to run miniprot")
        sys.exit(1)
    return miniprot_output


def lifton_miniprot_with_ref_protein(m_feature, m_feature_db, ref_db, ref_gene_id, ref_trans_id, tgt_fai, 
ref_proteins, ref_trans, tree_dict, ref_features_dict, args):
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
        - args: LiftOn arguments

        Returns:
        lifton_gene: LiftOn gene instance
        lifton_transcript_id: LiftOn transcript ID
        lifton_status: LiftOn status
    """
    mtrans_id = m_feature.attributes["ID"][0]
    # Create LifOn gene instance
    m_gene_feature = copy.deepcopy(m_feature)
    m_gene_feature.featuretype = "gene"
    lifton_gene = lifton_class.Lifton_GENE(ref_gene_id, m_gene_feature, copy.deepcopy(ref_db[ref_gene_id].attributes), tree_dict, ref_features_dict, args)
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
    m_lifton_aln = align.lifton_parasail_align(Lifton_trans, m_entry, tgt_fai, ref_proteins, ref_trans_id)
    lifton_status.annotation =  "miniprot"
    lifton_status.lifton_aa = m_lifton_aln.identity
    return lifton_gene, Lifton_trans, Lifton_trans.entry.id, lifton_status


def process_miniprot(mtrans, ref_db, m_feature_db, tree_dict, tgt_fai, ref_proteins, ref_trans, ref_features_dict, fw_score, m_id_2_ref_id_trans_dict, ref_features_len_dict, ref_trans_exon_num_dict, ref_features_reverse_dict, args):
    mtrans_id = mtrans.attributes["ID"][0]
    mtrans_interval = Interval(mtrans.start, mtrans.end, mtrans_id)
    is_overlapped = lifton_utils.check_ovps_ratio(mtrans, mtrans_interval, args.overlap, tree_dict)
    lifton_gene = None
    lifton_trans = None
    if not is_overlapped:
        ref_trans_id = m_id_2_ref_id_trans_dict[mtrans_id]            
        ref_gene_id, ref_trans_id = lifton_utils.get_ref_ids_miniprot(ref_features_reverse_dict, mtrans_id, m_id_2_ref_id_trans_dict)
        if ref_gene_id is None:
            return None
        if ref_trans_id in ref_proteins.keys() and ref_trans_id in ref_trans.keys():
            if ref_trans_id != None:
                # Check if the additional copy is valid
                # 1. Remove processed pseudogenes: 1 CDS in miniprot but >1 CDS in reference
                # 2. Check the trans ratio is in the range of min_minprot and max_minprot
                if len(list(m_feature_db.children(mtrans, featuretype='CDS'))) == 1 and ref_trans_exon_num_dict[ref_trans_id] > 1: # Processed pseudogene
                    return None
                miniprot_trans_ratio = (mtrans.end - mtrans.start + 1) / ref_features_len_dict[ref_gene_id]
                if miniprot_trans_ratio > args.min_miniprot and miniprot_trans_ratio < args.max_miniprot:
                    lifton_gene, lifton_trans, transcript_id, lifton_status = lifton_miniprot_with_ref_protein(mtrans, m_feature_db, ref_db.db_connection, ref_gene_id, ref_trans_id, tgt_fai, ref_proteins, ref_trans, tree_dict, ref_features_dict, args)
                    lifton_gene.transcripts[transcript_id].entry.attributes["miniprot_annotation_ratio"] = [f"{miniprot_trans_ratio:.3f}"]
                else: # Invalid miniprot transcript
                    return None
            else: # Skip those cannot be found in reference.
                return None
        lifton_trans_aln, lifton_aa_aln = lifton_gene.orf_search_protein(lifton_trans.entry.id, ref_trans_id, tgt_fai, ref_proteins, ref_trans, lifton_status)
        lifton_utils.print_lifton_status(transcript_id, mtrans, lifton_status, DEBUG=args.debug)
        lifton_gene.add_lifton_gene_status_attrs("miniprot")
        lifton_gene.add_lifton_trans_status_attrs(transcript_id, lifton_status)
        lifton_utils.write_lifton_status(fw_score, transcript_id, mtrans, lifton_status)
    return lifton_gene