import re, sys, os, copy
from Bio.Seq import Seq
from lifton import align, lifton_class, run_liftoff, run_miniprot, logger
from lifton.liftoff import liftoff_main

def check_miniprot_installed():
    """
        This function checks if miniprot is installed.

        Parameters:
        None

        Returns:
        Exist the program if miniprot is not installed.
    """
    miniprot_installed = run_miniprot.check_miniprot_installed()
    # print("miniprot_installed: ", miniprot_installed)
    if not miniprot_installed:
        if not miniprot_installed:
            print("miniprot is not properly installed.")
        return sys.exit(1)


def get_truncated_protein(ref_proteins):
    """
        This function gets the truncated proteins.

        Parameters:
        - ref_proteins: reference proteins dictionary

        Returns:
        truncated_proteins: truncated proteins dictionary
    """
    truncated_proteins = {}
    for record in ref_proteins.keys():
        protein = ref_proteins[record]
        if not check_protein_valid(str(protein)):
            truncated_proteins[record] = protein
    # print("truncated_proteins: ", len(truncated_proteins))
    # print("good_protein: ", good_protein)
    # print("bad_protein: ", bad_protein)
    return truncated_proteins


def write_seq_2_file(outdir, ref_seqs, target):
    """
        This function writes the reference sequences to a file.

        Parameters:
        - outdir: output directory
        - ref_seqs: reference sequences dictionary
        - target: target type ('truncated_proteins', 'proteins', or 'transcripts')

        Returns:
        ref_seqs_file: reference sequences file path
    """
    if target == "truncated_proteins":
        ref_seqs_file = outdir + "/proteins_truncated.fa"
    elif target == "proteins":
        ref_seqs_file = outdir + "/proteins.fa"
    elif target == "transcripts":
        ref_seqs_file = outdir + "/transcripts.fa"
    fw = open(ref_seqs_file, 'w')
    # Iterate through the original FASTA and write the records to the new FASTA file
    for record in ref_seqs.keys():
        seq = ref_seqs[record]
        fw.write(f'>{record}\n{seq}\n')
    fw.close()
    return ref_seqs_file


def check_protein_valid(protein):
    """
        This function checks if the protein is valid.

        Parameters:
        - protein: protein sequence

        Returns:
        True if the protein is valid, False otherwise.
    """
    # The length of the protein has to be greater than 0
    if len(protein) == 0:
        return False
    # Start with M
    if protein[0] != "M":
        return False
    # End with *
    if protein[-1] != "*":
        return False
    # Only 1 * in the string
    if protein.count("*") != 1:
        return False
    return True
    

def exec_liftoff(outdir, args):
    """
        This function executes liftoff.

        Parameters:
        - outdir: output directory
        - args: arguments

        Returns:
        liftoff_annotation: liftoff annotation file path
    """
    liftoff_annotation = args.liftoff
    if liftoff_annotation is None or not os.path.exists(liftoff_annotation):
        print("\n*********************")
        print("** Running Liftoff **")
        print("*********************")
        liftoff_annotation = run_liftoff.run_liftoff(outdir, args)
    return liftoff_annotation


def exec_miniprot(outdir, args, tgt_genome, ref_proteins_file):
    """
        This function executes miniprot.

        Parameters:
        - outdir: output directory
        - args: arguments
        - tgt_genome: target genome
        - ref_proteins_file: reference protein file

        Returns:
        miniprot_annotation: miniprot annotation file path
    """
    check_miniprot_installed()
    miniprot_annotation = args.miniprot
    if miniprot_annotation is None or not os.path.exists(miniprot_annotation):
        print("\n**********************")
        print("** Running miniprot **")
        print("**********************")
        miniprot_annotation = run_miniprot.run_miniprot(outdir, args, tgt_genome, ref_proteins_file)
    return miniprot_annotation


def custom_bisect_insert(sorted_list, element_to_insert):
    """
        This function bisects the sorted list and inserts the element.

        Parameters:
        - sorted_list: sorted list
        - element_to_insert: element to insert

        Returns:
        None
    """
    low = 0
    high = len(sorted_list)

    while low < high:
        mid = (low + high) // 2
        if sorted_list[mid].entry.end < element_to_insert.entry.end:
            low = mid + 1
        else:
            high = mid
    sorted_list.insert(low, element_to_insert)


def get_ID_base(id):
    """
        This function gets the ID base by simply removing the last substring after "_".

        Parameters:
        - id: ID

        Returns:    
        ID base
    """
    # Regular expression pattern to match the desired substrings
    splits = id.split("_")
    try:
        int(splits[-1])
        id_base = "_".join(splits[:-1])
    except:
        id_base = id
    return id_base


def get_ID(feature):
    """
        This function gets the ID.

        Parameters:
        - feature: feature

        Returns:    
        the original ID and the ID base
    """
    # print("feature: ", feature)
    # id_spec={"ID"}
    id = feature.id
    id_base = get_ID_base(id)
    return id, id_base


def get_parent_features_to_lift(feature_types_file):
    """
        This function gets the parent features to lift.

        Parameters:
        - feature_types_file: feature types file (from '-f' / '--features')

        Returns:
        feature_types: feature types
    """
    feature_types = ["gene"]
    if feature_types_file is not None:
        feature_types = []
        f = open(feature_types_file)
        for line in f.readlines():
            feature_types.append(line.rstrip())
    return feature_types


def LiftOn_eval_alignment(eval_trans, locus, tgt_fai, ref_proteins, ref_trans_id, lifton_status):
    eval_aln = align.lifton_parasail_align(eval_trans, locus, tgt_fai, ref_proteins, ref_trans_id)
    if eval_aln != None:
        lifton_status.lifton_aa = eval_aln.identity
    return eval_aln


def LiftOn_liftoff_alignment(lifton_trans, locus, tgt_fai, ref_proteins, ref_trans_id, lifton_status):
    liftoff_aln = align.lifton_parasail_align(lifton_trans, locus, tgt_fai, ref_proteins, ref_trans_id)
    if liftoff_aln != None:
        lifton_status.liftoff = liftoff_aln.identity
    return liftoff_aln


def LiftOn_miniprot_alignment(chromosome, transcript, m_id_dict, m_feature_db, tree_dict, fai, ref_proteins, ref_trans_id, lifton_status):
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
            m_entry = m_feature_db[m_id]
            _, overlap = segments_overlap_length((m_entry.start, m_entry.end), (transcript.start, transcript.end))
            if not overlap or m_entry.seqid != transcript.seqid:
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
            ovps_miniprot = tree_dict[chromosome].overlap(m_entry.start, m_entry.end)
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
            miniprot_trans = lifton_class.Lifton_TRANS(m_id, "", "", 0, m_entry, {})
            exons = m_feature_db.children(m_entry, featuretype=('CDS', 'stop_codon'), order_by='start')
            for exon in list(exons):
                miniprot_trans.add_exon(exon)
            cdss = m_feature_db.children(m_entry, featuretype=('CDS', 'stop_codon'), order_by='start') 
            cds_num = 0
            for cds in list(cdss):
                cds_num += 1
                miniprot_trans.add_cds(cds)
            tmp_m_lifton_aln = align.lifton_parasail_align(miniprot_trans, m_entry, fai, ref_proteins, ref_trans_id)
            if m_lifton_aln == None or tmp_m_lifton_aln.identity > lifton_status.miniprot:
                m_lifton_aln = tmp_m_lifton_aln
                lifton_status.miniprot = m_lifton_aln.identity
    return m_lifton_aln, has_valid_miniprot


def get_ref_liffover_features(features, ref_db, intermediate_dir, args):
    """
        This function gets the reference liftover features.

        Parameters:
        - features: list of features to liftover
        - ref_db: reference database

        Returns:
        ref_features_dict: reference features dictionary (gene id -> transcript id)
        ref_features_reverse_dict: reference features reverse dictionary (transcript id -> gene id)
    """

    fw_gene = open(f"{intermediate_dir}/ref_feature.txt", "w")
    fw_trans = open(f'{intermediate_dir}/ref_transcript.txt', 'w')    
    ref_features_dict = {}
    ref_features_len_dict = {}
    ref_features_reverse_dict = {}
    ref_trans_exon_num_dict = {}
    new_gene_feature = lifton_class.Lifton_feature("Lifton-gene")
    ref_features_dict["LiftOn-gene"] = new_gene_feature
    for f_itr in features:
        for locus in ref_db.db_connection.features_of_type(f_itr):
            CDS_children = list(ref_db.db_connection.children(locus, featuretype='CDS'))
            feature = lifton_class.Lifton_feature(locus.id)
            # Write out reference gene features IDs
            # Decide if its type
            gene_type_key = ""
            if args.annotation_database.upper() == "REFSEQ":
                gene_type_key = "gene_biotype"
            elif args.annotation_database.upper() == "GENCODE" or args.annotation_database.upper() == "ENSEMBL" or args.annotation_database.upper() == "CHESS":
                gene_type_key = "gene_type"

            if gene_type_key in locus.attributes.keys():
                if locus.attributes[gene_type_key][0] == "protein_coding" and len(CDS_children) > 0:
                    feature.is_protein_coding = True
                    fw_gene.write(f"{locus.id}\tcoding\n")
                elif (locus.attributes[gene_type_key][0] == "lncRNA" or locus.attributes[gene_type_key][0] == "ncRNA"):
                    feature.is_non_coding = True
                    fw_gene.write(f"{locus.id}\tnon-coding\n")
                else:
                    fw_gene.write(f"{locus.id}\tother\n")
            else:
                fw_gene.write(f"{locus.id}\tother\n")
            exon_children = list(ref_db.db_connection.children(locus, featuretype='exon', level=1, order_by='start'))
            if len(exon_children) > 0:
                __process_ref_liffover_features(locus, ref_db, None)
            else:
                transcripts = ref_db.db_connection.children(locus, level=1)
                for transcript in list(transcripts):
                    __process_ref_liffover_features(transcript, ref_db, feature)
                    ref_features_reverse_dict[transcript.id if not args.evaluation_liftoff_chm13 else locus.id[4:]] = locus.id
                    all_CDS_in_trans = list(ref_db.db_connection.children(transcript, featuretype='CDS', order_by='start'))
                    if len(all_CDS_in_trans) > 0:
                        ref_trans_exon_num_dict[transcript.id if not args.evaluation_liftoff_chm13 else locus.id[4:]] = len(all_CDS_in_trans)
                    else:
                        ref_trans_exon_num_dict[transcript.id if not args.evaluation_liftoff_chm13 else locus.id[4:]] = 0
                    # Write out reference trans feature IDs
                    if feature.is_protein_coding and transcript.featuretype == "mRNA":
                        fw_trans.write(f"{transcript.id}\tcoding\n")
                    elif feature.is_non_coding and (transcript.featuretype == "ncRNA" or transcript.featuretype == "nc_RNA" or transcript.featuretype == "lncRNA" or transcript.featuretype == "lnc_RNA"):
                        fw_trans.write(f"{transcript.id}\tnon-coding\n")
                    else:
                        fw_trans.write(f"{transcript.id}\tother\n")
            ref_features_dict[locus.id if not args.evaluation_liftoff_chm13 else locus.id[5:]] = feature
            all_CDS_children = list(ref_db.db_connection.children(locus, featuretype='CDS', order_by='start'))
            if len(all_CDS_children) > 0:
                ref_features_len_dict[locus.id if not args.evaluation_liftoff_chm13 else locus.id[5:]] = all_CDS_children[-1].end - all_CDS_children[0].start + 1
            else:
                ref_features_len_dict[locus.id if not args.evaluation_liftoff_chm13 else locus.id[5:]] = 0
    fw_gene.close()
    fw_trans.close()
    return ref_features_dict, ref_features_len_dict, ref_features_reverse_dict, ref_trans_exon_num_dict


def __process_ref_liffover_features(locus, ref_db, feature):
    if feature != None:
        feature.children.add(locus.id)


def miniprot_id_mapping(m_feature_db):
    """
        This function creates a dictionary of miniprot id to reference id.

        Parameters:
        - m_feature_db: miniprot feature database

        Returns:
        ref_id_2_m_id_trans_dict: reference id to miniprot transcript ids dictionary
        m_id_2_ref_id_trans_dict: miniprot transcript id to reference id dictionary
    """
    ref_id_2_m_id_trans_dict = {}
    m_id_2_ref_id_trans_dict = {}
    for feature in m_feature_db.features_of_type("mRNA"):
        miniprot_id = feature["ID"][0]
        aa_trans_id = str(feature.attributes["Target"][0]).split(" ")[0]
        if aa_trans_id in ref_id_2_m_id_trans_dict.keys():
            ref_id_2_m_id_trans_dict[aa_trans_id].append(miniprot_id)
        else:
            ref_id_2_m_id_trans_dict[aa_trans_id] = [miniprot_id]
        m_id_2_ref_id_trans_dict[miniprot_id] = aa_trans_id
    return ref_id_2_m_id_trans_dict, m_id_2_ref_id_trans_dict


def get_ref_ids_liftoff(ref_features_dict, liftoff_gene_id, liftoff_trans_id):
    """
        This function gets the reference IDs from Liftoff IDs.

        Parameters:
        - ref_features_dict: reference features dictionary
        - liftoff_gene_id: Liftoff gene ID
        - liftoff_trans_id: Liftoff transcript ID

        Returns:
        # Three cases that will be passed into this function
        1.  Liftoff_gene_id, None
        2.  None, Liftoff_trans_id
        3.  Liftoff_gene_id, Liftoff_trans_id
    """
    if liftoff_gene_id is None:
        ref_trans_id = __extract_ref_ids(ref_features_dict, liftoff_trans_id)
        return None, ref_trans_id 

    if liftoff_trans_id is None:
        ref_gene_id = __extract_ref_ids(ref_features_dict, liftoff_gene_id)
        return ref_gene_id, None

    ref_gene_id = __extract_ref_ids(ref_features_dict, liftoff_gene_id)
    if ref_gene_id is None:
        return None, None
    else:
        ref_trans_id = get_ID_base(liftoff_trans_id)
        return ref_gene_id, ref_trans_id


def __extract_ref_ids(ref_features_dict, liftoff_id):
    if liftoff_id in ref_features_dict.keys():
        return liftoff_id
    else:
        ref_id = get_ID_base(liftoff_id)
        if ref_id in ref_features_dict.keys():
            return ref_id
        else:
            return None


def get_ref_ids_miniprot(ref_features_reverse_dict, miniprot_trans_id, m_id_2_ref_id_trans_dict):
    if miniprot_trans_id not in m_id_2_ref_id_trans_dict.keys():
        return None, None
    ref_trans_id = m_id_2_ref_id_trans_dict[miniprot_trans_id]
    if ref_trans_id not in ref_features_reverse_dict.keys():
        return None, ref_trans_id
    return ref_features_reverse_dict[ref_trans_id], ref_trans_id


def write_lifton_eval_status(fw_score, transcript_id, transcript, lifton_status):
    final_status = ";".join(lifton_status.status)
    fw_score.write(f"{transcript_id}\t{lifton_status.eval_dna}\t{lifton_status.eval_aa}\t{lifton_status.annotation}\t{final_status}\t{transcript.seqid}:{transcript.start}-{transcript.end}\n")


def print_lifton_status(transcript_id, transcript, lifton_status, DEBUG=False):
    final_status = ";".join(lifton_status.status)
    logger.log(f"{transcript_id}\t{lifton_status.liftoff}\t{lifton_status.miniprot}\t{lifton_status.lifton_dna}\t{lifton_status.lifton_aa}\t{lifton_status.annotation}\t{final_status}\t{transcript.seqid}:{transcript.start}-{transcript.end}\n", debug=DEBUG)


def write_lifton_status(fw_score, transcript_id, transcript, lifton_status):
    final_status = ";".join(lifton_status.status)
    fw_score.write(f"{transcript_id}\t{lifton_status.liftoff}\t{lifton_status.miniprot}\t{lifton_status.lifton_dna}\t{lifton_status.lifton_aa}\t{lifton_status.annotation}\t{final_status}\t{transcript.seqid}:{transcript.start}-{transcript.end}\n")


def write_lifton_eval_status(fw_score, transcript_id, transcript, lifton_status):
    final_status = ";".join(lifton_status.status)
    fw_score.write(f"{transcript_id}\t{lifton_status.lifton_dna}\t{lifton_status.lifton_aa}\t{final_status}\t{transcript.seqid}:{transcript.start}-{transcript.end}\n")


def write_lifton_chains(fw_chain, transcript_id, chains):
    chain_ls = ";".join(chains)
    fw_chain.write(f"{transcript_id}\t{chain_ls}\n")


def segments_overlap_length(segment1, segment2):
    """
        This function gets the length of the overlapping segments.

        Parameters:
        - segment1: segment 1 in tuple (start, end)
        - segment2: segment 2 in tuple (start, end)

        Returns:
        The length of the overlapping segments.
    """
    if len(segment1) != 2 or len(segment2) != 2:
        raise ValueError("Segments must have exactly 2 endpoints")
    # Sort the segments by their left endpoints
    segment1, segment2 = sorted([segment1, segment2], key=lambda x: x[0])
    ovp_len = segment1[1] - segment2[0] + 1
    ovp = False
    if ovp_len > 0: ovp = True
    return ovp_len, ovp


def check_ovps_ratio(mtrans, mtrans_interval, overlap_ratio, tree_dict):
    """
        This function checks the overlap ratio.

        Parameters:
        - mtrans: miniprot transcript gffutils feature
        - mtrans_interval: miniprot transcript interval
        - overlap_ratio: overlap ratio
        - tree_dict: tree dictionary

        Returns:
        True if the overlap ratio is greater than the threshold, False otherwise.
    """
    is_overlapped = False
    if mtrans.seqid not in tree_dict.keys():
        return False
    ovps = tree_dict[mtrans.seqid].overlap(mtrans_interval)
    for ovp in ovps:
        ovp_len, _ = segments_overlap_length((mtrans_interval[0], mtrans_interval[1]), (ovp[0], ovp[1]))
        ref_len = ovp[1] - ovp[0] + 1
        target_len = mtrans_interval[1] - mtrans_interval[0] + 1
        # Overlapping does not extend the ratio of the reference
        if (ovp_len / min(ref_len, target_len)) > overlap_ratio:
            is_overlapped = True
            break
    return is_overlapped