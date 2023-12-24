import re, sys, os, copy
from Bio.Seq import Seq
from lifton import align, lifton_class, run_liftoff, run_miniprot, logger
from lifton.liftoff import liftoff_main

def check_miniprot_installed():
    ################################
    # Checkk if liftoff and miniprot are installed
    ################################
    miniprot_installed = run_miniprot.check_miniprot_installed()
    # print("miniprot_installed: ", miniprot_installed)
    if not miniprot_installed:
        if not miniprot_installed:
            print("miniprot is not properly installed.")
        return sys.exit(1)

def get_truncated_protein(ref_proteins):
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
    # Run liftoff with no extra-copies
    ################################
    # Check if liftoff and miniprot results are generated
    ################################
    # liftoff_annotation = outdir + "/" + "liftoff.gff3"
    # print("liftoff_annotation  : ", liftoff_annotation)
    
    ################################
    # Execute liftoff
    ################################
    liftoff_annotation = args.liftoff
    if liftoff_annotation is None or not os.path.exists(liftoff_annotation):
        print("\n*********************")
        print("** Running Liftoff **")
        print("*********************")
        liftoff_annotation = run_liftoff.run_liftoff(outdir, args)
    return liftoff_annotation


def exec_miniprot(outdir, args, tgt_genome, ref_proteins_file):
    ################################
    # Check if liftoff and miniprot results are generated
    ################################
    # miniprot_annotation = outdir + "/" + "miniprot.gff3"
    # print("miniprot_annotation : ", miniprot_annotation)
    ################################
    # Execute miniprot
    ################################
    check_miniprot_installed()
    miniprot_annotation = args.miniprot
    if miniprot_annotation is None or not os.path.exists(miniprot_annotation):
        print("\n**********************")
        print("** Running miniprot **")
        print("**********************")
        miniprot_annotation = run_miniprot.run_miniprot(outdir, args, tgt_genome, ref_proteins_file)
    return miniprot_annotation


def get_child_types(parent_types, db):
    child_types = set()
    for parent in parent_types:
        for feature in db.db_connection.features_of_type(featuretype=parent):
            child_count = 0
            for child in db.db_connection.children(feature):
                child_count += 1
                if db.is_lowest_child(child):
                    child_types.add(child.featuretype)
            if child_count == 0:
                child_types.add(feature.featuretype)
    return child_types


def segments_overlap(segment1, segment2):
    # Check if the segments have valid endpoints
    # print("Checking two segments overlapping.!")
    if len(segment1) != 2 or len(segment2) != 2:
        raise ValueError("Segments must have exactly 2 endpoints")
    
    # Sort the segments by their left endpoints
    segment1, segment2 = sorted([segment1, segment2], key=lambda x: x[0])
    # Check if the right endpoint of the first segment is greater than or equal to the left endpoint of the second segment
    return segment1[1] >= segment2[0]


def custom_bisect_insert(sorted_list, element_to_insert):
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
    # # Regular expression pattern to match the desired substrings
    splits = id.split("_")
    try:
        int(splits[-1])
        id_base = "_".join(splits[:-1])
    except:
        id_base = id
    return id_base

def get_ID(feature):
    # print("feature: ", feature)

    # id_spec={"ID"}
    id = feature.id
    id_base = get_ID_base(id)
    return id, id_base

def get_parent_features_to_lift(feature_types_file):
    feature_types = ["gene"]
    if feature_types_file is not None:
        feature_types = []
        f = open(feature_types_file)
        for line in f.readlines():
            feature_types.append(line.rstrip())
    return feature_types


def LiftOn_check_miniprot_alignment(lifton_trans, chromosome, transcript, lifton_status, m_id_dict, m_feature_db, tree_dict, fai, ref_proteins, ref_trans_id):
    m_lifton_aln = None
    has_valid_miniprot = False

    if (ref_trans_id in m_id_dict.keys()) and (ref_trans_id in ref_proteins.keys()):
        m_ids = m_id_dict[ref_trans_id]
        for m_id in m_ids:
            ##################################################
            # Check 1: Check if the miniprot transcript is overlapping the current gene locus
            ##################################################
            m_entry = m_feature_db[m_id]
            overlap = segments_overlap((m_entry.start, m_entry.end), (transcript.start, transcript.end))
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

            ################################
            # Step 3: valid miniprot transcript exists => check if the miniprot transcript is valid
            ################################
            has_valid_miniprot = True

            ###########################
            # Add LifOn transcript instance
            ###########################
            miniprot_trans = lifton_class.Lifton_TRANS(m_id, "", "", 0, m_entry, {})
            ###########################
            # Add LiftOn exons
            ###########################
            exons = m_feature_db.children(m_entry, featuretype=('CDS', 'stop_codon'), order_by='start')
            for exon in list(exons):
                miniprot_trans.add_exon(exon)
            # logger.log(f"\tAfter adding exons\n", debug=DEBUG)
            ###########################
            # Add LiftOn CDS
            ###########################
            cdss = m_feature_db.children(m_entry, featuretype=('CDS', 'stop_codon'), order_by='start') 
            cds_num = 0
            for cds in list(cdss):
                cds_num += 1
                miniprot_trans.add_cds(cds)

            tmp_m_lifton_aln = align.parasail_align("miniprot", miniprot_trans, m_entry, fai, ref_proteins, ref_trans_id, lifton_status)
            if m_lifton_aln == None or tmp_m_lifton_aln.identity > lifton_status.miniprot:
                m_lifton_aln = tmp_m_lifton_aln
                # SETTING miniprot identity score                
                lifton_status.miniprot = m_lifton_aln.identity

    return m_lifton_aln, has_valid_miniprot


def get_ref_liffover_features(features, ref_db):
    ref_features_dict = {}
    ref_features_reverse_dict = {}
    new_gene_feature = lifton_class.Lifton_feature("Lifton-gene")
    ref_features_dict["LiftOn-gene"] = new_gene_feature

    for f_itr in features:
        for locus in ref_db.db_connection.features_of_type(f_itr):

            CDS_children = list(ref_db.db_connection.children(locus, featuretype='CDS'))

            feature = lifton_class.Lifton_feature(locus.id)
            if len(CDS_children) > 0:
                # This is the protien-coding gene
                feature.is_protein_coding = True
                
            exon_children = list(ref_db.db_connection.children(locus, featuretype='exon', level=1, order_by='start'))

            if len(exon_children) > 0:
                process_ref_liffover_features(locus, ref_db, None)
            else:
                transcripts = ref_db.db_connection.children(locus, level=1)
                for transcript in list(transcripts):
                    # print(transcript)
                    process_ref_liffover_features(transcript, ref_db, feature)
                    ref_features_reverse_dict[transcript.id] = locus.id
            ref_features_dict[locus.id] = feature
    return ref_features_dict, ref_features_reverse_dict


def process_ref_liffover_features(locus, ref_db, feature):
    if feature != None:
        feature.children.add(locus.id)


def get_ref_ids_liftoff(ref_features_dict, liftoff_gene_id, liftoff_trans_id):
    # Three cases that will be passed into this function
    #  1.  Liftoff_gene_id, None
    #  2.  None, Liftoff_trans_id
    #  3.  Liftoff_gene_id, Liftoff_trans_id
    if liftoff_gene_id is None:
        ref_trans_id = extract_ref_ids(ref_features_dict, liftoff_trans_id)
        return None, ref_trans_id 

    if liftoff_trans_id is None:
        ref_gene_id = extract_ref_ids(ref_features_dict, liftoff_gene_id)
        return ref_gene_id, None

    ref_gene_id = extract_ref_ids(ref_features_dict, liftoff_gene_id)
    if ref_gene_id is None:
        return None, None
    else:
        ref_trans_id = get_ID_base(liftoff_trans_id)
        return ref_gene_id, ref_trans_id


def extract_ref_ids(ref_features_dict, liftoff_id):
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

def segments_overlap_length(segment1, segment2):
    # Check if the segments have valid endpoints
    # print("Checking two segments overlapping.!")
    if len(segment1) != 2 or len(segment2) != 2:
        raise ValueError("Segments must have exactly 2 endpoints")
    
    # Sort the segments by their left endpoints
    segment1, segment2 = sorted([segment1, segment2], key=lambda x: x[0])
    # print("Checking miniprot overlapped length: ", segment1[1] - segment2[0] + 1)

    return segment1[1] - segment2[0] + 1

def check_ovps_ratio(mtrans, mtrans_interval, overlap_ratio, tree_dict):
    is_overlapped = False
    if mtrans.seqid not in tree_dict.keys():
        return is_overlapped
    
    ovps = tree_dict[mtrans.seqid].overlap(mtrans_interval)
    for ovp in ovps:
        ovp_len = segments_overlap_length((mtrans_interval[0], mtrans_interval[1]), (ovp[0], ovp[1]))
        ref_len = ovp[1] - ovp[0] + 1
        # Overlapping does not extend the ratio of the reference
        if (ovp_len / ref_len) > overlap_ratio:
            # print("Overlapped!!: ", (ovp_len / ref_len))
            is_overlapped = True
            break
    return is_overlapped


