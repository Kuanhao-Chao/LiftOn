import re, sys, os, copy
from Bio.Seq import Seq
from lifton import align, lifton_class, run_liftoff, run_miniprot
from lifton.liftoff import liftoff_main

def check_liftoff_installed():
    ################################
    # Checkk if liftoff and miniprot are installed
    ################################
    liftoff_installed = run_liftoff.check_liftoff_installed()
    print("liftoff_installed : ", liftoff_installed)
    if not liftoff_installed:
        if not liftoff_installed:
            print("Liftoff is not properly installed.")
        return sys.exit(1)


def check_miniprot_installed():
    ################################
    # Checkk if liftoff and miniprot are installed
    ################################
    miniprot_installed = run_miniprot.check_miniprot_installed()
    print("miniprot_installed: ", miniprot_installed)
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
        print(">> Running Liftoff ...")
        liftoff_args = copy.deepcopy(args)
        liftoff_outdir = os.path.dirname(args.output) + "/liftoff/"
        os.makedirs(liftoff_outdir, exist_ok=True)
        liftoff_annotation = liftoff_outdir + "liftoff.gff3"
        liftoff_args.output = liftoff_annotation
        liftoff_main.run_all_liftoff_steps(liftoff_args)
        run_liftoff.run_liftoff(args)
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
        print(">> Running miniprot ...")
        miniprot_annotation = run_miniprot.run_miniprot(args, tgt_genome, ref_proteins_file)
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

    # print(segment1, segment2)

    if len(segment1) != 2 or len(segment2) != 2:
        raise ValueError("Segments must have exactly 2 endpoints")
    
    # Sort the segments by their left endpoints
    segment1, segment2 = sorted([segment1, segment2], key=lambda x: x[0])


    # Check if the right endpoint of the first segment is greater than or equal to the left endpoint of the second segment
    # print(segment1[1] >= segment2[0])

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

# def get_ID_base(id):
#     id_base = id.split("_")[0]

#     id_base.split("-")[:-1]
#     return id_base

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
    id = feature["ID"][0]
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


def update_copy(id_base, copy_num_dict):
    if id_base in copy_num_dict.keys():
        copy_num_dict[id_base] += 1
    else:
        copy_num_dict[id_base] = 0
    

def LiftOn_check_miniprot_alignment(chromosome, transcript, lifton_status, m_id_dict, m_feature_db, tree_dict, fai, ref_proteins, ref_trans_id):
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
                print("Not overlapped")
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

            if m_lifton_aln == None or m_lifton_aln.identity > lifton_status.miniprot:
                m_lifton_aln = align.parasail_align("miniprot", m_feature_db, m_entry, fai, ref_proteins, ref_trans_id, lifton_status)
                # SETTING miniprot identity score                
                lifton_status.miniprot = m_lifton_aln.identity
        
    return m_lifton_aln, has_valid_miniprot


def get_ref_liffover_features(features, ref_db):
    ref_features_dict = {}
    ref_trans_2_gene_dict = {}
    gene_info_dict = {}
    trans_info_dict = {}
    
    for feature in features:
        for gene in ref_db.db_connection.features_of_type(feature):#, limit=("CM033155.1", 0, 
            gene_id = gene.attributes["ID"][0]
            gene_info = copy.deepcopy(gene)
            lifton_gene_info = lifton_class.Lifton_GENE_info(gene_info.attributes, gene_id)
            gene_info_dict[gene_id] = lifton_gene_info
            
            if gene_id not in ref_features_dict.keys():
                ref_features_dict[gene_id] = {}

            transcripts = ref_db.db_connection.children(gene, level=1)
            for transcript in list(transcripts):
                transcript_id = transcript["ID"][0]
                transcript_info = copy.deepcopy(transcript)
                lifton_trans_info = lifton_class.Lifton_TRANS_info(transcript_info.attributes, transcript_id, gene_id)
                trans_info_dict[transcript_id] = lifton_trans_info

                ref_features_dict[gene_id][transcript_id] = False
                ref_trans_2_gene_dict[transcript_id] = gene_id

    print("ref gene count : ", len(ref_features_dict), "(", len(gene_info_dict), ")")
    print("ref trans count: ", len(trans_info_dict))
    return ref_features_dict, ref_trans_2_gene_dict, gene_info_dict, trans_info_dict


def write_lifton_status(fw_score, transcript_id, transcript, lifton_status):
    final_status = ";".join(lifton_status.status)
    fw_score.write(f"{transcript_id}\t{lifton_status.liftoff}\t{lifton_status.miniprot}\t{lifton_status.lifton}\t{lifton_status.annotation}\t{final_status}\t{transcript.seqid}:{transcript.start}-{transcript.end}\n")


def write_lifton_status(fw_score, transcript_id, transcript, lifton_status):
    final_status = ";".join(lifton_status.status)
    fw_score.write(f"{transcript_id}\t{lifton_status.liftoff}\t{lifton_status.miniprot}\t{lifton_status.lifton}\t{lifton_status.annotation}\t{final_status}\t{transcript.seqid}:{transcript.start}-{transcript.end}\n")


def segments_overlap_length(segment1, segment2):
    # Check if the segments have valid endpoints
    # print("Checking two segments overlapping.!")
    if len(segment1) != 2 or len(segment2) != 2:
        raise ValueError("Segments must have exactly 2 endpoints")
    
    # Sort the segments by their left endpoints
    segment1, segment2 = sorted([segment1, segment2], key=lambda x: x[0])
    # print("Checking miniprot overlapped length: ", segment1[1] - segment2[0] + 1)

    return segment1[1] - segment2[0] + 1

def check_ovps_ratio(mtrans_interval, ovps, overlap_ratio):
    is_overlapped = False
    for ovp in ovps:
        ovp_len = segments_overlap_length((mtrans_interval[0], mtrans_interval[1]), (ovp[0], ovp[1]))
        ref_len = ovp[1] - ovp[0] + 1
        # Overlapping does not extend the ratio of the reference
        if (ovp_len / ref_len) > overlap_ratio:
            # print("Overlapped!!: ", (ovp_len / ref_len))
            is_overlapped = True
            break
    return is_overlapped


