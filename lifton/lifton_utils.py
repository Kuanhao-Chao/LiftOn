import re, sys, os
from Bio.Seq import Seq
from lifton import align, lifton_class, run_liftoff, run_miniprot

def check_software_installed():
    ################################
    # Checkk if liftoff and miniprot are installed
    ################################
    # print(run_liftoff.run_all_liftoff_steps)
    liftoff_installed = run_liftoff.check_liftoff_installed()
    miniprot_installed = run_miniprot.check_miniprot_installed()
    print("liftoff_installed : ", liftoff_installed)
    print("miniprot_installed: ", miniprot_installed)
    if not liftoff_installed or not miniprot_installed:
        if not liftoff_installed:
            print("Liftoff is not properly installed.")
        if not miniprot_installed:
            print("Miniprot is not properly installed.")
        return sys.exit(1)

def exec_liftoff(outdir):
    ################################
    # Check if liftoff and miniprot results are generated
    ################################
    liftoff_annotation = outdir + "/" + "liftoff.gff3"
    print("liftoff_annotation  : ", liftoff_annotation)
    ################################
    # Execute liftoff and miniprot
    ################################
    if not os.path.exists(liftoff_annotation):
        run_liftoff.run_liftoff()


def exec_miniprot(outdir):
    ################################
    # Check if liftoff and miniprot results are generated
    ################################
    miniprot_annotation = outdir + "/" + "miniprot.gff3"
    print("miniprot_annotation : ", miniprot_annotation)
    ################################
    # Execute liftoff and miniprot
    ################################
    if not os.path.exists(miniprot_annotation):
        run_miniprot.run_miniprot()


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


def get_feature_types(feature_arg):
    feature_types = ['gene']
    if feature_arg is not None:
        with open(feature_arg) as fa:
            for line in fa:
                feature_types.append(line.strip())
    return feature_types


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

def get_parent_features_to_lift(feature_types_file):
    feature_types = ["gene"]
    if feature_types_file is not None:
        f = open(feature_types_file)
        for line in f.readlines():
            feature_types.append(line.rstrip())
    return feature_types


def update_copy(id_base, copy_num_dict):
    if id_base in copy_num_dict.keys():
        copy_num_dict[id_base] += 1
    else:
        copy_num_dict[id_base] = 0
    

def LiftOn_check_miniprot_alignment(chromosome, transcript, lifton_status, m_id_dict, m_feature_db, tree_dict, fai, fai_protein, transcript_id):
    m_lifton_aln = None
    has_valid_miniprot = False
    if (transcript_id in m_id_dict.keys()):
        #############################################
        # Step 3.6.1.1: Liftoff annotation is not perfect & miniprot annotation exists => Fix by protein information
        #############################################
        m_ids = m_id_dict[transcript_id]

        for m_id in m_ids:

            ##################################################
            # Check 1: Check if the miniprot transcript is overlapping the current gene locus
            ##################################################
            m_entry = m_feature_db[m_id]
            overlap = segments_overlap((m_entry.start, m_entry.end), (transcript.start, transcript.end))
            if not overlap or m_entry.seqid != transcript.seqid:
                print("Not overlapping")
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
                # print("\tovp_liftoff: ", ovp_liftoff)
            # print("liftoff_set : ", liftoff_set)
            
            for ovp_miniprot in ovps_miniprot:
                if ovp_miniprot[2] not in liftoff_set:
                    # Miniprot overlap to more genes
                    miniprot_cross_gene_loci = True
                    break
            if miniprot_cross_gene_loci:
                continue

            ################################
            # Step 3.6.2: Protein sequences are in both Liftoff and miniprot & overlap
            #   Fix & update CDS list
            ################################
            ################################
            # Step 3.6.3: miniprot transcript alignment
            ################################
            has_valid_miniprot = True

            if m_lifton_aln == None or m_lifton_aln.identity > lifton_status.miniprot:
                # # Writing out truncated miniprot annotation
                # m_lifton_aln.write_alignment(outdir, "miniprot", m_id)
                # SETTING miniprot identity score
                
                m_lifton_aln = align.parasail_align("miniprot", m_feature_db, m_entry, fai, fai_protein, transcript_id)
                lifton_status.miniprot = m_lifton_aln.identity
        
    return m_lifton_aln, has_valid_miniprot
            