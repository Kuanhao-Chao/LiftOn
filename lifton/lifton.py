from lifton import extract_features
import argparse
from pyfaidx import Fasta, Faidx
from Bio.Seq import Seq
import parasail

def get_id_fraction(reference, target):
    matches = 0
    for i, letter in enumerate(reference):
        if letter == target[i]:
            matches += 1
    return matches, max(len(reference), len(target))

def segments_overlap(segment1, segment2):
    # Check if the segments have valid endpoints
    if len(segment1) != 2 or len(segment2) != 2:
        raise ValueError("Segments must have exactly 2 endpoints")
    
    # Sort the segments by their left endpoints
    segment1, segment2 = sorted([segment1, segment2], key=lambda x: x[0])

    # Check if the right endpoint of the first segment is greater than or equal to the left endpoint of the second segment
    return segment1[1] >= segment2[0]


def parasail_align(tool, db, db_entry, fai, fai_protein, aa_trans_id):
    # Get the children of the entry
    # cds_children = [child for child in db.children(db_entry, featuretype='CDS')]
    # print("Before cds_children: ", cds_children)
    cds_children = []
    for child in db.children(db_entry, featuretype='CDS'):
        if len(cds_children) == 0:
            cds_children.append(child)
            continue
        idx_insert = 0
        for idx_c in range(len(cds_children)):
            itr_c = cds_children[idx_c]
            if child.start > itr_c.end:
                idx_insert += 1
        
        cds_children.insert(idx_insert, child)

    # print("After cds_children: ", cds_children)
    # Iterate through the children and print their attributes
    trans_seq = ""
    for cds_idx, cds in enumerate(cds_children):

        # Include the stop coding for the last CDS(+) / first CDS(-) for miniprot 
        if tool == "miniprot" and cds_idx == 0 and cds.strand == '-':
            cds.start = cds.start -3
        if tool == "miniprot" and cds_idx == len(cds_children)-1 and cds.strand == '+':
            cds.end = cds.end + 3

        p_seq = cds.sequence(fai)
        p_seq = Seq(p_seq)

        # Chaining the CDS features
        if cds.strand == '-':
            # if tool == "liftoff":
            trans_seq = p_seq + trans_seq
            # elif tool == "miniprot":
            #     trans_seq = trans_seq + p_seq
        elif cds.strand == '+':
            trans_seq = trans_seq + p_seq
        # print('>' + cds.id + '\n' + p_seq)
    # print(cds.strand+' trans_seq: ' + trans_seq)
    
    protein_seq = str(trans_seq.translate())
    ref_protein_seq = str(fai_protein[aa_trans_id])

    matrix = parasail.Matrix("blosum62")
    gap_open = 11
    gap_extend = 1

    reference_seq = str(ref_protein_seq) + "*"
    extracted_seq = str(protein_seq)
    # print("db_entry: ", db_entry)
    # print("\treference_seq: ", reference_seq)
    # print("\textracted_seq: ", extracted_seq)
    
    # (Query, Reference)
    extracted_parasail_res = parasail.nw_trace_scan_sat(extracted_seq, reference_seq, gap_open, gap_extend, matrix)

    # Extract the alignment information
    alignment_score = extracted_parasail_res.score
    alignment_query = extracted_parasail_res.traceback.query
    alignment_comp = extracted_parasail_res.traceback.comp
    alignment_ref = extracted_parasail_res.traceback.ref

    cigar = extracted_parasail_res.cigar
    # alignment_start_query = extracted_parasail_res.traceback.query_begin
    # alignment_end_query = extracted_parasail_res.traceback.query_end
    # alignment_start_comp = extracted_parasail_res.traceback.comp_begin
    # alignment_end_comp = extracted_parasail_res.traceback.comp_end

    # print("extracted_seq: ", extracted_seq)
    # # Print the additional alignment results
    # print("\t>> alignment_score: ", alignment_score)
    # print("\t>> alignment_query: ", alignment_query)
    # print("\t>> alignment_comp: ", alignment_comp)
    # # print("\t>> cigar.seq   : ", cigar.seq)
    # # # use decode attribute to return a decoded cigar string
    # # print("\t>> cigar.decode: ", cigar.decode)
    # print("\t alignment_ref : ", alignment_ref)
    
    extracted_matches, extracted_length = get_id_fraction(extracted_parasail_res.traceback.ref, extracted_parasail_res.traceback.query)
    extracted_identity = extracted_matches/extracted_length

    return extracted_identity, cds_children, alignment_query, alignment_comp, alignment_ref





def grouping_m_l_children():
    pass


def fix_transcript_annotation(m_children, m_aln_query, m_aln_comp, m_aln_ref, l_children, l_aln_query, l_aln_comp, l_aln_ref):
    print("number of children, m: ", len(m_children))
    print("number of children, l: ", len(l_children))

    m_len_sum = 0
    


    m_c_idx = 0
    l_c_idx = 0
    while m_c_idx != (len(m_children)-1) or l_c_idx != (len(l_children)-1):
        m_c = m_children[m_c_idx]
        l_c = l_children[l_c_idx]

        if m_c.end > l_c.end:
            l_c_idx += 1
        elif m_c.end < l_c.end:
            m_c_idx += 1
        elif m_c.end == l_c.end:
            l_c_idx += 1
            m_c_idx += 1
            print("\tGroup!")
        print("l_c_idx: ", l_c_idx)
        print("m_c_idx: ", m_c_idx)
    print("\tGroup!")

    # for child in m_children:
    #     child_len = child.end - child.start + 1
    #     m_len_sum += child_len
    #     print("start: ", child.start, ";  start: ", child.end, ";  len: ", child_len, ";  len_sum: ", m_len_sum, ";  aa_len: ", m_len_sum/3)
    # print("m_aln_query: ", len(m_aln_query))
    
    # print("number of children: ", len(l_children))
    # l_len_sum = 0
    # for child in l_children:
    #     child_len = child.end - child.start + 1
    #     l_len_sum += child_len
    #     print("start: ", child.start, ";  start: ", child.end, ";  len: ", child_len, ";  len_sum: ", l_len_sum, ";  aa_len: ", l_len_sum/3)
    # print("l_aln_query: ", len(l_aln_query))







def parse_args(arglist):
    print("arglist: ", arglist)
    parser = argparse.ArgumentParser(description='Lift features from one genome assembly to another')
    parser.add_argument('target', help='target fasta genome to lift genes to')

    # refrgrp = parser.add_argument_group('Required input (annotation)')
    # mxgrp = refrgrp.add_mutually_exclusive_group(required=True)
    # mxgrp.add_argument(
    #     '-g', metavar='GFF', help='annotation file to lift over in GFF or GTF format'
    # )
    # mxgrp.add_argument(
    #     '-db', metavar='DB', help='name of feature database; if not specified, the -g '
    #                               'argument must be provided and a database will be built automatically'
    # )

    outgrp = parser.add_argument_group('Output')
    outgrp.add_argument(
        '-o', default='stdout', metavar='FILE',
        help='write output to FILE in same format as input; by default, output is written to terminal (stdout)'
    )

    outgrp.add_argument(
        '-dir', default='intermediate_files', metavar='DIR',
        help='name of directory to save intermediate fasta and SAM files; default is "intermediate_files"',
    )

    parser.add_argument('-V', '--version', help='show program version', action='version', version='v1.6.3')
    parser.add_argument(
        '-p', default=1, type=int, metavar='P', help='use p parallel processes to accelerate alignment; by default p=1'
    )

    parser.add_argument(
        '-proteins', metavar='fasta', required=True,
        help='the reference protein sequences.'
    )

    liftoffrefrgrp = parser.add_argument_group('Required input (Liftoff annotation)')

    liftoffgrp = liftoffrefrgrp.add_mutually_exclusive_group(required=True)
    liftoffgrp.add_argument(
        '-liftoff', metavar='gff',
        help='the annotation generated by Liftoff'
    )
    liftoffgrp.add_argument(
        '-liftoffdb', metavar='gff-DB',
        help='name of Liftoff database; if not specified, the -liftoff '
                                  'argument must be provided and a database will be built automatically'
    )

    miniprotrefrgrp = parser.add_argument_group('Required input (miniprot annotation)')

    miniprotgrp = miniprotrefrgrp.add_mutually_exclusive_group(required=True)
    miniprotgrp.add_argument(
        '-miniprot', metavar='gff',
        help='the annotation generated by miniprot'
    )
    miniprotgrp.add_argument(
        '-miniprotdb', metavar='gff-DB',
        help='name of miniprot database; if not specified, the -miniprot '
                                  'argument must be provided and a database will be built automatically'
    )

    parser._positionals.title = 'Required input (sequences)'
    parser._optionals.title = 'Miscellaneous settings'

    parser._action_groups = [parser._positionals, liftoffrefrgrp, miniprotrefrgrp, outgrp, parser._optionals]
    args = parser.parse_args(arglist)
    return args



# def extract_features_to_lift(ref_chroms, liftover_type, parents_to_lift, args):
#     print("extracting features")
#     if os.path.exists(args.dir) is False:
#         os.mkdir(args.dir)
#     feature_db = create_feature_db_connections(args)
#     feature_hierarchy, parent_order = seperate_parents_and_children(feature_db, parents_to_lift)
#     get_gene_sequences(feature_hierarchy.parents, ref_chroms, args, liftover_type)
#     return feature_hierarchy, feature_db, parent_order


# def lift_original_annotation(ref_chroms, target_chroms, lifted_features_list, args, unmapped_features, parents_to_lift):
#     liftover_type = "chrm_by_chrm"
#     if target_chroms[0] == args.target and args.exclude_partial == False:
#         min_cov, min_seqid = 0.05, 0.05
#     else:
#         min_cov, min_seqid = args.a, args.s

#     feature_hierarchy, feature_db, ref_parent_order = extract_features.extract_features_to_lift(ref_chroms,
#                                                                                                 liftover_type,
#                                                                                                 parents_to_lift, args)

def run_all_liftoff_steps(args):
    print(">> run_all_lifton_steps")
    lifted_feature_list = {}
    unmapped_features = []

    liftover_type = "chrm_by_chrm"
    ref_chroms = []
    
    fai = Fasta(args.target)
    print("fai: ", fai.keys())

    fai_protein = Fasta(args.proteins)
    print("fai: ", fai_protein["rna-NM_001370185.1-2"])

    l_feature_db, m_feature_db = extract_features.extract_features_to_fix(ref_chroms, liftover_type, args)
    print("l_feature_db: ", l_feature_db)
    


    for feature in m_feature_db.features_of_type("mRNA"):
        # Print all attributes and their values for the feature
        # print(feature)

        aa_trans_id = str(feature.attributes["Target"][0]).split(" ")[0]
        # miniprot_identity = float(feature.attributes["Identity"][0])

        miniprot_trans_id = feature.attributes["ID"][0]
        m_entry = m_feature_db[miniprot_trans_id]

        # print(m_entry)
        miniprot_identity = 0.0
        miniprot_identity, m_children, m_aln_query, m_aln_comp, m_aln_ref = parasail_align("miniprot", m_feature_db, m_entry, fai, fai_protein, aa_trans_id)
        
        # for attr_name, attr_value in feature.attributes.items():
        #     print(f"{attr_name}: {attr_value}")

        liftoff_identity = 0.0
        try:
            l_entry = l_feature_db[aa_trans_id]
            liftoff_identity, l_children, l_aln_query, l_aln_comp, l_aln_ref = parasail_align("liftoff", l_feature_db, l_entry, fai, fai_protein, aa_trans_id)
        except:
            print("An exception occurred")


        if miniprot_identity > liftoff_identity and liftoff_identity > 0:

            overlap = segments_overlap((m_entry.start, m_entry.end), (l_entry.start, l_entry.end))
            if (overlap and m_entry.seqid == l_entry.seqid):

                fix_transcript_annotation(m_children, m_aln_query, m_aln_comp, m_aln_ref, l_children, l_aln_query, l_aln_comp, l_aln_ref)


                print(aa_trans_id)
                print(m_entry)
                print(l_entry)
                print("miniprot_identity: ", miniprot_identity, "; number of children: ", len(m_children))
                print("liftoff_identity: ", liftoff_identity, "; number of children: ", len(l_children))
                print("\n\n")
        elif liftoff_identity == 0:
            pass


    #     try:
    #         m_entry = m_feature_db[aa_trans_id]
    #         miniprot_identity = parasail_align(m_feature_db, m_entry, fai, fai_protein, aa_trans_id)
    #     except:
    #         print("An exception occurred")


        # if aa_trans_id in l_feature_db:
        #     print(aa_trans_id)
        # else:
        #     print(f"Feature not found for {aa_trans_id}")

        # if aa_trans_id in l_feature_db:
        #     print(aa_trans_id) 
            # Example 1: Extract sequence by feature ID
            # l_entry = l_feature_db[aa_trans_id]
            # print(l_entry)

    # chrom_seq = reference_fasta_idx[current_chrom][:].seq

    # if liftover_type == "unplaced":
    #     open(args.dir + "/unplaced_genes.fa", 'w')
    # for chrom in ref_chroms:
    #     fasta_out = get_fasta_out(chrom, args.reference, liftover_type, args.dir)
    #     sorted_parents = sorted(list(parent_dict.values()), key=lambda x: x.seqid)

    #     if len(sorted_parents) == 0:
    #         sys.exit(
    #             "GFF does not contain any gene features. Use -f to provide a list of other feature types to lift over.")
    #     write_gene_sequences_to_file(chrom, args.reference, fai, sorted_parents, fasta_out, args)
    #     fasta_out.close()


    # l_feature_hierarchy, l_feature_db, l_ref_parent_order, m_feature_hierarchy, m_feature_db, m_ref_parent_order = extract_features.extract_features_to_fix(ref_chroms, liftover_type, args)

        # self.parents = parents
        # self.intermediates = intermediates
        # self.children = children

    
    

    # m_entry = m_feature_db["MP000005"]
    # print(m_entry)


    # # print("l_feature_hierarchy: ", l_feature_hierarchy.parents)
    # # print("l_feature_hierarchy: ", l_feature_hierarchy.intermediates)
    
    # # print("l_feature_hierarchy: ", l_feature_hierarchy)
    # print("l_feature_db: ", l_feature_db)
    # # print("l_ref_parent_order: ", len(l_ref_parent_order))
    # # print("m_feature_hierarchy: ", m_feature_hierarchy)
    # print("m_feature_db: ", m_feature_db)
    # # print("m_ref_parent_order: ", len(m_ref_parent_order))



def main(arglist=None):
    args = parse_args(arglist)
    print("Run Lifton!!")
    print(args)
    run_all_liftoff_steps(args)