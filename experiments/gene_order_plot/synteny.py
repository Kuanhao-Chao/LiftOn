import numpy as np
from plot_gene_order import *
from filepaths import *
from annotation import *
import nltk
import warnings


def analyze_synteny(ref_db, target_db,  ref_fa, target_fa, args, ref_features_dict):
    print('Analyzing synteny')
    output_file = filepaths.build_filepath([args.dir, filepaths.SYNTENY_OUTPUTS['gene_order']])
    if filepaths.make_file(output_file, args.force):
        if args.r_sort:
            ref_chrom_order = parse_chrom_order_list(args.r_sort, ref_fa)
            target_chrom_order= parse_chrom_order_list(args.t_sort, target_fa)
            ref_gene_order = order_genes(ref_chrom_order, target_chrom_order, ref_db, target_db, args.ft)
            target_gene_order = order_genes(target_chrom_order, ref_chrom_order, target_db, ref_db, args.ft)
        else:
            ref_gene_order, target_gene_order = infer_chrom_and_gene_order(ref_fa, target_fa, ref_db, target_db, args.ft)
        output_rows = print_order_output(ref_gene_order, target_gene_order, ref_db,target_db, output_file, args, ref_features_dict)
        plot_gene_order(output_rows, args)


def infer_chrom_and_gene_order(ref_fa, target_fa, ref_db, target_db, feature_types):
        less_contigs_fa, more_contigs_fa, less_contigs_db, more_contigs_db = find_most_contiguous_genome(ref_fa, target_fa, ref_db, target_db)
        default_chrom_order = order_chroms(less_contigs_fa)
        default_gene_order = order_genes(default_chrom_order, more_contigs_fa.keys(), less_contigs_db, more_contigs_db, feature_types)
        matched_chrom_order = match_chrom_order(default_gene_order, more_contigs_db, feature_types)
        matched_gene_order = order_genes(matched_chrom_order , less_contigs_fa.keys(), more_contigs_db, less_contigs_db, feature_types)
        return get_ref_and_target_gene_order(default_gene_order, matched_gene_order,ref_fa == more_contigs_fa)


def find_most_contiguous_genome(ref_fa, target_fa, ref_db, target_db):
    if len([seq for seq in ref_fa]) < len([seq for seq in target_fa]):
        return ref_fa, target_fa, ref_db, target_db
    return target_fa, ref_fa, target_db, ref_db


def get_ref_and_target_gene_order(default_gene_order, matched_gene_order, ref_has_more_contigs):
    if ref_has_more_contigs:
        return matched_gene_order, default_gene_order
    return default_gene_order, matched_gene_order


def order_chroms(fa):
    chroms = []
    for chrom,seq in fa.items():
        chroms.append(chrom)
    ordered_chroms = np.array(chroms)
    return get_chrom_order_dict(ordered_chroms)


def parse_chrom_order_list(chrom_order_txt_file, fa):
    chrom_order = []
    with open(chrom_order_txt_file) as co:
        for line in co:
            chrom = line.strip()
            if chrom not in fa:
                chrom_not_found_str = chrom +" not in assembly. Skipping"
                warnings.warn(chrom_not_found_str)
            else:
                chrom_order.append(line.strip())
    return get_chrom_order_dict(np.array(chrom_order))

def get_chrom_order_dict(ordered_chroms):
    chrom_order_dict = {}
    for i in range(len(ordered_chroms)):
        chrom_order_dict[ordered_chroms[i]] = i
    return chrom_order_dict


def match_chrom_order(ref_gene_order, target_db, feature_types):
    chrom_to_ref_gene_order = get_default_order_of_target_genes(target_db, ref_gene_order, feature_types)
    target_positions = get_median_default_gene_order(chrom_to_ref_gene_order)
    target_positions.sort(key = lambda x: x[1])
    ordered_chroms = np.array(target_positions)[:,0]
    return get_chrom_order_dict(ordered_chroms)


def get_default_order_of_target_genes(target_db, ref_gene_order, feature_types):
    chrom_to_ref_gene_order = {}
    all_genes = target_db.get_features_of_type(feature_types)
    for gene in all_genes:
        if gene.seqid not in chrom_to_ref_gene_order:
            chrom_to_ref_gene_order[gene.seqid] = [0]
        if gene.id in ref_gene_order:
            chrom_to_ref_gene_order[gene.seqid].append(ref_gene_order[gene.id])
    return chrom_to_ref_gene_order


def get_median_default_gene_order(chrom_to_ref_gene_order):
    chrom_and_avg_gene_number = []
    for chrom, gene_order_numbers in chrom_to_ref_gene_order.items():
        chrom_and_avg_gene_number.append([chrom, np.median(gene_order_numbers)])
    return chrom_and_avg_gene_number


def order_genes(chrom_order1, chrom_order2, gene_db1, gene_db2, feature_types):
    gene_order_dict = {}
    all_genes1 = [gene for gene in gene_db1.get_features_of_type(feature_types) if gene.seqid in chrom_order1]
    all_genes2 = [gene.id for gene in gene_db2.get_features_of_type(feature_types) if gene.seqid in chrom_order2]
    filtered_genes = [gene for gene in all_genes1 if  gene.id in all_genes2]
    filtered_genes.sort(key=lambda x: [chrom_order1[x.seqid], x.start])
    for i in range(len(filtered_genes)):
        gene_order_dict[filtered_genes[i].id] =i
    return gene_order_dict


def print_order_output(ref_gene_order, target_gene_order, ref_db, target_db, output_file, args, ref_features_dict):
    if args.edit_distance:
        edit_distance = nltk.edit_distance(list(ref_gene_order.values()), [target_gene_order[gene] for gene in
                                                                           ref_gene_order],transpositions=True)
    else:
        edit_distance = 'not calculated'
    output_rows = []
    target_gene_dict = target_db.get_feature_dict(args.ft)

    score_dict = get_scores(args, ref_features_dict)
    add_ref_genes_to_print(output_rows, ref_gene_order, target_gene_order, ref_db, target_gene_dict, args.ft, score_dict)
    with open(output_file, 'w') as f:
        f.write("Edit distance=" + str(edit_distance) + "\n" + "\n")
        for row in output_rows:
            f.write("\t".join([ str(col) for col in row]) + "\n")
    if len(output_rows) == 0:
        mismatch_ids_warning = "No features with matching IDs to plot"
        warnings.warn(mismatch_ids_warning)
    return output_rows


def add_ref_genes_to_print(output_rows, ref_gene_order, target_gene_order, ref_db, target_gene_dict, feature_types, score_dict):
    ref_gene_dict = ref_db.get_feature_dict(feature_types)
    for gene in ref_gene_order:
        ref_seq_name = ref_gene_dict[gene].seqid
        target_order = target_gene_order[gene]
        target_seq_name = target_gene_dict[gene].seqid
        # seq_id = get_perc_id(target_gene_dict[gene])
        seq_id = get_protein_score(gene, score_dict)
        if seq_id == -1:
            continue
        
        output_rows.append([gene, ref_gene_order[gene], target_order, ref_seq_name, target_seq_name, seq_id])

def get_protein_score(gene_id, score_dict):
    if gene_id in score_dict.keys():
        return score_dict[gene_id]
    else:
        return -1

def get_scores(args, ref_features_dict):
    score_fname = os.path.join(os.path.dirname(args.dir[:-1]), "score.txt")
    trans_scores = {}
    with open(score_fname, "r") as fr:
        lines = fr.read().splitlines()
        for line in lines:
            eles = line.split("\t")           
            trans_id = eles[0]
            lifton_score = float(eles[3])
            trans_scores[trans_id] = lifton_score

    # print("trans_scores: ", trans_scores)
    gene_scores = {}

    for gene, transs in ref_features_dict.items():
        max_score = 0
        for trans in transs.children:
            if trans not in trans_scores.keys():
                continue
            if trans_scores[trans] > max_score:
                max_score = trans_scores[trans]
        gene_scores[gene] = max_score

    return gene_scores
    # print(gene_scores)
    # print(len(gene_scores))
