import numpy as np
from plot_gene_order import *
from filepaths import *
from annotation import *
import nltk
import warnings


def flip_gene_order_for_chromosome(chrom_to_gene_order, ref_gene_order, chrom_to_gene_dict):
    for chrom, gene_order_dict in chrom_to_gene_order.items():
        if len(gene_order_dict) > 1:
            gene_order_list = list(gene_order_dict.values())
            max_order = max(gene_order_list)
            min_order = min(gene_order_list)

            gene_order_list = list(gene_order_dict.values())
            # Check if almost all genes are in reverse order
            
            print("Check the reverse chromosome order: ", chrom)
            if is_reverse_order(gene_order_dict, max_order, min_order, ref_gene_order):
            
                # Flip the gene order for this chromosome
                print(f"Reverse {chrom}")
                for gene_id, order in gene_order_dict.items():
                    print("Before: ", gene_order_dict[gene_id])
                    gene_order_dict[gene_id] = max_order - order + min_order
                    print("\tref_gene_order: ", ref_gene_order[gene_id])
                    print("After: ", gene_order_dict[gene_id])
                    print("\n")
                

def remove_outliers(list1, list2, threshold=0.4):
    # Calculate the mean and standard deviation for each list
    mean1 = sum(list1) / len(list1)
    mean2 = sum(list2) / len(list2)
    std1 = (sum((x - mean1) ** 2 for x in list1) / len(list1)) ** 0.5
    std2 = (sum((x - mean2) ** 2 for x in list2) / len(list2)) ** 0.5

    # Identify outliers based on the threshold (e.g., 2 standard deviations)
    outliers1 = [index for index, value in enumerate(list1) if abs(value - mean1) > threshold * std1]
    outliers2 = [index for index, value in enumerate(list2) if abs(value - mean2) > threshold * std2]

    # Combine outlier indices from both lists
    outlier_indices = set(outliers1 + outliers2)

    # Remove outliers at the same indices from both lists
    filtered_list1 = [value for index, value in enumerate(list1) if index not in outlier_indices]
    filtered_list2 = [value for index, value in enumerate(list2) if index not in outlier_indices]

    return filtered_list1, filtered_list2

def is_reverse_order(order_dict, max_order, min_order, ref_gene_order):
    # Check if at least 50% of the order list is in descending order
    sum_fwd_diff = 0
    sum_rvs_diff = 0

    # min_ref_order = 0
    # gene_count = 0
    tgt_ordering = []
    ref_ordering = []
    for geneid, tgt_order in order_dict.items():
        ref_order = ref_gene_order[geneid]
        tgt_ordering.append(tgt_order)
        ref_ordering.append(ref_order)
        print(f"{tgt_order};  {ref_order}")
    
    sub_tgt_ordering, sub_ref_ordering = remove_outliers(tgt_ordering, ref_ordering)
    print("sub_tgt_ordering: ", sub_tgt_ordering)
    print("sub_ref_ordering: ", sub_ref_ordering)
    correlation_coefficient = np.corrcoef(sub_tgt_ordering, sub_ref_ordering)[0, 1]

    # Interpret the correlation coefficient
    if correlation_coefficient > 0:
        print("Positive correlation")
        return False
    elif correlation_coefficient < 0:
        print("Negative correlation")
        return True

    #     # # min_ref_order = min(min_ref_order, ref_order)
    #     # # gene_count += 1

    #     # sum_fwd_diff += (abs(tgt_order - min_order - ref_order))
    #     # sum_rvs_diff += (abs(abs(max_order - tgt_order) - min_order - ref_order))   
    #     # print(f"tgt_order: {tgt_order}; ref_order: {ref_order}")
    
    # if sum_fwd_diff > sum_rvs_diff:
    #     print("is_reverse_order: ", True)
    #     return True
    # else:
    #     print("is_reverse_order: ", False)
    #     return False

def get_chrom_to_gene_order(gene_order, gene_db, feature_types):
    chrom_to_gene_order = {}
    all_genes = gene_db.get_features_of_type(feature_types)
    for gene in all_genes:
        if gene.seqid not in chrom_to_gene_order:
            chrom_to_gene_order[gene.seqid] = {}
        if gene.id in gene_order:
            chrom_to_gene_order[gene.seqid][gene.id] = gene_order[gene.id]
    return chrom_to_gene_order

# def get_tgt_2_ref_gene_order_diff(chrom_to_target_gene_order, ref_gene_order):

def concatenate_genes(chrom_to_gene_order):
    concatenated_gene_order = {}
    for chrom, gene_order_dict in chrom_to_gene_order.items():
        # Concatenate genes for each chromosome
        concatenated_gene_order.update(gene_order_dict)
    return concatenated_gene_order

def calculate_chromosome_gene_ratio(chrom_to_gene_order, total_gene):
    chr_gene_ratio = {}
    for chrom, gene_order_dict in chrom_to_gene_order.items():
        gene_on_chrom = len(gene_order_dict)
        ratio = gene_on_chrom / total_gene
        chr_gene_ratio[chrom] = ratio
        # if total_genes > 0:
        #     reverse_genes = sum(1 for order in gene_order_dict.values() if order != 0)
        #     forward_genes = total_genes - reverse_genes
        #     ratio = forward_genes / total_genes

        #     print(f"Chromosome {chrom} - Forward Genes: {forward_genes}, Reverse Genes: {reverse_genes}, Ratio: {ratio:.2f}")
    return chr_gene_ratio




def analyze_synteny(ref_db, target_db,  ref_fa, target_fa, args, ref_features_dict):
    print('Analyzing synteny')
    output_file = filepaths.build_filepath([args.dir, filepaths.SYNTENY_OUTPUTS['gene_order']])
    if filepaths.make_file(output_file, args.force):
        score_dict = get_scores(args, ref_features_dict)
        if args.r_sort:
            ref_chrom_order = parse_chrom_order_list(args.r_sort, ref_fa)
            target_chrom_order= parse_chrom_order_list(args.t_sort, target_fa)
            ref_gene_order = order_genes(ref_chrom_order, target_chrom_order, ref_db, target_db, args.ft)
            target_gene_order = order_genes(target_chrom_order, ref_chrom_order, target_db, ref_db, args.ft)
        else:
            ref_gene_order, target_gene_order = infer_chrom_and_gene_order(ref_fa, target_fa, ref_db, target_db, args.ft, score_dict)




        # Create dictionaries to store gene order for each chromosome
        chrom_to_ref_gene_order = get_chrom_to_gene_order(ref_gene_order, ref_db, args.ft)
        chrom_to_target_gene_order = get_chrom_to_gene_order(target_gene_order, target_db, args.ft)

        ref_chr_gene_ratio = calculate_chromosome_gene_ratio(chrom_to_ref_gene_order, len(ref_gene_order))
        tgt_chr_gene_ratio = calculate_chromosome_gene_ratio(chrom_to_target_gene_order, len(target_gene_order))

        # get_tgt_2_ref_gene_order_diff(chrom_to_target_gene_order, ref_gene_order)

        print("chrom_to_ref_gene_order    : ", chrom_to_ref_gene_order)
        print("chrom_to_target_gene_order : ", chrom_to_target_gene_order)
        
        # Flip gene order for affected chromosomes
        flip_gene_order_for_chromosome(chrom_to_target_gene_order, ref_gene_order, target_db.get_feature_dict(args.ft))

        # Concatenate genes corresponding to the chromosome orders
        concatenated_target_gene_order = concatenate_genes(chrom_to_target_gene_order)

        print("concatenated_target_gene_order: ", concatenated_target_gene_order)
        output_rows = print_order_output(ref_gene_order, concatenated_target_gene_order, ref_db,target_db, output_file, args, ref_features_dict, score_dict)

        # print("ref_gene_order    : ", ref_gene_order)
        # print("target_gene_order : ", target_gene_order)
        # output_rows = print_order_output(ref_gene_order, target_gene_order, ref_db,target_db, output_file, args, ref_features_dict, score_dict)
        plot_gene_order(output_rows, args, ref_chr_gene_ratio, tgt_chr_gene_ratio)


def infer_chrom_and_gene_order(ref_fa, target_fa, ref_db, target_db, feature_types, score_dict):
    less_contigs_fa, more_contigs_fa, less_contigs_db, more_contigs_db = find_most_contiguous_genome(ref_fa, target_fa, ref_db, target_db)
    default_chrom_order = order_chroms(less_contigs_fa)

    default_gene_order = order_genes(default_chrom_order, more_contigs_fa.keys(), less_contigs_db, more_contigs_db, feature_types, score_dict)
    matched_chrom_order = match_chrom_order(default_gene_order, more_contigs_db, feature_types)

    print(f"Target chromosomes ordering  : {len(less_contigs_fa.keys())}\t{less_contigs_fa.keys()}" )
    print(f"Reference chromsome ordering : {len(matched_chrom_order)}\t{matched_chrom_order}")
    matched_gene_order = order_genes(matched_chrom_order , less_contigs_fa.keys(), more_contigs_db, less_contigs_db, feature_types, score_dict)
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


def order_genes(chrom_order1, chrom_order2, gene_db1, gene_db2, feature_types, score_dict):
    gene_order_dict = {}
    all_genes1 = [gene for gene in gene_db1.get_features_of_type(feature_types) if (gene.seqid in chrom_order1) and (gene.id in score_dict.keys()) and (score_dict[gene.id] != 0)]
    all_genes2 = [gene.id for gene in gene_db2.get_features_of_type(feature_types) if (gene.seqid in chrom_order2) and (gene.id in score_dict.keys()) and (score_dict[gene.id] != 0)]
    filtered_genes = [gene for gene in all_genes1 if  gene.id in all_genes2]
    filtered_genes.sort(key=lambda x: [chrom_order1[x.seqid], x.start])
    for i in range(len(filtered_genes)):
        gene_order_dict[filtered_genes[i].id] =i
    return gene_order_dict


def print_order_output(ref_gene_order, target_gene_order, ref_db, target_db, output_file, args, ref_features_dict, score_dict):
    if args.edit_distance:
        edit_distance = nltk.edit_distance(list(ref_gene_order.values()), [target_gene_order[gene] for gene in
                                                                           ref_gene_order],transpositions=True)
    else:
        edit_distance = 'not calculated'
    output_rows = []
    target_gene_dict = target_db.get_feature_dict(args.ft)

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
