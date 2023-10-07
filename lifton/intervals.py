from intervaltree import Interval, IntervalTree

def initialize_interval_tree(l_feature_db):
    tree_dict = {}
    chr_num_ls = [*range(1, 23)] 
    chr_num_ls += ['X', 'Y', 'M']
    for i in chr_num_ls:
        tree = IntervalTree()
        tree_dict["chr" + str(i)] = tree
    
    ################################
    # Adding gene intervals into intervaltree
    ################################
    for gene in l_feature_db.features_of_type('gene'):
        gene_interval = Interval(gene.start, gene.end, gene.attributes["ID"][0])
        chromosome = gene.seqid
        tree_dict[chromosome].add(gene_interval)

    return tree_dict