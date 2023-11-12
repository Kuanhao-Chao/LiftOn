from intervaltree import Interval, IntervalTree

def initialize_interval_tree(l_feature_db):
    tree_dict = {}

    ################################
    # Adding gene intervals into intervaltree
    ################################
    for gene in l_feature_db.features_of_type('gene'):
        gene_interval = Interval(gene.start, gene.end, gene.attributes["ID"][0])
        chromosome = gene.seqid

        if chromosome not in tree_dict.keys():
            tree_dict[chromosome] = IntervalTree()
        tree_dict[chromosome].add(gene_interval)
    return tree_dict