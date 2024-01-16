from intervaltree import Interval, IntervalTree

def initialize_interval_tree(l_feature_db, features):
    tree_dict = {}
    # Adding gene intervals into intervaltree
    for feature in features:
        for locus in l_feature_db.features_of_type(feature):#, limit=("chrX", 0, 200914237)):
            gene_interval = Interval(locus.start, locus.end, locus.id)
            chromosome = locus.seqid
            if chromosome not in tree_dict.keys():
                tree_dict[chromosome] = IntervalTree()
            tree_dict[chromosome].add(gene_interval)
    return tree_dict