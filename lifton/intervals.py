from intervaltree import Interval, IntervalTree


def _make_interval(start, end, locus_id):
    """V2.8 fix: intervaltree disallows zero-length intervals
    (start == end) but NCBI GFF3 permits single-base features. Widen
    the half-open end by 1 so the tree builds without crashing.
    """
    if end <= start:
        end = start + 1
    return Interval(start, end, locus_id)


def initialize_interval_tree(l_feature_db, features):
    tree_dict = {}
    # Adding gene intervals into intervaltree
    for feature in features:
        for locus in l_feature_db.features_of_type(feature):#, limit=("chrX", 0, 200914237)):
            gene_interval = _make_interval(locus.start, locus.end, locus.id)
            chromosome = locus.seqid
            if chromosome not in tree_dict.keys():
                tree_dict[chromosome] = IntervalTree()
            tree_dict[chromosome].add(gene_interval)
    return tree_dict