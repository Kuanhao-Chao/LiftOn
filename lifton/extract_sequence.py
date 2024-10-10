from Bio.Seq import Seq

def extract_features(ref_db, features, ref_fai):
    ref_trans = {}
    ref_proteins = {}
    for feature in features:
        for locus in ref_db.db_connection.features_of_type(feature):
            __inner_extract_feature(ref_db, locus, ref_fai, ref_trans, ref_proteins)
    return ref_trans, ref_proteins


# def __inner_extract_feature(ref_db, feature, ref_fai, ref_trans, ref_proteins):
#     # If exon is the first level children
#     children_exons = list(ref_db.db_connection.children(feature, featuretype='exon', level=1))
#     children_CDSs = list(ref_db.db_connection.children(feature, featuretype=('start_codon', 'CDS', 'stop_codon'), level=1))
#     if len(children_exons) > 0 or len(children_CDSs) > 0:
#         # print(f"Parent: {feature.id};  exon: {len(children_exons)}; CDS: {len(children_CDSs)}")
#         if len(children_exons) > 0:
#             trans_seq = get_dna_sequence(feature, ref_fai, children_exons)
#             trans_seq = trans_seq.upper()
#             ref_trans[feature.id] = trans_seq
#         if len(children_CDSs) > 0:
#             protein_seq = get_protein_sequence(feature, ref_fai, children_CDSs)
#             protein_seq = protein_seq.upper()
#             ref_proteins[feature.id] = protein_seq
#     else:
#         for child in ref_db.db_connection.children(feature, level=1, order_by='start'):
#             __inner_extract_feature(ref_db, child, ref_fai, ref_trans, ref_proteins)


def __inner_extract_feature(ref_db, feature, ref_fai, ref_trans, ref_proteins):
    if feature.featuretype == 'gene':
        # For GTF files, transcripts are children of genes
        transcripts = list(ref_db.db_connection.children(feature, featuretype='transcript', order_by='start'))
        for transcript in transcripts:
            __inner_extract_feature(ref_db, transcript, ref_fai, ref_trans, ref_proteins)
    elif feature.featuretype == 'transcript':
        # Get exons and CDSs from transcripts
        exons = list(ref_db.db_connection.children(feature, featuretype='exon', order_by='start'))
        CDSs = list(ref_db.db_connection.children(feature, featuretype='CDS', order_by='start'))
        if exons:
            trans_seq = get_dna_sequence(feature, ref_fai, exons)
            trans_seq = trans_seq.upper()
            ref_trans[feature.id] = trans_seq
        if CDSs:
            protein_seq = get_protein_sequence(feature, ref_fai, CDSs)
            protein_seq = protein_seq.upper()
            ref_proteins[feature.id] = protein_seq
    else:
        # For GFF files, exons and CDSs might be direct children of the feature
        exons = list(ref_db.db_connection.children(feature, featuretype='exon', level=1, order_by='start'))
        CDSs = list(ref_db.db_connection.children(feature, featuretype='CDS', level=1, order_by='start'))
        if exons or CDSs:
            if exons:
                trans_seq = get_dna_sequence(feature, ref_fai, exons)
                trans_seq = trans_seq.upper()
                ref_trans[feature.id] = trans_seq
            if CDSs:
                protein_seq = get_protein_sequence(feature, ref_fai, CDSs)
                protein_seq = protein_seq.upper()
                ref_proteins[feature.id] = protein_seq
        else:
            # Recursively check for children in case of nested features
            for child in ref_db.db_connection.children(feature, order_by='start'):
                __inner_extract_feature(ref_db, child, ref_fai, ref_trans, ref_proteins)


def merge_children_intervals(children):
    if not children:
        return []
    intervals = [[child.start, child.end] for child in children]
    intervals.sort(key=lambda interval: interval[0])
    merged = [intervals[0]]
    for current in intervals[1:]:
        previous = merged[-1]
        if current[0] <= previous[1]:
            previous[1] = max(previous[1], current[1])
        else:
            merged.append(current)
    return merged


def get_dna_sequence(parent_feature, fasta, features):
    chrom = parent_feature.seqid
    strand = parent_feature.strand
    merged_features = merge_children_intervals(features)
    sequence =''
    if chrom not in fasta.keys():
        return sequence
    for start, end in merged_features:
        sequence += str(fasta[chrom][start -1: end])
    if strand == "-":
        sequence = str(Seq(sequence).reverse_complement())
    sequence += 'N' * get_padding_length(len(sequence))
    return sequence


def get_padding_length(sequence_length):
    return (3-sequence_length%3)%3


def get_protein_sequence(parent_feature, fasta, features):
    dna = get_dna_sequence(parent_feature, fasta, features)
    return str(Seq(dna).translate())