def get_ref_liffover_features(features, ref_db):
    ref_features_dict = {}
    ref_features_reverse_dict = {}
    # ref_trans_2_gene_dict = {}
    # gene_info_dict = {}
    # trans_info_dict = {}
    new_gene_feature = Lifton_feature("Lifton-gene")
    ref_features_dict["LiftOn-gene"] = new_gene_feature

    for f_itr in features:
        for locus in ref_db.db_connection.features_of_type(f_itr):
            feature = Lifton_feature(locus.id)
            exon_children = list(ref_db.db_connection.children(locus, featuretype='exon', level=1, order_by='start'))

            if len(exon_children) > 0:
                process_ref_liffover_features(locus, ref_db, None)
            else:
                transcripts = ref_db.db_connection.children(locus, level=1)
                for transcript in list(transcripts):
                    # print(transcript)
                    process_ref_liffover_features(transcript, ref_db, feature)
                    ref_features_reverse_dict[transcript.id] = locus.id
            ref_features_dict[locus.id] = feature

    return ref_features_dict, ref_features_reverse_dict

def process_ref_liffover_features(locus, ref_db, feature):
    if feature != None:
        feature.children.add(locus.id)


class Lifton_feature:
    def __init__(self, id):
        self.id = id
        self.copy_num = 0
        self.children = set()

