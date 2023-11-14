from lifton import lifton_utils, logger

def miniprot_id_mapping(m_feature_db):
    m_id_dict = {}
    m_id_2_ref_id_trans_dict = {}
    for feature in m_feature_db.features_of_type("mRNA"):
        # Print all attributes and their values for the feature
        miniprot_id = feature["ID"][0]

        aa_trans_id = str(feature.attributes["Target"][0]).split(" ")[0]
        # print("aa_trans_id: ", aa_trans_id)
        if aa_trans_id in m_id_dict.keys():
            m_id_dict[aa_trans_id].append(miniprot_id)
        else:
            m_id_dict[aa_trans_id] = [miniprot_id]
        m_id_2_ref_id_trans_dict[miniprot_id] = aa_trans_id

    ###################################
    # Printing the miniprot dictionary
    ###################################
    # for key, vals in m_id_dict.items():
    #     print("key : ", key)
    #     print("vals: ", vals)


    # print(" m_feature_db.features_of_type('mRNA'):",  m_feature_db.all_features())
    # for feature in m_feature_db.features_of_type("mRNA"):
    #     print("feature ", feature)
    
    return m_id_dict, m_id_2_ref_id_trans_dict


def liftoff_id_mapping(l_feature_db, features):
    l_feature_id_2_ref_feature_id_dict = {}
    ref_feature_id_2_l_feature_id_dict = {}
    
    for feature in features:
        print("feature: ", feature)
        for locus in l_feature_db.features_of_type(feature):#, limit=("OX291666.1", 14237, 214237)):
            print("Locus: ", locus)
            liftoff_get_id_locus(locus, l_feature_db)

    print("End of liftoff_id_mapping!")

    # for feature in l_feature_db.features_of_type("mRNA"):
    #     # Print all attributes and their values for the feature
    #     miniprot_id = feature["ID"][0]

    #     aa_trans_id = str(feature.attributes["Target"][0]).split(" ")[0]
    #     # print("aa_trans_id: ", aa_trans_id)
    #     if aa_trans_id in l_id_dict.keys():
    #         l_id_dict[aa_trans_id].append(miniprot_id)
    #     else:
    #         l_id_dict[aa_trans_id] = [miniprot_id]
    #     l_id_2_ref_id_trans_dict[miniprot_id] = aa_trans_id

    ###################################
    # Printing the miniprot dictionary
    ###################################
    # for key, vals in m_id_dict.items():
    #     print("key : ", key)
    #     print("vals: ", vals)


    # print(" m_feature_db.features_of_type('mRNA'):",  m_feature_db.all_features())
    # for feature in m_feature_db.features_of_type("mRNA"):
    #     print("feature ", feature)
    
    # return m_id_dict, m_id_2_ref_id_trans_dict


def liftoff_get_id_locus(locus, l_feature_db):
    # Check if there are exons in the children
    children = list(l_feature_db.children(locus, featuretype='exon', level=1, order_by='start'))
    if len(children) == 0:
        # locus is in gene level: No exons in the children
        liftoff_gene_id, ref_gene_id = lifton_utils.get_ID(locus)
        logger.log(f"Liftoff: liftoff_gene_id\t: {liftoff_gene_id}\t{ref_gene_id}", debug=True)

        # ###########################
        # # 5.1 Create LifOn gene instance
        # ###########################
        # lifton_gene = lifton_class.Lifton_GENE(liftoff_gene_id, locus, ref_gene_info_dict[ref_gene_id], gene_copy_num_dict, tree_dict)
        # # Assign new gene ID with copy_number updated
        # liftoff_gene_id = lifton_gene.entry.id

        transcripts = l_feature_db.children(locus, level=1)
        for transcript in list(transcripts):
            liftoff_get_id_locus(transcript, l_feature_db)
        print("\n")

    else:
        # locus is in transcript level: exons are in its children
        print(f"\tTranscript level: {locus.id}")