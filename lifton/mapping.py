def id_mapping(m_feature_db):
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
