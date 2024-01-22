from lifton import lifton_utils, logger

def miniprot_id_mapping(m_feature_db):
    """
        This function creates a dictionary of miniprot id to reference id.

        Parameters:
        - m_feature_db: miniprot feature database

        Returns:
        ref_id_2_m_id_trans_dict: reference id to miniprot transcript ids dictionary
        m_id_2_ref_id_trans_dict: miniprot transcript id to reference id dictionary
    """
    ref_id_2_m_id_trans_dict = {}
    m_id_2_ref_id_trans_dict = {}
    for feature in m_feature_db.features_of_type("mRNA"):
        # Print all attributes and their values for the feature
        miniprot_id = feature["ID"][0]

        aa_trans_id = str(feature.attributes["Target"][0]).split(" ")[0]
        # print("aa_trans_id: ", aa_trans_id)
        if aa_trans_id in ref_id_2_m_id_trans_dict.keys():
            ref_id_2_m_id_trans_dict[aa_trans_id].append(miniprot_id)
        else:
            ref_id_2_m_id_trans_dict[aa_trans_id] = [miniprot_id]
        m_id_2_ref_id_trans_dict[miniprot_id] = aa_trans_id
    # Printing the miniprot dictionary
    # for key, vals in m_id_dict.items():
    #     print("key : ", key)
    #     print("vals: ", vals)
    return ref_id_2_m_id_trans_dict, m_id_2_ref_id_trans_dict