import sys, os
import annotation

def miniprot_id_mapping(m_feature_db):
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

    ###################################
    # Printing the miniprot dictionary
    ###################################
    # for key, vals in m_id_dict.items():
    #     print("key : ", key)
    #     print("vals: ", vals
    return ref_id_2_m_id_trans_dict, m_id_2_ref_id_trans_dict


def segments_overlap(segment1, segment2):
    # Check if the segments have valid endpoints
    # print("Checking two segments overlapping.!")

    # print(segment1, segment2)

    if len(segment1) != 2 or len(segment2) != 2:
        raise ValueError("Segments must have exactly 2 endpoints")
    
    # Sort the segments by their left endpoints
    segment1, segment2 = sorted([segment1, segment2], key=lambda x: x[0])


    # Check if the right endpoint of the first segment is greater than or equal to the left endpoint of the second segment
    # print(segment1[1] >= segment2[0])

    return segment1[1] >= segment2[0]


def main():

    ################################
    # Step 4: Run liftoff & miniprot
    ################################
    target = sys.argv[1]
    liftoff_annotation = f"/ccb/salz2/kh.chao/Lifton/results/{target}/liftoff/liftoff.gff3_db"
    miniprot_annotation = f"/ccb/salz2/kh.chao/Lifton/results/{target}/miniprot/miniprot.gff3_db"


    ################################
    # Step 5: Run LiftOn algorithm
    ################################
    ################################
    # Step 5.0: Create liftoff and miniprot database
    ################################
    l_feature_db = annotation.Annotation(liftoff_annotation, True).db_connection
    m_feature_db = annotation.Annotation(miniprot_annotation, True).db_connection

    print("l_feature_db: ", l_feature_db)
    print("m_feature_db: ", m_feature_db)


    ################################
    # Step 5.1: Creating miniprot 2 Liftoff ID mapping
    ################################
    ref_id_2_m_id_trans_dict, m_id_2_ref_id_trans_dict = miniprot_id_mapping(m_feature_db)

    # print("ref_id_2_m_id_trans_dict: ", ref_id_2_m_id_trans_dict)
    # print("m_id_2_ref_id_trans_dict: ", m_id_2_ref_id_trans_dict)
    
    out_dirname = f"/ccb/salz2/kh.chao/Lifton/results/{target}/visualization/circos/"
    os.makedirs(out_dirname, exist_ok=True)
    out_fname_circos = f"{out_dirname}/links_nonovp.csv"
    out_fname_circos_ovp = f"{out_dirname}/links_ovp.csv"

    out_f = open(out_fname_circos, "w")
    out_f_ovp = open(out_fname_circos_ovp, "w")

    out_f.write("chr1\tstart1\tend1\tchr2\tstart2\tend2\n")
    out_f_ovp.write("chr1\tstart1\tend1\tchr2\tstart2\tend2\n")

    # one_2_one_mapping = 0
    # one_2_multi_mapping = 0

    ovp_one_2_one_mapping_trans_count = 0
    ovp_multi__mapping_trans_count = 0
    nonovp_trans_count = 0

    for ref_id, m_ids in ref_id_2_m_id_trans_dict.items():
        try :
            ovp = False            
            ref_start = l_feature_db[ref_id].start
            ref_end = l_feature_db[ref_id].end
            ref_chr = l_feature_db[ref_id].seqid[3:] + "_t"
            ref_feature = f"{ref_chr},{ref_start},{ref_end}"


            for m_id in m_ids:
                m_start = m_feature_db[m_id].start
                m_end = m_feature_db[m_id].end
                m_chr = m_feature_db[m_id].seqid[3:] + "_r"
                m_feature = f"{m_chr},{m_start},{m_end}"

                if (l_feature_db[ref_id].seqid[3:] != m_feature_db[m_id].seqid[3:]) or not segments_overlap((int(ref_start), int(ref_end)), (int(m_start), int(m_end))):
                    entry = f"{ref_feature},{m_feature}"
                    print(entry)
                    out_f.write(f"{entry}\n")
                elif (l_feature_db[ref_id].seqid[3:] == m_feature_db[m_id].seqid[3:]) and segments_overlap((int(ref_start), int(ref_end)), (int(m_start), int(m_end))):
                    ovp = True
                    entry = f"{ref_feature},{m_feature}"
                    out_f_ovp.write(f"{entry}\n")


            if ovp and len(m_ids) == 1:
                ovp_one_2_one_mapping_trans_count += 1
            elif ovp and len(m_ids) > 1:
                ovp_multi__mapping_trans_count += 1
            else:
                nonovp_trans_count += 1
        except:
            print("ref_id: ", ref_id)
            print("ref_id not found in liftoff database")
            continue

    out_f.close()
    out_f_ovp.close()

    print(f"Liftoff & miniprot 1-1 mapping: {ovp_one_2_one_mapping_trans_count}")
    print(f"Liftoff & miniprot at least 1 mapping: {ovp_multi__mapping_trans_count}")
    print(f"Liftoff & miniprot non-mapping: {nonovp_trans_count}")






    # # out_fname_circos = sys.argv[1]
    # with open(out_fname_circos, "w") as out_f:
    #     # with open("/home/choh1/lifton_circos/human_refseq_test_links.csv") as f:
    #     #     f.readline()
    #     #     for line in f:
    #     #         line = line.rstrip().split(",")
    #     #         name1 = line[0].split('_')[0]


if __name__ == "__main__":
    main()