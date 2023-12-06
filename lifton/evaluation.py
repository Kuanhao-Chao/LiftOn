from lifton import lifton_class, lifton_utils, align, logger
def tgt_evaluate(locus, db, tgt_fai, ref_features_dict, ref_proteins, DEBUG):
        # lifton_gene, locus, ref_db, l_feature_db, ref_id_2_m_id_trans_dict, m_feature_db, tree_dict, tgt_fai, ref_proteins, ref_trans, ref_features_dict, fw_score, DEBUG):
    # Check if there are exons in the children
    print("Inside tgt_evaluate!")
    # print(locus)
    exon_children = list(db.children(locus, featuretype='exon', level=1, order_by='start'))
    if len(exon_children) == 0:
        ###########################
        # These are middle features without exons
        #   => lifton_gene is not None & there are no exon children
        ###########################
        features = db.children(locus, level=1)
        for feature in list(features):
            tgt_evaluate(feature, db, tgt_fai, ref_features_dict, ref_proteins, DEBUG)

    else:
    # else:
    #     if lifton_gene is None:
    #         ###########################
    #         # These are gene (1st) features with direct exons 
    #         #   => lifton_gene is None & there are direct exon children 
    #         ###########################
    #         ref_gene_id, ref_trans_id = lifton_utils.get_ref_ids_liftoff(ref_features_dict, locus.id, None)
    #         ref_trans_id = ref_gene_id

    #         lifton_gene = lifton_class.Lifton_GENE(ref_gene_id, copy.deepcopy(locus), copy.deepcopy(ref_db[ref_gene_id].attributes), tree_dict, ref_features_dict, tmp = True)


    #     else:
    #         ###########################
    #         # These are transcript features with direct exons 
    #         #   => lifton_gene is not None & there are direct exon children 
    #         #   => (transcript level)
    #         ###########################
    #         ref_gene_id, ref_trans_id = lifton_utils.get_ref_ids_liftoff(ref_features_dict, lifton_gene.entry.id, locus.id)

    #     logger.log(f"\tTranscript level ref_gene_id\t: {ref_gene_id}; ref_trans_id\t:{ref_trans_id}", debug=DEBUG)

    #     ###########################
    #     # Processing transcript (feature with exon features)
    #     ###########################

    #     ###########################
    #     # Add LifOn transcript instance
    #     ###########################
    #     lifton_trans = lifton_gene.add_transcript(ref_trans_id, copy.deepcopy(locus), copy.deepcopy(ref_db[ref_trans_id].attributes))

    #     ###########################
    #     # Add LiftOn exons
    #     ###########################
    #     exons = l_feature_db.children(locus, featuretype='exon', order_by='start')
    #     for exon in list(exons):
    #         lifton_gene.add_exon(lifton_trans.entry.id, exon)
    #     logger.log(f"\tAfter adding exons\n", debug=DEBUG)

        ###########################
        # Add LiftOn CDS
        ###########################
        cdss = db.children(locus, featuretype='CDS', level=1, order_by='start') 
        cds_num = 0
        for cds in list(cdss):
            cds_num += 1
            # lifton_gene.add_cds(lifton_trans.entry.id, cds)
        logger.log(f"\tAfter adding CDSs\n", debug=DEBUG)
        
        # miniprot_aln, has_valid_miniprot = lifton_utils.LiftOn_check_miniprot_alignment(locus.seqid, locus, lifton_status, ref_id_2_m_id_trans_dict, m_feature_db, tree_dict, tgt_fai, ref_proteins, ref_trans_id)

        if (cds_num > 0):
            #############################################
            # Liftoff has protein
            #############################################                
            ref_gene_id, ref_trans_id = lifton_utils.get_ref_ids_liftoff(ref_features_dict, locus.id, locus.id)
            lifton_status = lifton_class.Lifton_Status()                
            liftoff_aln = align.parasail_align("lifton", db, locus, tgt_fai, ref_proteins, ref_trans_id, lifton_status)

            logger.log(f"After liftoff parasail_align\n", debug=DEBUG)

            if liftoff_aln is None:
                #############################################
                # There is no reference protein -> just keep Liftoff annotation
                #############################################
                logger.log("\t* Has CDS but no ref protein", debug=DEBUG)
                lifton_status.annotation = "target_no_ref_protein"
                lifton_status.status = ["no_ref_protein"]

            elif liftoff_aln.identity < 1:
                #############################################
                # Liftoff annotation is not perfect
                #############################################
            
                #############################################
                # Running chaining algorithm if there are valid miniprot alignments
                #############################################
                # if has_valid_miniprot:
                #     logger.log("\t* Has CDS and valid miniprot", debug=DEBUG)
                #     lifton_status.annotation = "LiftOn_chaining_algorithm" 
                #     cds_list = fix_trans_annotation.chaining_algorithm(liftoff_aln, miniprot_aln, tgt_fai)
                #     lifton_gene.update_cds_list(lifton_trans.entry.id, cds_list)
                #     logger.log("\tHas cds & protein & valid miniprot annotation!", debug=DEBUG)
                # else:
                #     logger.log("\t* has CDS but invalid miniprot", debug=DEBUG)
                #     lifton_status.annotation = "Liftoff_truncated"
                #     logger.log("\tHas cds & protein & invalid miniprot annotation!", debug=DEBUG)
                    
                #############################################
                # Step 3.6.1.2: Check if there are mutations in the transcript
                #############################################
                print("Score not perfect!")
                lifton_trans_aln, lifton_aa_aln = lifton_gene.fix_truncated_protein(lifton_trans.entry.id, ref_trans_id, tgt_fai, ref_proteins, ref_trans, lifton_status)                

                for mutation in lifton_status.status:
                    if mutation != "synonymous" and mutation != "identical" and mutation != "nonsynonymous":
                        lifton_aa_aln.write_alignment(intermediate_dir, "lifton_AA", mutation, transcript_id)
                        lifton_trans_aln.write_alignment(intermediate_dir, "lifton_DNA", mutation, transcript_id)

    #         elif liftoff_aln.identity == 1:
    #             #############################################
    #             # Step 3.6.2: Liftoff annotation is perfect
    #             #############################################
    #             lifton_status.annotation = "Liftoff_identical"
    #             # SETTING LiftOn identity score => Same as Liftoff
    #             lifton_status.lifton = liftoff_aln.identity
    #             lifton_status.status = ["identical"]
    #     else:
    #         #############################################
    #         # Liftoff no protein
    #         #############################################
    #         if has_valid_miniprot:
    #             ###########################
    #             # Condition 3: LiftOn does not have proteins & miniprot has proteins
    #             ###########################
    #             logger.log("\t* No CDS; miniprot has ref protein", debug=DEBUG)
    #             print(f"ref_trans_id: {ref_trans_id}; \t ref_id_2_m_id_trans_dict[ref_trans_id]: {ref_id_2_m_id_trans_dict[ref_trans_id]}")

    #             if ref_trans_id in ref_id_2_m_id_trans_dict.keys():
    #                 for mtrans_id in ref_id_2_m_id_trans_dict[ref_trans_id]:
    #                     mtrans = m_feature_db[mtrans_id]

    #                     mtrans_interval = Interval(mtrans.start, mtrans.end, mtrans_id)
    #                     is_overlapped = lifton_utils.check_ovps_ratio(mtrans, mtrans_interval, 0.7, tree_dict)
    #                     print("is_overlapped: ", is_overlapped)

    #                     if is_overlapped:
    #                         lifton_gene, transcript_id, lifton_status = run_miniprot.lifton_miniprot_with_ref_protein(mtrans, m_feature_db, ref_db, ref_gene_id, ref_trans_id, tgt_fai, ref_proteins, ref_trans, tree_dict, ref_features_dict, DEBUG)
    #                         lifton_trans = lifton_gene.transcripts[transcript_id]
    #                     else:
    #                         lifton_status.annotation = "Liftoff_no_protein"
    #                         lifton_status.status = ["no_protein"]
    #         else:
    #             ###########################
    #             # Condition 4: LiftOn does not have proteins & miniprot does not have proteins
    #             ###########################
    #             if (ref_trans_id in ref_proteins.keys()) :
    #                 logger.log("\t* No CDS & no miniprot but have ref protein", debug=DEBUG)
    #                 lifton_status.annotation = "Liftoff_no_ref_protein"
    #                 lifton_status.status = ["no_ref_protein"]
    #             else:
    #                 logger.log("\t* No CDS & no miniprot & no ref protein", debug=DEBUG)
    #                 lifton_status.annotation = "Liftoff_nc_transcript"
    #                 lifton_status.status = ["nc_transcript"]

    #     lifton_utils.write_lifton_status(fw_score, lifton_trans.entry.id, locus, lifton_status)

    #     lifton_gene.add_lifton_status_attrs(lifton_trans.entry.id, lifton_status)
        
    #     ###########################
    #     # Truncated reference proteins
    #     ###########################
    #     # if transcript_id_base in trunc_ref_proteins.keys() or transcript_id in trunc_ref_proteins.keys():
    #     #     print("Reference proteins is truncated.")
    #     #     lifton_status.annotation = "Reference_protein_truncated"
    #     #     lifton_status.status = ["truncated_ref_protein"]
        
    # return lifton_gene