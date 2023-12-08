import copy
from lifton import lifton_class, lifton_utils, align, logger

def tgt_evaluate(lifton_gene, locus, ref_db, db, tree_dict, tgt_fai, ref_features_dict, ref_proteins, ref_trans, fw_score, DEBUG):
    # Check if there are exons in the children
    exon_children = list(db.children(locus, featuretype='exon', level=1, order_by='start'))


    if len(exon_children) == 0:
        if lifton_gene is None:    
            ###########################
            # These are gene (1st) features without direct exons 
            #   => Create LifOn gene instance
            ########################### 
            ref_gene_id, ref_trans_id = lifton_utils.get_ref_ids_liftoff(ref_features_dict, locus.id, None)
            
            # logger.log(f"Before Gene level ref_gene_id\t: {ref_gene_id}; ref_trans_id\t:{ref_trans_id};  lifton_gene.copy_number\t:", debug=DEBUG)
            lifton_gene = lifton_class.Lifton_GENE(ref_gene_id, copy.deepcopy(locus), copy.deepcopy(ref_db[ref_gene_id].attributes), tree_dict, ref_features_dict)
            # logger.log(f"Gene level ref_gene_id\t: {ref_gene_id}; ref_trans_id\t:{ref_trans_id};  lifton_gene.copy_number\t:{lifton_gene.copy_num}", debug=DEBUG)

            transcripts = db.children(locus, level=1)
            for transcript in list(transcripts):
                tgt_evaluate(lifton_gene, transcript, ref_db, db, tree_dict, tgt_fai, ref_features_dict, ref_proteins, ref_trans, fw_score, DEBUG)

        else:
            ###########################
            # These are middle features without exons
            #   => lifton_gene is not None & there are no exon children
            ###########################
            lifton_feature = lifton_gene.add_feature(copy.deepcopy(locus))                
            # logger.log(f"\tMiddle other feature level lifton_feature\t: {lifton_feature.entry.id};", debug=DEBUG)
            features = db.children(locus, level=1)
            for feature in list(features):
                tgt_evaluate(lifton_feature, feature, ref_db, db, tree_dict, tgt_fai, ref_features_dict, ref_proteins, ref_trans, fw_score, DEBUG)



    else:

        if lifton_gene is None:
            ###########################
            # These are gene (1st) features with direct exons 
            #   => lifton_gene is None & there are direct exon children 
            ###########################
            ref_gene_id, ref_trans_id = lifton_utils.get_ref_ids_liftoff(ref_features_dict, locus.id, None)
            ref_trans_id = ref_gene_id

            lifton_gene = lifton_class.Lifton_GENE(ref_gene_id, copy.deepcopy(locus), copy.deepcopy(ref_db[ref_gene_id].attributes), tree_dict, ref_features_dict, tmp = True)


        else:
            ###########################
            # These are transcript features with direct exons 
            #   => lifton_gene is not None & there are direct exon children 
            #   => (transcript level)
            ###########################
            ref_gene_id, ref_trans_id = lifton_utils.get_ref_ids_liftoff(ref_features_dict, lifton_gene.entry.id, locus.id)

        # logger.log(f"\tTranscript level ref_gene_id\t: {ref_gene_id}; ref_trans_id\t:{ref_trans_id}", debug=DEBUG)

        ###########################
        # Processing transcript (feature with exon features)
        ###########################
        lifton_status = lifton_class.Lifton_Status()                

        ###########################
        # Add LifOn transcript instance
        ###########################
        lifton_trans = lifton_gene.add_transcript(ref_trans_id, copy.deepcopy(locus), copy.deepcopy(ref_db[ref_trans_id].attributes))

        ###########################
        # Add LiftOn exons
        ###########################
        exons = db.children(locus, featuretype='exon', order_by='start')
        for exon in list(exons):
            lifton_gene.add_exon(lifton_trans.entry.id, exon)
        # logger.log(f"\tAfter adding exons\n", debug=DEBUG)
        ###########################
        # Add LiftOn CDS
        ###########################
        cdss = db.children(locus, featuretype='CDS', order_by='start') 
        cds_num = 0
        for cds in list(cdss):
            cds_num += 1
            lifton_gene.add_cds(lifton_trans.entry.id, cds)
        # logger.log(f"\tAfter adding CDSs\n", debug=DEBUG)
        
        if (cds_num > 0):
            #############################################
            # Liftoff has protein
            #############################################                
            lifton_aln = align.parasail_align("lifton", db, locus, tgt_fai, ref_proteins, ref_trans_id, lifton_status)

            # logger.log(f"After liftoff parasail_align\n", debug=DEBUG)

            if lifton_aln is None:
                #############################################
                # There is no reference protein -> just keep Liftoff annotation
                #############################################
                # logger.log("\t* Has CDS but no ref protein", debug=DEBUG)
                lifton_status.annotation = "Eval_no_ref_protein"
                lifton_status.status = ["no_ref_protein"]

            elif lifton_aln.identity <= 1:

                #############################################
                # Step 3.6.1.2: Check if there are mutations in the transcript
                #############################################
                lifton_trans_aln, lifton_aa_aln = lifton_gene.fix_truncated_protein(lifton_trans.entry.id, ref_trans_id, tgt_fai, ref_proteins, ref_trans, lifton_status)
                
                if lifton_trans_aln is not None:
                    lifton_status.eval_dna = lifton_trans_aln.identity
                else:
                    lifton_status.eval_dna = 0
                if lifton_aa_aln is not None:
                    lifton_status.eval_aa = lifton_aa_aln.identity

                # print("lifton_trans_aln: ", lifton_trans_aln)
                # print("lifton_aa_aln: ", lifton_aa_aln)
        else:
            lifton_status.annotation = "Eval_no_cdss"
            lifton_status.status = ["no_cdss"]


        print("lifton_trans.entry.id: ", lifton_trans.entry.id)
        lifton_utils.write_lifton_eval_status(fw_score, lifton_trans.entry.id, locus, lifton_status)



        # ref_gene_id, ref_trans_id = lifton_utils.get_ref_ids_liftoff(ref_features_dict, lifton_gene.entry.id, locus.id)


        # ###########################
        # # Add LiftOn CDS
        # ###########################
        # cdss = db.children(locus, featuretype='CDS', level=1, order_by='start') 
        # cds_num = 0
        # for cds in list(cdss):
        #     cds_num += 1
        # logger.log(f"\tAfter adding CDSs\n", debug=DEBUG)
        
        # if (cds_num > 0):
        #     #############################################
        #     # Liftoff has protein
        #     #############################################    
        #     # Need to fix this part
        #     ref_trans_id = lifton_utils.get_ID_base(locus.id)
        #     lifton_status = lifton_class.Lifton_Status()                
        #     print(f"ref_trans_id: {ref_trans_id}")
            
        #     lifton_aln = align.parasail_align("lifton", db, locus, tgt_fai, ref_proteins, ref_trans_id, lifton_status)
        #     logger.log(f"After liftoff parasail_align\n", debug=DEBUG)

        #     if lifton_aln is None:
        #         #############################################
        #         # There is no reference protein -> just keep Liftoff annotation
        #         #############################################
        #         logger.log("\t* Has CDS but no ref protein", debug=DEBUG)
        #         lifton_status.annotation = "target_no_ref_protein"
        #         lifton_status.status = ["no_ref_protein"]
                
        #     elif lifton_aln.identity <= 1:
        #         #############################################
        #         # Step 3.6.1.2: Check if there are mutations in the transcript
        #         #############################################

        #         if self.ref_gene_id in ref_features_dict.keys():
        #             copy_num = ref_features_dict[self.ref_gene_id].copy_num


        #         Lifton_trans = lifton_class.Lifton_TRANS(ref_trans_id, ref_trans_id, locus.id, 1, locus, {})

        #         ###########################
        #         # Add LiftOn exons
        #         ###########################
        #         print("locus: ", locus)
        #         exons = db.children(locus, featuretype='exon', order_by='start')
        #         for exon in list(exons):
        #             Lifton_trans.add_exon(exon)
        #             print("exon: ", exon)
        #         logger.log(f"\tAfter adding exons\n", debug=DEBUG)
        #         ###########################
        #         # Add LiftOn CDS
        #         ###########################
        #         cdss = db.children(locus, featuretype='CDS', order_by='start') 
        #         cds_num = 0
        #         for cds in list(cdss):
        #             cds_num += 1
        #             Lifton_trans.add_cds(cds)
        #             print("cds: ", cds)
        #         logger.log(f"\tAfter adding CDSs\n", debug=DEBUG)


        #         lifton_trans_aln, lifton_aa_aln = Lifton_trans.fix_truncated_protein(tgt_fai, ref_proteins, ref_trans, lifton_status)                


        #         lifton_utils.write_lifton_status(fw_score, Lifton_trans.entry.id, locus, lifton_status)


        #         print(f"lifton_trans_aln: {lifton_trans_aln}")
        #         print(f"lifton_aa_aln: {lifton_aa_aln}")
        #         print(f"lifton_status.status: {lifton_status.status}")


        #         # for mutation in lifton_status.status:
        #         #     if mutation != "synonymous" and mutation != "identical" and mutation != "nonsynonymous":
        #         #         lifton_aa_aln.write_alignment(intermediate_dir, "lifton_AA", mutation, transcript_id)
        #         #         lifton_trans_aln.write_alignment(intermediate_dir, "lifton_DNA", mutation, transcript_id)
