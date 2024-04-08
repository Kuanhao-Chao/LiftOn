import copy
from lifton import lifton_class, lifton_class_eval, lifton_utils, align, logger

def tgt_evaluate(lifton_gene, locus, ref_db, db, tree_dict, tgt_fai, ref_features_dict, ref_proteins, ref_trans, fw_score, DEBUG):
    # Check if there are exons in the children
    exon_children = list(db.children(locus, featuretype='exon', level=1, order_by='start'))
    if len(exon_children) == 0:
        if lifton_gene is None:    
            ###########################
            # These are gene (1st) features without direct exons 
            #   => Create LifOn gene instance
            ########################### 
            ref_gene_id, ref_trans_id = lifton_utils.get_ref_ids_liftoff(ref_features_dict, "gene-"+(locus.id),None)

            # ref_gene_id, ref_trans_id = lifton_utils.get_ref_ids_liftoff(ref_features_dict, (locus.id),None)

            # logger.log(f"Before Gene level ref_gene_id\t: {ref_gene_id}; ref_trans_id\t:{ref_trans_id};  lifton_gene.copy_number\t:\n", debug=DEBUG)
            if ref_gene_id is None:
                return 
            
            lifton_gene = lifton_class_eval.Lifton_GENE_EVAL(ref_gene_id, copy.deepcopy(locus))
            logger.log(f"Gene level ref_gene_id\t: {ref_gene_id}; ref_trans_id\t:{ref_trans_id};", debug=DEBUG)

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

        # logger.log(f"\t\tlocus.id: {locus.id}", debug=DEBUG)

        if lifton_gene is None:
            ###########################
            # These are gene (1st) features with direct exons 
            #   => lifton_gene is None & there are direct exon children 
            ###########################

            ref_gene_id, ref_trans_id = lifton_utils.get_ref_ids_liftoff(ref_features_dict, "gene-"+(locus.id), None)

            # ref_gene_id, ref_trans_id = lifton_utils.get_ref_ids_liftoff(ref_features_dict, (locus.id), None)

            ref_trans_id = ref_gene_id

            if ref_gene_id is None:
                return 
            
            lifton_gene = lifton_class_eval.Lifton_GENE_EVAL(ref_gene_id, copy.deepcopy(locus))


        else:
            ###########################
            # These are transcript features with direct exons 
            #   => lifton_gene is not None & there are direct exon children 
            #   => (transcript level)
            ###########################

            ref_gene_id, ref_trans_id = lifton_utils.get_ref_ids_liftoff(ref_features_dict, "gene-"+(lifton_gene.entry.id), "rna-"+(locus.id))

            # ref_gene_id, ref_trans_id = lifton_utils.get_ref_ids_liftoff(ref_features_dict, lifton_gene.entry.id, (locus.id))


        logger.log(f"\tTranscript level ref_gene_id\t: {ref_gene_id}; ref_trans_id\t:{ref_trans_id}", debug=DEBUG)

        if ref_trans_id is None:
            return 
        
        ###########################
        # Processing transcript (feature with exon features)
        ###########################
        lifton_status = lifton_class.Lifton_Status()                
        try:
            ###########################
            # Add LifOn transcript instance
            ###########################
            lifton_trans = lifton_gene.add_transcript(ref_trans_id, copy.deepcopy(locus))
        except:
            return

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
        cdss = db.children(locus, featuretype=('CDS', 'stop_codon'), order_by='start') 
        cds_num = 0
        for cds in list(cdss):
            cds_num += 1
            lifton_gene.add_cds(lifton_trans.entry.id, cds)
        # logger.log(f"\tAfter adding CDSs\n", debug=DEBUG)
        
        if (cds_num > 0):
            #############################################
            # Liftoff has protein
            #############################################             
            if ref_trans_id not in ref_proteins.keys():
                #############################################
                # There is no reference protein -> just keep Liftoff annotation
                #############################################
                # logger.log("\t* Has CDS but no ref protein", debug=DEBUG)
                lifton_status.annotation = "Eval_no_ref_protein"
                lifton_status.status = ["no_ref_protein"]

            else:
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

        # print(">> written lifton_trans.entry.id: ", lifton_trans.entry.id)
        lifton_utils.write_lifton_eval_status(fw_score, lifton_trans.entry.id, locus, lifton_status)
