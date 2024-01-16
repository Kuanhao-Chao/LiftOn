import subprocess
import os, copy
from lifton import lifton_class, logger, fix_trans_annotation, lifton_utils, run_miniprot, align
from lifton.liftoff import liftoff_main
from lifton.liftoff.tests import test_basic, test_advanced
from intervaltree import Interval, IntervalTree

def run_liftoff(output_dir, args):
    liftoff_args = copy.deepcopy(args)
    liftoff_outdir = output_dir + "liftoff/"    
    os.makedirs(liftoff_outdir, exist_ok=True)
    liftoff_annotation = liftoff_outdir + "liftoff.gff3"
    liftoff_args.output = liftoff_annotation
    liftoff_main.run_all_liftoff_steps(liftoff_args)
    # test_basic.test_yeast(liftoff_outdir + "test_basic/")
    # test_advanced.test_yeast(liftoff_outdir + "test_advance/")
    return liftoff_annotation


def process_liftoff(lifton_gene, locus, ref_db, l_feature_db, ref_id_2_m_id_trans_dict, m_feature_db, tree_dict, tgt_fai, ref_proteins, ref_trans, ref_features_dict, fw_score, fw_chain, write_chains, DEBUG):
    exon_children = list(l_feature_db.children(locus, featuretype='exon', level=1, order_by='start'))
    if len(exon_children) == 0:
        if lifton_gene is None:    
            # Gene (1st) features without direct exons 
            #   => Create LifOn gene instance
            ref_gene_id, ref_trans_id = lifton_utils.get_ref_ids_liftoff(ref_features_dict, locus.id, None)
            lifton_gene = lifton_class.Lifton_GENE(ref_gene_id, copy.deepcopy(locus), copy.deepcopy(ref_db[ref_gene_id].attributes), tree_dict, ref_features_dict)
            logger.log(f"Gene level ref_gene_id\t: {ref_gene_id}; ref_trans_id\t:{ref_trans_id};  lifton_gene.copy_number\t:{lifton_gene.copy_num}", debug=DEBUG)
            transcripts = l_feature_db.children(locus, level=1)
            for transcript in list(transcripts):
                process_liftoff(lifton_gene, transcript, ref_db, l_feature_db, ref_id_2_m_id_trans_dict, m_feature_db, tree_dict, tgt_fai, ref_proteins, ref_trans, ref_features_dict, fw_score, fw_chain, write_chains, DEBUG)
        else:
            # Middle features without exons
            #   => lifton_gene is not None & there are no exon children
            lifton_feature = lifton_gene.add_feature(copy.deepcopy(locus))                
            logger.log(f"\tother feature middle level lifton_feature\t: {lifton_feature.entry.id};", debug=DEBUG)
            features = l_feature_db.children(locus, level=1)
            for feature in list(features):
                process_liftoff(lifton_feature, feature, ref_db, l_feature_db, ref_id_2_m_id_trans_dict, m_feature_db, tree_dict, tgt_fai, ref_proteins, ref_trans, ref_features_dict, fw_score, fw_chain, write_chains, DEBUG)
    else:
        if lifton_gene is None:
            # Gene (1st) features with direct exons 
            #   => lifton_gene is None & there are direct exon children 
            ref_gene_id, ref_trans_id = lifton_utils.get_ref_ids_liftoff(ref_features_dict, locus.id, None)
            ref_trans_id = ref_gene_id
            lifton_gene = lifton_class.Lifton_GENE(ref_gene_id, copy.deepcopy(locus), copy.deepcopy(ref_db[ref_gene_id].attributes), tree_dict, ref_features_dict, tmp = True)
        else:
            # Transcript features with direct exons 
            #   => lifton_gene is not None & there are direct exon children 
            #   => (transcript level)
            ref_gene_id, ref_trans_id = lifton_utils.get_ref_ids_liftoff(ref_features_dict, lifton_gene.entry.id, locus.id)
        logger.log(f"\tTranscript level ref_gene_id\t: {ref_gene_id}; ref_trans_id\t:{ref_trans_id}", debug=DEBUG)
        # Processing transcript (feature with exon features)
        lifton_status = lifton_class.Lifton_Status()                
        # Add LifOn transcript instance
        lifton_trans = lifton_gene.add_transcript(ref_trans_id, copy.deepcopy(locus), copy.deepcopy(ref_db[ref_trans_id].attributes))
        exons = l_feature_db.children(locus, featuretype='exon', order_by='start')
        for exon in list(exons):
            lifton_gene.add_exon(lifton_trans.entry.id, exon)
        cdss = l_feature_db.children(locus, featuretype=('CDS', 'stop_codon'), order_by='start') 
        cds_num = 0
        for cds in list(cdss):
            cds_num += 1
            lifton_gene.add_cds(lifton_trans.entry.id, cds)
        # miniprot alignment
        miniprot_aln, has_valid_miniprot = lifton_utils.LiftOn_check_miniprot_alignment(lifton_trans, locus.seqid, locus, lifton_status, ref_id_2_m_id_trans_dict, m_feature_db, tree_dict, tgt_fai, ref_proteins, ref_trans_id)
        if (cds_num > 0):
            # Liftoff alignment                  
            liftoff_aln = align.parasail_align("liftoff", lifton_trans, locus, tgt_fai, ref_proteins, ref_trans_id, lifton_status)
            if liftoff_aln is None:
                # There is no reference protein -> just keep Liftoff annotation
                logger.log("\t* Has CDS but no ref protein", debug=DEBUG)
                lifton_status.annotation = "Liftoff_no_ref_protein"
                lifton_status.status = ["no_ref_protein"]
            elif liftoff_aln.identity == 1:
                # Liftoff protein annotation is perfect
                logger.log("\t* Liftoff protein identical", debug=DEBUG)
                lifton_status.lifton_aa = 1
                # DNA level alignment
                lifton_trans_aln = lifton_gene.align_trans_dna(lifton_trans.entry.id, ref_trans_id, tgt_fai, ref_trans, lifton_status)
            elif liftoff_aln.identity < 1:
                # Liftoff protein annotation is not perfect
                if has_valid_miniprot:
                    logger.log("\t* Has CDS and valid miniprot", debug=DEBUG)
                    lifton_status.annotation = "LiftOn_chaining_algorithm" 
                    cds_list, chains = fix_trans_annotation.chaining_algorithm(liftoff_aln, miniprot_aln, tgt_fai, debug=DEBUG)
                    if write_chains:
                        lifton_utils.write_lifton_chains(fw_chain, lifton_trans.entry.id, chains)
                    lifton_gene.update_cds_list(lifton_trans.entry.id, cds_list)
                else:
                    logger.log("\t* Has cds & protein but invalid miniprot annotation!", debug=DEBUG)
                    lifton_status.annotation = "Liftoff_truncated"                    
                # Check if there are mutations in the transcript
                lifton_trans_aln, lifton_aa_aln = lifton_gene.fix_truncated_protein(lifton_trans.entry.id, ref_trans_id, tgt_fai, ref_proteins, ref_trans, lifton_status)
        else:
            # Liftoff has no protein
            if has_valid_miniprot:
                # LiftOn does not have proteins & miniprot has proteins
                logger.log("\t* No CDS; miniprot has ref protein", debug=DEBUG)
                # logger.log(f"ref_trans_id: {ref_trans_id}; \t ref_id_2_m_id_trans_dict[ref_trans_id]: {ref_id_2_m_id_trans_dict[ref_trans_id]}", debug=DEBUG)
                if ref_trans_id in ref_id_2_m_id_trans_dict.keys():
                    for mtrans_id in ref_id_2_m_id_trans_dict[ref_trans_id]:
                        mtrans = m_feature_db[mtrans_id]
                        mtrans_interval = Interval(mtrans.start, mtrans.end, mtrans_id)
                        is_overlapped = lifton_utils.check_ovps_ratio(mtrans, mtrans_interval, 0.7, tree_dict)
                        if is_overlapped:
                            lifton_gene, transcript_id, lifton_status = run_miniprot.lifton_miniprot_with_ref_protein(mtrans, m_feature_db, ref_db, ref_gene_id, ref_trans_id, tgt_fai, ref_proteins, ref_trans, tree_dict, ref_features_dict, DEBUG)
                            lifton_trans = lifton_gene.transcripts[transcript_id]
                        else:
                            lifton_status.annotation = "Liftoff_no_protein"
                            lifton_status.status = ["no_protein"]
            else:
                # LiftOn does not have proteins & miniprot does not have proteins
                if (ref_trans_id in ref_proteins.keys()) :
                    logger.log("\t* No CDS & no miniprot but have ref protein", debug=DEBUG)
                    lifton_status.annotation = "Liftoff_no_ref_protein"
                    lifton_status.status = ["no_ref_protein"]
                else:
                    logger.log("\t* No CDS & no miniprot & no ref protein", debug=DEBUG)
                    lifton_status.annotation = "Liftoff_nc_transcript"
                    lifton_status.status = ["nc_transcript"]
        lifton_utils.print_lifton_status(lifton_trans.entry.id, locus, lifton_status, DEBUG=DEBUG)
        lifton_gene.add_lifton_status_attrs(lifton_trans.entry.id, lifton_status)
        # Writing out LiftOn entries & scores
        lifton_utils.write_lifton_status(fw_score, lifton_trans.entry.id, locus, lifton_status)
    return lifton_gene