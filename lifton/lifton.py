from lifton import extract_features, mapping, intervals, extra_copy, align, adjust_cds_boundaries, fix_trans_annotation, lifton_class, lifton_utils, annotation, sequence
from intervaltree import Interval, IntervalTree
import argparse
from argparse import Namespace
from pyfaidx import Fasta, Faidx
import copy, os
# from lifton import run_miniprot, sequence, annotation
# import subprocess
# import sys
# from Bio.Seq import Seq

# from liftoff.lifton import parse_chrm_files, get_parent_features_to_lift, liftover_types

import time


def args_outgrp(parser):
    outgrp = parser.add_argument_group('* Output settings')
    outgrp.add_argument(
        '-o', '--output', default='stdout', metavar='FILE',
        help='write output to FILE in same format as input; by default, output is written to terminal (stdout)'
    )
    outgrp.add_argument(
        '-u', default='unmapped_features.txt', metavar='FILE',
        help='write unmapped features to FILE; default is "unmapped_features.txt"',
    )
    outgrp.add_argument(
        '-exclude_partial', action='store_true',
        help='write partial mappings below -s and -a threshold to unmapped_features.txt; if true '
             'partial/low sequence identity mappings will be included in the gff file with '
             'partial_mapping=True, low_identity=True in comments'
    )
    outgrp.add_argument(
        '-dir', '--directory', default='intermediate_files', metavar='DIR',
        help='name of directory to save intermediate fasta and SAM files; default is "intermediate_files"',
    )
    return outgrp

def args_aligngrp(parser):
    aligngrp = parser.add_argument_group('Alignments')
    aligngrp.add_argument('-mm2_options', metavar='=STR', type=str, default='-a --end-bonus '
                                                                            '5 --eqx -N 50 '
                                                                            '-p 0.5',
                          help='space delimited minimap2 parameters. By default ="-a --end-bonus 5 --eqx -N 50 -p 0.5"')
    aligngrp.add_argument(
        '-a', default=0.5, metavar='A', type=float,
        help='designate a feature mapped only if it aligns with coverage ≥A; by default A=0.5',
    )
    aligngrp.add_argument(
        '-s', default=0.5, metavar='S', type=float,
        help='designate a feature mapped only if its child features (usually exons/CDS) align '
             'with sequence identity ≥S; by default S=0.5'
    )
    aligngrp.add_argument(
        '-d', metavar='D', default=2.0, type=float,
        help='distance scaling factor; alignment nodes separated by more than a factor of D in '
             'the target genome will not be connected in the graph; by default D=2.0'
    )
    aligngrp.add_argument(
        '-flank', default=0, metavar='F', type=float, help="amount of flanking sequence to align as a "
                                                           "fraction [0.0-1.0] of gene length. This can improve gene "
                                                           "alignment where gene structure  differs between "
                                                           "target and "
                                                           "reference; by default F=0.0")
    return aligngrp

def args_optional(parser):
    parser.add_argument('-V', '--version', help='show program version', action='version', version='v1.6.3')
    parser.add_argument(
        '-t', '--threads', default=1, type=int, metavar='THREADS', help='use t parallel processes to accelerate alignment; by default p=1'
    )
    parser.add_argument('-m', help='Minimap2 path', metavar='PATH')
    parser.add_argument('-f', '--features', metavar='TYPES', help='list of feature types to lift over')
    parser.add_argument(
        '-infer-genes', required=False, action='store_true',
        help='use if annotation file only includes transcripts, exon/CDS features'
    )
    parser.add_argument(
        '-infer_transcripts', action='store_true', required=False,
        help='use if annotation file only includes exon/CDS features and does not include transcripts/mRNA'
    )
    parser.add_argument(
        '-chroms', metavar='TXT', help='comma seperated file with corresponding chromosomes in '
                                       'the reference,target sequences',
    )
    parser.add_argument(
        '-unplaced', metavar='TXT',
        help='text file with name(s) of unplaced sequences to map genes from after genes from '
             'chromosomes in chroms.txt are mapped; default is "unplaced_seq_names.txt"',
    )
    parser.add_argument('-copies', action='store_true', help='look for extra gene copies in the target genome')
    parser.add_argument(
        '-sc', default=1.0, metavar='SC', type=float,
        help='with -copies, minimum sequence identity in exons/CDS for which a gene is considered '
             'a copy; must be greater than -s; default is 1.0',
    )
    parser.add_argument('-overlap', default=0.1, metavar='O', help="maximum fraction [0.0-1.0] of overlap allowed by 2 "
                                                                   "features; by default O=0.1", type=float)
    parser.add_argument('-mismatch', default=2, metavar='M', help="mismatch penalty in exons when finding best "
                                                                  "mapping; by default M=2", type=int)
    parser.add_argument('-gap_open', default=2, metavar='GO', help="gap open penalty in exons when finding best "
                                                                   "mapping; by default GO=2", type=int)
    parser.add_argument('-gap_extend', default=1, metavar='GE', help="gap extend penalty in exons when finding best "
                                                                     "mapping; by default GE=1", type=int)
    parser.add_argument('-subcommand', required=False,  help=argparse.SUPPRESS)
    parser.add_argument('-polish', required=False, action='store_true', default = False)
    parser.add_argument('-cds', required=False, action="store_true", default=True, help="annotate status of each CDS "
                                                                                        "(partial, missing start, "
                                                                                        "missing stop, inframe stop "
                                                                                        "codon)")



def parse_args(arglist):
    parser = argparse.ArgumentParser(description='Lift features from one genome assembly to another')
    parser.add_argument('target', help='target fasta genome to lift genes to')
    parser.add_argument('reference', help='reference fasta genome to lift genes from')
    
    parser_outgrp = args_outgrp(parser)
    parser_aligngrp = args_aligngrp(parser)
    args_optional(parser)

    referencegrp = parser.add_argument_group('* Required input (Reference annotation)')    
    referencegrp.add_argument(
        '-g', '--reference-annotation', metavar='GFF',  required=True,
        help='the reference annotation file to lift over in GFF or GTF format (or) '
                'name of feature database; if not specified, the -g '
                'argument must be provided and a database will be built automatically'
    )

    referencegrp.add_argument(
        '-P', '--proteins', metavar='FASTA', required=False, default=None,
        help='the reference protein sequences.'
    )

    referencegrp.add_argument(
        '-T', '--transcripts', metavar='FASTA', required=False, default=None,
        help='the reference transcript sequences.'
    )

    ###################################
    # START for the LiftOn algorithm
    ###################################
    liftoffrefrgrp = parser.add_argument_group('* Optional input (Liftoff annotation)')
    liftoffrefrgrp.add_argument(
        '-L', '--liftoff', metavar='gff', default=None,
        help='the annotation generated by Liftoff (or) '
                'name of Liftoff gffutils database; if not specified, the -liftoff '
                'argument must be provided and a database will be built automatically'
    )

    miniprotrefrgrp = parser.add_argument_group('* Optional input (miniprot annotation)')
    miniprotrefrgrp.add_argument(
        '-M', '--miniprot', metavar='gff', default=None,
        help='the annotation generated by miniprot (or) '
                'name of miniprot gffutils database; if not specified, the -miniprot '
                'argument must be provided and a database will be built automatically'
    )

    ###################################
    # END for the LiftOn algorithm
    ###################################

    parser._positionals.title = '* Required input (sequences)'
    parser._optionals.title = '* Miscellaneous settings'
    # parser._action_groups = [referencegrp, parser_outgrp, parser._optionals]

    parser._action_groups = [parser._positionals, referencegrp, liftoffrefrgrp, miniprotrefrgrp, parser_outgrp, parser._optionals, parser_aligngrp]
    args = parser.parse_args(arglist)

    if '-a' not in args.mm2_options:
        args.mm2_options += ' -a'
    if '--eqx' not in args.mm2_options:
        args.mm2_options += ' --eqx'
    if '-N' not in args.mm2_options:
        args.mm2_options += " -N 50"
    if '-p' not in args.mm2_options:
        args.mm2_options += " -p 0.5"
    if '--end-bonus' not in args.mm2_options:
        args.mm2_options += "--end-bonus 5"
    if (float(args.s) > float(args.sc)):
        parser.error("-sc must be greater than or equal to -s")
    if (args.chroms is None and args.unplaced is not None):
        parser.error("-unplaced must be used with -chroms")

    return args


def run_all_lifton_steps(args):
    ################################
    # Step 0: Reading target & reference genomes
    ################################
    tgt_genome = args.target
    ref_genome = args.reference

    # Setting output directory & intermediate_dir
    if args.output == "stdout": 
        outdir = "."
    else:
        outdir = os.path.dirname(args.output)
        outdir = outdir if outdir is not "" else "."
        os.makedirs(outdir, exist_ok=True)

    intermediate_dir = f"{outdir}/intermediate_files/"
    os.makedirs(intermediate_dir, exist_ok=True)
    args.directory = intermediate_dir
    
    print(f">> outdir           : {outdir}")
    print(f">> intermediate_dir : {intermediate_dir}")

    print(">> Reading target genome ...")
    tgt_fai = Fasta(tgt_genome)
    tgt_chromosomes = tgt_fai.keys()
    print(">> Reading reference genome ...")
    ref_fai = Fasta(ref_genome)
    ref_chromosomes = ref_fai.keys()
    print(args)

    ################################
    # Step 1: Building database from the reference annotation
    ################################
    ref_db = annotation.Annotation(args.reference_annotation, args.infer_genes)
    # feature_types = lifton_utils.get_feature_types(args.features)
    # This is for special annotation file where there are no exon type
    # child_types = lifton_utils.get_child_types(feature_types, ref_db)

    # Get all the parent features to liftover    
    features = lifton_utils.get_parent_features_to_lift(args.features)

    ################################
    # Step 1: Get all the parent features to liftover
    ################################
    ref_features_dict, ref_trans_2_gene_dict, ref_gene_info_dict, ref_trans_info_dict = lifton_utils.get_ref_liffover_features(features, ref_db)
        


    # ################################
    # # Step 2: Creating protein & DNA dictionaries from the reference annotation
    # ################################
    ref_proteins_file = args.proteins    
    if ref_proteins_file is None or not os.path.exists(ref_proteins_file):
        print(">> Creating transcript protein dictionary from the reference annotation ...")
        ref_proteins = sequence.SequenceDict(ref_db, ref_fai, ['CDS', 'start_codon', 'stop_codon'], True)        
        ref_proteins_file = lifton_utils.write_seq_2_file(intermediate_dir, ref_proteins, "proteins")

    else:
        print(">> Reading transcript protein dictionary from the reference fasta ...")
        ref_proteins = Fasta(ref_proteins_file)
    print("\t * number of proteins: ", len(ref_proteins.keys()))

    trunc_ref_proteins = lifton_utils.get_truncated_protein(ref_proteins)
    trunc_ref_proteins_file = lifton_utils.write_seq_2_file(intermediate_dir, trunc_ref_proteins, "truncated_proteins")

    ref_trans_file = args.transcripts    
    if ref_trans_file is None or not os.path.exists(ref_trans_file):
        print(">> Creating transcript DNA dictionary from the reference annotation ...")
        ref_trans = sequence.SequenceDict(ref_db, ref_fai, ['exon'], False)
        ref_trans_file = lifton_utils.write_seq_2_file(intermediate_dir, ref_trans, "transcripts")

    else:
        print(">> Reading transcript DNA dictionary from the reference fasta ...")
        ref_trans = Fasta(ref_trans_file)

    print(">> Creating transcript DNA dictionary from the reference annotation ...")
    print("\t * number of transcripts: ", len(ref_trans.keys()))






    # ################################
    # # Step TEST: integrate Liftoff inside LiftOn
    # ################################
    # # if args.chroms is not None:
    # #     ref_chroms, target_chroms = parse_chrm_files(args.chroms)
    # # else:
    # #     ref_chroms = [args.reference]
    # #     target_chroms = [args.target]
    # # parent_features_to_lift = get_parent_features_to_lift(args.features)
    # # lifted_feature_list = {}
    # # unmapped_features = []
    # # feature_db, feature_hierarchy, ref_parent_order = liftover_types.lift_original_annotation(ref_chroms, target_chroms,
    # #                                                                                           lifted_feature_list, args,
    # #                                                                                           unmapped_features,
    # #                                                                                           parent_features_to_lift)
    # # unmapped_features = map_unmapped_features(unmapped_features, target_chroms, lifted_feature_list, feature_db,
    # #                                           feature_hierarchy, ref_parent_order, args)
    # # map_features_from_unplaced_seq(unmapped_features, lifted_feature_list, feature_db, feature_hierarchy,
    # #                                ref_parent_order, args)
    # # write_unmapped_features_file(args.u, unmapped_features)
    # # map_extra_copies(args, lifted_feature_list, feature_hierarchy, feature_db, ref_parent_order)

    # # if args.cds and args.polish is False:
    # #     check_cds(lifted_feature_list, feature_hierarchy, args)
    # # if args.polish:
    # #      print("polishing annotations")
    # #      check_cds(lifted_feature_list, feature_hierarchy, args)
    # #      write_new_gff.write_new_gff(lifted_feature_list, args, feature_db)
    # #      find_and_polish_broken_cds(args, lifted_feature_list,feature_hierarchy, ref_chroms,
    # #                                                       target_chroms,
    # #                            unmapped_features, feature_db, ref_parent_order)
    # #      if args.output != 'stdout':
    # #          args.output += "_polished"
    # # write_new_gff.write_new_gff(lifted_feature_list, args, feature_db)











    ################################
    # Step 3: Run liftoff & miniprot
    ################################
    liftoff_annotation = lifton_utils.exec_liftoff(outdir, args)
    miniprot_annotation = lifton_utils.exec_miniprot(outdir, args, tgt_genome, ref_proteins_file)


    ################################
    # Step 4: Run LiftOn algorithm
    ################################
    ################################
    # Step 4.0: Create liftoff and miniprot database
    ################################
    print(">> Creating liftoff database : ", liftoff_annotation)
    l_feature_db = annotation.Annotation(liftoff_annotation, args.infer_genes).db_connection
    
    print(">> Creating miniprot database : ", miniprot_annotation)
    m_feature_db = annotation.Annotation(miniprot_annotation, args.infer_genes).db_connection

    # liftoff_liftover_count = 0
    # for feature in features:
    #     for gene in l_feature_db.features_of_type(feature):#, limit=("CM033155.1", 0, 
    #         gene_id = gene.attributes["ID"][0]
            
    #         if gene_id in ref_features_dict.keys():
    #             ref_features_dict[gene_id] = True

    #         liftoff_liftover_count += 1
    
    # print("liftoff_liftover_count: ", liftoff_liftover_count)


    # miniprot_liftover_count = 0
    # for feature in ["mRNA"]:
    #     for gene in m_feature_db.features_of_type(feature):#, limit=("CM033155.1", 0, 
    #         miniprot_liftover_count += 1
    
    # print("miniprot_liftover_count: ", miniprot_liftover_count)




    fw = open(args.output, "w")
    fw_lifton_fixed_perfect = open(args.output+".perfect.gff3", "w")
    fw_lifton_fixed_unperfect = open(args.output+".unperfect.gff3", "w")

    fw_score = open(outdir+"/score.txt", "w")

    ################################
    # Step 4.1: Creating miniprot 2 Liftoff ID mapping
    ################################
    m_id_dict, m_id_2_ref_id_trans_dict = mapping.id_mapping(m_feature_db)
    # print("m_id_dict; ", m_id_dict)
    # print("m_id_2_ref_id_trans_dict; ", m_id_2_ref_id_trans_dict)
    
    ################################
    # Step 4.2: Initializing intervaltree
    ################################
    tree_dict = intervals.initialize_interval_tree(l_feature_db)

    # # Dictionary for extra copy
    gene_copy_num_dict = {}
    trans_copy_num_dict = {}

    LIFTOFF_BAD_PROT_TRANS_COUNT = 0
    LIFTOFF_GOOD_PROT_TRANS_COUNT = 0
    LIFTOFF_NC_TRANS_COUNT = 0

    MINIPROT_BAD_PROT_TRANS_COUNT = 0
    MINIPROT_GOOD_PROT_TRANS_COUNT = 0

    LIFTON_TOTAL_GENE_COUNT = 0
    LIFTON_TOTAL_TRANS_COUNT = 0
    LIFTON_BAD_PROT_TRANS_COUNT = 0
    LIFTON_GOOD_PROT_TRANS_COUNT = 0
    LIFTON_C_TRANS_LOST_COUNT = 0
    LIFTON_NC_TRANS_COUNT = 0
    LIFTON_C_TRANS_NOREF_COUNT = 0

    LIFTON_MINIPROT_FIXED_GENE_COUNT = 0


    ################################
    # Step 4.3: Iterate gene entries & fixing CDS lists
    ################################
    # For missing transcripts.
    gene_copy_num_dict["gene-LiftOn"] = 0
    
    # For transcripts without CDSs
    tmp_outdir = intermediate_dir + "/tmp/"
    os.makedirs(tmp_outdir, exist_ok=True)
    fw_truncted_trans = open(tmp_outdir+"/truncated_protein_trans.txt", "w")
    fw_no_ref_trans = open(tmp_outdir+"/no_ref_protein_trans.txt", "w")
    fw_ref_trans_loss = open(tmp_outdir+"/ref_protein_trans_loss.txt", "w")
    fw_nc_trans = open(tmp_outdir+"/nc_trans.txt", "w")


    ################################
    # Step 5: Process Liftoff genes & transcripts
    ################################
    for feature in features:
        for gene in l_feature_db.features_of_type(feature):#, limit=("CM033155.1", 0, 877531)):
            LIFTON_TOTAL_GENE_COUNT += 1
            gene_id, gene_id_base = lifton_utils.get_ID(gene)
            chromosome = gene.seqid

            # Update gene copy number dictionary
            lifton_utils.update_copy(gene_id_base, gene_copy_num_dict)

            ###########################
            # Create LifOn gene instance
            ###########################
            lifton_gene = lifton_class.Lifton_GENE(gene)
            
            # iterate through Liftoff transcripts
            #   Assumption: all 1st level are transcripts
            transcripts = l_feature_db.children(gene, level=1)
            for transcript in list(transcripts):
                LIFTON_TOTAL_TRANS_COUNT += 1
                transcript_id, transcript_id_base = lifton_utils.get_ID(transcript)

                # create status for each LiftOn transcript
                lifton_status = lifton_class.Lifton_Status()
                
                ###########################
                # Add LifOn transcript instance
                ###########################
                lifton_gene.add_transcript(transcript)

                # Mark transcript as lifted-over
                ref_features_dict[gene_id_base][transcript_id_base] = True

                # Update transcript copy number dictionary
                lifton_utils.update_copy(transcript_id_base, trans_copy_num_dict)
                
                print("\ttranscript_id\t: ", transcript_id)
                if transcript_id != transcript_id_base:
                    print("&& transcript_id_base\t: ", transcript_id_base)


                ###########################
                # Add LiftOn exons
                ###########################
                exons = l_feature_db.children(transcript, featuretype='exon')  # Replace 'exon' with the desired child feature type
                for exon in list(exons):
                    lifton_gene.add_exon(transcript_id, exon)

                ###########################
                # Add LiftOn CDS
                ###########################
                cdss = l_feature_db.children(transcript, featuretype='CDS')  # Replace 'exon' with the desired child feature type
                cds_num = 0
                for cds in list(cdss):
                    cds_num += 1
                    lifton_gene.add_cds(transcript_id, cds)


                if (cds_num > 0) and (transcript_id_base in ref_proteins.keys()):
                    liftoff_aln = align.parasail_align("liftoff", l_feature_db, transcript, tgt_fai, ref_proteins, transcript_id_base)

                    # SETTING Liftoff identity score
                    lifton_status.liftoff = liftoff_aln.identity

                    if liftoff_aln.identity < 1:
                        LIFTOFF_BAD_PROT_TRANS_COUNT += 1
                        #############################################
                        # Step 3.6.1: Liftoff annotation is not perfect
                        #############################################
                        lifton_status.annotation = "Liftoff_truncated"
                    
                        miniprot_aln, has_valid_miniprot = lifton_utils.LiftOn_check_miniprot_alignment(chromosome, transcript, lifton_status, m_id_dict, m_feature_db, tree_dict, tgt_fai, ref_proteins, transcript_id_base)

                        #############################################
                        # Step 3.6.1.1: Running chaining algorithm if there are valid miniprot alignments
                        #############################################
                        if has_valid_miniprot:
                            LIFTON_MINIPROT_FIXED_GENE_COUNT += 1
                            cds_list = fix_trans_annotation.chaining_algorithm(liftoff_aln, miniprot_aln, tgt_fai, fw)
                            lifton_gene.update_cds_list(transcript_id, cds_list)
                            lifton_status.annotation = "Liftoff_miniprot_chaining_algorithm" 
                        else:
                            print("Has cds & protein & miniprot annotation; but no valid miniprot annotation!")
                            
                        #############################################
                        # Step 3.6.1.2: Check if there are mutations in the transcript
                        #############################################
                        lifton_trans_aln, lifton_aa_aln = lifton_gene.fix_truncated_protein(transcript_id, transcript_id_base, tgt_fai, ref_proteins, ref_trans, lifton_status)
                        # SETTING LiftOn identity score
                        
                        if lifton_aa_aln.identity == 1:
                            LIFTON_GOOD_PROT_TRANS_COUNT += 1
                            # Writing out perfect LiftOn annotation
                            lifton_gene.write_entry(fw_lifton_fixed_perfect)


                        elif lifton_aa_aln.identity < 1:
                            LIFTON_BAD_PROT_TRANS_COUNT += 1
                            # Writing out truncated LiftOn annotation
                            lifton_gene.write_entry(fw_lifton_fixed_unperfect)
                            
                            for mutation in lifton_status.status:
                                if mutation != "synonymous" and mutation != "identical" and mutation != "nonsynonymous":
                                    lifton_aa_aln.write_alignment(intermediate_dir, "lifton_AA", mutation, transcript_id)
                                    lifton_trans_aln.write_alignment(intermediate_dir, "lifton_DNA", mutation, transcript_id)

                    elif liftoff_aln.identity == 1:
                        #############################################
                        # Step 3.6.2: Liftoff annotation is perfect
                        #############################################
                        LIFTOFF_GOOD_PROT_TRANS_COUNT += 1
                        LIFTON_GOOD_PROT_TRANS_COUNT += 1

                        miniprot_aln, has_valid_miniprot = lifton_utils.LiftOn_check_miniprot_alignment(chromosome, transcript, lifton_status, m_id_dict, m_feature_db, tree_dict, tgt_fai, ref_proteins, transcript_id_base)

                        # SETTING LiftOn identity score => Same as Liftoff
                        lifton_status.lifton = liftoff_aln.identity
                        lifton_status.annotation = "Liftoff_identical"
                        lifton_status.status = ["identical"]


                    if lifton_status.miniprot == 1:
                        MINIPROT_GOOD_PROT_TRANS_COUNT += 1
                    elif lifton_status.miniprot < 1:
                        MINIPROT_BAD_PROT_TRANS_COUNT += 1                    


                    final_status = ";".join(lifton_status.status)
                    fw_score.write(f"{transcript_id}\t{lifton_status.liftoff}\t{lifton_status.miniprot}\t{lifton_status.lifton}\t{lifton_status.annotation}\t{final_status}\t{transcript.seqid}:{transcript.start}-{transcript.end}\n")

                elif transcript_id_base in trunc_ref_proteins.keys() or transcript_id in trunc_ref_proteins.keys():
                    print("Reference proteins is truncated.")
                    lifton_status.annotation = "Reference_protein_truncated"
                    lifton_status.status = ["truncated_ref_protein"]
                    fw_truncted_trans.write(transcript_id+"\n")

                else:
                    # Only liftoff annotation
                    if (cds_num > 0) and (transcript_id_base not in ref_proteins.keys()):
                        print("Have CDS but no ref protein")
                        LIFTON_C_TRANS_NOREF_COUNT += 1
                        lifton_status.annotation = "No_reference_protein"
                        lifton_status.status = ["no_ref_protein"]
                        fw_no_ref_trans.write(transcript_id+"\n")

                    elif (cds_num == 0) and (transcript_id_base in ref_proteins.keys() or transcript_id in ref_proteins.keys()) :
                        print("No CDS but have ref protein")
                        LIFTON_C_TRANS_LOST_COUNT += 1
                        lifton_status.annotation = "Liftoff_protein_transcript_loss"
                        lifton_status.status = ["protein_transcript_loss"]
                        fw_ref_trans_loss.write(transcript_id+"\n")
                    
                    else:
                        print("No CDS and no ref protein")
                        LIFTOFF_NC_TRANS_COUNT += 1
                        LIFTON_NC_TRANS_COUNT += 1
                        lifton_status.annotation = "Liftoff_nc_transcript"
                        lifton_status.status = ["nc_transcript"]
                        fw_nc_trans.write(transcript_id+"\n")


            ###########################
            # Step 4.7: Writing out LiftOn entries
            ###########################
            # if lifton_status.lifton == 1 and lifton_status.liftoff < 1 and lifton_status.miniprot < 1:
            lifton_gene.write_entry(fw)

            ###########################
            # Step 4.8: Adding LiftOn intervals
            ###########################
            gene_interval = Interval(lifton_gene.entry.start, lifton_gene.entry.end, gene_id)
            tree_dict[chromosome].add(gene_interval)


    ################################
    # Step 6: Process miniprot transcripts
    ################################
    # EXTRA_COPY_MINIPROT_COUNT = 0 
    # NEW_LOCUS_MINIPROT_COUNT = 0
    # EXTRA_COPY_MINIPROT_COUNT, NEW_LOCUS_MINIPROT_COUNT = extra_copy.find_extra_copy(m_feature_db, tree_dict, m_id_2_ref_id_trans_dict, ref_gene_info_dict, ref_trans_info_dict, ref_trans_2_gene_dict, gene_copy_num_dict, trans_copy_num_dict, fw, ref_features_dict)

    for mtrans in m_feature_db.features_of_type('mRNA'):
        chromosome = mtrans.seqid
        mtrans_id = mtrans.attributes["ID"][0]
        mtrans_interval = Interval(mtrans.start, mtrans.end, mtrans_id)
        ovps = tree_dict[chromosome].overlap(mtrans_interval)
    

        if len(ovps) == 0:

            ref_trans_id = m_id_2_ref_id_trans_dict[mtrans_id]
            gene_entry_base = copy.deepcopy(mtrans)
            trans_entry_base = copy.deepcopy(mtrans)

            if ref_trans_id in ref_trans_info_dict.keys():
                
                ref_gene_id = ref_trans_2_gene_dict[ref_trans_id] 
                # Update gene copy number dictionary
                lifton_utils.update_copy(ref_gene_id, gene_copy_num_dict)

                ###########################
                # Create LifOn gene instance
                ###########################
                gene_attrs = ref_gene_info_dict[ref_gene_id]
                lifton_gene = lifton_class.Lifton_GENE(gene_entry_base)
                new_extra_cp_gene_id = lifton_gene.update_gene_info(ref_gene_id, chromosome, mtrans.start, mtrans.end, gene_attrs, gene_copy_num_dict)

                ###########################
                # Create LifOn transcript instance
                ###########################
                trans_attrs = ref_trans_info_dict[ref_trans_id]
                new_extra_cp_trans_id = lifton_gene.create_new_transcript(False, ref_gene_id, ref_trans_id, trans_entry_base, chromosome, mtrans.start, mtrans.end, trans_attrs, gene_copy_num_dict, trans_copy_num_dict)

                ref_features_dict[ref_gene_id][ref_trans_id] = True

                #######################################
                # Step 5.3: Create exon / CDS entries
                #######################################
                cdss = m_feature_db.children(mtrans, featuretype='CDS')  # Replace 'exon' with the desired child feature type
                # print("cdss len: ", len(cdss))
                for cds in list(cdss):
                    # entry.attributes["ID"]
                    lifton_gene.add_exon(new_extra_cp_trans_id, cds)
                    cds_copy = copy.deepcopy(cds)
                    lifton_gene.add_cds(new_extra_cp_trans_id, cds_copy)


                # create status for each LiftOn transcript
                lifton_status = lifton_class.Lifton_Status()
                                

                miniprot_aln, has_valid_miniprot = lifton_utils.LiftOn_check_miniprot_alignment(chromosome, mtrans, lifton_status, m_id_dict, m_feature_db, tree_dict, tgt_fai, ref_proteins, transcript_id_base)

                lifton_status.liftoff = 0
                lifton_status.miniprot = miniprot_aln.identity
                lifton_status.lifton = lifton_status.miniprot
                lifton_status.annotation =  "miniprot"
                lifton_status.status = [""]

                # Write scores for each transcript



            else:
                # Reason it's missing => the mRNA does not belong to gene (vdj segments) || the mRNA is not in the reference annotation
                NEW_LOCUS_MINIPROT_COUNT += 1
                #######################################
                # Step 5.1: Create the gene entry
                #######################################
                ref_gene_id = "gene-LiftOn"
                lifton_gene = lifton_class.Lifton_GENE(gene_entry_base)
                gene_attrs = lifton_class.Lifton_GENE_info({}, ref_gene_id)
                new_extra_cp_gene_id = lifton_gene.update_gene_info(ref_gene_id, chromosome, mtrans.start, mtrans.end, gene_attrs, gene_copy_num_dict)

                #######################################
                # Step 5.2: Create the transcript entry
                #######################################
                trans_attrs = lifton_class.Lifton_TRANS_info({}, ref_trans_id, ref_gene_id)
                new_extra_cp_trans_id = lifton_gene.create_new_transcript(True, ref_gene_id, ref_trans_id, trans_entry_base, chromosome, mtrans.start, mtrans.end, trans_attrs, gene_copy_num_dict, trans_copy_num_dict)

                #######################################
                # Step 5.3: Create the exon entry
                #######################################
                #######################################
                # Step 5.4: Create the CDS entry
                #######################################
                cdss = m_feature_db.children(mtrans, featuretype='CDS')  # Replace 'exon' with the desired child feature type
                # print("cdss len: ", len(cdss))
                for cds in list(cdss):
                    lifton_gene.add_exon(ref_trans_id, cds)
                    cds_copy = copy.deepcopy(cds)
                    lifton_gene.add_cds(ref_trans_id, cds_copy)
            
            # Add new transcript loci into interval tree
            tree_dict[chromosome].add(mtrans_interval)
            lifton_gene.write_entry(fw)



    print("*********************************************************************")
    print("LiftOn total gene loci\t\t\t\t: ", LIFTON_TOTAL_GENE_COUNT)
    print("LiftOn total transcript\t\t\t\t: ", LIFTON_TOTAL_TRANS_COUNT)
    print("LiftOn good protein trans count\t\t\t: ", LIFTON_GOOD_PROT_TRANS_COUNT)
    print("LiftOn unperfect protein trans count\t\t: ", LIFTON_BAD_PROT_TRANS_COUNT)
    print("\tLiftOn miniprot fixed trans count\t\t: ", LIFTON_MINIPROT_FIXED_GENE_COUNT)
    print("LiftOn non-coding trans count\t\t\t: ", LIFTON_NC_TRANS_COUNT)
    print("LiftOn coding trans with no reference count\t\t\t: ", LIFTON_C_TRANS_NOREF_COUNT)
    print("LiftOn coding trans with protein lost count\t\t\t: ", LIFTON_C_TRANS_LOST_COUNT)

    print("\n")
    print("\tLiftoff good protein trans count\t\t\t: ", LIFTOFF_GOOD_PROT_TRANS_COUNT)
    print("\tLiftoff unperfect protein trans count\t\t\t: ", LIFTOFF_BAD_PROT_TRANS_COUNT)
    print("\tLiftoff non-coding trans count\t\t\t\t: ", LIFTOFF_NC_TRANS_COUNT)

    print("\n")
    print("\tminiprot good protein trans count\t\t\t: ", MINIPROT_GOOD_PROT_TRANS_COUNT)
    print("\tminiprot unperfect protein trans count\t\t\t: ", MINIPROT_BAD_PROT_TRANS_COUNT)
    print("*********************************************************************")

    fw.close()
    fw_truncted_trans.close()
    fw_no_ref_trans.close()
    fw_ref_trans_loss.close()
    fw_nc_trans.close()

def main(arglist=None):
    args = parse_args(arglist)
    print("Run Lifton!!")
    run_all_lifton_steps(args)