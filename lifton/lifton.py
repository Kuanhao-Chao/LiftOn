from lifton import extract_features, mapping, intervals, extra_copy, align, adjust_cds_boundaries, fix_trans_annotation, lifton_class, lifton_utils, annotation, sequence, stats, logger, run_miniprot, run_liftoff, __version__
from intervaltree import Interval, IntervalTree
import argparse
from argparse import Namespace
from pyfaidx import Fasta, Faidx
import copy, os
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
    parser.add_argument('-V', '--version', help='show program version', action='version', version=__version__)
    parser.add_argument('-D', '--debug', action='store_true', help='Run debug mode')
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
    DEBUG = args.debug
    # Setting output directory & intermediate_dir
    if args.output == "stdout": 
        outdir = "."
    else:
        outdir = os.path.dirname(args.output)
        outdir = outdir if outdir != "" else "."
        os.makedirs(outdir, exist_ok=True)

    intermediate_dir = f"{outdir}/intermediate_files/"
    os.makedirs(intermediate_dir, exist_ok=True)
    args.directory = intermediate_dir
    
    print(">> Reading target genome ...")
    tgt_fai = Fasta(tgt_genome)
    print(">> Reading reference genome ...")
    ref_fai = Fasta(ref_genome)

    ################################
    # Step 1: Building database from the reference annotation
    ################################
    ref_db = annotation.Annotation(args.reference_annotation, args.infer_genes)
    # Get all the parent features to liftover    
    features = lifton_utils.get_parent_features_to_lift(args.features)


    ################################
    # Step 1: Get all the parent features to liftover
    ################################
    ref_features_dict, ref_trans_2_gene_dict, ref_gene_info_dict, ref_trans_info_dict = lifton_utils.get_ref_liffover_features(features, ref_db)
    logger.log("ref_trans_2_gene_dict: ", len(ref_trans_2_gene_dict.keys()), debug=DEBUG)
    logger.log("ref_gene_info_dict: ", len(ref_gene_info_dict.keys()), debug=DEBUG)
    logger.log("ref_trans_info_dict: ", len(ref_trans_info_dict.keys()), debug=DEBUG)
        

    ################################
    # Step 2: Creating protein & DNA dictionaries from the reference annotation
    ################################
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
    fw = open(args.output, "w")
    fw_score = open(outdir+"/score.txt", "w")
    fw_unmapped = open(outdir+"/unmapped_features.txt", "w")
    fw_extra_copy = open(outdir+"/extra_copy_features.txt", "w")

    
    ################################
    # Step 4.1: Creating miniprot 2 Liftoff ID mapping
    ################################
    m_id_dict, m_id_2_ref_id_trans_dict = mapping.id_mapping(m_feature_db)

    
    ################################
    # Step 4.2: Initializing intervaltree
    ################################
    tree_dict = intervals.initialize_interval_tree(l_feature_db)
    # # Dictionary for extra copy
    gene_copy_num_dict = {}
    trans_copy_num_dict = {}


    ################################
    # Step 4.3: Iterate gene entries & fixing CDS lists
    ################################
    # For missing transcripts.
    gene_copy_num_dict["gene-LiftOn"] = 0


    ################################
    # Step 5: Process Liftoff genes & transcripts
    ################################
    for feature in features:
        for gene in l_feature_db.features_of_type(feature):#, limit=("CP101152.1", 22584243, 22587514)):
            liftoff_gene_id, ref_gene_id = lifton_utils.get_ID(gene)
            logger.log(f"Liftoff: liftoff_gene_id\t: {liftoff_gene_id}\t{ref_gene_id}\n", debug=DEBUG)
            ###########################
            # 5.1 Create LifOn gene instance
            ###########################
            lifton_gene = lifton_class.Lifton_GENE(liftoff_gene_id, gene, ref_gene_info_dict[ref_gene_id], gene_copy_num_dict, tree_dict)
            # Assign new gene ID with copy_number updated
            liftoff_gene_id = lifton_gene.entry.id
            
            ###########################
            # 5.2 iterate through Liftoff transcripts
            #   Assumption: all 1st level are transcripts
            ###########################
            transcripts = l_feature_db.children(gene, level=1)
            for transcript in list(transcripts):
                transcript_id, ref_trans_id = lifton_utils.get_ID(transcript)
                # Mark transcript as lifted-over
                ref_features_dict[ref_gene_id][ref_trans_id] = True

                logger.log(f"Liftoff: transcript_id\t: {transcript_id}\t{ref_trans_id}\n", debug=DEBUG)

                # create status for each LiftOn transcript
                lifton_status = lifton_class.Lifton_Status()                

                ###########################
                # Add LifOn transcript instance
                ###########################

                transcript_id = lifton_gene.add_transcript(transcript, ref_trans_info_dict[ref_trans_id], gene_copy_num_dict, trans_copy_num_dict)


                logger.log(f"Liftoff: transcript_id\t: {transcript_id}\t{ref_trans_id}\n", debug=DEBUG)

                ###########################
                # Add LiftOn exons
                ###########################
                exons = l_feature_db.children(transcript, featuretype='exon')  # Replace 'exon' with the desired child feature type
                for exon in list(exons):
                    lifton_gene.add_exon(transcript_id, exon)

                logger.log(f"After adding exons\n", debug=DEBUG)
                ###########################
                # Add LiftOn CDS
                ###########################
                cdss = l_feature_db.children(transcript, featuretype='CDS')  # Replace 'exon' with the desired child feature type
                cds_num = 0
                for cds in list(cdss):
                    cds_num += 1
                    lifton_gene.add_cds(transcript_id, cds)
                logger.log(f"After adding CDSs\n", debug=DEBUG)
                miniprot_aln, has_valid_miniprot = lifton_utils.LiftOn_check_miniprot_alignment(gene.seqid, transcript, lifton_status, m_id_dict, m_feature_db, tree_dict, tgt_fai, ref_proteins, ref_trans_id)
                logger.log(f"After LiftOn_check_miniprot_alignment\n", debug=DEBUG)

                if (cds_num > 0):
                    #############################################
                    # Liftoff has protein
                    #############################################                
                    liftoff_aln = align.parasail_align("liftoff", l_feature_db, transcript, tgt_fai, ref_proteins, ref_trans_id, lifton_status)

                    logger.log(f"After liftoff parasail_align\n", debug=DEBUG)

                    if liftoff_aln is None:
                        #############################################
                        # There is no reference protein -> just keep Liftoff annotation
                        #############################################
                        logger.log("\t* Has CDS but no ref protein", debug=DEBUG)
                        lifton_status.annotation = "Liftoff_no_ref_protein"
                        lifton_status.status = ["no_ref_protein"]

                    elif liftoff_aln.identity < 1:
                        #############################################
                        # Liftoff annotation is not perfect
                        #############################################
                    
                        #############################################
                        # Running chaining algorithm if there are valid miniprot alignments
                        #############################################
                        if has_valid_miniprot:
                            logger.log("\t* Has CDS and valid miniprot", debug=DEBUG)
                            lifton_status.annotation = "LiftOn_chaining_algorithm" 
                            cds_list = fix_trans_annotation.chaining_algorithm(liftoff_aln, miniprot_aln, tgt_fai, fw)
                            lifton_gene.update_cds_list(transcript_id, cds_list)
                            logger.log("\tHas cds & protein & valid miniprot annotation!", debug=DEBUG)
                        else:
                            logger.log("\t* has CDS but invalid miniprot", debug=DEBUG)
                            lifton_status.annotation = "Liftoff_truncated"
                            logger.log("\tHas cds & protein & invalid miniprot annotation!", debug=DEBUG)
                            
                        #############################################
                        # Step 3.6.1.2: Check if there are mutations in the transcript
                        #############################################
                        lifton_trans_aln, lifton_aa_aln = lifton_gene.fix_truncated_protein(transcript_id, ref_trans_id, tgt_fai, ref_proteins, ref_trans, lifton_status)
                        

                            # for mutation in lifton_status.status:
                            #     if mutation != "synonymous" and mutation != "identical" and mutation != "nonsynonymous":
                            #         lifton_aa_aln.write_alignment(intermediate_dir, "lifton_AA", mutation, transcript_id)
                            #         lifton_trans_aln.write_alignment(intermediate_dir, "lifton_DNA", mutation, transcript_id)

                    elif liftoff_aln.identity == 1:
                        #############################################
                        # Step 3.6.2: Liftoff annotation is perfect
                        #############################################
                        lifton_status.annotation = "Liftoff_identical"
                        # SETTING LiftOn identity score => Same as Liftoff
                        lifton_status.lifton = liftoff_aln.identity
                        lifton_status.status = ["identical"]
                else:
                    #############################################
                    # Liftoff no protein
                    #############################################
                    if has_valid_miniprot:
                        ###########################
                        # Condition 3: LiftOn does not have proteins & miniprot has proteins
                        ###########################
                        logger.log("\t* No CDS; miniprot has ref protein", debug=DEBUG)

                        lifton_gene, transcript_id, lifton_status = run_miniprot.lifton_miniprot_with_ref_protein(transcript, m_feature_db, ref_gene_id, ref_trans_id, ref_gene_info_dict, ref_trans_info_dict, gene_copy_num_dict, trans_copy_num_dict, tgt_fai, ref_proteins, ref_trans, tree_dict, DEBUG)

                    else:
                        ###########################
                        # Condition 4: LiftOn does not have proteins & miniprot does not have proteins
                        ###########################
                        if (ref_trans_id in ref_proteins.keys()) :
                            logger.log("\t* No CDS & valid miniprot but have ref protein", debug=DEBUG)
                            lifton_status.annotation = "Liftoff_no_ref_protein"
                            lifton_status.status = ["no_ref_protein"]
                        
                        else:
                            logger.log("\t* No CDS & valid miniprot & no ref protein", debug=DEBUG)
                            lifton_status.annotation = "Liftoff_nc_transcript"
                            lifton_status.status = ["nc_transcript"]

                lifton_utils.write_lifton_status(fw_score, transcript_id, transcript, lifton_status)


                ###########################
                # Truncated reference proteins
                ###########################
                # if transcript_id_base in trunc_ref_proteins.keys() or transcript_id in trunc_ref_proteins.keys():
                #     print("Reference proteins is truncated.")
                #     lifton_status.annotation = "Reference_protein_truncated"
                #     lifton_status.status = ["truncated_ref_protein"]

                lifton_gene.add_lifton_status_attrs(transcript_id, lifton_status)
            ###########################
            # Writing out LiftOn entries
            ###########################
            lifton_gene.write_entry(fw)


    ################################
    # Step 6: Process miniprot transcripts
    ################################
    for mtrans in m_feature_db.features_of_type('mRNA'):
        mtrans_id = mtrans.attributes["ID"][0]
        mtrans_interval = Interval(mtrans.start, mtrans.end, mtrans_id)

        is_overlapped = lifton_utils.check_ovps_ratio(mtrans, mtrans_interval, args.overlap, tree_dict)

        if not is_overlapped:
            ref_trans_id = m_id_2_ref_id_trans_dict[mtrans_id]
            logger.log("miniprot: transcript_id\t: ", ref_trans_id, debug=DEBUG)

            if ref_trans_id in ref_trans_info_dict.keys():
                ref_gene_id = ref_trans_2_gene_dict[ref_trans_id] 
                
                # Mark transcript as lifted-over
                ref_features_dict[ref_gene_id][ref_trans_id] = True

                lifton_gene, new_m_trans_id, lifton_status = run_miniprot.lifton_miniprot_with_ref_protein(mtrans, m_feature_db, ref_gene_id, ref_trans_id, ref_gene_info_dict, ref_trans_info_dict, gene_copy_num_dict, trans_copy_num_dict, tgt_fai, ref_proteins, ref_trans, tree_dict, DEBUG)
            else:
                ref_gene_id = "gene-LiftOn"
                lifton_gene, new_m_trans_id, lifton_status = run_miniprot.lifton_miniprot_no_ref_protein(mtrans, m_feature_db, ref_gene_id, ref_trans_id, gene_copy_num_dict, trans_copy_num_dict, tree_dict, DEBUG)

            ###########################
            # Write scores for each transcript
            ###########################
            lifton_utils.write_lifton_status(fw_score, new_m_trans_id, mtrans, lifton_status)

            ###########################
            # Writing out LiftOn entries
            ###########################
            lifton_gene.add_lifton_status_attrs(new_m_trans_id, lifton_status)
            lifton_gene.write_entry(fw)

    stats.print_report(ref_features_dict, gene_copy_num_dict, trans_copy_num_dict, fw_unmapped, fw_extra_copy, debug=DEBUG)

    fw.close()
    fw_score.close()
    fw_unmapped.close()
    fw_extra_copy.close()

def main(arglist=None):
    args = parse_args(arglist)
    print("Run Lifton!!")
    run_all_lifton_steps(args)