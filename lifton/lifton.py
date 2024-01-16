from lifton import mapping, intervals, lifton_utils, annotation, extract_sequence, stats, logger, run_miniprot, run_liftoff, evaluation, __version__
from intervaltree import Interval
import argparse
from argparse import Namespace
from pyfaidx import Fasta, Faidx
import os
from concurrent.futures import ProcessPoolExecutor

def args_outgrp(parser):
    outgrp = parser.add_argument_group('* Output settings')
    outgrp.add_argument(
        '-o', '--output', default='lifton.gff3', metavar='FILE',
        help='write output to FILE in same format as input; by default, output is written to "lifton.gff3"'
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
    parser.add_argument('-E', '--evaluation', help='Run LiftOn in evaluation mode', action='store_true', default = False)
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
    ###################################
    # START for the LiftOn params
    ###################################
    referencegrp.add_argument(
        '-P', '--proteins', metavar='FASTA', required=False, default=None,
        help='the reference protein sequences.'
    )
    referencegrp.add_argument(
        '-T', '--transcripts', metavar='FASTA', required=False, default=None,
        help='the reference transcript sequences.'
    )
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
    # END for the LiftOn params
    ###################################
    parser._positionals.title = '* Required input (sequences)'
    parser._optionals.title = '* Miscellaneous settings'
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
    lifton_outdir = f"{outdir}/lifton_output/"
    args.directory = "intermediate_files/"
    intermediate_dir = f"{outdir}/lifton_output/{args.directory}"
    os.makedirs(intermediate_dir, exist_ok=True)
    args.directory = intermediate_dir
    logger.log(">> Reading target genome ...", debug=True)
    tgt_fai = Fasta(tgt_genome)
    logger.log(">> Reading reference genome ...", debug=True)
    ref_fai = Fasta(ref_genome)

    ################################
    # Step 1: Building database from the reference annotation
    ################################
    logger.log("\n>> Creating reference annotation database : ", args.reference_annotation, debug=True)
    ref_db = annotation.Annotation(args.reference_annotation, args.infer_genes)

    ################################
    # Step 2: Get all reference features to liftover
    ################################
    features = lifton_utils.get_parent_features_to_lift(args.features)
    ref_features_dict, ref_features_reverse_dict = lifton_utils.get_ref_liffover_features(features, ref_db)

    ################################
    # Step 3: Extract protein & DNA dictionaries from the selected reference features
    ################################
    ref_trans_file = args.transcripts    
    ref_proteins_file = args.proteins    
    if (ref_proteins_file is None) or (not os.path.exists(ref_proteins_file)) or (ref_trans_file is None) or (not os.path.exists(ref_trans_file)):
        logger.log(">> Creating transcript DNA dictionary from the reference annotation ...", debug=True)
        logger.log(">> Creating transcript protein dictionary from the reference annotation ...", debug=True)
        ref_trans, ref_proteins = extract_sequence.extract_features(ref_db, features, ref_fai)
        ref_proteins_file = lifton_utils.write_seq_2_file(intermediate_dir, ref_proteins, "proteins")
        ref_trans_file = lifton_utils.write_seq_2_file(intermediate_dir, ref_trans, "transcripts")
    else:
        logger.log(">> Reading transcript DNA dictionary from the reference fasta ...", debug=True)
        logger.log(">> Reading transcript protein dictionary from the reference fasta ...", debug=True)
        ref_trans = Fasta(ref_trans_file)
        ref_proteins = Fasta(ref_proteins_file)
    logger.log("\t * number of transcripts: ", len(ref_trans.keys()), debug=True)
    logger.log("\t * number of proteins: ", len(ref_proteins.keys()), debug=True)
    trunc_ref_proteins = lifton_utils.get_truncated_protein(ref_proteins)
    trunc_ref_proteins_file = lifton_utils.write_seq_2_file(intermediate_dir, trunc_ref_proteins, "truncated_proteins")
    logger.log("\t\t * number of truncated proteins: ", len(trunc_ref_proteins.keys()), debug=True)

    # # Evaluation mode
    # if args.evaluation:
    #     tgt_annotation = args.output
    #     ref_annotation = args.reference_annotation
    #     print("Run LiftOn in evaluation mode")
    #     print("lifton_outdir     : ", lifton_outdir)
    #     print("Ref genome        : ", ref_genome)
    #     print("Target genome     : ", tgt_genome)
    #     print("Ref annotation    : ", args.reference_annotation)
    #     print("Target annotation : ", args.output)
    #     print("ref_trans_file    : ", ref_trans_file)
    #     print("ref_proteins_file : ", ref_proteins_file)
    #     logger.log(">> Creating target database : ", tgt_annotation, debug=True)
    #     tgt_feature_db = annotation.Annotation(tgt_annotation, args.infer_genes).db_connection
    #     fw_score = open(lifton_outdir+"/eval.txt", "w")
    #     tree_dict = intervals.initialize_interval_tree(tgt_feature_db, features)
    #     for feature in features:
    #         for locus in tgt_feature_db.features_of_type(feature):#, limit=("chr1", 146652669, 146708545)):
    #             evaluation.tgt_evaluate(None, locus, ref_db.db_connection, tgt_feature_db, tree_dict, tgt_fai, ref_features_dict, ref_proteins, ref_trans, fw_score, DEBUG)
    #     fw_score.close()
    #     return

    # LiftOn mode
    ################################
    # Step 4: Run liftoff & miniprot
    ################################
    liftoff_annotation = lifton_utils.exec_liftoff(lifton_outdir, args)
    miniprot_annotation = lifton_utils.exec_miniprot(lifton_outdir, args, tgt_genome, ref_proteins_file)

    ################################
    # Step 5: Create liftoff and miniprot database
    ################################
    logger.log("\n>> Creating liftoff annotation database : ", liftoff_annotation, debug=True)
    l_feature_db = annotation.Annotation(liftoff_annotation, args.infer_genes).db_connection
    logger.log("\n>> Creating miniprot annotation database : ", miniprot_annotation, debug=True)
    m_feature_db = annotation.Annotation(miniprot_annotation, args.infer_genes).db_connection

    # Open output files
    fw = open(args.output, "w")
    fw_score = open(f"{lifton_outdir}/score.txt", "w")
    fw_chain = open(f"{lifton_outdir}/chain.txt", "w")
    fw_unmapped = open(f"{lifton_outdir}/unmapped_features.txt", "w")
    fw_extra_copy = open(f"{lifton_outdir}/extra_copy_features.txt", "w")

    ################################
    # Step 6: Creating miniprot 2 Liftoff ID mapping
    ################################
    ref_id_2_m_id_trans_dict, m_id_2_ref_id_trans_dict = mapping.miniprot_id_mapping(m_feature_db)

    ################################
    # Step 7: Initializing intervaltree
    ################################
    tree_dict = intervals.initialize_interval_tree(l_feature_db, features)

    ################################
    # Step 8: Process Liftoff genes & transcripts
    #     structure 1: gene -> transcript -> exon
    #     structure 2: transcript -> exon
    ################################
    for feature in features:
        for locus in l_feature_db.features_of_type(feature):#, limit=("NW_020825194.1", 28072487, 28072684)):
            lifton_gene = run_liftoff.process_liftoff(None, locus, ref_db.db_connection, l_feature_db, ref_id_2_m_id_trans_dict, m_feature_db, tree_dict, tgt_fai, ref_proteins, ref_trans, ref_features_dict, fw_score, fw_chain, DEBUG)
            # Writing out LiftOn entries
            lifton_gene.write_entry(fw)

    ################################
    # Step 9: Process miniprot transcripts
    ################################
    for mtrans in m_feature_db.features_of_type('mRNA'):
        mtrans_id = mtrans.attributes["ID"][0]
        mtrans_interval = Interval(mtrans.start, mtrans.end, mtrans_id)
        is_overlapped = lifton_utils.check_ovps_ratio(mtrans, mtrans_interval, args.overlap, tree_dict)
        if not is_overlapped:
            ref_trans_id = m_id_2_ref_id_trans_dict[mtrans_id]            
            # Link the reference trans ID to feature
            ref_gene_id, ref_trans_id = lifton_utils.get_ref_ids_miniprot(ref_features_reverse_dict, mtrans_id, m_id_2_ref_id_trans_dict)
            logger.log(f"miniprot: ref_gene_id: {ref_gene_id};  ref_trans_id: {ref_trans_id}\t: ", debug=DEBUG)
            if ref_trans_id in ref_proteins.keys() and ref_trans_id in ref_trans.keys():
                # The transcript match the reference transcript
                if ref_trans_id != None:
                    lifton_gene, transcript_id, lifton_status = run_miniprot.lifton_miniprot_with_ref_protein(mtrans, m_feature_db, ref_db.db_connection, ref_gene_id, ref_trans_id, tgt_fai, ref_proteins, ref_trans, tree_dict, ref_features_dict, DEBUG)
                else:
                    ref_gene_id = "LiftOn-gene"
                    lifton_gene, transcript_id, lifton_status = run_miniprot.lifton_miniprot_no_ref_protein(mtrans, m_feature_db, ref_gene_id, ref_trans_id, ref_features_dict, tree_dict, DEBUG)
            lifton_gene.add_lifton_status_attrs(transcript_id, lifton_status)
            lifton_utils.write_lifton_status(fw_score, transcript_id, mtrans, lifton_status)
            lifton_gene.write_entry(fw)

    ################################
    # Step 10: Printing stats
    ################################
    stats.print_report(ref_features_dict, fw_unmapped, fw_extra_copy, debug=DEBUG)
    # Close output files
    fw.close()
    fw_score.close()
    fw_chain.close()
    fw_unmapped.close()
    fw_extra_copy.close()


def main(arglist=None):
    banner = '''
====================================================================
An accurate homology lift-over tool between assemblies
====================================================================


    ██╗     ██╗███████╗████████╗ ██████╗ ███╗   ██╗
    ██║     ██║██╔════╝╚══██╔══╝██╔═══██╗████╗  ██║
    ██║     ██║█████╗     ██║   ██║   ██║██╔██╗ ██║
    ██║     ██║██╔══╝     ██║   ██║   ██║██║╚██╗██║
    ███████╗██║██║        ██║   ╚██████╔╝██║ ╚████║
    ╚══════╝╚═╝╚═╝        ╚═╝    ╚═════╝ ╚═╝  ╚═══╝
    '''
    print(banner)
    print(f"{__version__}\n")
    args = parse_args(arglist)
    run_all_lifton_steps(args)