from lifton import intervals, lifton_utils, annotation, extract_sequence, stats, logger, run_liftoff, run_miniprot, run_evaluation, __version__
from intervaltree import Interval
import argparse
from pyfaidx import Fasta
import os, sys
import time

def args_gffutils(parser):
    gffutils_grp = parser.add_argument_group('* gffutils parameters')
    gffutils_grp.add_argument(
        '--merge-strategy', default='create_unique', choices=['create_unique', 'merge', 'error', 'warning', 'replace'],
        help='Strategy for merging features when building the database. Default is "create_unique".'
    )
    gffutils_grp.add_argument(
        '--id-spec', default=None,
        help='Attribute to use as feature ID. Default is "ID".'
    )
    gffutils_grp.add_argument(
        '--force', default=False, action='store_true', help='Overwrite existing database.'
    )
    gffutils_grp.add_argument(
        '--verbose', default=False, action='store_true', help='Enable verbose output.'
    )
    return gffutils_grp


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
    aligngrp.add_argument('-mp_options', metavar='=STR', type=str, default='',
                          help='space delimited miniprot parameters.')
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
        '-min_miniprot', default=0.9, metavar='MIN_MINIPROT', type=float,
        help='The minimum length ratio of a protein-coding transcript to the longest protein-coding transcript within a gene locus, as identified exclusively by miniprot in the target genome, is set by default to MIN_MINIPROT=0.9.'
    )
    aligngrp.add_argument(
        '-max_miniprot', default=1.5, metavar='MAX_MINIPROT', type=float,
        help='The maximum length ratio of a protein-coding transcript to the longest protein-coding transcript within a gene locus, as identified exclusively by miniprot in the target genome, is set by default to MIN_MINIPROT=1.5.'
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
    parser.add_argument(
        '-time', '--measure_time', required=False, action='store_true',
        help='Enable time measurement for each step'
    )


def parse_args(arglist):
    parser = argparse.ArgumentParser(description='Lift features from one genome assembly to another')
    parser.add_argument('target', help='target fasta genome to lift genes to')
    parser.add_argument('reference', help='reference fasta genome to lift genes from')
    parser.add_argument('-E', '--evaluation', help='Run LiftOn in evaluation mode', action='store_true', default = False)
    parser.add_argument('-EL', '--evaluation-liftoff-chm13', help='Run LiftOn in evaluation mode', action='store_true', default = False)
    parser.add_argument('-c', '--write_chains', help='Write chaining files', action='store_true', default = True)
    parser.add_argument('--no-orf-search', help='Do not perform open reading frame search', action='store_true', default = False)
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
    referenceseqgrp = parser.add_argument_group('* Optional input (Reference sequences)')
    referenceseqgrp.add_argument(
        '-P', '--proteins', metavar='FASTA', required=False, default=None,
        help='the reference protein sequences. If not provided, the protein sequences will be extracted from the reference annotation. The ID of the protein sequences must match the ID of the transcript sequences.'
    )
    referenceseqgrp.add_argument(
        '-T', '--transcripts', metavar='FASTA', required=False, default=None,
        help='the reference transcript sequences. If not provided, the transcript sequences will be extracted from the reference annotation.'
    )
    liftoffrefrgrp = parser.add_argument_group('* Optional input (Liftoff annotation)')
    liftoffrefrgrp.add_argument(
        '-L', '--liftoff', metavar='gff', default=None,
        help='the annotation generated by Liftoff (or) '
                'name of Liftoff gffutils database; if not specified, '
                'a Liftoff database will be built automatically'
    )
    miniprotrefrgrp = parser.add_argument_group('* Optional input (miniprot annotation)')
    miniprotrefrgrp.add_argument(
        '-M', '--miniprot', metavar='gff', default=None,
        help='the annotation generated by miniprot (or) '
                'name of miniprot gffutils database; if not specified, '
                'a miniprot database will be built automatically'
    )
    parser_gffutils_grp = args_gffutils(parser)
    parser.add_argument('-ad', '--annotation-database', metavar='SOURCE', help='The source of the reference annotation (RefSeq / GENCODE / others).', default = "RefSeq")
    ###################################
    # END for the LiftOn params
    ###################################
    parser._positionals.title = '* Required input (sequences)'
    parser._optionals.title = '* Miscellaneous settings'
    parser._action_groups = [parser._positionals, referencegrp, referenceseqgrp, 
                             liftoffrefrgrp, miniprotrefrgrp, parser_gffutils_grp, 
                             parser_outgrp, parser._optionals, parser_aligngrp]
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
        args.mm2_options += " --end-bonus 5"
    if (float(args.s) > float(args.sc)):
        parser.error("-sc must be greater than or equal to -s")
    if (args.chroms is None and args.unplaced is not None):
        parser.error("-unplaced must be used with -chroms")    
    return args
    

def run_all_lifton_steps(args):
    t1 = time.process_time()
    ################################
    # Step 0: Reading target & reference genomes
    ################################
    tgt_genome = args.target
    ref_genome = args.reference
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
    stats_dir= f"{outdir}/lifton_output/stats/"
    os.makedirs(stats_dir, exist_ok=True)
    args.directory = intermediate_dir
    logger.log(">> Reading target genome ...", debug=True)
    tgt_fai = Fasta(tgt_genome)
    logger.log(">> Reading reference genome ...", debug=True)
    ref_fai = Fasta(ref_genome)    

    t2 = time.process_time()
    ################################
    # Step 1: Building database from the reference annotation
    ################################
    logger.log("\n>> Creating reference annotation database : ", args.reference_annotation, debug=True)
    ref_db = annotation.Annotation(args.reference_annotation, args.infer_genes, args.infer_transcripts, args.merge_strategy, args.id_spec, args.force, args.verbose)

    t3 = time.process_time()
    ################################
    # Step 2: Get all reference features to liftover
    ################################
    features = lifton_utils.get_parent_features_to_lift(args.features)
    ref_features_dict, ref_features_len_dict, ref_features_reverse_dict, ref_trans_exon_num_dict = lifton_utils.get_ref_liffover_features(features, ref_db, intermediate_dir, args)

    t4 = time.process_time()
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
    logger.log("\t\t * number of truncated proteins: ", len(trunc_ref_proteins.keys()), debug=True)

    ################################
    # optional Step: Evaluation mode
    ################################
    if args.evaluation:
        tgt_annotation = args.output
        ref_annotation = args.reference_annotation
        print("Run LiftOn in evaluation mode")
        print("lifton_outdir     : ", lifton_outdir)
        print("Ref genome        : ", ref_genome)
        print("Target genome     : ", tgt_genome)
        print("Ref annotation    : ", args.reference_annotation)
        print("Target annotation : ", args.output)
        print("ref_trans_file    : ", ref_trans_file)
        print("ref_proteins_file : ", ref_proteins_file)
        logger.log(">> Creating target database : ", tgt_annotation, debug=True)
        os.makedirs(lifton_outdir, exist_ok=True)
        tgt_feature_db = annotation.Annotation(tgt_annotation, args.infer_genes).db_connection
        fw_score = open(lifton_outdir+"/eval.txt", "w")
        tree_dict = {}
        processed_features = 0
        for feature in features:
            for locus in tgt_feature_db.features_of_type(feature):#, limit=("chr1", 146652669, 146708545)):
                lifton_gene = run_evaluation.evaluation(None, locus, ref_db.db_connection, tgt_feature_db, tree_dict, tgt_fai, ref_proteins, ref_trans, ref_features_dict, fw_score, args, ENTRY_FEATURE=True)
                if processed_features % 20 == 0:
                    sys.stdout.write("\r>> LiftOn evaluated: %i features." % processed_features)
                processed_features += 1
        fw_score.close()
        return
    
    ################################
    # Step 4: Run liftoff & miniprot
    ################################
    t5 = time.process_time()
    liftoff_annotation = lifton_utils.exec_liftoff(lifton_outdir, args)
    t6 = time.process_time()
    miniprot_annotation = lifton_utils.exec_miniprot(lifton_outdir, args, tgt_genome, ref_proteins_file)

    t7 = time.process_time()
    ################################
    # Step 5: Create liftoff and miniprot database
    ################################
    logger.log(f"\n>> Creating liftoff annotation database : {liftoff_annotation}", debug=True)
    l_feature_db = annotation.Annotation(liftoff_annotation, args.infer_genes, args.infer_transcripts, args.merge_strategy, args.id_spec, args.force, args.verbose).db_connection
    t8 = time.process_time()
    logger.log(f">> Creating miniprot annotation database : {miniprot_annotation}", debug=True)
    m_feature_db = annotation.Annotation(miniprot_annotation, args.infer_genes, args.infer_transcripts, args.merge_strategy, args.id_spec, args.force, args.verbose).db_connection
    fw = open(args.output, "w")
    fw_score = open(f"{lifton_outdir}/score.txt", "w")
    fw_unmapped = open(f"{stats_dir}/unmapped_features.txt", "w")
    fw_extra_copy = open(f"{stats_dir}/extra_copy_features.txt", "w")
    fw_mapped_feature = open(f'{stats_dir}/mapped_feature.txt', 'w')
    fw_mapped_trans = open(f'{stats_dir}/mapped_transcript.txt', 'w')
    fw_chain = open(f"{lifton_outdir}/chain.txt", "w") if args.write_chains else None

    t9 = time.process_time()
    ################################
    # Step 6: Creating miniprot 2 Liftoff ID mapping & Initializing intervaltree
    ################################
    ref_id_2_m_id_trans_dict, m_id_2_ref_id_trans_dict = lifton_utils.miniprot_id_mapping(m_feature_db)
    tree_dict = intervals.initialize_interval_tree(l_feature_db, features)
    transcripts_stats_dict = {'coding': {}, 'non-coding': {}, 'other': {}}
    processed_features = 0
    
    t10 = time.process_time()
    ################################
    # Step 7: Process Liftoff genes & transcripts
    #     structure 1: gene -> transcript -> exon
    #     structure 2: transcript -> exon
    ################################
    for feature in features:#CP132235.1:34100723-34103135
        for locus in l_feature_db.features_of_type(feature):#, limit=("chr1", 206759465, 206890507)):
            lifton_gene = run_liftoff.process_liftoff(None, locus, ref_db.db_connection, l_feature_db, ref_id_2_m_id_trans_dict, m_feature_db, tree_dict, tgt_fai, ref_proteins, ref_trans, ref_features_dict, fw_score, fw_chain, args, ENTRY_FEATURE=True)
            if lifton_gene is None or lifton_gene.ref_gene_id is None:
                continue
            lifton_gene.write_entry(fw,transcripts_stats_dict)
            if processed_features % 20 == 0:
                sys.stdout.write("\r>> LiftOn processed: %i features." % processed_features)
            processed_features += 1

    t11 = time.process_time()
    ################################
    # Step 8: Process miniprot transcripts
    ################################
    for mtrans in m_feature_db.features_of_type('mRNA'):
        lifton_gene = run_miniprot.process_miniprot(mtrans, ref_db, m_feature_db, tree_dict, tgt_fai, ref_proteins, ref_trans, ref_features_dict, fw_score, m_id_2_ref_id_trans_dict, ref_features_len_dict, ref_trans_exon_num_dict, ref_features_reverse_dict, args)
        if lifton_gene is None or lifton_gene.ref_gene_id is None:
            continue
        lifton_gene.write_entry(fw, transcripts_stats_dict)
        if processed_features % 20 == 0:
            sys.stdout.write("\r>> LiftOn processed: %i features." % processed_features)
        processed_features += 1
    
    t12 = time.process_time()
    ################################
    # Step 9: Printing stats
    ################################
    stats.print_report(ref_features_dict, transcripts_stats_dict, fw_unmapped, fw_extra_copy, fw_mapped_feature, fw_mapped_trans, debug=args.debug)
    fw.close()
    fw_score.close()
    fw_unmapped.close()
    fw_extra_copy.close()
    fw_mapped_feature.close()
    fw_mapped_trans.close()
    if args.write_chains: fw_chain.close()
    t13 = time.process_time()

    if args.measure_time:
        reading_target_reference_genomes= t2 - t1
        creating_reference_annotation_database = t3 - t2
        get_all_reference_features_to_lift = t4 - t3
        extract_protein_dna_dictionaries = t5 - t4
        run_liftoff_miniprot = t6 - t5
        create_liftoff_database = t8 - t6
        create_miniprot_database = t9 - t8
        miniprot_2_liftoff_id_mapping = t10 - t9
        process_liftoff_genes_transcripts = t11 - t10
        process_miniprot_transcripts = t12 - t11
        report_stats = t13 - t12
        overall_time = t13 - t1
        print("Time taken for each step:")
        print(f"Reading target & reference genomes: {reading_target_reference_genomes}")
        print(f"Creating reference annotation database: {creating_reference_annotation_database}")
        print(f"Get all reference features to liftover: {get_all_reference_features_to_lift}")
        print(f"Extract protein & DNA dictionaries: {extract_protein_dna_dictionaries}")
        print(f"Run liftoff & miniprot: {run_liftoff_miniprot}")
        print(f"Create liftoff database: {create_liftoff_database}")
        print(f"Create miniprot database: {create_miniprot_database}")
        print(f"Miniprot 2 Liftoff ID mapping: {miniprot_2_liftoff_id_mapping}")
        print(f"Process Liftoff genes & transcripts: {process_liftoff_genes_transcripts}")
        print(f"Process miniprot transcripts: {process_miniprot_transcripts}")
        print(f"Report stats: {report_stats}")
        print(f"Overall time: {overall_time}")
        fw_time = open(f"{outdir}/time.txt", "w")
        fw_time.write(f"{reading_target_reference_genomes}\tReading target & reference genomes\n")
        fw_time.write(f"{creating_reference_annotation_database}\tCreating reference annotation database\n")
        fw_time.write(f"{get_all_reference_features_to_lift}\tGet all reference features to liftover\n")
        fw_time.write(f"{extract_protein_dna_dictionaries}\tExtract protein & DNA dictionaries\n")
        fw_time.write(f"{run_liftoff_miniprot}\tRun liftoff & miniprot\n")
        fw_time.write(f"{create_liftoff_database}\tCreate liftoff database\n")
        fw_time.write(f"{create_miniprot_database}\tCreate miniprot database\n")
        fw_time.write(f"{miniprot_2_liftoff_id_mapping}\tMiniprot 2 Liftoff ID mapping\n")
        fw_time.write(f"{process_liftoff_genes_transcripts}\tProcess Liftoff genes & transcripts\n")
        fw_time.write(f"{process_miniprot_transcripts}\tProcess miniprot transcripts\n")
        fw_time.write(f"{report_stats}\tReport stats\n")
        fw_time.write(f"{overall_time}\tOverall time\n")
        fw_time.close()


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
    print(banner, file=sys.stderr)
    args = parse_args(arglist)
    if not run_miniprot.check_miniprot_installed(): 
        sys.exit("miniprot is not installed. Please install miniprot before running LiftOn.")
    run_all_lifton_steps(args)
