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
        '-p', '--proteins', metavar='FASTA', required=False, default=None,
        help='the reference protein sequences.'
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
    ref_chroms = []

    ################################
    # Step 0: Reading target & reference genomes
    ################################
    tgt_genome = args.target
    ref_genome = args.reference
    outdir = os.path.dirname(args.output)
    os.makedirs(outdir, exist_ok=True)
    outdir = outdir if outdir is not "" else "."
    print(">> Reading target genome ...")
    tgt_fai = Fasta(tgt_genome)
    
    print(">> Reading reference genome ...")
    ref_fai = Fasta(ref_genome)
    print(args)

    ################################
    # Step 1: Building database from the reference annotation
    ################################
    feature_types = lifton_utils.get_feature_types(args.features)
    ref_db = annotation.Annotation(args.reference_annotation, args.infer_genes)
    # This is for special annotation file where there are no exon type
    # child_types = lifton_utils.get_child_types(feature_types, ref_db)

    ################################
    # Step 2: Creating protein & DNA dictionaries from the reference annotation
    ################################
    ref_proteins_file = args.proteins    
    if ref_proteins_file is None:
        print(">> Creating transcript protein dictionary from the reference annotation ...")
        ref_proteins = sequence.SequenceDict(ref_db, ref_fai, ['CDS', 'start_codon', 'stop_codon'], True)        
        ref_proteins_file = lifton_utils.write_protein_2_file(outdir, ref_proteins, False)

    else:
        print(">> Reading transcript protein dictionary from the reference fasta ...")
        ref_proteins = Fasta(ref_proteins_file)

    trunc_ref_proteins = lifton_utils.get_truncated_protein(ref_proteins)
    trunc_ref_proteins_file = lifton_utils.write_protein_2_file(outdir, trunc_ref_proteins, True)

    print(">> Creating transcript DNA dictionary from the reference annotation ...")
    # ref_trans = sequence.SequenceDict(ref_db, ref_fai, child_types, False)
    ref_trans = sequence.SequenceDict(ref_db, ref_fai, ['exon'], False)

    # print("; len(ref_proteins.keys()): ", len(ref_proteins.keys()))
    # print(";  len(ref_trans.keys()): ", len(ref_trans.keys()))

    ################################
    # Step 3: Run liftoff & miniprot
    ################################
    liftoff_annotation = lifton_utils.exec_liftoff(outdir, args)
    miniprot_annotation = lifton_utils.exec_miniprot(outdir, args, tgt_genome, ref_proteins_file)
    # print("liftoff_annotation : ", liftoff_annotation)
    # print("miniprot_annotation: ", miniprot_annotation)


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
    fw_truncated = open(args.output+".truncated", "w")
    fw_score = open(outdir+"/score.txt", "w")

    ################################
    # Step 4.1: Creating miniprot 2 Liftoff ID mapping
    ################################
    m_id_dict, aa_id_2_m_id_dict = mapping.id_mapping(m_feature_db)
    print("m_id_dict; ", m_id_dict)
    print("aa_id_2_m_id_dict; ", aa_id_2_m_id_dict)
    
    ################################
    # Step 4.2: Initializing intervaltree
    ################################
    tree_dict = intervals.initialize_interval_tree(l_feature_db)

    # Dictionary for extra copy
    gene_copy_num_dict = {}
    trans_copy_num_dict = {}
    gene_info_dict = {}
    trans_info_dict = {}
    trans_2_gene_dict = {}

    LIFTOFF_TOTAL_GENE_COUNT = 0
    LIFTOFF_TOTAL_TRANS_COUNT = 0
    LIFTOFF_BAD_PROT_TRANS_COUNT = 0
    LIFTOFF_GOOD_PROT_TRANS_COUNT = 0
    LIFTOFF_NC_TRANS_COUNT = 0
    LIFTOFF_OTHER_TRANS_COUNT = 0

    LIFTON_TOTAL_TRANS_COUNT = 0
    LIFTON_BAD_PROT_TRANS_COUNT = 0
    LIFTON_GOOD_PROT_TRANS_COUNT = 0
    LIFTON_NC_TRANS_COUNT = 0
    LIFTON_OTHER_TRANS_COUNT = 0
    LIFTON_MINIPROT_FIXED_GENE_COUNT = 0


    ################################
    # Step 4.3: Iterate gene entries & fixing CDS lists
    ################################
    # For missing transcripts.
    gene_copy_num_dict["gene-LiftOn"] = 0
    features = lifton_utils.get_parent_features_to_lift(args.features)
    
    # For transcripts without CDSs
    tmp_outdir = outdir + "/tmp/"
    os.makedirs(tmp_outdir, exist_ok=True)
    fw_other_trans = open(tmp_outdir+"/other_trans.txt", "w")
    fw_nc_trans = open(tmp_outdir+"/nc_trans.txt", "w")

    for feature in features:
        for gene in l_feature_db.features_of_type(feature):#, limit=("chr1", 0, 80478771)):
            LIFTOFF_TOTAL_GENE_COUNT += 1
            chromosome = gene.seqid
            gene_id = gene.attributes["ID"][0]
            gene_id_base = lifton_utils.get_ID_base(gene_id)

            ################################
            # Step 3.1: Creating gene copy number dictionary
            ################################
            lifton_utils.update_copy(gene_id_base, gene_copy_num_dict)

            ################################
            # Step 3.2: Creating LiftOn gene & gene_info
            ################################
            lifton_gene = lifton_class.Lifton_GENE(gene)
            # gene_info = copy.deepcopy(gene)
            # lifton_gene_info = lifton_class.Lifton_GENE_info(gene_info.attributes, gene_id_base)
            # gene_info_dict[gene_id_base] = lifton_gene_info
            
            ################################
            # Step 3.3: Adding LiftOn transcripts
            ################################
            # Assumption that all 1st level are transcripts
            transcripts = l_feature_db.children(gene, level=1)
            for transcript in list(transcripts):
                LIFTOFF_TOTAL_TRANS_COUNT += 1
                lifton_status = lifton_class.Lifton_Status()
                lifton_gene.add_transcript(transcript)
                transcript_id = transcript["ID"][0]
                transcript_id_base = lifton_utils.get_ID_base(transcript_id)

                ################################
                # Step 3.3.1: Creating trans copy number dictionary
                ################################
                lifton_utils.update_copy(transcript_id_base, trans_copy_num_dict)
                print("\ttranscript_id\t: ", transcript_id)
                if transcript_id != transcript_id_base:
                    print("&& transcript_id_base\t: ", transcript_id_base)
                # transcript_info = copy.deepcopy(transcript)
                # lifton_trans_info = lifton_class.Lifton_TRANS_info(transcript_info.attributes, transcript_id_base, gene_id_base)
                trans_2_gene_dict[transcript_id_base] = gene_id_base
                # trans_info_dict[transcript_id_base] = lifton_trans_info

                ###########################
                # Step 3.4: Adding exons
                ###########################
                exons = l_feature_db.children(transcript, featuretype='exon')  # Replace 'exon' with the desired child feature type
                for exon in list(exons):
                    lifton_gene.add_exon(transcript_id, exon)
                
                ###########################
                # Step 3.5: Adding CDS
                ###########################
                cdss = l_feature_db.children(transcript, featuretype='CDS')  # Replace 'exon' with the desired child feature type
                cds_num = 0
                for cds in list(cdss):
                    cds_num += 1
                    lifton_gene.add_cds(transcript_id, cds)


                #############################################
                # Step 3.6: Processing transcript
                #############################################
                if (cds_num > 0) and (transcript_id_base in ref_proteins.keys()):
                    liftoff_aln = align.parasail_align("liftoff", l_feature_db, transcript, tgt_fai, ref_proteins, transcript_id_base)

                    # SETTING Liftoff identity score
                    lifton_status.liftoff = liftoff_aln.identity

                    if liftoff_aln.identity < 1:
                        # This is for debugging purpose
                        # lifton_gene.remove_transcript(transcript_id)
                        
                        LIFTOFF_BAD_PROT_TRANS_COUNT += 1
                        #############################################
                        # Step 3.6.1: Liftoff annotation is not perfect
                        #############################################
                        lifton_status.annotation = "LiftOff_truncated"

                        # # Writing out truncated LiftOff annotation
                        # liftoff_aln.write_alignment(outdir, "liftoff", transcript_id)
                    
                        miniprot_aln, has_valid_miniprot = lifton_utils.LiftOn_check_miniprot_alignment(chromosome, transcript, lifton_status, m_id_dict, m_feature_db, tree_dict, tgt_fai, ref_proteins, transcript_id_base)


                        #############################################
                        # Step 3.6.1.1: Running chaining algorithm if there are valid miniprot alignments
                        #############################################
                        if has_valid_miniprot:
                            LIFTON_MINIPROT_FIXED_GENE_COUNT += 1
                            cds_list = fix_trans_annotation.chaining_algorithm(liftoff_aln, miniprot_aln, tgt_fai, fw)
                            lifton_gene.update_cds_list(transcript_id, cds_list)
                            lifton_status.annotation = "LiftOff_miniprot_chaining_algorithm" 
                        else:
                            print("Has cds & protein & miniprot annotation; but no valid miniprot annotation!")
                            
                        #############################################
                        # Step 3.6.1.2: Check if there are mutations in the transcript
                        #############################################
                        on_lifton_trans_aln, on_lifton_aa_aln = lifton_gene.fix_truncated_protein(transcript_id, transcript_id_base, tgt_fai, ref_proteins, ref_trans, lifton_status)
                        # SETTING LiftOn identity score
                        
                        if on_lifton_aa_aln.identity == 1:
                            LIFTON_GOOD_PROT_TRANS_COUNT += 1
                        elif on_lifton_aa_aln.identity < 1:
                            # Writing out truncated LiftOn annotation
                            LIFTON_BAD_PROT_TRANS_COUNT += 1
                            
                            for mutation in lifton_status.status:
                                if mutation != "synonymous" and mutation != "identical" and mutation != "nonsynonymous":
                                    on_lifton_aa_aln.write_alignment(outdir, "lifton_AA", mutation, transcript_id)
                                    on_lifton_trans_aln.write_alignment(outdir, "lifton_DNA", mutation, transcript_id)

                    elif liftoff_aln.identity == 1:
                        #############################################
                        # Step 3.6.2: Liftoff annotation is perfect
                        #############################################
                        LIFTOFF_GOOD_PROT_TRANS_COUNT += 1
                        LIFTON_GOOD_PROT_TRANS_COUNT += 1

                        miniprot_aln, has_valid_miniprot = lifton_utils.LiftOn_check_miniprot_alignment(chromosome, transcript, lifton_status, m_id_dict, m_feature_db, tree_dict, tgt_fai, ref_proteins, transcript_id_base)

                        # SETTING LiftOn identity score => Same as Liftoff
                        lifton_status.lifton = liftoff_aln.identity
                        lifton_status.annotation = "LiftOff_identical"
                        lifton_status.status = ["identical"]

                    final_status = ";".join(lifton_status.status)
                    fw_score.write(f"{transcript_id}\t{lifton_status.liftoff}\t{lifton_status.miniprot}\t{lifton_status.lifton}\t{lifton_status.annotation}\t{final_status}\t{transcript.seqid}:{transcript.start}-{transcript.end}\n")


                else:
                    # Only liftoff annotation
                    print("No cds || no protein")
                    if (cds_num > 0):
                        LIFTOFF_OTHER_TRANS_COUNT += 1
                        LIFTON_OTHER_TRANS_COUNT += 1
                        lifton_status.annotation = "LiftOff_no_ref_protein"
                        lifton_status.status = ["no_ref_protein"]
                        fw_other_trans.write(transcript_id+"\n")

                    else:
                        LIFTOFF_NC_TRANS_COUNT += 1
                        LIFTON_NC_TRANS_COUNT += 1
                        lifton_status.annotation = "LiftOff_nc_transcript"
                        lifton_status.status = ["nc_transcript"]
                        fw_nc_trans.write(transcript_id+"\n")



            ###########################
            # Step 4.7: Writing out LiftOn entries
            ###########################
            lifton_gene.write_entry(fw)
            LIFTON_TOTAL_TRANS_COUNT += 1
            # print("Final!!")
            # lifton_gene.print_gene()

            ###########################
            # Step 4.8: Adding LiftOn intervals
            ###########################
            gene_interval = Interval(lifton_gene.entry.start, lifton_gene.entry.end, gene_id)
            tree_dict[chromosome].add(gene_interval)


    ################################
    # Step 5: Finding extra copies
    ################################
    EXTRA_COPY_MINIPROT_COUNT = 0 
    NEW_LOCUS_MINIPROT_COUNT = 0
    # EXTRA_COPY_MINIPROT_COUNT, NEW_LOCUS_MINIPROT_COUNT = extra_copy.find_extra_copy(m_feature_db, tree_dict, aa_id_2_m_id_dict, gene_info_dict, trans_info_dict, trans_2_gene_dict, gene_copy_num_dict, trans_copy_num_dict, fw)

    print("Liftoff total gene loci\t\t\t: ", LIFTOFF_TOTAL_GENE_COUNT)
    print("Liftoff total transcript\t\t\t: ", LIFTOFF_TOTAL_TRANS_COUNT)
    print("Liftoff unperfect protein trans count\t\t\t: ", LIFTOFF_BAD_PROT_TRANS_COUNT)
    print("Liftoff good protein trans count\t\t\t: ", LIFTOFF_GOOD_PROT_TRANS_COUNT)
    print("Liftoff non-coding trans count\t\t\t: ", LIFTOFF_NC_TRANS_COUNT)
    print("Liftoff OTHER trans count\t\t\t: ", LIFTOFF_OTHER_TRANS_COUNT)
    print("\n\n")

    print("LiftOn total transcript\t\t\t: ", LIFTON_TOTAL_TRANS_COUNT)
    print("LiftOn unperfect protein trans count\t\t\t: ", LIFTON_BAD_PROT_TRANS_COUNT)
    print("LiftOn good protein trans count\t\t\t: ", LIFTON_GOOD_PROT_TRANS_COUNT)
    print("LiftOn non-coding trans count\t\t\t: ", LIFTON_NC_TRANS_COUNT)
    print("LiftOn OTHER trans count\t\t\t: ", LIFTON_OTHER_TRANS_COUNT)
    print("LiftOn miniprot fixed trans count\t\t\t: ", LIFTON_MINIPROT_FIXED_GENE_COUNT)

    fw.close()
    fw_truncated.close()
    fw_other_trans.close()
    fw_nc_trans.close()

def main(arglist=None):
    args = parse_args(arglist)
    print("Run Lifton!!")
    run_all_lifton_steps(args)