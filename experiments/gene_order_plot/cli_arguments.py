import argparse

def parse_args(arglist):
    parser = argparse.ArgumentParser(description='Compare gene annotations across genome assemblies')
    parser.add_argument('subcommand', choices=['clusters', 'variants', 'synteny', 'all'])
    parser.add_argument('-r', help='reference fasta', required=True)
    parser.add_argument('-t', help='target fasta', required=True)
    parser.add_argument('-rg', metavar='GFF/GTF or DB', help='reference annotation file to lift over in GFF or GTF '
                                                            'format or gffutils database created in previous liftoff '
                                                             'or liftofftools run', required=True)
    parser.add_argument('-tg', metavar='GFF/GTF or DB', help=' target annotation file to lift over in GFF or GTF '
                                                             'format or gffutils databased created in previous '
                                                             'liftoff or liftofftools run', required=True)
    parser.add_argument('-c', action='store_true', default=False, help='analyze protein coding gene clusters only',
                        required=False)
    parser.add_argument('-f', required=False, help='text file with additional feature types besides genes to analyze')
    parser.add_argument('-infer-genes', required=False, action='store_true', default=False)
    parser.add_argument('-dir', default="liftofftools_output", required=False, help="output directory")
    parser.add_argument('-force', required=False, default=False, action='store_true', help="force overwrite of "
                                                                                           "output/intermediate files in -dir")
    clusters_group = parser.add_argument_group("clusters arguments")
    clusters_group.add_argument('-mmseqs_path', help="mmseqs path if not in working directory or PATH", required=False)
    clusters_group.add_argument('-mmseqs_params', default="--min-seq-id 0.9 -c 0.9", metavar='=STR',
                        required=False,help='space delimited list of additional mmseqs parameters. Default="--min-seq-id 0.9 -c 0.9"')
    synteny_group = parser.add_argument_group("synteny arguments")
    synteny_group.add_argument('-edit-distance', help="calculate edit distance between reference gene order and "
                                                      "target gene order", required=False, action='store_true')
    synteny_group.add_argument('-r-sort', help="txt file with the order of the reference chromosomes to be plotted on the x-axis", required=False, default=None)
    synteny_group.add_argument('-t-sort', help="txt file with the order of the target chromosomes to be plotted on the y-axis", required=False, default=None)                                    
    parser.add_argument('-V', '--version', help='show program version', action='version', version='v0.4.3')
    parser._positionals.title = 'Subcommands'
    args = parser.parse_args(arglist)
    if bool(args.r_sort) ^ bool(args.t_sort):
        parser.error('-r-sort and -t-sort must be given together')
    return args
