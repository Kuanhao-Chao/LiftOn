import os
import shutil

CLSUTER_OUTPUTS = { 'clusters': 'clusters', 'unmapped': 'unmapped_closest_paralogs'}

MMSEQS_INTERMEDIATES = {'dir': "mmseqs_intermediates", 'ref_coding_prefix': 'reference_coding',
                        'target_coding_prefix':'target_coding', 'ref_noncoding_prefix': 'ref_noncoding',
                        'target_noncoding_prefix': 'target_noncoding'}

SYNTENY_OUTPUTS = {'gene_order': 'gene_order', 'plot': 'gene_order_plot.pdf' }

VARIANTS_OUTPUTS = { 'alignments': 'alignments.fa', 'variants': 'variant_effects'}


def build_filepath(path_parts):
    path = ""
    for part in path_parts:
        part.rstrip("/")
        path += part + "/"
    return path[:-1]


def make_directory(path):
    if os.path.exists(path) is False:
        os.mkdir(path)


def remove_directory(path):
    if os.path.exists(path):
        shutil.rmtree(path)


def make_file(file_name, force):
    return os.path.exists(file_name) is False or force


def remove_file(file_name):
    if os.path.exists(file_name):
        os.remove(file_name)