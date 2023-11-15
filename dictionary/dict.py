import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os
import gffutils
from Bio import SeqIO
import sys
from progress.bar import Bar

'''
One of the LiftOn modes is taking LiftOff and miniprot annotations along with protein sequences to generate LiftOn annotations.
Because some people prefer to directly download protein sequences from NCBI rather than using gffread to extract proteins, 
and NCBI assigns proteins with different IDs; for example, instead of starting with NM, they start with NP.
So, we need a Python function that can take LiftOff annotations and miniprot annotations and match miniprot proteins 
to Liftoff protein gene loci by creating a ID mapping dictionary.
The function should be able to handle various scenarios, including direct mapping with transcript IDs 
(applicable when protein sequences are extracted using gffread) and 
performing protein_ID to transcript_ID mapping (relevant for proteins directly downloaded from NCBI).
'''

def create_id_mapping(liftoff_annotations, miniprot_annotations):
    """
    Create an ID mapping dictionary between LiftOff and miniprot annotations.

    Parameters:
    - liftoff_annotations: Dictionary containing LiftOff annotations (gene_id to transcript_id mapping).
    - miniprot_annotations: Dictionary containing miniprot annotations (protein_id to gene_id mapping).

    Returns:
    - id_mapping: Dictionary mapping protein_ids to transcript_ids or gene_ids.
    """
    id_mapping = {}

    # direct mapping with transcript IDs
    for gene_id, transcript_id in liftoff_annotations.items():
        if transcript_id in miniprot_annotations:
            protein_id = miniprot_annotations[transcript_id]
            id_mapping[protein_id] = transcript_id

    # protein_ID to transcript_ID mapping
    for protein_id, gene_id in miniprot_annotations.items():
        if protein_id not in id_mapping:
            id_mapping[protein_id] = gene_id

    return id_mapping


def read_annotations_from_file(file_path):
    """
    Read annotations from a file and return a dictionary.

    Parameters:
    - file_path: Path to the annotation file.

    Returns:
    - annotations: Dictionary containing annotations.
    """

    annotations = {}
    with open(file_path, 'r') as file:
        for line in file:
            # assuming the file has two columns separated by a tab or space
            fields = line.strip().split('\t')  # change to '\t' or ' ' based on your file format
            if len(fields) == 2:
                annotations[fields[0]] = fields[1]
    return annotations


def output_mapping(id_mapping, output_file):
    """
    Write the ID mapping dictionary to a TSV file.

    Parameters:
    - id_mapping: Dictionary containing the ID mapping.
    - output_file: Path to the output TSV file.
    """
    with open(output_file, 'w') as file:
        for protein_id, mapped_id in id_mapping.items():
            file.write(f"{protein_id}\t{mapped_id}\n")


def preview_file(file_path, num_lines=5):
    """
    Preview the content of an annotation file.

    Parameters:
    - file_path: Path to the annotation file.
    - num_lines: Number of lines to preview (default is 5).
    """
    with open(file_path, 'r') as file:
        for _ in range(num_lines):
            line = file.readline().strip()
            if not line:
                break
            print(line)
    

def get_files_in_directory(directory_path):
    """
    Get a list of files in the specified directory.

    Parameters:
    - directory_path: Path to the directory.

    Returns:
    - file_list: List of file names in the directory.
    """
    file_list = []
    for file_name in os.listdir(directory_path):
        if os.path.isfile(os.path.join(directory_path, file_name)):
            file_list.append(file_name)
    return sorted(file_list)

def get_db(file_path, db_file=''):
    """
    Read data from a GFF file and return a list of features.

    Parameters:
    - file_path: Path to the GFF file.
    - db_file: Path to the SQLite database file (if none given, will create in a /db/ subdir)

    Returns:
    - db: The database generatd form gffutils
    - db_file: The path to the database

    """    
    # create an output for the database file
    base_name = os.path.splitext(os.path.basename(file_path))[0] + '.db'
    if db_file == '':
        db_dir = os.path.join(os.path.dirname(file_path), 'db/')
        db_name = base_name + '.db'
        os.makedirs(db_dir, exist_ok=True)
        db_file = db_dir+db_name
        print(f'[Info] No database name given, creating at: {db_file}')

    # get path of the input file
    fin = os.path.abspath(file_path)

    # generate database if empty (~5 mins / 3.42 mil features), 
    if not os.path.exists(db_file):
        try:
            db = gffutils.create_db(fin, db_file, merge_strategy="create_unique", force=True, \
                disable_infer_transcripts=True, disable_infer_genes=True, verbose=True)
        except:
            print('[Info] Error creating database, searching for line error...')
            find_problem_lines(file_path, f'./results/{db_name}_PROBLEM_LINES.out', quit=False) # DEBUGGING
    else: 
        db = gffutils.FeatureDB(db_file)

    return db, db_file


def find_problem_lines(gff_file, outfile='', quit=True):
    '''
    Finds the line in the GFF file causing error for gffutils.
    
    Parameters:
    - gff_file: filepath to the GFF file
    - output: output file with error, only if given
    - quit: if true, will stop running program when problem lines found 

    '''
    f = open(gff_file, 'r')
    lines = f.readlines()
    problem_lines = []
    pbar = Bar('Searching GFF file... ', max=len(lines))
    for i in range(len(lines)):
        line = lines[i]
        if line[0] != "#":
            try:
                gffutils.create_db(line, ":memory:", from_string=True, force=True)
            except:
                problem_lines.append(i+1)

        pbar.next()
    pbar.finish()

    if len(problem_lines) != 0:
        if outfile != '':
            with (open(outfile), 'w') as f:
                for l in problem_lines:
                    f.write(str(l) + '\n')
        if quit:
            raise ValueError("ERROR: Incorrect GFF/GTF syntax on line(s) " + str(list(problem_lines)))


''' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ MAIN DRIVER FUNCTION ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'''
def run():

    # get the annotation files in a paired list
    liftoff_dir = './annotations/liftoff/'
    miniprot_dir = './annotations/miniprot/'
    liftoff_files = get_files_in_directory(liftoff_dir)
    miniprot_files = get_files_in_directory(miniprot_dir)
    # paired_list = list(zip(liftoff_files, miniprot_files))
    paired_list = [(liftoff_dir+lf, miniprot_dir+mf) for lf, mf in zip(liftoff_files, miniprot_files)]
    # print(paired_list)

    for liftoff_f, miniprot_f in paired_list[:1]:
        print(f'previewing: {liftoff_f}, {miniprot_f}')
        preview_file(liftoff_f)
        print('-'*170)
        preview_file(miniprot_f)
        liftoff_db, _ = get_db(liftoff_f)
        miniprot_db, _ = get_db(miniprot_f)

    # preview files
    # for file in liftoff_files:
    #     preview_file()

    # read in the files
    # liftoff_file = './annotations/liftoff/'
    # miniprot_file = './annotations/miniprot/'
    
    # liftoff_annotations = read_annotations_from_file(liftoff_file)
    # miniprot_annotations = read_annotations_from_file(miniprot_file)

    # # create mapping
    # mapping = create_id_mapping(liftoff_annotations, miniprot_annotations)

    # # output the results in a tsv file
    # output_file = './results/output.tsv'
    # output_mapping(mapping, output_file)

    return

if __name__ == '__main__':
    run()