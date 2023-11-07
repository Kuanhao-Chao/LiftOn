from gffutils import FeatureDB
import pybedtools
import os, re
import json

import sys
from gffutils import FeatureDB
from Bio import SeqIO
import parasail

def read_fasta(file_path):
    sequences = {}
    current_sequence = ""
    current_header = ""

    with open(file_path, 'r') as fasta_file:
        for line in fasta_file:
            line = line.strip()
            if line.startswith('>'):
                if current_header != "":
                    sequences[current_header] = current_sequence
                current_header = line[1:]
                current_sequence = ""
            else:
                current_sequence += line

        if current_header != "":
            sequences[current_header] = current_sequence

    return sequences

def get_id_fraction(reference, target):
    matches = 0
    for i, letter in enumerate(reference):
        if letter == target[i]:
            matches += 1
        if target[i] == "*" or target[i] == ".":
            break
    # print("matches: ", matches)
    # print("len(reference): ", len(reference))
    # print("len(target): ", len(target))
    return matches, max(len(reference), len(target))



TARGET = "CHM13_MANE"

genome = "/ccb/salz3/kh.chao/ref_genome/homo_sapiens/T2T-CHM13/chm13v2.0.fa"

output_dir = f"/ccb/salz2/kh.chao/Lifton/results/{TARGET}/eviann_lifton_cmp/"
protein_fa = f"/ccb/salz2/kh.chao/PR_liftoff_protein_search/data/MANE_RefSeq/MANE.GRCh38.v1.2.refseq_genomic.AA.cleaned.fa"

os.makedirs(output_dir, exist_ok=True)


evionn_protein_fa = f"/ccb/salz2/kh.chao/Lifton/data/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna_AA.fa"

lifton_protein_fa = f"/ccb/salz2/kh.chao/Lifton/results/{TARGET}/{TARGET}_lifton_protein.fa"
lifton_protein_gff = f"/ccb/salz2/kh.chao/Lifton/results/{TARGET}/{TARGET}_lifton.gff3" 

gene_loci_dict_fn = "/ccb/salz2/kh.chao/Lifton/results/CHM13_MANE/eviann_lifton_cmp/gene_loci.json"


print("protein_fa: ", protein_fa)
print("evionn_protein_fa: ", evionn_protein_fa)
print("lifton_protein_fa: ", lifton_protein_fa)
print("lifton_protein_gff: ", lifton_protein_gff)

with open(gene_loci_dict_fn, 'r') as openfile:
    gene_loci_dict = json.load(openfile)
print("gene_loci_dict: ", gene_loci_dict)


evionn_sequences = read_fasta(evionn_protein_fa)
lifton_sequences = read_fasta(lifton_protein_fa)
MANE_sequences = read_fasta(protein_fa)

both_cnt = 0
evionn_only_cnt = 0
lifton_only_cnt = 0
miss_both_cnt = 0

fwname = output_dir + "identities.txt"
fwname_both = output_dir + "both.txt"
fwname_evionn_only = output_dir + "eviann_only.txt"
fwname_lifton_only = output_dir + "lifton_only.txt"
fwname_miss_both = output_dir + "miss_both.txt"


fw = open(fwname, "w")
fw_both = open(fwname_both, "w")
fw_evionn_only = open(fwname_evionn_only, "w")
fw_lifton_only = open(fwname_lifton_only, "w")
fw_miss_both = open(fwname_miss_both, "w")

EARLY_STOP_IN_REFERENCE = 0

print("lifton_sequences.keys(): ", lifton_sequences.keys())
for id, sequence in MANE_sequences.items():
    # print(id)
    process = id in gene_loci_dict.keys() and id in lifton_sequences.keys()
    print(process)
    if process:
        if len(gene_loci_dict[id]) == 0:
            lifton_only_cnt += 1
            fw_lifton_only.write(id + "\n")

        else:

            both_cnt += 1
            fw_both.write(id + "\n")
            both_cnt += 1
            
            matrix = parasail.Matrix("blosum62")
            gap_open = 11
            gap_extend = 1

            reference_seq = MANE_sequences[id]

            # if "*" in reference_seq or "." in reference_seq:
            #     EARLY_STOP_IN_REFERENCE += 1
            #     continue
            
            max_evionn_identity = 0
            for eviann_id in gene_loci_dict[id]:
                # Mapping eviann
                eviann_seq = evionn_sequences[eviann_id]
                
                eviann_parasail_res = parasail.nw_trace_scan_sat(eviann_seq, reference_seq, gap_open, gap_extend, matrix)
                eviann_matches, eviann_length = get_id_fraction(eviann_parasail_res.traceback.ref, eviann_parasail_res.traceback.query)
                max_evionn_identity = max(max_evionn_identity, eviann_matches/eviann_length)

            max_lifton_identity = 0
            # Mapping lifton
            lifton_seq = lifton_sequences[id]            
            lifton_parasail_res = parasail.nw_trace_scan_sat(lifton_seq, reference_seq, gap_open, gap_extend, matrix)

            lifton_matches, lifton_length = get_id_fraction(lifton_parasail_res.traceback.ref, lifton_parasail_res.traceback.query)
            max_lifton_identity = max(max_lifton_identity, lifton_matches/lifton_length)
            # #print(matches, length)
            fw.write(f"{id}\t{max_evionn_identity}\t{max_lifton_identity}\n")

    else:
        fw_miss_both.write(id + "\n")
        miss_both_cnt += 1

print("lifton_only_cnt: ", lifton_only_cnt)
print("miss_both_cnt: ", miss_both_cnt)
print("both_cnt: ", both_cnt)

print("EARLY_STOP_IN_REFERENCE: ", EARLY_STOP_IN_REFERENCE)