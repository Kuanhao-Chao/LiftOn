from gffutils import FeatureDB
import pybedtools
import os, re

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
    return matches, max(len(reference), len(target))

TARGET = sys.argv[1]

if TARGET == "CHM13_MANE" or TARGET == "CHM13_GRCh38" or TARGET == "GRCh38_RefSeq" or TARGET == "Mus_musculus":
    print("Running with ", TARGET)
    genome = ""
    if TARGET == "GRCh38_RefSeq":
        genome = "/ccb/salz3/kh.chao/ref_genome/homo_sapiens/NCBI_Refseq_chr_fixed/GCF_000001405.40_GRCh38.p14_genomic.fna"
    elif TARGET == "CHM13_MANE" or TARGET == "CHM13_GRCh38":
        genome = "/ccb/salz3/kh.chao/ref_genome/homo_sapiens/T2T-CHM13/chm13v2.0.fa"
    elif TARGET == "Mus_musculus":
        genome = "/ccb/salz3/kh.chao/ref_genome/mus_musculus/NCBI_Refseq/GCF_000001635.27_GRCm39_genomic.fna"
else:
    sys.exit(-1)



output_dir = f"/ccb/salz2/kh.chao/Lifton/results/CHM13_MANE/liftoff_lifton_cmp/"

os.makedirs(output_dir, exist_ok=True)

MANE_protein_fa = f"/ccb/salz2/kh.chao/PR_liftoff_protein_search/data/MANE_RefSeq/MANE.GRCh38.v1.2.refseq_genomic.AA.cleaned.fa"

liftoff_protein_fa = f"/ccb/salz2/kh.chao/PR_liftoff_protein_search/results/liftoff/{TARGET}/{TARGET}_AA.fa"

lifton_protein_fa = f"/ccb/salz2/kh.chao/Lifton/results/{TARGET}/{TARGET}_protein.fa"

lifton_protein_gff = f"/ccb/salz2/kh.chao/Lifton/results/{TARGET}/lifton.gff3" 


print("MANE_protein_fa: ", MANE_protein_fa)
print("liftoff_protein_fa: ", liftoff_protein_fa)
print("lifton_protein_fa: ", lifton_protein_fa)
print("lifton_protein_gff: ", lifton_protein_gff)

liftoff_sequences = read_fasta(liftoff_protein_fa)
lifton_sequences = read_fasta(lifton_protein_fa)
MANE_sequences = read_fasta(MANE_protein_fa)

both_cnt = 0
liftoff_only_cnt = 0
lifton_only_cnt = 0
miss_both_cnt = 0

fwname = output_dir + "identities.txt"
fwname_both = output_dir + "both.txt"
fwname_liftoff_only = output_dir + "liftoff_only.txt"
fwname_lifton_only = output_dir + "lifton_only.txt"
fwname_miss_both = output_dir + "miss_both.txt"

fw = open(fwname, "w")
fw_both = open(fwname_both, "w")
fw_liftoff_only = open(fwname_liftoff_only, "w")
fw_lifton_only = open(fwname_lifton_only, "w")
fw_miss_both = open(fwname_miss_both, "w")

# print("lifton_sequences.keys(): ", lifton_sequences.keys())
for id, sequence in MANE_sequences.items():

    in_liftoff = id in liftoff_sequences.keys()
    in_lifton = id in lifton_sequences.keys()
    
    # print([i for i in lifton_ids[id] if i in lifton_sequences.keys()])
    # in_lifton_ls = [i for i in lifton_ids[id] if i in lifton_sequences.keys()]
    # in_lifton = len(in_lifton_ls) > 0

    if in_liftoff and in_lifton:
        fw_both.write(id + "\n")
        both_cnt += 1

        matrix = parasail.Matrix("blosum62")
        gap_open = 11
        gap_extend = 1

        reference_seq = MANE_sequences[id]

        # Mapping liftoff
        liftoff_seq = liftoff_sequences[id]
        
        liftoff_parasail_res = parasail.nw_trace_scan_sat(reference_seq, liftoff_seq, gap_open, gap_extend, matrix)

        liftoff_matches, liftoff_length = get_id_fraction(liftoff_parasail_res.traceback.ref, liftoff_parasail_res.traceback.query)
        liftoff_identity = liftoff_matches/liftoff_length

        max_lifton_identity = 0
        # Mapping lifton
        lifton_seq = lifton_sequences[id]
        
        lifton_parasail_res = parasail.nw_trace_scan_sat(reference_seq, lifton_seq, gap_open, gap_extend, matrix)

        lifton_matches, lifton_length = get_id_fraction(lifton_parasail_res.traceback.ref, lifton_parasail_res.traceback.query)
        max_lifton_identity = max(max_lifton_identity, lifton_matches/lifton_length)
        # #print(matches, length)
        # #break
        fw.write(f"{id}\t{liftoff_identity}\t{max_lifton_identity}\n")

    elif in_liftoff and not in_lifton:
        fw_liftoff_only.write(id + "\n")
        liftoff_only_cnt += 1
    elif not in_liftoff and in_lifton:
        fw_lifton_only.write(id + "\n")
        lifton_only_cnt += 1
    elif not in_liftoff and not in_lifton:
        fw_miss_both.write(id + "\n")
        miss_both_cnt += 1