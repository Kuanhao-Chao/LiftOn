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

# using the definition of BLAST identity
#   defined as the number of matching bases over the number of alignment columns.
def get_id_fraction(reference, target):
    matches = 0
    cmp_column = 0
    for i, letter in enumerate(reference):
        if (target[i] == "*" or target[i] == ".") and (i != (len(reference) - 1)):
            break        
        
        cmp_column += 1
        if letter == target[i]:
            matches += 1
    return matches, len(reference)

TARGET = sys.argv[1]

protein_fa = ""
if TARGET == "CHM13_MANE" or TARGET == "CHM13_RefSeq" or TARGET == "GRCh38_RefSeq" or TARGET == "Han1" or TARGET == "Ash1" or TARGET == "PR1" or TARGET == "Mus_musculus_MANE":
    print("Running with ", TARGET)
    genome = ""
    if TARGET == "GRCh38_RefSeq":
        genome = "/ccb/salz3/kh.chao/ref_genome/homo_sapiens/NCBI_Refseq_chr_fixed/GCF_000001405.40_GRCh38.p14_genomic.fna"
        protein_fa = f"/ccb/salz2/kh.chao/PR_liftoff_protein_search/data/protein.fasta"

    elif TARGET == "CHM13_MANE": 
        genome = "/ccb/salz3/kh.chao/ref_genome/homo_sapiens/T2T-CHM13/chm13v2.0.fa"
        protein_fa = f"/ccb/salz2/kh.chao/PR_liftoff_protein_search/data/MANE_RefSeq/MANE.GRCh38.v1.2.refseq_genomic.AA.cleaned.fa"

    elif TARGET == "CHM13_RefSeq":
        genome = "/ccb/salz3/kh.chao/ref_genome/homo_sapiens/T2T-CHM13/chm13v2.0.fa"
        protein_fa = f"/ccb/salz2/kh.chao/PR_liftoff_protein_search/data/protein.fasta"

    elif TARGET == "Han1":
        genome = "/ccb/salz3/kh.chao/PR_liftoff_protein_search/data/Han1/v1.0/Assembly/Han1_v1.2.fasta"
        protein_fa = f"/ccb/salz2/kh.chao/PR_liftoff_protein_search/data/MANE_RefSeq/MANE.GRCh38.v1.2.refseq_genomic.AA.cleaned.fa"

    elif TARGET == "Ash1":
        genome = "/ccb/salz3/kh.chao/PR_liftoff_protein_search/data/Ash1/v2.2/Assembly/Ash1_v2.2.fa"
        protein_fa = f"/ccb/salz2/kh.chao/PR_liftoff_protein_search/data/MANE_RefSeq/MANE.GRCh38.v1.2.refseq_genomic.AA.cleaned.fa"

    elif TARGET == "PR1":
        genome = "/ccb/salz3/kh.chao/PR_liftoff_protein_search/data/PR1/v3.0/Assembly/PR1.fa"
        protein_fa = f"/ccb/salz2/kh.chao/PR_liftoff_protein_search/data/MANE_RefSeq/MANE.GRCh38.v1.2.refseq_genomic.AA.cleaned.fa"

    elif TARGET == "Mus_musculus_MANE":
        genome = "/ccb/salz3/kh.chao/ref_genome/mus_musculus/NCBI_Refseq/GCF_000001635.27_GRCm39_genomic.fna"
        protein_fa = f"/ccb/salz2/kh.chao/PR_liftoff_protein_search/data/MANE_RefSeq/MANE.GRCh38.v1.2.refseq_genomic.AA.cleaned.fa"
else:
    sys.exit(-1)



output_dir = f"/ccb/salz2/kh.chao/Lifton/results/{TARGET}/liftoff_lifton_cmp/"

os.makedirs(output_dir, exist_ok=True)

liftoff_protein_fa = f"/ccb/salz2/kh.chao/PR_liftoff_protein_search/results/liftoff/{TARGET}/{TARGET}_liftoff_AA.fa"
lifton_protein_fa = f"/ccb/salz2/kh.chao/Lifton/results/{TARGET}/{TARGET}_lifton_protein.fa"
lifton_protein_gff = f"/ccb/salz2/kh.chao/Lifton/results/{TARGET}/{TARGET}_lifton.gff3" 


print("protein_fa: ", protein_fa)
print("liftoff_protein_fa: ", liftoff_protein_fa)
print("lifton_protein_fa: ", lifton_protein_fa)
print("lifton_protein_gff: ", lifton_protein_gff)

liftoff_sequences = read_fasta(liftoff_protein_fa)
lifton_sequences = read_fasta(lifton_protein_fa)
MANE_sequences = read_fasta(protein_fa)

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

EARLY_STOP_IN_REFERENCE = 0
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

        # if "*" in reference_seq or "." in reference_seq:
        #     EARLY_STOP_IN_REFERENCE += 1
        #     continue

        # Mapping liftoff
        liftoff_seq = liftoff_sequences[id]
        
        liftoff_parasail_res = parasail.nw_trace_scan_sat(liftoff_seq+"*", reference_seq+"*", gap_open, gap_extend, matrix)

        # print("liftoff_seq: ", liftoff_seq+"*")
        # print("reference_seq: ", reference_seq+"*")

        liftoff_matches, liftoff_length = get_id_fraction(liftoff_parasail_res.traceback.ref, liftoff_parasail_res.traceback.query)
        liftoff_identity = liftoff_matches/liftoff_length

        # print("liftoff_matches: ", liftoff_matches)
        # print("liftoff_length: ", liftoff_length)
        # print("liftoff_identity: ", liftoff_identity)

        # if id == "rna-NM_012230.5" or id == "rna-NM_001005183.1" or id =="rna-NM_001365389.2" or id == "rna-NM_181612.3":
        #     print("Liftoff; ", id)
        #     print(f"{id}\t{liftoff_parasail_res.traceback.ref}\t{liftoff_parasail_res.traceback.query}\t{liftoff_matches}\t{liftoff_length}\t{liftoff_identity}\n")


        max_lifton_identity = 0
        # Mapping lifton
        lifton_seq = lifton_sequences[id]
        
        lifton_parasail_res = parasail.nw_trace_scan_sat(lifton_seq+"*", reference_seq+"*", gap_open, gap_extend, matrix)


        lifton_matches, lifton_length = get_id_fraction(lifton_parasail_res.traceback.ref, lifton_parasail_res.traceback.query)
        lifton_identity = lifton_matches/lifton_length

        if id == "rna-NM_012230.5" or id == "rna-NM_001005183.1" or id =="rna-NM_001365389.2" or id == "rna-NM_181612.3":
            print("Liftoff; ", id)
            print(f"{id}\t{lifton_parasail_res.traceback.ref}\t{lifton_parasail_res.traceback.query}\t{lifton_matches}\t{lifton_length}\t{lifton_identity}\n")

        fw.write(f"{id}\t{liftoff_identity}\t{lifton_identity}\n")

    elif in_liftoff and not in_lifton:
        fw_liftoff_only.write(id + "\n")
        liftoff_only_cnt += 1
    elif not in_liftoff and in_lifton:
        fw_lifton_only.write(id + "\n")
        lifton_only_cnt += 1
    elif not in_liftoff and not in_lifton:
        fw_miss_both.write(id + "\n")
        miss_both_cnt += 1

print("EARLY_STOP_IN_REFERENCE: ", EARLY_STOP_IN_REFERENCE)