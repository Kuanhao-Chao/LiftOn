from gffutils import FeatureDB
import os, sys
import subprocess

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

if TARGET == "Han1" or TARGET == "Ash1" or TARGET == "PR1" or TARGET == "CHM13_MANE" or TARGET == "CHM13_GRCh38" or TARGET == "GRCh38_RefSeq" or TARGET == "GRCh38_MANE" or TARGET == "Mus_musculus":
    print("Running with ", TARGET)
    genome = ""
    if TARGET == "Han1":
        genome = "/ccb/salz2/kh.chao/PR_liftoff_protein_search/data/Han1/v1.0/Assembly/Han1_v1.2.fasta"
    elif TARGET == "Ash1":
        genome = "/ccb/salz2/kh.chao/PR_liftoff_protein_search/data/Ash1/v2.2/Assembly/Ash1_v2.2.fa"
    elif TARGET == "PR1":
        genome = "/ccb/salz2/kh.chao/PR_liftoff_protein_search/data/PR1/v3.0/Assembly/PR1.fa"
    elif TARGET == "GRCh38_RefSeq" or TARGET == "GRCh38_MANE":
        genome = "/ccb/salz2/kh.chao/PR_liftoff_protein_search/data/NCBI_Refseq_chr_fixed/GCF_000001405.40_GRCh38.p14_genomic.fna"
    elif TARGET == "CHM13_MANE" or TARGET == "CHM13_GRCh38":
        genome = "/ccb/salz2/kh.chao/ref_genome/homo_sapiens/T2T-CHM13/chm13v2.0.fa"
    elif TARGET == "Mus_musculus":
        genome = "/ccb/salz2/kh.chao/PR_liftoff_protein_search/data/mus_musculus/NCBI_Refseq/GCF_000001635.27_GRCm39_genomic.fna"
else:
    sys.exit(-1)

print("Your target: ", TARGET)
      
output_dir = f"/ccb/salz2/kh.chao/PR_liftoff_protein_search/results/liftoff/{TARGET}/MANE_liftoff_miniprot_cmp/"

os.makedirs(output_dir, exist_ok=True)


frname_both = output_dir + "identities_coords.txt"

fwname_both_liftoff = output_dir + "identities_coords_liftoff_better.txt"
fwname_both_miniprot = output_dir + "identities_coords_miniprot_better.txt"

command = "awk '$2 > $3+0.2' " + frname_both + " > " + fwname_both_liftoff
print(command)
# fw = open(fwname_both_liftoff, "w")
# subprocess.call(commands, stdout=fw)
# fw.close()

command = "awk '$2+0.2 < $3' " + frname_both + " > " + fwname_both_miniprot
print(command)

# commands = ["awk", '\'$2+0.2 < $3\'', frname_both]
# print(commands)
# fw = open(fwname_both_miniprot, "w")
# subprocess.call(commands, stdout=fw)
# fw.close()
