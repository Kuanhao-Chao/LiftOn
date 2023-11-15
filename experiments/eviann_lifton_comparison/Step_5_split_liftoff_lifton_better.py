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


TARGET = "CHM13_MANE"

print("Your target: ", TARGET)

output_dir = f"/ccb/salz2/kh.chao/Lifton/results/{TARGET}/eviann_lifton_cmp/"

os.makedirs(output_dir, exist_ok=True)


frname_both = output_dir + "identities_coords.txt"

frname_both_truncated = output_dir + "identities_coords_truncated.txt"

fwname_both_eviann = output_dir + "identities_coords_eviann_better.txt"
fwname_both_lifton = output_dir + "identities_coords_lifton_better.txt"


command = "awk '{if($2 < 1 || $3 < 1) {print}}' " + frname_both + " > " + frname_both_truncated
print(command)

os.system(command=command)


command = "awk '$2 > $3' " + frname_both + " > " + fwname_both_eviann
print(command)

os.system(command=command)

# commands = ["awk", '\'$2 > $3+0.2\'', frname_both]

# fw = open(fwname_both_eviann, "w")
# subprocess.call(commands, stdout=fw)
# fw.close()

command = "awk '$2 < $3' " + frname_both + " > " + fwname_both_lifton
print(command)
os.system(command=command)



# commands = ["awk", '\'$2+0.2 < $3\'', frname_both]
# # print(commands)
# fw = open(fwname_both_lifton, "w")
# subprocess.call(commands, stdout=fw)
# fw.close()
