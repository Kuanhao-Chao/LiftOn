from Bio import Seq

sequence = Seq.Seq("ATGTTGCTCTGATGAGGGGTGAAGCAAGGTTACAGTGAGAAGGGCCTGGAGGGAGGAGGTCCTGGAGGAGGGGGG")

# Find ORFs manually
start_codon = "ATG"
stop_codons = ["TAA", "TAG", "TGA"]

orf_list = []

for frame in range(3):
    for i in range(frame, len(sequence), 3):
        codon = str(sequence[i:i+3])
        if codon == start_codon:
            orf = ""
            for j in range(i, len(sequence), 3):
                codon = str(sequence[j:j+3])
                if codon in stop_codons:
                    break
                orf += codon
            if orf:
                orf_list.append(orf)

# Print the ORFs
for i, orf in enumerate(orf_list):
    print(f"ORF {i+1}: {orf}")
