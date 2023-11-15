####################
# Goal of this script is to count the number of 
# 1. genes, 
# 2. transcripts, 
# 3. protein-coding genes,
# 4. protein-coding transcripts
####################

import argparse

def count_features(gff_file):
    gene_count = 0
    transcript_count = 0
    protein_coding_gene_count = 0
    protein_coding_transcript_count = 0

    with open(gff_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) < 3:
                continue

            feature_type = fields[2]
            attributes = fields[-1]

            if feature_type == 'gene':
                gene_count += 1
                if 'protein_coding' in attributes:
                    protein_coding_gene_count += 1

            if feature_type == 'mRNA' or feature_type == 'transcript':
                transcript_count += 1
                if feature_type == 'mRNA' or 'protein_coding' in attributes:
                    protein_coding_transcript_count += 1

    return gene_count, transcript_count, protein_coding_gene_count, protein_coding_transcript_count

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Count features in a GFF file')
    parser.add_argument('gff_file', help='Input GFF file')
    args = parser.parse_args()

    gene_count, transcript_count, protein_coding_gene_count, protein_coding_transcript_count = count_features(args.gff_file)

    print(f'Total number of genes: {gene_count}')
    print(f'Total number of transcripts: {transcript_count}')
    print(f'Total number of protein-coding genes: {protein_coding_gene_count}')
    print(f'Total number of protein-coding transcripts: {protein_coding_transcript_count}')