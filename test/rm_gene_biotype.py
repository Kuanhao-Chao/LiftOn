def remove_gene_biotype(gff_line):
    parts = gff_line.split('\t')
    attributes = parts[-1]
    attributes_dict = dict(item.split('=') for item in attributes.split(';') if '=' in item)
    attributes_dict.pop('gene_biotype', None)
    modified_attributes = ';'.join(f"{key}={value}" for key, value in attributes_dict.items())
    parts[-1] = modified_attributes
    return '\t'.join(parts)

def process_gff_file(input_file, output_file):
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            # Check if the line contains data to modify
            if line.strip() and not line.startswith('#'):  # Ignore empty and comment lines
                modified_line = remove_gene_biotype(line.strip())
                outfile.write(modified_line + '\n')
            else:
                # Write the line as is if it's a comment or empty
                outfile.write(line)

# Specify the input and output file names
input_gff = 'GRCh38_chr22.gff3'
output_gff = 'GRCh38_chr22_no_biotype.gff3'

# Process the file
process_gff_file(input_gff, output_gff)

print(f"Processed file saved as '{output_gff}'.")
