from Bio.Seq import Seq

def determine_file_format(file_path):
    """Detect if the input file is GTF or GFF3 format.
    
    GTF format uses space-separated key-value pairs with gene_id and transcript_id.
    GFF3 format uses semicolon-separated key=value pairs with ID and Parent.
    
    Default format is GFF3 (most common and recommended format).
    
    Returns:
        "GTF format" if GTF format is detected
        "GFF format" if GFF3 format is detected (default)
        "Unknown format" if format cannot be determined (should be treated as GFF3)
    """
    gff_keys = {"ID", "Parent", "Name"}
    gtf_keys = {"gene_id", "transcript_id"}
    gtf_count = 0
    gff_count = 0
    lines_checked = 0
    max_lines_to_check = 100  # Check first 100 non-comment lines
    
    try:
        with open(file_path, 'r') as file:
            for line in file:
                # Skip header or comment lines
                if line.startswith("#"):
                    continue
                
                columns = line.strip().split("\t")
                if len(columns) != 9:
                    # Not a valid GTF/GFF file line, skip
                    continue
                
                attributes_column = columns[8]
                
                # GTF format: key-value pairs separated by spaces and semicolons
                # e.g., "gene_id \"ENSG00000139618\"; transcript_id \"ENST00000380152\";"
                # Check for GTF-specific patterns (key followed by space and quoted value)
                # GTF uses: gene_id "value"; transcript_id "value";
                if any((key + " " in attributes_column or key + "\"" in attributes_column) and 
                       not (key + "=" in attributes_column) for key in gtf_keys):
                    gtf_count += 1
                
                # GFF3 format: key=value pairs separated by semicolons
                # e.g., "ID=gene1;Parent=;Name=Gene1"
                # Check for GFF3-specific patterns (key=value)
                if any(key + "=" in attributes_column for key in gff_keys):
                    gff_count += 1
                
                lines_checked += 1
                if lines_checked >= max_lines_to_check:
                    break
        
        # Determine format based on counts
        # GTF format is more specific, so if we see GTF patterns (without =), it's likely GTF
        # Default to GFF3 format (most common and recommended)
        if gtf_count > 0 and gtf_count > gff_count:
            return "GTF format"
        elif gff_count > 0:
            return "GFF format"
        else:
            # No clear indicators found, default to GFF3 (most common format)
            return "GFF format"
            
    except Exception as e:
        # If file cannot be read, default to GFF3 format
        return "GFF format"


def extract_features(ref_db, features, ref_fai):
    ref_trans = {}
    ref_proteins = {}
    counter = 0
    for feature in features:
        for locus in ref_db.db_connection.features_of_type(feature):
            # print(f"Extracting features for {locus.id}")
            counter += 1
            __inner_extract_feature(ref_db, locus, ref_fai, ref_trans, ref_proteins)
    print(f"Extracted features for {counter} features")
    return ref_trans, ref_proteins


def __inner_extract_feature(ref_db, feature, ref_fai, ref_trans, ref_proteins):
    # If exon is the first level children
    children_exons = list(ref_db.db_connection.children(feature, featuretype='exon', level=1))
    children_CDSs = list(ref_db.db_connection.children(feature, featuretype=('start_codon', 'CDS', 'stop_codon'), level=1))
    # print(f"Parent: {feature.id};  exon: {len(children_exons)}; CDS: {len(children_CDSs)}")
    if len(children_exons) > 0 or len(children_CDSs) > 0:
        # print(f"Parent: {feature.id};  exon: {len(children_exons)}; CDS: {len(children_CDSs)}")
        if len(children_exons) > 0:
            trans_seq = get_dna_sequence(feature, ref_fai, children_exons)
            trans_seq = trans_seq.upper()
            ref_trans[feature.id] = trans_seq
        if len(children_CDSs) > 0:
            protein_seq = get_protein_sequence(feature, ref_fai, children_CDSs)
            protein_seq = protein_seq.upper()
            ref_proteins[feature.id] = protein_seq
    else:
        for child in ref_db.db_connection.children(feature, level=1, order_by='start'):
            __inner_extract_feature(ref_db, child, ref_fai, ref_trans, ref_proteins)
            

def merge_children_intervals(children):
    if not children:
        return []
    intervals = [[child.start, child.end] for child in children]
    intervals.sort(key=lambda interval: interval[0])
    merged = [intervals[0]]
    for current in intervals[1:]:
        previous = merged[-1]
        if current[0] <= previous[1]:
            previous[1] = max(previous[1], current[1])
        else:
            merged.append(current)
    return merged


def get_dna_sequence(parent_feature, fasta, features):
    chrom = parent_feature.seqid
    strand = parent_feature.strand
    merged_features = merge_children_intervals(features)
    sequence =''
    if chrom not in fasta.keys():
        return sequence
    for start, end in merged_features:
        sequence += str(fasta[chrom][start -1: end])
    if strand == "-":
        sequence = str(Seq(sequence).reverse_complement())
    sequence += 'N' * get_padding_length(len(sequence))
    return sequence


def get_padding_length(sequence_length):
    return (3-sequence_length%3)%3


def get_protein_sequence(parent_feature, fasta, features):
    dna = get_dna_sequence(parent_feature, fasta, features)
    return str(Seq(dna).translate())