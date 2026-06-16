from Bio.Seq import Seq
from lifton import logger
from lifton.exceptions import LiftOnInputError

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


# ---------------------------------------------------------------------------
# Phase 15b — streaming extractor. Writes FASTA records as it goes,
# never materialising the full {id: seq} dict. Caller re-opens the
# resulting files via pyfaidx.Fasta for lazy mmap-style access.
# RAM ceiling: ~one feature at a time (mid-100s of bytes).
# ---------------------------------------------------------------------------

import os as _os


def extract_features_to_fasta(ref_db, features, ref_fai, out_dir):
    """Streaming alternative to :func:`extract_features`.

    Writes ``transcripts.fa`` and ``proteins.fa`` under ``out_dir`` as
    each reference feature's sequence is computed; nothing is held in
    Python dict form. Returns the two file paths so callers can hand
    them straight to ``pyfaidx.Fasta(path)``.

    The byte content of each FASTA record is identical to the legacy
    path's ``write_seq_2_file`` output (same record header, same
    in-order traversal, same uppercased sequence + N-padding) so
    downstream identity scores do not drift.
    """
    _os.makedirs(out_dir, exist_ok=True)
    trans_path = _os.path.join(out_dir, "transcripts.fa")
    prot_path = _os.path.join(out_dir, "proteins.fa")
    counter = 0
    with open(trans_path, "w") as ft, open(prot_path, "w") as fp:
        for feature in features:
            for locus in ref_db.db_connection.features_of_type(feature):
                counter += 1
                _stream_inner(ref_db, locus, ref_fai, ft, fp)
    print(f"Extracted features (streaming) for {counter} features")
    return trans_path, prot_path


def _stream_inner(ref_db, feature, ref_fai, ft, fp):
    # Phase 18: one ordered children() query + an in-Python featuretype
    # partition replaces the prior three per-feature queries (2x fewer
    # SQLite round-trips on the common exon/CDS-bearing path, 3x on a bare
    # gene parent). Byte-neutral: merge_children_intervals re-sorts the
    # exon/CDS lists by start before concatenation, so the SQL row order is
    # discarded for those branches; the recursive-descent branch below
    # keeps order_by='start' (the only order-sensitive path).
    all_children = list(ref_db.db_connection.children(
        feature, level=1, order_by='start'))
    children_exons = [c for c in all_children if c.featuretype == 'exon']
    children_CDSs = [c for c in all_children
                     if c.featuretype in ('start_codon', 'CDS', 'stop_codon')]
    if len(children_exons) > 0 or len(children_CDSs) > 0:
        if len(children_exons) > 0:
            try:
                trans_seq = get_dna_sequence(feature, ref_fai, children_exons)
                if trans_seq:
                    trans_seq = trans_seq.upper()
                    ft.write(f">{feature.id}\n{trans_seq}\n")
            except Exception as e:
                logger.log_warning(
                    f"extract_features_to_fasta: transcript {feature.id}: {e}"
                )
        if len(children_CDSs) > 0:
            try:
                protein_seq = get_protein_sequence(
                    feature, ref_fai, children_CDSs)
                if protein_seq:
                    protein_seq = protein_seq.upper()
                    fp.write(f">{feature.id}\n{protein_seq}\n")
            except Exception as e:
                logger.log_warning(
                    f"extract_features_to_fasta: protein {feature.id}: {e}"
                )
    else:
        for child in all_children:
            _stream_inner(ref_db, child, ref_fai, ft, fp)


def __inner_extract_feature(ref_db, feature, ref_fai, ref_trans, ref_proteins):
    # If exon is the first level children
    # Phase 18: collapse the prior 3 children() queries into one ordered
    # query + in-Python partition (see _stream_inner for the byte-neutral
    # argument; merge_children_intervals re-sorts, the else-branch keeps
    # order_by='start').
    all_children = list(ref_db.db_connection.children(feature, level=1, order_by='start'))
    children_exons = [c for c in all_children if c.featuretype == 'exon']
    children_CDSs = [c for c in all_children if c.featuretype in ('start_codon', 'CDS', 'stop_codon')]
    # print(f"Parent: {feature.id};  exon: {len(children_exons)}; CDS: {len(children_CDSs)}")
    if len(children_exons) > 0 or len(children_CDSs) > 0:
        if len(children_exons) > 0:
            # V1.2 fix: log instead of swallowing silently. A user whose
            # transcript silently disappears can't debug it; a [WARNING]
            # gives them the feature id and the underlying cause.
            try:
                trans_seq = get_dna_sequence(feature, ref_fai, children_exons)
                if trans_seq:
                    trans_seq = trans_seq.upper()
                    ref_trans[feature.id] = trans_seq
            except Exception as e:
                logger.log_warning(
                    f"extract_features: failed to extract transcript "
                    f"sequence for {feature.id}: {e}"
                )

        if len(children_CDSs) > 0:
            try:
                protein_seq = get_protein_sequence(feature, ref_fai, children_CDSs)
                if protein_seq:
                    protein_seq = protein_seq.upper()
                    ref_proteins[feature.id] = protein_seq
            except Exception as e:
                logger.log_warning(
                    f"extract_features: failed to extract protein "
                    f"sequence for {feature.id}: {e}"
                )
    else:
        for child in all_children:
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
        # V2.1 fix: GFF3 is 1-based inclusive; pyfaidx is 0-based half-open.
        # Negative `start - 1` silently wraps to the chromosome's tail and
        # corrupts the extracted sequence. Reject explicitly.
        if start < 1:
            raise LiftOnInputError(
                f"get_dna_sequence: feature on {chrom} has start={start}; "
                "GFF3 coordinates must be 1-based inclusive (start >= 1)."
            )
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