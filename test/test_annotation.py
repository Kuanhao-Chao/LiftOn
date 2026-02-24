import pytest
import os
from lifton import extract_sequence
from intervaltree import Interval

def test_determine_file_format_gff3(tmp_path):
    # Create a mock GFF3 file
    gff3_file = tmp_path / "test.gff3"
    gff3_file.write_text(
        "##gff-version 3\n"
        "chr1\tLiftOn\tgene\t100\t200\t.\t+\t.\tID=gene1;Name=Gene1\n"
        "chr1\tLiftOn\tmRNA\t100\t200\t.\t+\t.\tID=rna1;Parent=gene1\n"
    )
    format_type = extract_sequence.determine_file_format(str(gff3_file))
    assert format_type == "GFF format", f"Expected 'GFF format', got {format_type}"

def test_determine_file_format_gtf(tmp_path):
    # Create a mock GTF file
    gtf_file = tmp_path / "test.gtf"
    gtf_file.write_text(
        "chr1\tLiftOn\tgene\t100\t200\t.\t+\t.\tgene_id \"ENSG01\"; transcript_id \"ENST01\";\n"
        "chr1\tLiftOn\texon\t100\t200\t.\t+\t.\tgene_id \"ENSG01\"; transcript_id \"ENST01\"; exon_number \"1\";\n"
    )
    format_type = extract_sequence.determine_file_format(str(gtf_file))
    assert format_type == "GTF format", f"Expected 'GTF format', got {format_type}"

def test_determine_file_format_empty(tmp_path):
    # Empty file should default to GFF format
    empty_file = tmp_path / "empty.gff3"
    empty_file.write_text("")
    format_type = extract_sequence.determine_file_format(str(empty_file))
    assert format_type == "GFF format", f"Expected 'GFF format' for empty file, got {format_type}"

def test_determine_file_format_missing():
    # Missing file should default to GFF format
    format_type = extract_sequence.determine_file_format("nonexistent_file.gff3")
    assert format_type == "GFF format", f"Expected 'GFF format' for missing file, got {format_type}"
