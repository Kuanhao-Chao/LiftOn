from alignment import *
from filepaths import VARIANTS_OUTPUTS, build_filepath, make_file
from sequences import get_padding_length
import warnings

def analyze_variants(ref_proteins, target_proteins, ref_trans, target_trans, target_db, ref_db, args):
    print("Analyzing variants")
    variant_effects = build_filepath([args.dir, VARIANTS_OUTPUTS['variants']])
    if make_file(variant_effects, args.force):
        with open(variant_effects, 'w') as of_ve:
            for tran in target_trans:
                tran_feature = target_db.db_connection[tran]
                if tran_feature.featuretype == "gene" and len(target_trans[tran]) == tran_feature.end-tran_feature.start +1 + get_padding_length(tran_feature.end - tran_feature.start +1):
                    no_transcripts_warning = "Warning: Gene " + tran + " has no transcripts.Skipping"
                    warnings.warn(no_transcripts_warning)
                else:
                    source_name = target_db.get_source_name(tran)
                    if  source_name in ref_trans:
                        align_dna, seq_id_dna= align_and_evaluate_seqs(source_name, tran, ref_trans, target_trans)
                        if tran in target_proteins and source_name in ref_proteins:
                            align_protein, seq_id_protein = align_and_evaluate_seqs(source_name, tran,ref_proteins,target_proteins)
                            variant_effect = find_variants(align_dna, align_protein, seq_id_dna, seq_id_protein)
                        else:
                            variant_effect, seq_id_protein = 'NA' , 'NA'
                        write_output(source_name, tran, seq_id_dna, seq_id_protein, variant_effect, of_ve)
            add_unmapped_trans(ref_trans, target_trans, ref_db, of_ve)


def align_and_evaluate_seqs(ref_name, target_name, ref_db, target_db):
    ref_seq, target_seq = ref_db[ref_name].upper(), target_db[target_name].upper()
    if ref_seq != target_seq:
        if ref_db.is_protein:
            try:
                align = ProteinAlignment(ref_seq, target_seq)
            except AttributeError as error:
                print("This sequence encountered a failure with parasail and will be skipped: ")
                print(ref_name, target_name)
                return None, 1.0
        else:
            try:
                align = DNAAlignment(ref_seq, target_seq)
            except AttributeError as error:
                print("This sequence encountered a failure with parasail and will be skipped: ")
                print(ref_name, target_name)
                return None, 1.0
        try:
            seq_id = align.calculate_seq_id()
        except AttributeError as error:
            print("This sequence encountered a failure with parasail and will be skipped: ")
            print(ref_name, target_name)
            return None, 1.0
    else:
        return None, 1.0
    return align.alignment, seq_id


def find_variants(align_dna, align_protein, seq_id_dna, seq_id_protein):
    if seq_id_dna == 1.0:
        return 'identical'
    if seq_id_protein == 1.0:
        return 'synonymous'
    ref_align = align_dna.traceback.query
    target_align = align_dna.traceback.ref
    if ref_align[0] == "-":
        return "5'_truncated"
    if target_align[-1] == "-":
        return "3'_truncated"
    if target_align[0:3] != 'ATG' and ref_align[0:3] != target_align[0:3]:
        return 'start_lost'
    if is_frameshift(ref_align):
        return 'frameshift'
    if is_frameshift(target_align):
        return 'frameshift'
    if has_stop_codon(align_protein.traceback.query,align_protein.traceback.ref):
        return 'stop_gained'
    if "-" in ref_align:
        return 'inframe_insertion'
    if "-" in target_align:
        return 'inframe_deletion'
    return 'nonsynonymous'


def is_frameshift(dna):
    prev_letter=''
    count = 0
    for letter in dna:
        if letter!=prev_letter:
            if prev_letter == "-" and count %3 !=0:
                return True
            count = 1
        else:
            count += 1
        prev_letter = letter
    return False


def has_stop_codon(ref_align, target_align):
    for i, letter in enumerate(target_align):
        if letter == "*" and ref_align[i] != "*":
            return True
    return False


def write_output(source_name, tran, seq_id_dna, seq_id_protein, variant_effect,of_ve):
    out_values = [source_name,tran, str(seq_id_dna), str(seq_id_protein), variant_effect]
    of_ve.write("\t".join(out_values) + "\n")


def add_unmapped_trans(ref_trans, target_trans, ref_db, output_file):
    for tran in ref_trans:
        if tran not in target_trans and ref_db.db_connection[tran].featuretype != "gene":
            out_str= [tran, 'unmapped']
            output_file.write("\t".join(out_str) + "\n")
