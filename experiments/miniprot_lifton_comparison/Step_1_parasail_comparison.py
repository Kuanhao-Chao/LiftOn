from gffutils import FeatureDB, create_db
import pybedtools
import os, re

import sys
from gffutils import FeatureDB
from Bio import SeqIO
import parasail

def segments_overlap(segment1, segment2):
    if len(segment1) != 2 or len(segment2) != 2:
        raise ValueError("Segments must have exactly 2 endpoints")
    segment1, segment2 = sorted([segment1, segment2], key=lambda x: x[0])
    return segment1[1] >= segment2[0]


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
        if target[i] == "*" or target[i] == ".":
            break
    return matches, max(len(reference), len(target))


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


output_dir = f"/ccb/salz2/kh.chao/Lifton/results/{TARGET}/miniprot_lifton_cmp/"

os.makedirs(output_dir, exist_ok=True)

lifton_protein_fa = f"/ccb/salz2/kh.chao/Lifton/results/{TARGET}/{TARGET}_lifton_protein.fa"

lifton_gff_db = f"/ccb/salz2/kh.chao/Lifton/results/{TARGET}/{TARGET}_lifton.gff3_db"

miniprot_protein_fa = f"/ccb/salz2/kh.chao/PR_liftoff_protein_search/results/liftoff/{TARGET}/{TARGET}_miniprot_AA.fa"

miniprot_protein_gff = f"/ccb/salz2/kh.chao/PR_liftoff_protein_search/results/liftoff/{TARGET}/{TARGET}_miniprot.fix.gff"

print("protein_fa: ", protein_fa)
print("lifton_protein_fa: ", lifton_protein_fa)
print("lifton_gff_db: ", lifton_gff_db)

print("miniprot_protein_fa: ", miniprot_protein_fa)
print("miniprot_protein_gff: ", miniprot_protein_gff)


if not os.path.exists(lifton_gff_db):
    l_feature_db = create_db(
        lifton_gff_db[:-3],
        dbfn=lifton_gff_db,
        merge_strategy="merge",
        verbose = True,
        force=True,  # Set this to True to overwrite the database if it already exists    
    )
else:
    l_feature_db = FeatureDB(lifton_gff_db)

lifton_gene_loci = {}
for gene in l_feature_db.features_of_type('gene'):#, limit=("chr1", 0, 250000000)):
    chromosome = gene.seqid
    gene_id = gene.attributes["ID"][0]

    # gene_interval = Interval(lifton_gene.entry.start, lifton_gene.entry.end, gene_id)
    # tree_dict[chromosome].add(gene_interval)

    # Assumption that all 1st level are transcripts
    transcripts = l_feature_db.children(gene, level=1)  # Replace 'exon' with the desired 
    
    ###########################
    # Adding transcripts
    ###########################
    for transcript in list(transcripts):
        transcript_id = transcript["ID"][0]
        
        lifton_gene_loci[transcript_id] = (transcript.start, transcript.end)



miniprot_gene_loci = {}
# Step 1: create miniprot ID -> MANE ID dictionary
miniprot_ids = {}
with open(miniprot_protein_gff, 'r') as fr:
    lines = fr.read().splitlines()
    for line in lines:
        eles = line.split("\t")
        if len(eles) < 5 : 
            continue

        m_start = int(eles[3])
        m_end = int(eles[4])

        if len(eles) > 5 and eles[2] == "mRNA":
            props = eles[8].split(";")
            # print(props)
            # ID = props[0][3:]
            pattern = r'ID=([^; ]+)'
            match = re.search(pattern, eles[8])
            if match:
                extracted_text = match.group(1)
                ID = extracted_text
            else:
                pass

            pattern = r'Target=([^; ]+)'
            # Use re.search() to find the pattern in the input string
            match = re.search(pattern, eles[8])
            if match:
                extracted_text = match.group(1)
                TARGET = extracted_text
            else:
                pass

            # TARGET = props[4][7:]
            # TARGET = TARGET.split(" ")[0]
            if TARGET not in miniprot_ids.keys():
                miniprot_ids[TARGET] = [ID]
            else:
                miniprot_ids[TARGET].append(ID)
            
            miniprot_gene_loci[ID] = (m_start, m_end)





# print(miniprot_ids)
lifton_sequences = read_fasta(lifton_protein_fa)
miniprot_sequences = read_fasta(miniprot_protein_fa)
protein_sequences = read_fasta(protein_fa)

both_cnt = 0
lifton_only_cnt = 0
miniprot_only_cnt = 0
miss_both_cnt = 0

fwname = output_dir + "identities.txt"
fwname_both = output_dir + "both.txt"
fwname_lifton_only = output_dir + "lifton_only.txt"
fwname_miniprot_only = output_dir + "miniprot_only.txt"
fwname_miss_both = output_dir + "miss_both.txt"

fw = open(fwname, "w")
fw_both = open(fwname_both, "w")
fw_lifton_only = open(fwname_lifton_only, "w")
fw_miniprot_only = open(fwname_miniprot_only, "w")
fw_miss_both = open(fwname_miss_both, "w")

EARLY_STOP_IN_REFERENCE = 0
# print("miniprot_sequences.keys(): ", miniprot_sequences.keys())
for id, sequence in protein_sequences.items():

    in_lifton = id in lifton_sequences.keys() and id in lifton_gene_loci.keys()
    in_miniprot = id in miniprot_ids.keys()
    
    # print([i for i in miniprot_ids[id] if i in miniprot_sequences.keys()])
    # in_miniprot_ls = [i for i in miniprot_ids[id] if i in miniprot_sequences.keys()]
    # in_miniprot = len(in_miniprot_ls) > 0

    if in_lifton and in_miniprot:
        fw_both.write(id + "\n")
        both_cnt += 1

        matrix = parasail.Matrix("blosum62")
        gap_open = 11
        gap_extend = 1

        reference_seq = protein_sequences[id]

        if "*" in reference_seq or "." in reference_seq:
            EARLY_STOP_IN_REFERENCE += 1
            # continue


        lifton_seq = lifton_sequences[id]
        
        lifton_parasail_res = parasail.nw_trace_scan_sat(lifton_seq, reference_seq, gap_open, gap_extend, matrix)

        lifton_matches, lifton_length = get_id_fraction(lifton_parasail_res.traceback.ref, lifton_parasail_res.traceback.query)
        lifton_identity = lifton_matches/lifton_length


        max_miniprot_identity = 0
        for miniprotID in miniprot_ids[id]:

            miniprot_trans_locus = miniprot_gene_loci[miniprotID]
            lifton_trans_locus = lifton_gene_loci[id]

            overlap = segments_overlap(miniprot_trans_locus, lifton_trans_locus)


            if overlap:
                # Set the condition for what is considered overlap
                boundaries = []
                if miniprot_trans_locus[0] <= lifton_trans_locus[0]: 
                    boundaries.append(miniprot_trans_locus[0])
                    boundaries.append(lifton_trans_locus[0])
                else:
                    boundaries.append(lifton_trans_locus[0])
                    boundaries.append(miniprot_trans_locus[0])
                

                if miniprot_trans_locus[1] <= lifton_trans_locus[1]: 
                    boundaries.append(miniprot_trans_locus[1])
                    boundaries.append(lifton_trans_locus[1])
                else:
                    boundaries.append(lifton_trans_locus[1])
                    boundaries.append(miniprot_trans_locus[1])

                ovp_len = boundaries[2] - boundaries[1] + 1

                # Overlaping ratio
                print('boundaries: ', boundaries)
                if (ovp_len / (miniprot_trans_locus[1]-miniprot_trans_locus[0]+1) > 0.7) and (ovp_len / (lifton_trans_locus[1] - lifton_trans_locus[0]+1) > 0.7):
                    miniprot_seq = miniprot_sequences[miniprotID]
                    # print(miniprotID)
                    miniprot_parasail_res = parasail.nw_trace_scan_sat(miniprot_seq, reference_seq, gap_open, gap_extend, matrix)
                    miniprot_matches, miniprot_length = get_id_fraction(miniprot_parasail_res.traceback.ref, miniprot_parasail_res.traceback.query)
                    max_miniprot_identity = max(max_miniprot_identity, miniprot_matches/miniprot_length)

        fw.write(f"{id}\t{lifton_identity}\t{max_miniprot_identity}\n")

    elif in_lifton and not in_miniprot:
        fw_lifton_only.write(id + "\n")
        lifton_only_cnt += 1
    elif not in_lifton and in_miniprot:
        fw_miniprot_only.write(id + "\n")
        miniprot_only_cnt += 1
    elif not in_lifton and not in_miniprot:
        fw_miss_both.write(id + "\n")
        miss_both_cnt += 1

fw.close()
fw_both.close()
fw_lifton_only.close()
fw_miniprot_only.close()
fw_miss_both.close()

print("EARLY_STOP_IN_REFERENCE: ", EARLY_STOP_IN_REFERENCE)