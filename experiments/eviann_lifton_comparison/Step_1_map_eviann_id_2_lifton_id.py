import gffutils 
from intervaltree import Interval, IntervalTree
import json 
import os

output_dir = "/ccb/salz2/kh.chao/Lifton/results/CHM13_MANE/eviann_lifton_cmp" 
os.makedirs(output_dir, exist_ok=True)

ref_gff = f"/ccb/salz2/kh.chao/PR_liftoff_protein_search/data/MANE_RefSeq/MANE.GRCh38.v1.2.refseq_genomic.cleaned.gff"

lifton_gff = f"/ccb/salz2/kh.chao/Lifton/results/CHM13_MANE/CHM13_MANE_lifton.gff3" 

eviann_gff = "/ccb/salz2/kh.chao/Lifton/data/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna.chr_fix.gff"

chr_id_json = "/ccb/salz3/kh.chao/ref_genome/homo_sapiens/T2T-CHM13/ncbi/data/GCF_009914755.1/sequence_report.jsonl"

dbfn_ref_gff = ref_gff + "_db"
dbfn_eviann_gff = eviann_gff + "_db"
dbfn_lifton_gff = lifton_gff + "_db"

# Create the GFF database
# db_eviann = gffutils.create_db(
#     eviann_gff,
#     dbfn=dbfn_eviann_gff,
#     merge_strategy="merge",
#     force=True,  # Set this to True to overwrite the database if it already exists    
# )

db_eviann = gffutils.FeatureDB(dbfn_eviann_gff)
db_lifton = gffutils.FeatureDB(dbfn_lifton_gff)

# chr_dict = {}
# with open(chr_id_json, "r") as fr:
#     lines = fr.read().splitlines()
#     for line in lines:
#         tmp_dict = json.loads(line)
#         chr_dict[tmp_dict["refseqAccession"]] = tmp_dict["ucscStyleName"]

# print(chr_dict)





################################
# Step 2: Initializing intervaltree
################################
tree_dict = {}
chr_num_ls = [*range(1, 23)] 
chr_num_ls += ['X', 'Y', 'M']
for i in chr_num_ls:
    tree = IntervalTree()
    tree_dict["chr" + str(i)] = tree

for gene in db_eviann.features_of_type('mRNA'):#, limit=("chr4", 122344, 135031)):
    chromosome = gene.seqid
    gene_id = gene.attributes["ID"][0]
    gene_id_base = gene_id.split("_")[0]
    # print("&& gene_id      : ", gene_id)
    start = gene.start
    end = gene.end

    
    gene_interval = Interval(start, end, gene_id)
    tree_dict[chromosome].add(gene_interval)

# print(tree_dict.keys())

count_1_to_1 = 0
count_1_to_many = 0
count_none = 0

id_mapping_dic = {}

for gene in db_lifton.features_of_type('mRNA'):#, limit=
    chromosome = gene.seqid
    gene_id = gene.attributes["ID"][0]
    gene_id_base = gene_id.split("_")[0]
    # print("&& gene_id      : ", gene_id)
    start = gene.start
    end = gene.end

    lifton_gene_interval = Interval(start, end, gene_id)

    ovps = tree_dict[chromosome].overlap(lifton_gene_interval)
    print(len(ovps))

    eviann_ids = []
    boundaries = []
    for ovp in ovps:

        if lifton_gene_interval[0] <= ovp[0]: 
            boundaries.append(lifton_gene_interval[0])
            boundaries.append(ovp[0])
        else:
            boundaries.append(ovp[0])
            boundaries.append(lifton_gene_interval[0])
        

        if lifton_gene_interval[1] <= ovp[1]: 
            boundaries.append(lifton_gene_interval[1])
            boundaries.append(ovp[1])
        else:
            boundaries.append(ovp[1])
            boundaries.append(lifton_gene_interval[1])


        ovp_len = boundaries[2] - boundaries[1] + 1

        # Overlaping ratio
        print('boundaries: ', boundaries)
        if (ovp_len / (lifton_gene_interval[1] - lifton_gene_interval[0] + 1)  > 0.7) and (ovp_len / (ovp[1] - ovp[0] + 1)  > 0.7):
            eviann_ids.append(ovp[2])
        boundaries = []

    id_mapping_dic[gene_id] = eviann_ids
    
    if len(ovps) == 0:
        count_none += 1
    elif len(ovps) == 1:
        count_1_to_1 += 1
    elif len(ovps) > 1:
        count_1_to_many += 1

print("LiftOn map to none of Eviann: ", count_none)
print("LiftOn map to == 1 of Eviann: ", count_1_to_1)
print("LiftOn map to > 1 of Eviann: ", count_1_to_many)

# Serializing json
json_object = json.dumps(id_mapping_dic)
 
# Writing to sample.json
with open(f"{output_dir}/gene_loci.json", "w") as outfile:
    outfile.write(json_object)