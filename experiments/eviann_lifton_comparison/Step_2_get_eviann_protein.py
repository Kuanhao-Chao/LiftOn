import subprocess

output_protein_fa = "/ccb/salz2/kh.chao/Lifton/data/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna_AA.fa"

genome = "/ccb/salz3/kh.chao/ref_genome/homo_sapiens/T2T-CHM13/ncbi/data/GCF_009914755.1/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna"

refernece_gff3 = "/ccb/salz2/kh.chao/Lifton/data/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna.gff"

command = ['gffread', '-S', '-y', output_protein_fa, '-g', genome, refernece_gff3]

gffread_extract = subprocess.call(command)