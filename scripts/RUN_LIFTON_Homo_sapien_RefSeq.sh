python setup.py install

lifton -g /ccb/salz3/kh.chao/PR_liftoff_protein_search/data/NCBI_Refseq_chr_fixed/rRNA_removed/NCBI_RefSeq_no_rRNA.gff_db -dir results/CHM13_RefSeq/intermediate_files -o results/CHM13_RefSeq/CHM13_RefSeq_lifton.gff3 -chroms data/chroms_mapping.csv -f data/features.txt /ccb/salz3/kh.chao/ref_genome/homo_sapiens/T2T-CHM13/chm13v2.0.fa /ccb/salz2/kh.chao/PR_liftoff_protein_search/data/NCBI_Refseq_chr_fixed/GCF_000001405.40_GRCh38.p14_genomic.fna -t 10

