python setup.py install

lifton -g /ccb/salz3/kh.chao/ref_genome/homo_sapiens/chess/chess3.0.1.gff -dir results/CHM13_CHESS/intermediate_files -o results/CHM13_CHESS/CHM13_CHESS_lifton.gff3 -chroms data/chroms_mapping.csv -f data/features_chess.txt /ccb/salz3/kh.chao/ref_genome/homo_sapiens/T2T-CHM13/chm13v2.0.fa /ccb/salz2/kh.chao/PR_liftoff_protein_search/data/NCBI_Refseq_chr_fixed/GCF_000001405.40_GRCh38.p14_genomic.fna -t 20 --proteins /ccb/salz2/kh.chao/Lifton/results/CHM13_CHESS/proteins.fa

