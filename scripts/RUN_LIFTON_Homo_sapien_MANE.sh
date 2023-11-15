python setup.py install

#lifton --proteins /ccb/salz3/kh.chao/PR_liftoff_protein_search/data/MANE_RefSeq/MANE.GRCh38.v1.2.refseq_genomic.AA.cleaned.fa -g /ccb/salz3/kh.chao/PR_liftoff_protein_search/data/MANE_RefSeq/MANE.GRCh38.v1.2.refseq_genomic.cleaned.gff_db --miniprot /ccb/salz2/kh.chao/Lifton/results/CHM13_MANE/miniprot.gff3_db -o results/CHM13_MANE/CHM13_MANE_lifton.gff3 /ccb/salz3/kh.chao/ref_genome/homo_sapiens/T2T-CHM13/chm13v2.0.fa /ccb/salz2/kh.chao/PR_liftoff_protein_search/data/NCBI_Refseq_chr_fixed/GCF_000001405.40_GRCh38.p14_genomic.fna

lifton -g /ccb/salz3/kh.chao/PR_liftoff_protein_search/data/MANE_RefSeq/MANE.GRCh38.v1.2.refseq_genomic.cleaned.gff -dir results/CHM13_MANE/intermediate_files -o results/CHM13_MANE/CHM13_MANE_lifton.gff3 -chroms data/chroms_mapping.csv -f data/features.txt /ccb/salz3/kh.chao/ref_genome/homo_sapiens/T2T-CHM13/chm13v2.0.fa /ccb/salz2/kh.chao/PR_liftoff_protein_search/data/NCBI_Refseq_chr_fixed/GCF_000001405.40_GRCh38.p14_genomic.fna --proteins /ccb/salz2/kh.chao/Lifton/results/CHM13_MANE/proteins.fa -t 30 --liftoff results/CHM13_MANE/liftoff/liftoff.gff3_db --miniprot results/CHM13_MANE/miniprot/miniprot.gff3_db

#lifton --proteins /ccb/salz3/kh.chao/PR_liftoff_protein_search/data/MANE_RefSeq/MANE.GRCh38.v1.2.refseq_genomic.AA.cleaned.fa -g /ccb/salz3/kh.chao/PR_liftoff_protein_search/data/MANE_RefSeq/MANE.GRCh38.v1.2.refseq_genomic.cleaned.gff_db --liftoff /ccb/salz2/kh.chao/Lifton/results/CHM13_MANE/liftoff.gff3_db --miniprot /ccb/salz2/kh.chao/Lifton/results/CHM13_MANE/miniprot.gff3_db -o results/CHM13_MANE/CHM13_MANE_lifton.gff3 /ccb/salz3/kh.chao/ref_genome/homo_sapiens/T2T-CHM13/chm13v2.0.fa /ccb/salz2/kh.chao/PR_liftoff_protein_search/data/NCBI_Refseq_chr_fixed/GCF_000001405.40_GRCh38.p14_genomic.fna


