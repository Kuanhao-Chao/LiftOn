python setup.py install

#lifton --proteins /ccb/salz3/kh.chao/PR_liftoff_protein_search/data/MANE_RefSeq/MANE.GRCh38.v1.2.refseq_genomic.AA.cleaned.fa -g /ccb/salz3/kh.chao/PR_liftoff_protein_search/data/MANE_RefSeq/MANE.GRCh38.v1.2.refseq_genomic.cleaned.gff_db --miniprot /ccb/salz2/kh.chao/Lifton/results/CHM13_MANE/miniprot.gff3_db -o results/CHM13_MANE/CHM13_MANE_lifton.gff3 /ccb/salz3/kh.chao/ref_genome/homo_sapiens/T2T-CHM13/chm13v2.0.fa /ccb/salz2/kh.chao/PR_liftoff_protein_search/data/NCBI_Refseq_chr_fixed/GCF_000001405.40_GRCh38.p14_genomic.fna

lifton -g /ccb/salz2/kh.chao/Lifton/data/yeast/GCA_000146045.2_R64_genomic.gff.gz -o results/yeast_TEST/lifton.gff3 -dir results/yeast_TEST/intermediate_files/ /ccb/salz2/kh.chao/Lifton/data/yeast/S288C_reference_genome_R64-4-1_20230830/S288C_reference_sequence_R64-4-1_20230830.fsa /ccb/salz2/kh.chao/Lifton/data/yeast/GCA_000146045.2_R64_genomic.fna.gz

#lifton -g /ccb/salz2/kh.chao/Lifton/data/yeast/GCA_000146045.2_R64_genomic.gff.gz --liftoff /ccb/salz2/kh.chao/Lifton/results/yeast_TEST/liftoff/liftoff.gff3_db --miniprot /ccb/salz2/kh.chao/Lifton/data/yeast/S288C_reference_genome_R64-4-1_miniprot.gff3 -o results/yeast_TEST/lifton.gff3 --dir results/yeast_TEST/intermediate_files/ /ccb/salz2/kh.chao/Lifton/data/yeast/S288C_reference_genome_R64-4-1_20230830/S288C_reference_sequence_R64-4-1_20230830.fsa /ccb/salz2/kh.chao/Lifton/data/yeast/GCA_000146045.2_R64_genomic.fna.gz

#lifton --proteins /ccb/salz3/kh.chao/PR_liftoff_protein_search/data/MANE_RefSeq/MANE.GRCh38.v1.2.refseq_genomic.AA.cleaned.fa -g /ccb/salz2/kh.chao/PR_liftoff_protein_search/results/liftoff/CHM13_MANE/exp/test.gff3 --miniprot /ccb/salz2/kh.chao/Lifton/results/CHM13_MANE/miniprot.gff3_db -o results/TEST/lifton.gff3 /ccb/salz3/kh.chao/ref_genome/homo_sapiens/T2T-CHM13/chm13v2.0.fa /ccb/salz2/kh.chao/PR_liftoff_protein_search/data/NCBI_Refseq_chr_fixed/GCF_000001405.40_GRCh38.p14_genomic.fna

