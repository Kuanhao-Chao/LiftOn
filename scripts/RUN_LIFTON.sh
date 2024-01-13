TARGET=$1
EVALUATION=$2
if [[ "$TARGET" == "human_to_chimp" || "$TARGET" == "human_to_chimp_test" ]]; then
    echo "running LiftOn on human_to_chimp"
    REFERENCE_gff="/ccb/salz3/kh.chao/PR_liftoff_protein_search/data/NCBI_Refseq_chr_fixed/rRNA_removed/NCBI_RefSeq_no_rRNA.gff_db"
    # LIFTOFF_gff="/ccb/salz2/jheinz3/shared/lifton/cross_species/human_to_chimp/NHGRI_mPanTro3_liftoff_from_GRCh38_no_alt.gff_polished_db"
    # MINIPROT_gff="/ccb/salz2/jheinz3/shared/lifton/cross_species/human_to_chimp/NHGRI_mPanTro3_miniprot_from_GRCh38.gff_db"
    TARGET_genome="/ccb/salz2/jheinz3/shared/lifton/cross_species/human_to_chimp/NHGRI_mPanTro3-v1.1.fna"
    REFERENCE_genome="/ccb/salz2/kh.chao/PR_liftoff_protein_search/data/NCBI_Refseq_chr_fixed/GCF_000001405.40_GRCh38.p14_genomic.fna"

elif [[ "$TARGET" == "human_mane_to_mouse" || "$TARGET" == "human_mane_to_mouse_test" ]]; then
    echo "running LiftOn on human_mane_to_mouse"
    REFERENCE_gff="/ccb/salz3/kh.chao/PR_liftoff_protein_search/data/MANE_RefSeq/MANE.GRCh38.v1.2.refseq_genomic.cleaned.gff_db"
    # LIFTOFF_gff="/ccb/salz2/jheinz3/shared/lifton/cross_species/human_to_chimp/NHGRI_mPanTro3_liftoff_from_GRCh38_no_alt.gff_polished_db"
    # MINIPROT_gff="/ccb/salz2/jheinz3/shared/lifton/cross_species/human_to_chimp/NHGRI_mPanTro3_miniprot_from_GRCh38.gff_db"
    TARGET_genome="/ccb/salz2/jheinz3/shared/lifton/mouse/GRCm39_genomic.fna"
    REFERENCE_genome="/ccb/salz2/kh.chao/PR_liftoff_protein_search/data/NCBI_Refseq_chr_fixed/GCF_000001405.40_GRCh38.p14_genomic.fna"

elif [[ "$TARGET" == "human_refseq_to_mouse" || "$TARGET" == "human_refseq_to_mouse_test" ]]; then
    echo "running LiftOn on human_refseq_to_mouse"
    REFERENCE_gff="/ccb/salz2/jheinz3/shared/lifton/cross_species/human_to_chimp/GRCh38.p14_refseq_genomic.gff_db"
    # LIFTOFF_gff="/ccb/salz2/jheinz3/shared/lifton/cross_species/human_to_chimp/NHGRI_mPanTro3_liftoff_from_GRCh38_no_alt.gff_polished_db"
    # MINIPROT_gff="/ccb/salz2/jheinz3/shared/lifton/cross_species/human_to_chimp/NHGRI_mPanTro3_miniprot_from_GRCh38.gff_db"
    TARGET_genome="/ccb/salz2/jheinz3/shared/lifton/mouse/GRCm39_genomic.fna"
    REFERENCE_genome="/ccb/salz2/jheinz3/shared/lifton/cross_species/human_to_chimp/GRCh38.p14_refseq_genomic.fna"

elif [[ "$TARGET" == "mouse_to_rat" || "$TARGET" == "mouse_to_rat_test" ]]; then
    echo "running LiftOn on mouse_to_rat"
    REFERENCE_gff="/ccb/salz2/jheinz3/shared/lifton/mouse/GRCm39_genomic.gff_db"
    LIFTOFF_gff="/ccb/salz2/jheinz3/shared/lifton/cross_species/mouse_to_rat/mRatBN7.2_liftoff_from_GRCm39_no_alt.gff_polished_db"
    MINIPROT_gff="/ccb/salz2/jheinz3/shared/lifton/cross_species/mouse_to_rat/mRatBN7.2_miniprot_from_GRCm39.gff_db"
    TARGET_genome="/ccb/salz2/jheinz3/shared/lifton/cross_species/mouse_to_rat/mRatBN7.2_genomic.fna"
    REFERENCE_genome="/ccb/salz2/jheinz3/shared/lifton/mouse/GRCm39_genomic.fna"

elif [[ "$TARGET" == "yeast" || "$TARGET" == "yeast_test" ]]; then
    REFERENCE_gff="/ccb/salz2/kh.chao/LiftOn/data/yeast/GCA_000146045.2_R64_genomic.gff.gz_db"
    LIFTOFF_gff="/ccb/salz2/kh.chao/LiftOn/results/yeast/liftoff/liftoff.gff3_db"
    MINIPROT_gff="/ccb/salz2/kh.chao/LiftOn/results/yeast/miniprot/miniprot.gff3_db"
    TARGET_genome="/ccb/salz2/kh.chao/LiftOn/data/yeast/S288C_reference_genome_R64-4-1_20230830/S288C_reference_sequence_R64-4-1_20230830.fsa"
    REFERENCE_genome="/ccb/salz2/kh.chao/LiftOn/data/yeast/GCA_000146045.2_R64_genomic.fna.gz"

elif [[ "$TARGET" == "arabadop" || "$TARGET" == "arabadop_test" ]]; then
    echo "running LiftOn on arabadop"
    REFERENCE_gff="/ccb/salz2/jheinz3/shared/lifton/arabadop/TAIR10.gff_db"
    LIFTOFF_gff="/ccb/salz2/jheinz3/shared/lifton/arabadop/ASM2311539v1_liftoff_from_TAIR10.gff"
    MINIPROT_gff="/ccb/salz2/jheinz3/shared/lifton/arabadop/ASM2311539v1_miniprot_from_TAIR10.gff"
    TARGET_genome="/ccb/salz2/jheinz3/shared/lifton/arabadop/ASM2311539v1_genomic.fna"
    REFERENCE_genome="/ccb/salz2/jheinz3/shared/lifton/arabadop/TAIR10.fna"

elif [[ "$TARGET" == "bee" || "$TARGET" == "bee_test" ]]; then
    echo "running LiftOn on bee"
    REFERENCE_gff="/ccb/salz2/jheinz3/shared/lifton/bee/HAv3.1_genomic.gff_db"
    LIFTOFF_gff="/ccb/salz2/jheinz3/shared/lifton/bee/ASM1932182v1_liftoff_from_HAv3.1_no_alts.gff_polished_db"
    MINIPROT_gff="/ccb/salz2/jheinz3/shared/lifton/bee/ASM1932182v1_from_miniprot_HAv3.1.gff_db"
    TARGET_genome="/ccb/salz2/jheinz3/shared/lifton/bee/ASM1932182v1_genomic.fna"
    REFERENCE_genome="/ccb/salz2/jheinz3/shared/lifton/bee/HAv3.1_genomic.fna"

elif [[ "$TARGET" == "mouse" || "$TARGET" == "mouse_test" ]]; then
    echo "running LiftOn on mouse"
    REFERENCE_gff="/ccb/salz2/jheinz3/shared/lifton/mouse/GRCm39_genomic.gff_db"
    #LIFTOFF_gff="/ccb/salz2/jheinz3/shared/lifton/mouse/mMusMuc1.1_liftoff_from_GRCm39_no_alt.gff_polished_db"
    #MINIPROT_gff="/ccb/salz2/jheinz3/shared/lifton/mouse/mMusMuc1.1_from_miniprot_GRCm39.gff_db"
    TARGET_genome="/ccb/salz2/jheinz3/shared/lifton/mouse/NOD_SCID_genomic.fna"
    REFERENCE_genome="/ccb/salz2/jheinz3/shared/lifton/mouse/GRCm39_genomic.fna"

elif [[ "$TARGET" == "rice" || "$TARGET" == "rice_test" ]]; then
    echo "running LiftOn on rice"
    REFERENCE_gff="/ccb/salz2/jheinz3/shared/lifton/rice/IRGSP_genomic.gff_db"
    #LIFTOFF_gff="/ccb/salz2/jheinz3/shared/lifton/rice/ASM2616768v1_liftoff_from_IRGSP_no_alts.gff_polished"
    #MINIPROT_gff="/ccb/salz2/jheinz3/shared/lifton/rice/ASM2616768v1_from_miniprot_IRGSP.gff"
    TARGET_genome="/ccb/salz2/jheinz3/shared/lifton/rice/ASM3414082v1_genomic.fna"
    REFERENCE_genome="/ccb/salz2/jheinz3/shared/lifton/rice/IRGSP_genomic.fna"

elif [[ "$TARGET" == "human_chess" || "$TARGET" == "human_chess_test" ]]; then
    echo "running LiftOn on human_chess"
    REFERENCE_gff="/ccb/salz3/kh.chao/ref_genome/homo_sapiens/chess/chess3.0.1.gff_db"
    # LIFTOFF_gff="/ccb/salz2/jheinz3/shared/lifton/rice/ASM2616768v1_liftoff_from_IRGSP_no_alts.gff_polished_db"
    # MINIPROT_gff="/ccb/salz2/jheinz3/shared/lifton/rice/ASM2616768v1_from_miniprot_IRGSP.gff_db"
    TARGET_genome="/ccb/salz3/kh.chao/ref_genome/homo_sapiens/T2T-CHM13/chm13v2.0.fa"
    REFERENCE_genome="/ccb/salz3/kh.chao/ref_genome/homo_sapiens/chess/hg38_p12_ucsc.no_alts.no_fixs.fa"
    ADDITIONAL_ARG="-f /ccb/salz2/kh.chao/LiftOn/data/features_chess.txt"

elif [[ "$TARGET" == "human_mane" || "$TARGET" == "human_mane_test" ]]; then
    echo "running LiftOn on human_mane"
    REFERENCE_gff="/ccb/salz3/kh.chao/PR_liftoff_protein_search/data/MANE_RefSeq/MANE.GRCh38.v1.2.refseq_genomic.cleaned.gff_db"
    # LIFTOFF_gff="/ccb/salz2/jheinz3/shared/lifton/rice/ASM2616768v1_liftoff_from_IRGSP_no_alts.gff_polished_db"
    # MINIPROT_gff="/ccb/salz2/jheinz3/shared/lifton/rice/ASM2616768v1_from_miniprot_IRGSP.gff_db"
    TARGET_genome="/ccb/salz3/kh.chao/ref_genome/homo_sapiens/T2T-CHM13/chm13v2.0.fa"
    REFERENCE_genome="/ccb/salz2/kh.chao/PR_liftoff_protein_search/data/NCBI_Refseq_chr_fixed/GCF_000001405.40_GRCh38.p14_genomic.fna"

elif [[ "$TARGET" == "human_refseq" || "$TARGET" == "human_refseq_test" ]]; then
    echo "running LiftOn on human_refseq"
    REFERENCE_gff="/ccb/salz3/kh.chao/PR_liftoff_protein_search/data/NCBI_Refseq_chr_fixed/rRNA_removed/NCBI_RefSeq_no_rRNA.gff_db"
    # LIFTOFF_gff="/ccb/salz2/jheinz3/shared/lifton/rice/ASM2616768v1_liftoff_from_IRGSP_no_alts.gff_polished_db"
    # MINIPROT_gff="/ccb/salz2/jheinz3/shared/lifton/rice/ASM2616768v1_from_miniprot_IRGSP.gff_db"
    TARGET_genome="/ccb/salz3/kh.chao/ref_genome/homo_sapiens/T2T-CHM13/chm13v2.0.fa"
    REFERENCE_genome="/ccb/salz2/kh.chao/PR_liftoff_protein_search/data/NCBI_Refseq_chr_fixed/GCF_000001405.40_GRCh38.p14_genomic.fna"
elif [[ "$TARGET" == "drosophila" || "$TARGET" == "drosophila_test" ]]; then
    echo "running drosophila"
    REFERENCE_gff="/ccb/salz2/jheinz3/shared/lifton/cross_species/drosophila/d.menogaster_genomic.gff"
    TARGET_genome="/ccb/salz2/jheinz3/shared/lifton/cross_species/drosophila/d.simulans_genomic.fna"
    REFERENCE_genome="/ccb/salz2/jheinz3/shared/lifton/cross_species/drosophila/d.menogaster_genomic.fna"

elif [[ "$TARGET" == "drosophila_erecta" || "$TARGET" == "drosophila_erecta_test" ]]; then
    echo "running drosophila erecta"
    REFERENCE_gff="/ccb/salz2/jheinz3/shared/lifton/cross_species/drosophila/d.menogaster_genomic.gff_db"
    TARGET_genome="/ccb/salz2/jheinz3/shared/lifton/cross_species/drosophila/d.erecta_genomic.fna"
    REFERENCE_genome="/ccb/salz2/jheinz3/shared/lifton/cross_species/drosophila/d.menogaster_genomic.fna"

else
    echo "Cannot find the genome. Exit.."
    exit -1
fi

LIFTOFF_gff="/ccb/salz2/kh.chao/LiftOn/results/$TARGET/lifton_output/liftoff/liftoff.gff3_db"
MINIPROT_gff="/ccb/salz2/kh.chao/LiftOn/results/$TARGET/lifton_output/miniprot/miniprot.gff3"

intermediate_dir="/ccb/salz2/kh.chao/LiftOn/results/$TARGET/lifton_output/intermediate_files/"

output_LIFTON_gff="/ccb/salz2/kh.chao/LiftOn/results/$TARGET/lifton.gff3"

ref_proteins="${intermediate_dir}proteins.fa"
ref_trans="${intermediate_dir}transcripts.fa"
log_file="/ccb/salz2/kh.chao/LiftOn/results/$TARGET/output.log"

echo "lifoff annotation: $LIFTOFF_gff"
echo "miniprot annotation: $MINIPROT_gff"
echo "intermediate_dir: $intermediate_dir"
echo "lifton -D -g $REFERENCE_gff -dir $intermediate_dir -o $output_LIFTON_gff --liftoff $LIFTOFF_gff --miniprot $MINIPROT_gff --proteins $ref_proteins --transcripts $ref_trans -copies $TARGET_genome $REFERENCE_genome $ADDITIONAL_ARG $EVALUATION"

echo lifton -D -g $REFERENCE_gff -o $output_LIFTON_gff --liftoff $LIFTOFF_gff --miniprot $MINIPROT_gff --proteins $ref_proteins --transcripts $ref_trans -copies $TARGET_genome $REFERENCE_genome $ADDITIONAL_ARG $EVALUATION
