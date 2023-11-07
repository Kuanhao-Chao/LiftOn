TARGET=$1
if [ "$TARGET" == "human_to_chimp" ]; then
    echo "running LiftOn on human_to_chimp"
    REFERENCE_gff="/ccb/salz2/jheinz3/shared/lifton/cross_species/human_to_chimp/GRCh38.p14_refseq_genomic_no_alt.gff"
    LIFTOFF_gff="/ccb/salz2/jheinz3/shared/lifton/cross_species/human_to_chimp/NHGRI_mPanTro3_liftoff_from_GRCh38_no_alt.gff_polished_db"
    MINIPROT_gff="/ccb/salz2/jheinz3/shared/lifton/cross_species/human_to_chimp/NHGRI_mPanTro3_miniprot_from_GRCh38.gff_db"
    TARGET_genome="/ccb/salz2/jheinz3/shared/lifton/cross_species/human_to_chimp/NHGRI_mPanTro3-v1.1.fna"
    REFERENCE_genome="/ccb/salz2/jheinz3/shared/lifton/cross_species/human_to_chimp/GRCh38.p14_refseq_genomic.fna"

elif [ "$TARGET" == "mouse_to_rat" ]; then
    echo "running LiftOn on mouse_to_rat"
    REFERENCE_gff="/ccb/salz2/jheinz3/shared/lifton/mouse/GRCm39_no_alt_genomic.gff"
    LIFTOFF_gff="/ccb/salz2/jheinz3/shared/lifton/cross_species/mouse_to_rat/mRatBN7.2_liftoff_from_GRCm39_no_alt.gff_polished_db"
    MINIPROT_gff="/ccb/salz2/jheinz3/shared/lifton/cross_species/mouse_to_rat/mRatBN7.2_miniprot_from_GRCm39.gff_db"
    TARGET_genome="/ccb/salz2/jheinz3/shared/lifton/cross_species/mouse_to_rat/mRatBN7.2_genomic.fna"
    REFERENCE_genome="/ccb/salz2/jheinz3/shared/lifton/mouse/GRCm39_genomic.fna"

elif [ "$TARGET" == "yeast"  ]; then
    REFERENCE_gff="/ccb/salz2/kh.chao/Lifton/data/yeast/GCA_000146045.2_R64_genomic.gff.gz_db"
    LIFTOFF_gff="/ccb/salz2/kh.chao/Lifton/results/yeast/liftoff/liftoff.gff3_db"
    MINIPROT_gff="/ccb/salz2/kh.chao/Lifton/results/yeast/miniprot/miniprot.gff3_db"
    TARGET_genome="/ccb/salz2/kh.chao/Lifton/data/yeast/S288C_reference_genome_R64-4-1_20230830/S288C_reference_sequence_R64-4-1_20230830.fsa"
    REFERENCE_genome="/ccb/salz2/kh.chao/Lifton/data/yeast/GCA_000146045.2_R64_genomic.fna.gz"

elif [[ "$TARGET" == "arabadop" || "$TARGET" == "arabadop_test" ]]; then
    echo "running LiftOn on arabadop"
    REFERENCE_gff="/ccb/salz2/jheinz3/shared/lifton/arabadop/TAIR10.gff"
    LIFTOFF_gff="/ccb/salz2/jheinz3/shared/lifton/arabadop/Tanz-1_liftoff_from_TAIR10.gff_polished_db"
    MINIPROT_gff="/ccb/salz2/jheinz3/shared/lifton/arabadop/Tanz-1_from_miniprot_TAIR10.gff_db"
    TARGET_genome="/ccb/salz2/jheinz3/shared/lifton/arabadop/Tanz-1.fna"
    REFERENCE_genome="/ccb/salz2/jheinz3/shared/lifton/arabadop/TAIR10.fna"

elif [ "$TARGET" == "bee" ]; then
    echo "running LiftOn on bee"
    REFERENCE_gff="/ccb/salz2/jheinz3/shared/lifton/bee/HAv3.1_genomic_no_alts.gff"
    LIFTOFF_gff="/ccb/salz2/jheinz3/shared/lifton/bee/ASM1932182v1_liftoff_from_HAv3.1_no_alts.gff_polished_db"
    MINIPROT_gff="/ccb/salz2/jheinz3/shared/lifton/bee/ASM1932182v1_from_miniprot_HAv3.1.gff_db"
    TARGET_genome="/ccb/salz2/jheinz3/shared/lifton/bee/ASM1932182v1_genomic.fna"
    REFERENCE_genome="/ccb/salz2/jheinz3/shared/lifton/bee/HAv3.1_genomic.fna"

elif [ "$TARGET" == "bee_test" ]; then
    echo "running LiftOn on bee test"
    REFERENCE_gff="/ccb/salz2/jheinz3/shared/lifton/bee/HAv3.1_genomic_no_alts.gff"
    LIFTOFF_gff="/ccb/salz2/kh.chao/Lifton/results/$TARGET/liftoff/liftoff.gff3_db"
    MINIPROT_gff="/ccb/salz2/kh.chao/Lifton/results/$TARGET/miniprot/miniprot.gff3_db"
    TARGET_genome="/ccb/salz2/jheinz3/shared/lifton/bee/ASM1932182v1_genomic.fna"
    REFERENCE_genome="/ccb/salz2/jheinz3/shared/lifton/bee/HAv3.1_genomic.fna"

elif [[ "$TARGET" == "mouse" || "$TARGET" == "mouse_test" ]]; then
    echo "running LiftOn on mouse"
    REFERENCE_gff="/ccb/salz2/jheinz3/shared/lifton/mouse/GRCm39_no_alt_genomic.gff"
    LIFTOFF_gff="/ccb/salz2/jheinz3/shared/lifton/mouse/mMusMuc1.1_liftoff_from_GRCm39_no_alt.gff_polished_db"
    MINIPROT_gff="/ccb/salz2/jheinz3/shared/lifton/mouse/mMusMuc1.1_from_miniprot_GRCm39.gff_db"
    TARGET_genome="/ccb/salz2/jheinz3/shared/lifton/mouse/mMusMuc1.1_genomic.fna"
    REFERENCE_genome="/ccb/salz2/jheinz3/shared/lifton/mouse/GRCm39_genomic.fna"

elif [[ "$TARGET" == "rice" || "$TARGET" == "rice_test" ]]; then
    echo "running LiftOn on rice"
    REFERENCE_gff="/ccb/salz2/jheinz3/shared/lifton/rice/IRGSP_genomic_no_alts.gff"
    LIFTOFF_gff="/ccb/salz2/jheinz3/shared/lifton/rice/ASM2616768v1_liftoff_from_IRGSP_no_alts.gff_polished_db"
    MINIPROT_gff="/ccb/salz2/jheinz3/shared/lifton/rice/ASM2616768v1_from_miniprot_IRGSP.gff_db"
    TARGET_genome="/ccb/salz2/jheinz3/shared/lifton/rice/ASM2616768v1_genomic.fna"
    REFERENCE_genome="/ccb/salz2/jheinz3/shared/lifton/rice/IRGSP_genomic.fna"

else
    echo "Cannot find the genome. Exit.."
    exit -1
fi

LIFTOFF_gff="/ccb/salz2/kh.chao/Lifton/results/$TARGET/liftoff/liftoff.gff3_db"
MINIPROT_gff="/ccb/salz2/kh.chao/Lifton/results/$TARGET/miniprot/miniprot.gff3_db"

intermediate_dir="/ccb/salz2/kh.chao/Lifton/results/$TARGET/intermediate_files/"
output_LIFTON_gff="/ccb/salz2/kh.chao/Lifton/results/$TARGET/lifton.gff3"
ref_proteins="${intermediate_dir}proteins.fa"
ref_trans="${intermediate_dir}transcripts.fa"

echo "lifoff annotation: $LIFTOFF_gff"
echo "miniprot annotation: $MINIPROT_gff"
echo "intermediate_dir: $intermediate_dir"
echo "lifton -g $REFERENCE_gff -dir $intermediate_dir -o $output_LIFTON_gff --liftoff $LIFTOFF_gff --miniprot $MINIPROT_gff --proteins $ref_proteins --transcripts $ref_trans $TARGET_genome $REFERENCE_genome"

lifton -g $REFERENCE_gff -dir $intermediate_dir -o $output_LIFTON_gff --liftoff $LIFTOFF_gff --miniprot $MINIPROT_gff --proteins $ref_proteins --transcripts $ref_trans $TARGET_genome $REFERENCE_genome
