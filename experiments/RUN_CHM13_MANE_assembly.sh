TARGET=$1
EXPERIMENT=$2

if [ "$TARGET" == "Han1" ]; then
    echo "running liftoff on Han1 genome"
    TARGET_genome="/ccb/salz3/kh.chao/PR_liftoff_protein_search/data/Han1/v1.0/Assembly/Han1_v1.2.fasta"
elif [ "$TARGET" == "Ash1" ]; then
    echo "running liftoff on Ash1 genome"
    TARGET_genome="/ccb/salz3/kh.chao/PR_liftoff_protein_search/data/Ash1/v2.2/Assembly/Ash1_v2.2.fa"
elif [ "$TARGET" == "PR1" ]; then
    echo "running liftoff on PR1 genome"
    TARGET_genome="/ccb/salz3/kh.chao/PR_liftoff_protein_search/data/PR1/v3.0/Assembly/PR1.fa"
else
    echo "Cannot find the genome. Exit.."
    exit -1
fi


python Step_1_gffread_get_liftoff_protein.py $TARGET

if [ "$EXPERIMENT" == "liftoff" ]; then
    ################################
    # Liftoff & Lifton comparison
    ###############################
    echo python liftoff_lifton_comparison/Step_1_parasail_comparison.py $TARGET
    python liftoff_lifton_comparison/Step_1_parasail_comparison.py $TARGET

    echo python liftoff_lifton_comparison/Step_2_grep_transcript_coords.py $TARGET
    python liftoff_lifton_comparison/Step_2_grep_transcript_coords.py $TARGET

    echo python liftoff_lifton_comparison/Step_3_split_liftoff_lifton_better.py $TARGET
    python liftoff_lifton_comparison/Step_3_split_liftoff_lifton_better.py $TARGET

    echo python liftoff_lifton_comparison/Step_4_plot_parasail_liftoff_scores.py $TARGET
    python liftoff_lifton_comparison/Step_4_plot_parasail_liftoff_scores.py $TARGET

    echo python liftoff_lifton_comparison/Step_5_plot_frequency_plot.py $TARGET
    python liftoff_lifton_comparison/Step_5_plot_frequency_plot.py $TARGET

elif [ "$EXPERIMENT" == "miniprot" ]; then
    ################################
    # miniprot & Lifton comparison
    ################################
    echo python miniprot_lifton_comparison/Step_1_parasail_comparison.py $TARGET
    python miniprot_lifton_comparison/Step_1_parasail_comparison.py $TARGET

    echo python miniprot_lifton_comparison/Step_2_grep_transcript_coords.py $TARGET
    python miniprot_lifton_comparison/Step_2_grep_transcript_coords.py $TARGET

    echo python miniprot_lifton_comparison/Step_3_split_liftoff_miniprot_better.py $TARGET
    python miniprot_lifton_comparison/Step_3_split_liftoff_miniprot_better.py $TARGET

    echo python miniprot_lifton_comparison/Step_4_plot_parasail_liftoff_scores.py $TARGET
    python miniprot_lifton_comparison/Step_4_plot_parasail_liftoff_scores.py $TARGET

    echo python miniprot_lifton_comparison/Step_5_plot_frequency_plot.py $TARGET
    python miniprot_lifton_comparison/Step_5_plot_frequency_plot.py $TARGET
fi