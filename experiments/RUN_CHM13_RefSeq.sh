EXPERIMENT=$1

python Step_1_gffread_get_liftoff_protein.py CHM13_RefSeq

if [ "$EXPERIMENT" == "liftoff" ]; then

    ################################
    # Liftoff & Lifton comparison
    ################################
    python liftoff_lifton_comparison/Step_1_parasail_comparison.py CHM13_RefSeq
    python liftoff_lifton_comparison/Step_2_grep_transcript_coords.py CHM13_RefSeq
    python liftoff_lifton_comparison/Step_3_split_liftoff_lifton_better.py CHM13_RefSeq

    python liftoff_lifton_comparison/Step_4_plot_parasail_liftoff_scores.py CHM13_RefSeq
    python liftoff_lifton_comparison/Step_5_plot_frequency_plot.py CHM13_RefSeq

elif [ "$EXPERIMENT" == "miniprot" ]; then

    ################################
    # miniprot & Lifton comparison
    ################################
    python miniprot_lifton_comparison/Step_1_parasail_comparison.py CHM13_RefSeq
    python miniprot_lifton_comparison/Step_2_grep_transcript_coords.py CHM13_RefSeq
    python miniprot_lifton_comparison/Step_3_split_liftoff_miniprot_better.py CHM13_RefSeq
    python miniprot_lifton_comparison/Step_4_plot_parasail_liftoff_scores.py CHM13_RefSeq
    python miniprot_lifton_comparison/Step_5_plot_frequency_plot.py CHM13_RefSeq
fi