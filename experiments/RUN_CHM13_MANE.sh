EXPERIMENT=$1

python Step_1_gffread_get_liftoff_protein.py CHM13_MANE

if [ "$EXPERIMENT" == "liftoff" ]; then
    ################################
    # Liftoff & Lifton comparison
    ###############################
    python liftoff_lifton_comparison/Step_1_parasail_comparison.py CHM13_MANE
    python liftoff_lifton_comparison/Step_2_grep_transcript_coords.py CHM13_MANE
    python liftoff_lifton_comparison/Step_3_split_liftoff_lifton_better.py CHM13_MANE

    python liftoff_lifton_comparison/Step_4_plot_parasail_liftoff_scores.py CHM13_MANE
    python liftoff_lifton_comparison/Step_5_plot_frequency_plot.py CHM13_MANE

elif [ "$EXPERIMENT" == "miniprot" ]; then
    ################################
    # miniprot & Lifton comparison
    ################################
    python miniprot_lifton_comparison/Step_1_parasail_comparison.py CHM13_MANE
    python miniprot_lifton_comparison/Step_2_grep_transcript_coords.py CHM13_MANE
    python miniprot_lifton_comparison/Step_3_split_liftoff_miniprot_better.py CHM13_MANE
    python miniprot_lifton_comparison/Step_4_plot_parasail_liftoff_scores.py CHM13_MANE
    python miniprot_lifton_comparison/Step_5_plot_frequency_plot.py CHM13_MANE

elif [ "$EXPERIMENT" == "eviann" ]; then
    python eviann_lifton_comparison/Step_1_map_eviann_id_2_lifton_id.py
    python eviann_lifton_comparison/Step_2_get_eviann_protein.py 
    python eviann_lifton_comparison/Step_3_parasail_comparison.py
    python eviann_lifton_comparison/Step_4_grep_transcript_coords.py
    python eviann_lifton_comparison/Step_5_split_liftoff_lifton_better.py
    python eviann_lifton_comparison/Step_6_plot_parasail_liftoff_scores.py CHM13_MANE
    python eviann_lifton_comparison/Step_7_plot_frequency_plot.py CHM13_MANE
fi
