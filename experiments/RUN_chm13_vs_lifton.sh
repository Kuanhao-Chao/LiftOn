EXPERIMENT=$1

python Step_1_gffread_get_liftoff_protein.py CHM13_RefSeq

################################
# CHM13 & Lifton comparison
################################
python CHM13_lifton_comparison copy/Step_1_parasail_comparison.py CHM13_RefSeq
python CHM13_lifton_comparison copy/Step_2_grep_transcript_coords.py CHM13_RefSeq
python CHM13_lifton_comparison copy/Step_3_split_liftoff_lifton_better.py CHM13_RefSeq

python CHM13_lifton_comparison copy/Step_4_plot_parasail_liftoff_scores.py CHM13_RefSeq
python CHM13_lifton_comparison copy/Step_5_plot_frequency_plot.py CHM13_RefSeq