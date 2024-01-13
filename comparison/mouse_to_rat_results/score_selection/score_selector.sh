#!/bin/bash

## LIFTON BEST

input_file="/home/smao10/lifton/comparison/output/mouse_to_rat/lifton_best.out"
output_file="/home/smao10/lifton/comparison/mouse_to_rat_results/score_selection/best.out"

# lifton score is greater than 0.9 and should improve liftoff by at least 0.25
awk -F'\t' '$5 > 0.9 && $5 > $2 + 0.25 { print }' "$input_file" > "$output_file"
echo "Best lifton regions of interest have been written to $output_file"

# show number of lines filtered
input_lines=$(wc -l < "$input_file")
output_lines=$(wc -l < "$output_file")
echo -e "Input: $input_lines\tOutput: $output_lines\tLines filtered: $((input_lines - output_lines))"

## LIFTON EQ LIFTOFF

input_file="/home/smao10/lifton/comparison/output/mouse_to_rat/lifton_eq_liftoff.out"
output_file="/home/smao10/lifton/comparison/mouse_to_rat_results/score_selection/eq_liftoff.out"

# lifton score is greater than 0.95 and should improve miniprot score by at least 0.2
awk -F'\t' '$5 > 0.95 && $5 > $3 + 0.2 { print }' "$input_file" > "$output_file"
echo "Lifton equal to liftoff regions of interest have been written to $output_file"

# show number of lines filtered
input_lines=$(wc -l < "$input_file")
output_lines=$(wc -l < "$output_file")
echo -e "Input: $input_lines\tOutput: $output_lines\tLines filtered: $((input_lines - output_lines))"


## LIFTON EQ MINIPROT

input_file="/home/smao10/lifton/comparison/output/mouse_to_rat/lifton_eq_miniprot.out"
output_file="/home/smao10/lifton/comparison/mouse_to_rat_results/score_selection/eq_miniprot.out"

# lifton score is greater than 0.95 and should improve liftoff score by at least 0.5
awk -F'\t' '$5 > 0.95 && $5 > $2 + 0.5 { print }' "$input_file" > "$output_file"
echo "Lifton equal to miniprot regions of interest have been written to $output_file"

# show number of lines filtered
input_lines=$(wc -l < "$input_file")
output_lines=$(wc -l < "$output_file")
echo -e "Input: $input_lines\tOutput: $output_lines\tLines filtered: $((input_lines - output_lines))"