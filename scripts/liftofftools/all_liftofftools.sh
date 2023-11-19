# Define a vector of strings
# my_strings=("human_chess_test" "human_man" "mouse" "rice" "bee" "arabadop" "human_to_chimp" "mouse_to_rat" "drosophila")

my_strings=("human_refseq_test" "human_chess_test" "human_mane_test" "mouse_test" "rice_test" "bee_test" "arabadop_test" "human_to_chimp_test" "mouse_to_rat_test" "drosophila_test" "drosophila_erecta_test" "human_mane_to_mouse_test" "human_refseq_to_mouse_test")

# Iterate over the vector
for fruit in "${my_strings[@]}"; do
    echo "./RUN_LIFTOFFTOOLS.sh $fruit"
    ./RUN_LIFTOFFTOOLS.sh $fruit
done
