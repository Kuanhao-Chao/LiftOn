# Define a vector of strings
my_strings=("human_chess" "human_man" "mouse" "rice" "bee" "arabadop" "human_to_chimp" "mouse_to_rat" "drosophila")

# Iterate over the vector
for fruit in "${my_strings[@]}"; do
    echo "./RUN_LIFTOFFTOOLS.sh $fruit"
    ./RUN_LIFTOFFTOOLS.sh $fruit
done
