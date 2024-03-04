
# "human_mane_to_mouse" "yeast"
declare -a arr=("human_refseq" "human_mane" "human_chess" "rice" "arabadop" "bee" "mouse" "human_to_chimp" "drosophila_erecta" "mouse_to_rat")

for i in "${arr[@]}"
do 
    echo "running $i"
    # Cp all figures
    echo cp -r /ccb/salz2/kh.chao/LiftOn/results/$i/lifton_output/visualization/* $i/
    cp -r /ccb/salz2/kh.chao/LiftOn/results/$i/lifton_output/visualization/* $i/
    # rm -rf $i/Liftoff/
    # rm -rf $i/miniprot/
    rm -rf $i/liftoff_frequency.png
    rm -rf $i/lifton_frequency.png
    rm -rf $i/miniprot_frequency.png

    # # Cp gene order plot
    # echo cp -r /ccb/salz2/kh.chao/LiftOn/results/$i/liftofftools_output/gene_order_plot.pdf $i/
    # cp -r /ccb/salz2/kh.chao/LiftOn/results/$i/liftofftools_output/gene_order_plot.pdf $i/
    # # convert -density 200 /ccb/salz2/kh.chao/LiftOn/results/$i/liftofftools_output/gene_order_plot.pdf $i/gene_order_plot.png

    # rm -rf $i/gene_order_plot.pdf

    # convert -density 250 $i/circos_plot.pdf $i/circos_plot.png

done