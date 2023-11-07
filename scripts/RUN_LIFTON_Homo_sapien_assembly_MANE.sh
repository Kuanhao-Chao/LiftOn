python setup.py install

TARGET=$1
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

mkdir /ccb/salz2/kh.chao/Lifton/results/$TARGET

lifton --proteins /ccb/salz2/kh.chao/PR_liftoff_protein_search/data/protein.fasta --liftoffdb /ccb/salz2/kh.chao/PR_liftoff_protein_search/results/liftoff/${TARGET}/${TARGET}.sort.gff3_db --miniprotdb /ccb/salz2/kh.chao/PR_liftoff_protein_search/results/liftoff/${TARGET}/${TARGET}_miniprot.fix.gff_db -o results/$TARGET/${TARGET}_lifton.gff3 ${TARGET_genome}
