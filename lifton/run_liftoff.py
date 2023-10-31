import liftoff
import subprocess

def check_liftoff_installed():
    installed = False
    try:
        installed = liftoff.__version__ >= "1.6.3"
    except:
        pass
    return installed


def run_liftoff():

    # liftoff -g /ccb/salz2/kh.chao/PR_liftoff_protein_search/data/MANE_RefSeq/MANE.GRCh38.v1.2.refseq_genomic.cleaned.gff -o /ccb/salz2/kh.chao/PR_liftoff_protein_search/results/liftoff/$TARGET/${TARGET}_liftoff.gff3 -chroms /ccb/salz2/kh.chao/PR_liftoff_protein_search/data/chroms_mapping.csv -copies -polish $TARGET_genome /ccb/salz2/kh.chao/PR_liftoff_protein_search/data/NCBI_Refseq_chr_fixed/GCF_000001405.40_GRCh38.p14_genomic.fna

    # # command = ["liftoff", "-g", arg,  "--version"]

    # installed = False
    # try:
    #     res = subprocess.run(command)
    #     print(res)
    #     installed = True
    # except: 
    #     pass

    print(">> run_liftoff")
    pass