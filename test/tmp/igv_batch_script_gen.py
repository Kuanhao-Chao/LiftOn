import os 
def main():

    os.makedirs("images", exist_ok=True)
    fw = open("igv_batch.bat", "w")
    fw.write("new\n")
    fw.write("genome chm13v2.0\n")
    fw.write("snapshotDirectory /Users/chaokuan-hao/Documents/Projects/IGV-snapshot-automator/IGV_Snapshots\n")

    fw.write("load /Users/chaokuan-hao/Documents/Projects/PR_Liftoff/data/MANE__GRCh38_2_CHM13/CHM13_MANE.sort.sorted.gff3\n")

    fw.write("expand CHM13_MANE.sort.sorted.gff3\n")

    fw.write("load /Users/chaokuan-hao/Documents/Projects/PR_Liftoff/data/MANE__GRCh38_2_CHM13/CHM13_MANE_miniprot.fix.sorted.gff\n")

    fw.write("expand CHM13_MANE_miniprot.fix.sorted.gff\n")

    fw.write("load /Users/chaokuan-hao/Documents/Projects/PR_Liftoff/data/MANE__GRCh38_2_CHM13/Eviann.fna.sorted.gff\n")

    fw.write("expand Eviann.fna.sorted.gff\n")

    fw.write("load /Users/chaokuan-hao/Documents/Projects/PR_Liftoff/data/MANE__GRCh38_2_CHM13/CHM13_MANE.sort.sorted.gff3\n")

    fw.write("expand CHM13_MANE.sort.sorted.gff3\n")

    counter = 0
    with open("identities_coords_truncated.txt", "r") as fr:
        lines = fr.read().splitlines()
        for line in lines:
            eles = line.split("\t")
            print(eles[4])

            fw.write(f"goto {eles[4]}\n")
            fw.write(f"snapshot images/{eles[4]}.png\n")
            counter += 1

            # if counter > 10:
            #     break

    fw.write("exit\n")

if __name__ == "__main__":
    main()