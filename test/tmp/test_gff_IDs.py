import gffutils

def check_gff_id(gff_file):
    with open(gff_file, 'r') as file:
        for line in file:
            # Split the line into columns using tab as the delimiter
            columns = line.split('\t')

            # Get the last column
            last_column = columns[-1]

            # Count the occurrences of 'ID=' in the last column
            # id_count = last_column.count('ID=')

            # # Check if there is more than one 'ID=' in the last column
            # if id_count > 1:
            #     print(line.strip())

            counter = 0 
            key_val_pair = last_column.split(";")
            for key_val in key_val_pair:
                if key_val.startswith("ID="):
                    counter += 1
                    # print(key_val)

            if counter > 1:
                print(line)

if __name__ == "__main__":
    gff_file_path = "/ccb/salz2/kh.chao/Lifton/results/drosophila/liftoff/liftoff.gff3"
    check_gff_id(gff_file_path)