import gffutils

def check_gff_validity(gff_file):
    try:
        # Load GFF file using gffutils
        db = gffutils.create_db(gff_file, ':memory:')

        # Check for any issues
        for feature in db.features_of_type():
            # Check if the ID field has more than one value
            if ',' in feature.id:
                print(f"Error: Feature '{feature.id}' has multiple values in the ID field.")

        print("GFF file is valid.")

    except gffutils.FeatureNotFoundError as e:
        print(f"Error: {e}")
    # except gffutils.FeatureFormatError as e:
    #     print(f"Error: {e}")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")

if __name__ == "__main__":
    gff_file_path = "/ccb/salz2/kh.chao/Lifton/results/drosophila/liftoff/liftoff.gff3"
    check_gff_validity(gff_file_path)