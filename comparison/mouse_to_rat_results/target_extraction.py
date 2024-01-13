import numpy as np 
import pandas as pd 
import gffutils
from gffutils import gffwriter
from progress.bar import Bar
import sqlite3
import os

def read_targets(path):
    with open(path, 'r') as file:
        target_ids = [line.strip().split('\t')[0] for line in file]
        return target_ids

def extract_targets(db1, db2, db3, targets):
    # Create a new GFF file for the extracted transcripts
    output_gff_file = 'combined_extracted_targets.gff'
    with open(output_gff_file, 'w') as output_gff:
        # Iterate through target IDs
        for target_id in targets:
            # Search for the feature with the matching ID in the combined FeatureDB
            feature = combined_db[target_id]
            if feature:
                # Write the GFF line to the new GFF file
                output_gff.write(str(feature) + '\n')

                # Print relevant information about the target
                print(f"Target ID: {target_id}")
                print(f"Chromosome: {feature.seqid}")
                print(f"Start: {feature.start}")
                print(f"End: {feature.end}")
                print(f"Strand: {feature.strand}")
                print(f"Attributes: {feature.attributes}")
                print("-" * 20)

def main():

    liftoff_db = gffutils.FeatureDB('./lifton_output/liftoff/liftoff.db', keep_order=True)
    miniprot_db = gffutils.FeatureDB('./lifton_output/miniprot_target/miniprot_target.db', keep_order=True)
    lifton_db = gffutils.FeatureDB('./lifton_output/lifton/lifton.db', keep_order=True)
    targets = read_targets('./extracted_gff/targets.txt')
    
    print(str(miniprot_db['rna-NM_001305875.1']))

    #extract_targets(liftoff_db, miniprot_db, lifton_db, targets)


if __name__ == '__main__': 
    main()