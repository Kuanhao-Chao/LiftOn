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
   
    for target in targets:
        print(f'Extracting for target {target}...')

        feature1 = db1[target]
        print(str(feature1) + '\n')
        feature2 = db2[target]
        print(str(feature2) + '\n')
        feature3 = db3[target]
        print(str(feature3) + '\n')
        


def main():

    liftoff_db = gffutils.FeatureDB('./lifton_output/liftoff/liftoff.db', keep_order=True)
    miniprot_db = gffutils.FeatureDB('./lifton_output/miniprot_target/miniprot_target.db', keep_order=True)
    lifton_db = gffutils.FeatureDB('./lifton_output/lifton/lifton.db', keep_order=True)
    targets = read_targets('./extracted_gff/targets.txt')

    extract_targets(liftoff_db, miniprot_db, lifton_db, targets)


if __name__ == '__main__': 
    main()