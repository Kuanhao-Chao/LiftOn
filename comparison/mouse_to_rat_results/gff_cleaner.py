import numpy as np 
import pandas as pd 
import gffutils
from gffutils import gffwriter
from progress.bar import Bar
import sqlite3
import os

def create_database(gff_filepath, db_filepath, overwrite=False, change_miniprot_id=False):
    
    def get_orig_id(feature): 
        try:
            return feature["Target"][0].split()[0] # gffutils gets confused by block numbers at end      
        except:
            return feature["Parent"][0] # throwaway
    
    if not os.path.exists(db_filepath) or overwrite:
        try:
            if change_miniprot_id:
                # Use the standard IDs instead of miniprot's ID
                db = gffutils.create_db(gff_filepath, dbfn=db_filepath, id_spec=get_orig_id, force=True, keep_order=True, merge_strategy='create_unique', sort_attribute_values=True, verbose=True)
            else: 
                # Create the GFF database
                db = gffutils.create_db(gff_filepath, dbfn=db_filepath, force=True, keep_order=True, merge_strategy='create_unique', sort_attribute_values=True, verbose=True)
        except sqlite3.IntegrityError as e:
            print(f"Error during database creation: {e}.\nRetrying with skipping merge strategy.")
            # Use warning strategy to ignore all duplicate IDs during creation
            db = gffutils.create_db(gff_filepath, dbfn=db_filepath, force=True, keep_order=True, merge_strategy='warning', sort_attribute_values=True, verbose=True)
        
    else:
        try:
            # Load the GFF database
            db = gffutils.FeatureDB(db_filepath, keep_order=True)
        except TypeError as e:
            print(f"Error accessing database: {e}.\nRecreating database.")
            # Force update the database
            db = create_database(gff_filepath, db_filepath, overwrite=True)
    return db
    
def write_to_gff(db, gff_filepath, overwrite=False):
    if not os.path.exists(gff_filepath) or overwrite:
        print(f'Writing to gff at {gff_filepath}')
        gff_out = gffwriter.GFFWriter(gff_filepath)
        pbar = Bar('Writing genes...', max=db.count_features_of_type('gene'))
        for gene_rec in db.all_features(featuretype="gene"):
            gff_out.write_gene_recs(db, gene_rec.id)
            pbar.next()
        pbar.finish() 
        gff_out.close()

    else: 
        print(f'File exists at {gff_filepath}, skipping.')

def write_no_gene(db, gff_filepath, overwrite=False):
    if not os.path.exists(gff_filepath) or overwrite:
        print(f'Writing a gff version without gene-level transcripts at {gff_filepath}')
        gff_out = gffwriter.GFFWriter(gff_filepath)
        pbar = Bar('Writing features...', max=db.count_features_of_type(None))
        for feature in db.all_features():
            if feature.featuretype == 'gene':
                continue
            gff_out.write_rec(feature)
            pbar.next()
        pbar.finish() 
        gff_out.close()

    else: 
        print(f'File exists at {gff_filepath}, skipping.')

def write_pc_only(db, gff_filepath, overwrite=False):
    if not os.path.exists(gff_filepath) or overwrite:
        print(f'Writing a gff version with only protein-coding transcripts at {gff_filepath}')
        gff_out = gffwriter.GFFWriter(gff_filepath)
        pbar = Bar('Writing features...', max=db.count_features_of_type(None))
        for feature in db.all_features():
            if feature.featuretype == 'mRNA':
                gff_out.write_rec(feature) # Write mRNA
                gff_out.write_mRNA_children(db, feature.id) # Write children of mRNA in order
            pbar.next()
        pbar.finish() 
        gff_out.close()

    else: 
        print(f'File exists at {gff_filepath}, skipping.')  

def main():

    METHODS = ['liftoff', 'miniprot', 'lifton']

    for method in METHODS:
        print(f'Currently executing gff_cleaner for {method}...')
        # Creating sqlite3 database with gffutils
        gff_filepath = f'./lifton_output/{method}/{method}.sorted.gff3'
        db_filepath = f'./lifton_output/{method}/{method}.db'
        db = create_database(gff_filepath, db_filepath)
        feature_counts = {feature: db.count_features_of_type(feature) for feature in db.featuretypes()}
        print(f'Database contains the following feature types with counts: {feature_counts}')

        # Cleaning out the long transcripts (gene-level)
        new_gff_filepath = os.path.splitext(gff_filepath)[0] + '.cleaned.gff3'
        no_gene_filepath = os.path.splitext(gff_filepath)[0] + '.nogene.gff3'
        pc_only_filepath = os.path.splitext(gff_filepath)[0] + '.pconly.gff3'
        write_to_gff(db, new_gff_filepath) 
        write_no_gene(db, no_gene_filepath)
        write_pc_only(db, pc_only_filepath)


if __name__ == '__main__':
    # main()

    db = create_database('./lifton_output/miniprot/miniprot.sorted.gff3', './lifton_output/miniprot_target/miniprot_target.db', overwrite=True, change_miniprot_id=True)