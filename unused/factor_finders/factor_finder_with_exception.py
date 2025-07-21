import pandas as pd
import sys

def print_algn_ids_for_tf_and_exclude_cell_id(file_path, target_tf_uniprot_id, exclude_cell_id):
    df = pd.read_csv(file_path, sep='\t')
    filtered_df = df[(df['tf_uniprot_id'] == target_tf_uniprot_id) & (df['cell_id'] != exclude_cell_id)]
    for algn_id in filtered_df['algn_id']:
        print(algn_id)

if __name__ == "__main__":
    file_path = '/home/subpolare/gtrd/meta_6_may.tsv'
    target_tf_uniprot_id = 'P49711'
    exclude_cell_id = 514
    print_algn_ids_for_tf_and_exclude_cell_id(file_path, target_tf_uniprot_id, exclude_cell_id)
