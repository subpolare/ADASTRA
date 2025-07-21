import pandas as pd
import sys

def print_algn_ids(file_path, target_tf_uniprot_id):
    df = pd.read_csv(file_path, sep='\t')
    filtered_df = df[df['tf_uniprot_id'] == target_tf_uniprot_id]
    for algn_id in filtered_df['algn_id']:
        print(algn_id)

file_path = '/home/subpolare/gtrd/meta_6_may.tsv'
target_tf_uniprot_id = sys.argv[1]
print_algn_ids(file_path, target_tf_uniprot_id)
