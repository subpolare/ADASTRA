import pandas as pd
import sys

def print_top_tf_uniprot_ids(file_path, top_n=10):
    df = pd.read_csv(file_path, sep='\t')
    top_tf_uniprot_ids = df['tf_uniprot_id'].value_counts().head(top_n)
    print(f"Top {top_n} most frequent tf_uniprot_id:")
    for tf_uniprot_id, count in top_tf_uniprot_ids.items():
        print(f"{tf_uniprot_id}: {count}")

if __name__ == "__main__":
    file_path = '/home/subpolare/gtrd/meta_6_may.tsv'
    print_top_tf_uniprot_ids(file_path)
