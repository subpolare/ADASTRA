import os
import pandas as pd
import requests
import argparse

UNIPROT_URL = 'https://rest.uniprot.org/uniprotkb/{}.json'

def fetch_uniprot_identifier(tf_uniprot_id):
    url = UNIPROT_URL.format(tf_uniprot_id)
    response = requests.get(url)
    
    if response.status_code == 200:
        data = response.json()
        protein_name = data.get('uniProtkbId', None)
        return protein_name
    else:
        print(f'Error: UniProt request failed for {tf_uniprot_id}, status code {response.status_code}')
        return None

def rename_file(input_filepath, meta_filename):
    input_filepath = os.path.abspath(input_filepath)
    directory, filename = os.path.split(input_filepath)
    
    parts = filename.split('_')
    if len(parts) < 3 or not parts[2].startswith('ALIGNS'):
        print(f'Error: Unable to extract algn_id from filename {filename}')
        return
    
    algn_id = parts[2].split('.')[0]  
    
    meta_df = pd.read_csv(meta_filename, sep = '\t')
    
    row = meta_df[meta_df['algn_id'] == algn_id]
    if row.empty:
        print(f'Error: algn_id {algn_id} not found in metadata')
        return
    
    tf_uniprot_id = row.iloc[0]['tf_uniprot_id']  
    
    uniprot_name = fetch_uniprot_identifier(tf_uniprot_id)
    if not uniprot_name:
        print(f'Error: Unable to retrieve UniProt ID for {tf_uniprot_id}')
        return
    
    new_prefix = uniprot_name.split('_')[0]  
    new_filename = f'{new_prefix}_{"_".join(parts[1:])}'  
    new_filepath = os.path.join(directory, new_filename)
    
    os.rename(input_filepath, new_filepath)
    print(f'File {input_filepath} renamed to {new_filepath}')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = 'Script for renaming files based on UniProt data.')
    parser.add_argument('--input', required = True, help = 'Absolute path to input file')
    parser.add_argument('--meta', required = True, help = 'Metadata file')

    args = parser.parse_args()
    rename_file(args.input, args.meta)
