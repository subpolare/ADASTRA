from tqdm import tqdm
import pandas as pd
import requests
import os

def get_protein_name(uniprot_id):
    url = f'https://rest.uniprot.org/uniprotkb/{uniprot_id}.json'
    response = requests.get(url)
    if response.status_code == 200:
        data = response.json()
        return data.get('genes', [])[0].get('geneName', {}).get('value', '')
    else:
        return 'Unknown'

df = pd.read_csv('/home/subpolare/gtrd/meta/meta_6_may.tsv', delimiter = '\t')
for file in tqdm(os.listdir('/sandbox/subpolare/gtrd/allele_counts/')): 
    if file.startswith('ALIGNS'): 
        algn_id = file.split('_')[0]
        match_row = df[df['algn_id'] == algn_id]
        
        try: 
            if not match_row.empty:
                uniprot_id = match_row['tf_uniprot_id'].values[0]
                name = get_protein_name(uniprot_id).split(' ')[0]
                smth, _, _ = file.split('.')
                n1, n2 = smth.split('_')
                os.rename(f'/sandbox/subpolare/gtrd/allele_counts/{file}', f'/sandbox/subpolare/gtrd/allele_counts/{name}_{n2}_{n1}.vcf.gz')
            else:
                uniprot_id = 'UnknownTF'
                smth, _, _ = file.split('.')
                n1, n2 = smth.split('_')
                os.rename(f'/sandbox/subpolare/gtrd/allele_counts/{file}', f'/sandbox/subpolare/gtrd/allele_counts/UnknownTF_{n2}_{n1}.vcf.gz')
        except FileNotFoundError:
            smth, _, _ = file.split('.')
            n1, n2 = smth.split('_')
            os.rename(f'/sandbox/subpolare/gtrd/allele_counts/{file}', f'/sandbox/subpolare/gtrd/allele_counts/UniProtID_{uniprot_id}_{n2}_{n1}.vcf.gz')
