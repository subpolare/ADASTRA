from warnings import simplefilter
from tqdm import tqdm
import pandas as pd
import numpy as np
import requests

simplefilter(action = 'ignore')

df = pd.read_csv('/home/subpolare/gtrd/processed.tsv', sep = '\t')
meta = pd.read_csv('/home/subpolare/gtrd/meta/meta_6_may.tsv', sep = '\t') 
cells = pd.read_csv('/home/subpolare/gtrd/meta/meta_cells_and_tissues.tsv', sep = '\t') 

for idx in tqdm(range(len(df)), colour = 'green', desc = 'Processing with metadata'):
    align_id = df.at[idx, 'ALIGN']
    meta_row = meta[meta['algn_id'] == align_id]
    
    if not meta_row.empty:
        meta_row = meta_row.iloc[0]
        
        if pd.isna(df.at[idx, 'TF_UNIPROT_AC']):
            df.at[idx, 'TF_UNIPROT_AC'] = meta_row['tf_uniprot_id']

        if pd.isna(df.at[idx, 'GTRD_CELL_TYPE_NAME']):
            cell_id = df.at[idx, 'GTRD_CELL_TYPE_ID']
            cell_row = cells[cells['id'] == cell_id]

            if not cell_row.empty:
                cell_row = cell_row.iloc[0]
                df.at[idx, 'GTRD_CELL_TYPE_NAME'] = cell_row['title']
            else: 
                df.at[idx, 'GTRD_CELL_TYPE_NAME'] = 'NULL'
        
        if pd.isna(df.at[idx, 'GEO_GSE']):
            df.at[idx, 'GEO_GSE'] = meta_row['geo_gsm']

        if pd.isna(df.at[idx, 'GEO_GSE']):
            geo_gse = meta_row['geo_gsm']
            if pd.isna(geo_gse) or geo_gse == 'NULL':
                df.at[idx, 'GEO_GSE'] = 'NULL'
            else:
                df.at[idx, 'GEO_GSE'] = meta_row['geo_gsm']
        
        if pd.isna(df.at[idx, 'ENCODE']):
            encode = meta_row['encode']
            if pd.isna(encode) or encode == 'NULL':
                df.at[idx, 'ENCODE'] = 'NULL'
            else:
                df.at[idx, 'ENCODE'] = meta_row['encode']
        
        if pd.isna(df.at[idx, 'EXP']):
            df.at[idx, 'EXP'] = meta_row['id']

        if pd.isna(df.at[idx, 'IS_CONTROL']):
            control_id = meta_row['control_id']
            if pd.isna(control_id) or control_id == 'NULL':
                df.at[idx, 'IS_CONTROL'] = 'False'
            else:
                df.at[idx, 'IS_CONTROL'] = meta_row['control_id']

df = df.sort_values(by = ['TF_UNIPROT_ID', 'GTRD_CELL_TYPE_ID', 'EXP'])
df.to_csv('/home/subpolare/gtrd/processed.tsv', sep = '\t', index = False)
