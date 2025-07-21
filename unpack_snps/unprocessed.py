#!/usr/bin/env python3

# python3 scripts/unpack_snps/unprocessed.py --meta meta/meta_6_may.tsv --directory allele_counts --output unprocessed.txt

import argparse
import csv
import os

def extract_algn_id(filename):
    parts = filename.split('_')
    if len(parts) >= 3:
        return parts[2].split('.')[0]
    return None

def main():
    parser = argparse.ArgumentParser(description = 'Extract missing algn_ids from meta table that are not present in the filenames in the directory.')
    parser.add_argument('-m', '--meta', required = True, help = 'Path to the meta table file (TSV format).')
    parser.add_argument('-d', '--directory', required = True, help = 'Path to the directory containing the files.')
    parser.add_argument('-o', '--output', required = True, help = 'Output file name to save the missing algn_ids.')
    args = parser.parse_args()
    
    meta_algn_ids = set()
    with open(args.meta, 'r', newline = '') as meta_file:
        reader = csv.DictReader(meta_file, delimiter = '\t')
        for row in reader:
            algn_id = row.get('algn_id')
            if algn_id:
                meta_algn_ids.add(algn_id.strip())
    
    file_algn_ids = set()
    for fname in os.listdir(args.directory):
        full_path = os.path.join(args.directory, fname)
        if os.path.isfile(full_path):
            algn_id = extract_algn_id(fname)
            if algn_id:
                file_algn_ids.add(algn_id.strip())
    
    missing_ids = sorted(meta_algn_ids - file_algn_ids)
    
    with open(args.output, 'w') as out_file:
        for algn_id in missing_ids:
            out_file.write(algn_id + '\n')

if __name__ == '__main__':
    main()
