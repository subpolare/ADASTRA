import os

allele_counts_dir = '/home/subpolare/gtrd/allele_counts'

allele_files = os.listdir(allele_counts_dir)
allele_ids = set(f.split('_')[0] for f in allele_files if f.startswith('ALIGNS'))
allele_ids.update(f.split('_')[1] for f in allele_files if not f.startswith('ALIGNS'))

with open('/home/subpolare/gtrd/factors_gtrd.txt', 'r') as f:
    lines = f.readlines()

filtered_lines = [line for line in lines if not any(align_id in line for align_id in allele_ids)]

with open('/home/subpolare/gtrd/factors_gtrd_filtered.txt', 'w') as f:
    f.writelines(filtered_lines)

