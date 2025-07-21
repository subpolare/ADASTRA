import os

vcf_dir = os.path.expanduser('/home/subpolare/gtrd/allele_counts/')
vcf_files = [f for f in os.listdir(vcf_dir) if f.endswith('.vcf.gz')]
vcf_ids = [f.split('_')[1] for f in vcf_files]

with open('/home/subpolare/gtrd/CTCF.txt') as file:
    for line in file:
        align_id = line.strip().split('/')[-1].split('.')[0]
        if align_id not in vcf_ids:
            print(f'{align_id}')
