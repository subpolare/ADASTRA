import sys

map_dir    = sys.argv[1]
index      = sys.argv[2]
inp        = open(f'/mnt/flash/v.nachatoy/gtrd/batch{index}/tmp/filelist_for_s2.txt', 'r')
meta_table = open('/home/v.nachatoy/gtrd/meta/meta_6_may.tsv', 'r')
meta       = meta_table.readlines()
sample     = open(f'/mnt/flash/v.nachatoy/gtrd/batch{index}/filter_variants/sample_file.txt', 'w')
sample.write('sample\tindiv_id\tbam_file\tfiltered_sites_file\tsample_name\n')

for file in inp:
    cell_type = ''
    id = file.split('/')[-1].split('.')[0]
    bam_file            = f'/mnt/flash/v.nachatoy/gtrd/batch{index}/{map_dir}/' + id + '.marked.filtered.sorted.bam'
    filtered_sites_file = f'/mnt/flash/v.nachatoy/gtrd/batch{index}/filter_variants/output/{id}.bed.gz'

    for row in meta: 
        if id in row: 
            cell_type = row.split('\t')[5]
            break

    if cell_type == '': 
        break

    sample.write(f'{id}_{cell_type}\t{id}\t{bam_file}\t{filtered_sites_file}\t{id}\n')

inp.close()
sample.close()
meta_table.close()
