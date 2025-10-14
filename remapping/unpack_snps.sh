#!/usr/bin/bash

allele_counts=/mnt/flash/v.nachatoy/gtrd/batch${index}/remapping/vcf/allele_counts.vcf.gz
results=/mnt/data/v.nachatoy/gtrd/allele_counts

for sample in `bcftools query -l $allele_counts`; do
    bcftools view -c1 --output-type v -s ${sample} ${allele_counts} | bcftools filter -i 'AD[0:1] >= 5 && AD[0:0] >= 5' --output-type z -o ${results}/${sample}.vcf.gz
done
