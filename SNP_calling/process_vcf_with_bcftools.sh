# WRITE ONLY FULL ADDREESSES

index=$2
genome=/mnt/data/v.nachatoy/genome/GRCh38.primary_assembly.genome.fa
vcfs=/mnt/flash/v.nachatoy/gtrd/batch${index}/VCFs/
temp=/mnt/flash/v.nachatoy/gtrd/batch${index}/tmp/

# DO NOT EDIT BELOW

bcftools mpileup -A -r $1 --fasta-ref $genome --redo-BAQ --adjust-MQ 50 --gap-frac 0.05 --max-depth 10000 --max-idepth 200000 --annotate FORMAT/DP,FORMAT/AD --bam-list ${temp}filelist_for_s2.txt --output-type u \
| bcftools call --keep-alts --multiallelic-caller --format-fields GQ --output-type v \
| bcftools filter -i "INFO/DP>=10" --output-type z - \
| bcftools norm --check-ref x -m - --fasta-ref $genome \
| bcftools filter -i "QUAL>=10 & FORMAT/GQ>=20 & FORMAT/DP>=10" --SnpGap 3 --IndelGap 10 --set-GTs . \
| bcftools view -i 'GT="alt"' --trim-alt-alleles \
| bcftools annotate -x ^INFO/DP \
| bcftools +fill-tags -- -t all > ${vcfs}${1}.filtered.vcf

bcftools view ${vcfs}${1}.filtered.vcf -Oz -o ${vcfs}${1}.filtered.vcf.gz
bcftools index ${vcfs}${1}.filtered.vcf.gz

bcftools annotate -r $1 -a /sandbox/buyanchik/genome/common_all_20180418_compressed.vcf.gz --columns ID,CAF,TOPMED --output-type z ${vcfs}${1}.filtered.vcf.gz -o ${vcfs}${1}.filtered.annotated.vcf
bcftools view ${vcfs}${1}.filtered.annotated.vcf -Oz -o ${vcfs}${1}.filtered.annotated.vcf.gz
bcftools index ${vcfs}${1}.filtered.annotated.vcf.gz
