for file in $(ls -1 /sandbox/subpolare/gtrd/allele_counts); do
    name=$(echo $file | cut -d'_' -f3 | cut -d'.' -f1)
    cell=$(echo $file | cut -d'_' -f2)
    unip=$(echo $file | cut -d'_' -f1)'_HUMAN'
    echo -e "NA\t$name\tNA\t$unip\t$cell\tNA\tNA\tNA\tNA\tNA\tNA" >> ~/gtrd/processed.tsv
done

for file in /sandbox/subpolare/gtrd/allele_counts/*.vcf.gz; do
    python3 /home/subpolare/gtrd/scripts/unpack_snps/uniprot_parser.py \
        --input $file \
        --meta /home/subpolare/gtrd/meta/meta_6_may.tsv
done

mv /sandbox/subpolare/gtrd/allele_counts/* /home/subpolare/gtrd/allele_counts/
