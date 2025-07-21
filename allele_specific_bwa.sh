index=$1
threads=40
main=/home/v.nachatoy/gtrd/
temp=/mnt/flash/v.nachatoy/gtrd/batch${index}/tmp/
home=/mnt/flash/v.nachatoy/gtrd/batch${index}/
genome_dir=/mnt/data/v.nachatoy/genome/
reads_dir=/mnt/data/v.nachatoy/gtrd/reads/
scripts=/home/v.nachatoy/gtrd/scripts/
filter_reads=/home/v.nachatoy/bin/bwa_script/filter_reads.py
genome=/mnt/data/v.nachatoy/genome/GRCh38.primary_assembly.genome.fa

filter_se_reads() {
    index=$2
    home=/mnt/flash/v.nachatoy/gtrd/batch${index}/
    genome_dir=/mnt/data/v.nachatoy/genome/

    name=$(echo $1 | cut -d '.' -f1)
    map_dir=${home}mapped_se/
    python3 /home/v.nachatoy/bin/bwa_script/filter_reads.py ${map_dir}${name}.bam ${map_dir}${name}.marked.bam ${genome_dir}GRCh38.primary_assembly.genome.nuclear.txt
    rm ${map_dir}${name}.bam
}

filter_pe_reads() {
    index=$2
    home=/mnt/flash/v.nachatoy/gtrd/batch${index}/
    genome_dir=/mnt/data/v.nachatoy/genome/

    name=$(echo $1 | cut -d '.' -f1)
    map_dir=${home}mapped_pe/
    python3 /home/v.nachatoy/bin/bwa_script/filter_reads.py ${map_dir}${name}.bam ${map_dir}${name}.marked.bam ${genome_dir}GRCh38.primary_assembly.genome.nuclear.txt
    rm ${map_dir}${name}.bam
}

post_processing() {
    temp=/mnt/flash/v.nachatoy/gtrd/batch${index}/tmp/
    threads=40
    map_dir=$1
    index=$2

    for file in $(ls -1 ${map_dir}*.marked.bam); do
            name=$(echo $file | rev | cut -d'/' -f1 | rev | cut -d'.' -f1)
        samtools sort -T ${temp} -@ $threads -l0 ${map_dir}${name}.marked.bam | \
        samtools view -@ $threads -b -F 512 - > ${map_dir}${name}.marked.filtered.bam
        samtools sort -T ${temp} -@ $threads -l0 ${map_dir}${name}.marked.filtered.bam > ${map_dir}${name}.marked.filtered.sorted.bam
        samtools index ${map_dir}${name}.marked.filtered.sorted.bam
        rm ${map_dir}${name}.marked.bam ${map_dir}${name}.marked.filtered.bam
    done
}

SNP_calling() {
    index=$3
    map_dir=$1
    threads=40
    scripts=/home/v.nachatoy/gtrd/scripts/
    genome_dir=/mnt/data/v.nachatoy/genome/
    temp=/mnt/flash/v.nachatoy/gtrd/batch${index}/tmp/
    vcfs=/mnt/flash/v.nachatoy/gtrd/batch${index}/VCFs/

    ls -1 ${map_dir}*.sorted.bam > ${temp}batch${index}_filelist_for_s2.txt
    parallel -j $threads -a ${genome_dir}hg38.to_parallel_bcftools ${scripts}SNP_calling/process_vcf_with_bcftools.sh {} $index
    sed "s/^/\/mnt\/flash\/v.nachatoy\/gtrd\/batch${index}\/VCFs\//g" ${genome_dir}hg38.to_parallel_bcftools | sed 's/$/\.filtered\.annotated\.vcf\.gz/g' > ${temp}batch${index}.sorted.txt

    bcftools concat --threads $threads --output-type z -f ${temp}batch${index}.sorted.txt > ${vcfs}all.filtered.annotated.vcf.gz
    bcftools index --threads $threads ${vcfs}all.filtered.annotated.vcf.gz
    bcftools view --threads $threads -m2 -M2 -v snps --output-type z ${vcfs}all.filtered.annotated.vcf.gz > ${vcfs}all.filtered.annotated.snps.vcf.gz
    bcftools index --threads $threads ${vcfs}all.filtered.annotated.snps.vcf.gz
}

export -f filter_se_reads
export -f filter_pe_reads
export -f post_processing
export -f SNP_calling

for reads in $(ls -1 ${main}batch${index}_reads?.txt); do
    while IFS= read -r filename; do
        file=${reads_dir}${filename}
        name=$(echo $file | rev | cut -d'/' -f1 | rev | cut -d'.' -f1)
        samtools sort -T ${temp} -@ $threads -n $file -o ${temp}${name}.sorted.bam
        samtools bam2fq -@ $threads ${temp}${name}.sorted.bam \
            -0 ${temp}${name}_unpaired.fastq.gz \
            -1 ${temp}${name}_R1.fastq.gz \
            -2 ${temp}${name}_R2.fastq.gz

        if [ $(stat -c%s ${temp}${name}_unpaired.fastq.gz) -lt 1048576 ]; then
            map_dir=${home}mapped_pe/
            bwa mem -t $threads -R "@RG\tID:${name}\tSM:${name}" $genome \
                ${temp}${name}_R1.fastq.gz ${temp}${name}_R2.fastq.gz \
                | samtools view -@ $threads -b > ${map_dir}${name}.bam
        else
            map_dir=${home}mapped_se/
            bwa mem -t $threads -R "@RG\tID:${name}\tSM:${name}" $genome \
                ${temp}${name}_unpaired.fastq.gz \
                | samtools view -@ $threads -b > ${map_dir}${name}.bam
        fi

        rm -r ${temp}${name}.sorted.bam ${temp}${name}_R1.fastq.gz ${temp}${name}_R2.fastq.gz ${temp}${name}_unpaired.fastq.gz
    done < $reads

    if [ "$(find ${home}mapped_se/ -type f | wc -l)" -gt 0 ]; then
        parallel -j $threads filter_se_reads ::: $(ls -1 ${home}mapped_se/) ::: $index
        post_processing ${home}mapped_se/ $index
        SNP_calling ${home}mapped_se/ ${home}VCFs/ $index
        python3 ${scripts}filter_variants/sample_file.py mapped_se $index
        nextflow run -dsl2 ${scripts}filter_variants/filter_variants.nf \
            -w ${home}work --outdir ${home}filter_variants/output \
            --samples_file ${home}filter_variants/sample_file.txt \
            --genotype_file ${home}VCFs/all.filtered.annotated.snps.vcf.gz
        nextflow run ${scripts}remapping/allelic_mapping_se.nf \
            -w ${home}work --outdir ${home}remapping \
            --samples_file ${home}filter_variants/sample_file.txt \
            --genotype_file ${home}VCFs/all.filtered.annotated.snps.vcf.gz \
            --scratch $temp
        ${scripts}remapping/unpack_snps.sh $index
        python3 ${scripts}unpack_snps/factor_changer.py
    fi

    if [ "$(find ${home}mapped_pe/ -type f | wc -l)" -gt 0 ]; then
        parallel -j $threads filter_pe_reads ::: $(ls -1 ${home}mapped_pe/) ::: $index
        post_processing ${home}mapped_pe/ $index
        SNP_calling ${home}mapped_pe/ ${home}VCFs/ $index
        python3 ${scripts}filter_variants/sample_file.py mapped_pe $index
        nextflow run -dsl2 ${scripts}filter_variants/filter_variants.nf \
            -w ${home}work --outdir ${home}filter_variants/output \
            --samples_file ${home}filter_variants/sample_file.txt \
            --genotype_file ${home}VCFs/all.filtered.annotated.snps.vcf.gz
        nextflow run ${scripts}remapping/allelic_mapping_pe.nf \
            -w ${home}work --outdir ${home}remapping \
            --samples_file ${home}filter_variants/sample_file.txt \
            --genotype_file ${home}VCFs/all.filtered.annotated.snps.vcf.gz \
            --scratch $temp
        ${scripts}remapping/unpack_snps.sh $index
        python3 ${scripts}unpack_snps/factor_changer.py
    fi

    rm -r ${home}mapped_se/* ${home}mapped_pe/* ${temp}* ${home}VCFs/* ${home}filter_variants/output/* ${home}remapping/* ${home}work/*
done

rm ${main}reads?.txt
head -n 50 ${main}batch${index}_download.txt > ${main}batch${index}_reads1.txt
head -n 100 ${main}batch${index}_download.txt | tail -n 50 > ${main}batch${index}_reads2.txt
head -n 150 ${main}batch${index}_download.txt | tail -n 50 > ${main}batch${index}_reads3.txt
head -n 200 ${main}batch${index}_download.txt | tail -n 50 > ${main}batch${index}_reads4.txt
head -n 250 ${main}batch${index}_download.txt | tail -n 50 > ${main}batch${index}_reads5.txt