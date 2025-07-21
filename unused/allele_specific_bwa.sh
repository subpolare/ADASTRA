genome=/home/subpolare/genome/GRCh38.primary_assembly.genome.fa
reads=$1
bwa=/home/subpolare/bin/bwa/bwa
threads=50

# МНОГО ФУНКЦИЙ
# Фильтрация одноконцевых чтений после BWA
filter_se_reads() {
    name=$(echo $1 | cut -d '.' -f1)
    map_dir='/sandbox/subpolare/gtrd/mapped_se/'
    python3 /home/subpolare/bin/bwa_script/filter_reads.py ${map_dir}${name}.bam /sandbox/subpolare/tmp/${name}.marked.bam /home/subpolare/genome/GRCh38.primary_assembly.genome.nuclear.txt
}

# Фильтрация парноконцевых чтений после BWA 
filter_pe_reads() {
    name=$(echo $1 | cut -d '.' -f1)
    map_dir='/sandbox/subpolare/gtrd/mapped_pe/'
    python3 /home/subpolare/bin/bwa_script/filter_reads.py ${map_dir}${name}.bam /sandbox/subpolare/tmp/${name}.marked.bam /home/subpolare/genome/GRCh38.primary_assembly.genome.nuclear.txt
}

# Фильтрация и обработка BAM после фильтрации BWA скриптом. Единственное отличие от Андрея пометил комментарием
post_processing() {
    map_dir=$1
    files=(/sandbox/subpolare/tmp/*.marked.bam)
    for file in "${files[@]}"; do
	name=$(echo $file | rev | cut -d'/' -f1 | rev | cut -d'.' -f1)
        samtools sort -@ $threads -l0 /sandbox/subpolare/tmp/${name}.marked.bam | samtools view -b -F 512 - > ${map_dir}${name}.marked.filtered.bam
        samtools sort -@ $threads -l0 ${map_dir}${name}.marked.filtered.bam > ${map_dir}${name}.marked.filtered.sorted.bam
        # Единственный скрипт, которого не было в пайплайне Андрея. Почему-то после BWA исчезает тэг RG внутри BAM, этот скрипт возвращает его на место
        /home/subpolare/gtrd/scripts/BAM_processing/add_RG.sh ${map_dir}${name}.marked.filtered.sorted.bam  # на выходе .marked.filtered.sorted.edited.bam
        samtools sort -@ $threads -l0 ${map_dir}${name}.marked.filtered.sorted.edited.bam > ${map_dir}${name}.marked.filtered.edited.sorted.bam
        samtools index ${map_dir}${name}.marked.filtered.edited.sorted.bam
        rm ${map_dir}${name}.bam ${map_dir}${name}.marked.filtered.bam ${map_dir}${name}.marked.filtered.sorted.edited.bam ${map_dir}${name}.marked.filtered.sorted.bam
    	rm /sandbox/subpolare/tmp/${name}.marked.bam 
    done  
}

# Снипокллинг. Аналогичен тому, что было у Андрея, за небольшими исключениями, пометил их
SNP_calling() {
    vcfs=$2
	ls -1 $1*.sorted.bam > /sandbox/subpolare/tmp/filelist_for_s2.txt
	awk '{print $1}' /sandbox/buyanchik/genome/GRCh38.primary_assembly.genome.chrom_sizes > /home/subpolare/genome/hg38.to_parallel_bcftools
    # Разница между двумя скриптами process_vcf_with_bcftools_[se/pe].sh только в том, в какую директорию они складываюи результаты (VCFs_se или VCFs_pe)
	parallel -j $threads -a /home/subpolare/genome/hg38.to_parallel_bcftools /home/subpolare/gtrd/scripts/SNP_calling/process_vcf_with_bcftools_${3}.sh
	
    # Почему-то впихивание внутрь выражения sed переменной $3 не привело ни к чему, кроме еррора, так что сделал такой костыль
	if [ "$3" = "se" ]; then
		sed 's/^/\/sandbox\/subpolare\/gtrd\/VCFs_se\//g' /home/subpolare/genome/hg38.to_parallel_bcftools | sed 's/$/\.filtered\.annotated\.vcf\.gz/g' > /sandbox/subpolare/tmp/files.sorted.txt
	elif [ "$3" = "pe" ]; then
		sed 's/^/\/sandbox\/subpolare\/gtrd\/VCFs_pe\//g' /home/subpolare/genome/hg38.to_parallel_bcftools | sed 's/$/\.filtered\.annotated\.vcf\.gz/g' > /sandbox/subpolare/tmp/files.sorted.txt
	fi
	
	bcftools concat --output-type z -f /sandbox/subpolare/tmp/files.sorted.txt > ${vcfs}all.filtered.annotated.vcf.gz
	bcftools index ${vcfs}all.filtered.annotated.vcf.gz
	bcftools view -m2 -M2 -v snps --output-type z ${vcfs}all.filtered.annotated.vcf.gz > ${vcfs}all.filtered.annotated.snps.vcf.gz
	bcftools index ${vcfs}all.filtered.annotated.snps.vcf.gz
    # Складирование промежуточных VCF (до WASP) внутрь директории ~/gtrd/vcfs_before_wasp
	/home/subpolare/gtrd/scripts/SNP_calling/unpack_snps.sh $3
	rm /sandbox/subpolare/tmp/files.sorted.txt
}

export -f filter_se_reads
export -f filter_pe_reads
export -f post_processing
export -f SNP_calling

# ОСНОВНОЙ ПАЙПЛАЙН
# Ниже все BAM из ~/gtrd/gtrd_files сначала конвертируются в fastq, а потом заново картируются при помощи BWA
while IFS= read -r filename; do
    file='/home/subpolare/gtrd/reads'$filename
    name=$(echo $file | rev | cut -d'/' -f1 | rev | cut -d'.' -f1)
    samtools sort -@ $threads -n $file -o /sandbox/subpolare/tmp/${name}.sorted.bam
    samtools bam2fq -@ $threads /sandbox/subpolare/tmp/${name}.sorted.bam -0 /sandbox/subpolare/tmp/${name}_unpaired.fastq.gz -1 /sandbox/subpolare/tmp/${name}_R1.fastq.gz -2 /sandbox/subpolare/tmp/${name}_R2.fastq.gz

    # Проверяет, есть ли в непарных чтениях что-то (пара ридов там проскальзывает, поэтому поставил порог по байтам чуть больше, чем ноль)
    if [ $(stat -c%s "/sandbox/subpolare/tmp/${name}_unpaired.fastq.gz") -lt 1048576 ]; then     # Если плюс-минус ничего нет, то картируем парные чтения внутрь директории mapped_pe
        map_dir=/sandbox/subpolare/gtrd/mapped_pe/
        $bwa mem -t $threads $genome /sandbox/subpolare/tmp/${name}_R1.fastq.gz /sandbox/subpolare/tmp/${name}_R2.fastq.gz | samtools view -@ $threads -b > ${map_dir}${name}.bam
    else                                                                       # Если что-то есть, то картируем непарные чтения внутрь директории mapped_se
        map_dir=/sandbox/subpolare/gtrd/mapped_se/                                
        $bwa mem -t $threads $genome /sandbox/subpolare/tmp/${name}_unpaired.fastq.gz | samtools view -@ $threads -b > ${map_dir}${name}.bam
    fi

    rm -r /sandbox/subpolare/tmp/${name}.sorted.bam /sandbox/subpolare/tmp/${name}_R1.fastq.gz /sandbox/subpolare/tmp/${name}_R2.fastq.gz /sandbox/subpolare/tmp/${name}_unpaired.fastq.gz
done < $reads

# Обработка одноконцевых чтений
if [ "$(find "/sandbox/subpolare/gtrd/mapped_se/" -type f | wc -l)" -gt 0 ]; then            # Если в папке mapped_se есть хоть что-то, то... 
    parallel -j $threads filter_se_reads ::: $(ls -1 '/sandbox/subpolare/gtrd/mapped_se/')   # Фильтрация после выравнивания при помощи bwa_script/filter_reads.py
    post_processing '/sandbox/subpolare/gtrd/mapped_se/'                                     # Сортировка, фильтрация и тд полученных BAM файлов
    SNP_calling /sandbox/subpolare/gtrd/mapped_se/ /sandbox/subpolare/gtrd/VCFs_se/ se       # Снипколлинг
    ls -1 /sandbox/subpolare/gtrd/mapped_se/*.sorted.bam > /sandbox/subpolare/tmp/filelist_for_s2.txt                               
    python3 /home/subpolare/gtrd/scripts/filter_variants/sample_file.py mapped_se            # Создание sample_file для nextflow скриптов
    nextflow run -dsl2 /home/subpolare/gtrd/scripts/filter_variants/filter_variants.nf -w /sandbox/subpolare/gtrd/work --genotype_file '/sandbox/subpolare/gtrd/VCFs_se/all.filtered.annotated.snps.vcf.gz'
    nextflow run scripts/remapping/allelic_mapping_se.nf -w /sandbox/subpolare/gtrd/work     # Перекартирование и WASP 
    ~/gtrd/scripts/remapping/unpack_snps.sh                                                  # Распаковка отдельных VCF внутрь ~/gtrd/allele_counts/
    python3 ~/gtrd/scripts/unpack_snps/factor_changer.py
fi

# Обработка парноконцевых чтений, полностью аналогична
if [ "$(find "/sandbox/subpolare/gtrd/mapped_pe/" -type f | wc -l)" -gt 0 ]; then  
    parallel -j $threads filter_pe_reads ::: $(ls -1 '/sandbox/subpolare/gtrd/mapped_pe/')
    post_processing '/sandbox/subpolare/gtrd/mapped_pe/'
    SNP_calling /sandbox/subpolare/gtrd/mapped_pe/ /sandbox/subpolare/gtrd/VCFs_pe/ pe
    ls -1 /sandbox/subpolare/gtrd/mapped_pe/*.sorted.bam > /sandbox/subpolare/tmp/filelist_for_s2.txt
    python3 /home/subpolare/gtrd/scripts/filter_variants/sample_file.py mapped_pe
    nextflow run -dsl2 /home/subpolare/gtrd/scripts/filter_variants/filter_variants.nf -w /sandbox/subpolare/gtrd/work --genotype_file '/sandbox/subpolare/gtrd/VCFs_pe/all.filtered.annotated.snps.vcf.gz'
    nextflow run scripts/remapping/allelic_mapping_pe.nf -w /sandbox/subpolare/gtrd/work
    ~/gtrd/scripts/remapping/unpack_snps.sh
    python3 ~/gtrd/scripts/unpack_snps/factor_changer.py
fi

# Удаление всего ненужного
# rm -r /sandbox/subpolare/gtrd/mapped_se/* /sandbox/subpolare/gtrd/mapped_pe/* /sandbox/subpolare/tmp/* 
# rm -r /sandbox/subpolare/gtrd/VCFs_se/* /sandbox/subpolare/gtrd/VCFs_pe/* scripts/filter_variants/output/* scripts/remapping/output/* 
# rm -r /sandbox/subpolare/gtrd/work/* /home/subpolare/gtrd/.nextflow*
