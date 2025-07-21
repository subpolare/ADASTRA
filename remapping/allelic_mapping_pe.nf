#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

params.outdir=''
params.scratch=''
params.samples_file=''
params.genotype_file=''
params.genome='/mnt/data/v.nachatoy/genome/main_nextflow/GRCh38.primary_assembly.genome'
params.indexes='/mnt/data/v.nachatoy/genome/GRCh38.primary_assembly.genome.fa'

nuclear_chroms = "$params.genome" + ".nuclear.txt"
genome_chrom_sizes_file= "$params.genome"  + ".chrom_sizes"
genome_fai = "$params.genome" + ".fa.fai"

process generate_h5_tables {

	scratch true

	input:
		path vcf_file
                path chrom_sizes
		path csi_file

	output:
		path '*.h5'

	script:
	"""
	chroms=("\$(tabix -l ${vcf_file})")
	for chrom in \${chroms[@]}; do
		echo \${chrom}
		bcftools view -r \${chrom} -Oz ${vcf_file} > \${chrom}.chrom.vcf.gz
		bcftools index \${chrom}.chrom.vcf.gz
	done

	gzip -c ${chrom_sizes} > chrom_sizes.txt.gz

	~/bin/WASP/snp2h5/snp2h5 --chrom chrom_sizes.txt.gz \
		--format vcf \
		--haplotype haplotypes.h5 \
		--snp_index snp_index.h5 \
		--snp_tab snp_tab.h5 \
		*.chrom.vcf.gz
	"""
}

process remap_bamfiles {

	tag "${sample}"
	scratch "${params.scratch}"
	publishDir params.outdir + "/remapped", mode: 'symlink'
	cpus 50

	input:
	tuple val(indiv_id), val(sample), val(bam_file), val(filtered_sites_file), val(sample_name)
	path '*'
	path '*'

	output:
	tuple val(indiv_id), val(sample), val(filtered_sites_file), path("${sample}.bam"), path("${sample}.bam.bai"), path("${sample}.passing.bam"), path("${sample}.passing.bam.bai")

	script:
	"""
	merge_files=""
	remapped_merge_files=""

	samtools sort \
		-@ ${task.cpus} \
		-o sorted.bam \
		-O bam \
		${bam_file}

	python3 ~/bin/WASP/mapping/find_intersecting_snps.py \
		--is_paired_end \
		--is_sorted \
		--output_dir \${PWD} \
		--snp_tab snp_tab.h5 \
		--snp_index snp_index.h5  \
		--haplotype haplotypes.h5 \
		--samples ${indiv_id} \
		sorted.bam

	bwa mem -t ${task.cpus} ${params.indexes} sorted.remap.fq1.gz sorted.remap.fq2.gz | samtools view -b > remapped.bam
	python3 ~/bin/bwa_script/filter_reads.py remapped.bam remapped.marked.bam ${nuclear_chroms}

	samtools sort -@ ${task.cpus} -l0 remapped.marked.bam | samtools view -b -F 512 - > remapped.marked.filtered.bam
	python3 ~/bin/WASP/mapping/filter_remapped_reads.py sorted.to.remap.bam remapped.marked.filtered.bam remapped.keep.bam
	
	merge_files="\${merge_files} sorted.bam"
	remapped_merge_files="\${remapped_merge_files} remapped.keep.bam sorted.keep.bam"

	samtools merge -f \
		reads.before.bam \
		\${merge_files}
	samtools sort \
		-@${task.cpus} \
		-o reads.before.sorted.bam  \
		reads.before.bam
	samtools merge -f \
		reads.passing.bam \
		\${remapped_merge_files}
	samtools sort \
		-@${task.cpus} \
		-o reads.passing.sorted.bam  \
		reads.passing.bam 
	mv reads.before.sorted.bam ${sample}.bam
	samtools index ${sample}.bam
	mv reads.passing.sorted.bam ${sample}.passing.bam
	samtools index ${sample}.passing.bam
	"""
}

process count_reads {

	tag "${sample}"
	publishDir params.outdir + "/count_reads", mode: 'symlink'

	input:
	tuple val(indiv_id), val(sample), val(filtered_sites_file), path(bam_all_file), path(bam_all_index_file), path(bam_passing_file), path(bam_passing_index_file)

	output:
	tuple val(indiv_id), val(sample), path("${indiv_id}.bed.gz")
	tuple path("${indiv_id}.bed.gz"), path("${indiv_id}.bed.gz.tbi")

	script:
	"""
	~/gtrd/scripts/remapping/count_tags_pileup.py ${filtered_sites_file} ${bam_all_file} ${bam_passing_file} \
	| sort-bed - | bgzip -c > ${indiv_id}.bed.gz

	tabix -p bed ${indiv_id}.bed.gz
	"""
}

process recode_vcf {
	publishDir params.outdir + "/vcf", mode: 'move'

	input:
	path samples_file
	path 'sample_map.tsv'
        path vcf_file
	path '*'
	path '*'

	output:
	path 'allele_counts.vcf.gz*'

	script:
	"""
	~/gtrd/scripts/remapping/recode_vcf.py ${vcf_file} sample_map.tsv allele_counts.vcf
	bgzip -c allele_counts.vcf > allele_counts.vcf.gz
	bcftools index allele_counts.vcf.gz
	"""
}

workflow {
    samples = Channel.fromPath(params.samples_file).splitCsv(header:true, sep:'\t').map{ row -> tuple( row.indiv_id, row.sample, row.bam_file, row.filtered_sites_file, row.sample_name) }
    genotype_h5 = generate_h5_tables(params.genotype_file,genome_chrom_sizes_file,"${params.genotype_file}.csi")
    remapped = remap_bamfiles(samples,nuclear_chroms,genotype_h5.collect())
    (counts,reads) = count_reads(remapped)
    map_file = counts.map{ [it[1], it[0], it[2].name].join("\t") }.collectFile(name: 'sample_map.tsv', newLine: true).first()
    read_file = reads.flatMap{ [file(it[0]), file(it[1])] }
    recode_vcf(params.samples_file,map_file,"${params.genotype_file}",read_file.collect(),"${params.genotype_file}.csi")
}
