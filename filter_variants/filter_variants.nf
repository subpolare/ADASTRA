#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

params.samples_file=''
params.genotype_file=''
params.outdir=''

params.min_GQ=50 // Minimum genotype quality
params.min_DP=10 // Minimum read depth over SNP
params.min_AD=5  // Minimum reads per alleles

process filter_variants {

	tag "${indiv_id}"
	
    publishDir "${params.outdir}", mode: 'move'

	input:
	val indiv_id
	path genotype_file
	val min_DP
	val min_AD
	val min_GQ

	output:
	path "${indiv_id}.bed.gz"
        path "${indiv_id}.bed.gz.tbi"

	script:
	"""
        bcftools query \
		-s ${indiv_id} \
		-i'GT="1/1" || GT="1/0" || GT="0/1"' \
		-f'%CHROM\\t%POS0\\t%POS\\t%CHROM:%POS:%REF:%ALT\\t%REF\\t%ALT\\t[%GT\\t%GQ\\t%DP\\t%AD{0}\\t%AD{1}]\\n' \
		${genotype_file} \
	| awk -v OFS="\\t" \
		-v min_GQ=${min_GQ} -v min_AD=${min_AD} -v min_DP=${min_DP}\
		'\$8<min_GQ { next; } \$9<min_DP { next; }\
		(\$7=="0/1" || \$7=="1/0" || \$7=="0|1" || \$7=="1|0") && (\$10<min_AD || \$11<min_AD) { next; } \
		{ print; }' \
	| sort-bed - \
	| grep -v chrX | grep -v chrY | grep -v chrM | grep -v _random | grep -v _alt | grep -v chrUn | grep -v KI | grep -v GL \
	| bgzip -c > ${indiv_id}.bed.gz

	tabix -p bed ${indiv_id}.bed.gz
	"""
}

workflow {
    indiv = Channel.fromPath(params.samples_file).splitCsv(header:true, sep:'\t').map{ row -> row.indiv_id }.unique()
    filter_variants(indiv,params.genotype_file,params.min_DP,params.min_AD,params.min_GQ)
}

        
