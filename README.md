<div align="center">
  <h1>ADASTRA Remapping Pipeline</h1>
</div>

This repository contains a reproducible **Nextflow** pipeline that rebuilds the early stages of the **[ADASTRA workflow](https://www.nature.com/articles/s41467-021-23007-0)** — from mapped **ChIP-seq** experiments in **[GTRD database](https://gtrd.biouml.org/)** to the generation of **VCF files with allele-specific read counts**.

It implements all preprocessing steps required before downstream statistical analysis of allele-specific binding (ASB).

## 🔮 Overview

The pipeline performs alignment and SNP calling, corrects mapping bias with an updated NumPy-compatible version of **[WASP](https://github.com/bmvdgeijn/WASP)**, merges and filters BAM files, counts allele-specific reads, and produces a unified `allele_counts.vcf.gz` suitable for further ASB analysis.

### ⚙️ Main steps

1. **Input preparation.** Reads metadata from a sample table (`samples.tsv`) and genotype file (`.vcf.gz + .tbi`).

2. **HDF5 SNP tables.** Splits genotype VCF by chromosome and builds `snp_tab.h5`, `snp_index.h5`, and `haplotypes.h5` for WASP.

3. **Remapping (WASP).** Sorts and filters BAM files, rewrites reads overlapping SNPs, remaps them with `bwa mem`, and removes discordant alignments.

4. **Merging and indexing.** Combines pre- and post-remapping reads into final BAMs: `sample.bam`, `sample.passing.bam`, and their indexes.

5. **Counting allele-specific reads.** Produces compressed BED tables (`sample.bed.gz`) with read counts per heterozygous SNP.

6. **VCF recoding.** Assembles all results into a single `allele_counts.vcf.gz` with `sample_map.tsv`.

### ✨ Tech stack with sleek badge showcase

_**TL;DR: Nextflow**, Shell, Python_ 

![nextflow](https://img.shields.io/badge/nextflow-%2344CC11.svg?&style=for-the-badge&logo=nextflow&logoColor=white) ![shell](https://img.shields.io/badge/shell-%234EAA25.svg?&style=for-the-badge&logo=gnu-bash&logoColor=white) ![python](https://img.shields.io/badge/python%20-%234584B6.svg?&style=for-the-badge&logo=python&logoColor=white) ![numpy](https://img.shields.io/badge/numpy-%23013243.svg?&style=for-the-badge&logo=numpy&logoColor=white) 

- **Nextflow** ≥ 24  
- **bwa**, **samtools**, **bcftools**, **tabix**  
- **Python** ≥ 3.10 with **numpy** and **pytables**  
- **WASP** (fork [bmvdgeijn/WASP](https://github.com/bmvdgeijn/WASP), patched for modern NumPy)
