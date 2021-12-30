#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --time=24:00:00
#SBATCH --mem=30g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=haasx092@umn.edu
#SBATCH -p amdsmall
#SBATCH --account=jkimball
#SBATCH -o concat_filtered_vcfs.out
#SBATCH -e concat_filtered_vcfs.err

# Include path to the working directory
cd /scratch.global/haasx092/reneth_gwas/211227_snp_calling_results

module load bcftools

chr1='211227_snp_calling_results_ZPchr0001_filtered.recode.vcf'
chr2='211227_snp_calling_results_ZPchr0002_filtered.recode.vcf'
chr3='211227_snp_calling_results_ZPchr0003_filtered.recode.vcf'
chr4='211227_snp_calling_results_ZPchr0004_filtered.recode.vcf'
chr5='211227_snp_calling_results_ZPchr0005_filtered.recode.vcf'
chr6='211227_snp_calling_results_ZPchr0006_filtered.recode.vcf'
chr7='211227_snp_calling_results_ZPchr0007_filtered.recode.vcf'
chr8='211227_snp_calling_results_ZPchr0008_filtered.recode.vcf'
chr9='211227_snp_calling_results_ZPchr0009_filtered.recode.vcf'
chr10='211227_snp_calling_results_ZPchr0010_filtered.recode.vcf'
chr11='211227_snp_calling_results_ZPchr0011_filtered.recode.vcf'
chr12='211227_snp_calling_results_ZPchr0012_filtered.recode.vcf'
chr13='211227_snp_calling_results_ZPchr0013_filtered.recode.vcf'
chr14='211227_snp_calling_results_ZPchr0014_filtered.recode.vcf'
chr15='211227_snp_calling_results_ZPchr0015_filtered.recode.vcf'
scf16='211227_snp_calling_results_ZPchr0016_filtered.recode.vcf'
scf458='211227_snp_calling_results_ZPchr0458_filtered.recode.vcf'

bcftools concat $chr1 $chr2 $chr3 $chr4 $chr5 $chr6 $chr7 $chr8 $chr9 $chr10 $chr11 $chr12 $chr13 $chr14 $chr15 $scf16 $scf458 > merged_vcf_files.vcf
