#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --time=24:00:00
#SBATCH --mem=30g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=haasx092@umn.edu
#SBATCH -p amdsmall
#SBATCH --account=jkimball
#SBATCH -o run_normalize_awk.out
#SBATCH -e run_normalize_awk.err

cd /scratch.global/haasx092/reneth_gwas/211227_snp_calling_results

prefix="211227_snp_calling_results"
normalize_prefix="211227_snp_calling"
mktemp | read tmp
cat ${prefix}_*_filtered.recode.vcf | awk -f ./normalize.awk > '$tmp' 2> $normalize_prefix.err
cp '$tmp' 211227_snp_calling_normalize.tsv
