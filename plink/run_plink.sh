#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --time=24:00:00
#SBATCH --mem=30g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=haasx092@umn.edu
#SBATCH -p amdsmall
#SBATCH --account=jkimball
#SBATCH -o run_plink.out
#SBATCH -e run_plink.err

cd /scratch.global/haasx092/reneth_gwas/211227_snp_calling_results

module load plink

plink --vcf merged_vcf_files.vcf --mind 0.99 --double-id --allow-extra-chr --recode --out reneth_gwas

# PCA calculation
plink --pca --file reneth_gwas --allow-extra-chr -out reneth_gwas_pca
