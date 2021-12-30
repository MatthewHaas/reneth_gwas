#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --time=24:00:00
#SBATCH --mem=30g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=haasx092@umn.edu
#SBATCH -p amdsmall
#SBATCH --account=jkimball
#SBATCH -o plot_plink_pca.out
#SBATCH -e plot_plink_pca.err

cd /scratch.global/haasx092/reneth_gwas/211227_snp_calling_results

module load R/3.6.0

Rscript plot_plink_pca.R reneth_gwas_pca.eigenvec reneth_gwas_pca.eigenval 211227_reneth_gwas.pdf 211227_reneth_gwas.Rdata
