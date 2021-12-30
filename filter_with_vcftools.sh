#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --time=24:00:00
#SBATCH --mem=30g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=haasx092@umn.edu
#SBATCH -p amdsmall
#SBATCH --account=jkimball
#SBATCH -o filter_with_vcftools.out
#SBATCH -e filter_with_vcftools.err

# Include path to the working directory
cd /scratch.global/haasx092/reneth_gwas/211227_snp_calling_results

# Filter VCF files for missing data (20% missing/80% present), biallelic sites only, no indels, and a depth of 4
for i in $(cat vcf_file_list.txt)
do
STEM=$(echo ${i} | cut -f 1 -d ".")
~/vcftools/bin/vcftools --gzvcf  $i --max-missing 0.90 --min-alleles 2 --max-alleles 2 --maf 0.03 --remove-indels --minDP 8 --recode --recode-INFO-all --out ${STEM}_filtered
done
