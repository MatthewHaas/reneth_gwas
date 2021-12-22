#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --time=24:00:00
#SBATCH --mem=30g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=haasx092@umn.edu
#SBATCH -p amdsmall
#SBATCH --account=jkimball
#SBATCH -o run_cutadapt.out
#SBATCH -e run_cutadapt.err

cd /scratch.global/haasx092/reneth_gwas

module load cutadapt
module load parallel

# Note: the trimmed fastq files are going to the scratch directory where they will be deleted after 30 days
# Adapter sequence corresponds to the Nextera Transposase Sequence
cat 211222_reneth_gwas_sample_directory_list.txt | parallel -k --will-cite -j 10 \
        sh ./cutadapt_wrapper_script.sh
