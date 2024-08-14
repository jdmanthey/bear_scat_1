#!/bin/sh
#SBATCH --chdir=./
#SBATCH --job-name=filter
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --partition=nocona
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=8G
#SBATCH --array=1-97

source activate bcftools

# define main working directory
workdir=/lustre/scratch/jmanthey/19_bears

# define variables
region_array=$( head -n${SLURM_ARRAY_TASK_ID} ${workdir}/scaffolds.txt | tail -n1 )

# filter for structure
vcftools --vcf ${workdir}/04_vcf/${region_array}.vcf --max-missing 1.0 \
--mac 2 --max-alleles 2 --max-maf 0.49 --recode --recode-INFO-all \
--out ${workdir}/05_filtered/structure_${region_array}

# filter for kinship

vcftools --vcf ${workdir}/04_vcf/${region_array}.vcf --max-missing 0.2 \
--mac 1 --max-alleles 2 --max-maf 0.49 --recode --recode-INFO-all \
--out ${workdir}/05_filtered/kinship_${region_array}

