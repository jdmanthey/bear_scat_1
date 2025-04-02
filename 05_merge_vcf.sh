#!/bin/sh
#SBATCH --chdir=./
#SBATCH --job-name=merge
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --partition=nocona
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=8G
#SBATCH --array=1-97

source activate bcftools

# define main working directory
workdir=/lustre/scratch/jmanthey/19_bears

# define variables
region_array=$( head -n${SLURM_ARRAY_TASK_ID} ${workdir}/scaffolds.txt | tail -n1 )

# run bcftools to merge the vcf files
bcftools merge -m id --regions ${region_array} ${workdir}/03_vcf/*vcf.gz > \
${workdir}/04_vcf/${region_array}.vcf

# filter for structure (only subset with good genotyping coverage)
vcftools --vcf ${workdir}/04_vcf/${region_array}.vcf --keep bear_keep.txt --max-missing 1.0 \
--mac 2 --max-alleles 2 --max-maf 0.49 --thin 1000 --recode --recode-INFO-all \
--out ${workdir}/05_pca/structure_${region_array}

# filter for kinship (all individuals)
vcftools --vcf ${workdir}/04_vcf/${region_array}.vcf --max-missing 0.2 \
--mac 1 --max-alleles 2 --max-maf 0.49 --thin 1000 --recode --recode-INFO-all \
--out ${workdir}/07_kinship/kinship_${region_array}
