#!/bin/bash
#SBATCH --chdir=./
#SBATCH --job-name=heterozygosity
#SBATCH --partition nocona
#SBATCH --nodes=1 --ntasks=2
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=4G
#SBATCH --array=1-15

workdir=/lustre/scratch/jmanthey/19_bears

basename_array=$( head -n${SLURM_ARRAY_TASK_ID} ${workdir}/basenames_bears.txt | tail -n1 )

gzip -cd ${workdir}/03_vcf/${basename_array}.vcf.gz | grep -v "^#" | cut -f10 > ${basename_array}.genotypes

echo ${basename_array} > ${basename_array}.diversity

echo "0/0 count" >> ${basename_array}.diversity
grep "0/0" ${basename_array}.genotypes | wc -l >> ${basename_array}.diversity

echo "0/1 count" >> ${basename_array}.diversity
grep "0/1" ${basename_array}.genotypes | wc -l >> ${basename_array}.diversity

echo "1/1 count" >> ${basename_array}.diversity
grep "1/1" ${basename_array}.genotypes | wc -l >> ${basename_array}.diversity

echo "total count" >> ${basename_array}.diversity
grep -v "\\./\\." ${basename_array}.genotypes | wc -l >> ${basename_array}.diversity

rm ${basename_array}.genotypes


