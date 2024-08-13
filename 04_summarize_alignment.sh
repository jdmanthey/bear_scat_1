#!/bin/bash
#SBATCH --chdir=./
#SBATCH --job-name=align_stats
#SBATCH --partition nocona
#SBATCH --nodes=1 --ntasks=2
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=4G
#SBATCH --array=1-15

source activate bcftools

workdir=/lustre/scratch/jmanthey/19_bears

basename_array=$( head -n${SLURM_ARRAY_TASK_ID} ${workdir}/basenames_bears.txt | tail -n1 )

echo ${basename_array} > ${basename_array}.stats

# samtools depth sum of aligned sites
echo "samtools depth sum of aligned sites" >> ${basename_array}.stats
samtools depth  ${workdir}/01_bam_files/${basename_array}_final.bam  |  awk '{sum+=$3} END { print "Sum = ",sum}' >> ${basename_array}.stats

# proportion dupes
echo "proportion duplicates" >> ${basename_array}.stats
head -n8 ${workdir}/01_bam_files/${basename_array}_markdups_metric_file.txt | tail -n1 | cut -f9 >> ${basename_array}.stats

# number of genotyped sites passing minimum depth filter
echo "sites genotyped" >> ${basename_array}.stats
gzip -cd ${workdir}/03_vcf/${basename_array}.vcf.gz | grep -v "^#" | wc -l >> ${basename_array}.stats
