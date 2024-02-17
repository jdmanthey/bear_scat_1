#!/bin/bash
#SBATCH --chdir=./
#SBATCH --job-name=bam2
#SBATCH --partition quanah
#SBATCH --nodes=1 --ntasks=12
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=8G
#SBATCH --array=1-26

module load intel java singularity bowtie2

export SINGULARITY_CACHEDIR="/lustre/work/jmanthey/singularity-cachedir"

# define main working directory
workdir=/lustre/scratch/jmanthey/02_bear_scat

basename_array=$( head -n${SLURM_ARRAY_TASK_ID} ${workdir}/basenames_bears.txt | tail -n1 )

# define the reference genome
refgenome=/home/jmanthey/references/GCA_024610735.1_mUrsAme1.0.p_genomic.fna

# run bowtie2
bowtie2 --threads 12 -x /home/jmanthey/references/bear -1 ${workdir}/01_cleaned/${basename_array}_R1.fastq.gz -2 ${workdir}/01_cleaned/${basename_array}_R2.fastq.gz -S ${workdir}/01b_bam_files/${basename_array}.sam

# filter for mapped reads
samtools view -b -f 2 -@ 12 -o ${workdir}/01b_bam_files/${basename_array}.bam ${workdir}/01b_bam_files/${basename_array}.sam

# remove sam
rm ${workdir}/01b_bam_files/${basename_array}.sam

# clean up the bam file
singularity exec $SINGULARITY_CACHEDIR/gatk_4.2.3.0.sif gatk CleanSam -I ${workdir}/01b_bam_files/${basename_array}.bam -O ${workdir}/01b_bam_files/${basename_array}_cleaned.bam

# remove the raw bam
rm ${workdir}/01b_bam_files/${basename_array}.bam

# sort the cleaned bam file
singularity exec $SINGULARITY_CACHEDIR/gatk_4.2.3.0.sif gatk SortSam -I ${workdir}/01b_bam_files/${basename_array}_cleaned.bam -O ${workdir}/01b_bam_files/${basename_array}_cleaned_sorted.bam --SORT_ORDER coordinate

# remove the cleaned bam file
rm ${workdir}/01b_bam_files/${basename_array}_cleaned.bam

# add read groups to sorted and cleaned bam file
singularity exec $SINGULARITY_CACHEDIR/gatk_4.2.3.0.sif gatk AddOrReplaceReadGroups -I ${workdir}/01b_bam_files/${basename_array}_cleaned_sorted.bam -O ${workdir}/01b_bam_files/${basename_array}_cleaned_sorted_rg.bam --RGLB 1 --RGPL illumina --RGPU unit1 --RGSM ${basename_array}

# remove cleaned and sorted bam file
rm ${workdir}/01b_bam_files/${basename_array}_cleaned_sorted.bam

# remove duplicates to sorted, cleaned, and read grouped bam file (creates final bam file)
singularity exec $SINGULARITY_CACHEDIR/gatk_4.2.3.0.sif gatk MarkDuplicates --REMOVE_DUPLICATES true --MAX_FILE_HANDLES_FOR_READ_ENDS_MAP 100 -M ${workdir}/01b_bam_files/${basename_array}_markdups_metric_file.txt -I ${workdir}/01b_bam_files/${basename_array}_cleaned_sorted_rg.bam -O ${workdir}/01b_bam_files/${basename_array}_final.bam

# remove sorted, cleaned, and read grouped bam file
rm ${workdir}/01b_bam_files/${basename_array}_cleaned_sorted_rg.bam

# index the final bam file
samtools index ${workdir}/01b_bam_files/${basename_array}_final.bam
