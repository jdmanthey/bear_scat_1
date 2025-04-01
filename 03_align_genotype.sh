#!/bin/bash
#SBATCH --chdir=./
#SBATCH --job-name=genotype
#SBATCH --partition nocona
#SBATCH --nodes=1 --ntasks=24
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=4G
#SBATCH --array=1-45

source activate bowtie2

threads=24

# define main working directory
workdir=/lustre/scratch/jmanthey/19_bears

basename_array=$( head -n${SLURM_ARRAY_TASK_ID} ${workdir}/basenames_bears.txt | tail -n1 )

# define the reference genome
refgenome=/home/jmanthey/references/GCA_024610735.1_mUrsAme1.0.p_genomic.fna

# run bbduk
/lustre/work/jmanthey/bbmap/bbduk.sh in1=${workdir}/00_fastq/${basename_array}_R1.fastq.gz \
in2=${workdir}/00_fastq/${basename_array}_R2.fastq.gz out1=${workdir}/01_cleaned/${basename_array}_R1.fastq.gz \
out2=${workdir}/01_cleaned/${basename_array}_R2.fastq.gz minlen=50 ftl=10 qtrim=rl trimq=10 \
ktrim=r k=25 mink=7 ref=/lustre/work/jmanthey/bbmap/resources/adapters.fa hdist=1 tbo tpe

# run bowtie2
bowtie2 --threads ${threads} -x /home/jmanthey/references/bear \
-1 ${workdir}/01_cleaned/${basename_array}_R1.fastq.gz \
-2 ${workdir}/01_cleaned/${basename_array}_R2.fastq.gz \
-S ${workdir}/01_bam_files/${basename_array}.sam

source activate bcftools

# filter for mapped reads
samtools view -b -f 2 -@ ${threads} -o ${workdir}/01_bam_files/${basename_array}.bam \
${workdir}/01_bam_files/${basename_array}.sam

# remove sam
rm ${workdir}/01_bam_files/${basename_array}.sam

# clean up the bam file
java -jar /home/jmanthey/picard.jar CleanSam I=${workdir}/01_bam_files/${basename_array}.bam \
O=${workdir}/01_bam_files/${basename_array}_cleaned.bam

# remove the raw bam
rm ${workdir}/01_bam_files/${basename_array}.bam

# sort the cleaned bam file
java -jar /home/jmanthey/picard.jar SortSam I=${workdir}/01_bam_files/${basename_array}_cleaned.bam \
O=${workdir}/01_bam_files/${basename_array}_cleaned_sorted.bam SORT_ORDER=coordinate

# remove the cleaned bam file
rm ${workdir}/01_bam_files/${basename_array}_cleaned.bam

# add read groups to sorted and cleaned bam file
java -jar /home/jmanthey/picard.jar AddOrReplaceReadGroups \
I=${workdir}/01_bam_files/${basename_array}_cleaned_sorted.bam \
O=${workdir}/01_bam_files/${basename_array}_cleaned_sorted_rg.bam \
RGLB=1 RGPL=illumina RGPU=unit1 RGSM=${basename_array}

# remove cleaned and sorted bam file
rm ${workdir}/01_bam_files/${basename_array}_cleaned_sorted.bam

# remove duplicates to sorted, cleaned, and read grouped bam file (creates final bam file)
java -jar /home/jmanthey/picard.jar MarkDuplicates \
REMOVE_DUPLICATES=true MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=100 \
M=${workdir}/01_bam_files/${basename_array}_markdups_metric_file.txt \
I=${workdir}/01_bam_files/${basename_array}_cleaned_sorted_rg.bam \
O=${workdir}/01_bam_files/${basename_array}_final.bam

# remove sorted, cleaned, and read grouped bam file
rm ${workdir}/01_bam_files/${basename_array}_cleaned_sorted_rg.bam

# index the final bam file
samtools index ${workdir}/01_bam_files/${basename_array}_final.bam

# run bcftools to genotype
bcftools mpileup --skip-indels -C 0 -d 200 --min-MQ 10 --threads ${threads} \
-f ${refgenome} ${workdir}/01_bam_files/${basename_array}_final.bam | bcftools \
call -m --threads ${threads} -o ${workdir}/02_vcf/${basename_array}.vcf

# bgzip
bgzip ${workdir}/02_vcf/${basename_array}.vcf

#tabix
tabix ${workdir}/02_vcf/${basename_array}.vcf.gz

# filter individual vcf files
bcftools view -i 'MIN(DP)>5' ${workdir}/02_vcf/${basename_array}.vcf.gz > ${workdir}/03_vcf/${basename_array}.vcf

# bgzip
bgzip ${workdir}/03_vcf/${basename_array}.vcf

#tabix
tabix ${workdir}/03_vcf/${basename_array}.vcf.gz

# stats
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

# genotype counts per individual

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




