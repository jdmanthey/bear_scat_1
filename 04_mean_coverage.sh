#!/bin/bash
#SBATCH --chdir=./
#SBATCH --job-name=depth1
#SBATCH --partition quanah
#SBATCH --nodes=1 --ntasks=4
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=8G

for i in $( ls *final.bam ); do 
  echo $i >> mean_depth.txt
  samtools depth -a file.bam | awk '{c++;s+=$3}END{print s/c}' >> mean_depth.txt
done

