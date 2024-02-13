while read -r basename; do
	cat *${basename}*_R1_* >> ../${basename}_R1.fastq.gz
	cat *${basename}*_R2_* >> ../${basename}_R2.fastq.gz
done < basenames_bears.txt
