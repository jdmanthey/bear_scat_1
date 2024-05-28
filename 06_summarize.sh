source activate bcftools

for i in $( ls *final.bam ); do
echo $i

# samtools depth sum of aligned sites
samtools depth  $i  |  awk '{sum+=$3} END { print "Sum = ",sum}'

done


for i in $( ls *markdups* ); do
echo $i

head -n8 $i | tail -n1 | cut -f9

done

