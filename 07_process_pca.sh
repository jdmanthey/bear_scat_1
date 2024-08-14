source activate bcftools

cd /lustre/scratch/jmanthey/19_bears/05_filtered


##########################################
##########################################
### cat all individual scaffolds together
##########################################
##########################################

grep "^#" structure_JANIGQ010000097.1.recode.vcf > structure.vcf

for i in $( ls structure_*.recode.vcf ); do grep -v "^#" $i >> structure.vcf; done

grep "^#" kinship_JANIGQ010000097.1.recode.vcf > kinship.vcf

for i in $( ls kinship_*.recode.vcf ); do grep -v "^#" $i >> kinship.vcf; done

##########################################
##########################################
### remove all individual scaffold files
##########################################
##########################################

rm *recode*

##########################################
##########################################
### thin for 5 kbp separation between SNPs for PCA
### --max-missing 0.4 for kinship with a bit of thinning
##########################################
##########################################

vcftools --vcf structure.vcf --thin 5000 --recode --recode-INFO-all --out structure

vcftools --vcf kinship.vcf --max-missing 0.4 --thin 1000 --recode --recode-INFO-all --out kinship

##########################################
##########################################
### run plink for pca
##########################################
##########################################

# make chromosome map for this vcf
grep -v "#" structure.vcf | cut -f 1 | uniq | awk '{print $0"\t"$0}' > chrom_map.txt

#plink output format
vcftools --vcf structure.vcf  --plink --chrom-map chrom_map.txt --out structure 

# convert  with plink
plink --file structure --recode12 --allow-extra-chr \
--out structure_plink

# run pca on dataset
plink --file structure_plink --pca --allow-extra-chr \
--out structure_plink_pca

