interactive -p nocona

cd references

module load gcc/9.2.0 bwa

bwa index GCA_024610735.1_mUrsAme1.0.p_genomic.fna

samtools faidx GCA_024610735.1_mUrsAme1.0.p_genomic.fna

cd ..

java -jar picard.jar CreateSequenceDictionary R=/home/jmanthey/references/GCA_024610735.1_mUrsAme1.0.p_genomic.fna O=/home/jmanthey/references/GCA_024610735.1_mUrsAme1.0.p_genomic.dict

