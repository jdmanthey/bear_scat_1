options(scipen=999)

library(vcfR)
library(adegenet)
library(StAMPP)


# read in the vcf, create a distance matrix, and output in a
# format compatible with splitstree

# read vcf
x <- read.vcfR("structure.recode.vcf")
	
# convert to genlight
x <- vcfR2genlight(x)
	
# give fake population names (not used anyway)
pop(x) <- a_rep@ind.names
	
# calculate Nei's distance 
x <- stamppNeisD(x, pop = FALSE)
	
# output name
outname <- "bear_distmat.phy"
	
# write output distance matrix in phylip format
stamppPhylip(distance.mat=x, file=outname)
	












