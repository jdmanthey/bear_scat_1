cat *stats >> ../_all_stats.txt

# in R:
options(scipen=999)

x <- scan("_all_stats.txt", what="character", sep="\t")

# number of individuals
n_inds <- 45

id <- x[1:n_inds * 7 - 6]
sum_aligned <- as.numeric(sapply(strsplit(x[1:n_inds * 7 - 4], "=  "), "[[", 2))
duplicates <- as.numeric(x[1:n_inds * 7 - 2])
sites <- as.numeric(x[1:n_inds * 7])


output <- data.frame(id=as.character(id), sum_aligned=as.numeric(sum_aligned), duplicates=as.numeric(duplicates), sites=as.numeric(sites))

output2 <- output[order(output$sum_aligned, decreasing=T),]
plot(output2$sites)

output3 <- output[order(output$sites, decreasing=T),]
plot(output3$sites)


output_keep <- output[output$sites >= 900000000,]

write(output_keep$id, file="bear_keep.txt", ncolumns=1)


