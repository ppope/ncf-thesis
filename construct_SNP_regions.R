#This is a script to create the template for the list of regions from which to extract SNPs.
#Genomic coordinates will be filled in by hand, there could be a script to scrape the data from the web.


#define template
genes <- c("ACP5", "C2orf32", "CCDC50", "CPEB1", "CST6", "DKK3", "GJB2", "GSTT2", "HTRA3", "HTRA4", "IGDCC4", "LAYN", "LOX", "MARVELD3", "MSX1",
           "MT1M", "PAK6-2", "RBP1", "ROBO3", "TMEM22", "TPM1", "ZNF365", "BNIP3", "EMILIN2", "MSRB3", "VAV3", "AKAP12", "SYNE1", "IL20R",
           "p16", "RASSF1A", "MGMT-O6", "SFRP2", "SFRP1", "BETA3", "NOVEL2", "IGFBP3", "RPRM", "DCR1", "HLHP", "AK5", "TSLC1", "TUBB4",
           "XT3", "PCDH10", "MMP2", "APC2", "GPNMB", "BOLL", "ICAM5", "JPH3", "GNB4")

N <- length(genes)
dfrm <- data.frame("genes" = genes, "left.coord"= rep(" ", N), "right.coord"= rep(" ", N), "left.coord.minus.100k" = rep(" ", N), "right.coord.plus.100K"= rep(" ", N) )

write.csv(dfrm, file="D:/Phil/SNPs/data/SNP_regions.csv", row.names=FALSE)