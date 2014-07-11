


#retrieve list of SNPs. The SNPs list is the same for both ADC and SCC.
load("D:/Phil/SNPs/less_big_ADC_SNP_data.RData")
SNP.list <- all.SNPs
rm(list= ls()[-which(ls() == "SNP.list")]) #warning this will clear the workspace



setwd("D:/Phil/betareg/data/")
SNP.anno <- read.csv("GenomeWideSNP_6.na32.annot.Allene.csv", header=TRUE)

ind <- match(SNP.list, SNP.anno[,1])
ind <- na.omit(ind)
SNP.anno.extract <- SNP.anno[ind,]

write.csv(SNP.anno.extract, "SNP_annotation.csv", row.names=FALSE)