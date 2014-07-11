
load("D:/Phil/SNPs/less_big_ADC_SNP_data.RData")
ADC.gene.SNP.list <- gene.SNPs.list

rm( list=ls()[-c(which(ls() == "ADC.gene.SNP.list")])

load("D:/Phil/SNPs/less_big_SCC_SNP_data.RData")
SCC.gene.SNP.list <- gene.SNPs.list

rm( list=ls()[-c(which(ls() == "ADC.gene.SNP.list"), which(ls() == "SCC.gene.SNP.list"))])

which((unlist(ADC.gene.SNP.list) == unlist(SCC.gene.SNP.list)) == TRUE)

#same list of SNPs