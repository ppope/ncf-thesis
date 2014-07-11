



#read in list of sig probes
setwd("D:/Phil/significant_probes/output/")
sig.results <-  read.csv("sig_reg_results_no_NA_afro.csv", stringsAsFactors=FALSE)
sig.SNPs <- unique(sig.results$SNP.ID)




load("D:/Phil/SNPs/less_big_ADC_SNP_data.RData")
ADC.SNP.data.no.IDs <- data[, c(1, match(sig.SNPs, names(data)))]



load("D:/Phil/SNPs/less_big_SCC_SNP_data.RData")
ADC.SNP.data.no.IDs <- data[, c(1, match(sig.SNPs, names(data)))]




#This chunk reads in the 450k methylation data.
setwd("D:/Phil/find_methylation_status/output")
ADC.450k.meth.data <- read.csv("find_methylation_status_ADC_methylated_data.csv", header=TRUE)
SCC.450k.meth.data <- read.csv("find_methylation_status_SCC_methylated_data.csv", header=TRUE)


#read in list of sig probes
setwd("D:/Phil/significant_probes/output/")
sig.results <-  read.csv("sig_reg_results_no_NA_afro.csv", stringsAsFactors=FALSE)
sig.probes <- unique(sig.results$Meth.Probe.Name)
sig.SNPs <- unique(sig.results$SNP.ID)
sig.SNPs <- SNP_list_parse(sig.SNPs)

#read in SNP data




#this chunk shave down meth data to include only significant probes
ADC.450k.meth.data <- ADC.450k.meth.data[, c(3,4,match(sig.probes, names(ADC.450k.meth.data)))]
SCC.450k.meth.data <- SCC.450k.meth.data[, c(3,4,match(sig.probes, names(SCC.450k.meth.data)))]

setwd("D:/Phil/betareg/data")
write.csv(ADC.450k.meth.data, "ADC_total_sig_meth_data.csv", row.names=FALSE)
write.csv(SCC.450k.meth.data, "SCC_total_sig_meth_data.csv", row.names=FALSE)


