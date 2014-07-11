#Combined ADC/SCC into "Lung" tumor type
#Pool with GBM and HNT

SNP_list_parse <- function(string){
  
  x <- strsplit(string, split = "-")[[1]]
  y <- paste(x[1],x[2], sep=".")
  return(y)
  
}


#This chunk reads in the unadjusted significant results to get the list of sig probes.
setwd("D:/Phil/validate_HNT_GBM_COAD/data")
sig.results <- read.csv("sig_reg_results_no_NA_afro.csv", stringsAsFactors=FALSE)
sig.SNPs <- unique(sig.results$SNP.ID)
sig.SNPs2 <- sapply(sig.SNPs, SNP_list_parse)

#This chunk reads in the SNP data
ADC.SNP.data <- read.csv("ADC_SNP_data_no_27k.csv", stringsAsFactors=FALSE)
SCC.SNP.data <- read.csv("SCC_SNP_data_no_27k.csv", stringsAsFactors=FALSE)
ADC.SNP.data <- ADC.SNP.data[, c(1,2,match(sig.SNPs2, names(ADC.SNP.data)))]
SCC.SNP.data <- SCC.SNP.data[, c(1,2,match(sig.SNPs2, names(SCC.SNP.data)))]
GBM.SNP.data <- read.csv("GBM_SNP_data.csv", stringsAsFactors=FALSE)
HNT.SNP.data <- read.csv("HNT_SNP_data.csv", stringsAsFactors=FALSE)

#This chunk combines the SNP data
ADC.SNP.data <- cbind(ADC.SNP.data, tumor=rep("LUAD", nrow(ADC.SNP.data)))
SCC.SNP.data <- cbind(SCC.SNP.data, tumor=rep("LUSC", nrow(SCC.SNP.data)))
all.SNP.data <- rbind(ADC.SNP.data, SCC.SNP.data, cbind(GBM.SNP.data, tumor=rep("GBM", nrow(GBM.SNP.data))), cbind(HNT.SNP.data, tumor=rep("HNT", nrow(HNT.SNP.data))))


#This chunk reads in methylation data and combines it.
GBM.meth.data <- read.csv("GBM_meth_data.csv", stringsAsFactors=FALSE)
HNT.meth.data <- read.csv("HNT_meth_data.csv", stringsAsFactors=FALSE)
ADC.meth.data <- read.csv("ADC_total_sig_meth_data.csv", stringsAsFactors=FALSE)
SCC.meth.data <- read.csv("SCC_total_sig_meth_data.csv", stringsAsFactors=FALSE)
names(ADC.meth.data)[c(1,2)] <- c("sample.ID", "patient.ID")
names(SCC.meth.data)[c(1,2)] <- c("sample.ID", "patient.ID")
all.meth.data <- rbind(ADC.meth.data, SCC.meth.data, GBM.meth.data, HNT.meth.data)

all.meth.data <- all.meth.data[-which(duplicated(all.meth.data$patient.ID) == TRUE), ]
all.SNP.data <- all.SNP.data[-which(duplicated(all.SNP.data$patient.ID) == TRUE), ]

write.csv(all.SNP.data, "pooled_SNP_data.csv", row.names=FALSE)
write.csv(all.meth.data, "pooled_meth_data.csv", row.names=FALSE)
