
library(betareg)




#This chunk finds probes common between those identified as significant and those in the 27k dataset.
setwd("D:/Phil/adjusted_betareg/data")
unadjusted.reg.results <- read.csv("sig_reg_results.csv", stringsAsFactors=FALSE)
sig.probes <- unique(unadjusted.reg.results$Meth.Probe.Name)
setwd("D:/Phil/methylation_27k/output/")
ADC.meth.27k.data <- read.csv("ADC_27k_meth_data.csv", stringsAsFactors=FALSE)
SCC.meth.27k.data <- read.csv("SCC_27k_meth_data.csv", stringsAsFactors=FALSE)
ADC.meth.27k.probes <- names(ADC.meth.27k.data)[-(1:5)]
SCC.meth.27k.probes <- names(SCC.meth.27k.data)[-(1:5)]
ADC.common.probe <- sig.probes[which( is.na(match(sig.probes, ADC.meth.27k.probes)) == FALSE)]
SCC.common.probe <- sig.probes[which( is.na(match(sig.probes, SCC.meth.27k.probes)) == FALSE)]
#same probe for both ADC and SCC.
# cg12317456 is the only probe common between those identified as significant and those in the 27k dataset.
common.probe <- ADC.common.probe
ADC.meth.27k.data <- ADC.meth.27k.data[ , c(1:5, which(names(ADC.meth.27k.data) == common.probe))]
SCC.meth.27k.data <- SCC.meth.27k.data[ , c(1:5, which(names(SCC.meth.27k.data) == common.probe))]


#This chunk preallocates the table for regression results.
new.reg.results.table <- unadjusted.reg.results[which(unadjusted.reg.results$Meth.Probe.Name == common.probe), ]
new.reg.results.table$P.Value <- 0
new.reg.results.table$Prevalence <- 0
new.reg.results.table$Estimate <- 0

#This chunk finds the SNPs to be tested corresponding to the common probe and parses it so it can be matched to SNPs in datasets.
SNP_list_parse <- function(string){
  
  x <- strsplit(string, split = "-")[[1]]
  y <- paste(x[1],x[2], sep=".")
  return(y)
  
}
common.SNP <- unique(new.reg.results.table$SNP.ID)
common.SNP <- SNP_list_parse(common.SNP)


#This chunk reads in the 27k data.
ADC.SNP.27k.data <- read.csv("ADC_total_SNP_data.csv", stringsAsFactors=FALSE)
ADC.SNP.27k.data <- ADC.SNP.27k.data[which(ADC.SNP.27k.data$platform.indicator == 0), c(1:3,which(names(ADC.SNP.27k.data) == common.SNP))]
SCC.SNP.27k.data <- read.csv("SCC_total_SNP_data.csv", stringsAsFactors=FALSE)
SCC.SNP.27k.data <- SCC.SNP.27k.data[which(SCC.SNP.27k.data$platform.indicator == 0), c(1:3,which(names(SCC.SNP.27k.data) == common.SNP))]



setwd("D:/Phil/betareg/data")

#This chunk reads in 450k meth data and combines with 27k meth data.
ADC.meth.450k.data <- read.csv("D:/Phil/betareg/data/meth_genes/ADC_meth_genes/find_methylation_status_ADC_MMP2.csv", stringsAsFactors=FALSE)
ADC.meth.450k.data <- ADC.meth.450k.data[, c(1:4,which(names(ADC.meth.450k.data) == common.probe))]
ADC.meth.450k.data$platform.indicator <- rep(1, nrow(ADC.meth.450k.data))
ADC.meth.450k.data$Adeno.or.SCC <- rep("ADC", nrow(ADC.meth.450k.data))
ADC.meth.450k.data[c(5,6)] <- ADC.meth.450k.data[c(6,5)]
names(ADC.meth.450k.data) <- names(ADC.meth.27k.data)
ADC.meth.data <- rbind(ADC.meth.450k.data, ADC.meth.27k.data)
SCC.meth.450k.data <- read.csv("D:/Phil/betareg/data/meth_genes/SCC_meth_genes/find_methylation_status_SCC_MMP2.csv", stringsAsFactors=FALSE)
SCC.meth.450k.data <- SCC.meth.450k.data[, c(1:4,which(names(SCC.meth.450k.data) == common.probe))]
SCC.meth.450k.data$platform.indicator <- rep(1, nrow(SCC.meth.450k.data))
SCC.meth.450k.data$Adeno.or.SCC <- rep("SCC", nrow(SCC.meth.450k.data))
SCC.meth.450k.data[c(5,6)] <- SCC.meth.450k.data[c(6,5)]
names(SCC.meth.450k.data) <- names(SCC.meth.27k.data)
SCC.meth.data <- rbind(SCC.meth.450k.data, SCC.meth.27k.data)


#This chunk reads in SNP data corresponding to 450k samples and combines it with the 27k SNP data.
ADC.SNP.450k.data <- read.csv("D:/Phil/betareg/data/ADC_SNP_data_clean/MMP2_ADC_SNP_data.csv", stringsAsFactors=FALSE)
ADC.SNP.450k.data <- ADC.SNP.450k.data[ , c(1:3, match(common.SNP, names(ADC.SNP.450k.data)))]
names(ADC.SNP.27k.data) <- names(ADC.SNP.450k.data)
ADC.SNP.data <- rbind(ADC.SNP.450k.data, ADC.SNP.27k.data)
SCC.SNP.450k.data <- read.csv("D:/Phil/betareg/data/SCC_SNP_data_clean/MMP2_SCC_SNP_data.csv", stringsAsFactors=FALSE)
SCC.SNP.450k.data <- SCC.SNP.450k.data[ , c(1:3, match(common.SNP, names(SCC.SNP.450k.data)))]
names(SCC.SNP.27k.data) <- names(SCC.SNP.450k.data)
SCC.SNP.data <- rbind(SCC.SNP.450k.data, SCC.SNP.27k.data)


#This chunk parses the sample IDs into patient IDs (for matching SNP to meth data) and then combines SNP and meth data.
ADC.TCGA.IDs <- ADC.SNP.data$Patient.ID
ADC.patient.IDs <- substr(ADC.TCGA.IDs, 1, 12)
ADC.SNP.data$Patient.ID <- ADC.patient.IDs
names(ADC.meth.data)[4] <- "Patient.ID"
ADC.reg.table <- merge(ADC.meth.data[c(4,5,6)], ADC.SNP.data[c(2,4)])
SCC.TCGA.IDs <- SCC.SNP.data$Patient.ID
SCC.patient.IDs <- substr(SCC.TCGA.IDs, 1, 12)
SCC.SNP.data$Patient.ID <- SCC.patient.IDs
names(SCC.meth.data)[4] <- "Patient.ID"
SCC.reg.table <- merge(SCC.meth.data[c(4,5,6)], SCC.SNP.data[c(2,4)])


calc_allele_freq <- function(SNP){
  
  x <- length(which(SNP == 0)) #AA
  y <- length(which(SNP == 1)) #Aa
  z <- length(which(SNP == 2)) #aa
  
  a.freq <- (2*z + y)/(2*(x + y + z))
  return(a.freq)
  
}


ADC.450k.reg.table <- ADC.reg.table[ which(ADC.reg.table$platform.indicator == 1), ]
ADC.27k.reg.table <- ADC.reg.table[ which(ADC.reg.table$platform.indicator == 0), ]


SCC.450k.reg.table <- SCC.reg.table[ which(SCC.reg.table$platform.indicator == 1), ]
SCC.27k.reg.table <- SCC.reg.table[ which(SCC.reg.table$platform.indicator == 0), ]

summary(SCC.450k.reg.table$cg12317456[ which(SCC.450k.reg.table$SNP_A.8319905 == 1)])
summary(SCC.450k.reg.table$cg12317456[ which(SCC.450k.reg.table$SNP_A.8319905 == 2)])

summary(SCC.27k.reg.table$cg12317456[ which(SCC.27k.reg.table$SNP_A.8319905 == 1)])
summary(SCC.27k.reg.table$cg12317456[ which(SCC.27k.reg.table$SNP_A.8319905 == 2)])



#repeat unadjusted model for that probe and SNP pair.
beta_reg_unadjusted_27k <- function(reg.table, histo){
  
  Beta_Value <- reg.table[,3]
  SNP <- reg.table[,4]
  
  model <- betareg(Beta_Value ~ SNP, link="logit")
  
  k <- which(new.reg.results.table$ADC.or.SCC == histo)
  
  new.reg.results.table$Estimate[k] <- model$coefficients$mean[2]
  new.reg.results.table$P.Value[k] <- summary(model)$coefficients$mean[2,4]
  new.reg.results.table$Prevalence[k] <- calc_allele_freq(SNP)
  
  return(new.reg.results.table)

}


new.reg.results.table <- beta_reg_unadjusted_27k(ADC.reg.table, "ADC")
new.reg.results.table <- beta_reg_unadjusted_27k(SCC.reg.table, "SCC")

new.reg.results.table <- cbind(new.reg.results.table, "platform" = rep("450k+27k",2))
old.reg.results.table <- cbind(unadjusted.reg.results[unadjusted.reg.results$Meth.Probe.Name == common.probe,], "platform" = rep("450k",2))


combined.reg.results <- rbind(new.reg.results.table, old.reg.results.table)

setwd("D:/Phil/betareg/output")
write.csv(combined.reg.results, "unadjusted_reg_results_with_27k.csv", row.names=FALSE)











