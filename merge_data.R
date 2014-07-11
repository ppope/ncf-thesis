

process_450k_data <- function(meth.450k.data, hist){
  
  
  #Add 450k platform indicator
  N <- nrow(meth.450k.data)
  
  platform.indicator <- rep("1", N)
  ADC.or.SCC <- rep(hist, N)
  meth.450k.data[[1]] <- ADC.or.SCC
  names(meth.450k.data)[1] <- "ADC.or.SCC"
  
  probe.data <- meth.450k.data[, -(1:4)]
  clinical.data <- meth.450k.data[, (1:4)]
  clinical.data <- cbind(clinical.data, platform.indicator)
  merged.450k.data <- cbind(clinical.data, probe.data)
  return(merged.450k.data)
  
}


merge_meth_dfrms <- function(meth.27k.data, meth.450k.data){
  
  probes.450k <- names(meth.450k.data)[-(1:5)]
  probe.data.27k <- meth.27k.data[, -(1:5)]
  clinical.data.27k <- meth.27k.data[, (1:5)]
  N <- length(probes.450k)
  M <- nrow(meth.27k.data)
  new.27k.probe.dfrm <- as.data.frame(matrix(NA, M, N))
  names(new.27k.probe.dfrm) <- probes.450k
  ind <- match(names(probe.data.27k), probes.450k)
  new.27k.probe.dfrm[,ind] <- probe.data.27k
  meth.27k.data <- cbind(clinical.data.27k, new.27k.probe.dfrm)  
  meth.27k.data[,5] <- as.factor(meth.27k.data[,5])
  total.meth.data <- rbind(meth.450k.data, meth.27k.data)
  
  return(total.meth.data)
}


merge_SNP_anno <- function(SNP.450k.anno, SNP.27k.anno){
  
  SNP.450k.anno <- cbind(SNP.450k.anno, rep(1, nrow(SNP.450k.anno)))
  names(SNP.450k.anno)[4] <- "platform.indicator"
  SNP.27k.anno <- cbind(SNP.27k.anno, rep(0, nrow(SNP.27k.anno)))
  names(SNP.27k.anno)[4] <- "platform.indicator"
  
  total.anno <- rbind(SNP.450k.anno, SNP.27k.anno)
  return(total.anno)
}

match_SNP_to_meth <- function(meth.data, SNP.anno){
  
  ind <- match(as.character(SNP.anno$TCGA.ID), as.character(meth.data$sample.name))
  meth.data <- meth.data[ind,]
  return(meth.data)
  
}

SNP_filename_parse <- function(string){
  
  x <- strsplit(string, split = "\\.")[[1]][1]
  return(x)
  
}



#This chunk reads in and merges the methylation data for 27k and 450k.
setwd("D:/Phil/find_methylation_status/output")
ADC.450k.meth.data <- read.csv("find_methylation_status_ADC_methylated_data.csv", header=TRUE)
SCC.450k.meth.data <- read.csv("find_methylation_status_SCC_methylated_data.csv", header=TRUE)
setwd("D:/Phil/methylation_27k/output")
ADC.27k.meth.data <- read.csv("ADC_27k_meth_data.csv", header=TRUE)
ADC.450k.meth.data <- process_450k_data(ADC.450k.meth.data, "ADC")
ADC.total.meth.data <- merge_meth_dfrms(ADC.27k.meth.data, ADC.450k.meth.data)
SCC.27k.meth.data <- read.csv("SCC_27k_meth_data.csv", header=TRUE)
SCC.450k.meth.data <- process_450k_data(SCC.450k.meth.data, "SCC")
SCC.total.meth.data <- merge_meth_dfrms(SCC.27k.meth.data, SCC.450k.meth.data)


#this chunk merges the 27k and 450k SNP to TCGA ID files.
setwd("D:/Phil/methylation_27k/data")
ADC.450k.SNP.anno <- read.csv("ADC_SNP_files_to_TCGA_ID.csv", header=TRUE)
ADC.27k.SNP.anno <- read.csv("ADC_27k_SNP_files_to_TCGA_ID.csv", header=TRUE)
ADC.SNP.anno <- merge_SNP_anno(ADC.450k.SNP.anno, ADC.27k.SNP.anno)
SCC.450k.SNP.anno <- read.csv("SCC_SNP_files_to_TCGA_ID.csv", header=TRUE)
SCC.27k.SNP.anno <- read.csv("SCC_27k_SNP_files_to_TCGA_ID.csv", header=TRUE)
SCC.SNP.anno <- merge_SNP_anno(SCC.450k.SNP.anno, SCC.27k.SNP.anno)



#This chunk shaves down the total.meth.data files to include only those samples that have SNP data.
ADC.total.meth.data2 <- match_SNP_to_meth(ADC.total.meth.data, ADC.SNP.anno)
SCC.total.meth.data2 <- match_SNP_to_meth(SCC.total.meth.data, SCC.SNP.anno)


setwd("D:/Phil/methylation_27k/output")
write.csv(ADC.total.meth.data2, "ADC_total_meth_data.csv", row.names=FALSE)
write.csv(SCC.total.meth.data2, "SCC_total_meth_data.csv", row.names=FALSE)
write.csv(ADC.SNP.anno, "ADC_total_SNP_files_to_TCGA_ID.csv", row.names=FALSE)
write.csv(SCC.SNP.anno, "SCC_total_SNP_files_to_TCGA_ID.csv", row.names=FALSE)



#This chunk writes all SNP data to one table.
load("D:/Phil/SNPs/less_big_ADC_SNP_data.RData")
ADC.SNP.data <- data
ADC.SNP.data$Name <- as.character(ADC.SNP.data$Name)
ADC.SNP.data$Name <- sapply(ADC.SNP.data$Name, SNP_filename_parse)
load("D:/Phil/SNPs/less_big_SCC_SNP_data.RData")
SCC.SNP.data <- data
SCC.SNP.data <- SCC.SNP.data[,-2]
SCC.SNP.data$Name <- as.character(SCC.SNP.data$Name)
SCC.SNP.data$Name <- sapply(SCC.SNP.data$Name, SNP_filename_parse)

#This chunk adds annotation to the SNP data.
names(ADC.SNP.anno)[2] <- "SNP.file.ref"
names(ADC.SNP.data)[1] <- "SNP.file.ref"
names(SCC.SNP.anno)[2] <- "SNP.file.ref"
names(SCC.SNP.data)[1] <- "SNP.file.ref"
ADC.SNP.anno.data <- merge(ADC.SNP.anno, ADC.SNP.data)
SCC.SNP.anno.data <- merge(SCC.SNP.anno, SCC.SNP.data)

#This chunk removes SNP files from tumor in the SCC dataset.
tumor.ind <- which(SCC.SNP.anno.data$SNP.file.type == "Tumor")
SCC.SNP.anno.data <- SCC.SNP.anno.data[-tumor.ind,]

write.csv(ADC.SNP.anno.data, "ADC_total_SNP_data.csv", row.names=FALSE)
write.csv(SCC.SNP.anno.data, "SCC_total_SNP_data.csv", row.names=FALSE)
