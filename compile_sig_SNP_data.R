
TCGA_patient_ID_parser <- function(ID){
  
  x <- substr(ID, 1, 12)
  return(x)
  
}

setwd("D:/Phil/methylation_27k/output")

ADC.SNP.data <- read.csv("ADC_total_SNP_data.csv", stringsAsFactors=FALSE)
ADC.SNP.data <- ADC.SNP.data[which(ADC.SNP.data$platform.indicator == 1), -4]
names(ADC.SNP.data)[1:2] <- c("SNP.file.name", "sample.ID")
ADC.SNP.data$patient.ID <- sapply(ADC.SNP.data$sample.ID, TCGA_patient_ID_parser, USE.NAMES=FALSE)
ADC.SNP.data <- ADC.SNP.data[c(2, 4443, 1, 3, 4:4442)]




SCC.SNP.data <- read.csv("SCC_total_SNP_data.csv", stringsAsFactors=FALSE)
SCC.SNP.data <- SCC.SNP.data[which(SCC.SNP.data$platform.indicator == 1), -4]
names(SCC.SNP.data)[1:2] <- c("SNP.file.name", "sample.ID")
SCC.SNP.data$patient.ID <- sapply(SCC.SNP.data$sample.ID, TCGA_patient_ID_parser, USE.NAMES=FALSE)
SCC.SNP.data <- SCC.SNP.data[c(2, 4443, 1, 3, 4:4442)]




setwd("D:/Phil/adjusted_betareg/data")
write.csv(ADC.SNP.data, "ADC_SNP_data_no_27k.csv", row.names=FALSE)
write.csv(SCC.SNP.data, "SCC_SNP_data_no_27k.csv", row.names=FALSE)