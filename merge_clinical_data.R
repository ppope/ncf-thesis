

setwd("/home/delores/Academic/LRRI2013/Phil/methylation_27k/output")
ADC.SNP.data <- read.csv("ADC_total_SNP_data.csv")
SCC.SNP.data <- read.csv("SCC_total_SNP_data.csv")


setwd("/home/delores/Academic/LRRI2013/Phil/find_sex/output")
ADC.anno <- read.csv("ADC_anno_clinical.csv")
SCC.anno <- read.csv("SCC_anno_clinical.csv")

names(ADC.SNP.data)[1] <- "SNP.file.name"
names(SCC.SNP.data)[1] <- "SNP.file.name"

ADC.SNP.clinical.data <- merge(ADC.anno, ADC.SNP.data)
SCC.SNP.clinical.data <- merge(SCC.anno, SCC.SNP.data)


write.csv(ADC.SNP.clinical.data, "ADC_SNP_clinical_data.csv", row.names=FALSE)
write.csv(SCC.SNP.clinical.data, "SCC_SNP_clinical_data.csv", row.names=FALSE)

x1 <- which(ADC.SNP.clinical.data$tumor_stage == "null")
y1 <- which(ADC.SNP.clinical.data$gender == "null")
z1 <- which(ADC.SNP.clinical.data$age == "null")
u1 <- which(ADC.SNP.clinical.data$tobacco_smoking_history_indicator == "null")

ADC.null.ind <- unique(c(x1,y1,z1,u1))
#length(ADC.null.ind)
#Lose 33 samples (all 450k) when adding above covariates.

x2 <- which(SCC.SNP.clinical.data$tumor_stage == "null")
y2 <- which(SCC.SNP.clinical.data$gender == "null")
z2 <- which(SCC.SNP.clinical.data$age == "null")
u2 <- which(SCC.SNP.clinical.data$tobacco_smoking_history_indicator == "null")

SCC.null.ind <- unique(c(x2,y2,z2,u2))
#length(SCC.null.ind)
#Lose 23 samples (21/23 450k, 2/23 27k) when adding above covariates.
#That is OK.


