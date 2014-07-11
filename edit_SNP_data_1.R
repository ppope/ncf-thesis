
setwd("D:/Phil/SNPs/data/")

ADC.file <- read.table("broad.mit.edu_LUAD.Genome_Wide_SNP_6.sdrf.txt",sep="\t",header=T)
ADC.SNP.files.to.TCGA.ID <- ADC.file[,c(1,7)]
barcodes <- as.character(ADC.SNP.files.to.TCGA.ID[,1])
ADC.SNP.files.to.TCGA.ID[,1]  <- substring(barcodes,1,16)
# Need to make list of correspondence between TCGA ID and SNP names
# Match Extract.Names with Hybridization.Name
# There are repeated samples from some patients.

ADC.potentially.methylated.data <- read.csv("D:/Phil/find_methylation_status/output/find_methylation_status_ADC_methylated_data.csv", header=TRUE)

ADC.meth.sample.names <- as.character(ADC.potentially.methylated.data$sample.name)

#Following contains samples that have both SNP and Methylation data.
ADC.SNP.meth.samples <- ADC.meth.sample.names[which(is.na(match(ADC.meth.sample.names,ADC.SNP.files.to.TCGA.ID[,1])) == FALSE)]

ADC.SNP.to.TCGA.meth <- ADC.SNP.files.to.TCGA.ID[match(ADC.SNP.meth.samples, ADC.SNP.files.to.TCGA.ID[,1]),]

names(ADC.SNP.to.TCGA.meth) <- c("Patient.ID", "SNP.ID")

write.csv(ADC.SNP.to.TCGA.meth, "D:/Phil/SNPs/data/ADC_SNP_files_to_TCGA_ID.csv", row.names=FALSE)



SCC.file <- read.table("broad.mit.edu_LUSC.Genome_Wide_SNP_6.sdrf.txt",sep="\t",header=T)
SCC.SNP.files.to.TCGA.ID <- SCC.file[,c(2,9)]
barcodes <- as.character(SCC.SNP.files.to.TCGA.ID[,1])
SCC.SNP.files.to.TCGA.ID[,1]  <- substring(barcodes,1,16)
# Need to make list of correspondence between TCGA ID and SNP names
# Match Extract.Names with Hybridization.Name

SCC.potentially.methylated.data <- read.csv("D:/Phil/find_methylation_status/output/find_methylation_status_SCC_methylated_data.csv", header=TRUE)
SCC.meth.sample.names <- as.character(SCC.potentially.methylated.data$sample.name)

#Following contains samples that have both SNP and Methylation data.
SCC.SNP.meth.samples <- SCC.meth.sample.names[which(is.na(match(SCC.meth.sample.names,SCC.SNP.files.to.TCGA.ID[,1])) == FALSE)]
SCC.SNP.to.TCGA.meth <- SCC.SNP.files.to.TCGA.ID[match(SCC.SNP.meth.samples, SCC.SNP.files.to.TCGA.ID[,1]),]
names(SCC.SNP.to.TCGA.meth) <- c("Patient.ID", "SNP.ID")
write.csv(SCC.SNP.to.TCGA.meth, "D:/Phil/SNPs/data/SCC_SNP_files_to_TCGA_ID.csv", row.names=FALSE)