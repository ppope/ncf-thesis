
setwd("D:/Phil/find_methylation_status/data")
meth.data <- read.csv("Methylation450K_TCGA_SCCvsAdeno_2_total_ordered.csv", header=TRUE)

probenames <- read.csv("probenames_TCGA_SCCvsAdeno_2_Adeno_1.csv", header=TRUE)
#probes listed in "probenames_TCGA_SCCvsAdeno_2_Adeno_1.csv" and "probenames_TCGA_SCCvsAdeno_2_SCC_1.csv" are identical.


ind1 <- which(probenames[,2] == "SYNE1")
ind2 <- which(probenames[,2] == "SFRP2")
ind3 <- which(probenames[,2] == "TPM1")

SYNE1.probes <- as.character(probenames[ind1,1])
SFRP2.probes <- as.character(probenames[ind2,1])
TPM1.probes <- as.character(probenames[ind3,1])

meth.data.probenames <- names(meth.data)[-(1:4)]


anno <- read.csv("annotation_TCGA_SCCvsAdeno_2.csv", header=TRUE) 
annoind <- c(which(anno$gene.symbol == "SYNE1"), which(anno$gene.symbol == "SFRP2"), which(anno$gene.symbol == "TPM1"))
cut.anno <- anno[annoind,]

setwd("D:/Phil/betareg/output")

SYNE1.data <-meth.data[, c(1:4, match(SYNE1.probes, meth.data.probenames))]
SFRP2.data <-meth.data[, c(1:4, match(SFRP2.probes, meth.data.probenames))]
TPM1.data <-meth.data[, c(1:4, match(TPM1.probes, meth.data.probenames))]

write.csv(SYNE1.data, "SYNE1_total_data.csv", row.names=FALSE)
write.csv(SFRP2.data, "SFRP2_total_data.csv", row.names=FALSE)
write.csv(TPM1.data, "TPM1_total_data.csv", row.names=FALSE)
write.csv(cut.anno, "SYNE1_SFRP2_TPM1_annotation.csv", row.names=FALSE)

