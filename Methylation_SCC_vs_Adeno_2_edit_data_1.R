
setwd("Z:/Leng/Phil/Methylation_SCC_vs_Adeno_2/output")
filenames <- c("Methylation450K_TCGA_SCCvsAdeno_2_Adeno_1.csv","Methylation450K_TCGA_SCCvsAdeno_2_SCC_1.csv")

#Adenoma
data1 <- read.csv(filenames[1],sep=",",header=TRUE)
data1[[1]] <- as.character(data1[[1]])
barcodes1 <- data1[[1]][3:256]
barcodes1.edited <- substr(barcodes1,35,62)
data1[[1]][3:256] <- barcodes1.edited
colnames(data1) <- NULL

#remove duplicate probes
column.names1 <- as.vector(as.matrix(data1[1,]))
duplicates1 <- which(duplicated(column.names1) == TRUE)
data1 <- data1[,-duplicates1]

#create table of probe names and corresponding gene symbols
probe.name1 <- t(data1[1,])
gene.symbol1 <- t(data1[2,])
gene.and.probe.names1 <- data.frame(probe.name1,gene.symbol1)
colnames(gene.and.probe.names1) <- NULL
#remove gene symbols from data
data1 <- data1[-2,]
sample.names1 <- c("sample.name", substr(barcodes1.edited,1,16))
patient.names1 <- c("patient.name", substr(barcodes1.edited,1,12))
data11 <- cbind(sample.names1,patient.names1,data1)
data11[1,3] <- "barcode"
names(data11) <- NULL
write.csv(gene.and.probe.names1,"Z:/Leng/Phil/Methylation_SCC_vs_Adeno_2/output/probenames_TCGA_SCCvsAdeno_2_Adeno_1.csv", row.names=FALSE)
write.csv(data11,"Z:/Leng/Phil/Methylation_SCC_vs_Adeno_2/output/Methylation450K_TCGA_SCCvsAdeno_2_Adeno_2.csv", row.names=FALSE)

#Squaemous
data2 <- read.csv(filenames[2],sep=",",header=TRUE)
data2[[1]] <- as.character(data2[[1]])
barcodes2 <- data2[[1]][3:202]
barcodes2.edited <- substr(barcodes2,46,73)
data2[[1]][3:202] <- barcodes2.edited
colnames(data2) <- NULL

#remove duplicate probes
column.names2 <- as.vector(as.matrix(data2[1,]))
duplicates2 <- which(duplicated(column.names2) == TRUE)
data2 <- data2[,-duplicates2]

#create table of probe names and corresponding gene symbols
probe.name2 <- t(data2[1,])
gene.symbol2 <- t(data2[2,])
gene.and.probe.names2 <- data.frame(probe.name2,gene.symbol2)
colnames(gene.and.probe.names2) <- NULL
#remove gene symbols from data
data2 <- data2[-2,]
sample.names2 <- c("sample.name", substr(barcodes2.edited,1,16))
patient.names2 <- c("patient.name", substr(barcodes2.edited,1,12))
data22 <- cbind(sample.names2,patient.names2,data2)
data22[1,3] <- "barcode"
names(data22) <- NULL
write.csv(gene.and.probe.names2,"Z:/Leng/Phil/Methylation_SCC_vs_Adeno_2/output/probenames_TCGA_SCCvsAdeno_2_SCC_1.csv", row.names=FALSE)
write.csv(data22,"Z:/Leng/Phil/Methylation_SCC_vs_Adeno_2/output/Methylation450K_TCGA_SCCvsAdeno_2_SCC_2.csv", row.names=FALSE)
