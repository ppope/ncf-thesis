
setwd("D:/Phil/find_methylation_status/data")

data <- read.csv("Methylation450K_TCGA_SCCvsAdeno_2_total_ordered.csv", header=TRUE)


tumor.dfrm <- data[which(data$sample.type == "Primary Tumor"),]
STN.dfrm <- data[which(data$sample.type == "Solid Tissue Normal"),]


tumor.data <- unlist(tumor.dfrm[,5:1680])

STN.data <- unlist(STN.dfrm[,5:1680])


library(ggplot2)

STN.hist <- qplot(STN.data, geom="histogram", bin=0.01, main="Histogram of BVs for STN data")
tumor.hist <- qplot(tumor.data, geom="histogram", bin=0.01, main="Histogram of BVs for tumor data")

