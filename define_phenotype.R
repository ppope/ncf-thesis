
setwd("D:/Phil/find_methylation_status/data")



data <- read.csv("Methylation450K_TCGA_SCCvsAdeno_2_total_ordered.csv", header=TRUE)

#split data on histology
ADC.data <- data[which(data$Adeno.or.SCC == "Adeno"), ]
SCC.data <- data[which(data$Adeno.or.SCC == "SCC"), ]

#ADC!
#split data on sample type
ADC.tumor.data <- ADC.data[which(ADC.data$sample.type == "Primary Tumor"), ]
ADC.STN.data <- ADC.data[which(ADC.data$sample.type == "Solid Tissue Normal"), ]

# find probes in Solid Tissue Normal (STN) with low methylation. These probes are candidates for being methylated in tumor.
# To do so, calculate means of each column, find which are low methylation.
# For x a methylation value, we consider x <= .2 to be low methylation, and x >= .4 to be high methylation.

ADC.STN.means <- c(rep(NA,4), apply(ADC.STN.data[,-(1:4)], 2, function(x) mean(x, na.rm=TRUE)))
ADC.STN.low.meth.ind <- which(ADC.STN.means <= .2)
ADC.tumor.high.meth.count <- apply(ADC.tumor.data[, ADC.STN.low.meth.ind], 2, function(x) length(which(x >= .40)))

#We select probes that have at least five samples with high methylation.
ADC.potentially.methylated.probes <- names(which(ADC.tumor.high.meth.count >= 5))
ADC.potentially.methylated.data <- ADC.tumor.data[,c(1:4, match(ADC.potentially.methylated.probes, names(ADC.data)))]

write.csv(ADC.potentially.methylated.data, file="D:/Phil/find_methylation_status/output/find_methylation_status_ADC_methylated_data.csv", row.names=FALSE)

#SCC!
#split data on sample type
SCC.tumor.data <- SCC.data[which(SCC.data$sample.type == "Primary Tumor"), ]
SCC.STN.data <- SCC.data[which(SCC.data$sample.type == "Solid Tissue Normal"), ]

# find probes in Solid Tissue Normal (STN) with low methylation. These probes are candidates for being methylated in tumor.
# To do so, calculate means of each column, find which are low methylation.
# For x a methylation value, we consider x <= .2 to be low methylation, and x >= .4 to be high methylation.

SCC.STN.means <- c(rep(NA,4), apply(SCC.STN.data[,-(1:4)], 2, function(x) mean(x, na.rm=TRUE)))
SCC.STN.low.meth.ind <- which(SCC.STN.means <= .2)
SCC.tumor.high.meth.count <- apply(SCC.tumor.data[, SCC.STN.low.meth.ind], 2, function(x) length(which(x >= .40)))

#We select probes that have at least five samples with high methylation.
SCC.potentially.methylated.probes <- names(which(SCC.tumor.high.meth.count >= 5))
SCC.potentially.methylated.data <- SCC.tumor.data[,c(1:4, match(SCC.potentially.methylated.probes, names(SCC.data)))]

write.csv(SCC.potentially.methylated.data, file="D:/Phil/find_methylation_status/output/find_methylation_status_SCC_methylated_data.csv", row.names=FALSE)


#Combine STN data from both histologies.
STN.probe.names <- union(names(SCC.potentially.methylated.data),names(ADC.potentially.methylated.data))
STN.data <- data[which(data$sample.type == "Solid Tissue Normal"), match(STN.probe.names, names(data))]
write.csv(STN.data, file="D:/Phil/find_methylation_status/output/find_methylation_status_STN_data.csv", row.names=FALSE)




