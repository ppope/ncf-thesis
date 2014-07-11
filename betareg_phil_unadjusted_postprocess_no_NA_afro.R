

setwd("D:/Phil/betareg/data/")

SNP.anno <- read.csv("SNP_annotation.csv", header=TRUE)
SNP.anno <- SNP.anno[,1:2]
names(SNP.anno) <- c("SNP.ID", "SNP.Probe.Name")
 

reg.results <- read.csv("D:/Phil/betareg/output/total_reg_results_unadjusted_no_NA_afro.csv", header=TRUE)
reg.results[,3] <- sub("\\.", "-", reg.results[,3])
names(reg.results)[c(3,2)] <- c("SNP.ID", "Meth.Probe.Name")
reg.results <- merge(SNP.anno, reg.results)

meth.anno <- read.csv("annotation_TCGA_SCCvsAdeno_2.csv", header=TRUE)
meth.anno <- meth.anno[c("Name", "MAPINFO", "UCSC_RefGene_Group")]
names(meth.anno)[1] <- "Meth.Probe.Name"
reg.results <- merge(meth.anno, reg.results)

write.csv(reg.results, "D:/Phil/betareg/output/total_reg_results_unadjusted_postprocessed_no_NA_afro.csv", row.names=FALSE)