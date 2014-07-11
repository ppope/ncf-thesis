

setwd("Z:/Leng/Phil/Methylation_SCC_vs_Adeno_2/output")

probenames <- read.csv("probenames_TCGA_SCCvsAdeno_2_Adeno_1.csv",sep=",",header=TRUE)

annotation.file <- read.csv("Z:/Leng/Phil/Methylation_SCC_vs_Adeno_2/data/HumanMethylation450_15017482_v.1.2.csv", 
                            sep=",", header=TRUE, skip=7)

names(probenames)[1] <- "Name"

merged.annotation <- merge(probenames, annotation.file)

#sort on first Chromosome, and then on MAPINFO
permute <- order(as.numeric(as.character(merged.annotation$CHR)),merged.annotation$MAPINFO)

merged.annotation2 <- merged.annotation[permute,]

write.csv(merged.annotation2,"Z:/Leng/Phil/Methylation_SCC_vs_Adeno_2/output/annotation_TCGA_SCCvsAdeno_2.csv", row.names=FALSE)

total.data <- read.csv("Z:/Leng/Phil/Methylation_SCC_vs_Adeno_2/output/Methylation450K_TCGA_SCCvsAdeno_2_total.csv", header=T)


merged.annotation2.probenames <- as.character(merged.annotation2$Name)
idx=match(merged.annotation2.probenames,colnames(total.data))
a=total.data[,c(1:5,idx)]

#check if ordering matches up
colnames(a) == c("Adeno.or.SCC", "sample.type",  "sample.name",  "bar.code", "patient.name", merged.annotation2.probenames)
# => Yes.

write.csv(a,"Z:/Leng/Phil/Methylation_SCC_vs_Adeno_2/output/Methylation450K_TCGA_SCCvsAdeno_2_total_ordered.csv", 
          row.names=FALSE)