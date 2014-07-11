setwd("D:/Phil/SNPs/data/")


annotation <- read.csv("GenomeWideSNP_6.na32.annot.csv",header=T,skip=21)
gene.info <- read.csv("SNP_regions.csv", header=TRUE)

colnames(annotation)

N <- nrow(gene.info)

for (i in 1:N){
  
  annotation2 <- subset(annotation, subset = (Chromosome == gene.info$chromosome[i])) 
  SNP.data <- subset(annotation2, subset = (as.numeric(as.character(Physical.Position)) >= gene.info$left.coord.minus.100k[i]  & 
                                            as.numeric(as.character(Physical.Position)) <= gene.info$right.coord.plus.100k[i]) )
  
  gene <- as.character(gene.info$genes[i])
  filename <- paste("D:/Phil/SNPs/output/SNP_Info/", gene, "_SNP_Info.csv", sep = "")
  write.csv(SNP.data, filename, row.names=FALSE)

}


