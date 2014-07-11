
setwd("Z:/Leng/Phil/find_methylation_status/data")
filenames <- c("Methylation450K_TCGA_SCCvsAdeno_2_total_ordered.csv", "probenames_TCGA_SCCvsAdeno_2_Adeno_1.csv")


meth.data <- read.csv(filenames[1], header=TRUE)
symbols.to.probes <- read.csv(filenames[2], header=TRUE)


setwd("D:/Phil/SNPs/")
genes.to.symbols <- read.csv("D:/Phil/SNPs/genes_to_symbols.csv", header=TRUE)

N <- nrow(genes.to.symbols)

genes.to.symbols.list <- vector("list", length=N)
names(genes.to.symbols.list) <- as.character(genes.to.symbols[[1]])

#iterate on gene list.
#find all corresponding gene.symbols

for(i in 1:N){
  for(j in 2:4){
    if (genes.to.symbols[i,j] != "")
      genes.to.symbols.list[[i]][j-1] <- as.character(genes.to.symbols[i,j])
    else
      break
  }
}
  

#find corresponding probes to gene symbols.

genes.to.probes.list <- vector("list", length=N)
names(genes.to.probes.list) <- as.character(genes.to.symbols[[1]])

for(i in 1:N){
  M <- length(genes.to.symbols.list[[i]])
  probe.list <- NULL
  for(j in 1:M){
    probe.ind <- which(as.character(symbols.to.probes$gene.symbol) == genes.to.symbols.list[[i]][j])
    probe.list <- c(probe.list, as.character(symbols.to.probes[probe.ind,1]))
  }
  genes.to.probes.list[[i]] <- probe.list
}


#print out separate file.








