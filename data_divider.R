



#Need to cut down list of probes to genes to include only methylated probes.

ADC.genes.to.meth.probes.list <- vector("list", N)
names(ADC.genes.to.meth.probes.list) <- names(genes.to.probes.list)

for(i in 1:N){
  
  matches <- match(genes.to.probes.list[[i]], ADC.potentially.methylated.probes)
  ind <- which(is.na(matches) == FALSE)
  ADC.genes.to.meth.probes.list[[i]] <- genes.to.probes.list[[i]][ind]
                                  
}

#Create a separate file of methylation data for each gene.

for (i in 1:N){
  
  probe.ind <- match(ADC.genes.to.meth.probes.list[[i]], names(ADC.potentially.methylated.data))
  data <- ADC.potentially.methylated.data[, c(1:4,probe.ind)]
  filename <- paste("D:/Phil/find_methylation_status/output/meth_genes/ADC_meth_genes/find_methylation_status_ADC_", names(ADC.genes.to.meth.probes.list)[i], ".csv", sep="")
  
  write.csv(data, filename, row.names=FALSE)
}





SCC.genes.to.meth.probes.list <- vector("list", N)
names(SCC.genes.to.meth.probes.list) <- names(genes.to.probes.list)

for(i in 1:N){
  
  matches <- match(genes.to.probes.list[[i]], SCC.potentially.methylated.probes)
  ind <- which(is.na(matches) == FALSE)
  SCC.genes.to.meth.probes.list[[i]] <- genes.to.probes.list[[i]][ind]
  
}

#Create a separate file of methylation data for each gene.

for (i in 1:N){
  
  probe.ind <- match(SCC.genes.to.meth.probes.list[[i]], names(SCC.potentially.methylated.data))
  data <- SCC.potentially.methylated.data[, c(1:4,probe.ind)]
  filename <- paste("D:/Phil/find_methylation_status/output/meth_genes/SCC_meth_genes/find_methylation_status_SCC_", names(SCC.genes.to.meth.probes.list)[i], ".csv", sep="")
  
  write.csv(data, filename, row.names=FALSE)
}
