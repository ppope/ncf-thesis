

setwd("D:/Phil/significant_probes/data")
reg.results <- read.csv("total_reg_results_unadjusted_postprocessed_no_NA_afro.csv", header=TRUE)

#Only want significant ADC and SCC pairs.

#Only want probe-SNP pairs with significance at 0.05 and prevalence in [0.05,0.95]
ind1 <- which(reg.results$P.Value <= 0.05)
ind2 <- which(reg.results$Prevalence >= 0.05)
ind3 <- which(reg.results$Prevalence <= 0.95)

indie <- intersect(intersect(ind1, ind2), ind3)
sig.reg.results <- reg.results[indie,]
N <- nrow(sig.reg.results)
sig.ind <- rep(FALSE, N)


#This chunk finds which probe-SNP pairs have significant results for both histologies,
#and those whose SNP estimate is in the same direction for both histologies.
for (i in 1:N){
  
  meth.probe <- as.character(sig.reg.results$Meth.Probe.Name[i])
  SNP <- as.character(sig.reg.results$SNP.ID[i])
  ind <- intersect(which(sig.reg.results$Meth.Probe.Name == meth.probe), which(sig.reg.results$SNP.ID == SNP))
  if (length(ind) == 2){
    if ( sign(sig.reg.results$Estimate[ind[1]]) == sign(sig.reg.results$Estimate[ind[2]]) ) sig.ind[i] <- TRUE
  } 
}

sig.reg.results <- sig.reg.results[which(sig.ind == TRUE), ]
sig.reg.results <- sig.reg.results[order(sig.reg.results$Meth.Probe.Name, sig.reg.results$SNP.ID),]
setwd("D:/Phil/significant_probes/output")
write.csv(sig.reg.results, "sig_reg_results_no_NA_afro.csv", row.names=FALSE)



