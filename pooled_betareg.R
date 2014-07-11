#Combined ADC/SCC into "Lung" tumor type
#Pool with GBM and HNT

#Run betareg, adjusting for tumor type.


setwd("D:/Phil/adjusted_betareg/data")

library(betareg)

SNP_list_parse <- function(string){
  
  x <- strsplit(string, split = "-")[[1]]
  y <- paste(x[1],x[2], sep=".")
  return(y)
  
}

calc_allele_freq <- function(SNP){
  
  x <- length(which(SNP == 0)) #AA
  y <- length(which(SNP == 1)) #Aa
  z <- length(which(SNP == 2)) #aa
  
  a.freq <- (2*z + y)/(2*(x + y + z))
  return(a.freq)
  
}


setwd("D:/Phil/adjusted_betareg/data")
unadjusted.reg.results <- read.csv("sig_reg_results_no_NA_afro.csv", stringsAsFactors=FALSE)
results.table <- unadjusted.reg.results[which(unadjusted.reg.results$ADC.or.SCC == "ADC"), -c(7,8,11)]
N <- nrow(results.table)
results.table$Prevalence <- rep(0, N)
results.table$Sample.Size <- rep(0, N)

setwd("D:/Phil/validate_HNT_GBM_COAD/data")
all.meth.data <- read.csv("pooled_meth_data.csv", stringsAsFactors=FALSE)
all.SNP.data <- read.csv("pooled_SNP_data.csv", stringsAsFactors=FALSE)

meth.data <- all.meth.data
SNP.data <- all.SNP.data
reg.results.table <- results.table


pooled_betareg <- function(meth.data, SNP.data, reg.results.table){
  
  
  #This chunk sets up results table.    
  N <- nrow(reg.results.table) 
  
  
  for (i in 1:N){
    
    meth.probe <- reg.results.table$Meth.Probe[i]
    SNP.ID <- reg.results.table$SNP.ID[i]    
    SNP.ID <- SNP_list_parse(SNP.ID)
    test.pair.ind <- list(match(meth.probe, names(meth.data)), match(SNP.ID, names(SNP.data)))
    
    
    if (any(is.na(unlist(test.pair.ind) == TRUE))){
      
      reg.results.table$Estimate[i] <- NA
      reg.results.table$P.Value[i] <- NA          
      next
    }
    
    
    BV <- meth.data[test.pair.ind[[1]]]
    row.names(BV) <- meth.data$patient.ID
    
    cov <- SNP.data[c(test.pair.ind[[2]], 29)]
    row.names(cov) <- SNP.data$patient.ID
    reg.table <- merge(BV, cov, by="row.names") 
    
    reg.results.table$Sample.Size[i] <- nrow(reg.table)
    a.freq <- calc_allele_freq(reg.table[3])
    reg.results.table$Prevalence[i] <- a.freq
    
    
    
    if ((a.freq == 1)|(a.freq == 0)){
      
      reg.results.table$Estimate[i] <- NA
      reg.results.table$P.Value[i] <- NA          
      next
    }
    
    
    Beta_Value <- unlist(reg.table[2], use.names=FALSE)
    SNP <- unlist(reg.table[3], use.names=FALSE)
    Tumor <- as.factor(unlist(reg.table[4], use.names=FALSE))
    model <- betareg(Beta_Value ~ SNP + Tumor, link="logit")
    reg.results.table$SNP.Estimate[i] <- model$coefficients$mean[2]
    reg.results.table$SNP.P.Value[i] <- summary(model)$coefficients$mean[2,4]

  }
  
  
  return(reg.results.table)   
}


pooled.betareg.results <- pooled_betareg(all.meth.data, all.SNP.data, results.table)

setwd("D:/Phil/validate_HNT_GBM_COAD/output/")
write.csv(pooled.betareg.results, "pooled_betareg_results.csv", row.names=FALSE)

