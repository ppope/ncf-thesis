



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


#setwd("/home/delores/Academic/LRRI2013/Phil/betareg/output")
#setwd("D:/Phil/betareg/output")
#unadjusted.reg.results <- read.csv("total_reg_results_unadjusted_postprocessed.csv", stringsAsFactors=FALSE)

setwd("D:/Phil/adjusted_betareg/data")
unadjusted.reg.results <- read.csv("sig_reg_results.csv", stringsAsFactors=FALSE)
STN.meth.data <- read.csv("find_methylation_status_STN_data.csv", stringsAsFactors=FALSE)
ADC.covariates <- read.csv("ADC_total_covariates.csv", stringsAsFactors=FALSE)
SCC.covariates <- read.csv("SCC_total_covariates.csv", stringsAsFactors=FALSE)
SCC.covariates <- SCC.covariates[-anyDuplicated(SCC.covariates$patient.ID), ]

all.covariates <- rbind(ADC.covariates, SCC.covariates)


results.table <- unadjusted.reg.results[ which(unadjusted.reg.results$ADC.or.SCC == "ADC"), c(1,2,3,4,5,6,10)]

M <- nrow(results.table)
results.table <- cbind(results.table, "Estimate" = rep(0, M), "P.Value" = rep(0, M), "Prevalence" = rep(0, M), "sample.size" = rep(0,M))


#The function below runs beta regression on data subsets.
#Takes as parameters a list of genes, and a histology. 
#The outer loops iterates over the methylation probes and the inner loop iterates ovet the SNP probes.
#Outputs a file for each iteration of the loop.


meth.data <- STN.meth.data
covariates <- ADC.covariates
hist <- "ADC"
reg.results.table <- results.table



beta_reg <- function(meth.data, covariates, reg.results.table){
  
  
  #This chunk sets up results table.    
  N <- nrow(reg.results.table) 

    
  for (i in 1:N){
    
    meth.probe <- reg.results.table$Meth.Probe[i]
    SNP.ID <- reg.results.table$SNP.ID[i]    
    SNP.ID <- SNP_list_parse(SNP.ID)
    test.pair.ind <- list(match(meth.probe, names(meth.data)), match(SNP.ID, names(covariates)))
        
    
    if (any(is.na(unlist(test.pair.ind) == TRUE))){
      
      reg.results.table$Estimate[i] <- NA
      reg.results.table$P.Value[i] <- NA          
      next
    }
    
    
    BV <- meth.data[c(test.pair.ind[[1]])]
    row.names(BV) <- meth.data$patient.name
    
    SNPx <- covariates[test.pair.ind[[2]]]
    row.names(SNPx) <- covariates$patient.ID
    

    reg.table <- merge(BV, SNPx, by="row.names") 
    a.freq <- calc_allele_freq(reg.table[3])
    reg.results.table$Prevalence[i] <- a.freq
    reg.results.table$sample.size[i] <- nrow(reg.table)
    
    #Calculate allele freq with SNP data being tested.
    
    
    if ((a.freq == 1)|(a.freq == 0)){
      
      reg.results.table$Estimate[i] <- NA
      reg.results.table$P.Value[i] <- NA          
      next
    }
    
    
    Beta_Value <- unlist(reg.table[2], use.names=FALSE)
    SNP <- unlist(reg.table[3], use.names=FALSE)
    
    
    model <- betareg(Beta_Value ~ SNP, link="logit")
    reg.results.table$Estimate[i] <- model$coefficients$mean[2]
    reg.results.table$P.Value[i] <- summary(model)$coefficients$mean[2,4]   
      
  }
  
  return(reg.results.table)   
}




STN.reg.results <- beta_reg(STN.meth.data, all.covariates, results.table)


#setwd("/home/delores/Academic/LRRI2013/Phil/adjusted_betareg/output")
setwd("D:/Phil/betareg/output")

write.csv(STN.reg.results, "STN_unadjusted_reg_results.csv", row.names=FALSE)












