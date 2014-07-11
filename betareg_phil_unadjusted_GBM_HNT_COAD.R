



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
results.table <- unadjusted.reg.results[which(unadjusted.reg.results$ADC.or.SCC == "ADC"), -11]
N <- nrow(results.table)
results.table$Estimate <- rep(0, N)
results.table$P.Value <- rep(0, N)
results.table$Prevalence <- rep(0, N)
results.table$Sample.Size <- rep(0, N)
results.table$Tumor <- rep("", N)


setwd("D:/Phil/validate_HNT_GBM_COAD/data")

#This chunk reads in methylation data.
COAD.meth.data <- read.csv("COAD_meth_data.csv", stringsAsFactors=FALSE)
GBM.meth.data <- read.csv("GBM_meth_data.csv", stringsAsFactors=FALSE)
HNT.meth.data <- read.csv("HNT_meth_data.csv", stringsAsFactors=FALSE)


#This chunk reads in SNP data.
COAD.SNP.data <- read.csv("COAD_SNP_data.csv", stringsAsFactors=FALSE)
GBM.SNP.data <- read.csv("GBM_SNP_data.csv", stringsAsFactors=FALSE)
HNT.SNP.data <- read.csv("HNT_SNP_data.csv", stringsAsFactors=FALSE)


#Some patients have multiple samples.

#Remove duplicated data.
COAD.meth.data <- COAD.meth.data[-which(duplicated(COAD.meth.data$patient.ID) == TRUE),]
GBM.meth.data <- GBM.meth.data[-which(duplicated(GBM.meth.data$patient.ID) == TRUE),]
HNT.meth.data <- HNT.meth.data[-which(duplicated(HNT.meth.data$patient.ID) == TRUE),]
COAD.SNP.data <- COAD.SNP.data[-which(duplicated(COAD.SNP.data$patient.ID) == TRUE),]
GBM.SNP.data <- GBM.SNP.data[-which(duplicated(GBM.SNP.data$patient.ID) == TRUE),]
HNT.SNP.data <- HNT.SNP.data[-which(duplicated(HNT.SNP.data$patient.ID) == TRUE),]


#The function below runs beta regression on data subsets.
#Takes as parameters a list of genes, and a histology. 
#The outer loops iterates over the methylation probes and the inner loop iterates ovet the SNP probes.
#Outputs a file for each iteration of the loop.


#meth.data <- COAD.meth.data
#covariates <- COAD.SNP.data
#hist <- "COAD"
#reg.results.table <- results.table



beta_reg <- function(meth.data, covariates, reg.results.table, hist){
  
  
  #This chunk sets up results table.    
  N <- nrow(reg.results.table) 
  reg.results.table$Tumor<- rep(hist, N)
  
    
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
    
    
    BV <- meth.data[test.pair.ind[[1]]]
    row.names(BV) <- meth.data$patient.ID
    
    SNPx <- covariates[test.pair.ind[[2]]]
    row.names(SNPx) <- covariates$patient.ID
    reg.table <- merge(BV, SNPx, by="row.names") 
    
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
    model <- betareg(Beta_Value ~ SNP, link="logit")
    reg.results.table$Estimate[i] <- model$coefficients$mean[2]
    reg.results.table$P.Value[i] <- summary(model)$coefficients$mean[2,4]
      
       
  }
    
  
  return(reg.results.table)   
}


ptm <- proc.time()

COAD.reg.results <- beta_reg(COAD.meth.data, COAD.SNP.data, results.table, "COAD")
GBM.reg.results <- beta_reg(GBM.meth.data, GBM.SNP.data, results.table, "GBM")
HNT.reg.results <- beta_reg(HNT.meth.data, HNT.SNP.data, results.table, "HNT")

time <- proc.time() - ptm

setwd("D:/Phil/validate_HNT_GBM_COAD/output")


x <- rbind(GBM.reg.results, HNT.reg.results, COAD.reg.results)
write.csv(x, "unadjusted_GBM_HNT_COAD_reg_results.csv", row.names=FALSE)












