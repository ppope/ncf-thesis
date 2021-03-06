



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
unadjusted.reg.results <- read.csv("sig_reg_results_no_NA_afro.csv", stringsAsFactors=FALSE)


ADC.meth.data <- read.csv("ADC_total_meth_data.csv", stringsAsFactors=FALSE)
SCC.meth.data <- read.csv("SCC_total_meth_data.csv", stringsAsFactors=FALSE)
SCC.meth.data <- SCC.meth.data[-anyDuplicated(SCC.meth.data$patient.name), ]


ADC.covariates <- read.csv("ADC_total_covariates.csv", stringsAsFactors=FALSE)
SCC.covariates <- read.csv("SCC_total_covariates.csv", stringsAsFactors=FALSE)
SCC.covariates <- SCC.covariates[-anyDuplicated(SCC.covariates$patient.ID), ]

results.table <- unadjusted.reg.results[ which(unadjusted.reg.results$ADC.or.SCC == "ADC"), c(1,2,3,4,5,6,10)]

M <- nrow(results.table)
results.table <- cbind(results.table, "ADC.or.SCC" = rep("", M), "Estimate" = rep(0, M), "P.Value" = rep(0, M), "Prevalence" = rep(0, M))



#remove 27k.

remove_27k <- function(data){
  
  data <- data[which(data$platform.indicator == 1), -which(names(data) == "platform.indicator")]
  return(data)
               
}

#remove patients with NA ancestry > 0.10 and african ancestry > 0.10
remove_NA_afro <- function(covariates){
    
  data <- data[unique(c(which(covariates$proportion.african > 0.10), which(covariates$proportion.native.amer > 0.10))), ]
  return(data)
  
}


ADC.meth.data <- remove_27k(ADC.meth.data)
ADC.covariates <- remove_27k(ADC.covariates)
SCC.meth.data <- remove_27k(SCC.meth.data)
SCC.covariates <- remove_27k(SCC.covariates)

ADC.covariates <- remove_NA_afro(ADC.covariates)
SCC.covariates <- remove_NA_afro(SCC.covariates)


SCC.covariates[which(SCC.covariates$tumor_stage == "Stage II"), which(names(SCC.covariates) == "tumor_stage")] <- 1

#The function below runs beta regression on data subsets.
#Takes as parameters a list of genes, and a histology. 
#The outer loops iterates over the methylation probes and the inner loop iterates ovet the SNP probes.
#Outputs a file for each iteration of the loop.


meth.data <- SCC.meth.data
covariates <- SCC.covariates
reg.results.table <- results.table
hist <- "SCC"


beta_reg <- function(meth.data, covariates, reg.results.table, hist){
  
  
  #This chunk sets up results table.    
  N <- nrow(reg.results.table) 
  reg.results.table$ADC.or.SCC <- rep(hist, N)
  
    
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
    row.names(BV) <- meth.data$patient.name
    
    cov <- covariates[c(5, 6, 8, 9, 10, 11, test.pair.ind[[2]])]
    row.names(cov) <- covariates$patient.ID
    
    reg.table <- merge(BV, cov, by="row.names") 
    a.freq <- calc_allele_freq(reg.table[9])
    reg.results.table$Prevalence[i] <- a.freq
    reg.results.table$Sample.Size <- nrow(reg.table)
    
    #Calculate allele freq with SNP data being used in model.
    
    
    if ((a.freq == 1)|(a.freq == 0)){
      
      reg.results.table$Estimate[i] <- NA
      reg.results.table$P.Value[i] <- NA          
      next
    }
    
    
    Beta_Value <- unlist(reg.table[2], use.names=FALSE)
    Ethnicity_Afro <- unlist(reg.table[3], use.names=FALSE)
    Ethnicity_NA <- unlist(reg.table[4], use.names=FALSE)
    Age <- unlist(reg.table[5], use.names=FALSE)
    Gender <- unlist(reg.table[6], use.names=FALSE)
    Smoking_Status <- as.factor(unlist(reg.table[7], use.names=FALSE))
    Stage <- as.numeric(unlist(reg.table[8], use.names=FALSE))
    SNP <- unlist(reg.table[9], use.names=FALSE)
    
    model <- betareg(Beta_Value ~ SNP + Ethnicity_NA + Ethnicity_Afro  + Age + Gender + Smoking_Status + Stage, link="logit")
    reg.results.table$Estimate[i] <- model$coefficients$mean[2]
    reg.results.table$P.Value[i] <- summary(model)$coefficients$mean[2,4]
          
  }
    
  
  return(reg.results.table)   
}




ADC.reg.results <- beta_reg(ADC.meth.data, ADC.covariates, results.table, "ADC")
#ADC.reg.results <- ADC.reg.results[-which(is.na(ADC.reg.results$P.Value) == TRUE), ]
SCC.reg.results <- beta_reg(SCC.meth.data, SCC.covariates, results.table, "SCC")
#SCC.reg.results <- SCC.reg.results[-which(is.na(SCC.reg.results$P.Value) == TRUE), ]


#setwd("/home/delores/Academic/LRRI2013/Phil/adjusted_betareg/output")
setwd("D:/Phil/adjusted_betareg/output")


names(ADC.reg.results) == names(SCC.reg.results)

x <- rbind(ADC.reg.results, SCC.reg.results)
write.csv(x, "adjusted_reg_results_no_27k_no_NA_afro.csv", row.names=FALSE)












