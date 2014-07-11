



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
  
  data <- data[which(data$platform.indicator == 1), -which(names(data) == platform.indicator)]
  return(data)
               
}


#The function below runs beta regression on data subsets.
#Takes as parameters a list of genes, and a histology. 
#The outer loops iterates over the methylation probes and the inner loop iterates ovet the SNP probes.
#Outputs a file for each iteration of the loop.


meth.data <- ADC.meth.data
covariates <- ADC.covariates
reg.results.table <- results.table
hist <- "ADC"


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
    
    
    BV <- meth.data[c(test.pair.ind[[1]],5)]
    row.names(BV) <- meth.data$patient.name
    
    cov <- covariates[c(6, 7, 9, 10, 11, 12, test.pair.ind[[2]])]
    row.names(cov) <- covariates$patient.ID
    
    #Remove NAs, much of the 27k data is missing for probes of interest.
    BV <- BV[which(is.na(BV[,1]) == FALSE), 1:2, drop=FALSE]
    reg.table <- merge(BV, cov, by="row.names") 
    a.freq <- calc_allele_freq(reg.table[10])
    reg.results.table$Prevalence[i] <- a.freq
    
    #Calculate allele freq with SNP data being used in model.
    
    
    if ((a.freq == 1)|(a.freq == 0)){
      
      reg.results.table$Estimate[i] <- NA
      reg.results.table$P.Value[i] <- NA          
      next
    }
    
    
    Beta_Value <- unlist(reg.table[2], use.names=FALSE)
    Platform.Indicator <- unlist(reg.table[3], use.names=FALSE)
    Ethnicity_NA <- unlist(reg.table[4], use.names=FALSE)
    Ethnicity_Afro <- unlist(reg.table[5], use.names=FALSE)
    Age <- unlist(reg.table[6], use.names=FALSE)
    Gender <- unlist(reg.table[7], use.names=FALSE)
    Smoking_Status <- as.factor(unlist(reg.table[8], use.names=FALSE))
    Smoking_Status <- unlist(reg.table[8])
    Stage <- unlist(reg.table[9], use.names=FALSE)
    SNP <- unlist(reg.table[10], use.names=FALSE)
    
    
    
    #Data with identical platform.indicator values throws error in regression.
    if (all(Platform.Indicator == Platform.Indicator[1]) == TRUE){
      
      model <- betareg(Beta_Value ~ SNP + Ethnicity_NA + Ethnicity_Afro + Age + Gender + Smoking_Status + Stage, link="logit")
      reg.results.table$Estimate[i] <- model$coefficients$mean[2]
      reg.results.table$P.Value[i] <- summary(model)$coefficients$mean[2,4]   
      
    } else{
      
      model <- betareg(Beta_Value ~ SNP + Ethnicity_NA + Ethnicity_Afro  + Age + Gender + Smoking_Status + Stage + Platform.Indicator, link="logit")
      reg.results.table$Estimate[i] <- model$coefficients$mean[2]
      reg.results.table$P.Value[i] <- summary(model)$coefficients$mean[2,4]
      
    }
       
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
write.csv(x, "adjusted_reg_results.csv", row.names=FALSE)












