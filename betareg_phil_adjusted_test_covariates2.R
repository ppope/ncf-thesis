

#Test association between significant CpG sites and Histology, Stage, Survival.
#Also test association between significant SNPs and Histology, Stage, Survival.


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
results.table <- unadjusted.reg.results[which(unadjusted.reg.results$ADC.or.SCC == "ADC"), -c(7,8,11) ]
N <- nrow(results.table)
results.table$Histology.Estimate <- rep(0, N)
results.table$Stage.Estimate <- rep(0, N)
results.table$Survival.Estimate <- rep(0, N)
results.table$Histology.P.Value <- rep(0, N)
results.table$Stage.P.Value <- rep(0, N)
results.table$Survival.P.Value <- rep(0, N)
results.table$Prevalence <- rep(0, N)
results.table$Sample.Size <- rep(0, N)


sig.probes <- unique(unadjusted.reg.results$Meth.Probe.Name)
sig.SNPs <- unique(unadjusted.reg.results$SNP.ID)
sig.SNPs2 <- sapply(sig.SNPs, SNP_list_parse, USE.NAMES=FALSE)


ADC.meth.data <- read.csv("ADC_total_meth_data.csv", stringsAsFactors=FALSE)
SCC.meth.data <- read.csv("SCC_total_meth_data.csv", stringsAsFactors=FALSE)
SCC.meth.data <- SCC.meth.data[-anyDuplicated(SCC.meth.data$patient.name), ]
ADC.meth.data <- ADC.meth.data[which(ADC.meth.data$platform.indicator == 1), c(3,4,match(sig.probes, names(ADC.meth.data)))]
SCC.meth.data <- SCC.meth.data[which(SCC.meth.data$platform.indicator == 1), c(3,4,match(sig.probes, names(SCC.meth.data)))]
all.meth.data <- rbind(ADC.meth.data, SCC.meth.data)
names(all.meth.data)[2] <- "patient.ID"

ADC.covariates <- read.csv("ADC_total_covariates.csv", stringsAsFactors=FALSE)
SCC.covariates <- read.csv("SCC_total_covariates.csv", stringsAsFactors=FALSE)
SCC.covariates <- SCC.covariates[-anyDuplicated(SCC.covariates$patient.ID), ]
ADC.covariates <- ADC.covariates[which(ADC.covariates$platform.indicator == 1), -c(3,4,5)]
SCC.covariates <- SCC.covariates[which(SCC.covariates$platform.indicator == 1), -c(3,4,5)]

ADC.SNP.data <- ADC.covariates[, c(1, 2, match(sig.SNPs2, names(ADC.covariates)))]
SCC.SNP.data <- SCC.covariates[, c(1, 2, match(sig.SNPs2, names(SCC.covariates)))]
all.SNP.data <- rbind(ADC.SNP.data, SCC.SNP.data)


#Recode Histology
#0 : ADC
#1 : SCC

ADC.covariates2 <- cbind(ADC.covariates[, c(1,2,9)], "Histol" = rep("0", nrow(ADC.covariates)))
SCC.covariates2 <- cbind(SCC.covariates[, c(1,2,9)], "Histol" = rep("1", nrow(SCC.covariates)))
all.covariates <- rbind(ADC.covariates2, SCC.covariates2)


#Add survival data to covariates
setwd("D:/Phil/find_sex/data/")
ADC.clinical.data <- read.csv("clinical_patient_public_LUAD.csv")
ADC.survival.data <- ADC.clinical.data[ which(ADC.clinical.data$days_to_death != "null"), c(1,5)]
SCC.clinical.data <- read.csv("clinical_patient_all_LUSC.csv")
SCC.survival.data <- SCC.clinical.data[ which(SCC.clinical.data$days_to_death != "null"), c(1,7)]
survival.data <- rbind(ADC.survival.data, SCC.survival.data)
names(survival.data)[1] <- "patient.ID"

all.covariates <- merge(all.covariates, survival.data)
all.covariates[which(all.covariates$tumor_stage == "Stage II"), 3] <- 1


meth.results.table <- results.table[match(sig.probes, results.table$Meth.Probe.Name), -c(4,5,7,8)]
SNP.results.table <- results.table[match(sig.SNPs, results.table$SNP.ID), -c(1,2,3)]



reg_SNP_clinical <- function(SNP.data, covariates, reg.results.table){
  
  
  #This chunk sets up results table.    
  N <- nrow(reg.results.table) 
   
  for (i in 1:N){
    
    SNP <- SNP_list_parse(reg.results.table$SNP.ID[i])
    SNP.values <- SNP.data[ , which(names(SNP.data) == SNP), drop=FALSE]
    row.names(SNP.values) <- SNP.data$patient.ID
    cov <- covariates[c(3,4,5)]
    row.names(cov) <- covariates$patient.ID
    reg.table <- merge(SNP.values, cov, by="row.names") 
    
    
    a.freq <- calc_allele_freq(SNP.values)
    reg.results.table$Prevalence[i] <- a.freq
    reg.results.table$Sample.Size[i] <- nrow(reg.table)
    
    if ((a.freq == 1)|(a.freq == 0)){
      
      reg.results.table$Stage.Estimate[i] <- NA
      reg.results.table$Histology.Estimate[i] <- NA
      reg.results.table$Surival.Estimate[i] <- NA
      reg.results.table$Stage.P.Value[i] <- NA
      reg.results.table$Histology.P.Value[i] <- NA
      reg.results.table$Surival.P.Value[i] <- NA       
      next
    }
    
    SNP <- unlist(reg.table[2], use.names=FALSE)
    Stage <- as.numeric(unlist(reg.table[3], use.names=FALSE))
    Histology <- as.numeric(unlist(reg.table[4], use.names=FALSE))
    Survival<- as.numeric(unlist(reg.table[5], use.names=FALSE))
    
    model <- lm(SNP ~ Stage + Histology + Survival) 
       
    reg.results.table$Stage.Estimate[i] <- model$coefficients[2]
    reg.results.table$Histology.Estimate[i] <- model$coefficients[3]
    reg.results.table$Survival.Estimate[i] <- model$coefficients[4]
    
    reg.results.table$Stage.P.Value[i] <- summary(model)$coefficients[2,4]
    reg.results.table$Histology.P.Value[i] <- summary(model)$coefficients[3,4]
    reg.results.table$Survival.P.Value[i] <- summary(model)$coefficients[4,4]
    
    
    reg.results.table$r.squared[i] <- summary(model)$r.squared
    
  }
  
  return(reg.results.table)   
}



beta_reg_meth_clinical<- function(meth.data, covariates, reg.results.table){
  
  
  #This chunk sets up results table.    
  N <- nrow(reg.results.table) 
  
  for (i in 1:N){
    
    meth.probe <- reg.results.table$Meth.Probe[i]
    BV <- meth.data[match(meth.probe, names(meth.data))]
    row.names(BV) <- meth.data$patient.ID 
    cov <- covariates[c(3,4,5)]
    row.names(cov) <- covariates$patient.ID
    reg.table <- merge(BV, cov, by="row.names") 
    reg.results.table$Sample.Size[i] <- nrow(reg.table)
  
    
    Beta_Value <- unlist(reg.table[2], use.names=FALSE)
    Stage <- as.numeric(unlist(reg.table[3], use.names=FALSE))
    Histology <- as.numeric(unlist(reg.table[4], use.names=FALSE))
    Survival<- as.numeric(unlist(reg.table[5], use.names=FALSE))
    
    model <- betareg(Beta_Value ~ Stage + Histology + Survival , link="logit") 
    
    reg.results.table$Stage.Estimate[i] <- model$coefficients$mean[2]
    reg.results.table$Histology.Estimate[i] <- model$coefficients$mean[3]
    reg.results.table$Survival.Estimate[i] <- model$coefficients$mean[4]
    
    reg.results.table$Stage.P.Value[i] <- summary(model)$coefficients$mean[2,4]
    reg.results.table$Histology.P.Value[i] <- summary(model)$coefficients$mean[3,4]
    reg.results.table$Survival.P.Value[i] <- summary(model)$coefficients$mean[4,4]
    
  }
  
  return(reg.results.table)   
}



SNP.clinical.results <- reg_SNP_clinical(all.SNP.data, all.covariates, SNP.results.table)
meth.clinical.results <- beta_reg_meth_clinical(all.meth.data, all.covariates, meth.results.table)



setwd("D:/Phil/adjusted_betareg/output")

write.csv(SNP.clinical.results, "SNP_clinical_linreg_results.csv", row.names=FALSE)
write.csv(meth.clinical.results, "meth_clinical_betareg_results.csv", row.names=FALSE)








