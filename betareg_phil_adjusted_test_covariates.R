



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
  
  data <- data[which(data$platform.indicator == 1), -which(names(data) == "platform.indicator")]
  return(data)
               
}

ADC.meth.data <- remove_27k(ADC.meth.data)
ADC.covariates <- remove_27k(ADC.covariates)
SCC.meth.data <- remove_27k(SCC.meth.data)
SCC.covariates <- remove_27k(SCC.covariates)




sig.SNPs <- results.table$SNP.ID[which(results.table$Meth.Probe.Name == "cg12160578")]
sig.SNPs2 <- sapply(sig.SNPs, SNP_list_parse, USE.NAMES=FALSE)

match(sig.SNPs, names(ADC.covariates))

ADC.SNP.data2 <- ADC.covariates[, c(1, 2, match(sig.SNPs2, names(ADC.covariates)))]
ADC.meth.data2 <- ADC.meth.data[, c(1:4,which(names(ADC.meth.data) == "cg12160578"))]
ADC.covariates2 <- ADC.covariates[, 1:11]
meth.results.table <- results.table[match("cg12160578", results.table$Meth.Probe.Name), -c(4,5,7)]
SNP.results.table <- results.table[match(sig.SNPs, results.table$SNP.ID), -c(1,2,3, 9,10,11)]
SNP.results.table <- cbind(SNP.results.table, "r.squared" = rep(0, nrow(SNP.results.table)), "P.Value" = rep(0, nrow(SNP.results.table)))


SNP.data <- ADC.SNP.data2
covariates <- ADC.covariates2
reg.results.table <- SNP.results.table
hist <- "SCC"



reg_no_meth <- function(SNP.data, covariates, reg.results.table, hist){
  
  
  #This chunk sets up results table.    
  N <- nrow(reg.results.table) 
  reg.results.table$ADC.or.SCC <- rep(hist, N)
  
  
  for (i in 1:N){
    
    SNP <- SNP_list_parse(reg.results.table$SNP.ID[i])
    SNP.values <- SNP.data[ , which(names(SNP.data) == SNP), drop=FALSE]
    row.names(SNP.values) <- SNP.data$patient.ID
    cov <- covariates[c(5, 6, 8, 9, 10, 11)]
    row.names(cov) <- covariates$patient.ID
    reg.table <- merge(SNP.values, cov, by="row.names") 
    
    SNP <- unlist(reg.table[2], use.names=FALSE)
    Ethnicity_NA <- unlist(reg.table[4], use.names=FALSE)
    
    model <- lm(SNP ~ Ethnicity_NA) 
       
    reg.results.table$Estimate[i] <- model$coefficients[2]
    reg.results.table$P.Value[i] <- summary(model)$coefficients[2,4]
    reg.results.table$r.squared[i] <- summary(model)$r.squared
    
  }
  
  return(reg.results.table)   
}


beta_reg_no_SNP<- function(meth.data, covariates, reg.results.table, hist){
  
  
  #This chunk sets up results table.    
  N <- nrow(reg.results.table) 
  reg.results.table$ADC.or.SCC <- rep(hist, N)
  
  
  for (i in 1:N){
    
    meth.probe <- reg.results.table$Meth.Probe[i]
    
    BV <- meth.data[5]
    row.names(BV) <- meth.data$patient.name
    
    cov <- covariates[c(5, 6, 8, 9, 10, 11)]
    row.names(cov) <- covariates$patient.ID
    
    reg.table <- merge(BV, cov, by="row.names") 
    
    Beta_Value <- unlist(reg.table[2], use.names=FALSE)
    Ethnicity_Afro <- unlist(reg.table[3], use.names=FALSE)
    Ethnicity_NA <- unlist(reg.table[4], use.names=FALSE)
    Age <- unlist(reg.table[5], use.names=FALSE)
    Gender <- unlist(reg.table[6], use.names=FALSE)
    Smoking_Status <- as.factor(unlist(reg.table[7], use.names=FALSE))
    Stage <- unlist(reg.table[8], use.names=FALSE)
    
    model <- betareg(Beta_Value ~ Ethnicity_NA + Ethnicity_Afro  + Age + Gender + Smoking_Status + Stage, link="logit")
    
    
  }
  
  return(model)   
}


ADC.no.SNP.model <- beta_reg_no_SNP(ADC.meth.data2, ADC.covariates2, meth.results.table, "ADC")
ADC.cor.SNP.NA <- reg_no_meth(ADC.SNP.data2, ADC.covariates2, SNP.results.table, "ADC")



setwd("D:/Phil/adjusted_betareg/output")

write.csv(ADC.cor.SNP.NA, "cor_SNP_NA.csv", row.names=FALSE)








