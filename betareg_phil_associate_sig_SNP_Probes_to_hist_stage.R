

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
results.table$Estimate <- rep(0, N)
results.table$P.Value <- rep(0, N)
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
all.covariates[which(all.covariates$tumor_stage == "Stage II"), 3] <- 1


meth.results.table <- results.table[match(sig.probes, results.table$Meth.Probe.Name), -c(4,5,7,8)]
SNP.results.table <- results.table[match(sig.SNPs, results.table$SNP.ID), -c(1,2,3)]

histol.response <- all.covariates[4]
stage.response <- all.covariates[3]

row.names(histol.response) <- all.covariates[,2]
row.names(stage.response) <- all.covariates[,2]



#response <- histol.response
#covariates <- all.meth.data
#reg.results.table <- meth.results.table
#SNP.or.Meth <- "meth"

logreg_histo_stage <- function(response, covariates, reg.results.table, SNP.or.meth){
  
  
  #This chunk sets up results table.    
  N <- nrow(reg.results.table) 
   
  for (i in 1:N){
    
    
    if (SNP.or.meth == "SNP"){
      
      SNP <- SNP_list_parse(reg.results.table$SNP.ID[i])
      cov <- covariates[ , which(names(covariates) == SNP), drop=FALSE]
      row.names(cov) <- covariates$patient.ID    
      a.freq <- calc_allele_freq(cov)
      reg.results.table$Prevalence[i] <- a.freq
      
      
      if ((a.freq == 1)|(a.freq == 0)){
        
        reg.results.table$Estimate[i] <- NA
        reg.results.table$P.Value[i] <- NA       
        next
      }
      
    } else if (SNP.or.meth == "meth"){
      
      probe <- reg.results.table$Meth.Probe.Name[i]
      cov <- covariates[ , which(names(covariates) == probe), drop=FALSE]
      row.names(cov) <- covariates$patient.ID     
    }
    
    
    reg.table <- merge(response, cov, by="row.names") 
    reg.results.table$Sample.Size[i] <- nrow(reg.table)
    
    Response <- as.factor(unlist(reg.table[2], use.names=FALSE))
    Covariates <- unlist(reg.table[3], use.names=FALSE)
    
    model <- glm(Response ~ Covariates, family=binomial("logit"))
       
    reg.results.table$Estimate[i] <- model$coefficients[2]
    reg.results.table$P.Value[i] <- summary(model)$coefficients[2,4]    
  }
  
  return(reg.results.table)   
}



SNP.stage.results <- logreg_histo_stage(stage.response, all.SNP.data, SNP.results.table, "SNP")
SNP.histol.results <- logreg_histo_stage(histol.response, all.SNP.data, SNP.results.table, "SNP")

meth.stage.results <- logreg_histo_stage(stage.response, all.meth.data, meth.results.table, "meth")
meth.histol.results <- logreg_histo_stage(histol.response, all.meth.data, meth.results.table, "meth")



setwd("D:/Phil/adjusted_betareg/output")

write.csv(SNP.stage.results , "SNP_stage_association_results.csv", row.names=FALSE)
write.csv(SNP.histol.results, "SNP_histol_association_results.csv", row.names=FALSE)

write.csv(meth.stage.results , "meth_stage_association_results.csv", row.names=FALSE)
write.csv(meth.histol.results, "meth_histol_association_results.csv", row.names=FALSE)






