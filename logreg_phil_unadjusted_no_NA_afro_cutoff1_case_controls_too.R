


#Calculate cutoff by taking average of beta values in pooled normal + 0.2
#Less than cutoff is unmethylated, greater than is methylated.
#logistic regression

#If you have these SNPs what are the odds a particular probe is methylated?

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
results.table <- unadjusted.reg.results[which(unadjusted.reg.results$ADC.or.SCC == "ADC"), ]
names(results.table)[7] <- "OR.Estimate"
N <- nrow(results.table)
results.table$OR.Estimate <- rep(0, N)
results.table$P.Value <- rep(0, N)
results.table$Prevalence <- rep(0, N)
results.table$Sample.Size <- rep(0, N)
results.table$ADC.or.SCC <- rep("", N)

sig.probes <- unique(unadjusted.reg.results$Meth.Probe.Name)
sig.SNPs <- unique(unadjusted.reg.results$SNP.ID)
sig.SNPs <- sapply(sig.SNPs, SNP_list_parse, USE.NAMES=FALSE)

setwd("D:/Phil/betareg/data")
ADC.meth.data <- read.csv("ADC_total_sig_meth_data.csv", stringsAsFactors=FALSE)
SCC.meth.data <- read.csv("SCC_total_sig_meth_data.csv", stringsAsFactors=FALSE)


#Calculate cutoffs
STN.meth.data <- read.csv("find_methylation_status_STN_data.csv", stringsAsFactors=FALSE)
STN.meth.data <- STN.meth.data[ , match(sig.probes, names(STN.meth.data))]

meth.cutoffs1 <- colMeans(STN.meth.data) + 0.2

#These cutoffs were also tested, the one above was actually used.
meth.cutoffs2 <- rep(0.25, length(sig.probes))
names(meth.cutoffs2) <- sig.probes
ADC.meth.cutoffs3 <- apply(ADC.meth.data[-(1:2)], 2, median)
SCC.meth.cutoffs3 <- apply(SCC.meth.data[-(1:2)], 2, median)





#Recode methylation data
dichotomize_meth_data <- function(meth.data, cutoffs){
  
  N <- ncol(meth.data)
  for (i in 3:N){
    
    cutoff <- cutoffs[match(names(meth.data)[i], names(cutoffs))]
    meth.data[which(meth.data[i] < cutoff), i] <- 0
    meth.data[which(meth.data[i] >= cutoff), i] <- 1 
  }

  return(meth.data)
  
}

ADC.meth.data.dicho <- dichotomize_meth_data(ADC.meth.data, meth.cutoffs1)
SCC.meth.data.dicho <- dichotomize_meth_data(SCC.meth.data, meth.cutoffs1)


#This chunk creates case control tables.
create_case_control_table <- function(reg.table){
  
  N <- nrow(reg.table)

  case.control.dfrm$Cases <- colSums(reg.table)
  case.control.dfrm$Controls <- N - case.control.dfrm$Cases
  
}



setwd("D:/Phil/adjusted_betareg/data")
ADC.SNP.data <- read.csv("ADC_SNP_data_no_27k.csv", stringsAsFactors=FALSE)
SCC.SNP.data <- read.csv("SCC_SNP_data_no_27k.csv", stringsAsFactors=FALSE)
ADC.SNP.data <- ADC.SNP.data[, c(1,2,match(sig.SNPs, names(ADC.SNP.data)))]
SCC.SNP.data <- SCC.SNP.data[, c(1,2,match(sig.SNPs, names(SCC.SNP.data)))]

M <- length(sig.probes)
case.control.table <- data.frame("Meth.Probe" = sig.probes, "Controls"=rep(0,M), "Cases" = rep(0,M))

#meth.data <- ADC.meth.data.dicho
#SNP.data <- ADC.SNP.data
#reg.results.table <- results.table
#hist <- "ADC"

logistic_reg <- function(meth.data, SNP.data, reg.results.table, hist, case.control.dfrm){
  
  
  #This chunk sets up results table.    
  N <- nrow(reg.results.table) 
  reg.results.table$ADC.or.SCC<- rep(hist, N)
  
  
  for (i in 1:N){
    
    meth.probe <- reg.results.table$Meth.Probe.Name[i]
    SNP.ID <- reg.results.table$SNP.ID[i]    
    SNP.ID <- SNP_list_parse(SNP.ID)
    test.pair.ind <- list(match(meth.probe, names(meth.data)), match(SNP.ID, names(SNP.data)))
    
    
    if (any(is.na(unlist(test.pair.ind) == TRUE))){
      
      reg.results.table$OR.Estimate[i] <- NA
      reg.results.table$P.Value[i] <- NA          
      next
    }
    
    
    meth.binary <- meth.data[test.pair.ind[[1]]]
    row.names(meth.binary) <- meth.data$patient.name
    
    SNPx <- SNP.data[test.pair.ind[[2]]]
    row.names(SNPx) <- SNP.data$patient.ID
    reg.table <- merge(meth.binary, SNPx, by="row.names") 
    
    probe.ind <- match(meth.probe, case.control.dfrm$Meth.Probe)
    case.control.dfrm$Cases[probe.ind] <- sum(reg.table[2])
    case.control.dfrm$Controls[probe.ind] <- nrow(reg.table) - sum(reg.table[2])
    
    reg.results.table$Sample.Size[i] <- nrow(reg.table)
    a.freq <- calc_allele_freq(reg.table[3])
    reg.results.table$Prevalence[i] <- a.freq
    
    
    if ((a.freq == 1)|(a.freq == 0)){
      
      reg.results.table$OR.Estimate[i] <- NA
      reg.results.table$P.Value[i] <- NA          
      next
    }
    
    
    Meth_Binary <- unlist(reg.table[2], use.names=FALSE)
    SNP <- unlist(reg.table[3], use.names=FALSE)
    model <- glm(Meth_Binary ~ SNP, family=binomial("logit"))
    reg.results.table$OR.Estimate[i] <- exp(coef(model))[2]
    reg.results.table$P.Value[i] <- summary(model)$coefficients[2,4]
    
  }
  
  
  return(list(reg.results.table, case.control.dfrm))
}





x <- logistic_reg(dichotomize_meth_data(ADC.meth.data, meth.cutoffs1), ADC.SNP.data, results.table, "ADC", case.control.table)
ADC.logreg.results <- x[[1]]
ADC.case.controls <- x[[2]]
y <- logistic_reg(dichotomize_meth_data(SCC.meth.data, meth.cutoffs1), SCC.SNP.data, results.table, "SCC", case.control.table)
SCC.logreg.results <- y[[1]]
SCC.case.controls <- y[[2]]



setwd("D:/Phil/betareg/output/")
z <- rbind(ADC.logreg.results, SCC.logreg.results)
w <- cbind(rbind(ADC.case.controls, SCC.case.controls), "Histology" = c(rep("ADC", M), rep("SCC", M)))
u <- as.data.frame(t(meth.cutoffs1))

write.csv(z, "unadjusted_logreg_results_cutoff1.csv", row.names=FALSE)
write.csv(w, "unadjusted_logreg_case_controls_cutoff1.csv", row.names=FALSE)
write.csv(u, "logreg_cutoffs1.csv", row.names=FALSE)

