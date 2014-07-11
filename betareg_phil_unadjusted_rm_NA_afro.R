#setwd("/media/PEP_USB/LRRI2013/Phil/betareg/data/ADC_SNP_data")
#ADC.meth.filenames <- list.files("/media/PEP_USB/LRRI2013/Phil/betareg/data/meth_genes/ADC_meth_genes/", full.names=TRUE)
#ADC.SNP.filenames <- list.files("/media/PEP_USB/LRRI2013/Phil/betareg/data/ADC_SNP_data", full.names=TRUE)
#myGG <- qplot(x=reg.table[,5], y=reg.table[,1], xlab="SNP_A.1790605", ylab="cg27542552", main="BVs vs. SNP values for cg27542552 and SNP_A.1790605  ")

setwd("D:/Phil/betareg/data/")
library(betareg)


meth_filename_split <- function(string){
  
  x <- strsplit(string, split = "_")[[1]][8]
  y <- strsplit(x, split = "[.]")[[1]][1] 
  return(y)
  
}

SNP_filename_split <- function(string){
  
  x <- strsplit(string, split = "_")[[1]][4]
  y <- strsplit(x, split = "/")[[1]][2]
  return(y)
  
}

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

#contains probe_ID, probe_name, Chromosome, Allele_A, and Allele_B info.
SNP.anno <- read.csv("SNP_annotation.csv", header=TRUE)
SNP.list <- as.character(SNP.anno[,1])
SNP.anno[,1] <- sapply(SNP.list, SNP_list_parse)


ADC.meth.filenames <- list.files("D:/Phil/betareg/data/meth_genes/ADC_meth_genes/", full.names=TRUE)
ADC.SNP.filenames <- list.files("D:/Phil/betareg/data/ADC_SNP_data_clean", full.names=TRUE)

SCC.meth.filenames <- list.files("D:/Phil/betareg/data/meth_genes/SCC_meth_genes/", full.names=TRUE)
SCC.SNP.filenames <- list.files("D:/Phil/betareg/data/SCC_SNP_data_clean", full.names=TRUE)

#This chunk creates a list of genes whose corresponding data is used in the regression loop.
#Some genes have no methylation data (i.e. no probes were identifed as methylated) and must be removed.

ADC.meth.genes <- sapply(ADC.meth.filenames, meth_filename_split)
ADC.SNP.genes <- sapply(ADC.SNP.filenames, SNP_filename_split)
names(ADC.meth.genes) <- NULL
names(ADC.SNP.genes) <- NULL
ADC.genes <- intersect(ADC.meth.genes, ADC.SNP.genes)
#37 genes in the ADC dataset have both methylation and SNP data.

SCC.meth.genes <- sapply(SCC.meth.filenames, meth_filename_split)
SCC.SNP.genes <- sapply(SCC.SNP.filenames, SNP_filename_split)
names(SCC.meth.genes) <- NULL
names(SCC.SNP.genes) <- NULL
SCC.genes <- intersect(SCC.meth.genes, SCC.SNP.genes)
#35 genes in the ADC dataset have both methylation and SNP data.


#Remove patients with NA > 0.10, Afro > 0.10.

setwd("D:/Phil/adjusted_betareg/data")
ADC.covariates <- read.csv("ADC_total_covariates.csv", stringsAsFactors=FALSE)
SCC.covariates <- read.csv("SCC_total_covariates.csv", stringsAsFactors=FALSE)

ADC.rm.IDs <- ADC.covariates$patient.ID[unique(c(which(ADC.covariates$proportion.african > 0.10), which(ADC.covariates$proportion.native.amer > 0.10)))]
SCC.rm.IDs <- SCC.covariates$patient.ID[unique(c(which(SCC.covariates$proportion.african > 0.10), which(SCC.covariates$proportion.native.amer > 0.10)))]


#genes <- ADC.genes
#hist <- "ADC"
#rm.IDs <- ADC.rm.IDs

#The function below runs beta regression on data subsets.
#Takes as parameters a list of genes, and a histology. 
#The outer loops iterates over the methylation probes and the inner loop iterates ovet the SNP probes.
#Outputs a file for each iteration of the loop.
beta_reg_rm_NA_afro <- function(genes, hist, rm.IDs){
  
  for (l in 1:length(genes)){
    
    #This chunk sets up the regression table.
    current.gene <- genes[l]
    current.hist <- hist
    
    if (current.hist == "ADC"){
      meth.file <- ADC.meth.filenames[which(ADC.meth.genes == current.gene)]
      SNP.file <- ADC.SNP.filenames[which(ADC.SNP.genes == current.gene)]
    }else if (current.hist == "SCC"){
      meth.file <- SCC.meth.filenames[which(SCC.meth.genes == current.gene)]
      SNP.file <- SCC.SNP.filenames[which(SCC.SNP.genes == current.gene)]
    }
    
    meth.data <- read.csv(meth.file, header=TRUE) #first four columns are not BVs
    names(meth.data)[3] <- "Patient.ID"
    #remove IDs
    rm.ind <- match(rm.IDs, meth.data$patient.name)
    rm.ind <- rm.ind[which(is.na(rm.ind) == FALSE)]
    meth.data <- meth.data[-rm.ind, ]
    
    
    SNP.data <- read.csv(SNP.file, header=TRUE)
    
    
    #Merge on Patient.ID
    reg.table <- merge(meth.data[,-c(1,2,4)], SNP.data[,-c(1,3)])  
    rownames(reg.table) <- reg.table[,1]
    reg.table <- as.matrix(reg.table[,-1])
    meth.probes <- names(meth.data)[-(1:4)] 
    SNP.probes <-names(SNP.data)[-(1:3)]
    N <- length(meth.probes)
    M <- length(SNP.probes)  
    k <- 1
    
    reg.results.table <- as.data.frame(matrix(0,N*M,9))
    names(reg.results.table) <- c("Gene", "Probe", "SNP", "Estimate", "P.Value", "Prevalence", "Test.Allele", "ADC.or.SCC", "Sample.Size")
    
    for (i in 1:N){
      
      meth.probe.name <- meth.probes[i]
      for (j in 1:M){
        
        SNP.probe.name <- SNP.probes[j]
        Beta_Value <- reg.table[,i]
        SNP <- reg.table[,j+N] 
        a.freq <- calc_allele_freq(SNP)
        test.allele <- as.character(SNP.anno$Allele_B[match(SNP.probe.name, SNP.anno$probe_ID)])
        reg.results.table$Prevalence[k] <- a.freq
        reg.results.table$Sample.Size[k] <- nrow(reg.table)
        reg.results.table$Test.Allele[k] <- test.allele
        reg.results.table$Gene[k] <- current.gene
        reg.results.table$ADC.or.SCC[k] <- current.hist
        reg.results.table$Probe[k] <- meth.probe.name
        reg.results.table$SNP[k] <- SNP.probe.name
        
        if ((a.freq == 1)|(a.freq == 0)){
          
          reg.results.table$Estimate[k] <- NA
          reg.results.table$P.Value[k] <- NA
          k <- k+1
          
          next
        }
        
        model <- betareg(Beta_Value ~ SNP, link="logit")
        reg.results.table$Estimate[k] <- model$coefficients$mean[2]
        reg.results.table$P.Value[k] <- summary(model)$coefficients$mean[2,4]
        k <- k+1
      }
    }
    
    directory <- paste0("D:/Phil/betareg/output/", current.hist, "_", "reg_results_unadjusted_no_NA_afro/")
    filename <- paste0(current.hist, "_reg_results_", current.gene, ".csv")
    write.csv(reg.results.table, paste0(directory, filename), row.names=FALSE)
  }
}


ptm <- proc.time()
beta_reg_rm_NA_afro(ADC.genes, hist="ADC", rm.IDs=ADC.rm.IDs)
ADC.runtime <- proc.time() - ptm
ADC.runtime

#ADC 
#user  system elapsed 
#1322.65    0.10 1323.84

ptm <- proc.time()
beta_reg_rm_NA_afro(SCC.genes, hist="SCC", rm.IDs=SCC.rm.IDs)
SCC.runtime <- proc.time() - ptm
SCC.runtime

#SCC.runtime
#user  system elapsed 







