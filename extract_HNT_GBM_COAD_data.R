#Three different cancers
#Two different types of data



TCGA_parser <- function(ID){
  
   x <- substr(ID, 1, 16)
  return(x)
  
}

TCGA_patient_parser <- function(TCGA.ID){
  
  x <- substr(TCGA.ID, 1, 12)
  return(x)
  
}

SNP_filename_split <- function(string){
  
  x <- strsplit(string, split = "/")[[1]][5]
  y <- strsplit(x, split = "\\.")[[1]][1]
  return(y)
  
}

SNP_filename_split2 <- function(string){
  
  x <- strsplit(string, split = "/")[[1]][6]
  y <- strsplit(x, split = "\\.")[[1]][1]
  return(y)
  
}

extract_meth_data <- function(directory){
  
  filenames <- list.files(directory, full.names=TRUE)
  N <- length(filenames)
  M <- length(sig.probes)
  data.table <- data.frame(matrix(0,N,(M+2)))
  names(data.table) <- c("sample.ID", "patient.ID", sig.probes)
  
  
  for (i in 1:N){
    
    current.file <- read.delim(filenames[i], stringsAsFactors=FALSE)
    data.table$sample.ID[i] <- TCGA_parser(current.file$barcode[1])
    data.table$patient.ID[i] <- TCGA_patient_parser(data.table$sample.ID[i])
    data.table[i,-(1:2)] <- t(current.file[match(sig.probes, current.file$probe.name), 3])
    
  }
  
  return(data.table)
  
}


extract_SNP_data <- function(directory, path, type){
  
  filenames <- list.files(directory, full.names=TRUE)
  if (type == "HNT") ref.names <- sapply(filenames, SNP_filename_split, USE.NAMES=FALSE)
  if (type == "COAD") ref.names <- sapply(filenames, SNP_filename_split2, USE.NAMES=FALSE)
  if (type == "GBM") ref.names <- sapply(filenames, SNP_filename_split, USE.NAMES=FALSE)
  SNP.anno <- read.delim(path, stringsAsFactors=FALSE)
  SNP.anno <- SNP.anno[ ,c(2,8)]
  N <- length(filenames)
  M <- length(sig.SNPs)
  data.table <- data.frame(matrix(0,N,(M+2)))
  names(data.table) <- c("sample.ID", "patient.ID", sig.SNPs)
  
  
  for (i in 1:N){
    
    current.file <- read.delim(filenames[i], stringsAsFactors=FALSE, skip=1)
    TCGA.ID <- SNP.anno[match(ref.names[i], SNP.anno[,2]), 1]
    TCGA.ID <- TCGA_parser(TCGA.ID)
    data.table$sample.ID[i] <- TCGA.ID
    data.table$patient.ID[i] <- TCGA_patient_parser(TCGA.ID)
    data.table[i,-(1:2)] <- t(current.file[match(sig.SNPs, current.file$Composite.Element.REF), 2])
    
  }
  
  return(data.table)
  
}

#This chunk constructs the lists of significant probes and SNPs to be extracted.
setwd("D:/Phil/adjusted_betareg/data")
unadjusted.reg.results <- read.csv("sig_reg_results_no_NA_afro.csv", stringsAsFactors=FALSE)
sig.probes <- unique(unadjusted.reg.results$Meth.Probe.Name)
sig.SNPs <- unique(unadjusted.reg.results$SNP.ID)

#These chunks extracts the methylation data.
HNT.meth.directory <- "Z:/TCGA/HNT/Methylation_450K/DNA_Methylation/JHU_USC__HumanMethylation450/Level_3/"
COAD.meth.directory <- "Z:/TCGA/COAD/Methylation_450K/DNA_Methylation/JHU_USC__HumanMethylation450/Level_3/"
GBM.meth.directory <- "G:/GBM/GBM_Methylation/DNA_Methylation/JHU_USC__HumanMethylation450/Level_3"

ptm <- proc.time()
HNT.meth.data <- extract_meth_data(HNT.meth.directory)
HNT.meth.time <- proc.time() - ptm
HNT.meth.time
#user  system elapsed 
#1379.44   28.73 1535.84 

ptm <- proc.time()
COAD.meth.data <- extract_meth_data(COAD.meth.directory)
COAD.meth.time <- proc.time() - ptm
COAD.meth.time
#user  system elapsed 
#2058.86   40.31 2272.17 


ptm <- proc.time()
GBM.meth.data <- extract_meth_data(GBM.meth.directory)
GBM.meth.time <- proc.time() - ptm
GBM.meth.time
#user  system elapsed 
#859.31    5.93  866.31


#These chunks extracts the SNP data
HNT.SNP.directory <- "Z:/TCGA/HNT/SNP_HNSC"
HNT.SNP.anno.path <- "Z:/TCGA/HNT/METADATA/BI__Genome_Wide_SNP_6/broad.mit.edu_HNSC.Genome_Wide_SNP_6.sdrf.txt"
ptm <- proc.time()
HNT.SNP.data <- extract_SNP_data(HNT.SNP.directory, HNT.SNP.anno.path, "HNT")
HNT.SNP.time <- proc.time() - ptm
HNT.SNP.time
#user  system elapsed 
#1518.95   36.99 1701.02 


COAD.SNP.directory <- "D:/Phil/data/COAD/SNP"
COAD.SNP.anno.path <- "D:/Phil/data/COAD/broad.mit.edu_COAD.Genome_Wide_SNP_6.sdrf.txt"
ptm <- proc.time()
COAD.SNP.data <- extract_SNP_data(COAD.SNP.directory, COAD.SNP.anno.path, "COAD")
COAD.SNP.time <- proc.time() - ptm
COAD.SNP.time
#user  system elapsed 
#901.51    9.62  909.77 


GBM.SNP.directory <- "G:/GBM/SNP/GBM_SNP"
GBM.SNP.anno.path <- "G:/GBM/SNP/broad.mit.edu_GBM.Genome_Wide_SNP_6.sdrf.txt"
ptm <- proc.time()
GBM.SNP.data <- extract_SNP_data(GBM.SNP.directory, GBM.SNP.anno.path, "GBM")
GBM.SNP.time <- proc.time() - ptm
GBM.SNP.time
#user  system elapsed 
#3956.06   47.79 4046.93 


setwd("D:/Phil/validate_HNT_GBM_COAD/data/")

write.csv(HNT.SNP.data, "HNT_SNP_data.csv", row.names=FALSE)
write.csv(HNT.meth.data, "HNT_meth_data.csv", row.names=FALSE)
write.csv(COAD.SNP.data, "COAD_SNP_data.csv", row.names=FALSE)
write.csv(COAD.meth.data, "COAD_meth_data.csv", row.names=FALSE)
write.csv(GBM.SNP.data, "GBM_SNP_data.csv", row.names=FALSE)
write.csv(GBM.meth.data, "GBM_meth_data.csv", row.names=FALSE)



