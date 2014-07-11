
#Find which probes are common between those identifed as methylation on 450k and the 27k.
#read in list of 27k probes.



ADC_filename_27k_parser <- function(filename){
  
  x <- strsplit(filename, "__")[[1]][3]
  return(x)
  
}

SCC_filename_27k_parser <- function(filename){
  
  x <- strsplit(filename, "\\.")[[1]][6]
  return(x)
}


TCGA_parser <- function(ID){
  
  y <- substr(ID, 1, 16)
  return(y)
  
}


TCGA_patient_parser <- function(ID){
  
  y <- substr(ID, 1, 12)
  return(y)
  
}


extract_27k_data <- function(hist){
  
  setwd("D:/Phil/methylation_27k/data")
  
  common.probes.file <- paste0(hist, "_common_probes_450k_27k.csv")
  common.probes <- read.csv(common.probes.file, header=TRUE)
  common.probes <- as.character(common.probes[,1])
  
  SNP.anno.file <- paste0(hist, "_27k_SNP_files_to_TCGA_ID.csv")
  SNP.anno <- read.csv(SNP.anno.file, header=TRUE)
  TCGA.IDs.27k <- as.character(SNP.anno[,1])
  
  if (hist == "ADC") {
    file.directory.27k <- "Z:/Leng/LUAD/JHU_USC__HumanMethylation27/Level_3/"
    file.parser <- ADC_filename_27k_parser
    hist1 <- "Adeno"
  }
  
  if (hist == "SCC") {
    file.directory.27k <- "Z:/Leng/LUSC/LUSC.HumanMethylation27.Level_3/"
    file.parser <- SCC_filename_27k_parser
    hist1 <- "SCC"
  }
  
  meth.27k.filenames <- list.files(file.directory.27k)
  meth.27k.filenames.parsed <- sapply(meth.27k.filenames, file.parser)
  meth.27k.filenames.parsed <- sapply(meth.27k.filenames.parsed, TCGA_parser)
  names(meth.27k.filenames.parsed) <- NULL
  
  
  ind <- match(TCGA.IDs.27k, meth.27k.filenames.parsed)
  meth.27k.files <- meth.27k.filenames[ind]
  
  
  N <- length(meth.27k.files)
  
  meth.27k.data <- data.frame("ADC.or.SCC" = rep(hist, N), "sample.type"= rep("Primary Tumor", N), 
                              "sample.name" = rep("", N), "patient.name" = rep("", N), "platform.indicator" = rep(0,N))
  meth.27k.data$sample.name <- as.character(meth.27k.data$sample.name)
  meth.27k.data$patient.name <- as.character(meth.27k.data$sample.name)
  
  probe.mat <- matrix(0, N, length(common.probes))
  probe.mat <- as.data.frame(probe.mat)
  names(probe.mat) <- common.probes
  
  for (i in 1:N){
    
    meth.27k.data$sample.name[i] <- TCGA_parser(file.parser(meth.27k.files[i]))
    meth.27k.data$patient.name[i] <- TCGA_patient_parser(meth.27k.data$sample.name[i])
    
    x <- paste0(file.directory.27k, meth.27k.files[i])
    
    if (hist == "ADC") {
      current.meth.data <- read.delim(x, header=TRUE) 
      probe.col <- 2
      beta.col <- 3
    }
       
    if (hist == "SCC") {
      current.meth.data <- read.delim(x, header=TRUE, skip=1)
      probe.col <- 1
      beta.col <- 2
    }
     
    
    probe.col
    
    probe.ind <- match(common.probes, current.meth.data[, probe.col])
    current.meth.data <- current.meth.data[probe.ind, ]
    probes <- as.character(current.meth.data[, probe.col])
    current.meth.data <- t(current.meth.data[, beta.col])
    colnames(current.meth.data) <- probes
    ind2 <- match(colnames(current.meth.data), names(probe.mat))
    probe.mat[i,] <- current.meth.data[ind2]
    
  }
  meth.27k.data <- cbind(meth.27k.data, probe.mat)
  return(meth.27k.data)
}



ADC.27k.meth.data <- extract_27k_data("ADC")
SCC.27k.meth.data <- extract_27k_data("SCC")

setwd("D:/Phil/methylation_27k/output")

write.csv(ADC.27k.meth.data, "ADC_27k_meth_data.csv", row.names=FALSE )
write.csv(SCC.27k.meth.data, "SCC_27k_meth_data.csv", row.names=FALSE )


