

#read in list of AIMs
setwd("D:/Phil/SNPs2/data/AIMs/")
AIM.file <- read.csv("genome1000.selected.csv")
AIMs <- names(AIM.file)[-1]

#Match rs name with SNP.ID
setwd("D:/Phil/SNPs2/data/")
SNP.anno <- read.csv("GenomeWideSNP_6.na32.annot.Allene.csv")
SNP.IDs <- as.character(SNP.anno$probe_ID[match(AIMs, SNP.anno$Probe_name)])
directory <- "Z:/Leng/LUSC/LUSC.Genome_Wide_SNP_6.Level_2"

extract_SNP_data <- function(SNP.IDs, directory){
  
  setwd(directory)
  filenames <- list.files()
  #preallocate matrix for SNP data.
  N <- length(filenames)
  M <- length(SNP.IDs)
  SNP.data <- matrix(0, N, M)
  colnames(SNP.data) <- SNP.IDs
  
  for (i in 1:N){
    
    #Columns of files are: Composite Element REF, Call, Confidence
    current.file <- read.table(filenames[i], sep = "\t", header = F, skip = 2)
    trans.file <- t(current.file[,c(1,2)])
    SNP.data[i,] <- as.numeric(trans.file[2, match(SNP.IDs, trans.file[1,])])
    
  }
  
  return(SNP.data)
  
}

SNP_file_split <- function(filename){
  
  x <- strsplit(filename, "\\.")[[1]][1]
  
}

ptm <- proc.time()
SCC.AIM.SNP.data <- extract_SNP_data(SNP.IDs, directory)
time <- proc.time() - ptm
time

#ADC
#    user   system  elapsed 
#12178.50    77.53 12409.93 

#SCC
#user   system  elapsed 
#16748.70    73.74 17138.48 

filenames <- list.files(directory)
filenames.parsed <- sapply(filenames,SNP_file_split)
names(filenames.parsed) <- NULL


SCC.AIM.SNP.data2 <- as.data.frame(cbind("Name" = filenames.parsed, SCC.AIM.SNP.data))

setwd("D:/Phil/SNPs2/output/")
write.csv(SCC.AIM.SNP.data2, "SCC_AIM_SNP_data.csv", row.names=FALSE)




