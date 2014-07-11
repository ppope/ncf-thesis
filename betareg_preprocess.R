
setwd("D:/Phil/betareg/data/")
ADC.SNP.to.patient.ID <- read.csv("ADC_SNP_files_to_TCGA_ID.csv", header=TRUE) #contains list of TCGA IDs that have both methylation and SNP data.
names(ADC.SNP.to.patient.ID) <- c("Patient.ID", "Name", "SNP.file.type")
SCC.SNP.to.patient.ID <- read.csv("SCC_SNP_files_to_TCGA_ID.csv", header=TRUE) #contains list of TCGA IDs that have both methylation and SNP data.
names(SCC.SNP.to.patient.ID) <- c("Patient.ID", "Name", "SNP.file.type")

#this block parses the ADC SNP.ID names
setwd("D:/Phil/betareg/data/ADC_SNP_data/")
ADC.SNP.data.filenames <- list.files()
ADC.SNP.file <-  read.csv(ADC.SNP.data.filenames[1], header=TRUE)
ADC.SNP.filenames <- as.character(ADC.SNP.file$Name)
SNP.split <- function(x){
  
  y <- strsplit(x, split = "\\.")[[1]][1]
  return(y)
  
}
ADC.SNP.filenames.parsed <- sapply(ADC.SNP.filenames, SNP.split)

#this block attaches patient IDs to each file of ADC SNP data.
for (i in 1:length(ADC.SNP.data.filenames)){
  
  ADC.SNP.file <-  read.csv(ADC.SNP.data.filenames[i], header=TRUE)
  ADC.SNP.file$Name <- ADC.SNP.filenames.parsed
  x <- merge(ADC.SNP.to.patient.ID, ADC.SNP.file)
  filename <- paste0("D:/Phil/betareg/data/ADC_SNP_data_clean/", ADC.SNP.data.filenames[i])
  write.csv(x, filename, row.names=FALSE)
  
}


#the following two blocks do the same for SCC

#this block parses the SCC SNP.ID names
setwd("D:/Phil/betareg/data/SCC_SNP_data/")
SCC.SNP.data.filenames <- list.files()
SCC.SNP.file <-  read.csv(SCC.SNP.data.filenames[1], header=TRUE)
SCC.SNP.filenames <- as.character(SCC.SNP.file$Name)
SCC.SNP.filenames.parsed <- sapply(SCC.SNP.filenames, SNP.split)

#this block attaches patient IDs to each file of SCC SNP data.
for (i in 1:length(SCC.SNP.data.filenames)){
  
  SCC.SNP.file <-  read.csv(SCC.SNP.data.filenames[i], header=TRUE)
  SCC.SNP.file$Name <- SCC.SNP.filenames.parsed
  x <- merge(SCC.SNP.to.patient.ID, SCC.SNP.file)
  filename <- paste0("D:/Phil/betareg/data/SCC_SNP_data_clean/", SCC.SNP.data.filenames[i])
  write.csv(x, filename, row.names=FALSE)
  
}