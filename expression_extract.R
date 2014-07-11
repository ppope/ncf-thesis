####Read all the expression arrays and merge them into one dataset



TCGA_parser <- function(ID){
  
  x <- substr(ID, 1, 16)
  return(x)
  
}

TCGA_patient_parser <- function(TCGA.ID){
  
  x <- substr(TCGA.ID, 1, 12)
  return(x)
  
}

extract_expression_data <- function(directory, histol){
  
  filenames <- list.files(directory, full.names=TRUE)
  N <- length(filenames)
  M <- length(sig.genes)
  data.table <- data.frame(matrix(0,N,(M+2)))
  names(data.table) <- c("sample.ID", "patient.ID", sig.genes)
  
  
  for (i in 1:N){
    
    
    current.file <- read.table(filenames[i], sep="\t", header=FALSE, na.strings ="NULL", quote=" ")
    
    if (histol == "SCC") barcode <- as.character(unlist(read.table(filenames[i], sep="\t", header=FALSE, na.strings ="NULL", nrows=1)[2]))
    if (histol == "ADC")  {
      
      barcode <- as.character(current.file[1,1])
      current.file <- current.file[,c(2,3)]
      
    }
    
    
    data.table$sample.ID[i] <- TCGA_parser(barcode)
    data.table$patient.ID[i] <- TCGA_patient_parser(data.table$sample.ID[i])
    gene.ind <- sapply(sig.genes, function(x) { which(current.file[,1] == x)} )
    data.table[i,-(1:2)] <- t(current.file[gene.ind, 2])
    
  }
  
  return(data.table)
  
}

#ADC.express.directory <- "Z:/Leng/LUAD/LUAD_UNC__AgilentG4502A_07_3_Level_3"
#ADC.express.directory <- "Z:/Leng/LUSC/LUSC.AgilentG4502A_07_3.Level_3"


SCC.express.directory <- "/home/delores/Academic/LRRI2013/data/LUSC.AgilentG4502A_07_3.Level_3"
ADC.express.directory <- "/home/delores/Academic/LRRI2013/data/LUAD_UNC__AgilentG4502A_07_3_Level_3"
sig.genes <- c("MMP2", "SYNE1", "TPM1", "SFRP2", "DKK3")

ADC.expression <- extract_expression_data(ADC.express.directory, "ADC")
SCC.expression <- extract_expression_data(SCC.express.directory, "SCC")

#setwd("D:/Phil/correlate_expression/data")
setwd("/home/delores/Academic/LRRI2013/Phil/correlate_expression/data/")

write.csv(ADC.expression,"ADC_expression_data.csv", row.names=FALSE)
write.csv(SCC.expression,"SCC_expression_data.csv", row.names=FALSE)