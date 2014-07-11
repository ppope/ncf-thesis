setwd("/home/delores/Academic/LRRI2013/Phil/SNPs/data")


ADC.anno.data <- read.delim("broad.mit.edu_LUAD.Genome_Wide_SNP_6.sdrf.txt")
ADC.anno.IDs <- as.character(ADC.anno.data[,1])
ADC.anno.SNP.files <- as.character(ADC.anno.data[,7])
ADC.anno.IDs <- substr(ADC.anno.IDs, 1, 16)

#This chunk removes samples that are duplicates. Same patient, multiple files.
dup.ind <- which(duplicated(ADC.anno.IDs) == TRUE)
ADC.anno.IDs <- ADC.anno.IDs[-dup.ind]
ADC.anno.SNP.files <- ADC.anno.SNP.files[-dup.ind]
ADC.meth.data <- read.csv("/home/delores/Academic/LRRI2013/Phil/betareg/data/find_methylation_status_ADC_methylated_data.csv")
#below three all have same indices.
ADC.meth.tumor.IDs <- as.character(ADC.meth.data[,3]) #All Sample IDs are unique
ADC.meth.patient.IDs <- substr(ADC.meth.tumor.IDs,1,12) #All Patient IDs are unique
ADC.meth.STN.IDs <- paste0(ADC.meth.patient.IDs, "-11A")

#There are three types of SNP data, STN, BTN, and tumor. 
#In order of priority we want to use STN first, then BTN, then tumor SNP data. 

N <- length(ADC.meth.tumor.IDs)
ADC.SNP.files <- rep(0,N)
ADC.SNP.file.types <- rep(0,N)

for (i in 1:N){
  
  ADC.tumor.ID <- ADC.meth.tumor.IDs[i]
  ADC.patient.ID <- substr(ADC.tumor.ID,1,12)
  ADC.STN.ID <- paste0(ADC.patient.ID, "-11A")
  ADC.BTN.ID <- paste0(ADC.patient.ID, "-10A")
  
  x <- match(ADC.STN.ID, ADC.anno.IDs)
  y <- match(ADC.BTN.ID, ADC.anno.IDs)
  z <- match(ADC.tumor.ID, ADC.anno.IDs)
  
  if (is.na(x) == FALSE){
    ADC.SNP.files[i] <- ADC.anno.SNP.files[x]
    ADC.SNP.file.types[i] <- "STN"
  }
  else{
    if (is.na(y) == FALSE){
      ADC.SNP.files[i] <- ADC.anno.SNP.files[y]
      ADC.SNP.file.types[i] <- "BTN"
    }
    else{
      if (is.na(z) == FALSE){
        ADC.SNP.files[i] <- ADC.anno.SNP.files[z]
        ADC.SNP.file.types[i] <- "tumor"
      }
      else{
        ADC.SNP.files[i] <- NA
        ADC.SNP.file.types[i] <- NA
      }
    }
  }  
}

na.ind <- which(is.na(ADC.SNP.files) == TRUE)
ADC.meth.tumor.IDs <- ADC.meth.tumor.IDs[-na.ind]
ADC.SNP.files <- ADC.SNP.files[-na.ind]
ADC.SNP.file.types <- ADC.SNP.file.types[-na.ind]
ADC.dfrm <- data.frame("TCGA.ID" = ADC.meth.tumor.IDs, "SNP.files" = ADC.SNP.files, "SNP.file.type" = ADC.SNP.file.types )

length(which(ADC.dfrm$SNP.file.type == "STN"))
length(which(ADC.dfrm$SNP.file.type == "BTN"))

#69 ADC samples have STN SNP data
#68 ADC samples have BTN SNP data



#Below essentially does the same for SCC as above.

SCC.anno.data <- read.delim("broad.mit.edu_LUSC.Genome_Wide_SNP_6.sdrf.txt")
SCC.anno.IDs <- as.character(SCC.anno.data[,2])
SCC.anno.SNP.files <- as.character(SCC.anno.data[,9])
SCC.anno.IDs <- substr(SCC.anno.IDs, 1, 16)

#This chunk removes samples that are duplicates. Same patient, multiple files.
dup.ind <- which(duplicated(SCC.anno.IDs) == TRUE)
SCC.anno.IDs <- SCC.anno.IDs[-dup.ind]
SCC.anno.SNP.files <- SCC.anno.SNP.files[-dup.ind]
SCC.meth.data <- read.csv("/home/delores/Academic/LRRI2013/Phil/betareg/data/find_methylation_status_SCC_methylated_data.csv")
#below three all have same indices.
SCC.meth.tumor.IDs <- as.character(SCC.meth.data[,3]) 
#All Sample IDs are unique
SCC.meth.patient.IDs <- substr(SCC.meth.tumor.IDs,1,12) 
#All Patient IDs are unique
SCC.meth.STN.IDs <- paste0(SCC.meth.patient.IDs, "-11A")

#There are three types of SNP data, STN, BTN, and tumor. 
#In order of priority we want to use STN first, then BTN, then tumor SNP data. 

M <- length(SCC.meth.tumor.IDs)
SCC.SNP.files <- rep(0,M)
SCC.SNP.file.types <- rep(0,M)

for (i in 1:M){
  
  SCC.tumor.ID <- SCC.meth.tumor.IDs[i]
  SCC.patient.ID <- substr(SCC.tumor.ID,1,12)
  SCC.STN.ID <- paste0(SCC.patient.ID, "-11A")
  SCC.BTN.ID <- paste0(SCC.patient.ID, "-10A")
  
  x <- match(SCC.STN.ID, SCC.anno.IDs)
  y <- match(SCC.BTN.ID, SCC.anno.IDs)
  z <- match(SCC.tumor.ID, SCC.anno.IDs)
  
  if (is.na(x) == FALSE){
    SCC.SNP.files[i] <- SCC.anno.SNP.files[x]
    SCC.SNP.file.types[i] <- "STN"
  }
  else{
    if (is.na(y) == FALSE){
      SCC.SNP.files[i] <- SCC.anno.SNP.files[y]
      SCC.SNP.file.types[i] <- "BTN"
    }
    else{
      if (is.na(z) == FALSE){
        SCC.SNP.files[i] <- SCC.anno.SNP.files[z]
        SCC.SNP.file.types[i] <- "tumor"
      }
      else{
        SCC.SNP.files[i] <- NA
        SCC.SNP.file.types[i] <- NA
      }
    }
  }  
}

na.ind <- which(is.na(SCC.SNP.files) == TRUE)

if (na.ind != 0){
  SCC.meth.tumor.IDs <- SCC.meth.tumor.IDs[-na.ind]
  SCC.SNP.files <- SCC.SNP.files[-na.ind]
  SCC.SNP.file.types <- SCC.SNP.file.types[-na.ind]
}

SCC.dfrm <- data.frame("TCGA.ID" = SCC.meth.tumor.IDs, "SNP.files" = SCC.SNP.files, "SNP.file.type" = SCC.SNP.file.types )

#94 SCC samples have STN SNP data
#55 SCC samples have BTN SNP data

length(which(SCC.dfrm$SNP.file.type == "STN"))
length(which(SCC.dfrm$SNP.file.type == "BTN"))

setwd("/home/delores/Academic/LRRI2013/Phil/betareg/data")

write.csv(ADC.dfrm, "ADC_SNP_files_to_TCGA_ID.csv", row.names=FALSE)
write.csv(SCC.dfrm, "SCC_SNP_files_to_TCGA_ID.csv", row.names=FALSE)


patient.IDs <- c(ADC.meth.patient.IDs, SCC.meth.patient.IDs)
hist <- c(rep("ADC", length(ADC.meth.patient.IDs)), rep("SCC", length(SCC.meth.patient.IDs)))

patient.IDs.dfrm <- data.frame(patient.IDs, hist)
write.csv(patient.IDs.dfrm, "patient_IDs.csv", row.names=FALSE)



