
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


setwd("Z:/Leng/LUAD/JHU_USC__HumanMethylation27/Level_3/")
ADC.27k.filenames <- list.files()
ADC.27k.meth.data <- read.delim(ADC.27k.filenames[1], header=TRUE )
ADC.27k.meth.probes <- as.character(ADC.27k.meth.data$probe.name)
ADC.27k.sample.IDs <- sapply(ADC.27k.filenames, ADC_filename_27k_parser)
ADC.27k.sample.IDs <- sapply(ADC.27k.sample.IDs , TCGA_parser)
names(ADC.27k.sample.IDs) <- NULL

setwd("D:/Phil/find_methylation_status/output")
ADC.meth.data <- read.csv("find_methylation_status_ADC_methylated_data.csv", header=TRUE)
ADC.meth.probes <- names(ADC.meth.data)[-(1:4)]
ADC.common.probes <- ADC.27k.meth.probes[match(ADC.meth.probes, ADC.27k.meth.probes)]
ADC.common.probes <- ADC.common.probes[which(is.na(ADC.common.probes) == FALSE)]


setwd("Z:/Leng/LUSC/LUSC.HumanMethylation27.Level_3/")
SCC.27k.filenames <- list.files()
SCC.27k.meth.data <- read.delim(SCC.27k.filenames[1], header=TRUE, skip=1 )
SCC.27k.meth.probes <- as.character(SCC.27k.meth.data$Composite.Element.REF)
SCC.27k.sample.IDs <- sapply(SCC.27k.filenames, SCC_filename_27k_parser)
SCC.27k.sample.IDs <- sapply(SCC.27k.sample.IDs , TCGA_parser)
names(SCC.27k.sample.IDs) <- NULL

setwd("D:/Phil/find_methylation_status/output")
SCC.meth.data <- read.csv("find_methylation_status_SCC_methylated_data.csv", header=TRUE)
SCC.meth.probes <- names(SCC.meth.data)[-(1:4)]
SCC.common.probes <- SCC.27k.meth.probes[match(SCC.meth.probes, SCC.27k.meth.probes)]
SCC.common.probes <- SCC.common.probes[which(is.na(SCC.common.probes) == FALSE)]

#67/440 probes common between 450k and 27k for ADC
#63/317 probes common between 450k and 27k for SCC

#Check which have SNP data.

setwd("D:/Phil/betareg/data")
ADC.SNP.anno <- read.delim("broad.mit.edu_LUAD.Genome_Wide_SNP_6.sdrf.txt", header=TRUE)
ADC.SNP.anno <- ADC.SNP.anno[ ,c(1,7)]
ADC.SNP.anno[,1] <- sapply(ADC.SNP.anno[,1], TCGA_parser)

SCC.SNP.anno <- read.delim("broad.mit.edu_LUSC.Genome_Wide_SNP_6.sdrf.txt", header=TRUE)
SCC.SNP.anno <- SCC.SNP.anno[ ,c(2,9)]
SCC.SNP.anno[,1] <- sapply(SCC.SNP.anno[,1], TCGA_parser)



#Need to match ADC 27k sample IDs to STN, then BTN, then tumor SNP data.
sample_identifier <- function(TCGA.ID){
  x <- strsplit(TCGA.ID, "-")
  if(x[[1]][4] == "01A"){
    y <- "Tumor"
  }
  else if(x[[1]][4] == "11A"){
    y <- "STN"
  }
  else if(x[[1]][4] == "10A"){
    y <- "BTN"
  } 
  return(y)
}

find_27k_SNP_match <- function(tumor.ID, SNP.anno){
  
  patient.ID <- substring(tumor.ID, 1, 12)
  STN.ID <- paste0(patient.ID, "-11A")
  BTN.ID <- paste0(patient.ID, "-10A")
  
  x <- match(STN.ID, SNP.anno[,1])
  y <- match(BTN.ID, SNP.anno[,1])
  z <- match(tumor.ID, SNP.anno[,1])
  
  if (is.na(x) == FALSE){
    SNP.file <- SNP.anno[x,2]
    SNP.file.type <- "STN"
  }else{
      if (is.na(y) == FALSE){
        SNP.file <- SNP.anno[y,2]
        SNP.file.type <- "BTN"
      }else{
        if (is.na(z) == FALSE){
          SNP.file <- SNP.anno[z,2]
          SNP.file.type <- "Tumor"
        }else{
          SNP.file <- NA
          SNP.file.type <- NA
         }
      } 
    }
  
  return(c(tumor.ID, as.character(SNP.file), SNP.file.type))
  
}


setwd("D:/Phil/methylation_27k/output")

ADC.27k.sample.types <- sapply(ADC.27k.sample.IDs, sample_identifier)
ADC.27k.tumor.IDs <- ADC.27k.sample.IDs[which(ADC.27k.sample.types == "Tumor")]
#Take all 27k tumor IDs, and find a corresponding BTN or STN SNP file.
ADC.27k.SNP.anno <- as.data.frame(t(sapply(ADC.27k.tumor.IDs, find_27k_SNP_match, SNP.anno=ADC.SNP.anno)))
row.names(ADC.27k.SNP.anno) <- NULL
names(ADC.27k.SNP.anno) <- c("TCGA.ID", "SNP.files", "SNP.file.type")
ADC.27k.SNP.anno <- ADC.27k.SNP.anno[which(is.na(as.character(ADC.27k.SNP.anno[,2])) == FALSE),]
write.csv(ADC.27k.SNP.anno, "ADC_27k_SNP_files_to_TCGA_ID.csv", row.names=FALSE)


SCC.27k.sample.types <- sapply(SCC.27k.sample.IDs, sample_identifier)
SCC.27k.tumor.IDs <- SCC.27k.sample.IDs[which(SCC.27k.sample.types == "Tumor")]
SCC.27k.SNP.anno <- as.data.frame(t(sapply(SCC.27k.tumor.IDs, find_27k_SNP_match, SNP.anno=SCC.SNP.anno)))
row.names(SCC.27k.SNP.anno) <- NULL
names(SCC.27k.SNP.anno) <- c("TCGA.ID", "SNP.files", "SNP.file.type")
SCC.27k.SNP.anno <- SCC.27k.SNP.anno[which(is.na(as.character(SCC.27k.SNP.anno[,2])) == FALSE),]
write.csv(SCC.27k.SNP.anno, "SCC_27k_SNP_files_to_TCGA_ID.csv", row.names=FALSE)


#compare SNP files for 450k and 27k data to see if any overlap between patients.

ADC.450k.SNP.anno <- read.csv("ADC_SNP_files_to_TCGA_ID.csv")
match(ADC.27k.SNP.anno[,1], as.character(ADC.450k.SNP.anno[,1]))
#No overlap between ADC 450k and 27k SNP files.

SCC.450k.SNP.anno <- read.csv("SCC_SNP_files_to_TCGA_ID.csv")
match(SCC.27k.SNP.anno[,1], as.character(SCC.450k.SNP.anno[,1]))
#No overlap between SCC 450k and 27k SNP files.


write.csv(ADC.common.probes, "ADC_common_probes_450k_27k.csv", row.names=FALSE)
write.csv(SCC.common.probes, "SCC_common_probes_450k_27k.csv", row.names=FALSE)
