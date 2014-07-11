extract_X_chromosome_data <- function(data){
  
  current.SNP.file <- data
  X.SNP.ind <- match(X.SNPs, current.SNP.file$Composite.Element.REF)
  #X.SNP.ind <- X.SNP.ind[which(is.na(X.SNP.ind) == FALSE)]
  X.SNP.data <- current.SNP.file[X.SNP.ind, ]
  X.SNP.data <- X.SNP.data[,2] 
  return( X.SNP.data)
}


extract_Y_chromosome_data <- function(data){
  
  current.SNP.file <- data
  Y.SNP.ind <- match(Y.SNPs, current.SNP.file$Composite.Element.REF)
  #Y.SNP.ind <- Y.SNP.ind[which(is.na(Y.SNP.ind) == FALSE)]
  Y.SNP.data <- current.SNP.file[Y.SNP.ind, ]
  Y.SNP.data <- Y.SNP.data[,2] 
  return( Y.SNP.data)
}



get_Y_chromosome_data <- function(directory){
  
  setwd(directory) 
  SNP.files <- list.files(directory)
  N <- length(Y.SNPs)
  M <- length(SNP.files)
  Y <- as.data.frame(matrix(0, N, M+1))
  Y[,1] <- Y.SNPs
  
  for (i in 1:M){
    
    data <- read.delim(SNP.files[i], header=TRUE, skip=1)
    y <- extract_Y_chromosome_data(data)
    Y[,i+1] <- y
  }
  
  return(Y)
}




SNP.anno <- read.csv("GenomeWideSNP_6.na32.annot.csv", header=TRUE, skip=21)
SNP.sex.anno <- SNP.anno[, c(6,22)]


Y.ind <- which(SNP.anno$Chromosome == "Y")
Y.SNPs <- as.character(SNP.anno$Probe.Set.ID[Y.ind])

ADC.SNP.directory <- "/home/delores/Academic/LRRI2013/data/broad.mit.edu_LUAD.Genome_Wide_SNP_6.Level_2"
SCC.SNP.directory <- "/home/delores/Academic/LRRI2013/data/LUSC.Genome_Wide_SNP_6.Level_2"

#ADC.SNP.directory <- "Z:/Leng/LUAD/broad.mit.edu_LUAD.Genome_Wide_SNP_6.Level_2/"
#SCC.SNP.directory <- "Z:/Leng/LUSC/LUSC.Genome_Wide_SNP_6.Level_2/"


ADC.Y <- get_Y_chromosome_data(ADC.SNP.directory)
ADC.Y <- na.omit(ADC.Y)

SCC.Y <- get_Y_chromosome_data(SCC.SNP.directory)
SCC.Y <- na.omit(SCC.Y)

setwd("D:/Phil/find_sex/")
save.image( file="xy.Rdata")

#X.pseudo.autosomal.data.ind1 <- which(SNP.sex.anno[,1] != 0) 
#X.pseudo.autosomal.data.ind2 <- which(SNP.sex.anno[,2] != 0)
#X.pseudo.autosomal1.SNPs <- as.character(SNP.anno[X.pseudo.autosomal.data.ind1,1])
#X.pseudo.autosomal2.SNPs <- as.character(SNP.anno[X.pseudo.autosomal.data.ind2,1])
#match(X.pseudo.autosomal1.SNPs, ADC.Y[,1])
#match(X.pseudo.autosomal2.SNPs, SCC.Y[,1])
#match(X.pseudo.autosomal1.SNPs, SCC.Y[,1])
#match(X.pseudo.autosomal2.SNPs, SCC.Y[,1])

#None of the SNPs in the X.pseudo.autosomal.region match to the SNPs in the Y chromosome.

TCGA_parser <- function(TCGA.ID){
  
  x <-substr(as.character(TCGA.ID), 1, 16)
  return(x)
  
}

TCGA_parser2 <- function(TCGA.ID){
  
  x <-substr(as.character(TCGA.ID), 1, 12)
  return(x)
  
}


SNP_file_parser <- function(file){
  
  x <- strsplit(file, "\\.")[[1]][1]
  return(x)
  
}


#setwd("D:/Phil/betareg/data")
setwd("/home/delores/Academic/LRRI2013/Phil/betareg/data")

ADC.SNP.anno <- read.delim("broad.mit.edu_LUAD.Genome_Wide_SNP_6.sdrf.txt")
ADC.SNP.TCGA.IDs.and.files <- ADC.SNP.anno[, c(1,7)]
SCC.SNP.anno <- read.delim("broad.mit.edu_LUSC.Genome_Wide_SNP_6.sdrf.txt")
SCC.SNP.TCGA.IDs.and.files <- SCC.SNP.anno[, c(2,9)]
#No sex data in either file.


ADC.SNP.TCGA.IDs.and.files[,1] <- sapply(ADC.SNP.TCGA.IDs.and.files[,1], TCGA_parser)
SCC.SNP.TCGA.IDs.and.files[,1] <- sapply(SCC.SNP.TCGA.IDs.and.files[,1], TCGA_parser)

ADC.SNP.filenames <- list.files(ADC.SNP.directory)
ADC.SNP.filenames <- sapply(ADC.SNP.filenames, SNP_file_parser, USE.NAMES=FALSE)

SCC.SNP.filenames <- list.files(SCC.SNP.directory)
SCC.SNP.filenames <- sapply(SCC.SNP.filenames, SNP_file_parser, USE.NAMES=FALSE)


#Read in SNP file to TCGA ID data
#Match to SNP.files

ADC.TCGA.IDs <- ADC.SNP.TCGA.IDs.and.files[match(ADC.SNP.filenames, ADC.SNP.TCGA.IDs.and.files[,2]), 1]
SCC.TCGA.IDs <- SCC.SNP.TCGA.IDs.and.files[match(SCC.SNP.filenames, SCC.SNP.TCGA.IDs.and.files[,2]), 1]

ADC.TCGA.patient.IDs <- sapply(ADC.TCGA.IDs, TCGA_parser2)
SCC.TCGA.patient.IDs <- sapply(SCC.TCGA.IDs, TCGA_parser2)

rownames(ADC.Y) <- ADC.Y[,1]
ADC.Y <- ADC.Y[,-1]
rownames(SCC.Y) <- SCC.Y[,1]
SCC.Y <- SCC.Y[,-1]

names(ADC.Y) <- c(ADC.TCGA.IDs)
names(SCC.Y) <- c(SCC.TCGA.IDs)

ADC.Y <- t(ADC.Y)
SCC.Y <- t(SCC.Y)

ADC.gender <- rep(NA, nrow(ADC.Y))
SCC.gender <- rep(NA, nrow(SCC.Y))
 
ADC.Y <- as.data.frame(cbind(ADC.gender, ADC.Y))
SCC.Y <- as.data.frame(cbind(SCC.gender, SCC.Y))

#Read in clinical data
setwd("/home/delores/Academic/LRRI2013/Phil/find_sex/data/")
ADC.clinical.data <- read.csv("clinical_patient_public_LUAD.csv")
ADC.clinical.data <- ADC.clinical.data[, c(1,2,19,20,38,39)]

SCC.clinical.data <- read.csv("clinical_patient_all_LUSC.csv")
SCC.clinical.data <- SCC.clinical.data[, c(1,2,22,24,54,55)]


#Match sex and ethnicity clinical data to SNP data, also SNP file
#See if pattern holds.

ADC.Y[,1] <- as.character(ADC.clinical.data$gender[match(row.names(ADC.Y), ADC.clinical.data[,1])])
SCC.Y[,1] <- as.character(SCC.clinical.data$gender[match(row.names(SCC.Y), SCC.clinical.data[,1])])

#Most samples have gender data already!


#Merge clinical data and SNP annotation.
ADC.all.patient.IDs <- sapply(ADC.SNP.TCGA.IDs.and.files[,1], TCGA_parser2,USE.NAMES=FALSE)
ADC.new.anno <- cbind(ADC.all.patient.IDs, ADC.SNP.TCGA.IDs.and.files)
names(ADC.clinical.data)[1] <- "patient.ID"
names(ADC.new.anno) <- c("patient.ID", "sample.ID", "SNP.file.name")
ADC.new.anno <-merge(ADC.new.anno, ADC.clinical.data)

SCC.all.patient.IDs <- sapply(SCC.SNP.TCGA.IDs.and.files[,1], TCGA_parser2,USE.NAMES=FALSE)
SCC.new.anno <- cbind(SCC.all.patient.IDs, SCC.SNP.TCGA.IDs.and.files)
names(SCC.clinical.data)[1] <- "patient.ID"
names(SCC.new.anno) <- c("patient.ID", "sample.ID", "SNP.file.name")
SCC.new.anno <-merge(SCC.new.anno, SCC.clinical.data)


setwd( "/home/delores/Academic/LRRI2013/Phil/find_sex/output")
write.csv(ADC.new.anno, "ADC_anno_clinical.csv", row.names=FALSE)
write.csv(SCC.new.anno, "SCC_anno_clinical.csv", row.names=FALSE)


                     