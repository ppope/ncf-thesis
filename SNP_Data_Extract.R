

setwd("Z:/Leng/LUAD/broad.mit.edu_LUAD.Genome_Wide_SNP_6.Level_2")
#read in SNP data, check if they are birdseed (SNP genotyping algorithm)
filenames <- list.files()
idx <- NULL
for (i in c(1:length(filenames))){
  if ("birdseed" == strsplit(filenames[i],"\\.")[[1]][2]) idx <- c(idx,filenames[i])
}

N <- length(filenames)
file1 <- read.table(idx[1],sep="\t",header=F,skip=2)
M <- nrow(file1)
#Preallocate data matrix
data <- data.frame(matrix(ncol = M, nrow = N))


#Columns of files are: Composite Element REF, Call, Confidence
#SNP IDs are same for each file. 
#Extract out SNP IDs and first file's data. Data from remaining files will be added.
trans.file1 <- t(file1[,c(1,2)])
data[1:2,] <- trans.file1

#Make one big 'ole data file.
for (i in 2:N){
  file <- read.table(idx[i],sep="\t",header=F,skip=2)
  trans.file <- t(file[,2])
  data[i+1,] <- trans.file
}

# this takes a real long time.
# write.csv(data,"D:/Phil/SNPs/output/SNP_data.csv", row.names=c("Name",idx))


# extract SNPs listed in genes
# want a vector containing all SNPs for the genes, to cut down the data file
# and list of vectors containing SNPs for each gene.
# setwd("/home/delores/Academic/LRRI2013/Phil/SNPs/output/SNP_Info/")

setwd("D:/Phil/SNPs/output/SNP_Info/")

info.filenames <- list.files()
L <- length(info.filenames)
gene.SNPs.list <- vector(mode = "list", length = L)
all.SNPs <- c()

for (i in 1:L){
  filename <- info.filenames[i]
  gene <- strsplit(filename, "_")[[1]][1]
  gene.file <- read.csv(filename, header=TRUE)
  gene.SNPs.list[[i]] <- as.character(gene.file[,1])
  names(gene.SNPs.list)[[i]] <- gene
  all.SNPs <- c(all.SNPs, as.character(gene.file[,1]))
}

all.SNPs <- unique(all.SNPs)


#cut down data
SNPs.ind <- match(all.SNPs, data[1,])
#Not all SNPs in gene regions have corresponding entry in SNP-array files
SNPs.ind2 <- SNPs.ind[-which(is.na(SNPs.ind) == TRUE)]
data <- data[, SNPs.ind2]
names(data) <- data[1,]
data <- data[-1,]

Name <- idx

data <-  cbind(Name,data)


#make a separate file of SNP data for each gene.

setwd("D:/Phil/SNPs/output/ADC_SNP_data/")
O <- length(gene.SNPs.list)

for (i in 1:O){
  
  gene <- names(gene.SNPs.list)[i]
  match.idx <- match(gene.SNPs.list[[i]], names(data))
  match.idx <- match.idx[-(which(is.na(match.idx) == TRUE))]
  filename <- paste(gene, "_ADC_SNP_data.csv", sep="")
  write.csv(data[,c(1,match.idx)], filename, row.names=FALSE)
  
}


