

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


setwd("Z:/Leng/LUSC/LUSC.Genome_Wide_SNP_6.Level_2")
#read in SNP data, check if they are birdseed (SNP genotyping algorithm)
filenames <- list.files()
idx <- NULL
for (i in c(1:length(filenames))){
  if ("birdseed" == strsplit(filenames[i],"\\.")[[1]][2]) idx <- c(idx,filenames[i])
}


N <- length(filenames)
file1 <- read.table(idx[1],sep="\t",header=F,skip=2)
trans.file1 <- t(file1[,c(1,2)])

SNPs.ind <- match(all.SNPs, trans.file1[1,])
SNPs.ind2 <- SNPs.ind[-which(is.na(SNPs.ind) == TRUE)]


#Preallocate data matrix
M <- length(SNPs.ind2)
data <- data.frame(matrix(ncol = M, nrow = N))


names(data) <- trans.file1[1,SNPs.ind2]
data[1,] <- trans.file1[2,SNPs.ind2]


#Columns of files are: Composite Element REF, Call, Confidence
#SNP IDs are same for each file. 
#Extract out SNP IDs and first file's data. Data from remaining files will be added.

#Construct data file.

ptm <- proc.time()

for (i in 2:N){
  file <- read.table(idx[i],sep="\t",header=F,skip=2)
  trans.file <- t(file[,2])
  trans.file <- trans.file[SNPs.ind2]
  data[i,] <- trans.file
}

#end timer
proc.time() - ptm


#Above loop took 15643.83 seconds (4.35 hours)


Name <- idx
data <-  cbind(Name,data)


#make a separate file of SNP data for each gene.

setwd("D:/Phil/SNPs/output/SCC_SNP_data/")
O <- length(gene.SNPs.list)

for (i in 1:O){
  
  gene <- names(gene.SNPs.list)[i]
  match.idx <- match(gene.SNPs.list[[i]], names(data))
  match.idx <- match.idx[-(which(is.na(match.idx) == TRUE))]
  filename <- paste(gene, "_SCC_SNP_data.csv", sep="")
  write.csv(data[,c(1,match.idx)], filename, row.names=FALSE)
  
}


