
setwd("Z:/Leng/Phil/Methylation_SCC_vs_Adeno_1/output")
filenames <- c("Methylation450K_TCGA_SCCvsAdeno_Adeno_2.csv","Methylation450K_TCGA_SCCvsAdeno_SCC_2.csv")
histnames <- c("clinical_sample_all_LUAD_1.csv","clinical_sample_all_LUSC_1.csv")




#read in data
data1 <- read.csv(filenames[1],sep=",",header=TRUE)
#read in histology
hist1 <- read.csv(histnames[1],sep=",",header=TRUE)
names(hist1)[1] <- "sample.name"
data.and.hist1 <- merge(data1,hist1)


#read in data
data2 <- read.csv(filenames[2],sep=",",header=TRUE)
#read in histology
hist2 <- read.csv(histnames[2],sep=",",header=TRUE)
names(hist2)[1] <- "sample.name"

#why were there 9 files lost here?
#some patients have multiple samples
data.and.hist2 <- merge(data2,hist2)

x1 <- sort(names(data.and.hist1), decreasing=TRUE, index.return=TRUE)
permute1 <- x1$ix
data.and.hist1.ordered <- data.and.hist1[,permute1]

x2 <- sort(names(data.and.hist2), decreasing=TRUE, index.return=TRUE)
permute2 <- x2$ix
data.and.hist2.ordered <- data.and.hist2[,permute2]

#check if names match
names(data.and.hist1.ordered) == names(data.and.hist2.ordered)

#merge Adeno with SCC
merged.data <- rbind(data.and.hist1.ordered, data.and.hist2.ordered)

#create vector designating adeno or SCC
Adeno.or.SCC <- c(rep("Adeno",254), rep("SCC",191))

total.data <- cbind(Adeno.or.SCC, merged.data)

write.csv(total.data,"Z:/Leng/Phil/Methylation_SCC_vs_Adeno_1/output/Methylation450K_TCGA_SCCvsAdeno_total.csv", row.names=FALSE)

