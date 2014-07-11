
#Want to visualize the data.
#Heatmap is good choice.


#What data to use?
#Data corresponding to 52 genes of interest
#Data corresponding the methylated probes

#How to order the data?
#By chromosome location
#By genes
#By similarity

#Which distance metric to use?
#Euclidean
#Pearson

#How to add in color bar?
#Distinguish ADC, SCC, and STN.


library(gplots)

#Read in data
#Methylation450K_TCGA_SCCvsAdeno_2_total_ordered.csv contains all methylation data for the 52 probes
setwd("/home/delores/Academic/LRRI2013/Phil/find_methylation_status/data")
meth.data <- read.csv("Methylation450K_TCGA_SCCvsAdeno_2_total_ordered.csv" , stringsAsFactors=FALSE)
meth.data.numeric <- as.matrix(meth.data[-c(1,2,3,4)])
#remove NA columns
meth.data.numeric.no.na <- meth.data.numeric[, -which(is.na(meth.data.numeric[1,]) == TRUE)]
#create label for histology
histology <- ifelse(meth.data$Adeno.or.SCC == "Adeno", "ADC", "SCC")
for (i in 1:length(histology)) { if (meth.data$sample.type[i] == "Solid Tissue Normal") histology[i] <- "STN"}
#assign a color to each histology
color.map <- rep("", length(histology))
for (i in 1:length(histology)) { 
  
  if (histology[i] == "STN") color.map[i] <- "#00FF00" #STN is colored green.
  if (histology[i] == "ADC") color.map[i] <- "#FF0000" #ADC is colored red.
  if (histology[i] == "SCC") color.map[i] <- "#0000FF" #SCC is colored blue.

}
heatmap.2(meth.data.numeric.no.na, RowSideColors=color.map, trace="none", scale="none", xlab="Probes", ylab="Samples", labRow="", labCol="", hclustfun = hclust)


#Below is an unsuccessful attempt to cut the dendrogram height. 
clusf <- function(c) {
  
  hcl <- hclust(c)
  hcl <- cut(as.dendrogram(hcl), h=8)
  return(hcl$upper)
  
}
dendro <- cut(as.dendrogram(hclust(dist(meth.data.numeric.no.na))), h=8)[[1]]
heatmap.2(meth.data.numeric.no.na, Rowv=dendro, RowSideColors=color.map, trace="none", scale="none", xlab="Probes", ylab="Samples")


#This chunk creates a heatmap for the methylated probes only.
setwd("/home/delores/Academic/LRRI2013/Phil/find_methylation_status/output")
ADC.meth.data <- read.csv("find_methylation_status_ADC_methylated_data.csv" , stringsAsFactors=FALSE)
SCC.meth.data <- read.csv("find_methylation_status_SCC_methylated_data.csv" , stringsAsFactors=FALSE)
STN.meth.data <- read.csv("find_methylation_status_STN_data.csv" , stringsAsFactors=FALSE)
intersection <- intersect(names(ADC.meth.data), names(SCC.meth.data))
x <- rbind(ADC.meth.data[match(intersection, names(ADC.meth.data))], SCC.meth.data[match(intersection, names(SCC.meth.data))])
meth.data <- rbind(x, STN.meth.data[match(intersection, names(STN.meth.data))])
meth.data.numeric <- as.matrix(meth.data[-c(1,2,3,4)])
histology <- ifelse(meth.data$Adeno.or.SCC == "Adeno", "ADC", "SCC") #create label for histology
for (i in 1:length(histology)) { if (meth.data$sample.type[i] == "Solid Tissue Normal") histology[i] <- "STN"}
color.map <- rep("", length(histology))
for (i in 1:length(histology)) { 
  
  if (histology[i] == "STN") color.map[i] <- "#00FF00" #STN is colored green.
  if (histology[i] == "ADC") color.map[i] <- "#FF0000" #ADC is colored red.
  if (histology[i] == "SCC") color.map[i] <- "#0000FF" #SCC is colored blue.
  
}
heatmap.2(meth.data.numeric, RowSideColors=color.map, trace="none", scale="none", xlab="Probes", ylab="Samples", labRow="", labCol="", hclustfun = hclust)

#This chunk creates a histogram of all the 
