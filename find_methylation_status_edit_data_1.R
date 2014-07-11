
setwd("D:/Leng/Phil/find_methylation_status/data")
filenames <- c("Methylation450K_TCGA_SCCvsAdeno_2_total_ordered.csv", "probenames_TCGA_SCCvsAdeno_2_Adeno_1.csv")


gene.symbol.list <- c( "ACP5"       ,       "AK5"        ,       "AKAP12"      ,      "APC2"         ,     "BHLHE22"   ,       
                       "BHLHE22;LOC401463", "BHLHE41"      ,     "BNIP3"         ,    "BOLL"          ,    "CADM1"        ,    
                       "CCDC50"       ,     "CCDC50;UTS2D"   ,   "CDKN2A"     ,       "CDKN2A;CDKN2BAS"  , "CDKN2BAS;CDKN2A" , 
                       "CPEB1"      ,       "CRIP1"     ,        "CST6"      ,        "DKK3"      ,        "EMILIN2"     ,     
                       "GJB2"         ,     "GNB4"         ,     "GPNMB"        ,     "GSTT2B;GSTT2" ,     "HTRA3"    ,        
                       "HTRA4"         ,    "HTRA4;PLEKHA2" ,    "ICAM5"  ,           "ICAM5;ICAM4"   ,    "IGDCC4"  ,         
                       "IGFBP3"     ,       "IL20RA"    ,        "JPH3"         ,     "LAYN"      ,        "LOX"      ,        
                       "MARVELD3"      ,    "MGMT"        ,      "MIR220B;TUBB4"   ,  "MMP2"        ,      "MSRB3"   ,         
                       "MSX1"        ,      "MT1M"     ,         "PAK6;C15orf56" ,    "PCDH10"     ,       "RASSF1"    ,       
                       "RBP1"        ,      "ROBO3"     ,        "RPRM"     ,         "SFRP1"      ,       "SFRP2"     ,       
                       "SLC6A20"      ,     "SYNE1"     ,        "TMEM22"      ,      "TNFRSF10C"    ,     "TOX2"     ,        
                       "TPM1"         ,     "TUBB4"       ,      "TUBB4;MIR220B"   ,  "UTS2D;CCDC50"   ,   "VAV3"   ,          
                       "ZNF365")


data <- read.csv(filenames[1],sep=",",header=TRUE)

probenames <- read.csv(filenames[2],sep=",",header=TRUE)


#Need to attach corresponding gene symbols to probes, so data can be separated into files by gene symbols.

#put gene names in convenient form for attaching.
t.probenames1 <- t(probenames)
t.probenames2 <- t.probenames1[2,]
names(t.probenames2) <- t.probenames1[1,]



probe.data <- data[,-(1:4)]
sample.data <- data[,1:4]

#matches order of second argument to first.
order <- match(names(probe.data),names(t.probenames2))

names(probe.data)== names(t.probenames2)[order]

x <- rbind(t.probenames2[order],probe.data)


#want to split on the levels of gene.symbol


f <- as.factor(as.character(x[1,]))

N <- length(levels(f))


#make a list of indices in data corresponding to which gene symbol is being extracted (in the levels of gene.symbol.list)
ind <- vector("list",N)
j <- 1
for (i in levels(f)){
  ind[[j]] <- which(x[1,] == i)
  j <- j + 1
}

for (i in 1:N){
  
  probe.data.subset <- probe.data[,ind[[i]]]
  data.subset <- cbind(sample.data,probe.data.subset)
  
  #split again on adeno vs. SCC
  adeno.data.subset <- data.subset[which(data$Adeno.or.SCC == "Adeno"),]
  SCC.data.subset <- data.subset[which(data$Adeno.or.SCC == "SCC"),]
  
  adeno.filename <- paste("Z:/Leng/Phil/find_methylation_status/output/genes/find_methylation_status_adeno_",as.character(levels(f)[i]),".csv", sep="")
  SCC.filename <- paste("Z:/Leng/Phil/find_methylation_status/output/genes/find_methylation_status_SCC_",as.character(levels(f)[i]),".csv", sep="")
  
  write.csv(adeno.data.subset, file = adeno.filename, row.names=FALSE)
  write.csv(SCC.data.subset, file = SCC.filename, row.names=FALSE)
  
}

# Didn't end up needed the following:
# x <- data.frame("gene.symbols" = gene.symbol.list, "adeno.tumor.methylated" = rep(" ", 61), "SCC.tumor.methylated" = rep(" ", 61))
# write.csv(x, "Z:/Leng/Phil/find_methylation_status/output/RESULTS_find_methylation_status.csv", row.names=FALSE)

