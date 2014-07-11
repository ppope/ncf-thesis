###C2cd4d: C2CD4D;LOC100132111   LOC100132111;C2CD4D
###Clic6:CLIC6
###Didol:DIDO1 DIDO1;C20orf11 C20orf11;DIDO1
###Gdnf: GDNF
###Hnf1b:HNF1B
###Me3:ME3
###Msc:MSC
###Thbs2:THBS2
###IL_6:IL6
setwd("Z:/Leng/LUSC/LUSC.HumanMethylation450.Level_3")
filenames=list.files()
file1=read.table(filenames[1],sep="\t",header=T,skip=1)
dimnames(file1)[[2]]
x=table(file1$Gene_Symbol)
idx1=which(file1$Gene_Symbol=="C2CD4D;LOC100132111"|
file1$Gene_Symbol=="LOC100132111;C2CD4D")
idx2=which(file1$Gene_Symbol=="CLIC6")
idx3=which(file1$Gene_Symbol=="DIDO1;C20orf11" |
file1$Gene_Symbol=="C20orf11;DIDO1")
idx4=which(file1$Gene_Symbol=="GDNF")
idx5=which(file1$Gene_Symbol=="HNF1B")
idx6=which(file1$Gene_Symbol=="ME3")
idx7=which(file1$Gene_Symbol=="MSC")
idx8=which(file1$Gene_Symbol=="THBS2")
idx9=which(file1$Gene_Symbol=="IL6")
Extract1=t(file1[c(idx1,idx2,idx3,idx4,idx5,idx6,idx7,idx8,idx9),c(1,3,2)])
#6/5/2013: changed below from dimnames to row.names for my own ease of understanding (Phil) 
row.names(Extract1)<-c("probe.name","gene.symbol",as.character(filenames[1]))
data=NULL
data=rbind(data,Extract1) 

#added timer to determine run time of for loop
ptm <- proc.time()

for (i in c(2:length(filenames))){
 file=read.table(filenames[i],sep="\t",header=T,skip=1)
 trans.file=t(file[c(idx1,idx2,idx3,idx4,idx5,idx6,idx7,idx8,idx9),2])
 row.names(trans.file) <- as.character(filenames[i])
 data=rbind(data,trans.file)
}

#end timer
proc.time() - ptm

write.csv(data,"Z:/Leng/Phil/Methylation_SCC_vs_Adeno_1/data/Methylation450K_TCGA_SCCvsAdeno_SCC_1.csv")

#user  system elapsed 
#1671.60   29.92 2133.56
