###C2cd4d: C2CD4D;LOC100132111   LOC100132111;C2CD4D
###Clic6:CLIC6
###Didol:DIDO1 DIDO1;C20orf11 C20orf11;DIDO1
###Gdnf: GDNF
###Hnf1b:HNF1B
###Me3:ME3
###Msc:MSC
###Thbs2:THBS2
###IL_6:IL6
setwd("Z:/Leng/LUAD/JHU_USC__HumanMethylation450/Level_3")
filenames=list.files()
file1=read.table(filenames[1],sep="\t",header=T)
dimnames(file1)[[2]]
x=table(file1$gene.symbol)
idx1=which(file1$gene.symbol=="C2CD4D;LOC100132111"|
file1$gene.symbol=="LOC100132111;C2CD4D")
idx2=which(file1$gene.symbol=="CLIC6")
idx3=which(file1$gene.symbol=="DIDO1;C20orf11" |
file1$gene.symbol=="C20orf11;DIDO1")
idx4=which(file1$gene.symbol=="GDNF")
idx5=which(file1$gene.symbol=="HNF1B")
idx6=which(file1$gene.symbol=="ME3")
idx7=which(file1$gene.symbol=="MSC")
idx8=which(file1$gene.symbol=="THBS2")
idx9=which(file1$gene.symbol=="IL6")
Extract1=t(file1[c(idx1,idx2,idx3,idx4,idx5,idx6,idx7,idx8,idx9),c(2,4,3)])
dimnames(Extract1)[[1]]<-c("probe.name","gene.symbol",as.character(filenames[1]))
data=NULL
data=rbind(data,Extract1) 
for (i in c(2:length(filenames))){
 file=read.table(filenames[i],sep="\t",header=T)
 trans.file=t(file[c(idx1,idx2,idx3,idx4,idx5,idx6,idx7,idx8,idx9),3])
 dimnames(trans.file)[[1]]<-list(as.character(filenames[i]))
 data=rbind(data,trans.file)
}
write.csv(data,"Methylation450K_TCGA_SCCvsAdeno_Adeno.csv")


