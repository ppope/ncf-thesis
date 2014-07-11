setwd("Z:/Leng/LUAD/IL6")
il6.file=read.csv("IL6_SNP_Info.csv",header=T)
il6.snp=as.character(il6.file[,2])
setwd("Z:/Leng/LUAD/broad.mit.edu_LUAD.Genome_Wide_SNP_6.Level_2")
filenames=list.files()
idx=NULL
for (i in c(1:length(filenames))){
if ("birdseed" == strsplit(filenames[i],"\\.")[[1]][2]) idx=c(idx,filenames[i])
           } 
data=NULL
file1=read.table(idx[1],sep="\t",header=F,skip=2)
trans.file1=t(file1[,c(1,2)])
data=rbind(data,trans.file1)
snp.names=data[1,]
il6.idx=NULL
for (i in c(1:length(il6.snp))){
  match.idx=match(il6.snp[i],snp.names)
  il6.idx=cbind(il6.idx, match.idx)
}

data.il6=data[,il6.idx]

for (i in c(2:length(idx))){
 file=read.table(idx[i],sep="\t",header=F,skip=2)
 trans.file=t(file[,2])
 data.il6=rbind(data.il6,trans.file[,il6.idx])
}
setwd("Z:/Leng/LUAD/IL6")
write.csv(data.il6,"data.il6.csv",row.names=c("Name",idx))


