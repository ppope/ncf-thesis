setwd("Y:/Xiequn Zhang/AIMS/GenomeWideSNP_6.na32.annot.csv")
memory.limit(size=4000)
anote=read.csv("GenomeWideSNP_6.na32.annot.csv",header=T,skip=21)
colnames(anote)
##SNP to be extracted in chr7:22,764,766-22,767,766 NCBI Build37
anote2=subset(anote,subset=(Chromosome==7)) 
IL6= subset(anote2,subset=(as.numeric(as.character(Physical.Position))>=22764766  & 
as.numeric(as.character(Physical.Position))<=22767766) )
write.csv(IL6,"IL6_SNP_Info.csv")