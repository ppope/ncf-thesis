


###given symbol : symbol(s) in file
###ACP5: ACP5
###C2orf32: CRIP1
###CCDC50: CCDC50 CCDC50;UTS2D UTS2D;CCDC50
###CPEB1: CPEB1
###CST6: CST6
###DKK3: DKK3
###GJB2: GJB2
###GSTT2: GSTT2B;GSTT2 GSTT2:GSTT2B DTT;GSTT2 GSTT2;DTT
###HTRA3: HTRA3 
###HTRA4: HTRA4 HTRA4;PLEKHA2 PLEKHA2;HTRA4
###IGDCC4 (NOPE): IGDCC4
###LAYN: LAYN
###LOX: LOX
###MARVELD3: MARVELD3
###MSX1: MSX1
###MT1M: MT1M
###PAK6-2: PAK6;C15orf56 C15orf56;PAK6 
###(Note: the gene is PAK6 but the methylated CpG island is the one in the middle of the gene that is close to the TSS for the smaller transcripts. Only the probes corresponding to this spot were considered, with corresponding gene symbol PAK6;C15orf56. See http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&position=chr15:40509629-40569688&hgsid=338380169&knownGene=full)
###RBP1: RBP1
###ROBO3: ROBO3
###TMEM22: TMEM22
###TPM1: TPM1
###ZNF365: ZNF365
###BNIP3: BNIP3
###EMILIN2: EMILIN2
###MSRB3: MSRB3
###VAV3: VAV3
###AKAP12: AKAP12
###SYNE1: SYNE1
###IL20R: IL20RA
###p16: CDKN2A CDKN2A;CDKN2BAS CDKN2BAS;CDKN2A
###RASSF1A: RASSF1 RASSF1;ZMYND10 ZMYND10;RASSF1
###MGMT-O6: MGMT
###SFRP2: SFRP2
###SFRP1: SFRP1
###BETA3, BHLHB5: BHLHE22, BHLHE22;LOC401463 LOC401463;BHLHE22
###NOVEL2: TOX2
###IGFBP3: IGFBP3
###Reprimo (RPRM): RPRM
###DCR1: TNFRSF10C
###HLHP: BHLHE41
###AK5: AK5
###TSLC1: CADM1
###TUBB4: TUBB4 MIR220B;TUBB4 TUBB4;MIR220B
###XT3: SLC6A20
###PCDH10: PCDH10
###MMP2: MMP2
###APC2: APC2
###GPNMB: GPNMB
###BOLL: BOLL
###ICAM5: ICAM5 ICAM5;ICAM4
###JPH3: JPH3
###GNB4: GNB4
###SYNE1: SYNE1


setwd("Z:/Leng/LUAD/JHU_USC__HumanMethylation450/Level_3")
filenames=list.files()
file1=read.table(filenames[1],sep="\t",header=T)
dimnames(file1)[[2]]
x=table(file1$gene.symbol)


idx1=which(file1$gene.symbol=="ACP5")
idx2=which(file1$gene.symbol=="CRIP1")
idx3=which(file1$gene.symbol=="CCDC50" | file1$gene.symbol=="CCDC50;UTS2D" | file1$gene.symbol=="UTS2D;CCDC50")
idx4=which(file1$gene.symbol=="CPEB1")
idx5=which(file1$gene.symbol=="CST6")
idx6=which(file1$gene.symbol=="DKK3")
idx7=which(file1$gene.symbol=="GJB2")
idx8=which(file1$gene.symbol=="GSTT2B;GSTT2" | file1$gene.symbol=="GSTT2:GSTT2B")
idx9=which(file1$gene.symbol=="HTRA3")
idx10=which(file1$gene.symbol=="HTRA4" | file1$gene.symbol=="HTRA4;PLEKHA2" | file1$gene.symbol=="PLEKHA2;HTRA4")
idx11=which(file1$gene.symbol=="IGDCC4")
idx12=which(file1$gene.symbol=="LAYN")
idx13=which(file1$gene.symbol=="LOX")
idx14=which(file1$gene.symbol=="MARVELD3")
idx15=which(file1$gene.symbol=="MSX1")
idx16=which(file1$gene.symbol=="MT1M")
idx17=which(file1$gene.symbol=="PAK6;C15orf56" | file1$gene.symbol=="C15orf56;PAK6" )
idx18=which(file1$gene.symbol=="RBP1")
idx19=which(file1$gene.symbol=="ROBO3")
idx20=which(file1$gene.symbol=="TMEM22")
idx21=which(file1$gene.symbol=="TPM1")
idx22=which(file1$gene.symbol=="ZNF365")
idx23=which(file1$gene.symbol=="BNIP3")
idx24=which(file1$gene.symbol=="EMILIN2")
idx25=which(file1$gene.symbol=="MSRB3")
idx26=which(file1$gene.symbol=="VAV3")
idx27=which(file1$gene.symbol=="AKAP12")
idx28=which(file1$gene.symbol=="SYNE1")
idx29=which(file1$gene.symbol=="IL20RA")
idx30=which(file1$gene.symbol=="CDKN2A" | file1$gene.symbol=="CDKN2A;CDKN2BAS" | file1$gene.symbol=="CDKN2BAS;CDKN2A")
idx31=which(file1$gene.symbol=="RASSF1" | file1$gene.symbol=="RASSF1;ZYMND10" | file1$gene.symbol=="ZYMND10;RASSF1")
idx32=which(file1$gene.symbol=="MGMT")
idx33=which(file1$gene.symbol=="SFRP2")
idx34=which(file1$gene.symbol=="SFRP1")
idx35=which(file1$gene.symbol=="BHLHE22" | file1$gene.symbol=="BHLHE22;LOC401463" | file1$gene.symbol=="LOC401463;BHLHE22")
idx36=which(file1$gene.symbol=="TOX2")
idx37=which(file1$gene.symbol=="IGFBP3")
idx38=which(file1$gene.symbol=="RPRM")
idx39=which(file1$gene.symbol=="TNFRSF10C")
idx40=which(file1$gene.symbol=="BHLHE41")
idx41=which(file1$gene.symbol=="AK5")
idx42=which(file1$gene.symbol=="CADM1")
idx43=which(file1$gene.symbol=="TUBB4" | file1$gene.symbol=="MIR220B;TUBB4" | file1$gene.symbol=="TUBB4;MIR220B")
idx44=which(file1$gene.symbol=="SLC6A20")
idx45=which(file1$gene.symbol=="PCDH10")
idx46=which(file1$gene.symbol=="MMP2")
idx47=which(file1$gene.symbol=="APC2")
idx48=which(file1$gene.symbol=="GPNMB")
idx49=which(file1$gene.symbol=="BOLL")
idx50=which(file1$gene.symbol=="ICAM5" | file1$gene.symbol=="ICAM5;ICAM4" | file1$gene.symbol=="ICAM4;ICAM5")
idx51=which(file1$gene.symbol=="JPH3")
idx52=which(file1$gene.symbol=="GNB4")
idx53=which(file1$gene.symbol=="SYNE1")


idx.all <- c(idx1, idx2, idx3, idx4, idx5, idx6, idx7, idx8, idx9,
             idx10, idx11, idx12, idx13, idx14, idx15, idx16, idx17, idx18, idx19,
             idx20, idx21, idx22, idx23, idx24, idx25, idx26, idx27, idx28, idx29,
             idx30, idx31, idx32, idx33, idx34, idx35, idx36, idx37, idx38, idx39,
             idx40, idx41, idx42, idx43, idx44, idx45, idx46, idx47, idx48, idx49, 
             idx50, idx51, idx52, idx53)
  
Extract1=t(file1[idx.all,c(2,4,3)])
dimnames(Extract1)[[1]]<-c("probe.name","gene.symbol",as.character(filenames[1]))
data=NULL
data=rbind(data,Extract1) 

#added timer to determine run time of for loop
ptm <- proc.time()

for (i in c(2:length(filenames))){
 file=read.table(filenames[i],sep="\t",header=T)
 trans.file=t(file[idx.all,3])
 dimnames(trans.file)[[1]]<-list(as.character(filenames[i]))
 data=rbind(data,trans.file)
}

#end timer
proc.time() - ptm

write.csv(data,"Z:/Leng/Phil/Methylation_SCC_vs_Adeno_2/output/Methylation450K_TCGA_SCCvsAdeno_2_Adeno_1.csv")

#user  system elapsed 
#4326.71   54.62 5239.17 

