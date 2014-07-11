genes <- c("ACP5", "C2orf32", "CCDC50", "CPEB1", "CST6", "DKK3", "GJB2", "GSTT2", "HTRA3", "HTRA4", "IGDCC4", "LAYN", "LOX", "MARVELD3", "MSX1",
           "MT1M", "PAK6-2", "RBP1", "ROBO3", "TMEM22", "TPM1", "ZNF365", "BNIP3", "EMILIN2", "MSRB3", "VAV3", "AKAP12", "SYNE1", "IL20R",
           "p16", "RASSF1A", "MGMT-O6", "SFRP2", "SFRP1", "BETA3", "NOVEL2", "IGFBP3", "RPRM", "DCR1", "HLHP", "AK5", "TSLC1", "TUBB4",
           "XT3", "PCDH10", "MMP2", "APC2", "GPNMB", "BOLL", "ICAM5", "JPH3", "GNB4")

symbols1 <- c("ACP5", "CRIP1", "CCDC50", "CPEB1", "CST6", "DKK3", "GJB2",
                 "GSTT2B;GSTT2", "HTRA3", "HTRA4", "IGDCC4", "LAYN", "LOX", "MARVELD3", "MSX1", "MT1M", "PAK6;C15orf56", "RBP1",
                 "ROBO3", "TMEM22", "TPM1", "ZNF365", "BNIP3", "EMILIN2", "MSRB3", "VAV3", "AKAP12", "SYNE1",
                 "IL20RA", "CDKN2A", "RASSF1", "MGMT",
                 "SFRP2", "SFRP1", "BHLHE22", "TOX2", "IGFBP3", "RPRM", "TNFRSF10C",
                 "BHLHE41", "AK5", "CADM1", "TUBB4", "SLC6A20", "PCDH10", "MMP2", "APC2",
                 "GPNMB", "BOLL","ICAM5", "JPH3", "GNB4")

symbols2 <- c( "", "", "CCDC50;UTS2D", "", "", "", "", "GSTT2;GSTT2B", "", "HTRA4;PLEKHA2", "", "", "", "", "", "", "C15orf56;PAK6", "", "", "", "", "", "", "", 
               "", "", "", "", "", "", "", "", "CDKN2A;CDKN2BAS", "RASSF1;ZMYND10", "", "", "", "BHLHE22;LOC401463", "", "", "", "", "", "", "MIR220B;TUBB4", "", "", "", "", "ICAM5;ICAM4", "", "")


symbols3 <- c( "", "", "UTS2D;CCDC50", "", "", "", "", "DTT;GSTT2", "", "PLEKHA2;HTRA4", "", "", "", "", "", "", "", "", "", "", "", "", "", "", 
               "", "", "", "", "", "CDKN2BAS;CDKN2A", "ZMYND10;RASSF1", "", "", "", "LOC401463;BHLHE22", "", "", "", "", "", "", "TUBB4;MIR220B", "", "", "", "", "", "", "", "", "", "")

symbols4 <- c( "", "", "", "", "", "", "", "GSTT2;DTT", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", 
               "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "")



dfrm <- data.frame(cbind(genes, symbols1, symbols2, symbols3, symbols4))

dfrm

write.csv(dfrm, "D:/Phil/SNPs/genes_to_symbols.csv", row.names=FALSE)