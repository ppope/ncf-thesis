

#Preprocess structure results.
#Attach TCGA IDs to each pop.id
#Only need two columns, the third is dependent.

#Add Ethnicity to covariates.
#Run adjusted model.

#Run unadjusted model with 27k data.


setwd("D:/Phil/structure/data")

structure.results.raw <- read.table("results.txt", sep = "\t", skip=2)
N <- nrow(structure.results.raw)
structure.results <- data.frame("PopID" = rep(0, N), "proportion.african" = rep(0, N), 
                                "proportion.native.amer"= rep(0, N), "proportion.euro" = rep(0, N))


structure.results.raw <- as.vector(structure.results.raw)

structure_results_parser <- function(raw.data){
  
  x <- as.character(raw.data)
  x <- strsplit(x, " {1,}") #breaks at one or more " ".
  pop.id <- x[[1]][2]
  prop.afr <- x[[1]][6]
  prop.euro <- x[[1]][7]
  prop.na <- x[[1]][8]
  
  return(c(pop.id, prop.afr, prop.euro, prop.na))
  
}

for (i in 1:N) structure.results[i, ] <- structure_results_parser(structure.results.raw[i,])
head(structure.results)

POP.ID.to.TCGA.ID <- read.csv("LUAD.GBM.COAD.LUSC.POPID.ethnicity.csv", stringsAsFactors=FALSE)
POP.ID.to.TCGA.ID <- POP.ID.to.TCGA.ID[ , c(1,2,3,394)]

ethnicity.anno <- merge(structure.results, POP.ID.to.TCGA.ID)

ethnicity.anno <- ethnicity.anno[,c(2,3,4,5)]
head(ethnicity.anno)
setwd("D:/Phil/structure/output")
write.csv(ethnicity.anno, "ethnicity_anno.csv", row.names=FALSE)
