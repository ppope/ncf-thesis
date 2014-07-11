

setwd("D:/Phil/adjusted_betareg/data")

#Add new ethnicity data to covariates.

ethnicity.anno <- read.csv("ethnicity_anno.csv", stringsAsFactors=FALSE)
ADC.covariates <- read.csv("ADC_SNP_clinical_data.csv", stringsAsFactors=FALSE)
SCC.covariates <- read.csv("SCC_SNP_clinical_data.csv", stringsAsFactors=FALSE)
names(ethnicity.anno)[4] <- c("sample.ID")

#remove redundant TCGA ID column.
ADC.covariates <- ADC.covariates[-which(names(ADC.covariates) == "TCGA.ID")]
SCC.covariates <- SCC.covariates[-which(names(SCC.covariates) == "TCGA.ID")]

merge_ethnicity <- function(ethnicity, covariates){
  
  #remove existing ethnicity label
  
  covariates <- covariates[ , -which(names(covariates) == "ethnicity")]
  covariates <- merge(ethnicity, covariates)
  return(covariates)
   
}


#Recode tumor stage as
#Stage IA, Stage IB: 0 (Early Stage)
#Stage IIA, Stage IIB, Stage IIIA, Stage IIIB, Stage IV: 1 (Late Stage) 

#Recode tobacco history as
#"Current smoker":   0
#"Current reformed smoker for < or = 15 years": 1 
#"Current reformed smoker for > 15 years" : 2
#"Lifelong Non-smoker":  3

#Recode Gender as:
#Female: 0
#Male: 1

recode_covariates <- function(covariates){
  
  N <- nrow(covariates)
  
  for (i in 1:N){
    
    gender <- covariates[i,8]
    history <- covariates[i,9]
    stage <- covariates[i,10]
    
    if (gender == "FEMALE") { gender <- 0
    } else if (gender == "MALE") gender <- 1
    
    if (history == "Current smoker") { history <- 0
    } else if (history == "Current reformed smoker for < or = 15 years")  { history <- 1
    } else if (history == "Current reformed smoker for > 15 years")  { history <- 2  
    } else if (history == "Lifelong Non-smoker" ) history <- 3  
    
    if ((stage == "Stage IA") | (stage == "Stage IB") | (stage == "Stage I")) { stage <- 0
    } else if ((stage == "Stage IIA") | (stage == "Stage IIA") | (stage == "Stage IIB") | (stage == "Stage IIIA") | (stage == "Stage IIIB") | (stage == "Stage IV")) {
      stage <- 1
    }
    
   covariates[i,c(8,9,10)] <- c(gender, history, stage)
    
    
  }
  
  return(covariates)
  
}


remove_nulls <- function(covariates){
  
  
  x <- which(covariates$tumor_stage == "null")
  y <- which(covariates$gender == "null")
  z <- which(covariates$age == "null")
  u <- which(covariates$tobacco_smoking_history_indicator == "null")
  
  null.ind <- unique(c(x,y,z,u))
  
  covariates <- covariates[-null.ind,]
  
  return(covariates)
  
}



#This chunk merges ethnicity results with ADC.covariates.
ADC.covariates <- merge_ethnicity(ethnicity.anno, ADC.covariates)
SCC.covariates <- merge_ethnicity(ethnicity.anno, SCC.covariates)


#This chunk recodes the covariates.
ADC.covariates.recoded <- recode_covariates(ADC.covariates)
SCC.covariates.recoded <- recode_covariates(SCC.covariates)

#This chunk reorders columns for convenience.
ADC.covariates.recoded[1:12] <- ADC.covariates.recoded[c(1,6,5,11,12,2,3,4,7,8,9,10)]
names(ADC.covariates.recoded)[1:12] <- names(ADC.covariates.recoded)[c(1,6,5,11,12,2,3,4,7,8,9,10)]
SCC.covariates.recoded[1:12] <- SCC.covariates.recoded[c(1,6,5,11,12,2,3,4,7,8,9,10)]
names(SCC.covariates.recoded)[1:12] <- names(SCC.covariates.recoded)[c(1,6,5,11,12,2,3,4,7,8,9,10)]


#This chunk removes null values
ADC.covariates.recoded <- remove_nulls(ADC.covariates.recoded)
SCC.covariates.recoded <- remove_nulls(SCC.covariates.recoded)


write.csv(ADC.covariates.recoded, "ADC_total_covariates.csv", row.names=FALSE)
write.csv(SCC.covariates.recoded, "SCC_total_covariates.csv", row.names=FALSE)
