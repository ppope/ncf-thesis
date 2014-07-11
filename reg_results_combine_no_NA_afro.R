

count_results <- function(hist){
  
  directory <- paste0("D:/Phil/betareg/output/", hist, "_reg_results_unadjusted_no_NA_afro/")
  setwd(directory)
  filenames <- list.files()
  count <- rep(0, length=length(filenames))

  for (i in 1:length(filenames)){
    current.file <- read.csv(filenames[i], header=TRUE)
    count[i] <- nrow(current.file)
  }
  
  return(count)

}


combine_results <- function(hist, count){
  
  directory <- paste0("D:/Phil/betareg/output/", hist, "_reg_results_unadjusted_no_NA_afro/")
  setwd(directory)
  filenames <- list.files()  
  v <- vector(mode="character", length=sum(count))
  w <- vector(mode="numeric", length=sum(count))
  combined.results <- data.frame( "Gene"= v, "Probe" = v, "SNP" = v, "Estimate" = w, "P.Value" = w, 
                                  "Prevalence" = w, "Test.Allele" = v, "ADC.or.SCC" = v, "Sample.Size" = w)
  
  for (i in c(1,2,3,7,8)) combined.results[,i] <- as.character( combined.results[ ,i])

  ind1 <- 1
  ind2 <- count[1]
  
  
  
  for (i in 1:length(filenames)){
    
    current.file <- read.csv(filenames[i], header=TRUE)
    
    for (j in c(1,2,3,7,8)) current.file[ ,j] <- as.character( current.file[, j])
    
    combined.results[ind1:ind2, ] <- current.file
    
    combined.results[ind1:ind2, ] <- as.data.frame(current.file, stringsAsFactors = FALSE)
    ind1 <- ind2 + 1
    ind2 <- ind1 + count[i+1] - 1
    
  }
  
  return(combined.results)
}

ADC.count <- count_results("ADC") 
ADC.combined.results <- combine_results("ADC", ADC.count)
SCC.count <- count_results("SCC") 
SCC.combined.results <- combine_results("SCC", SCC.count)
all.results <- rbind(ADC.combined.results, SCC.combined.results)
write.csv(all.results, "D:/Phil/betareg/output/total_reg_results_unadjusted_no_NA_afro.csv", row.names=FALSE)


