

count_results <- function(hist){
  
  directory <- paste0("/home/delores/Academic/LRRI2013/Phil/betareg/output/", hist, "_reg_results/")
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
  
  directory <- paste0("/home/delores/Academic/LRRI2013/Phil/betareg/output/", hist, "_reg_results/")
  setwd(directory)
  filenames <- list.files()  
  combined.results <- as.data.frame(matrix(0, sum(count), 8))
  names(combined.results) <- c("Gene", "Probe", "SNP", "Estimate", "P.Value", 
                                "Prevalence", "Test.Allele", "ADC.or.SCC")
  ind1 <- 1
  ind2 <- count[1]
  
  for (i in 1:length(filenames)){
    
    current.file <- read.csv(filenames[i], header=TRUE)
    combined.results[ind1:ind2, ] <- current.file
    ind1 <- ind2 + 1
    ind2 <- ind1 + count[i+1] - 1
    
  }
  
  return(combined.results)
}

ADC.count <- count_results("ADC") 



SCC.count <- count_results("SCC") 
ADC.combined.results <- combine_results("ADC", ADC.count)
SCC.combined.results <- combine_results("SCC", SCC.count)
all.results <- rbind(ADC.combined.results, SCC.combined.results)
write.csv(all.results, /home/delores/Academic/LRRI2013/Phil/betareg/output/total_reg_results_unadjusted.csv", row.names=FALSE)
