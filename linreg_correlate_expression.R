

setwd("D:/Phil/correlate_expression/data")

ADC.ex.data <- read.csv("ADC_expression_data.csv", stringsAsFactors=FALSE) 
SCC.ex.data <- read.csv("SCC_expression_data.csv", stringsAsFactors=FALSE) 
SCC.ex.data <- SCC.ex.data[-which(duplicated(SCC.ex.data$patient.ID) == TRUE), ]

ADC.meth.data <- read.csv("D:/Phil/adjusted_betareg/data/ADC_total_meth_data.csv", stringsAsFactors=FALSE)
SCC.meth.data <- read.csv("D:/Phil/adjusted_betareg/data/SCC_total_meth_data.csv", stringsAsFactors=FALSE)
SCC.meth.data <- SCC.meth.data[-anyDuplicated(SCC.meth.data$patient.name), ]

results.table <- read.csv("D:/Phil/adjusted_betareg/data/sig_reg_results.csv", stringsAsFactors=FALSE)
results.table <- results.table[which(results.table$ADC.or.SCC == "ADC"), c(1,2,3,6,7,8,11)]
results.table$Estimate <- 0
results.table$P.Value <- 0
results.table$ADC.or.SCC <- 0
results.table <- cbind(results.table, "r.squared" = rep(0, nrow(results.table)))

ex.data <- SCC.ex.data
meth.data <- SCC.meth.data
linreg.results.table <- results.table
histol <- "SCC"



linreg <- function(ex.data, meth.data, linreg.results.table, histol){
  
  N <- nrow(linreg.results.table)
  linreg.results.table$ADC.or.SCC <- rep(histol, N)
  
  for (i in 1:N){
    
    current.probe <- linreg.results.table$Meth.Probe.Name[i]
    current.gene <- linreg.results.table$Gene[i]
    
    EX <- ex.data[match(current.gene, names(ex.data))]
    row.names(EX) <- ex.data$patient.ID
    
    BV <- meth.data[match(current.probe, names(meth.data))]
    row.names(BV) <- meth.data$patient.name
    BV <- BV[which(is.na(BV[,1]) == FALSE), 1, drop=FALSE]
  
    reg.table <- merge(BV, EX, by="row.names")
    
    if (nrow(reg.table) == 0){
      linreg.results.table$Estimate[i] <- NA
      linreg.results.table$P.Value[i] <- NA  
    } else{
      
      Beta_Value <- unlist(reg.table[2], use.names=FALSE)
      Expression_Value <- unlist(reg.table[3], use.names=FALSE)
      
      model <- lm(Expression_Value ~ Beta_Value)    
      linreg.results.table$Estimate[i] <- model$coefficients[2]
      linreg.results.table$P.Value[i] <- summary(model)$coefficients[2,4]
      linreg.results.table$r.squared[i] <- summary(model)$r.squared
    }      
  }
  
  return(linreg.results.table)
  
}

ADC.linreg.results <- linreg(ADC.ex.data, ADC.meth.data, results.table, "ADC")
SCC.linreg.results <- linreg(SCC.ex.data, SCC.meth.data, results.table, "SCC")

ADC.linreg.results <- ADC.linreg.results[which(is.na(ADC.linreg.results$P.Value) == FALSE), ]
SCC.linreg.results <- SCC.linreg.results[which(is.na(SCC.linreg.results$P.Value) == FALSE), ]

x <- rbind(ADC.linreg.results, SCC.linreg.results)

setwd("D:/Phil/correlate_expression/output/")

write.csv(x, "linreg_expresssion_results.csv", row.names=FALSE)

#SCC has expression data for all sig probes.
#ADC has expression data for only one probe, cg12317456

