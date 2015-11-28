#
# 3_Processing_Output.R
# Susan Johnston
# November 2015
#
#


#~~ Load libraries and functions

library(data.table)
library(plyr)
library(ggplot2)
library(reshape)

mod.freq <- function(vec) {sum(vec - 1)/(2*length(vec))}


#~~ Load simulation description data

RunID <- "a"
runOnEddie <- T

test.frame <- expand.grid(list(
  pop.size            = c(20, 50, 100, 150, 200, 500),
  loci.number         = c(20, 50, 100, 150, 200, 500),
  sel.thresh          = seq(0, 1, 0.2),
  modifier.found.freq = seq(0, 1, 0.2),
  generations         = 200,
  n.found.hap         = 100,
  iterations          = 100,
  modifier            = c("RNF212", "PRDM9")))

test.frame <- unique(test.frame)
test.frame$Model.Name <- paste0(RunID, 
                                "_n", test.frame$pop.size,
                                "_l", test.frame$loci.number,
                                "_s", test.frame$sel.thresh,
                                "_m", test.frame$modifier.found.freq,
                                "_g", test.frame$generations,
                                "_i", test.frame$iterations,
                                "_" , test.frame$modifier)


model.results <- dir()
model.results <- gsub(".Rdata", "", model.results)

test.frame$Ran <- ifelse(test.frame$Model.Name %in% model.results, "yes", "no")
table(test.frame$Ran)


test.frame <- subset(test.frame, Ran == "yes")
head(test.frame)

rm(model.results)

#~~ Load data and summarise by differences in modifier locus frequencies at RNF212


for(i in 1:nrow(test.frame)){
  
  if(i %in% seq(1, nrow(test.frame), 10)) print(paste("Processing line", i, "of", nrow(test.frame)))
  
  load(paste0(test.frame$Model.Name[i], ".Rdata"))
  sim.1 <- rbindlist(sim.1.list)
  head(sim.1)
  
  
  #~~ Genreate summary tables
  
  temp1 <- tapply(sim.1$PHENO, list(sim.1$GEN, sim.1$Simulation), mean)
  temp2 <- tapply(sim.1$modifier, list(sim.1$GEN, sim.1$Simulation), mod.freq)
  temp3 <- tapply(sim.1$PHENO, list(sim.1$GEN, sim.1$Simulation), var)
  
  temp1 <- melt(temp1)
  temp2 <- melt(temp2)
  temp3 <- melt(temp3)
  names(temp1) <- c("GEN", "Simulation", "MeanPHENO")
  names(temp2) <- c("GEN", "Simulation", "ModifierFreq")
  names(temp3) <- c("GEN", "Simulation", "VarPHENO")
  
  
  df1 <- join(temp1, temp2)
  df1 <- join(df1, temp3)
  
  rm(temp1, temp2, temp3)
  
  df1$ModifierStart <- test.frame$modifier.found.freq[i]
  
  
  temp1 <- tapply(sim.1$PHENO, sim.1$GEN, mean)
  temp2 <- tapply(sim.1$modifier, sim.1$GEN, mod.freq)
  temp3 <- tapply(sim.1$PHENO, sim.1$GEN, var)
  
  
  temp1 <- melt(temp1)
  temp2 <- melt(temp2)
  temp3 <- melt(temp3)
  
  names(temp1) <- c("GEN", "MeanPHENO")
  names(temp2) <- c("GEN", "ModifierFreq")
  names(temp3) <- c("GEN", "VarPHENO")
  
  df2 <- join(temp1, temp2)
  df2 <- join(df2, temp3)
  
  rm(temp1, temp2, temp3)
  
  df2$ModifierStart <- test.frame$modifier.found.freq[i]
  
  save(df1, df2, file = paste0(test.frame$Model.Name[i], "processed.Rdata"))
  
}
