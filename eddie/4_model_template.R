


library(data.table)
library(plyr)
library(ggplot2)
library(reshape)

mod.freq <- function(vec) {sum(vec - 1)/(2*length(vec))}
  

load(paste0(x, ".Rdata"))
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

df1$ModifierStart <- mod.freq


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

df2$ModifierStart <- mod.freq

save(df1, df2, file = paste0(x, "processed.Rdata"))