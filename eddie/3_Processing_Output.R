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

source("C:/Users/Susan Johnston/Desktop/R Functions/multiplot.R")

mod.freq <- function(vec) {sum(vec - 1)/(2*length(vec))}

vec <- subset(sim.1, GEN == 1 & Simulation == 1)$modifier

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


group_frame <- subset(test.frame, pop.size == 50 & loci.number == 50 & sel.thresh == 0.8 & modifier == "RNF212")

head(group_frame)

df1_cat <- NULL
df2_cat <- NULL

for(i in 1:nrow(group_frame)){
  
  load(paste0(group_frame$Model.Name[i], ".Rdata"))
  sim.1 <- rbindlist(sim.1.list)
  head(sim.1)
  
  
  #~~ Genreate summary tables
  
  temp1 <- tapply(sim.1$PHENO, list(sim.1$GEN, sim.1$Simulation), mean)
  temp2 <- tapply(sim.1$modifier, list(sim.1$GEN, sim.1$Simulation), mod.freq)
  temp1 <- melt(temp1)
  temp2 <- melt(temp2)
  names(temp1) <- c("GEN", "Simulation", "MeanPHENO")
  names(temp2) <- c("GEN", "Simulation", "ModifierFreq")
  
  df1 <- join(temp1, temp2)
  rm(temp1, temp2)
  
  df1$ModifierStart <- group_frame$modifier.found.freq[i]
    
  df1_cat <- rbind(df1_cat, df1)
  
  temp1 <- tapply(sim.1$PHENO, sim.1$GEN, mean)
  temp2 <- tapply(sim.1$modifier, sim.1$GEN, mod.freq)
  temp1 <- melt(temp1)
  temp2 <- melt(temp2)
  names(temp1) <- c("GEN", "MeanPHENO")
  names(temp2) <- c("GEN", "ModifierFreq")
  
  df2 <- join(temp1, temp2)
  rm(temp1, temp2)
  
  df2$ModifierStart <- group_frame$modifier.found.freq[i]
  
  df2_cat <- rbind(df2_cat, df2)

}

df1_cat$ModifierStart <- as.character(df1_cat$ModifierStart)

ggplot(df1_cat, aes(GEN, MeanPHENO, colour = as.factor(ModifierStart), group = interaction(ModifierStart, Simulation))) +
  geom_line() +
  scale_colour_brewer(palette = "Set1") +
  labs(title = "Mean phenotype per simulation")

ggplot(df1_cat, aes(GEN, ModifierFreq, group = Simulation)) +
  geom_line() +
  scale_colour_brewer(palette = "Set1") +
  labs(title = "Mean modifier frequency per simulation") +
  facet_wrap(~ModifierStart)


ggplot(df2_cat, aes(GEN, MeanPHENO, colour = as.factor(ModifierStart))) +
  geom_line() +
  scale_colour_brewer(palette = "Set1") +
  labs(title = "Mean phenotype per simulation")

ggplot(df2_cat, aes(GEN, ModifierFreq, colour = as.factor(ModifierStart))) +
  geom_line() +
  scale_colour_brewer(palette = "Set1") +
  labs(title = "Mean modifier frequency per simulation")



