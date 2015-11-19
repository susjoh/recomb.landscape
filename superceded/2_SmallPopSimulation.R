#
# Simulation of recombination landscape
# Susan E. Johnston
# Started: 5th November 2015
#
#

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 0. Set up working environment                                #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


library(ggplot2)
library(plyr)
library(data.table)

map <- read.table("data/soay_map.txt", header = T)

sapply(paste0("R/", dir("R")), source)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 1. Define simulation parameters and sampling distributions   #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


generations <- 200
no.females <- 100
no.males <- 100
no.loci <- 1000
no.offspring <- 2
no.founder.hap <- 100
male.sel.thresh <- 0.4
iterations <- 50
sampling.int <- 1


#~~ Create the sampling distrubutions

map.dist <- diff(map$cM.Position.Sex.Averaged[seq(1, nrow(map), sampling.int)])
map.dist <- map.dist[which(map.dist >= 0 & map.dist < 2)]

maf.info <- map$MAF

res.list.prdm9.present <- list()
res.list.prdm9.absent  <- list()

system.time({
  for(i in 1:iterations){
    print(paste("Running simulation", i))
    
    x1 <- createFounderObject(map.dist = map.dist, maf.info = maf.info, n.found.hap = no.founder.hap,
                              n.loci = no.loci, n.f = no.females, n.m = no.males, f.RS = no.offspring, f.RS.Pr = NULL)
    
    res.list.prdm9.present[[i]] <- rbindlist(
      simPopulationLandscape(map.dist = map.dist, maf.info = maf.info,
                                              n.found.hap = no.founder.hap, n.loci = no.loci, 
                                              n.f = no.females, n.m = no.males, f.RS = no.offspring, 
                                              sel.thresh.f = 1, sel.thresh.m = male.sel.thresh,
                                              prdm9.found.maf = 0.4, n.generations = generations, FounderObject = x1)$results)
    
    res.list.prdm9.present[[i]]$Simulation <- i
    
    res.list.prdm9.absent[[i]] <- rbindlist(
      simPopulationLandscape(map.dist = map.dist, maf.info = maf.info,
                                              n.found.hap = no.founder.hap, n.loci = no.loci, 
                                              n.f = no.females, n.m = no.males, f.RS = no.offspring, 
                                              sel.thresh.f = 1, sel.thresh.m = male.sel.thresh,
                                              prdm9.found.maf = 0, n.generations = generations, FounderObject = x1)$results)
    
    res.list.prdm9.absent[[i]]$Simulation <- i
  }
})

save(res.list.prdm9.present, res.list.prdm9.absent, 
     file = paste0("results/g", generations, "_it", iterations,  "_f", no.females, "_m", no.males, "_o", no.offspring,
                  "_l", no.loci, "_msel", male.sel.thresh, ".Rdata"))


     
     
load(paste0("results/g", generations, "_it", iterations,  "_f", no.females, "_m", no.males, "_o", no.offspring,
            "_l", no.loci, "_msel", male.sel.thresh, ".Rdata"))

test <- rbind(cbind(PRDM9.Start = 0.4, rbindlist(res.list.prdm9.present)),
              cbind(PRDM9.Start = 0  , rbindlist(res.list.prdm9.absent)))

test$Simulation <- as.factor(test$Simulation)
test$PRDM9.Start <- as.factor(test$PRDM9.Start)
head(test)

df1 <- ddply(test, .(GEN, PRDM9.Start, Simulation), summarise, MeanPHENO = mean(PHENO))

df3 <- ddply(test, .(GEN, PRDM9.Start, Simulation), summarise, VarPHENO = var(PHENO))

setkey(test, GEN, PRDM9.Start, Simulation)

testfunc <- function(vec) {sum(vec - 1)/(2*length(vec))}

df2 <- test[,list(MeanPHENO=mean(PHENO),
                  VarPheno = var(PHENO),
                  PopSize = length(PHENO),
                  PRDM9.Freq = testfunc(PRDM9)),
            by=list(GEN, PRDM9.Start, Simulation)] 

head(df2)
test



ggplot(df2, aes(GEN, MeanPHENO, col = PRDM9.Start, group = interaction(Simulation, PRDM9.Start))) + 
  geom_line()

ggplot(df2, aes(GEN, VarPheno, col = PRDM9.Start, group = interaction(Simulation, PRDM9.Start))) + 
  geom_line()

ggplot(df2, aes(GEN, PopSize, col = PRDM9.Start, group = interaction(Simulation, PRDM9.Start))) + 
  geom_line()

ggplot(df2, aes(GEN, PRDM9.Freq, col = PRDM9.Start, group = interaction(Simulation, PRDM9.Start))) + 
  geom_line()



ggplot(df1, aes(GEN, MeanPHENO, col = PRDM9.Start)) +
  geom_point(alpha = 0) +
  stat_smooth() +
  theme(axis.text.x  = element_text (size = 16, vjust = 1),
        axis.text.y  = element_text (size = 16, hjust = 1),
        strip.text.x = element_text (size = 16, vjust = 0.7),
        axis.title.y = element_text (size = 16, angle = 90, vjust = 1),
        axis.title.x = element_text (size = 16, vjust = 0.2),
        strip.background = element_blank()) +
  scale_colour_brewer(palette = "Set1") +
  labs(title = paste0("g", generations, "_it", iterations,  "_f", no.females, "_m", no.males, "_o", no.offspring,
                      "_l", no.loci, "_msel", male.sel.thresh))

beepr::beep()