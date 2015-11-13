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

#~~ Create the sampling distrubutions

map.dist <- diff(map$cM.Position.Sex.Averaged[seq(1, nrow(map), 10)])
map.dist <- map.dist[which(map.dist >= 0 & map.dist < 2)]

maf.info <- map$MAF


#~~ Testing the basic functions

x <- restartOnExtinct(simPopulationLandscape(map.dist = map.dist, maf.info = maf.info,
                                             n.found.hap = 100, n.loci = 100, n.f = 100, 
                                             n.m = 100, f.RS = 2, sel.thresh.f = 1, sel.thresh.m = 0.2,
                                             prdm9.found.maf = 0.4, n.generations = 100)$results)

#~~ Create a founder population

x1 <- createFounderObject(map.dist = map.dist, maf.info = maf.info, n.found.hap = 100, n.loci = 100, n.f = 100, n.m = 100, f.RS = 2, f.RS.Pr = NULL)


#~~ Run with a founder population to test effects of PRDM9

x2 <- restartOnExtinct(simPopulationLandscape(map.dist = map.dist, maf.info = maf.info,
                                             n.found.hap = 100, n.loci = 100, n.f = 100, 
                                             n.m = 100, f.RS = 2, sel.thresh.f = 1, sel.thresh.m = 0.2,
                                             prdm9.found.maf = 0.4, n.generations = 100, FounderObject = x1)$results)


res.list.prdm9.present <- list()
res.list.prdm9.absent  <- list()

system.time({
  for(i in 1:100){
    print(paste("Running simulation", i))
    
    x1 <- createFounderObject(map.dist = map.dist, maf.info = maf.info, n.found.hap = 100,
                              n.loci = 100, n.f = 100, n.m = 100, f.RS = 2, f.RS.Pr = NULL)
    
    res.list.prdm9.present[[i]] <- rbindlist(
      restartOnExtinct(simPopulationLandscape(map.dist = map.dist, maf.info = maf.info,
                             n.found.hap = 100, n.loci = 100, n.f = 100, 
                             n.m = 100, f.RS = 2, sel.thresh.f = 1, sel.thresh.m = 0.2,
                             prdm9.found.maf = 0.4, n.generations = 100, FounderObject = x1)$results)
    )
    res.list.prdm9.present[[i]]$Simulation <- i
    
    res.list.prdm9.absent[[i]] <- rbindlist(
      restartOnExtinct(simPopulationLandscape(map.dist = map.dist, maf.info = maf.info,
                             n.found.hap = 100, n.loci = 100, n.f = 100, 
                             n.m = 100, f.RS = 2, sel.thresh.f = 1, sel.thresh.m = 0.2,
                             prdm9.found.maf = 0, n.generations = 100, FounderObject = x1)$results)
    )
    res.list.prdm9.absent[[i]]$Simulation <- i
  }
})

save(res.list.prdm9.present, res.list.prdm9.absent, file = "test.Rdata")


load("test.Rdata")
test <- rbind(cbind(PRDM9.Start = 0.4, rbindlist(res.list.prdm9.present)),
              cbind(PRDM9.Start = 0  , rbindlist(res.list.prdm9.absent)))
test$Simulation <- as.factor(test$Simulation)
test$PRDM9.Start <- as.factor(test$PRDM9.Start)
head(test)

df1 <- ddply(test, .(GEN, PRDM9.Start, Simulation), summarise, MeanPHENO = mean(PHENO))

setkey(test, GEN, PRDM9.Start, Simulation)

df2 <- test[,list(MeanPHENO=mean(PHENO), PopSize = length(PHENO)),
            by=list(GEN, PRDM9.Start, Simulation)] 



ggplot(df2, aes(GEN, PopSize, col = PRDM9.Start, group = interaction(Simulation))) + geom_line()


ggplot(df1, aes(GEN, MeanPHENO, col = PRDM9.Start, group = interaction(Simulation))) +
  geom_point(alpha = 0) +
  stat_smooth() +
  scale_colour_brewer(palette = "Set1") +
  facet_wrap(~PRDM9.Start)

