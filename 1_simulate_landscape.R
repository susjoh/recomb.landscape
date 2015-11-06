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

map <- read.table("data/soay_map.txt", header = T)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 1. Define simulation parameters and sampling distributions   #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~ Create the sampling distrubutions

map.dist <- diff(map$cM.Position.Sex.Averaged[seq(1, nrow(map), 10)])
map.dist <- map.dist[which(map.dist >= 0 & map.dist < 2)]

maf.info <- map$MAF


simPopulationLandscape(map.dist = map.dist, maf.info = maf.info, n.found.hap = 100, n.loci = 100,
                       n.f = 100, n.m = 100, f.RS = 2, sel.thresh.f = 1, sel.thresh.m = 0.2,
                       prdm9.found.maf = 0.4, n.generations = 100)



require(data.table)
test <- rbindlist(results.list)

ggplot(test, aes(as.factor(GEN), PHENO)) + geom_boxplot()
ggplot(test, aes(GEN, PRDM9)) + geom_point() + stat_smooth()

