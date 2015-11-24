library(plyr)
library(data.table)

sapply(paste0("R/", dir("R", pattern = ".R")), source)


n.females         <- 50  # Number of females in the founder generation
n.males           <- 50  # Number of males in the founder generation
n.offspring       <- 2   # Number of offspring per female
male.sel.thresh   <- 0.4 # Top proportion of males that will be selected
female.sel.thresh <- 1   # Top proportion of females that will be selected
generations       <- 100 # Number of generations to run the simulation.

n.loci <- 100      # The number of loci contributing to the trait
n.found.hap <- 100 # Number of founder haplotypes


map.1  <- rep(1/n.loci, n.loci)
freq.1 <- rep(0.5, n.loci)

sim.1 <- simPopulationResponse(n.found.hap = n.found.hap,
                               n.loci = n.loci,
                               n.f = n.females, n.m = n.males,
                               map.list = list(map.1),
                               allele.freqs = freq.1,
                               f.RS = n.offspring,
                               sel.thresh.f = female.sel.thresh,
                               sel.thresh.m = male.sel.thresh,
                               modifier.found.freq = 0,
                               n.generations = generations,
                               force.equal.male.success = T,
                               force.equal.sex = T,
                               progressBar = F)

sim.1.res <- rbindlist(sim.1$results)

sim.1.res

df1 <- ddply(sim.1.res, .(GEN), summarise, MeanPHENO = mean(PHENO))
