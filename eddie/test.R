n.females           <- 50  
n.males             <- 50  
n.offspring         <- 2   
male.sel.thresh     <- 0.4 
female.sel.thresh   <- 1
generations         <- 100 
modifier.found.freq <- 0
n.loci <- 100
n.found.hap <- 100 
map.1  <- rep(1/n.loci, n.loci)
freq.1 <- rep(0.5, n.loci)



library(plyr)
library(data.table)

sapply(paste0("R/", dir("R", pattern = ".R")), source)



sim.1 <- simPopulationResponse(n.found.hap = n.found.hap,
                               n.loci = n.loci,
                               n.f = n.females, n.m = n.males,
                               map.list = list(map.1),
                               allele.freqs = freq.1,
                               f.RS = n.offspring,
                               sel.thresh.f = female.sel.thresh,
                               sel.thresh.m = male.sel.thresh,
                               modifier.found.freq = modifier.found.freq,
                               n.generations = generations,
                               force.equal.male.success = T,
                               force.equal.sex = T,
                               progressBar = F)

sim.1.res <- rbindlist(sim.1$results)

save(sim.1.res, "test.Rdata")