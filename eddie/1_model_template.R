

library(plyr)
library(data.table)


sapply(paste0("R/", dir("R", pattern = ".R")), source)


#~~ iterations

sim.1.list <- list()

for(it in 1:iterations){
  
  if(it %in% seq(1, iterations, 10)) print(paste("Running iteration", it, "of", iterations))
  #~~ generate map
  
  if(modifier == "RNF212"){
    x <- rnbinom(n.loci, 1, 0.1) # sample allele frequencies from negative binomial distribution/
    x <- x/sum(x)             # standardised to 1 Morgan
    
    het <- mean(c(1, rate.increase))
    
    map.list <- list(x, x*het, x*rate.increase)
    rm(x, het)  
  }
  
  if(modifier == "PRDM9"){
    x <- rnbinom(n.loci, 1, 0.1) # sample allele frequencies from negative binomial distribution/
    x <- x/sum(x)             # standardised to 1 Morgan
    
    y <- rnbinom(n.loci, 1, 0.1) # sample allele frequencies from negative binomial distribution/
    y <- y/sum(y)             # standardised to 1 Morgan
    
    map.list <- list(x, apply(data.frame(x, y), 1, mean), y)
    rm(x, y)  
  }
  
  #~~ generate allele frequencies
  
  freq.vec <- runif(n.loci, 0, 1)
  
  
  sim.1 <- simPopulationResponse(n.found.hap = n.found.hap,
                                 n.loci = n.loci,
                                 n.f = n.females, n.m = n.males,
                                 map.list = map.list,
                                 allele.freqs = freq.vec,
                                 f.RS = n.offspring,
                                 sel.thresh.f = female.sel.thresh,
                                 sel.thresh.m = male.sel.thresh,
                                 modifier.found.freq = modifier.found.freq,
                                 n.generations = generations,
                                 force.equal.male.success = T,
                                 force.equal.sex = T,
                                 progressBar = F)
  sim.1.list[[it]] <- sim.1
}

save(sim.1, file = paste0(model.name, ".Rdata"))
