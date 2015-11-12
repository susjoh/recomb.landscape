
library(rbenchmark)

map <- read.table("data/soay_map.txt", header = T)
map.dist <- diff(map$cM.Position.Sex.Averaged[seq(1, nrow(map), 10)])
map.dist <- map.dist[which(map.dist >= 0 & map.dist < 2)]
maf.info <- map$MAF


# test values for optimisation:
n.found.hap     <- 100    # Number of founder haplotypes generated
n.loci          <- 100    # Number of loci underlying the trait
n.f             <- 100    # Number of females
n.m             <- 100    # Number of males
f.RS            <- 2      # Number of offspring per female
f.RS.Pr         <- 1      # Probability of number of offspring per female.
sel.thresh.f    <- 1      # Selection threshold
sel.thresh.m    <- 0.2    # Selection threshold
prdm9.found.maf <- 0.4    # Starting freq of PRDM9
n.generations   <- 100
n.iterations    <- 100
return.haplos   <- TRUE
restart.on.extinction <- TRUE
progressBar = TRUE
SaveOnExtinction = FALSE


# simPopulationLandscape<- function(
#   map.dist,
#   maf.info,
#   n.found.hap = 100,
#   n.loci,
#   n.f,
#   n.m,
#   f.RS,
#   f.RS.Pr = NULL,
#   sel.thresh.f,
#   sel.thresh.m,
#   prdm9.found.maf, 
#   n.generations, 
#   return.haplos = FALSE,
#   progressBar = TRUE,
#   SaveOnExtinction = FALSE){

#~~ sample two landscapes and initial frequencies of alleles in founders


benchmark(r1   <- sample(map.dist/100, n.loci))  #0.02
benchmark(r2   <- sample(map.dist/100, n.loci))  #0.02
rhet <- (r1 + r2)/2

map.list <- list(r1, rhet, r2)

#~~ fix MAFs so that they show a range of frequencies between 0 & 1

benchmark(mafs.found <- sample(maf.info, n.loci))  #0.04
mafs.found <- mafs.found + (runif(n.loci) < 0.5)/2

#~~ determine PRDM9 genotype frequencies based on HWE

p <- prdm9.found.maf
q <- 1 - p
prdm9.found.prs <- c(p^2, 2*p*q, q^2)

#~~ Run Simulation

results.list <- list()
haplo.list <- list()



# NB Sex 1 = male, 2 = female (Xy, XX)

#~~ generate founder haplotypes

benchmark(founder.haplos <- lapply(1:n.found.hap, function (x) (runif(n.loci) < mafs.found) + 0L)) #0.25

#~~ generate diplotypes for n.f females and n.m males

gen.0 <- list()
gen.0[1:(n.f + n.m)] <- list(list(MOTHER = NA, FATHER = NA))


benchmark({       # 0.64  
  gen.0 <- list()
  gen.0[1:(n.f + n.m)] <- list(list(MOTHER = NA, FATHER = NA))
  
  for(i in 1:(n.f + n.m)){
    gen.0[[i]]["MOTHER"] <- sample(founder.haplos, size = 1)
    gen.0[[i]]["FATHER"] <- sample(founder.haplos, size = 1)
  }
})



#~~ create reference table
# 3.93 to 0.89
benchmark(ref.0 <- data.frame(GEN         = 0,
                              ID          = 1:length(gen.0),
                              MOTHER      = NA,
                              FATHER      = NA,
                              SEX         = sapply(1:length(gen.0), function(x) (runif(1) < 0.5) + 1L),
                              PRDM9       = sapply(1:length(gen.0), function(x) sample(1:3, size = 1, prob = prdm9.found.prs)),
                              PHENO       = sapply(1:length(gen.0), function(x) sum(gen.0[[x]][[1]]) + sum(gen.0[[x]][[2]]))))

# 0.11 and 0.11
benchmark(m.thresh <- sort(subset(ref.0, SEX == 1)$PHENO)[(1-sel.thresh.m)*length(subset(ref.0, SEX == 1)$PHENO)])
benchmark(f.thresh <- sort(subset(ref.0, SEX == 2)$PHENO)[(1-sel.thresh.f)*length(subset(ref.0, SEX == 2)$PHENO)])
if(length(m.thresh) == 0) m.thresh <- 0
if(length(f.thresh) == 0) f.thresh <- 0

#~~ remove IDs that will not be selected

ref.0$Bred <- 0
benchmark(ref.0$Bred[sort(c(which(ref.0$SEX == 1 & ref.0$PHENO >= m.thresh), #0.03
                            which(ref.0$SEX == 2 & ref.0$PHENO >= f.thresh)))] <- 1)


results.list[[1]] <- ref.0
haplo.list[[1]] <- gen.0

#~~ generate spaces for diplotypes of two offspring per female and sample best fathers

benchmark(if(progressBar == TRUE) pb = txtProgressBar(min = 1, max = n.generations, style = 3) )


for(gen in 1:n.generations){
  
  if(progressBar == TRUE) setTxtProgressBar(pb,gen)
  
  length.out <- f.RS*length(which(ref.0$SEX == 2 & ref.0$Bred == 1))
  
  if(any(c(length.out == 0, length(unique(ref.0$SEX)) == 1) == TRUE)) {
    if(SaveOnExtinction == TRUE){
      if(return.haplos == TRUE){
        return(list(results = results.list, haplos = haplo.list, maps = map.list))
      }
      if(return.haplos == FALSE){
        return(list(results = results.list, maps = map.list))
      }
    }
    stop(paste("Population has gone extinct at generation", gen))
  }
  
  ref.1 <- data.frame(GEN         = gen, #0.34
                      ID          = 1:length.out,
                      MOTHER      = rep(ref.0$ID[which(ref.0$SEX == 2 & ref.0$Bred == 1)], each = f.RS),
                      FATHER      = sample(ref.0$ID[which(ref.0$SEX == 1 & ref.0$Bred == 1)], size = length.out, replace = T),
                      SEX         = (runif(length.out) < 0.5) + 1L,
                      PRDM9       = NA,
                      PHENO       = NA)
  
  #~~ Transmit a gamete from parents to offspring
  
  gen.1 <- list()
  gen.1[1:length.out] <- list(list(MOTHER = NA, FATHER = NA))
  
  for(i in 1:length.out){
    
    # 3.98 seconds for 10K replications up to PRDM9 part
    #~~ MOTHER ~~#
    # 0.27 and 0.65 at 10K replicates
    haplos <- gen.0   [[ref.1$MOTHER[i]]]
    rmap   <- map.list[[ref.0$PRDM9 [which(ref.0$ID == ref.1$MOTHER[i])]]]
    
    #~~ sample crossover positions
    
    
    rec.pos <- which(((runif(length(rmap)) < rmap) + 0L) == 1) #0.25 @ 10K
    if(length(rmap) %in% rec.pos) rec.pos <- rec.pos[-which(rec.pos == length(rmap))]
    if(length(rec.pos) == 0) gen.1[[i]]["MOTHER"] <- haplos[sample.int(2, 1)] #0.25 @ 10K
    
    
    if(length(rec.pos) > 0){   #0.58 @ 10K
      
      haplos <- haplos[sample.int(2, 2, replace = F)] #0.14
      
      start.pos <- c(1, rec.pos[1:(length(rec.pos))] + 1)  #0.08
      stop.pos <- c(rec.pos, length(rmap)) #0.05
      
      fragments <- list()
      
      for(k in 1:length(start.pos)){  #0.03
        if(k %% 2 != 0) fragments[[k]] <- haplos[[1]][start.pos[k]:stop.pos[k]]
        if(k %% 2 == 0) fragments[[k]] <- haplos[[2]][start.pos[k]:stop.pos[k]]
      }
      
      gen.1[[i]]["MOTHER"] <- list(unlist(fragments)) #0.14
      
    }
    
    #~~ FATHER ~~#
    
    haplos <- gen.0   [[ref.1$FATHER[i]]]
    rmap   <- map.list[[ref.0$PRDM9 [which(ref.0$ID == ref.1$FATHER[i])]]]
    
    #~~ sample crossover positions
    
    rec.pos <- which(((runif(length(rmap)) < rmap) + 0L) == 1)
    if(length(rmap) %in% rec.pos) rec.pos <- rec.pos[-which(rec.pos == length(rmap))]
    if(length(rec.pos) == 0) gen.1[[i]]["FATHER"] <- haplos[sample.int(2, 1)]
    
    
    if(length(rec.pos) > 0){
      
      haplos <- haplos[sample.int(2, 2, replace = F)]
      
      start.pos <- c(1, rec.pos[1:(length(rec.pos))] + 1)
      stop.pos <- c(rec.pos, length(rmap))
      
      fragments <- list()
      
      for(k in 1:length(start.pos)){
        if(k %% 2 != 0) fragments[[k]] <- haplos[[1]][start.pos[k]:stop.pos[k]]
        if(k %% 2 == 0) fragments[[k]] <- haplos[[2]][start.pos[k]:stop.pos[k]]
      }
      
      gen.1[[i]]["FATHER"] <- list(unlist(fragments))
      
    }
    
    #~~ Deal with PRDM9
    
    
    prdm9.mum <- ref.0$PRDM9[which(ref.0$ID == ref.1$MOTHER[i])] #0.71
    prdm9.mum.2 <- ifelse(prdm9.mum == 1, 0,
                          ifelse(prdm9.mum == 3, 1,
                                 ifelse(prdm9.mum == 2, (runif(1) < 0.5) + 0L, NA)))  #0.14
    
    prdm9.dad <- ref.0$PRDM9[which(ref.0$ID == ref.1$FATHER[i])]
    prdm9.dad.2 <- ifelse(prdm9.dad == 1, 0,
                          ifelse(prdm9.dad == 3, 1,
                                 ifelse(prdm9.dad == 2, (runif(1) < 0.5) + 0L, NA)))
    
    ref.1$PRDM9[i] <- prdm9.mum.2 + prdm9.dad.2 + 1
    
  }
  
  
  rm(haplos, rmap, rec.pos, start.pos, stop.pos, fragments, k, prdm9.mum, prdm9.mum.2, prdm9.dad, prdm9.dad.2)
  
  ref.1$PHENO <- sapply(1:length(gen.1), function(x) sum(gen.1[[x]][[1]]) + sum(gen.1[[x]][[2]]))  #2.92
  
  #~~ Deal with IDs that will be selected
  
  #4.27
  m.thresh <- sort(ref.1$PHENO[which(ref.1$SEX == 1)])[(1-sel.thresh.m)*length(ref.1$PHENO[which(ref.1$SEX == 1)])]
  f.thresh <- sort(ref.1$PHENO[which(ref.1$SEX == 2)])[(1-sel.thresh.m)*length(ref.1$PHENO[which(ref.1$SEX == 2)])]
  
  if(length(m.thresh) == 0) m.thresh <- 0
  if(length(f.thresh) == 0) f.thresh <- 0
  
  #~~ remove IDs that will not be selected
  
  ref.1$Bred <- 0
  ref.1$Bred[sort(c(which(ref.1$SEX == 1 & ref.1$PHENO >= m.thresh),
                    which(ref.1$SEX == 2 & ref.1$PHENO >= f.thresh)))] <- 1
  
  results.list[[(gen + 1)]] <- ref.1
  if(return.haplos == TRUE) haplo.list[[(gen + 1)]] <- gen.1
  
  gen.0 <- gen.1
  ref.0 <- ref.1
  
}


#~~ Parse output

if(return.haplos == TRUE){
  list(results = results.list, haplos = haplo.list, maps = map.list)
} else {  
  list(results = results.list, maps = map.list)
}
}

