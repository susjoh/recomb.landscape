#' Simulate a population with directional selection, recombination and landscape
#' variation
#' 
#' @param map.dist A vector of cM distances between adjacent loci that are
#'   sampled to create the recombination landscape. Can be any length greater
#'   than the number of loci.
#' @param maf.info A vector of minor allele frequencies of trait loci. Can be
#'   any length greater than the number of loci.
#' @param n.found.hap The number of haplotypes in the founder population.
#'   Default is 100.
#' @param n.loci The number of loci contributing to phenotype. Assumes loci of
#'   equal and additive effect on phenotype.
#' @param n.f The number of females in the founder population.
#' @param n.m The number of males in the founder population.
#' @param f.RS Female reproductive success. At present simulations only run with
#'   a specific number of offspring for each female. This will be modified in
#'   future.
#' @param f.RS.Pr Not currently used
#' @param sel.thresh.f Selection threshold in females (value between 0 and 1,
#'   where 1 is all females selected)
#' @param sel.thresh.m Selection threshold in males (value between 0 and 1,
#'   where 1 is all males selected)
#' @param prdm9.found.maf The minor allele frequence at the PRDM9 locus, which
#'   will modify recombination landscape.
#' @param n.generations Number of generations to run the simulation
#' @param return.haplos Should the haplotype information be returned? Default =
#'   FALSE as will return a large amount of data!
#' @param verbose (Default = TRUE) Should information on progress be printed?


#test values for optimisation:
# n.found.hap     <- 100    # Number of founder haplotypes generated
# n.loci          <- 100    # Number of loci underlying the trait
# n.f             <- 100    # Number of females
# n.m             <- 100    # Number of males
# f.RS            <- 2      # Number of offspring per female
# f.RS.Pr         <- 1      # Probability of number of offspring per female.
# sel.thresh.f    <- 1      # Selection threshold
# sel.thresh.m    <- 0.2    # Selection threshold
# prdm9.found.maf <- 0.4    # Starting freq of PRDM9
# n.generations   <- 100
# n.iterations    <- 100
# return.haplos   <- TRUE
# restart.on.extinction <- TRUE
# progressBar = TRUE
# 
# map <- read.table("data/soay_map.txt", header = T)
# map.dist <- diff(map$cM.Position.Sex.Averaged[seq(1, nrow(map), 10)])
# map.dist <- map.dist[which(map.dist >= 0 & map.dist < 2)]
# maf.info <- map$MAF
# SaveOnExtinction = TRUE


simPopulationLandscape<- function(
  map.dist,
  maf.info,
  n.found.hap = 100,
  n.loci,
  n.f,
  n.m,
  f.RS,
  f.RS.Pr = NULL,
  sel.thresh.f,
  sel.thresh.m,
  prdm9.found.maf, 
  n.generations, 
  return.haplos = FALSE,
  progressBar = TRUE,
  SaveOnExtinction = FALSE,
  FounderObject = NULL){
  
  #~~ sample two landscapes and initial frequencies of alleles in founders
  
  if(is.null(FounderObject)){
    r1   <- sample(map.dist/100, replace = T, n.loci)
    r2   <- sample(map.dist/100, replace = T, n.loci)
    rhet <- (r1 + r2)/2 
    
      ggplot(data.frame(r = c(r1, rhet, r2), map = rep(1:3, each = length(r1)), x = rep(1:length(r1), times = 3)),
             aes(x, r, colour = factor(map))) +
        geom_line() +
        scale_colour_brewer(palette = "Set1")
    
    map.list <- list(r1, rhet, r2)
    
    #~~ fix MAFs so that they show a range of frequencies between 0 & 1
    
    mafs.found <- sample(maf.info, n.loci)
    mafs.found <- mafs.found + (runif(n.loci) < 0.5)/2
    
    #~~ determine PRDM9 genotype frequencies based on HWE
    
    p <- prdm9.found.maf
    q <- 1 - p
    prdm9.found.prs <- c(p^2, 2*p*q, q^2)
    
    #~~ Run Simulation
    
    
    
    # NB Sex 1 = male, 2 = female (Xy, XX)
    
    #~~ generate founder haplotypes
    
    founder.haplos <- lapply(1:n.found.hap, function (x) (runif(n.loci) < mafs.found) + 0L)
    
    #~~ generate diplotypes for n.f females and n.m males
    
    gen.0 <- list()
    gen.0[1:(n.f + n.m)] <- list(list(MOTHER = NA, FATHER = NA))
    
    for(i in 1:(n.f + n.m)){
      gen.0[[i]]["MOTHER"] <- sample(founder.haplos, size = 1)
      gen.0[[i]]["FATHER"] <- sample(founder.haplos, size = 1)
    }
    
    #~~ create reference table
    
    ref.0 <- data.frame(GEN         = 0,
                        ID          = 1:length(gen.0),
                        MOTHER      = NA,
                        FATHER      = NA,
                        #SEX         = sapply(1:length(gen.0), function(x) (runif(1) < 0.5) + 1L),
                        SEX         = rep(1:2, length.out = length(gen.0)),
                        PHENO       = sapply(1:length(gen.0), function(x) sum(gen.0[[x]][[1]]) + sum(gen.0[[x]][[2]])),
                        PRDM9       = sapply(1:length(gen.0), function(x) sample(1:3, size = 1, prob = prdm9.found.prs)))
    
  } else {
    
    if(any(n.f != FounderObject$n.f,
           n.m != FounderObject$n.m,
           n.loci != FounderObject$n.loci)) stop("Founder Object parameters do not match those of current simulation")
    
    ref.0 <- FounderObject$ref.0
    founder.haplos <- founder.haplos <- FounderObject$founder.haplos
    gen.0 <- FounderObject$gen.0
    map.list<- FounderObject$map.list

    
    #~~ determine PRDM9 genotype frequencies based on HWE
    
    p <- prdm9.found.maf
    q <- 1 - p
    prdm9.found.prs <- c(p^2, 2*p*q, q^2)
    
    ref.0$PRDM9 = sapply(1:length(gen.0), function(x) sample(1:3, size = 1, prob = prdm9.found.prs))        
  }
  
  m.thresh <- sort(ref.0$PHENO[which(ref.0$SEX == 1)])[(1-sel.thresh.m)*length(ref.0$PHENO[which(ref.0$SEX == 1)])]
  f.thresh <- sort(ref.0$PHENO[which(ref.0$SEX == 2)])[(1-sel.thresh.f)*length(ref.0$PHENO[which(ref.0$SEX == 2)])]
  
  if(length(m.thresh) == 0) m.thresh <- 0
  if(length(f.thresh) == 0) f.thresh <- 0
  
  #~~ remove IDs that will not be selected
  
  ref.0$Bred <- 0
  ref.0$Bred[sort(c(which(ref.0$SEX == 1 & ref.0$PHENO >= m.thresh),
                    which(ref.0$SEX == 2 & ref.0$PHENO >= f.thresh)))] <- 1
  
  
  results.list <- list()
  haplo.list <- list()
  
  results.list[[1]] <- ref.0
  haplo.list[[1]] <- gen.0
  
  #~~ generate spaces for diplotypes of two offspring per female and sample best fathers
  
  if(progressBar == TRUE) pb = txtProgressBar(min = 1, max = n.generations, style = 3) 
  
  for(gen in 1:n.generations){
    
    if(progressBar == TRUE) setTxtProgressBar(pb,gen)
    
    length.out <- f.RS*length(which(ref.0$SEX == 2 & ref.0$Bred == 1))
    
    if(any(c(length.out == 0, length(unique(ref.0$SEX)) == 1) == TRUE)) {
      if(SaveOnExtinction == TRUE){
        if(return.haplos == TRUE){
          return(list(results = results.list, haplos = haplo.list, maps = map.list))
        } else {  
          return(list(results = results.list, maps = map.list))
        }
      }
      
      break(paste("Population has gone extinct at generation", gen))
    }
    
    
    ref.1 <- data.frame(GEN         = gen,
                        ID          = 1:length.out,
                        MOTHER      = rep(ref.0$ID[which(ref.0$SEX == 2 & ref.0$Bred == 1)], each = f.RS),
                        FATHER      = sample(ref.0$ID[which(ref.0$SEX == 1 & ref.0$Bred == 1)], size = length.out, replace = T),
                        SEX         = (runif(length.out) < 0.5) + 1L,
                        PHENO       = NA,
                        PRDM9       = NA)
    
    #~~ Transmit a gamete from parents to offspring
    
    gen.1 <- list()
    gen.1[1:length.out] <- list(list(MOTHER = NA, FATHER = NA))
    
    for(i in 1:length.out){
      
      #~~ MOTHER ~~#
      
      haplos <- gen.0   [[ref.1$MOTHER[i]]]
      rmap   <- map.list[[ref.0$PRDM9 [which(ref.0$ID == ref.1$MOTHER[i])]]]
      
      #~~ sample crossover positions
      
      rec.pos <- which(((runif(length(rmap)) < rmap) + 0L) == 1)
      if(length(rmap) %in% rec.pos) rec.pos <- rec.pos[-which(rec.pos == length(rmap))]
      if(length(rec.pos) == 0) gen.1[[i]]["MOTHER"] <- haplos[sample.int(2, 1)]
      
      
      if(length(rec.pos) > 0){
        
        haplos <- haplos[sample.int(2, 2, replace = F)]
        
        start.pos <- c(1, rec.pos[1:(length(rec.pos))] + 1)
        stop.pos <- c(rec.pos, length(rmap))
        
        fragments <- list()
        
        for(k in 1:length(start.pos)){
          if(k %% 2 != 0) fragments[[k]] <- haplos[[1]][start.pos[k]:stop.pos[k]]
          if(k %% 2 == 0) fragments[[k]] <- haplos[[2]][start.pos[k]:stop.pos[k]]
        }
        
        gen.1[[i]]["MOTHER"] <- list(unlist(fragments))
        
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
      
      prdm9.mum <- ref.0$PRDM9[which(ref.0$ID == ref.1$MOTHER[i])]
      prdm9.mum.2 <- ifelse(prdm9.mum == 1, 0,
                            ifelse(prdm9.mum == 3, 1,
                                   ifelse(prdm9.mum == 2, (runif(1) < 0.5) + 0L, NA)))
      
      prdm9.dad <- ref.0$PRDM9[which(ref.0$ID == ref.1$FATHER[i])]
      prdm9.dad.2 <- ifelse(prdm9.dad == 1, 0,
                            ifelse(prdm9.dad == 3, 1,
                                   ifelse(prdm9.dad == 2, (runif(1) < 0.5) + 0L, NA)))
      
      ref.1$PRDM9[i] <- prdm9.mum.2 + prdm9.dad.2 + 1
      
      
    }
    
    rm(haplos, rmap, rec.pos, start.pos, stop.pos, fragments, k, prdm9.mum, prdm9.mum.2, prdm9.dad, prdm9.dad.2)
    
    ref.1$PHENO <- sapply(1:length(gen.1), function(x) sum(gen.1[[x]][[1]]) + sum(gen.1[[x]][[2]]))
    
    
    #~~ Deal with IDs that will be selected
    
    m.thresh <- sort(ref.1$PHENO[which(ref.1$SEX == 1)])[(1-sel.thresh.m)*length(ref.1$PHENO[which(ref.1$SEX == 1)])]
    f.thresh <- sort(ref.1$PHENO[which(ref.1$SEX == 2)])[(1-sel.thresh.f)*length(ref.1$PHENO[which(ref.1$SEX == 2)])]
    
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




