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


# test values for optimisation:
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


simPopulationLandscape <- function(
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
  verbose = TRUE){
  
  #~~ sample two landscapes and initial frequencies of alleles in founders
  
  r1   <- sample(map.dist/100, n.loci)
  r2   <- sample(map.dist/100, n.loci)
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
  
  results.list <- list()
  haplo.list <- list()
  
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
                      SEX         = sapply(1:length(gen.0), function(x) (runif(1) < 0.5) + 1L),
                      PRDM9       = sapply(1:length(gen.0), function(x) sample(1:3, size = 1, prob = prdm9.found.prs)),
                      PHENO       = sapply(1:length(gen.0), function(x) sum(unlist(gen.0[[x]]))))
  
  m.thresh <- sort(subset(ref.0, SEX == 1)$PHENO)[(1-sel.thresh.m)*length(subset(ref.0, SEX == 1)$PHENO)]
  f.thresh <- sort(subset(ref.0, SEX == 2)$PHENO)[(1-sel.thresh.f)*length(subset(ref.0, SEX == 2)$PHENO)]
  if(length(m.thresh) == 0) m.thresh <- 0
  if(length(f.thresh) == 0) f.thresh <- 0
  
  #~~ remove IDs that will not be selected
  
  ref.0$Bred <- 0
  ref.0$Bred[sort(c(which(ref.0$SEX == 1 & ref.0$PHENO >= m.thresh),
                    which(ref.0$SEX == 2 & ref.0$PHENO >= f.thresh)))] <- 1
  
  
  results.list[[1]] <- ref.0
  haplo.list[[1]] <- gen.0
  
  #~~ generate spaces for diplotypes of two offspring per female and sample best fathers
  
  for(gen in 1:n.generations){
    
    if(verbose == TRUE) print(paste("Generation", gen))
    
    length.out <- f.RS*nrow(subset(ref.0, SEX == 2 & Bred == 1))
    
    ref.1 <- data.frame(GEN         = gen,
                        ID          = 1:length.out,
                        MOTHER      = rep(subset(ref.0, SEX == 2 & Bred == 1)$ID, each = f.RS),
                        FATHER      = sample(subset(ref.0, SEX == 1 & Bred == 1)$ID, size = length.out, replace = T),
                        SEX         = sapply(1:length.out, function(x) (runif(1) < 0.5) + 1L),
                        PRDM9       = NA,
                        PHENO       = NA)
    
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
    
    ref.1$PHENO <- sapply(1:length(gen.1), function(x) sum(unlist(gen.1[[x]])))
    
    
    #~~ Deal with IDs that will be selected
    
    m.thresh <- sort(subset(ref.1, SEX == 1)$PHENO)[(1-sel.thresh.m)*length(subset(ref.1, SEX == 1)$PHENO)]
    f.thresh <- sort(subset(ref.1, SEX == 2)$PHENO)[(1-sel.thresh.f)*length(subset(ref.1, SEX == 2)$PHENO)]
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
    list(results = results.list, haplos = haplo.list)
  } else {  
    results.list
  }
}



