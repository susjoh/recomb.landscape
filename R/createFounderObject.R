#' Create a Founder object that can be used to start simulations with different PRDM9 values
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


# test values for optimisation:
# n.found.hap     <- 100    # Number of founder haplotypes generated
# n.loci          <- 100    # Number of loci underlying the trait
# n.f             <- 100    # Number of females
# n.m             <- 100    # Number of males
# f.RS            <- 2      # Number of offspring per female
# f.RS.Pr         <- 1      # Probability of number of offspring per female.
# 
# 
# map <- read.table("data/soay_map.txt", header = T)
# map.dist <- diff(map$cM.Position.Sex.Averaged[seq(1, nrow(map), 10)])
# map.dist <- map.dist[which(map.dist >= 0 & map.dist < 2)]
# maf.info <- map$MAF


createFounderObject<- function(
  map.dist,
  maf.info,
  n.found.hap = 100,
  n.loci,
  n.f,
  n.m,
  f.RS,
  f.RS.Pr = NULL){
  
  #~~ sample two landscapes and initial frequencies of alleles in founders
  
  r1   <- sample(map.dist/100, replace = T, n.loci)
  r2   <- sample(map.dist/100, replace = T, n.loci)
  rhet <- (r1 + r2)/2 
  
  #   ggplot(data.frame(r = c(r1, rhet, r2), map = rep(1:3, each = length(r1)), x = rep(1:length(r1), times = 3)),
  #          aes(x, r, colour = factor(map))) +
  #     geom_line() +
  #     scale_colour_brewer(palette = "Set1")
  
  map.list <- list(r1, rhet, r2)
  
  #~~ fix MAFs so that they show a range of frequencies between 0 & 1
  
  mafs.found <- sample(maf.info, n.loci)
  mafs.found <- mafs.found + (runif(n.loci) < 0.5)/2
  
  
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
                      #######                      PRDM9       = sapply(1:length(gen.0), function(x) sample(1:3, size = 1, prob = prdm9.found.prs)),
                      PHENO       = sapply(1:length(gen.0), function(x) sum(gen.0[[x]][[1]]) + sum(gen.0[[x]][[2]])))
  
  
  return(list(ref.0 = ref.0,
              founder.haplos = founder.haplos,
              gen.0 = gen.0,
              map.list = map.list,
              n.f = n.f,
              n.m = n.m,
              n.loci = n.loci))
  
  ##### START AT m.thresh
  
}



