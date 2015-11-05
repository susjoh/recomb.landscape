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

#~~ define starting values

n.found.hap     <- 100    # Number of founder haplotypes generated
n.loci          <- 100    # Number of loci underlying the trait
n.f             <- 10   # Number of females
n.m             <- 10   # Number of males
f.RS            <- 2      # Number of offspring per female
f.RS.Pr         <- 1      # Probability of number of offspring per female.
sel.thresh.f    <- 1      # Selection threshold
sel.thresh.m    <- 0.2    # Selection threshold
prdm9.found.maf <- 0.1    # Starting freq of PRDM9


#~~ sample two landscapes and initial frequencies of alleles in founders

r1   <- sample(map.dist/100, n.loci)
r2   <- sample(map.dist/100, n.loci)
rhet <- (r1 + r2)/2 

ggplot(data.frame(r = c(r1, rhet, r2), map = rep(1:3, each = length(r1)), x = rep(1:length(r1), times = 3)),
       aes(x, r, colour = factor(map))) +
  geom_line() +
  scale_colour_brewer(palette = "Set1")

#~~ fix MAFs so that they show a range of frequencies between 0 & 1

mafs.found <- sample(maf.info, n.loci)
mafs.found <- mafs.found + (runif(n.loci) < 0.5)/2
  
#~~ determine PRDM9 genotype frequencies based on HWE
  
p <- prdm9.found.maf
q <- 1 - p
prdm9.found.prs <- c(p^2, 2*p*q, q^2)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 3. Generate information in founder generation                #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# NB Sex 1 = male, 2 = female (Xy, XX)

#~~ generate founder haplotypes

founder.haplos <- lapply(1:n.found.hap, function (x) (runif(n.loci) < mafs.found) + 0L)

#~~ generate diplotypes for n.f females and n.m males

gen.0 <- list()
gen.0[1:(n.f + n.m)] <- list(list(MOTHER = NA, FATHER = NA))

for(i in 1:(n.f + n.m)){
  gen.0[[i]][[1]] <- sample(founder.haplos, size = 1)
  gen.0[[i]][[2]] <- sample(founder.haplos, size = 1)
}

#~~ create reference table

ref.0 <- data.frame(List.Key    = 1:length(gen.0),
                    GEN         = 0,
                    ID          = 1:(n.f + n.m),
                    MOTHER      = NA,
                    FATHER      = NA,
                    SEX         = sapply(1:length(gen.0), function(x) (runif(1) < 0.5) + 1L),
                    PRDM9       = sapply(1:length(gen.0), function(x) sample(1:3, size = 1, prob = prdm9.found.prs)),
                    PHENO       = sapply(1:length(gen.0), function(x) sum(unlist(gen.0[[x]]))))



#~~ create phenotype table for selection

ref.0.pheno <- data.frame(ID.Pheno = tapply(ref.0$Haplo.Pheno, ref.0$ID, sum))
ref.0.pheno$ID <- row.names(ref.0.pheno)
ref.0.pheno <- join(ref.0.pheno, unique(ref.0[,c("ID", "SEX")]))

ref.0.pheno.m <- subset(ref.0.pheno, SEX == 1)                    
ref.0.pheno.m <- ref.0.pheno.m[sample(1:nrow(ref.0.pheno.m), replace = F, size = nrow(ref.0.pheno.m)),]
ref.0.pheno.m <- arrange(ref.0.pheno.m, ID.Pheno)
ref.0.pheno.m <- ref.0.pheno.m[((1-sel.thresh.m)*nrow(ref.0.pheno.m) + 1):nrow(ref.0.pheno.m),]

ref.0.pheno.f <- subset(ref.0.pheno, SEX == 2)                    
ref.0.pheno.f <- ref.0.pheno.f[sample(1:nrow(ref.0.pheno.f), replace = F, size = nrow(ref.0.pheno.f)),]
ref.0.pheno.f <- arrange(ref.0.pheno.f, ID.Pheno)
ref.0.pheno.f <- ref.0.pheno.f[((1-sel.thresh.f)*nrow(ref.0.pheno.f) + 1):nrow(ref.0.pheno.f),]

ref.0 <- subset(ref.0, ID %in% c(ref.0.pheno.m$ID, ref.0.pheno.f$ID))

#~~ generate spaces for diplotypes of two offspring per female and sample best fathers

length.out <- f.RS*nrow(subset(ref.0, SEX == 2))

ref.1 <- data.frame(List.Key    = 1:length.out,
                    ID          = rep(1:(length.out/2), each = 2),
                    MOTHER      = rep(subset(ref.0, SEX == 2)$ID, each = f.RS),
                    FATHER      = rep(sample(subset(ref.0, SEX == 1)$ID,
                                             size = length.out/2,
                                             replace = T),
                                      each = 2),
                    SEX         = rep((runif(1) < 0.5) + 1L,
                                      each = 2,
                                      length.out = length.out),
                    Haplo       = rep(1:2, length.out = length.out),
                    Haplo.Pheno = NA)

#~~ sample fathers





rec.pos <- which(((runif(length(r1)) < r1) + 0L) == 1)

#~~ If there are no recombination events, sample one of the parental haplotypes at random

if(length(rec.pos) == 0){
  haplo.list[as.character(transped$Offspring.ID[j])][[1]][transped$Parent.ID.SEX[j]] <- haplo.list[as.character(transped$Parent.ID[j])][[1]][sample.int(2, 1)]
}

#~~ If there are recombination events, sample the order of the parental haplotypes,
#   exchange haplotypes and then transmit haplotype to offspring

if(length(rec.pos) > 0){
  
  parental.haplotypes <- haplo.list[as.character(transped$Parent.ID[j])][[1]][sample.int(2, 2, replace = F)]
  
  start.pos <- c(1, rec.pos[1:(length(rec.pos))] + 1)
  stop.pos <- c(rec.pos, length(r.female) + 1)
  
  fragments <- list()
  
  for(k in 1:length(start.pos)){
    if(k %% 2 != 0) fragments[[k]] <- parental.haplotypes[1][[1]][start.pos[k]:stop.pos[k]]
    if(k %% 2 == 0) fragments[[k]] <- parental.haplotypes[2][[1]][start.pos[k]:stop.pos[k]]
  }
  
  haplo.list[as.character(transped$Offspring.ID[j])][[1]][transped$Parent.ID.SEX[j]] <- list(unlist(fragments))
  
}
}
}


head(m.0.breedids)

gen.0.m




