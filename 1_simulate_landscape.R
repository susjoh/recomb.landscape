
library(ggplot2)
library(plyr)

#~~ Create the sampling distrubutions

map <- read.table("data/soay_map.txt", sep = ",", header = T)

map.dist <- diff(map$cM.Position.Sex.Averaged)
map.dist <- map.dist[which(map.dist >= 0 & map.dist < 2)]

maf.info <- maf$MAF

rm(maf, map)

#~~ define starting values

n.found.hap <- 100
n.loci <- 1000
n.f <- 1000
n.m <- 1000
f.RS <- 2
f.RS.Pr <- 1
sel.thresh <- 0.2

#~~ sample a landscape and initial frequencies of alleles in founders

r1 <- sample(map.dist/100, n.loci)
r2 <- sample(map.dist/100, n.loci)
mafs.found <- sample(maf.info, n.loci)


ggplot(data.frame(r = c(r1, r2), map = rep(1:2, each = length(r1)), x = rep(1:length(r1), times = 2)),
       aes(x, r, colour = factor(map))) +
  geom_line() +
  scale_colour_brewer(palette = "Set1")

#~~ generate founder haplotypes

founder.haplos <- lapply(1:n.found.hap, function (x) (runif(n.loci) < mafs.found) + 0L)

### FOUNDER GENERATION

#~~ generate diplotypes for n.f females and n.m males

gen.0.f <- sample(founder.haplos, size = 2 * n.f, replace = T)
gen.0.m <- sample(founder.haplos, size = 2 * n.m, replace = T)

#~~ selection from males - create a reference frame

m.0.info <- data.frame(ID = 1:n.m,
                       Pheno = unlist(lapply(gen.0.m[seq(1, 2 * n.m, 2)], sum)) +
                         unlist(lapply(gen.0.m[seq(2, 2 * n.m, 2)], sum)))

m.0.info <- m.0.info[sample(1:nrow(m.0.info), replace = F, size = nrow(m.0.info)),]
m.0.info <- arrange(m.0.info, Pheno)

m.0.breedids <- data.frame(ID = rep(sort(m.0.info$ID[((1-sel.thresh)*n.m + 1):(n.m)]), each = 2))
m.0.breedids$Haplo = rep(1:2, length.out = nrow(m.0.breedids))
m.0.breedids$Reference <- (m.0.breedids$ID*2) + (m.0.breedids$Haplo - 2)

#~~ each female has two offspring and is mated with one of the highest males at random


m.1.population


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




