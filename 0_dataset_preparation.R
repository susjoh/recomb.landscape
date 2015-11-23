library(plyr)

#~~ Create the sampling distrubutions

map <- read.table("../Recombination Rate Manuscript/SupplementaryInformation/S1_LinkageMapTable.csv", sep = ",", header = T)
maf <- read.table("../Recombination Rate Manuscript/SupplementaryInformation/S4_GWASresultsASREML.csv", sep = ",", header = T)

maf <- subset(maf, Type == "cistrans" & Model == "All" & Chr != 0)
maf <- subset(maf, select = c(SNP.Name, MAF))
map <- join(map, maf)

write.table(map, "data/soay_map.txt", row.names = F, sep = "\t", quote = F)

system("cmd", input = "git config --global user.name \"susjoh\"")

system("cmd", input = "git remote add origin https://github.com/susjoh/recomb_landscape.git")
system("cmd", input = "git push -u origin master")
