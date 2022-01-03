library(tidyverse)
library(reshape2)

setwd("Desktop")
#blastori = read.csv("blast.taxonomy_ori.txt", header = F, sep="\t")
blastori = read.csv("blast.taxonomy.SSU.txt", header = F, sep="\t")
taxembl = read.csv("embl.taxonomy.txt", header = F, sep="\t")

tax_pruned=taxembl[which(taxembl$V1 %in% blastori$V1),]
orileft=blastori[which(!(blastori$V1 %in% tax_pruned$V1)),]

combo = rbind(tax_pruned,orileft)
#write.csv(combo, file="blast.taxonomy.combo.txt", col.names = F, row.names = F)
write.csv(combo, file="blast.taxonomy.SSUcombo.txt", col.names = NULL, row.names = F)
