library(tidyverse)
library(reshape2)
library(scales)
library(RColorBrewer)
library(forcats)
library(vegan)
library(adegenet)
library(magrittr)

runid = "RUN2"
barcodes=c("134", "152","166", "17nine", "261", "272")
setwd("~")
#Listing up IDtables and their molten, pruned counterparts
IDLIST=list()
prunedlist=list()
IDLIST_def=list()
prunedlist_def=list()
for(i in 1:length(barcodes)){
  istring = paste("Downloads/",runid, barcodes[i], "IDtable.tsv", sep="")
  istring_def = paste("Downloads/",runid, barcodes[i], "IDtable_default.tsv", sep="")
  IDLIST[[i]] <- read.csv(file = istring, header = TRUE, sep = " ", stringsAsFactors = F)
  IDLIST_def[[i]] <- read.csv(file = istring_def, header = TRUE, sep = " ", stringsAsFactors = F)
  names(IDLIST)[i]<-paste(runid, barcodes[i],sep="")
  names(IDLIST_def)[i]<-paste(runid, barcodes[i],sep="")
  prunedlist[[i]] <- melt(IDLIST[[i]], id.vars = c("Feature_name", "OTU"))
  prunedlist_def[[i]] <- melt(IDLIST_def[[i]], id.vars = c("Feature_name", "OTU"))
  names(prunedlist)[i]<-paste(runid, barcodes[i],sep="")
  names(prunedlist_def)[i]<-paste(runid, barcodes[i],sep="")
  prunedlist[[i]] <- prunedlist[[i]][which(prunedlist[[i]]$value>0 & str_count(prunedlist[[i]]$OTU,"Metazoa")>0 & str_count(prunedlist[[i]]$OTU,"Hydroidolina")==0),]
  prunedlist_def[[i]] <- prunedlist_def[[i]][which(prunedlist_def[[i]]$value>0 & str_count(prunedlist_def[[i]]$OTU,"Metazoa")>0 & str_count(prunedlist_def[[i]]$OTU,"Hydrozoa")==0 ),]
}

comblist = lapply(seq_along(prunedlist), function(x) rbind(prunedlist[[x]], prunedlist_def[[x]]))

#Import reliability
rel_list = list()
for(i in 1:length(barcodes)){
  istring = paste("Downloads/METAXA2_",runid, barcodes[i], ".taxonomy-reliability.txt", sep="")
  istring_def = paste("Downloads/METAXA2_default_",runid, barcodes[i], ".taxonomy-reliability.txt", sep="")
  rel_i <- read.csv(file = istring, header = FALSE, sep = "\t", stringsAsFactors = F)
  rel_i_def <- read.csv(file = istring_def, header = FALSE, sep = "\t", stringsAsFactors = F)
  rel_i_comb <- rbind(rel_i, rel_i_def)
  names(rel_i_comb) <- c("Feature_name", "Relative_IDs")
  rel_list[[i]] <- rel_i_comb
}

#Combine reliability onto taxonomy
comblist_rel = comblist
for(i in 1:length(barcodes)){
  combi = comblist[[i]]
  reli = rel_list[[i]]
  reli_pruned = reli[which(reli$Feature_name %in% combi$Feature_name),]
  comblist_rel[[i]] <- cbind(combi, reli_pruned[match(combi$Feature_name, reli_pruned$Feature_name),"Relative_IDs"])
  names(comblist_rel[[i]])[5]<-"Relative_IDs"
}

#Eyeballing numbers
for(i in 1:length(barcodes)){
  print(barcodes[i])
  idtable = comblist_rel[[i]]
  idtable[which(idtable$Relative_ID %>% str_count("Isopo")>0),1:2] %>% print()
}

i = 1
plotlist2 = list()
for(t in comblist){
  #pdf(paste("Desktop/",barcodes[i],"_OTUplot.pdf",sep=""), width = 10, height = 10)
  plotlist2[[i]] <- ggplot(t, aes(x = variable, y = log(value), fill = OTU)) + geom_bar(position = "fill",stat = "identity") + theme(legend.key.size = unit(1,"line"), legend.text = element_text(size=4),legend.position="bottom", axis.text.x = element_text(angle = 90, hjust = 1)) + scale_fill_manual(values=rep(c("#89C5DA", "#DA5724", "#74D944", "#CE50CA", "#3F4921", "#C0717C", "#CBD588", "#5F7FC7", "#673770", "#D3D93E", "#38333E", "#508578", "#D7C1B1", "#689030", "#AD6F3B", "#CD9BCD", "#D14285", "#6DDE88", "#652926", "#7FDCC0", "#C84248", "#8569D5", "#5E738F", "#D1A33D", "#8A7C64", "#599861"),10)) + ggtitle(names(prunedlist)[i]) + guides(shape = guide_legend(override.aes = list(size = 1)))
  plotlist2[[i]] %>% print()
  #dev.off()
  i=i+1
}

combinetab = rbind(comblist_rel[[1]],comblist_rel[[2]],comblist_rel[[3]],comblist_rel[[4]],comblist_rel[[5]],comblist_rel[[6]])
combinetab$variable <- str_remove(combinetab$variable,"^X\\w\\w\\w_")
combinetab$OTU <-  str_remove(combinetab$OTU,"__.+$")
combinetab$OTU <-  str_remove(combinetab$OTU,"_NA.+$")
#combinetab$OTU <-  str_remove(combinetab$OTU,"uncultured eukaryote$")

#Eyeballing numbers
combinetab[which(combinetab$variable %>% str_count("124_Nano")>0), c(1,5)]


#ggplot(combinetab, aes(x = variable, y = log(value), fill = OTU)) + geom_bar(position = "fill",stat = "identity") + theme(legend.key.size = unit(1,"line"), legend.text = element_text(size=5),legend.position="bottom", axis.text.x = element_text(angle = 90, hjust = 1)) + scale_fill_manual(values=rep(c("#89C5DA", "#DA5724", "#74D944", "#CE50CA", "#3F4921", "#C0717C", "#CBD588", "#5F7FC7", "#673770", "#D3D93E", "#38333E", "#508578", "#D7C1B1", "#689030", "#AD6F3B", "#CD9BCD", "#D14285", "#6DDE88", "#652926", "#7FDCC0", "#C84248", "#8569D5", "#5E738F", "#D1A33D", "#8A7C64", "#599861"),10)) + ggtitle(names(prunedlist)[i]) + guides(shape = guide_legend(override.aes = list(size = 1)))
#Color alternative
combinetab[which(combinetab[,2] == unique(combinetab[,2] == unique(combinetab[,2]) | combinetab[,3] == unique(combinetab[,3]))),]
ggplot(combinetab, aes(x = variable, y = log(value), fill = Relative_IDs)) + geom_bar(position = "fill",stat = "identity") + theme(legend.key.size = unit(0.5,"line"), legend.text = element_text(size=3),legend.position="bottom", axis.text.x = element_text(angle = 90, hjust = 1)) + scale_fill_manual(values=c("#89C5DA", "#DA5724", "#74D944", "#CE50CA", "#3F4921", "#C0717C", "#CBD588", "#5F7FC7", "#673770", "#D3D93E", "#38333E", "#508578", "#D7C1B1", "#689030", "#AD6F3B", "#CD9BCD", "#D14285", "#6DDE88", "#652926", "#7FDCC0", "#C84248", "#8569D5", "#5E738F", "#D1A33D", "#8A7C64", "#599861",colors()[sample(length(colors()),length(colors())/4)])) + ggtitle("All Barcodes") + guides(shape = guide_legend(override.aes = list(size = 1)))
