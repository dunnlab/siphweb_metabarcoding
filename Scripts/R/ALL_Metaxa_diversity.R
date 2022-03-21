library(tidyverse)
library(reshape2)
library(scales)
library(RColorBrewer)
library(forcats)
library(vegan)
library(adegenet)
library(magrittr)
library(plyr)
library(paleotree)
library(phytools)
library(patchwork)
library(ggtree)
library(ggplotify)
library(scatterpie)
library(googledrive)

barcodes=c("134", "152","166", "179", "261", "272")
drive_get("/")
setwd("~/siphweb_metabarcoding/Data/IDtables_ALL/")

runids = c("RUN0","RUN1","RUN2","RUN3", "RUN4","RUN5")
metacomblist=list()
metarel_list=list()
for(run in 1:length(runids)){
  runid=runids[run]
  IDLIST=list()
  prunedlist=list()
  IDLIST_def=list()
  prunedlist_def=list()
  rel_list = list()
  if(runid=="RUN2"){
    barcodes=c("134", "152","166", "17nine", "261", "272")
  }
  else barcodes = c("134", "152","166", "179", "261", "272")
  for(i in 1:length(barcodes)){
    istring = paste(runid, barcodes[i], "IDtable.tsv", sep="")
    istring_def = paste(runid, barcodes[i], "IDtable_splus.tsv", sep="")
    IDLIST[[i]] <- read.csv(file = istring, header = TRUE, sep = " ", stringsAsFactors = F)
    IDLIST_def[[i]] <- read.csv(file = istring_def, header = TRUE, sep = " ", stringsAsFactors = F)
    names(IDLIST)[i]<-paste(runid, barcodes[i],sep="")
    names(IDLIST_def)[i]<-paste(runid, barcodes[i],sep="")
    prunedlist[[i]] <- cbind(melt(IDLIST[[i]], id.vars = c("Feature_name", "OTU")), rep(runids[run], nrow(melt(IDLIST[[i]], id.vars = c("Feature_name", "OTU")))), rep(barcodes[i], nrow(melt(IDLIST[[i]], id.vars = c("Feature_name", "OTU")))), rep("SILVA", nrow(melt(IDLIST[[i]], id.vars = c("Feature_name", "OTU")))))
    names(prunedlist[[i]]) <- c("Feature_name", "OTU","sample", "abundance","run","barcode","database")
    prunedlist[[i]]$sample <- str_replace_all(prunedlist[[i]]$sample, paste("X",barcodes[i],"_", sep = ""),"")
    if(runid=="RUN0"){
      prunedlist[[i]]$sample <- str_replace_all(prunedlist[[i]]$sample, "_[CGTA]+.[CGTA]+_L001","")
    }
    else prunedlist[[i]]$sample <- str_replace_all(prunedlist[[i]]$sample, "_\\w\\w\\w_\\w\\w\\w_S\\w+_L001","")
    def_melt <- melt(IDLIST_def[[i]], id.vars = c("Feature_name", "OTU"))
    def_melt_nrow <- nrow(def_melt)
    prunedlist_def[[i]] <- cbind(def_melt, rep(runids[run],def_melt_nrow), rep(barcodes[i], def_melt_nrow), rep("SILVA138plus", def_melt_nrow))
    names(prunedlist_def[[i]]) <- c("Feature_name", "OTU","sample", "abundance","run","barcode","database")
    prunedlist_def[[i]]$sample <- str_replace_all(prunedlist_def[[i]]$sample, paste("X",barcodes[i],"_", sep = ""),"")
    if(runid=="RUN0"){
      prunedlist_def[[i]]$sample <- str_replace_all(prunedlist_def[[i]]$sample, "_[CGTA]+.[CGTA]+_L001","")
    }
    else prunedlist_def[[i]]$sample <- str_replace_all(prunedlist_def[[i]]$sample, "_\\w\\w\\w_\\w\\w\\w_S\\w+_L001","")
    names(prunedlist)[i]<-paste(runid, barcodes[i],sep="")
    names(prunedlist_def)[i]<-paste(runid, barcodes[i],sep="")
    #prunedlist[[i]] <- prunedlist[[i]][which(prunedlist[[i]]$abundance>0 & str_count(prunedlist[[i]]$OTU,"Metazoa")>0 & str_count(prunedlist[[i]]$OTU,"Hydroidolina")==0 & str_count(prunedlist[[i]]$OTU,"Mammalia")==0 & str_count(prunedlist[[i]]$OTU,"Metazoa__")==0 & str_count(prunedlist[[i]]$OTU,"Tubificin")==0 & str_count(prunedlist[[i]]$OTU,"Catenulid")==0 & str_count(prunedlist[[i]]$OTU,"Pentastomida")==0 & str_count(prunedlist[[i]]$OTU,"Myxozoa")==0 & str_count(prunedlist[[i]]$OTU,"Insect")==0),]
    prunedlist[[i]] <- prunedlist[[i]][which(prunedlist[[i]]$abundance>0),]
    #prunedlist_def[[i]] <- prunedlist_def[[i]][which(prunedlist_def[[i]]$abundance>0 & str_count(prunedlist_def[[i]]$OTU,"Hydroidolina")==0 & str_count(prunedlist_def[[i]]$OTU,"Mammalia")==0 & str_count(prunedlist_def[[i]]$OTU,"Tubificin")==0 & str_count(prunedlist_def[[i]]$OTU,"Catenulid")==0 & str_count(prunedlist_def[[i]]$OTU,"Pentastomida")==0 & str_count(prunedlist_def[[i]]$OTU,"Myxozoa")==0 & str_count(prunedlist_def[[i]]$OTU,"Insect")==0),]
    prunedlist_def[[i]] <- prunedlist_def[[i]][which(prunedlist_def[[i]]$abundance>0),]
    irel = paste("METAXA2_",runid, barcodes[i], ".taxonomy-reliability.txt", sep="")
    irel_def = paste("METAXA2_plus_",runid, barcodes[i], ".taxonomy-reliability.txt", sep="")
    rel_i <- read.csv(file = irel, header = FALSE, sep = "\t", stringsAsFactors = F)
    rel_i <- cbind(rel_i, rep(runids[run],nrow(rel_i)), rep(barcodes[i], nrow(rel_i)), rep("SILVA", nrow(rel_i)))
    if (!file.size(irel_def) == 0) {
      rel_i_def <- read.csv(file = irel_def, header = FALSE, sep = "\t", stringsAsFactors = F)
      rel_i_def <- cbind(rel_i_def, rep(runids[run],nrow(rel_i_def)), rep(barcodes[i], nrow(rel_i_def)), rep("SILVA138plus", nrow(rel_i_def)))
      names(rel_i_def) <- c("Feature_name", "Relative_IDs", "run","barcode","database")
    }
    else rel_i_def <- data.frame()
    names(rel_i) <- c("Feature_name", "Relative_IDs", "run","barcode","database")
    rel_i_comb <- rbind(rel_i, rel_i_def)
    rel_list[[i]] <- rel_i_comb
  }
  comblist = lapply(seq_along(prunedlist), function(x){rbind(prunedlist[[x]], prunedlist_def[[x]])})
  metacomblist[[run]] <- comblist
  metarel_list[[run]] <- rel_list
}
names(metacomblist) <- runids
names(metarel_list) <- runids

## HERE the issues start ##

#Integrating reliability and taxonomy tables
metacomblist_rel=list()
for(run in 1:length(runids)){
  #Combine reliability onto OTUs
  comblist = metacomblist[[run]]
  comblist_rel = metacomblist[[run]]
  rel_list <- metarel_list[[run]]
  for(i in 1:length(barcodes)){
    combi = comblist[[i]]
    combi_132 = combi[which(combi$database == "SILVA"),]
    combi_plus = combi[which(combi$database == "SILVA138plus"),]
    reli = rel_list[[i]]
    reli_132 = reli[which(reli$database == "SILVA"),]
    reli_plus = reli[which(reli$database == "SILVA138plus"),] 
    reli_132_pruned = reli_132[which(reli_132$Feature_name %in% combi_132$Feature_name),]
    reli_plus_pruned = reli_plus[which(reli_plus$Feature_name %in% combi_plus$Feature_name),]
    merge_132 <- cbind(combi_132, reli_132_pruned[match(combi_132$Feature_name, reli_132_pruned$Feature_name),"Relative_IDs"])
    names(merge_132)[8] <- "Relative_IDs"
    merge_plus <- cbind(combi_plus, reli_plus_pruned[match(combi_plus$Feature_name, reli_plus_pruned$Feature_name),"Relative_IDs"])
    names(merge_plus)[8] <- "Relative_IDs"
    comblist_rel[[i]] <- rbind(merge_132,merge_plus) 
    names(comblist_rel[[i]])[8]<-"Relative_IDs"
  }
  metacomblist_rel[[run]] <- comblist_rel
}

#Import taxonomy
metatablist=list()
for(run in 1:length(runids)){
  tab_list = list()
  runid=runids[run]
  if(runid=="RUN2"){
    barcodes=c("134", "152","166", "17nine", "261", "272")
  }
  else barcodes = c("134", "152","166", "179", "261", "272")
  for(i in 1:length(barcodes)){
    istring_def = paste("METAXA2_plus_",runid, barcodes[i], ".taxonomy.txt", sep="")
    istring = paste("METAXA2_",runid, barcodes[i], ".taxonomy-table.tsv", sep="")
    tab_i <- read.csv(file = istring, header = FALSE, sep = ",", stringsAsFactors = F)
    if (!file.size(istring_def) == 0) {
    print("custom entries present")
    tab_i_def <- read.csv(file = istring_def, header = FALSE, sep = "\t", stringsAsFactors = F)[,1:3]
    }
    else tab_i_def <- data.frame()
    if(run==1 | run>4){
      print("processing RUN0 or RUN4/5")
      istring = paste("METAXA2_",runid, barcodes[i], ".taxonomy.txt", sep="")
      tab_i <- read.csv(file = istring, header = FALSE, sep = "\t", stringsAsFactors = F)[,1:3]
      istring_def = paste("METAXA2_plus_",runid, barcodes[i], ".taxonomy.txt", sep="")
      if (!file.size(istring_def) == 0) {
      print("custom entries present in RUN0/RUN4")
      tab_i_def <- read.csv(file = istring_def, header = FALSE, sep = "\t", stringsAsFactors = F)[,1:3]
      }
      else tab_i_def <- data.frame()
    }
    tab_i <- cbind(tab_i, rep(runids[run],nrow(tab_i)), rep(barcodes[i], nrow(tab_i)), rep("SILVA", nrow(tab_i)))
    names(tab_i) <- c("Feature_name", "IDs", "scores", "run","barcode","database")
    if (nrow(tab_i_def) > 0) {
    print("custom entries parsed")
    tab_i_def <- cbind(tab_i_def, rep(runids[run],nrow(tab_i_def)), rep(barcodes[i], nrow(tab_i_def)), rep("SILVA138plus", nrow(tab_i_def)))
    names(tab_i_def) <- c("Feature_name", "IDs", "scores", "run","barcode","database")
    tab_i_comb <- rbind(tab_i[,1:6], tab_i_def)
    }
    else tab_i_comb <- tab_i
    tab_i_comb$Feature_name <- str_remove_all(tab_i_comb$Feature_name,"\001")
    tab_list[[i]] <- tab_i_comb
  }
  metatablist[[run]] <- tab_list
}

#Combine taxonomy onto OTUs
metacomblist_tab=list()
for(run in 1:length(runids)){
  comblist_rel = metacomblist_rel[[run]]
  comblist_tab = metacomblist_rel[[run]]
  tab_list = metatablist[[run]]
  for(i in 1:length(barcodes)){
    combi = comblist_rel[[i]]
    combi_132 = combi[which(combi$database == "SILVA"),]
    combi_plus = combi[which(combi$database == "SILVA138plus"),]
    tabi = tab_list[[i]]
    tabi_132 = tabi[which(tabi$database == "SILVA"),]
    tabi_132_pruned <- tabi_132[which(tabi_132$Feature_name %in% combi_132$Feature_name),]
    tabi_plus = tabi[which(tabi$database == "SILVA138plus"),] 
    tabi_plus_pruned <- tabi_plus[which(tabi_plus$Feature_name %in% combi_plus$Feature_name),]
    merge_132 <- cbind(combi_132, tabi_132_pruned[match(combi_132$Feature_name, tabi_132_pruned$Feature_name),"IDs"])
    names(merge_132)[9] <- "IDs"
    merge_plus <- cbind(combi_plus, tabi_plus_pruned[match(combi_plus$Feature_name, tabi_plus_pruned$Feature_name),"IDs"])
    names(merge_plus)[9] <- "IDs"
    comblist_tab[[i]] <- rbind(merge_132,merge_plus)
    names(comblist_tab[[i]])[9]<-"IDs"
  }
  metacomblist_tab[[run]] <- comblist_tab
}

metacomblist <- reverseList(metacomblist)
metacomblist_rel <- reverseList(metacomblist_rel)
metacomblist_tab <- reverseList(metacomblist_tab)
metarel_list <- reverseList(metarel_list)
metatablist <- reverseList(metatablist)

comblist <- lapply(metacomblist, function(L){do.call(rbind, L)})
comblist_rel <- lapply(metacomblist_rel, function(L){do.call(rbind, L)})
comblist_tab <- lapply(metacomblist_tab, function(L){do.call(rbind, L)})
rel_list <- lapply(metarel_list, function(L){do.call(rbind, L)})
tablist <- lapply(metatablist, function(L){do.call(rbind, L)})

combinetab = rbind(comblist_tab[[1]],comblist_tab[[2]],comblist_tab[[3]],comblist_tab[[4]],comblist_tab[[5]],comblist_tab[[6]])
combinetab$OTU <-  str_remove(combinetab$OTU,"__.+$")
combinetab$OTU <-  str_remove(combinetab$OTU,"_NA.+$")

#The END#

combinetab[which(combinetab$database=="SILVA138plus"),] -> splus
splus <- splus[which(!(splus$sample %in% c("21_mix1","22_mix2","23_mix3","26_Forsk"))),]
#write.table(splus, "SILVA138plus_new.tsv", sep="\t", row.names = F, col.names = T)

combinetab[which(combinetab$database=="Custom" & combinetab$run == "RUN2" & combinetab$barcode == "17nine"),] -> RUN2_17nine
#write.table(RUN2_17nine, "customRUN2_17nine.tsv", sep="\t", row.names = F, col.names = T)

combinetab[which(combinetab$run == "RUN5" & combinetab$database=="SILVA"),] -> RUN5
#write.table(RUN5, "RUN5_SILVA.tsv", sep="\t", row.names = F, col.names = T)

combinetab[which(combinetab$database=="Custom" & str_count(combinetab$IDs,"Metazoa")>0 ),] -> customtab
#write.table(customtab, "customDBparsed.tsv", sep="\t", row.names = F, col.names = T)

combinetab[which(!(combinetab$Feature_name %in% allsamples$Feature_name)),] -> nonmetazoan
nonmetazoan <- nonmetazoan[which(str_count(nonmetazoan$sample,"mix")==0 ),]
#write.table(nonmetazoan, "nonmetazoan.tsv", sep="\t", row.names = F, col.names = T)

ggplot(customtab, aes(x = sample, y = log(abundance), fill = OTU)) + geom_bar(position = "fill",stat = "identity") + theme(legend.key.size = unit(1,"line"), legend.text = element_text(size=5),legend.position="bottom", axis.text.x = element_text(angle = 90, hjust = 1)) + scale_fill_manual(values=rep(c("#89C5DA", "#DA5724", "#74D944", "#CE50CA", "#3F4921", "#C0717C", "#CBD588", "#5F7FC7", "#673770", "#D3D93E", "#38333E", "#508578", "#D7C1B1", "#689030", "#AD6F3B", "#CD9BCD", "#D14285", "#6DDE88", "#652926", "#7FDCC0", "#C84248", "#8569D5", "#5E738F", "#D1A33D", "#8A7C64", "#599861"),10)) + ggtitle("Custom") + guides(shape = guide_legend(override.aes = list(size = 1)))

### Analyze Metaxa & Trawl Data
#Read metabarcoding data
setwd("~/siphweb_metabarcoding")
allsamples <- read.csv("Manual_curation/AllSamplesParsed.tsv", sep="\t", header = T, stringsAsFactors = F)
allsamples <- allsamples[,-66]
#Parse sequences/feature_names
allseqs <- read.csv("Data/seqs-feats.csv", sep="\t", header = F, stringsAsFactors = F)
names(allseqs) <- c("Feature_name","Sequence")
allseqs <- allseqs[which(allseqs$Feature_name %in% unique(allsamples$Feature_name)),]
seqs_samples <- full_join(allsamples,allseqs, by="Feature_name")
#write.table(seqs_samples, "SequencesParsed.tsv", sep="\t", col.names = T, row.names = F)

#Read trawl data
trawls = read.csv("Data/trawls_data_Nov2020.tsv", header=T, sep="\t", stringsAsFactors = F)
trawls$Broad.group[which(trawls$Broad.group == "Krill")] <- "Euphausid"
trawls$Broad.group[which(trawls$Broad.group == "Gymnosome pteropod")] <- "Gastropod"
trawls$Broad.group[which(trawls$Broad.group == "Appendicularia")] <- "Larvacean"
trawls$Broad.group[which(trawls$Broad.group == "Appendicularian")] <- "Larvacean"
trawls$Broad.group[which(trawls$Broad.group == "Benthic snail")] <- "Gastropod"
trawls$Broad.group[which(trawls$Broad.group == "Copepod ")] <- "Copepod"
samples <- split(trawls, trawls$Sample)
#samples.bygroup <- lapply(samples, function(x){aggregate(x$Representative.count, by = list(x$Broad.group), FUN = sum)})
samples.bygroup <- lapply(samples, function(x){aggregate(x$Count, by = list(x$Broad.group), FUN = sum)})
for(i in 1:length(samples.bygroup)){
  samples.bygroup[[i]] <- cbind(rep(names(samples.bygroup)[i], nrow(samples.bygroup[[i]])), samples.bygroup[[i]])
}
broad.taxa <- bind_rows(samples.bygroup)
names(broad.taxa) <- c("Sample", "Broad.group", "Representative.count")
broad.taxa$Representative.count[is.na(broad.taxa$Representative.count)]<-1

#Encode extraction-to-trawl match
ext <- read.csv("Data/extractions_GC_Nov2020.tsv", header=T, stringsAsFactors = F, sep="\t" )
trawl_to_specimen <- ext[,c(1,3,13)]
RUN0trawls = data.frame(Extraction..=allsamples$Extraction, Specimen.. = allsamples$Specimen, Corrresponding.prey.field.sample = rep(NA, nrow(allsamples)))[which(allsamples$run == "RUN0"),]
names(RUN0trawls) <- names(trawl_to_specimen)
trawl_to_specimen <- rbind(trawl_to_specimen, RUN0trawls)

dyads_by_specimen = list()
it=0
for(i in unique(allsamples$Specimen)){
  dyad_i = list()
  dyad_i[[1]] <- allsamples[which(allsamples$Specimen == i),]
  tows <- trawl_to_specimen$Corresponding.prey.field.sample[which(trawl_to_specimen$Specimen.. == i)]
  trawls_i <- trawls[which(trawls$Sample == tows),]
  if(nrow(trawls_i)>0){dyad_i[[2]] <- trawls_i}
  else dyad_i[[2]] <- NA
  names(dyad_i)=c("Diet", "Preyfield")
  it=it+1
  dyads_by_specimen[[it]] <- dyad_i
  names(dyads_by_specimen)[[it]] <- i
}

dyads_by_specimen$`BIOS19-D2-P1`$Preyfield[,c(2,3,4,7)]
dyads_by_specimen$`BIOS19-D2-P1`$Diet[,c(8:11,48)]

## FINAL FIGURE ###
#row by specimen - #broad GC group set of columns - #broad prey groups in trawl set of columns

#Widen format of GC data
diets_broad <- allsamples[which(allsamples$Interpretation == "Prey"),c("Specimen","Species","Broad.group")] %>% unique()
cast_broadiets <- dcast(diets_broad, Specimen + Species ~ Broad.group)
cast_broadiets[is.na(cast_broadiets)] <- 0
cast_broadiets[,c(3:ncol(cast_broadiets))] <- sapply(cast_broadiets[,c(3:ncol(cast_broadiets))], as.numeric)
cast_broadiets[is.na(cast_broadiets)] <- 1

#Widen format of trawl data
cast_trawl <- dcast(broad.taxa, Sample ~ Broad.group)
cast_trawl[is.na(cast_trawl)] <- 0

#Match trawlcast to specimen
trawlbyspecimen <- cast_trawl[match(trawl_to_specimen$Corresponding.prey.field.sample[match(cast_broadiets$Specimen, trawl_to_specimen$Specimen..)], cast_trawl$Sample),]
names(trawlbyspecimen)[1]<-"Trawl"
trawlbyspecimen <- data.frame(trawlbyspecimen, Phyllocarid = trawlbyspecimen$Stomatopod*0) #add phyllocarid zeroes

#Bind GC and trawls
GC_Trawl_broad <- cbind(cast_broadiets, trawlbyspecimen)
rownames(GC_Trawl_broad)<-GC_Trawl_broad$Specimen
GC_Trawl_broad$Species[which(GC_Trawl_broad$Species=="Apolemia undescribed sp")] <- "Apolemia sp"
GC_Trawl_broad$Species[which(GC_Trawl_broad$Species=="Nanomia (deep)")] <- "Nanomia sp. deep"
GC_Trawl_broad[,c("Cephalopod", "Cumacean")]

#Phylogeny of the species sampled
tree = read.newick(text="(Physalia physalis,((Apolemia sp,Apolemia rubriversa,Apolemia lanosa),((Bargmannia amoena,Bargmannia elongata,Bargmannia lata),(((Undescribed physonect L),(Resomia dunni,(Forskalia sp.,(Lychnagalma utricularia,(Halistemma rubrum,(Nanomia sp. Atlantic,Nanomia sp. shallow,Nanomia sp. deep)))))),(Vogtia serrata,(Desmophyes annectens,(Chuniphyes multidentata,((Sphaeronectes koellikeri,Sphaeronectes christiansonae),((Muggiaea atlantica,(Lensia conoidea,Sulculeolaria chuni)),Diphyes dispar)))))))))));")

rotateNodes<-function(tree,nodes,polytom=c(1,2,3),...){
  n<-length(tree$tip.label)
  if(nodes[1]=="all") nodes<-1:tree$Nnode+n
  for(i in 1:length(nodes))
    tree<-ape::rotate(tree,nodes[i],polytom)
  if(hasArg(reversible)) reversible<-list(...)$reversible
  else reversible<-TRUE
  if(reversible){
    ii<-which(tree$edge[,2]<=n)
    jj<-tree$edge[ii,2]
    tree$edge[ii,2]<-1:n
    tree$tip.label<-tree$tip.label[jj]
  }
  return(tree)
}

tree2 <- tree
#tree2 <- rotateNodes(tree,nodes="all")
plotTree(tree2)
is_tip <- tree2$edge[,2] <= length(tree2$tip.label)
ordered_tips <- tree2$edge[is_tip, 2]
SPorder <- tree2$tip.label[ordered_tips]

#plot GC
meltedGC = melt(GC_Trawl_broad[,1:16])
names(meltedGC)[3:4] <- c("Gut Content Prey", "Presence")
meltedGC$Species <- factor(meltedGC$Species, levels=SPorder)
IDorder <- unique(arrange(meltedGC, Species)$Specimen)
meltedGC$Specimen <- factor(meltedGC$Specimen, levels=IDorder)
PreyOrder <- c("Copepod", "Decapod", "Euphausid", "Mysid", "Lophogastrid","Stomatopod", "Amphipod", "Ostracod", "Bivalve", "Gastropod", "Larvacean", "Salp", "Scyphomedusa", "Ctenophore", "Actinopteri")
meltedGC$`Gut Content Prey` <- factor(meltedGC$`Gut Content Prey`, levels=PreyOrder)

GC <- ggplot(meltedGC, aes(x = `Gut Content Prey`, y = Specimen, fill = Presence)) + geom_tile(color="grey") +scale_fill_gradient(low = "white", high = "black") + theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = "none") + ylab("Specimen")

## by species
GCspecies <- aggregate(x = GC_Trawl_broad[,3:16], by = list(GC_Trawl_broad$Species), FUN = "sum")
names(GCspecies)[1] <- "Species"
GCspecies_melt = melt(GCspecies)
names(GCspecies_melt)[2:3] <- c("Gut Content Prey","N")
GCspecies_melt$Species <- factor(GCspecies_melt$Species, levels=SPorder)
GCspecies_melt$`Gut Content Prey` <- factor(GCspecies_melt$`Gut Content Prey`, levels=PreyOrder)
GC_spp <- ggplot(GCspecies_melt, aes(x = `Gut Content Prey`, y = Species, fill = log(N+1))) + geom_tile(color="grey") + scale_fill_gradient(low = "white", high = "black") + geom_text(aes(label=N), col="red") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ylab("Species")
##

#plot trawls
meltedTrawls = melt(GC_Trawl_broad[,c(1,2,17:ncol(GC_Trawl_broad))])
names(meltedTrawls)[4:5] <- c("Prey", "Abundance")
meltedTrawls <- ddply(meltedTrawls, .(Specimen), transform, Relative.Abundance = scale(log(Abundance+1)))
meltedTrawls <- meltedTrawls[which(!(meltedTrawls$Prey %in% c("Animal", "Egg"))),]
names(meltedTrawls)[4]<-"Ambient Prey"
PreyOrderTr <- c("Copepod", "Decapod", "Euphausid", "Mysid", "Lophogastrid", "Stomatopod",  "Ostracod", "Bivalve","Larvacean", "Salp",  "Scyphomedusa", "Ctenophore", "Actinopteri", "Cladoceran", "Amphipod", "Isopod", "Barnacle", "Cumacean", "Phyllocarid", "Brachiopod", "Echinoderm","Polychaete", "Chaetognath", "Gastropod", "Cephalopod", "Doliolid", "Pyrosome", "Hydromedusa", "Siphonophore", "Anthozoa")

meltedTrawls$`Ambient Prey` <- factor(meltedTrawls$`Ambient Prey`, levels=PreyOrderTr)
meltedTrawls$Species <- factor(meltedTrawls$Species, levels=SPorder)
IDorder <- unique(arrange(meltedTrawls, Species)$Specimen)
meltedTrawls$Specimen <- factor(meltedTrawls$Specimen, levels=IDorder)

TR <- ggplot(meltedTrawls, aes(x = `Ambient Prey`, y = Specimen, fill = Relative.Abundance)) + geom_tile() + scale_fill_gradient(low = "white", high = "dark blue", na.value = "grey") + theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), axis.title.y=element_blank(),axis.text.y=element_blank(),legend.position = "none", axis.ticks = element_blank())

#Selectivity?
LI <- function(feeding, overlap, output){
  lapply(as.list(1:nrow(overlap)), function(X){
    lapply(as.list(1:ncol(feeding)), function(Y){
      output[X,Y] <<- feeding[X,Y] - overlap[X,Y]
    })
  })
  return(output)
}


rownames(trawlbyspecimen) = cast_broadiets$Specimen
rownames(cast_broadiets) = cast_broadiets$Specimen

rel_trawls <- t(apply(trawlbyspecimen[,c(-1,-4,-18)],1,function(x){x <- x/sum(x, na.rm = T); return(x)}))  #-4 -18 to get rid of "Animal" and "Egg" concepts

colnames(rel_trawls)[which(!(colnames(rel_trawls) %in% colnames(cast_broadiets)))] #CHECK EVERY TIME!

expanded_broadiets <- data.frame(cast_broadiets, Amphipod = rep(0, nrow(cast_broadiets)),Anthozoa = rep(0, nrow(cast_broadiets)),Barnacle = rep(0, nrow(cast_broadiets)),Brachiopod = rep(0, nrow(cast_broadiets)),Cephalopod = rep(0, nrow(cast_broadiets)),Chaetognath = rep(0, nrow(cast_broadiets)),Cladoceran = rep(0, nrow(cast_broadiets)),Cumacean = rep(0, nrow(cast_broadiets)),Doliolid = rep(0, nrow(cast_broadiets)),Echinoderm = rep(0, nrow(cast_broadiets)), Gastropod = rep(0, nrow(cast_broadiets)), Pyrosome = rep(0, nrow(cast_broadiets)),Siphonophore = rep(0, nrow(cast_broadiets)), Hydromedusa = rep(0, nrow(cast_broadiets)), Isopod = rep(0, nrow(cast_broadiets)), Polychaete = rep(0, nrow(cast_broadiets)), Phyllocarid = rep(0, nrow(cast_broadiets)))

expanded_broadiets <- expanded_broadiets[,match(colnames(rel_trawls), colnames(expanded_broadiets))]

rel_GC <- t(apply(expanded_broadiets,1,function(x){x <- x/sum(x); return(x)}))

selectivity <- as.data.frame(matrix(nrow=nrow(rel_trawls), ncol=ncol(rel_trawls)))
names(selectivity) <- colnames(rel_trawls); rownames(selectivity) <- rownames(rel_GC)
selectivity <- LI(rel_GC, rel_trawls,selectivity)

cast_broadiets$Species[which(cast_broadiets$Species=="Apolemia undescribed sp")] <- "Apolemia sp"
cast_broadiets$Species[which(cast_broadiets$Species=="Nanomia (deep)")] <- "Nanomia sp. deep"
selectivity <- data.frame(selectivity, Species = cast_broadiets$Species, Specimen = cast_broadiets$Specimen)
selectivity_melt <- melt(selectivity, id.vars = c("Specimen", "Species"))
names(selectivity_melt) <- c("Specimen","Species", "Prey", "L.I.")

selectivity_melt$Species <- as.character(selectivity_melt$Species)
selectivity_melt$Species <- factor(selectivity_melt$Species, levels=SPorder[which(SPorder %in% selectivity_melt$Species)])
IDorder <- unique(arrange(selectivity_melt, Species)$Specimen)
selectivity_melt$Specimen <- factor(selectivity_melt$Specimen, levels=IDorder)
selectivity_melt$Prey <- factor(selectivity_melt$Prey, levels=PreyOrderTr)
names(selectivity_melt)[3]<-"Selectivity"
selectivity_melt[,4] <- as.numeric(selectivity_melt[,4])

SEL <- ggplot(selectivity_melt, aes(x = Selectivity, y = Specimen, fill = `L.I.`)) + geom_tile() + scale_fill_gradient2(low = "red", high = "blue", mid = "white", na.value = "grey") + theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), axis.title.y=element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank())

## by species selectivity
selectivity_spp <- aggregate(selectivity_melt$L.I.~selectivity_melt$Selectivity+selectivity_melt$Species, FUN=function(x){mean(x,na.rm = T)})
names(selectivity_spp) <- c("Prey","Species","L.I.")
GCspecies_melt$Species[which(!(GCspecies_melt$Species %in% selectivity_spp$Species))] %>% unique() -> missing_spp
expand.grid(unique(selectivity_spp$Prey),missing_spp)->missing_entries
missing_entries <- cbind(missing_entries, rep(NA, nrow(missing_entries)))
names(missing_entries) <- names(selectivity_spp)
selectivity_spp <- rbind(selectivity_spp, missing_entries)

SEL_spp <- ggplot(selectivity_spp, aes(x = Prey, y = Species, fill = `L.I.`)) + geom_tile() + scale_fill_gradient2(low = "red", high = "blue", mid = "white", na.value = "grey") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), axis.title.y=element_blank())

## 

#plot combined
pdf("GC-Trawls-SEL.pdf", height = 10, width = 12)
wrap_plots(GC,TR,SEL, widths = c(1,1.7,1.7))
dev.off()

pdf("GC-SEL_spp.pdf", height = 10, width = 12)
wrap_plots(GC_spp,SEL_spp, widths = c(1,1.7))
dev.off()

pdf("PhyloTree.pdf", height = 10, width = 12)
plotTree(tree2)
dev.off()

#Specimen tree
data.frame(GC_Trawl_broad$Species,GC_Trawl_broad$Specimen) %>% arrange(GC_Trawl_broad.Species)

ind_tree = read.newick(text="((P1,P2,P3),(((D10,SS4,D2),(D1,D2),(D1,D6)),(((SS12,D12),SS6,(SS7,D9)),(((D5,(D6,SS8)),(D8,((D9,D9),(SS11,(D8,((36,38,24,26,28,30,34),(40.1,BW7-4),(SS8,SS10))))))),((D8,SS10),(SS12,(SS9,(((35.2,BW60),D6),((BW65,(SS12,2_4)),(39,BW25,BW1,BW24))))))))))));")

pdf("Specimen_Figure/indTree.pdf", height = 10, width = 5)
plot(ind_tree %>% rotateNodes("all"))
dev.off()

ggtree(ind_tree %>% rotateNodes("all"))+ geom_tiplab()
#plot combined
#pdf("GC-Trawls-SEL.pdf", height = 10, width = 12)

#plot(ind_tree %>% rotateNodes("all"))
wrap_plots(GC,TR,SEL, widths = c(1,1,1.7,1.7))
#dev.off()

#### COMBINED PLOT WITH LITERATURE+VARS+DAPC_predictions
lit <- read.csv("Data/Literature-data/allinteractions.tsv", sep="\t", stringsAsFactors = F)

ROV_lit <- lit[which(lit$source %in% c("Pugh2009cResomiidae", "Choy_2017", "Hissman2005_newrho")),]
ROV_lit$Predator[which(ROV_lit$Predator=="Apolemia")] <- "Apolemia sp" 
ROV_lit$Predator[which(ROV_lit$Predator=="Bargmannia")] <- "Bargmannia elongata" #according to Steve June 02, 4.39PM
ROV_lit$Predator[which(ROV_lit$Predator=="Forskalia")] <- "Forskalia sp." 
ROV_lit$Predator[which(ROV_lit$Predator=="Nanomia bijuga")] <- "Nanomia sp. deep"
ROV_pruned <- ROV_lit[which(ROV_lit$Predator %in% unique(meltedGC$Species)),]

ROV_pruned$Predator <- factor(ROV_pruned$Predator, levels=SPorder)

ggplot(ROV_pruned, aes(x = Prey_broad, y = Predator, fill = log(Num_interactions))) + geom_tile() + scale_fill_gradient(low = "white", high = "black") + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), axis.title.y=element_blank()) 

GC_lit <- lit[which(!(lit$source %in% c("Pugh2009cResomiidae", "Choy_2017", "Hissman2005_newrho"))),]
GC_lit$Predator[which(GC_lit$Predator=="Nanomia bijuga" & GC_lit$source=="Biggs1977b_Feeding")] <- "Nanomia sp. Atlantic"
GC_lit$Prey_broad[which(GC_lit$Predator=="Nanomia sp. Atlantic" & GC_lit$source=="Biggs1977b_Feeding")] <- "Mysid"
GC_lit$Predator[which(GC_lit$Predator=="Nanomia bijuga" & GC_lit$source=="Purcell1981b_feeding")] <- "Nanomia sp. shallow"
GC_lit$Predator[which(GC_lit$Predator=="Apolemia uvaria")] <- "Apolemia sp" #questionable choice...
GC_lit$Predator[which(GC_lit$Predator=="Forskalia edwardsi")] <- "Forskalia sp." 
GC_lit$Predator[which(GC_lit$Predator=="Forskalia tholoides")] <- "Forskalia sp." 
GC_lit <- GC_lit[which(GC_lit$Prey != "Crustacea"),]
GC_pruned <- GC_lit[which(GC_lit$Predator %in% unique(meltedGC$Species)),]

ggplot(GC_pruned, aes(x = Prey_broad, y = Predator, fill = Num_interactions)) + geom_tile(color="grey") + scale_fill_gradient(low = "black", high = "black") + theme_bw() + theme(panel.grid.major = element_blank()) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), axis.title.y=element_blank())

MorphDAPC <- read.csv("Data/Literature-data/DAPCdiet_posteriors.tsv", sep="\t", stringsAsFactors = F)
#MorphDAPC <- MorphDAPC[which(MorphDAPC$Species %in% unique(meltedGC$Species)),]
meltedMorph <- melt(MorphDAPC)
names(meltedMorph)[2:3] <- c("Prey","Predicted")
meltedMorph$Species[which(meltedMorph$Species=="Physonect sp")] <- "Undescribed physonect L"
meltedMorph$Species[which(meltedMorph$Species=="Forskalia formosa")] <- "Forskalia sp."
meltedMorph$Species[which(meltedMorph$Species=="Forskalia asymmetrica")] <- "Forskalia sp."
meltedMorph$Species[which(meltedMorph$Species=="Forskalia edwardsii")] <- "Forskalia sp."

ggplot(meltedMorph, aes(x = Prey, y = Species, fill = Predicted)) + geom_tile() + scale_fill_gradient(low = "white", high = "black")+ theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), axis.title.y=element_blank())

morph <- rbind(mutate(meltedMorph[which(meltedMorph$Prey == "Large.crustacean"),], Prey_extrapolate="Decapod"),
      mutate(meltedMorph[which(meltedMorph$Prey == "Large.crustacean"),], Prey_extrapolate="Euphausid"),
      mutate(meltedMorph[which(meltedMorph$Prey == "Large.crustacean"),], Prey_extrapolate="Mysid"),
      mutate(meltedMorph[which(meltedMorph$Prey == "Large.crustacean"),], Prey_extrapolate="Lophogastrid"),
      mutate(meltedMorph[which(meltedMorph$Prey == "Large.crustacean"),], Prey_extrapolate="Stomatopod"),
      mutate(meltedMorph[which(meltedMorph$Prey == "Large.crustacean"),], Prey_extrapolate="Amphipod"),
      mutate(meltedMorph[which(meltedMorph$Prey == "Small.crustacean"),], Prey_extrapolate="Copepod"),
      mutate(meltedMorph[which(meltedMorph$Prey == "Small.crustacean"),], Prey_extrapolate="Cladoceran"),
      mutate(meltedMorph[which(meltedMorph$Prey == "Small.crustacean"),], Prey_extrapolate="Ostracod"),
      mutate(meltedMorph[which(meltedMorph$Prey == "Fish"),], Prey_extrapolate="Actinopteri"),
      mutate(meltedMorph[which(meltedMorph$Prey == "Mixed"),], Prey_extrapolate="Actinopteri"),
      mutate(meltedMorph[which(meltedMorph$Prey == "Mixed"),], Prey_extrapolate="Decapod"),
      mutate(meltedMorph[which(meltedMorph$Prey == "Mixed"),], Prey_extrapolate="Euphausid"),
      mutate(meltedMorph[which(meltedMorph$Prey == "Mixed"),], Prey_extrapolate="Mysid"),
      mutate(meltedMorph[which(meltedMorph$Prey == "Mixed"),], Prey_extrapolate="Lophogastrid"),
      mutate(meltedMorph[which(meltedMorph$Prey == "Mixed"),], Prey_extrapolate="Stomatopod"),
      mutate(meltedMorph[which(meltedMorph$Prey == "Mixed"),], Prey_extrapolate="Amphipod"),
      mutate(meltedMorph[which(meltedMorph$Prey == "Mixed"),], Prey_extrapolate="Copepod"),
      mutate(meltedMorph[which(meltedMorph$Prey == "Mixed"),], Prey_extrapolate="Cladoceran"),
      mutate(meltedMorph[which(meltedMorph$Prey == "Mixed"),], Prey_extrapolate="Ostracod"),
      mutate(meltedMorph[which(meltedMorph$Prey == "Mixed"),], Prey_extrapolate="Mollusc"),
      mutate(meltedMorph[which(meltedMorph$Prey == "Mixed"),], Prey_extrapolate="Polychaete"),
      mutate(meltedMorph[which(meltedMorph$Prey == "Mixed"),], Prey_extrapolate="Chaetognath"),
      mutate(meltedMorph[which(meltedMorph$Prey == "Mixed"),], Prey_extrapolate="Larvacean"),
      mutate(meltedMorph[which(meltedMorph$Prey == "Gelatinous"),], Prey_extrapolate="Gelatinous")
)
morph <- cbind(morph[,c(1,4,3)],rep("Morphology",nrow(morph)))
names(morph) <- c("Species","Prey","N","Source")

metabar <- cbind(GCspecies_melt, rep("Metabarcoding",nrow(GCspecies_melt)))
names(metabar) <- c("Species","Prey","N","Source")
lit_deep <- cbind(ROV_pruned[,c(1,4,6)], rep("Literature_Deep",nrow(ROV_pruned)))
names(lit_deep) <- c("Species","Prey","N","Source")
lit_shallow <- cbind(GC_pruned[,c(1,4,6)], rep("Literature_Shallow",nrow(GC_pruned)))
names(lit_shallow) <- c("Species","Prey","N","Source")
sources_combo <- rbind(metabar, lit_deep, lit_shallow)
sources_combo$Prey[which(sources_combo$Prey == "Fish")]<-"Actinopteri"
sources_combo$Prey[which(sources_combo$Prey == "Euphausiid")]<-"Euphausid"
#sources_combo$N[which(sources_combo$N == 0)]<-NA
sources_combo=sources_combo[which(sources_combo$Prey != "Crustacea"),]
sources_combo$Prey[which(sources_combo$Prey == "Salp")]<-"Gelatinous"
sources_combo$Prey[which(sources_combo$Prey == "Hydrozoa")]<-"Gelatinous"
sources_combo$Prey[which(sources_combo$Prey == "Medusae")]<-"Gelatinous"
sources_combo$Prey[which(sources_combo$Prey == "Scyphomedusa")]<-"Gelatinous"
sources_combo$Prey[which(sources_combo$Prey == "Scyphozoa")]<-"Gelatinous"
sources_combo$Prey[which(sources_combo$Prey == "Ctenophore")]<-"Gelatinous"
sources_combo$Prey[which(sources_combo$Prey == "Bivalve")]<-"Mollusc"
sources_combo$Prey[which(sources_combo$Prey == "Gastropod")]<-"Mollusc"
sources_combo$Prey[which(sources_combo$Prey == "Cephalopod")]<-"Mollusc"
sources_combo$N[which(sources_combo$N != 0)]<-1

PreyOrderComb <- c("Copepod","Decapod","Euphausid","Mysid","Lophogastrid","Stomatopod","Amphipod","Cladoceran","Ostracod","Mollusc","Polychaete","Chaetognath","Gelatinous","Larvacean","Actinopteri")
sources_combo$Prey <- factor(sources_combo$Prey,levels=PreyOrderComb)

sources_combo_nulls_M <- sources_combo
sources_combo_nulls_M$N <- 0
sources_combo_nulls_M$Source <- "Metabarcoding"
sources_combo_nulls_D <- sources_combo_nulls_M
sources_combo_nulls_D$Source <- "Literature_Deep"
sources_combo_nulls_S <- sources_combo_nulls_D
sources_combo_nulls_S$Source <- "Literature_Shallow"
sources_combo <- rbind(sources_combo,sources_combo_nulls_M,sources_combo_nulls_D,sources_combo_nulls_S)
#sources_combo <- unique(sources_combo)

morph <- morph[which(morph$Species %in% sources_combo$Species),]


# nudgesX <- function(x){
#   if(x=="Metabarcoding"){return(0.1)}
#   else if(x=="Litetarure_Deep"){return(0)}
#   else return(-0.1)
# }
# 
# nudgesY <- function(y){
#   if(y=="Metabarcoding"){return(-0.1)}
#   else if(x=="Litetarure_Deep"){return(0.1)}
#   else return(-0.1)
# }
# 
# sources_combo <- mutate(sources_combo, nudgeX = nudgesX(Source), nudgeY = nudgesY(Source))

M <- ggplot(sources_combo[which(sources_combo$Source=="Metabarcoding"),], aes(x=Prey,y=Species)) + geom_point(shape=15,position = position_nudge(0.08,-0.16), col=ifelse(sources_combo[which(sources_combo$Source=="Metabarcoding"),]$N == 0, NA, "red"),cex=2.5) + theme_bw()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), axis.title.y=element_blank())
MD <- M + geom_point(data = sources_combo[which(sources_combo$Source=="Literature_Deep"),], shape=15, position = position_nudge(-0.08,0.16), col=ifelse(sources_combo[which(sources_combo$Source=="Literature_Deep"),]$N == 0, NA, "blue"),cex=2.5) + theme_bw()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), axis.title.y=element_blank())
MDS <- MD + geom_point(data = sources_combo[which(sources_combo$Source=="Literature_Shallow"),],shape=15, position = position_nudge(-0.08,-0.16), col=ifelse(sources_combo[which(sources_combo$Source=="Literature_Shallow"),]$N == 0, NA, "green"),cex=2.5) + theme_bw()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), axis.title.y=element_blank())
MDS + geom_point(data = morph, position = position_nudge(0.08,0.16),shape=15, col=ifelse(morph$N<0.01, NA, "purple"),cex=2.5, aes(alpha=log(N+1))) + theme_bw()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), axis.title.y=element_blank())

#Supplementary compoplots ##### AND NOW SUPPLEMENTARY TABLES 

ref_samples <- allsamples[which(!(allsamples$Interpretation %in% c("Short","Unknown"))),]
ref_samples$Interpretation[which(ref_samples$Interpretation == "Cross")] <- "Contamination"

ref_samples$Broad.group[which(ref_samples$Broad.group %in% c("Hydroid","Hydrozoan","Hydromedusa"))] <- "Hydrozoan"

BG_order <- c("Siphonophore", PreyOrderTr[which(PreyOrderTr %in% unique(ref_samples$Broad.group) & PreyOrderTr != "Siphonophore")],"Sponge", "Ascidian", "Hydrozoan","Anemone", "Branchiopod", "Bryozoan",  "Nemertean", "Sipunculan", "Echiurid", "Rotifer", "Gastrotrich", "Platyctenid", "Orthonectid", "Trematode", "Cestode", "Nematode", "Myxozoan", "Ichthyiophonid", "Mite", "Insect", "Oligochaete", "Shark", "Tetrapod", "Pollen", "Eukaryote", "Bacteria")
#group_samples <- ref_samples[which(!(ref_samples$Interpretation %in% c("Environmental","Contamination"))),]

#species by interpretation

ref_samples[,c("Species","Interpretation","abundance")] %>% dcast(Species~Interpretation, value.var="abundance", fun.aggregate = sum) -> Isp_table
Isp_table[,c("Species", "Predator", "Prey", "Secondary", "Parasite", "Environmental", "Contamination")] -> Isp_table
write.table(Isp_table, "PLOSOne_Revisions/Supplement_tables/SM_Table1_Isp.tsv", sep="\t", row.names = F)

pdf("SM_spp_interpretation.pdf", width = 20, height = 10)
ggplot(ref_samples, aes(x = Species, y = log(abundance), fill = Interpretation)) + geom_bar(position = "fill",stat = "identity") + theme_bw() + theme(legend.key.size = unit(1,"line"), legend.text = element_text(size=5),legend.position="bottom", axis.text.x = element_text(angle = 90, hjust = 1)) + scale_fill_manual(values=c(colorRampPalette(brewer.pal(8, "Set2"))(8))) + guides(shape = guide_legend(override.aes = list(size = 1))) 
dev.off()

#species by interpretation and barcode

barIsp_list <- lapply(as.list(barcodes), function(B){ref_samples[which(ref_samples$barcode == B),c("Species","Interpretation","abundance")] %>% as.data.frame() %>% dcast(Species~Interpretation, value.var="abundance", fun.aggregate = sum) -> Ti; Ti[,c("Species", "Predator", "Prey", "Secondary", "Parasite", "Environmental", "Contamination")] -> Ti_ord;return(mutate(Ti_ord, Barcode = B))})
barIsp_table <- bind_rows(barIsp_list)
barIsp_table[is.na(barIsp_table)] <- 0
barIsp_table <- barIsp_table[,c("Species", "Barcode", "Predator", "Prey", "Secondary", "Parasite", "Environmental", "Contamination")]
write.table(barIsp_table, "PLOSOne_Revisions/Supplement_tables/SMTable2_barIsp.tsv", sep="\t", row.names = F)


pdf("SM_spp_interpretation_barcode.pdf", width = 20, height = 10)
ggplot(ref_samples, aes(x = Species, y = log(abundance), fill = Interpretation)) + geom_bar(position = "fill",stat = "identity") + theme_bw() + theme(legend.key.size = unit(1,"line"), legend.text = element_text(size=5),legend.position="bottom", axis.text.x = element_text(angle = 90, hjust = 1)) + scale_fill_manual(values=c(colorRampPalette(brewer.pal(8, "Set2"))(8))) + guides(shape = guide_legend(override.aes = list(size = 1))) + facet_grid(.~barcode)
dev.off()

#specimens by interpretation

ref_samples[,c("Specimen","Species","Interpretation","abundance")] %>% dcast(Specimen+Species~Interpretation, value.var="abundance", fun.aggregate = sum) -> Iind_table
Iind_table[order(Iind_table$Species),c("Specimen","Species", "Predator", "Prey", "Secondary", "Parasite", "Environmental", "Contamination")] -> Iind_table
write.table(Iind_table, "PLOSOne_Revisions/Supplement_tables/SMTable3_Iind.tsv", sep="\t", row.names = F)

pdf("SM_ind_interpretation.pdf", width = 20, height = 10)
ggplot(ref_samples, aes(x = paste(Species, Specimen), y = log(abundance), fill = Interpretation)) + geom_bar(position = "fill",stat = "identity") + theme_bw() + theme(legend.key.size = unit(1,"line"), legend.text = element_text(size=5),legend.position="bottom", axis.text.x = element_text(angle = 90, hjust = 1)) + scale_fill_manual(values=c(colorRampPalette(brewer.pal(8, "Set2"))(8))) + guides(shape = guide_legend(override.aes = list(size = 1))) + guides(shape = guide_legend(override.aes = list(size = 1)))
dev.off()

#specimens by interpretation and barcode

barIind_list <- lapply(as.list(barcodes), function(B){ref_samples[which(ref_samples$barcode == B),c("Specimen","Species","Interpretation","abundance")] %>% dcast(Specimen+Species~Interpretation, value.var="abundance", fun.aggregate = sum) -> Ti; Ti[order(Ti$Species),c("Specimen","Species", "Predator", "Prey", "Secondary", "Parasite", "Environmental", "Contamination")] -> Ti_ord;return(mutate(Ti_ord, Barcode = B))})
barIind_table <- bind_rows(barIind_list)
barIind_table[is.na(barIind_table)] <- 0
barIind_table <- barIind_table[order(barIind_table$Barcode),c("Specimen","Species", "Barcode", "Predator", "Prey", "Secondary", "Parasite", "Environmental", "Contamination")]
write.table(barIind_table, "PLOSOne_Revisions/Supplement_tables/SMTable4_barIind.tsv", sep="\t", row.names = F)

pdf("SM_ind_interpretation_barcode.pdf", width = 20, height = 30)
ggplot(ref_samples, aes(x = paste(Species, Specimen), y = log(abundance), fill = Interpretation)) + geom_bar(position = "fill",stat = "identity") + theme_bw() + theme(legend.key.size = unit(1,"line"), legend.text = element_text(size=5),legend.position="bottom", axis.text.x = element_text(angle = 90, hjust = 1)) + scale_fill_manual(values=c(colorRampPalette(brewer.pal(8, "Set2"))(8))) + guides(shape = guide_legend(override.aes = list(size = 1))) + guides(shape = guide_legend(override.aes = list(size = 1))) + facet_wrap(.~barcode, nrow=6, ncol=1)
dev.off()

#species by broad group

ref_samples[,c("Species","Broad.group", "abundance")] %>% dcast(Species~Broad.group, value.var="abundance", fun.aggregate = sum) -> BGsp_table
BGsp_table[,c(1,match(BG_order, colnames(BGsp_table)))] -> BGsp_table
write.table(BGsp_table, "PLOSOne_Revisions/Supplement_tables/SMTable5_BGsp.tsv", sep="\t", row.names = F)

pdf("SM_spp_broadgroup.pdf", width = 20, height = 10)
ggplot(ref_samples, aes(x = Species, y = log(abundance), fill = Broad.group)) + geom_bar(position = "fill",stat = "identity") + theme_bw() + theme(legend.key.size = unit(1,"line"), legend.text = element_text(size=5),legend.position="bottom", axis.text.x = element_text(angle = 90, hjust = 1)) + scale_fill_manual(values=c(colorRampPalette(brewer.pal(8, "Set2"))(16), colorRampPalette(brewer.pal(8, "Set1"))(16), colorRampPalette(brewer.pal(8, "Dark2"))(23))) + guides(shape = guide_legend(override.aes = list(size = 1))) 
dev.off()

#species by broad group and barcode
barBGsp_list <- lapply(as.list(barcodes), function(B){ref_samples[which(ref_samples$barcode == B),c("Species","Broad.group","abundance")] %>% dcast(Species~Broad.group, value.var="abundance",fun.aggregate = sum) -> Ti; Ti[,c(1,match(BG_order[which(BG_order %in% colnames(Ti))], colnames(Ti)))] -> Ti_ord;return(mutate(Ti_ord, Barcode = B))})
barBGsp_table <- bind_rows(barBGsp_list)
barBGsp_table[is.na(barBGsp_table)] <- 0
barBGsp_table <- barBGsp_table[,c(1,38,match(BG_order, colnames(barBGsp_table)))]
write.table(barBGsp_table, "PLOSOne_Revisions/Supplement_tables/SMTable6_barBGsp.tsv", sep="\t", row.names = F)

pdf("SM_spp_broadgroup_barcode.pdf", width = 20, height = 10)
ggplot(ref_samples, aes(x = Species, y = log(abundance), fill = Broad.group)) + geom_bar(position = "fill",stat = "identity") + theme_bw() + theme(legend.key.size = unit(1,"line"), legend.text = element_text(size=5),legend.position="bottom", axis.text.x = element_text(angle = 90, hjust = 1)) + scale_fill_manual(values=c(colorRampPalette(brewer.pal(8, "Set2"))(16), colorRampPalette(brewer.pal(8, "Set1"))(16), colorRampPalette(brewer.pal(8, "Dark2"))(23))) + guides(shape = guide_legend(override.aes = list(size = 1))) +  facet_grid(.~barcode)
dev.off()


#specimens by broad group

ref_samples[,c("Specimen","Species","Broad.group","abundance")] %>% dcast(Specimen+Species~Broad.group, value.var="abundance", fun.aggregate = sum) -> BGind_table
BGind_table[order(BGind_table$Species),c(1,2,match(BG_order, colnames(BGind_table)))] -> BGind_table
write.table(BGind_table, "PLOSOne_Revisions/Supplement_tables/SMTable7_BGind.tsv", sep="\t", row.names = F)

pdf("SM_ind_broadgroup.pdf", width = 20, height = 10)
ggplot(ref_samples, aes(x = paste(Species, Specimen), y = log(abundance), fill = Broad.group)) + geom_bar(position = "fill",stat = "identity") + theme_bw() + theme(legend.key.size = unit(1,"line"), legend.text = element_text(size=5),legend.position="bottom", axis.text.x = element_text(angle = 90, hjust = 1)) + scale_fill_manual(values=c(colorRampPalette(brewer.pal(8, "Set2"))(16), colorRampPalette(brewer.pal(8, "Set1"))(16), colorRampPalette(brewer.pal(8, "Dark2"))(23))) + guides(shape = guide_legend(override.aes = list(size = 1)))
dev.off()

#specimens by broad group and barcode
barBGind_list <- lapply(as.list(barcodes), function(B){ref_samples[which(ref_samples$barcode == B),c("Specimen","Species","Broad.group","abundance")] %>% dcast(Specimen+Species~Broad.group, value.var="abundance", fun.aggregate = sum) -> Ti; Ti[order(Ti$Species),c(1,2,match(BG_order[which(BG_order %in% colnames(Ti))], colnames(Ti)))] -> Ti_ord;return(mutate(Ti_ord, Barcode = B))})
barBGind_table <- bind_rows(barBGind_list)
barBGind_table[is.na(barBGind_table)] <- 0
barBGind_table <- barBGind_table[order(barBGind_table$Barcode),c(1,2,39,match(BG_order, colnames(barBGind_table)))]
write.table(barBGind_table, "PLOSOne_Revisions/Supplement_tables/SMTable8_barBGind.tsv", sep="\t", row.names = F)

pdf("SM_ind_broadgroup_barcode.pdf", width = 20, height = 30)
ggplot(ref_samples, aes(x = paste(Species, Specimen), y = log(abundance), fill = Broad.group)) + geom_bar(position = "fill",stat = "identity") + theme_bw() + theme(legend.key.size = unit(1,"line"), legend.text = element_text(size=5),legend.position="bottom", axis.text.x = element_text(angle = 90, hjust = 1)) + scale_fill_manual(values=c(colorRampPalette(brewer.pal(8, "Set2"))(16), colorRampPalette(brewer.pal(8, "Set1"))(16), colorRampPalette(brewer.pal(8, "Dark2"))(23))) + guides(shape = guide_legend(override.aes = list(size = 1))) + facet_wrap(.~barcode, nrow=6, ncol=1)
dev.off()

#############
# ONLY PREY -- Without siph signal #
prey_samples <- ref_samples[which(ref_samples$Interpretation == "Prey"),]

#species by broad group

prey_samples[,c("Species","Broad.group","abundance")] %>% dcast(Species~Broad.group, value.var="abundance", fun.aggregate = sum) -> BGspPrey_table
BGspPrey_table[,c(1,match(PreyOrder[which(PreyOrder %in% colnames(BGspPrey_table))], colnames(BGspPrey_table)))] -> BGspPrey_table
write.table(BGspPrey_table, "PLOSOne_Revisions/Supplement_tables/SMTable9_BGspPrey.tsv", sep="\t", row.names = F)

pdf("SM_spp_broadgroup_prey.pdf", width = 20, height = 10)
ggplot(prey_samples, aes(x = Species, y = log(abundance), fill = Broad.group)) + geom_bar(position = "fill",stat = "identity") + theme_bw() + theme(legend.key.size = unit(1,"line"), legend.text = element_text(size=5),legend.position="bottom", axis.text.x = element_text(angle = 90, hjust = 1)) + scale_fill_manual(values=c(colorRampPalette(brewer.pal(6, "Set1"))(14))) + guides(shape = guide_legend(override.aes = list(size = 1))) 
dev.off()

#specimens by broad group
prey_samples[,c("Specimen","Species","Broad.group","abundance")] %>% dcast(Specimen+Species~Broad.group, value.var="abundance", fun.aggregate = sum) -> BGindPrey_table
BGindPrey_table[order(BGindPrey_table$Species),c(1,2,match(PreyOrder[which(PreyOrder %in% colnames(BGindPrey_table))], colnames(BGindPrey_table)))] -> BGindPrey_table
write.table(BGindPrey_table, "PLOSOne_Revisions/Supplement_tables/SMTable10_BGindPrey.tsv", sep="\t", row.names = F)

pdf("SM_ind_broadgroup_prey.pdf", width = 20, height = 10)
ggplot(prey_samples, aes(x = paste(Species, Specimen), y = log(abundance), fill = Broad.group)) + geom_bar(position = "fill",stat = "identity") + theme_bw() + theme(legend.key.size = unit(1,"line"), legend.text = element_text(size=5),legend.position="bottom", axis.text.x = element_text(angle = 90, hjust = 1)) + scale_fill_manual(values=c(colorRampPalette(brewer.pal(6, "Set1"))(14))) + guides(shape = guide_legend(override.aes = list(size = 1))) 
dev.off()

#species by broad group and barcode

barBGspPrey_list <- lapply(as.list(barcodes), function(B){prey_samples[which(prey_samples$barcode == B),c("Species","Broad.group","abundance")] %>% dcast(Species~Broad.group, value.var="abundance", fun.aggregate = sum) -> Ti; Ti[,c(1,match(PreyOrder[which(PreyOrder %in% colnames(Ti))], colnames(Ti)))] -> Ti_ord;return(mutate(Ti_ord, Barcode = B))})
barBGspPrey_table <- bind_rows(barBGspPrey_list)
barBGspPrey_table[is.na(barBGspPrey_table)] <- 0
barBGspPrey_table <- barBGspPrey_table[,c(1,10,match(PreyOrder[which(PreyOrder %in% colnames(BGspPrey_table))], colnames(barBGspPrey_table)))]
write.table(barBGsp_table, "PLOSOne_Revisions/Supplement_tables/SMTable11_barBGspPrey.tsv", sep="\t", row.names = F)

pdf("SM_spp_broadgroup_barcode_prey.pdf", width = 20, height = 6)
ggplot(prey_samples, aes(x = Species, y = log(abundance), fill = Broad.group)) + geom_bar(position = "fill",stat = "identity") + theme_bw() + theme(legend.key.size = unit(1,"line"), legend.text = element_text(size=5),legend.position="bottom", axis.text.x = element_text(angle = 90, hjust = 1)) + scale_fill_manual(values=c(colorRampPalette(brewer.pal(6, "Set1"))(14))) + guides(shape = guide_legend(override.aes = list(size = 1))) +  facet_grid(.~barcode)
dev.off()

#specimens by broad group and barcode
barBGindPrey_list <- lapply(as.list(barcodes), function(B){prey_samples[which(prey_samples$barcode == B),c("Specimen","Species","Broad.group","abundance")] %>% dcast(Specimen+Species~Broad.group, value.var="abundance", fun.aggregate = sum) -> Ti; Ti[order(Ti$Species),c(1,2,match(PreyOrder[which(PreyOrder %in% colnames(Ti))], colnames(Ti)))] -> Ti_ord;return(mutate(Ti_ord, Barcode = B))})
barBGindPrey_table <- bind_rows(barBGindPrey_list)
barBGindPrey_table[is.na(barBGindPrey_table)] <- 0
barBGindPrey_table <- barBGindPrey_table[,c(1,2,11,match(PreyOrder[which(PreyOrder %in% colnames(barBGspPrey_table))], colnames(barBGspPrey_table)))] %>% .[,-12]
write.table(barBGindPrey_table, "PLOSOne_Revisions/Supplement_tables/SMTable12_barBGindPrey.tsv", sep="\t", row.names = F)

pdf("SM_ind_broadgroup_barcode_prey.pdf", width = 20, height = 10)
ggplot(prey_samples, aes(x = paste(Species, Specimen), y = log(abundance), fill = Broad.group)) + geom_bar(position = "fill",stat = "identity") + theme_bw() + theme(legend.key.size = unit(1,"line"), legend.text = element_text(size=5),legend.position="bottom", axis.text.x = element_text(angle = 90, hjust = 1)) + scale_fill_manual(values=c(colorRampPalette(brewer.pal(6, "Set1"))(14))) + guides(shape = guide_legend(override.aes = list(size = 1)))  + facet_wrap(.~barcode, nrow=6, ncol=1)
dev.off()

######
#Unique sequences 
#######
sum_samples <- allsamples
sum_samples$Interpretation[which(sum_samples$Interpretation == "Cross")] <- "Contamination"
sum_samples$Interpretation[which(sum_samples$Interpretation == "Short")] <- "Unknown"
sum_samples %>% distinct(Sequence, .keep_all = TRUE) -> uniseq_samples

# interpretation by specimen and barcode
type_order <- c("Specimen","Species", "Predator", "Prey", "Secondary", "Parasite", "Environmental", "Contamination", "Unknown")

FEATbarIind_list <- lapply(as.list(barcodes), function(B){uniseq_samples[which(uniseq_samples$barcode == B),c("Specimen","Species","Interpretation","Sequence")] %>% dcast(Specimen+Species~Interpretation, fun.aggregate = function(X){length(unique(X))}, value.var = "Sequence") -> Ti;print(Ti %>% head()); Ti[order(Ti$Species),type_order[which(type_order %in% colnames(Ti))]] -> Ti_ord;return(mutate(Ti_ord, Barcode = B))})
FEATbarIind_table <- bind_rows(FEATbarIind_list)
FEATbarIind_table[is.na(FEATbarIind_table)] <- 0
FEATbarIind_table <- FEATbarIind_table[order(FEATbarIind_table$Barcode),c("Specimen","Species", "Barcode", "Predator", "Prey", "Secondary", "Parasite", "Environmental", "Contamination", "Unknown")]

group_by(FEATbarIind_table[,c(-1,-2)],Barcode) %>% 
  dplyr::summarise(Predator=sum(Predator), Prey = sum(Prey), Secondary = sum(Secondary), Parasite = sum(Parasite), Environmental = sum(Environmental), Contamination = sum(Contamination), Unknown = sum(Unknown)) %>% 
  mutate(Total = Predator+Prey+Secondary+Parasite+Environmental+Contamination+Unknown) %>% 
  add_row(Barcode="TOTAL", Predator=sum(.$Predator), Prey=sum(.$Prey), Secondary=sum(.$Secondary), Parasite=sum(.$Parasite), Environmental=sum(.$Environmental), Contamination=sum(.$Contamination), Unknown=sum(.$Unknown), Total=sum(.$Total)) %>% 
  write.table("PLOSOne_Revisions/Supplement_tables/SMTable13_FEATbarIind.tsv", sep="\t", row.names = F)

#broad group by specimen & barcode
BG_order <- c(BG_order,"Unknown")
FEATbarBGind_list <- lapply(as.list(barcodes), function(B){
  uniseq_samples[which(uniseq_samples$barcode == B),c("Specimen","Species","Broad.group","Sequence")] %>% 
    dcast(Specimen+Species~Broad.group, fun.aggregate = function(X){length(unique(X))}, value.var = "Sequence") -> Ti;
  Ti[order(Ti$Species),c(1,2,match(BG_order[which(BG_order %in% colnames(Ti))], colnames(Ti)))] -> Ti_ord;return(mutate(Ti_ord, Barcode = B))})
FEATbarBGind_table <- bind_rows(FEATbarBGind_list)
FEATbarBGind_table[is.na(FEATbarBGind_table)] <- 0
FEATbarBGind_table <- FEATbarBGind_table[order(FEATbarBGind_table$Barcode),c("Specimen","Species","Barcode", BG_order[which(BG_order %in% colnames(FEATbarBGind_table))])]

group_by(FEATbarBGind_table[,c(-1,-2)],Barcode) %>% 
  dplyr::summarise_all(sum) %>% 
  mutate(Total = Siphonophore+Copepod+Decapod+Euphausid+Mysid+Lophogastrid+Stomatopod+Ostracod+Bivalve+Larvacean++Salp+Scyphomedusa+Ctenophore+Actinopteri+Amphipod+Isopod+Echinoderm+Polychaete+Chaetognath+Gastropod+Doliolid+Pyrosome+Sponge+Ascidian+Hydrozoan+Anemone+Branchiopod+Bryozoan+Nemertean+Sipunculan+Echiurid+Rotifer+Platyctenid+Orthonectid+Trematode+Cestode+Nematode+Myxozoan+Ichthyiophonid+Mite+Insect+Oligochaete+Tetrapod+Pollen+Eukaryote+Bacteria+Unknown) %>% 
  add_row(Barcode="TOTAL", Siphonophore=sum(.$Siphonophore), Copepod=sum(.$Copepod), Decapod=sum(.$Decapod), Euphausid=sum(.$Euphausid), Mysid=sum(.$Mysid), Lophogastrid=sum(.$Lophogastrid), Stomatopod=sum(.$Stomatopod), Ostracod=sum(.$Ostracod), Bivalve=sum(.$Bivalve), Larvacean=sum(.$Larvacean), Salp=sum(.$Salp), Scyphomedusa=sum(.$Scyphomedusa), Ctenophore=sum(.$Ctenophore), Actinopteri=sum(.$Actinopteri), Amphipod=sum(.$Amphipod), Isopod=sum(.$Isopod), Echinoderm=sum(.$Echinoderm), Polychaete=sum(.$Polychaete), Chaetognath=sum(.$Chaetognath), Gastropod=sum(.$Gastropod), Doliolid=sum(.$Doliolid), Pyrosome=sum(.$Pyrosome), Sponge=sum(.$Sponge), Ascidian=sum(.$Ascidian), Hydrozoan=sum(.$Hydrozoan), Anemone=sum(.$Anemone), Branchiopod=sum(.$Branchiopod), Bryozoan=sum(.$Bryozoan), Nemertean=sum(.$Nemertean), Sipunculan=sum(.$Sipunculan), Echiurid=sum(.$Echiurid), Rotifer=sum(.$Rotifer), Platyctenid=sum(.$Platyctenid), Orthonectid=sum(.$Orthonectid), Trematode=sum(.$Trematode), Cestode=sum(.$Cestode), Nematode=sum(.$Nematode), Myxozoan=sum(.$Myxozoan), Ichthyiophonid=sum(.$Ichthyiophonid), Mite=sum(.$Mite), Insect=sum(.$Insect), Oligochaete=sum(.$Oligochaete), Tetrapod=sum(.$Tetrapod), Pollen=sum(.$Pollen), Eukaryote=sum(.$Eukaryote), Bacteria=sum(.$Bacteria), Unknown=sum(.$Unknown), Total=sum(.$Total)) %>% 
  write.table("PLOSOne_Revisions/Supplement_tables/SMTable14_FEATbarBGind.tsv", sep="\t", row.names = F)

############

#Summary.table
split(allsamples[,2:4],f=allsamples$Species) %>% lapply(function(x){unique(x$Specimen) %>% length()}) %>% unlist() -> Nsampled
Nsampled <- data.frame(Species = names(Nsampled), N_sampled = Nsampled)
dcast(meltedGC, formula=Species~Specimen, fun = sum, value.var="Presence") -> sps
data.frame(Species=sps$Species,N=rowSums(sps[,-1]))->sps
NS <- left_join(Nsampled,sps,by="Species")
spMorphGuilds <- do.call(rbind,split(meltedMorph,f=meltedMorph$Species) %>% lapply(function(x){return(x[which(x$Predicted==max(x$Predicted)),])}))[,1:2]
names(spMorphGuilds)[2] <- "MorphGuild"
NS <- left_join(NS,spMorphGuilds, by="Species")
PNASguilds <- c("Small crustacean", "Small crustacean", "Small crustacean", "Small crustacean", "Small crustacean", "Large crustacean", "Mixed", "Mixed", "Mixed", "Large crustacean", "Large crustacean", "Mixed", "Small crustacean", "Large crustacean", "Fish", "Fish", "Fish", "Large crustacean", "Gelatinous", "Fish", "Fish", "Fish", "Large crustacean", "Large crustacean")
names(PNASguilds) <- c("Sulculeolaria quadrivalvis","Chelophyes appendiculata","Diphyes dispar","Sphaeronectes koellikeri","Hippopodius hippopus","Praya dubia","Agalma okenii","Athorybia rosacea","Agalma elegans","Nanomia sp. deep","Lychnagalma utricularia","Forskalia sp.","Cordagalma ordinatum","Resomia ornicephala","Erenna richardi","Erenna sirena","Stephanomia amphitrydis","Bargmannia amoena","Apolemia rubriversa","Rhizophysa eysenhardtii","Rhizophysa filiformis","Physalia physalis","Nanomia sp. Atlantic", "Nanomia sp. shallow")
NS <- left_join(NS,data.frame(Species=names(PNASguilds), PNAS_guilds = PNASguilds), by="Species")

NS$N[is.na(NS$N)] <- 0
write.csv(NS,"NStable.tsv",col.names = T,row.names = F,sep='\t')
