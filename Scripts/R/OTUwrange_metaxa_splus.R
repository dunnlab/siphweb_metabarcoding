library(stringr)

#Import arguments
myargs = commandArgs(trailingOnly=TRUE)

#Declare variables for each argument
runid = myargs[1]
barcode = myargs[2]
Assignments <- read.csv(file = myargs[3], header = FALSE, sep = "\t")
Features <- read.csv(file = myargs[4], header = FALSE, sep = "\t",stringsAsFactors=F)
outdir = myargs[5]

#Order the OTUs as in the features table
Features_pruned = Features[c(1,which(Features$V1 %in% Assignments$V1)),]
Taxonomy <- Assignments[match(Features_pruned[2:nrow(Features_pruned),1], Assignments[,1]),2:ncol(Assignments)]
print(class(Features_pruned[,1]))

#Paste the ranks into a single string
Ranks_pasted = apply(Taxonomy, 1, function(x){paste(x, collapse="_")})

#Combine the taxonomic ranks ID with the feature statistics
IDtable <- cbind(c("OTU",Ranks_pasted), c("Feature_name",Features_pruned[-1,1]), Features_pruned[,2:ncol(Features_pruned)])
sapply(IDtable,as.character)->IDtable
IDtable[,1]<-str_replace_all(IDtable[,1], "\\(\\d+\\.\\d+\\)", replacement="")
IDtable[,1]<-str_replace_all(IDtable[,1], "\\(\\d+\\)", replacement="")

#Write an outfile
IDtable[1,] %>% as.list() %>% as.character() -> colnames(IDtable)
IDtable<-IDtable[-1,]
IDtable[,3:ncol(IDtable)]<-sapply(IDtable[,3:ncol(IDtable)],as.numeric)
write.table(x=IDtable, file=paste(outdir, runid, barcode, "IDtable_custom.tsv", sep=""), col.names=TRUE,  row.names=FALSE)
