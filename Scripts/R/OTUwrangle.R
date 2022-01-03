#Import arguments
myargs = commandArgs(trailingOnly=TRUE)

#Declare variables for each argument
runid = myargs[1]
barcode = myargs[2]
Assignments <- read.csv(file = myargs[3], header = FALSE, sep = ",")[,c(4:10)]
Features <- read.csv(file = myargs[4], header = TRUE, sep = "\t")
outdir = myargs[5]

#Order the OTUs as in the features table
Features_pruned = Features[c(1,which(Features$ID %in% Assignments$V4)),]
Taxonomy <- Assignments[match(Features_pruned[2:nrow(Features_pruned),1], Assignments[,1]),2:7]

#Paste the ranks into a single string
Ranks_pasted = apply(Taxonomy, 1, function(x){paste(x, collapse="_")})

#Combine the taxonomic ranks ID with the feature statistics
IDtable <- cbind(c("OTU",Ranks_pasted), Features_pruned[,2:ncol(Features_pruned)])
IDtable<-apply(IDtable, 2, as.character)[-1,]

namevector <- as.character(colnames(IDtable))
namevector[1]<-"OTU"
tabtab<-rbind(namevector,IDtable)

#Write an outfile
write.table(x=tabtab, file=paste(outdir, runid, barcode, "IDtable.tsv", sep=""), col.names=FALSE,  row.names=FALSE)
