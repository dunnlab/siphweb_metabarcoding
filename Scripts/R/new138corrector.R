#MERGER

annotated_138 <- read.csv("~/Desktop/SILVA138_extractMar16.tsv", header = T, sep='\t', stringsAsFactors = F)
new_138 <- read.csv("~/Desktop/SILVA138plus_newparsed.tsv", header = T, sep='\t', stringsAsFactors = F)

sorted_ann138 <- annotated_138[order(annotated_138$run, annotated_138$barcode, annotated_138$sample, annotated_138$Feature_name),]
sorted_new138 <- new_138[order(new_138$run, new_138$barcode, new_138$sample, new_138$Feature_name),]


sorted_new138[,2:4]<-sorted_ann138[,2:4]
sorted_new138$Broad.group<-sorted_ann138$Broad.group
sorted_new138$Interpretation<-sorted_ann138$Interpretation
sorted_new138$Sequence <- sorted_ann138$Sequence

write.table(sorted_new138, "SILVA138plus_fixed.tsv", sep="\t", row.names = F, col.names = T)
