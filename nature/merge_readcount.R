sample<-dir("./","*_counts.txt")
sample1<-paste0(sort(sub("_counts.txt","",sample)),"_counts.txt")
samplename<-sort(sub("_count.txt","",sample))
count <- do.call("cbind", lapply(sample1, read.table, sep='\t', row.names=1))
colnames(count)<-samplename
write.table(count,"htseq_reads_count.txt",row.names = T,col.names = T,sep="\t",quote = FALSE)
