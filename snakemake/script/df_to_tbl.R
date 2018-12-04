args=(commandArgs(TRUE))
file<-args[1]
outfile<-args[2]
a<-read.table(file)
datb<-xtabs(V3 ~ V1 + V2, a)
datb<-t(datb)
write.table(datb,file=outfile,sep="\t")
