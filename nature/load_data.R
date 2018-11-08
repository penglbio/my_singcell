########QC##########
rm(list=ls())
library(scater); library(destiny) ;library("ROCR");library(GenomicRanges)
library(M3Drop);  library(M3DExampleData); library(SC3); .pardefault <- par()
library("beanplot")



argv <- commandArgs(T)
argv[2] <- "Mus_musculus.GRCm38.84.renameChr.bed"
argv[3] <- "htseq_reads_count.txt"
argv[4] <- "0.4-400000"


cellinfo_table <- argv[1]
gene_ano <- argv[2]
ht_count <- argv[3]
align_ro_R<-as.numeric(strsplit(argv[4],"-")[[1]])




setwd("/var/data/07.single_cell/mydata/")
conpath <- "/var/data/07.single_cell/Areg/data/"

load(paste0(conpath,"CellInfo_MergedASC.RData"))
cellinfo=bigtable.ok
cellinfo[which(cellinfo[,"Cell"]%in%"ND"),"Cell"]="NE"; cellinfo[which(cellinfo[,"Cell"]%in%c("2RD","RDE")),"Cell"]="RD"
cellinfo[which(cellinfo[,"Cell"]%in%"R+N"),"Cell"]="RN"; cellinfo[which(cellinfo[,"Cell"]%in%"2R+N"),"Cell"]="RN"
cellinfo[which(cellinfo[,"Cell"]%in%"3N+R"),"Cell"]="RN"; cellinfo[which(cellinfo[,"Cell"]%in%"POS"),"Cell"]="R"
cellinfo[which(cellinfo[,"Cell"]%in%"NEG"),"Cell"]="N"
exclude_doublet=unique(c(which(!cellinfo[,"Qual"]%in%c("","ok","super","ok (weird profile; possibly degraded)")),
                         which(!cellinfo[,"Cell"]%in%c("N","R","NE","R?","RD")),grep("N[.]2[.]",rownames(cellinfo)),grep("N[.]1d[.]",rownames(cellinfo)),
                         grep("P[.]d[.]",rownames(cellinfo)),grep("P[.]0[.]",rownames(cellinfo))))

doublet_cellinfo=rownames(cellinfo)[exclude_doublet] #### unique names substr(cellinfo.ok[,"Name"],1,8)
doublet_id<-sub("_[0|1]_$","",gsub("\\.","_",doublet_cellinfo[grep("^D",doublet_cellinfo)]))
doublet_id <- gsub("D_","BD_",sub("01|00","",doublet_id))
mynames=substr(mysplit(filemap[rownames(cellinfo.ok)],"09-",2,2),1,39)
names(mynames)=gsub("[.][.]","",paste(gsub("BD","D",gsub("_",".",substr(mynames,1,9))), substr(mynames,38,39),sep=""))


###load gene_anno table
gene84<- read.table(gene_ano,header = T,sep="\t",stringsAsFactors = F,row.names = 1)
colnames(gene84)[4:6] <- c("strand","geneName","geneBiotype")
head(gene84)
gene84.GR<- makeGRangesFromDataFrame(gene84,keep.extra.columns = T,ignore.strand = F,seqnames.field = "seqnames")
gene84.GR["ENSMUSG00000040856",]
###load ht_count table
counts <- read.table(ht_count,header = T,sep="\t",row.names = 1,stringsAsFactors = F)
colnames(counts) <- gsub("_counts.txt","",colnames(counts))
colSums(counts)
###load mapped
maplog <- read.table("maplog.txt",header = T,sep="\t",stringsAsFactors = F,row.names = 1)
mapped_ratio <- 1-rowSums(apply(maplog[,c(7,8,9)],2,function(x){as.numeric(sub("%","",x))}))/100
sample_0.4 <- rownames(maplog[mapped_ratio<=0.4,])
sample_400000 <- rownames(maplog[maplog[,1]*mapped_ratio <=400000,])
exclude <- union(union(sample_400000,sample_0.4),doublet_id)
##ok sample logfile and counts
head(counts)
countsMatrix <- counts[,!(colnames(counts) %in% exclude)]

maplog <- maplog[!(colnames(counts) %in% exclude),]
##QC figures data
genesPerCell <- apply(countsMatrix,2,function(x) length(which(x>0)))
genereadPerCell <- apply(countsMatrix,2,function(x) sum(x))/1000000
readsPerCell <- as.numeric(apply(maplog,1,function(x) x[1]))/1000000
names(readsPerCell) <- names(genesPercell)
countMatrix.cpm=sapply(1:length(genereadPerCell),function(x) countsMatrix[,x]/genereadPerCell[x])
countMatrix.cpm=log1p(countMatrix.cpm)