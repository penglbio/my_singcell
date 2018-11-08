rm(list=ls())
library(RColorBrewer)
library(scater)
library(destiny)
library("ROCR")
library(M3Drop)
library(M3DExampleData)
library(SC3)
library("beanplot")
library(Rtsne)
library(biomaRt)
library(affycoretools)
library(GenomicRanges)
library(MASS)
library(DESeq)
library(gplots)
library(topGO)
library(cummeRbund)
argv <- commandArgs(T)
argv[2] <- "Mus_musculus.GRCm38.84.renameChr.bed"
argv[3] <- "htseq_reads_count.txt"
argv[4] <- "0.4-400000"
argv[5] <- "markers.txt"
argv[6] <- "/var/data/07.single_cell/figout/"
argv[7] <- "spearman"

cellinfo_table <- argv[1]
gene_ano <- argv[2]
ht_count <- argv[3]
align_ro_R<-as.numeric(strsplit(argv[4],"-")[[1]])
maker_gene <- argv[5]
figout <- argv[6]
mymethod <- argv[7]


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
cellinfo_C1 <- cellinfo[rownames(cellinfo)[grep("^D",rownames(cellinfo))],]
rownames(cellinfo_C1) <- sub("D_","BD_",sub("01|00$","",sub("_[0|1]_$","",gsub("\\.","_",rownames(cellinfo_C1)))))
cellinfo_C1.ok <- cellinfo_C1[colnames(countMatrix.cpm),]
dim(cellinfo_C1.ok)
doublet_id <- gsub("D_","BD_",sub("01|00$","",doublet_id))


mynames=substr(mysplit(filemap[rownames(cellinfo.ok)],"09-",2,2),1,39)
names(mynames)=gsub("[.][.]","",paste(gsub("BD","D",gsub("_",".",substr(mynames,1,9))), substr(mynames,38,39),sep=""))

###########-------------------------------loaddata-----------------------------------------------------------------------------
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
###load mapped logfile
maplog <- read.table("maplog.txt",header = T,sep="\t",stringsAsFactors = F,row.names = 1)
mapped_ratio <- 1-rowSums(apply(maplog[,c(7,8,9)],2,function(x){as.numeric(sub("%","",x))}))/100
sample_0.4 <- rownames(maplog[mapped_ratio<=0.4,])
sample_400000 <- rownames(maplog[maplog[,1]*mapped_ratio <=400000,])
exclude <- union(union(sample_400000,sample_0.4),doublet_id)
##ok sample logfile and counts
head(counts)
dim(counts)
##rm the last 3 ht_count "__too_low_aQual;__not_aligned;__alignment_not_unique"
counts<-counts[!grepl("_",rownames(counts)),]
countsMatrix <- counts[,!(colnames(counts) %in% exclude)]
write.table(countsMatrix,"Clean_sample_count_reads.txt",quote = F,sep="\t")
maplog <- maplog[!(colnames(counts) %in% exclude),]
###load marker genes
markers <- gsub("\\s","",readLines(maker_gene))
markername <- sapply(strsplit(markers,":"),function(x){x[1]})
markergene <- strsplit(sapply(strsplit(markers,":"),function(x){x[2]}),",")
names(markergene) <- markername
allmarkers <- markergene
myMarkerAll<- names(read.csv("myMakerALL.txt"))
mychoiceMarkerAll <- names(read.csv("mychoiceMakerAll.txt"))

##load pop counts
load("SSAll-STARHTSeqTags.-161121.Counts.RData")
dim(htseqgenes.pop)
genesPerCell.pop=apply(htseqgenes.pop,2,function(x) length(which(x>0)))
genereadPerCell.pop=apply(htseqgenes.pop,2,function(x) sum(x))/1000000
countmatrix.pop.cpm=sapply(1:length(genereadPerCell.pop),function(x) htseqgenes.pop[,x]/genereadPerCell.pop[x])
colnames(countmatrix.pop.cpm)=colnames(htseqgenes.pop)
countmatrix.pop.cpm <- countmatrix.pop.cpm[grep("EN",rownames(countmatrix.pop.cpm)),]
dim(countmatrix.pop.cpm)
#########-----------------------------------------------------------------analysis--------------------------------------------------------------------

##QC figures data
genesPerCell <- apply(countsMatrix,2,function(x) length(which(x>0)))
genereadPerCell <- apply(countsMatrix,2,function(x) sum(x))/1000000
readsPerCell <- as.numeric(apply(maplog,1,function(x) x[1]))/1000000
names(readsPerCell) <- names(genesPerCell)
countMatrix.cpm=sapply(1:length(genereadPerCell),function(x) countsMatrix[,x]/genereadPerCell[x])
log_countMatrix.cpm=log1p(countMatrix.cpm)
colnames(log_countMatrix.cpm) <-colnames(countMatrix.cpm) <- colnames(countsMatrix)
rownames(log_countMatrix.cpm) <-rownames(countMatrix.cpm) <- rownames(countsMatrix)
##QC figures
pdf("QC_figure.pdf")
#genePerCell
beanplot(lapply(c("P1","P2","P3"),function(x) genesPerCell[grep(x,names(genesPerCell))]),names=c("P1","P2","P3"),main="Nr. of genes/sample")
##genereadsPerCell
beanplot(lapply(c("P1","P2","P3"),function(x) genereadPerCell[grep(x,names(genereadPerCell))]),names=c("P1","P2","P3"),main="Genereads/Cell",ylab="Reads annotated in genes/cell(      )")
text(-0.32,2.95,expression(10^6),xpd=TRUE,srt=90)
beanplot(lapply(c("P1","P2","P3"),function(x) readsPerCell[grep(x,names(readsPerCell))]),names=c("P1","P2","P3"),main="Reads/Cell",ylab="Aligned reads/cell(       )")
text(-0.05,2.5,expression(10^6),xpd=TRUE,srt=90)
dev.off()

###count Maker genes scores
head(countMatrix.cpm)
choosecols=c("ivory3","orange", "red")
choosecolsgenes=c("steelblue3","yellow", "red")
##bulk colour
filcols=brewer.pal(9,"Pastel1")
batchcols=colnames(countMatrix.cpm); names(batchcols)=colnames(countMatrix.cpm)
batchcols[grep("P1",batchcols)]=filcols[1]; batchcols[grep("P2",batchcols)]=filcols[2]
batchcols[grep("P3",batchcols)]=filcols[3]
####count marker_score and plot the density figures(Extand file--figj)
names(markergene)
marker_score_list <- list()
marker_cols_list <- list()
for (i in  c(names(markergene))){
  # i <- "immuno"
  marker_score_list[[i]] <- marker_score(markergene[[i]],gene84.GR,choosecols,0,30)[[1]]
  marker_cols_list[[i]] <- marker_score(markergene[[i]],gene84.GR,choosecols,0,30)[[2]]
  }

beanplot(lapply(c("P1","P2","P3"),function(x){marker_score_list[["adipo"]][grep(x,names(marker_score_list[["adipo"]]))]}))
draw_marker(samp1="adipo",samp2="stem",cols=batchcols,fdir=figout,figname = "adipo_vs_stem.pdf",limx=50)


####Spearman corr.
rep_cpm<-sapply(c("P1","P2","P3"),function(x)apply(countMatrix.cpm[,grep(x,colnames(countMatrix.cpm))],1,function(y)mean(y)))
Spearman_cor<- function(a_cpm,b_cpm,cornames){
  # a_cpm <- rep_cpm
  # b_cpm <- countmatrix.pop.cpm[,19:22]  
  # cornames <- c("P1","P2","P3")
   allcors <- lapply(1:ncol(a_cpm),function(x) sapply(1:ncol(b_cpm), function(y)scatterSmoothPlot(log1p(a_cpm[,x]),log1p(b_cpm[,y]),x,y,"spearman")))
   names(allcors)=cornames
   return(allcors)
}

cor_list <- Spearman_cor(rep_cpm,countmatrix.pop.cpm[,19:22],c("P1","P2","P3"))

####Spearman figure
pdf(paste0(figout,"Correlation_SC_vs_Pop.pdf"))
beanplot(cor_list,main="SC vs. Pop",ylab="Spearman cor")
dev.off()
##############--------------------------------------M3Dplot--------------------












#---------------mysplit-----------------------------------------
mysplit <- function (vector,splitsign,numbers,position) {
  return(unlist(strsplit(as.character(vector),splitsign))[seq(position,length(vector)*numbers,by=numbers)])
}
#--------------------------------------------------------------

##--------------getEnsgID-------------------------
getEnsgID <- function(genes.GR,genenames) {
  return(names(genes.GR[which(elementMetadata(genes.GR)$geneName%in%genenames)]))
}
#---------------------------------------------------


#--------------------markerscore and plot------------------------------------------
##count score
marker_score <- function(marker_list,genes_anno,choosecols,mymin,mymax){
  temp=log1p(countMatrix.cpm[getEnsgID(genes_anno,marker_list),])
  if(!is.matrix(temp)==TRUE){
     marker_cols=marker_score=temp
   }else{
    marker_cols <- marker_score <- colSums(temp)
    names(marker_cols) <- colnames(countsMatrix)
   }
  marker_cols <- valuesToColorsAbs(marker_cols,choosecols,names(marker_cols),mymin,mymax)
  list(marker_score,marker_cols)
}
##plot figure
draw_marker <- function(samp1,samp2,cols,fdir,figname,limx) {
  pdf(paste0(figout,figname))
  plot(marker_score_list[[samp1]],marker_score_list[[samp2]],xlab=samp1,ylab=samp2,main=paste0(samp1,"vs.",samp1),pch=19,col=cols,cex=0.4)
  text(marker_score_list[[samp1]],marker_score_list[[samp2]]+max(marker_score_list[[samp2]])/100,labels=substr(names(marker_score_list[[samp1]]),4,8),col=cols,cex=0.5)
  plot(density(marker_score_list[[samp1]]),main=samp1,xlim=c(0,limx))
  plot(density(marker_score_list[[samp2]]),main=samp2,xlim=c(0,limx))
  dev.off()
}

###valuesToColorsAbs
valuesToColorsAbs=function (myvector,choosecols,mynames,mymin,mymax) {
  myvector=c(mymin,mymax,myvector); names(myvector)[1:2]=c("A1","A2"); mynames=c("A1","A2",mynames)
  myvector[which(is.na(as.numeric(as.vector(myvector))))]=min(myvector[which(!is.na(as.numeric(as.vector(myvector))))])
  myvector=round(myvector,digits=2)
  colormap=seq(min(myvector),max(myvector),by=0.01)
  colormap=colorRampPalette(choosecols)(length(colormap))
  names(colormap)=round(seq(min(myvector),max(myvector),by=0.01),digits=2)
  rpfmycol=as.vector(colormap[match(myvector,names(colormap))])
  names(rpfmycol)=mynames
  return(rpfmycol[which(!names(rpfmycol)%in%c("A1","A2"))])
}

valuesToColors=function (myvector,choosecols,mynames) {
  myvector[which(is.na(as.numeric(as.vector(myvector))))]=min(myvector[which(!is.na(as.numeric(as.vector(myvector))))])
  myvector=round(myvector,digits=2)
  colormap=seq(min(myvector),max(myvector),by=0.01)
  colormap=colorRampPalette(choosecols)(length(colormap))
  names(colormap)=round(seq(min(myvector),max(myvector),by=0.01),digits=2)
  rpfmycol=as.vector(colormap[match(myvector,names(colormap))])
  names(rpfmycol)=mynames
  return(rpfmycol)
}



#----------------------scatterSmoothPlot---------------------------------
scatterSmoothPlot <- function(a,b,namea,nameb,mymethod) {
  which((b%in%c(NA,"NaN","-Inf","Inf")))
  myind=which((!a%in%c(NA,"NaN","-Inf","Inf")))&which((!b%in%c(NA,"NaN","-Inf","Inf")))
  mycor=cor(a[myind],b[myind],method=mymethod)
  mytitle=paste(paste(namea,nameb,collapse="|",sep="|"),prettyNum(mycor,digits=4),sep=":",collapse=":")
  smoothScatter(a[myind],b[myind],main=mytitle,xlab=namea,ylab=nameb)
  return(mycor)
}
#----------------------makeprefheatmaprlog----------------------------
makeprefheatmaprlog = function(genetable,mygenes.GR,colvals="grey",mainname="row.scaled",
                               rowvals="grey") {
  genetable <- mymatrix.withlab3
  mymatrix=genetable
  if (length(grep("ENS",rownames(mymatrix)))>0) {
    rownames(mymatrix)[grep("ENS",rownames(mymatrix))]=elementMetadata(mygenes.GR[rownames(mymatrix)[grep("ENS",rownames(mymatrix))]])$geneName
  }
  rowvals=valuesToColors(apply(mymatrix,1,function(x) median(x)),c("ivory","darkred"),names(apply(mymatrix,1,function(x) median(x))))
  mymatrix=mymatrix[which(!duplicated(rownames(mymatrix))),]
  dim(mymatrix)
  heatmap.2(mymatrix,scale="none",col=rev(colorRampPalette(c("red","orange","ivory","lightblue","darkblue"))(216)),
            main=mainname,density.info="none",notecex=0.8,notecol="black",keysize = 1,trace="none",cexRow=0.6,cexCol=0.4,
            ColSideColors=colvals[colnames(mymatrix)],RowSideColors=rowvals[rownames(mymatrix)])
  heatmap.2(mymatrix,scale="row",col=rev(colorRampPalette(c("red","ivory","lightblue"))(128)),
            main=mainname,density.info="none",notecex=0.8,notecol="black",keysize = 1,trace="none",cexRow=0.6,cexCol=0.4,
            ColSideColors=colvals[colnames(mymatrix)],RowSideColors=rowvals[rownames(mymatrix)])
  #### prettier heatmap #####
  pheatmap(mymatrix,scale="none",clustering_distance_rows="correlation",cluster_cols = F,main=mainname,
           clustering_method="average",fontsize_col=4,fontsize_row=7,show_rownames=T,cexRow=0.1,
           cexCol=0.1,ColSideColors=colvals[colnames(mymatrix)],RowSideColors=rowvals[rownames(mymatrix)])
  pheatmap(mymatrix,clustering_distance_rows="correlation",cluster_cols = F,main=mainname,
           clustering_method="average",scale = "row",fontsize_col=4,fontsize_row=7,show_rownames=T,cexRow=0.1,
           cexCol=0.1,ColSideColors=colvals[colnames(mymatrix)],RowSideColors=rowvals[rownames(mymatrix)])
  
}


#---------------------------------------------------GO Enrichment--------------------------------------------
#runGO.84(myIDs,allgenes.background,plotpath,ensembl84,"BP","elimCount",myp,paste(mainname,"GO-pcut-",i,sep="_"),ontoGOs.84,showgraph)
runGO.84 <- function(genes,mybackground,path,ensemblversion,ont,algorithm,myp,set,showgraph) {
  ensembl=useEnsembl(biomart="ensembl", dataset="mmusculus_gene_ensembl", version=ensemblversion)
  mygenes <-getBM(attributes=c("ensembl_gene_id","external_gene_name","description"), filters="ensembl_gene_id", values=genes, mart=ensembl)
  gene2GO<-getBM(attributes=c("ensembl_gene_id","go_id"), mart=ensembl)
  mm10GO <- gene2GO[gene2GO[,2]!="",]
  head(gene2GO)
  mm10.gene2GO <- by(mm10GO$go_id,mm10GO$ensembl_gene_id,function(x)as.character(x))
  hRes <- sigGOTable(mygenes$ensembl_gene_id,names(mm10.gene2GO),mm10.gene2GO,ont,myp,algorithm,set,showgraph)
  write.table(hRes,paste(path,"GO/",set,"_topGO_",myp,".",ont,".",algorithm,".txt",sep=""),append=FALSE,sep="\t",quote=FALSE, eol="\n",na="NA",dec=".",row.names=TRUE,col.names=TRUE)
  
  myGO=hRes[,1]
  mygenes=cbind(mygenes,"")
  mygenes=cbind(mygenes,"")
  colnames(mygenes)[4:5]=c("GO.ID","GO.Description")
  mygenes=as.matrix(mygenes)
  for(i in 1:dim(hRes)[1]) {
    x=hRes[i,1]
    x.info=gene2GO[which(gene2GO[,2]%in%x),]
    mygenes[which(mygenes[,1]%in%x.info[,1]),4]=paste(mygenes[which(mygenes[,1]%in%x.info[,1]),4],hRes[i,1],sep=";")
    mygenes[which(mygenes[,1]%in%x.info[,1]),5]=paste(mygenes[which(mygenes[,1]%in%x.info[,1]),5],hRes[i,2],sep=";")
  }
  return(mygenes)
}

sigGOTable <- function(selGenes, GOgenes, gene2GO,ont, maxP,mymethod,set,showgraph){
  inGenes <- factor(as.integer(GOgenes %in% selGenes))
  names(inGenes) <- GOgenes
  GOdata <- new("topGOdata", ontology=ont, allGenes=inGenes,annot=annFUN.gene2GO, gene2GO=gene2GO)
  myTestStat <- new(mymethod, testStatistic=GOFisherTest,name="Fisher test", cutOff=maxP)
  mySigGroups <- getSigGroups(GOdata, myTestStat)
  sTab <- GenTable(GOdata, mySigGroups, topNodes=length(usedGO(GOdata)))
  names(sTab)[length(sTab)] <- "p.value"
  sTab <- subset(sTab, as.numeric(p.value) < maxP)
  sTab$Term <- sapply(mget(sTab$GO.ID, env=GOTERM), Term)
  if ((showgraph%in%c(4:10))&(dim(sTab)[1]>0)) {
    showSigOfNodes(GOdata, topGO::score(mySigGroups), firstSigNodes = showgraph, useInfo = "def")
    textplot( sTab, valign="top"  )
    title(paste(set,ont,sep="_"))
  }
  return(sTab)}

