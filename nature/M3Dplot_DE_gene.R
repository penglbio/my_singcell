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
argv[7] <- "/var/data/07.single_cell/rawdata_out/"
argv[8] <- "spearman"

cellinfo_table <- argv[1]
gene_ano <- argv[2]
ht_count <- argv[3]
align_ro_R<-as.numeric(strsplit(argv[4],"-")[[1]])
maker_gene <- argv[5]
figout <- argv[6]
rawdata_out<- argv[7]
mymethod <- argv[8]

setwd("/var/data/07.single_cell/mydata/")
conpath <- "/var/data/07.single_cell/Areg/data/"

###load clean reads count
countsMatrix <- read.table("Clean_sample_count_reads.txt",header = T)
head(countsMatrix)

###all marker gene--------------------------
markers <- gsub("\\s","",readLines(maker_gene))
markername <- sapply(strsplit(markers,":"),function(x){x[1]})
markergene <- strsplit(sapply(strsplit(markers,":"),function(x){x[2]}),",")
names(markergene) <- markername
plotMarker=unique(c("cdtdRFP","Dlk1","Ppara","Pdgfra","Pdgfrb","Cd24a","Zfp423","Adipoq","Wt1","Ebf2",unlist(markergene[3:10])))

###M3Dplot to identify the DEgene
Normalized_data <- M3DropCleanData(countsMatrix, is.counts=T, min_detected_genes=3000)
dim(Normalized_data$data)
fits <- M3DropDropoutModels(Normalized_data$data)
data.frame(MM=fits$MMFit$SAr, Logistic=fits$LogiFit$SAr,DoubleExpo=fits$ExpoFit$SAr) #MM best fit
data.frame(MM=fits$MMFit$SSr, Logistic=fits$LogiFit$SSr, DoubleExpo=fits$ExpoFit$SSr) #MM best fit
DE_genes <- M3DropDifferentialExpression(Normalized_data$data,mt_method="fdr", mt_threshold=0.05)
dim(DE_genes)

gene84.GR[rownames(DE_genes[which(rownames(DE_genes)%in%getEnsgID(gene84.GR,plotMarker)),])]### Fabp4 as only one among lineage genes
heat_out <- M3DropExpressionHeatmap(DE_genes$Gene, Normalized_data$data, batchcols)
###M3Drop get the marker gene
cell_populations_4 = M3DropGetHeatmapCellClusters(heat_out, k=4); cell_populations_3 = M3DropGetHeatmapCellClusters(heat_out, k=3)
marker_genes_4 = M3DropGetMarkers(Normalized_data$data, cell_populations_4); marker_genes_3 = M3DropGetMarkers(Normalized_data$data, cell_populations_3)

####SC3
mymatrix.log <- log1p(Normalized_data$data[rownames(DE_genes),]); 
mysceset <- SingleCellExperiment(assays = list(counts = Normalized_data$data[rownames(DE_genes),])); 
mysceset <- calculateQCMetrics(mysceset); 
rowData(mysceset)$feature_symbol <- rownames(mysceset)
logcounts(mysceset) <- mymatrix.log
mysceset <- sc3(mysceset, ks = 3:4,n_cores=1,biology = TRUE)
#### run SC3 on entire matrix
countsMatrix <- as.matrix(countsMatrix)
mysceset.all<-SingleCellExperiment(assays = list(counts = countsMatrix[rownames(Normalized_data$data),colnames(Normalized_data$data)])); 
mysceset.all <- calculateQCMetrics(mysceset.all);
rowData(mysceset.all)$feature_symbol <- rownames(mysceset.all)
logcounts(mysceset.all) <- log1p(Normalized_data$data)
mysceset.all <- sc3(mysceset.all, ks = 2:6,n_cores=1,biology = TRUE)
#scclusters.3=sc3_summarise_results(mysceset, k = 3); scclusters.4=sc3_summarise_results(mysceset, k = 4)
scclusters.3=colData(mysceset)[,"sc3_3_clusters"]; scclusters.4=colData(mysceset)[,"sc3_4_clusters"]
names(scclusters.3)=names(scclusters.4)=rownames(colData(mysceset))
scclusters.3.all=colData(mysceset.all)[,"sc3_3_clusters"]; scclusters.4.all=colData(mysceset.all)[,"sc3_4_clusters"]
names(scclusters.3.all)=names(scclusters.4.all)=rownames(colData(mysceset.all))
#### explore clustering #####
pdf(paste0(figout,"SC3_fig_est.pdf"))
sc3_plot_cluster_stability(mysceset, k = 3); sc3_plot_cluster_stability(mysceset, k = 4)
sc3_plot_silhouette(mysceset, k = 3); sc3_plot_silhouette(mysceset, k = 4)
dev.off()
mysceset.est <- sc3_estimate_k(mysceset); mysceset.est2 <- sc3_estimate_k(mysceset.all)
metadata(mysceset.est)$sc3$k_estimation 
sc3_plot_consensus( mysceset, k = 3, show_pdata = c("log10_total_features","sc3_3_clusters", "sc3_3_log2_outlier_score"))
sc3_plot_expression(mysceset, k = 3, show_pdata = c("cell_id", "log10_total_features","sc3_3_clusters", "sc3_3_log2_outlier_score"))
sc3_plot_cluster_stability(mysceset, k = 3);
sc3_plot_cluster_stability(mysceset, k = 4)
metadata(mysceset.all)


##### get the clusters - either from hierarchical clustering or from SC3 + get the markers for each cluster ####
heat.3=filcols[cell_populations_3]; names(heat.3)=names(cell_populations_3)
heat.4=filcols[cell_populations_4]; names(heat.4)=names(cell_populations_4)
sc3.3=filcols[scclusters.3]; names(sc3.3)=names(cell_populations_3); cID=scclusters.3;names(cID)=names(cell_populations_3)
sc3.4=filcols[scclusters.4]; names(sc3.4)=names(cell_populations_4)
marker_genes_sc_3 <- M3DropGetMarkers(Normalized_data$data, scclusters.3)
marker_genes_sc_4 <- M3DropGetMarkers(Normalized_data$data, scclusters.4)

##### cluster colors 
ccols=list(heat.3, heat.4, sc3.3,sc3.4); names(ccols)=c("heat.3", "heat.4","sc3.3","sc3.4")
##### Heatmap & SC3 markers, top 30 and top 50
P3.top=head(marker_genes_3[marker_genes_3$Group==1,],30)
P2.top=head(marker_genes_3[marker_genes_3$Group==2,],30)
P1.top=head(marker_genes_3[marker_genes_3$Group==3,],30)

### SC3 
P1.top.sc3=head(marker_genes_sc_3[marker_genes_sc_3$Group==1,],30); P3.top.sc3=head(marker_genes_sc_3[marker_genes_sc_3$Group==2,],30)
P2.top.sc3=head(marker_genes_sc_3[marker_genes_sc_3$Group==3,],30)

### top 50
P1.top.sc3.50=head(marker_genes_sc_3[marker_genes_sc_3$Group==1,],50); P3.top.sc3.50=head(marker_genes_sc_3[marker_genes_sc_3$Group==2,],50)
P2.top.sc3.50=head(marker_genes_sc_3[marker_genes_sc_3$Group==3,],50)
load("../Areg/data/GeneTypes84.R")
names(GENETYPES84)
ofinterest= unlist(c(GENETYPES84[[1]],GENETYPES84[[19]],GENETYPES84[[10]],GENETYPES84[[11]],GENETYPES84[[16]],GENETYPES84[[18]],GENETYPES84[[20]]))
P3.top.sc3.50=P3.top.sc3.50[which(rownames(P3.top.sc3.50)%in%ofinterest),]
P2.top.sc3.50=P2.top.sc3.50[which(rownames(P2.top.sc3.50)%in%ofinterest),]
P1.top.sc3.50=P1.top.sc3.50[which(rownames(P1.top.sc3.50)%in%ofinterest),]

#### Generate different matrices for visualisation
mymatrix.withlab.big=log1p(countMatrix.cpm)[,colnames(mymatrix.log)]
mymatrix.withlab=log1p(Normalized_data$data[c(rownames(P3.top),rownames(P2.top),rownames(P1.top)),])
mymatrix.withlab2=log1p(Normalized_data$data[c(rownames(P3.top.sc3),rownames(P2.top.sc3),rownames(P1.top.sc3)),])
mymatrix.withlab3=log1p(Normalized_data$data[c(rownames(P3.top.sc3.50),rownames(P2.top.sc3.50),rownames(P1.top.sc3.50)),])
rownames(mymatrix.withlab)[grep("ENS",rownames(mymatrix.withlab))]=elementMetadata(gene84.GR[rownames(mymatrix.withlab)[grep("ENS",rownames(mymatrix.withlab))]])$geneName
rownames(mymatrix.withlab2)[grep("ENS",rownames(mymatrix.withlab2))]=elementMetadata(gene84.GR[rownames(mymatrix.withlab2)[grep("ENS",rownames(mymatrix.withlab2))]])$geneName
rownames(mymatrix.withlab3)[grep("ENS",rownames(mymatrix.withlab3))]=elementMetadata(gene84.GR[rownames(mymatrix.withlab3)[grep("ENS",rownames(mymatrix.withlab3))]])$geneName
rownames(mymatrix.withlab.big)[grep("ENS",rownames(mymatrix.withlab.big))]=elementMetadata(gene84.GR[rownames(mymatrix.withlab.big)[grep("ENS",rownames(mymatrix.withlab.big))]])$geneName




###Plot - diffusion map and t-SNE (used in manuscript); color according to clusters & markers
##1.Diffusion map
Cluster_and_marker_plot <- function(dimre_out,dimre_name){
  toplot=mymatrix.withlab.big
  pdf(paste0(figout,dimre_name,"marker_Cat.pdf"))
  for (i in 1:length(marker_cols_list)) {
    plot(dimre_out, pch = 19, col = marker_cols_list[[i]][colnames(mymatrix.log)],main =paste("Cats - ",names(marker_cols_list)[i]))
}
  dev.off()
  pdf(paste0(figout,dimre_name,"Clusters.pdf"))
  for (i in 1:length(ccols)) 
  {
    plot(dimre_out, pch = 19, col = ccols[[i]],main =paste("Clusters - ",names(ccols)[i]))
  }
  dev.off()
  mymax=max(mymatrix.withlab.big)*2/3; plotMarker=unique(c(myMarkerAll,mychoiceMarkerAll,unlist(allmarkers),rownames(mymatrix.withlab3)))
  plotMarker=sort(plotMarker[which(plotMarker%in%rownames(toplot))])
  choosecollist=list(c("ivory2","darkgrey", "black"),c("ivory2","orange", "red"),c("ivory2","lightblue", "blue"),c("ivory2","lightgreen", "darkgreen"))
  pdf(paste0(figout,dimre_name,"Markergene.pdf"))
   for (j in 1:length(choosecollist)) {
    for (i in 1:length(plotMarker)) {
      choosecols=choosecollist[[j]]; plotcol=valuesToColorsAbs(toplot[plotMarker[i],],choosecols,colnames(mymatrix.withlab),0,mymax)
      plot(dimre_out, pch = 19, col = plotcol,main =plotMarker[i])
    }}
  dev.off()
}

plotDiffusionMap(mysceset)
dif_out <- DiffusionMap(t(mymatrix.log))
plot(eigenvalues(dif_out), ylim = c(0,0.3), pch = 20, xlab = "Diffusion component (DC)", ylab = "Eigenvalue")
Cluster_and_marker_plot(dif_out,"dif")

##2.t-SNE
rtsne_out <- Rtsne(as.matrix(t(mymatrix.log)),theta=0.0000001,perplexity=as.integer(dim(mymatrix.log)[2]/6))
Cluster_and_marker_plot(rtsne_out$Y,"rtsne")

####GO annotation
P1.top.sc3.10=head(marker_genes_sc_3[marker_genes_sc_3$Group==1,],10); P3.top.sc3.10=head(marker_genes_sc_3[marker_genes_sc_3$Group==2,],10)
P2.top.sc3.10=head(marker_genes_sc_3[marker_genes_sc_3$Group==3,],10)
P1.top.sc3.TFs=head(marker_genes_sc_3[marker_genes_sc_3$Group==1,],100); P3.top.sc3.TFs=head(marker_genes_sc_3[marker_genes_sc_3$Group==2,],100)
P2.top.sc3.TFs=head(marker_genes_sc_3[marker_genes_sc_3$Group==3,],100)
names(GENETYPES84)

ofinterest= unlist(c(GENETYPES84[[1]],GENETYPES84[[20]]))
P3.top.sc3.TFs=P3.top.sc3.TFs[which(rownames(P3.top.sc3.TFs)%in%ofinterest),]
P2.top.sc3.TFs=P2.top.sc3.TFs[which(rownames(P2.top.sc3.TFs)%in%ofinterest),]
P1.top.sc3.TFs=P1.top.sc3.TFs[which(rownames(P1.top.sc3.TFs)%in%ofinterest),]
mymatrix.withlabTFs=log1p(Normalized_data$data[unique(c(rownames(P3.top.sc3.10),rownames(P2.top.sc3.10),rownames(P1.top.sc3.10),
                                                        rownames(P3.top.sc3.TFs),rownames(P2.top.sc3.TFs),rownames(P1.top.sc3.TFs))),])
rownames(mymatrix.withlabTFs)[grep("ENS",rownames(mymatrix.withlabTFs))]=elementMetadata(spikedgenes84.GR[rownames(mymatrix.withlabTFs)[grep("ENS",rownames(mymatrix.withlabTFs))]])$geneName

###cluster heatmap include top10 and top100 TFs(fig1.c)
distance = as.dist(1-cor(t(mymatrix.withlabTFs),method="spearman")); cluster = hclust(distance, method = "ward.D2")
distance.Col = dist(t(mymatrix.withlabTFs)); cluster.Col = hclust(distance, method = "ward.D2")
rowvals=valuesToColors(apply(mymatrix.withlabTFs,1,function(x) median(x)),c("ivory","darkred"),names(apply(mymatrix.withlabTFs,1,function(x) median(x))))
heatmap.2(mymatrix.withlabTFs[,names(sort(cID))],scale="none",col=rev(colorRampPalette(c("red","orange","ivory","lightblue","darkblue"))(214)),
          Colv = FALSE, Rowv = as.dendrogram(cluster), density.info="none",notecex=0.8,notecol="black",keysize = 1,trace="none",cexRow=0.8,cexCol=0.4,
          ColSideColors=ccols[[3]][names(sort(cID))],RowSideColors=rowvals[rownames(mymatrix.withlabTFs)],dendrogram ="row")
heatmap.2(mymatrix.withlabTFs[,names(sort(cID))],scale="none",col=rev(colorRampPalette(c("tomato4","sienna2","ivory","steelblue1","steelblue3"))(214)),
          Colv = FALSE, Rowv = as.dendrogram(cluster), density.info="none",notecex=0.8,notecol="black",keysize = 1,trace="none",cexRow=0.8,cexCol=0.4,
          ColSideColors=ccols[[3]][names(sort(cID))],RowSideColors=rowvals[rownames(mymatrix.withlabTFs)],dendrogram ="row")
heatmap.2(mymatrix.withlabTFs[,names(sort(cID))],scale="row",col=rev(colorRampPalette(c("tomato4","sienna2","ivory","steelblue1","steelblue3"))(214)),
          Colv = FALSE, Rowv = as.dendrogram(cluster), density.info="none",notecex=0.8,notecol="black",keysize = 1,trace="none",cexRow=0.8,cexCol=0.4,
          ColSideColors=ccols[[3]][names(sort(cID))],RowSideColors=rowvals[rownames(mymatrix.withlabTFs)],dendrogram ="row")
write.table(mymatrix.withlabTFs[,names(sort(cID))],paste(rawdata_out,"Fig1c-raw.C1Heatmap.txt",sep=""),sep="\t",row.names=T,col.names=T,quote=F)

####Extend marker gene FigS1.e

Markers_CMP_Plot <- function(markers_in_plot,barcols){
  temp=log1p(countMatrix.cpm[getEnsgID(gene84.GR,markers_in_plot),])
  rownames(temp)[grep("ENS",rownames(temp))]=elementMetadata(gene84.GR[rownames(temp)[grep("ENS",rownames(temp))]])$geneName
  temp=temp[markers_in_plot,]
  pdf(paste0(figout,"barplot_marker_cpm.pdf"))
  boxplot(lapply(1:6,function(x) temp[x,]),names=rownames(temp), col=barcols,las=2,ylab=c("Log. CPM"))
  write.table(t(temp),paste(rawdata_out,"FigS1e-raw.C1Marker.txt",sep=""),sep="\t",row.names=T,col.names=T,quote=F)
  dev.off()
  }
Markers_CMP_Plot(markers_in_plot=c("Actb","Ptprc","Pecam1","Cd34","Ly6a","Itgb1"),barcols=c("snow3","black","black","salmon2","salmon2","salmon2"))


####fig.Extend.k Fabp4 with marker genes
#### Correlations with Fabp4 and CD34
###figS1.k
adipo_genes <- c("Pdgfra","Pdgfrb","Pparg","Zfp423","Cd24a","Dlk1","Prdm16","Ppargc1a","Ebf2","Adipoq","Wt1")

key_cor_adipo_genes <- function(key,adipo_genes){
  mycors.p=c(sapply(adipo_genes,function(x) cor.test(log1p(countMatrix.cpm[getEnsgID(gene84.GR,key),]),log1p(countMatrix.cpm[getEnsgID(gene84.GR,c(x)),]),method="spearman")$p.val))
  mycors=c(sapply(adipo_genes,function(x) cor(log1p(countMatrix.cpm[getEnsgID(gene84.GR,key),]),log1p(countMatrix.cpm[getEnsgID(gene84.GR,c(x)),]),method="spearman")))
  cols=rep("grey",length(mycors.p))
  cols[which(mycors.p<0.01)]="black"
  names(cols)=names(mycors.p)
  pdf(paste0(figout,key,"_adipo_cortable.pdf"))
  barplot(sort(-log10(mycors.p),decreasing=T),main=paste0("-logp.Cor.",key),las=2,col=cols[names(sort(-log10(mycors.p),decreasing = T))])
  barplot(mycors[names(sort(-log10(mycors.p),decreasing=T))],main=paste0("Cor.",key),las=2,col=cols[names(mycors[names(sort(-log10(mycors.p),decreasing = T))])])
  dev.off()
  cor.table=cbind(mycors.p,mycors)
  write.table(cor.table,file=paste0(rawdata_out,key,"_adipo_corTable.txt"),row.names=T,col.names=T,sep="\t",quote=F)
}

key_cor_adipo_genes("Fabp4",adipo_genes)
key_cor_adipo_genes("Cd34",adipo_genes)
###Adipogenic, stemness in 3 cluster etc. scor
##### Adipogenic, stemness etc. score FigS1.i,S3.a
Catscore <- function(Cat){
  if(Cat %in% names(markergene)){
    Catscore <- marker_score(markergene[[Cat]],gene84.GR,choosecols,0,30)[[1]] 
  }else
  {
    Catscore <- marker_score(Cat,gene84.GR,choosecols,0,30)[[1]]
    
  }
  Cat.score <- lapply(c(1,2,3),function(x) Catscore[names(cID[which(cID%in%x)])])
  Cat.pval.c1=c(wilcox.test(Cat.score[[1]],Cat.score[[2]])$p.val,wilcox.test(Cat.score[[1]],Cat.score[[3]])$p.val,
                  wilcox.test(Cat.score[[2]],Cat.score[[3]])$p.val)
  names(Cat.pval.c1)=c("P1-P2","P1-P3","P2-P3")
  beancols=list(c("darkgreen",rep("black",3)),c("salmon3",rep("black",3)),c("blue",rep("black",3)))
  boxcols=c("darkgreen","red","blue")
  pdf(paste0(figout,Cat,"_score.pdf"))
  beanplot(Cat.score,col=beancols,main=Cat, bw="nrd0" )
  boxplot(Cat.score,col=boxcols,main=Cat)
  dev.off()
  write.table(Cat.pval.c1,paste(rawdata_out,"FigS1i-raw.",Cat,"ScorePvals.txt",sep=""),sep="\t",row.names=T,col.names=T,quote=F)
  return(Cat.score)
  }

adiposcore.cats=Catscore("adipo")
osteoscore.cats=Catscore("osteo")
endoscore.cats=Catscore("endo")
immunoscore.cats=Catscore("immuno")
clu_markerscore=list()
for(i in c("Cd55","Il13ra1","F3","Adam12","Aoc3","Adipoq","Retn","Cd34","Ly6a","Itgb1","Cd44","Il13ra1","Il33","Il13ra2","Postn","Fap","Col1a1","S100a4")){
  clu_markerscore[[i]]=Catscore(i)
}

######## positive vs. negative (data not shown in manuscript); Dlk1/RFP analysis
par(mfrow=c(3,6))
PvsNscore.cats=sapply(c(1,3,2),
        function(x) c(length(names(cID[which(cID%in%x)])[which(names(cID[which(cID%in%x)])%in%colnames(countMatrix.cpm[,which(cellinfo_C1.ok[,"Cell"]%in%"R")]))]),
                    length(names(cID[which(cID%in%x)])[which(names(cID[which(cID%in%x)])%in%colnames(countMatrix.cpm[,which(cellinfo_C1.ok[,"Cell"]%in%"N")]))])))
#### make beanplot with the expression of different genes per RED vs. SNOW
lapply(c("Fabp4","Pdgfra","Pdgfrb","Pparg","Zfp423","Cd24a","Dlk1","Prdm16","Ppargc1a","Ebf2","Adipoq","Wt1"),
       function(y) boxplot(lapply(c("R","N"),
             function(x) log1p(countMatrix.cpm[getEnsgID(gene84.GR,y),
             which(cellinfo_C1.ok[,"Cell"]%in%x)])),col=c("red","snow3"),names=c("R","N"),main=y))
### and RFP
sepvals=lapply(adipo_genes,function(y) lapply(c("R","N"),function(x) as.numeric(log1p(countMatrix.cpm[getEnsgID(spikedgenes84.GR,y)[1],
                                                                                                      which(cellinfo_C1.ok[,"Cell"]%in%x)]))))
pvals=sapply(sepvals,function(x) wilcox.test(x[[1]],x[[2]], alternative="greater")$p.val)
snames=sapply(adipo_genes,function(y) getEnsgID(gene84.GR,y))
names(pvals)=elementMetadata(gene84.GR[snames])$geneName
pvals=sapply(pvals,function(x) min(x,1-x))
barplot(apply(PvsN.cats,2,function(x)x/sum(x)),beside=F,names=c("P1","P2","P3"),col=c("red","snow3"),las=2)

###GO enrichment
###### run the Gene Ontology analysis on genes specific for one population (not included in .man)--BP-Ensembl84
allgenes.background=rownames(Normalized_data$data)
bigcut=0.000006
P1.top.sc3.200=head(marker_genes_sc_3[marker_genes_sc_3$Group==1,],200);#P1.top.sc3.GO= P1.top.sc3.GO[which( P1.top.sc3.GO[,3]<bigcut),]
P2.top.sc3.200=head(marker_genes_sc_3[marker_genes_sc_3$Group==2,],200);#P3.top.sc3.GO= P3.top.sc3.GO[which( P3.top.sc3.GO[,3]<bigcut),]
P3.top.sc3.200=head(marker_genes_sc_3[marker_genes_sc_3$Group==3,],200);#P2.top.sc3.GO= P2.top.sc3.GO[which( P2.top.sc3.GO[,3]<bigcut),]
golist.sign=list(P1.top.sc3.200,P2.top.sc3.200,P3.top.sc3.200)
names(golist.sign)=c("P1","P2","P3")
SC3.golist=list(NA)
myp=0.001;showgraph=10; 
mainname="ASC"
###keypoint do twice;run twice;run again
for (i in 1:length(golist.sign)){
  myIDs=rownames(golist.sign[[i]])#names(spikedgenes84.GR[which(elementMetadata(spikedgenes84.GR)$geneName%in%rownames(golist.sign[[i]][which(golist.sign[[i]][,2]>1),]))])
  if (length(myIDs)>10) {
    SC3.golist[[i]]=runGO.84(myIDs,allgenes.background,rawdata_out,84,"BP","elimCount",myp,paste(mainname,"GO-pcut-",i,sep="_"),showgraph)
  }
}





