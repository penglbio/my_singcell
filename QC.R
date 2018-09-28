##QC figures
#genePerCell
beanplot(lapply(c("P1","P2","P3"),function(x) genesPerCell[grep(x,names(genesPerCell))]),names=c("P1","P2","P3"),main="Nr. of genes/sample")
##genereadsPerCell
beanplot(lapply(c("P1","P2","P3"),function(x) genereadPerCell[grep(x,names(genereadPerCell))]),names=c("P1","P2","P3"),main="Genereads/Cell",ylab="Reads annotated in genes/cell(      )")
text(-0.32,2.95,expression(10^6),xpd=TRUE,srt=90)
beanplot(lapply(c("P1","P2","P3"),function(x) readsPerCell[grep(x,names(readsPerCell))]),names=c("P1","P2","P3"),main="Reads/Cell",ylab="Aligned reads/cell(       )")
text(-0.05,2.5,expression(10^6),xpd=TRUE,srt=90)