library(beanplot)
library(yaml)
##load data
argv <- commandArgs(T)
config_file <- argv[1]
cg_mat_file <- argv[2]
config_file<-config.yaml
# setwd("/var/data/03.jsk_SC_20181010/count_mat/")
# cg_mat_file<-"cell_gene_mat.tsv"
# config_file<-"config.yaml"
config <- yaml.load_file(config_file)

####UMINCell
cell_gene_mat<-read.table(cg_mat_file,header = T)
dim(cell_gene_mat)
UMI_hum_mouse<-cbind(human=apply(cell_gene_mat[grep("^ENSG",rownames(cell_gene_mat)),],2,sum),
                     mouse=apply(cell_gene_mat[grep("^ENSMUSG",rownames(cell_gene_mat)),],2,sum))

pdf("figure/UMI_Hg_vs_MM.pdf")
plot(UMI_hum_mouse,ylab="umiCounts/Cell",main=paste0(config$sample,"UmiN/Cell"),ylim=c(0,10000))
dev.off()

###GeneNCell
GeneNCell<-apply(cell_gene_mat,2,function(x){length(x[x>0])})
Sorted_GeneNCell<-sort(GeneNCell,decreasing = T)
pdf("figure/GeneNCell.pdf")
par(xpd=T)
beanplot(Sorted_GeneNCell,ylab="number of Gene/Cell",main=paste0(config$sample,"GeN/Cell"))
dev.off()

##Hg_vs_MM cells bar
h_id<-rownames(UMI_hum_mouse[UMI_hum_mouse[,1]>UMI_hum_mouse[,2],])
m_id<-rownames(UMI_hum_mouse[UMI_hum_mouse[,1]<UMI_hum_mouse[,2],])

pdf("figure/Hg_vs_MM_cell.pdf")
par(mai=c(1,3,1,3))
barplot(c(human=length(Human_cell),mouse=length(Mouse_cell)),ylab="Number of Cell",col=c("green","red"),width = 0.1,ylim = c(0,5000),main=paste0(config$sample,"Cells/h_v_m"))
dev.off()
###gene and umi in cell (Hg and MM separeted)
hce<-cell_gene_mat[,h_id]
mce<-cell_gene_mat[,m_id]
Hg_gn_cell<-apply(hce,2,function(x){length(x[x>0])})
MM_gn_cell<-apply(mce,2,function(x){length(x[x>0])})
Human_cell_umi<-UMI_hum_mouse[h_id,]
Mouse_cell_umi<-UMI_hum_mouse[m_id,]

pdf("figure/Hg_vs_MM_sp_gene_umi.pdf")
beanplot(list(human=Human_cell_umi,mouse=Mouse_cell_umi),ylab="UMIscounts/Cell",main=paste0(config$sample,"UmiN/Cell_h_v_m"))
beanplot(list(human=Hg_gn_cell,mouse=MM_gn_cell),ylab="Genescounts/Cell",main=paste0(config$sample,"GN/Cell_h_v_m"))
dev.off()



###
library("rJava")
library("xlsx")
df<-list(UMI_hum_mouse,Sorted_GeneNCell,h_id,m_id,Hg_gn_cell,MM_gn_cell,Human_umi_cell,Mouse_umi_cell)
df_name<-c("UMI_hum_mouse","Sorted_GeneNCell","h_id","m_id","Hg_gn_cell","MM_gn_cell","Human_umi_cell","Mouse_umi_cell")
wb<-createWorkbook()
for (i in 1:length(df)){
  sheet<-createSheet(wb,sheetName = as.character(df_name[i]) )
  addDataFrame(df[[i]],sheet,row.names = TRUE)
}
saveWorkbook(wb,"stat/cell_gene_umi_stat.xlsx")





