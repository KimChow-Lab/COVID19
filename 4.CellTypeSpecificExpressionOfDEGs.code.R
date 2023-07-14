library(Seurat)
library(harmony)
library(dplyr)
library(ggplot2)
library(presto)
library(tibble)
library(fgsea)
####Dealing with the scRNAseq dataset of mouse brain #####
setwd("/data2/deng/MouseCortexReference/Bhattacherjee_GSE124952")

#https://www.nature.com/articles/s41467-019-12054-3
GeneExpr=read.csv("rawData/GSE124952_expression_matrix.csv",row.names=1,check.names=F)
metaData=read.csv("rawData/GSE124952_meta_data.csv",row.names=1)
all(rownames(metaData)==colnames(GeneExpr))
MouseCortex <- CreateSeuratObject(counts = GeneExpr, project ="MousePFC")
MouseCortex <- NormalizeData(MouseCortex, normalization.method = "LogNormalize", scale.factor = 10000)
MouseCortex <- FindVariableFeatures(MouseCortex, selection.method = "vst", nfeatures = 2000)
MouseCortex <- ScaleData(MouseCortex)
MouseCortex <- RunPCA(MouseCortex, features = VariableFeatures(object = MouseCortex))
MouseCortex <- RunHarmony(MouseCortex,group.by.vars="orig.ident",reduction.save = "harmony")
MouseCortex<- RunUMAP(MouseCortex, reduction = "harmony", dims = 1:50)
MouseCortex<- FindNeighbors(MouseCortex, reduction = "harmony", dims = 1:50)
MouseCortex<- FindClusters(MouseCortex, resolution = 0.5)

table(MouseCortex$seurat_clusters,MouseCortex$cellType)
#check the relationship between our clusters and the cell type provided by meta.data

classicMarker=c("Rbfox1","Gad1","Mag","Vcan","Gfap","P2ry12","Flt1")
pdf("MarkerGene.pdf",height=10,width=10)
FeaturePlot(MouseCortex, features=classicMarker,raster=TRUE)&NoLegend()&NoAxes()&theme(panel.border = element_rect(fill=NA,color="black", linewidth=1.5, linetype="solid"))
dev.off()


Idents(MouseCortex)=MouseCortex$seurat_clusters
new.cluster.ids <- c("Excitatory","Excitatory","Excitatory","Endo","Excitatory","Excitatory","Endo","Oligo","Astro","Excitatory","Microglia","Excitatory","Inhibitory","OPC","NF Oligo","Inhibitory","Inhibitory","Excitatory","Microglia","Endo")
names(new.cluster.ids) <- levels(MouseCortex)
MouseCortex <- RenameIdents(MouseCortex, new.cluster.ids)
MouseCortex$MyCellType=Idents(MouseCortex)
mycolors <- colorRampPalette(brewer.pal(8, "Dark2"))(length(unique(MouseCortex$MyCellType)))

all(rownames(MouseCortex@meta.data)==rownames(metaData)) #TRUE
MouseCortex$cellType=metaData$CellType

pdf("MouseCortexCluster.pdf")
DimPlot(MouseCortex,label=T)&NoLegend()&NoAxes()
dev.off()

#fig.1 H
pdf("MouseCortexCellType.pdf")
DimPlot(MouseCortex,label=T,group.by="cellType")&NoLegend()&NoAxes()&theme(plot.title=element_blank())
dev.off()

saveRDS(MouseCortex,file="MouseCortexFromGSE124952.rds")



#### cell type specific expression of the DEGs #####

setwd("/data2/deng/COVID19B2/DEG")
H1N1DEG=read.table("Sig_DEG_D30_H1N1vsCtrl.txt",header=T)
BA2DEG=read.table("Sig_DEG_D30_BA2vsCtrl.txt",header=T)
COVIDWTDEG=read.table("Sig_DEG_D30_COVIDWTvsCtrl.txt",header=T)
sigDEG=rbind(H1N1DEG,BA2DEG,COVIDWTDEG)
dim(sigDEG)
#2578
geneList=intersect(sigDEG$Symbol,rownames(MouseCortex))
sigDEG=sigDEG[sigDEG$Symbol%in%geneList,]
dim(sigDEG)
#1902
sigDEGList=split(sigDEG$Symbol,paste0(sigDEG$Pattern,"_",sigDEG$Infection,sep=""))

fgsea_sets=sigDEGList
DefaultAssay(MouseCortex)="RNA"
Healthy.integratedMarker <- wilcoxauc(MouseCortex, 'cellType')
table(Healthy.integratedMarker$group)
for(cluster in unique(MouseCortex$cellType)){
print (cluster)
clusterCell<- Healthy.integratedMarker %>% dplyr::filter(group == cluster) %>% arrange(desc(logFC)) %>% dplyr::select(feature, logFC)
ranks<- deframe(clusterCell)
fgseaRes<- fgseaMultilevel(fgsea_sets, stats = ranks,eps=0, nPermSimple = 10000)
ranks=na.omit(ranks)
fwrite(fgseaRes, file=paste0("cellTypeExpForSigDEG/",cluster,".txt",sep=""), sep="\t", sep2=c("", " ", ""))
}

#fig.1 I
data=read.table("DEG/cellTypeExpForSigDEG.txt",header=T,sep="\t") #cellTypeExpForSigDEG.txt is a combine file of the cell types under cellTypeExpForSigDEG folder
data=data[order(data$NES,decreasing=T),]
data=data[,1:8]
data$sig=ifelse(data$pval<0.05,"Sig","NonSig")
data=data[order(data$pval),]
pN=table(data[data$pval<0.05&data$NES>0,"pathway"])
PathwayList=factor(data$pathway,levels=names(sort(pN)))
t=ggplot(data,aes(Cluster,pathway,size=-1*log(pval),colour=NES,shape=sig))+geom_point()+
scale_color_gradient2(low="blue",mid="white",high = "red")+
scale_shape_manual(values=c(3,19))+
theme_bw()+#theme(legend.position="bottom")+
theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x = element_text(size=rel(1.0),colour = "black",angle=90, vjust = 0.5, hjust=1),axis.text.y = element_text(size=rel(1.0),colour = "black"))
pdf("DEG/cellTypeExpForSigDEG.pdf",width=6,height=4)
print(t)
dev.off()
